#!/usr/bin/env python

"""
Kcounts -> Kfilter --[ kmatrix -> Kgwas ]--> Kextract

Extract fastq reads containing target kmers to produce new 
fastq files with subset of matching read(pair)s. The fastq
files that this is applied to do not need to have been 
processed earlier by kcount or any other tool. The target kmers
are identified in kfilter, ktree, or kgwas...

TODO: option to not filter invariant kmers in Kgroup? Since excluding
these may limit our ability to construct contigs? Not sure about this...

CLI usage:
----------
kmerkit kextract \
    --json /tmp/test.json \
    --samples ... \
    --min-kmers-per-read 5

"""

import os
import re
import gzip
import time
import subprocess
import concurrent.futures
import numpy as np
from loguru import logger
from kmerkit.kmctools import KMTBIN
from kmerkit.utils import KmerkitError, num_cpus
from kmerkit.kschema import Project, KextractData, KextractParams, KextractBase
# from kmerkit.parallel import Cluster

# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


class Kextract:
    """
    Extract fastq reads containing target kmers and write to new files.

    Files 

    Parameters
    ==========
    json_file (str):
        ...
    db_samples (list):
        ...
    fastq_dict (dict):
        ...
    pairs_union (bool):
        Keep pairs as long as one or the other (the union) contains the 
        target min number of kmers. If False then only one OR the other
        must contain min kmers (the intersection).
    """
    def __init__(self, json_file, samples=None, min_kmers_per_read=1, paired_union=True):

        # load project
        self.json_file = json_file
        self.project = Project.parse_file(json_file).dict()

        # type-check params
        self._params = KextractParams(
            min_kmers_per_read=min_kmers_per_read,
            paired_union=True,
        )
        self.params = self._params.dict()

        # output prefix
        self.prefix = os.path.join(
            self.project['workdir'], f"{self.project['name']}_kextract")

        # get input data from database samples AND newly entered samples.
        self.fastq_dict = {}
        self.select_samples(samples)


    def select_samples(self, samples):
        """
        Checks that db_samples are in kcount.data and that fastq_dict
        samples are not repeated in db_samples. ...
        """
        # ...
        init_set = set(self.project['kinit']['data'])
        input_set = set(samples)

        # ----------------------------
        # if input is None then use database ktrim or kinit data
        if not input_set:
            if self.project.get('ktrim'):
                self.fastq_dict = {
                    sname: self.project['ktrim']['data'][sname]['data_out']
                    for sname in self.project['ktrim']['data']
                }
            else:
                self.fastq_dict = self.project['kinit']['data']
            return

        # ----------------------------
        # select samples by group (0,1), filepath, or sname/regex, allowing
        # multiple of these options to be combined. Builds a list of sample
        # names and then uses utils.get_fastq_dict() to get {name: [path, path]}
        samples = set()
        for iset in sorted(input_set):

            # special keyword to select kfilter group 0 or 1
            if iset in ["0", "1"]:
                fsamps = set(self.project['kfilter']['params'][f'trait_{iset}'])
                for fsamp in fsamps:
                    if self.project.get("ktrim"):
                        ifiles = [str(i) for i in self.project['ktrim']['data'][fsamp]['data_out']]
                    else:
                        ifiles = [str(i) for i in self.project['kinit']['data'][fsamp]]
                    self.fastq_dict[fsamp] = ifiles

            # a sample name in the database, or regex pattern of sample names.
            else:
                # get name and files from database
                if iset in init_set:
                    if self.project.get('ktrim'):
                        self.fastq_dict[iset] = (
                            self.project['ktrim']['data'][iset]['data_out']
                        )
                    else:
                        self.fastq_dict[iset] = self.project['kinit']['data'][iset]

                # check if name is a regex and expand to get names
                else:
                    search = (re.match(iset, i) for i in init_set)
                    matches = [i.string for i in search if i]
                    if matches:
                        logger.debug(f"{iset} expanded to {matches}")
                        for sname in matches:
                            if self.project.get('ktrim'):
                                self.fastq_dict[sname] = (
                                    self.project['ktrim']['data'][sname]['data_out']
                                )
                            else:
                                self.fastq_dict[sname] = self.project['kinit']['data'][sname]
            logger.info(f"{len(self.fastq_dict)} samples selected.")
            logger.debug(f"FASTQ_DICT: {self.fastq_dict}")


    def get_reads_with_kmers(self, fastq, sname, readnum, threads):
        """
        Generate a directory full of fastq files for each sample
        where reads are only kept if they contain kmers from the
        specified kmer database.

        CMD: kmc_tools filter database input.fastq -ci10 -cx100 out.fastq
        """
        logger.info(f"extracting reads: {fastq}")
        # Here broken into [kmc_tools filter database <options>]
        cmd = [
            KMTBIN, "-hp", "-t{}".format(threads),
            "filter", 
            str(self.project['kfilter']['data']['database_passed']),
        ]

        # insert options to filter on kmers:
        # -ci<val> : exclude kmers occurring < val times.
        # -cx<val> : exclude kmers occurring > val times.
        # cmd.extend([f'-ci{self.params['min_depth}'])

        # add the input.fastq
        cmd.extend([str(fastq)])

        # insert options to filter on reads:
        # -ci<val> : exclude reads containing < val kmers.
        # -cx<val> : exclude reads containing > val kmers
        cmd.extend([f"-ci{self.params['min_kmers_per_read']}"])

        # add the output fastq path
        cmd.extend([self.prefix + f"_{sname}_R{readnum}_tmp.fastq"])

        # log and call the kmc_tool command
        logger.debug(" ".join(cmd))
        subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            # cwd=self.workdir,
        )


    def check_overwrite(self):
        """
        Warn user of overwriting.
        """
        if self.project['kextract']:
            logger.error(
                "\nKextract results exist, use force to overwrite, or consider "
                "using branching to produce new results on a separate named "
                "branch without overwriting previous results."
            )
            raise KmerkitError("Preventing data overwrite")        



    def run(self, force=False, threads=None, workers=None):
        """
        Iterate over all fastq files to call filter funcs.
        """
        # check for current step
        if not force:
            self.check_overwrite()

        # set cores values to limit njobs to ncores / 4
        if workers in [0, None]:
            workers = max(1, int(np.ceil(num_cpus() / 4)))
            threads = 4

        # if user set workers, then scale threads to match
        else:
            workers = int(workers)
            threads = (threads if threads else int(num_cpus() / workers))
            # (threads if threads else 4)
        logger.debug(f"workers={workers}; threads={threads};")

        # wrapped execution
        stats = {}
        writing = {}
        with concurrent.futures.ProcessPoolExecutor(workers) as lbview:
            for sname in self.fastq_dict:
                args = (self.fastq_dict[sname][0], sname, 1, threads)
                future_1 = lbview.submit(self.get_reads_with_kmers, *args)
                args = (self.fastq_dict[sname][1], sname, 2, threads)
                future_2 = lbview.submit(self.get_reads_with_kmers, *args)
                writing[sname] = [future_1, future_2]

            # track jobs, start follower jobs, and get results
            matching = {}
            while 1:

                # get finished pair writing jobs and send pair-matching job
                finished = [
                    i for i in writing if all(j.done() for j in writing[i])
                ]
                for sname in finished:
                    if 1:  # TODO: support SE data.
                        matching[sname] = lbview.submit(
                            # self.match_paired_reads, sname
                            self.new_match_paired_reads, sname, self.params["paired_union"],
                        )
                    # else:
                        # lbview.submit(self.concat_reads, sname)
                    writing.pop(sname)

                # get finished pair-matching jobs and enter to stats
                finished = [i for i in matching if matching[i].done()]
                for sname in finished:
                    nreads, data_out = matching[sname].result()
                    # if data_out is not None:
                    stats[sname] = KextractData(
                        data_in=self.fastq_dict[sname],
                        data_out=data_out,
                        kmer_matched_reads=nreads
                    )
                    logger.debug(stats[sname].json(indent=4))
                    matching.pop(sname)

                # end looping
                if (not writing) and (not matching):
                    break
                time.sleep(1)

        # save to project
        self.project['kextract'] = KextractBase(
            params=self._params,
            data=stats,
        )
        with open(self.json_file, 'w') as out:
            out.write(Project(**self.project).json(indent=4))



    def new_match_paired_reads(self, sname, union=True):
        """
        Get read pairs for all reads in which EITHER has a kmer match.
        Current approach may be memory crushing... TODO: write chunk 
        files.
        """
        logger.info(f"re-matching paired reads in {sname}")

        # get fastq file paths of original data
        old_fastqs = [str(i) for i in self.fastq_dict[sname]]
        new_fastqs = [
            f"{self.prefix}_{sname}_R1_tmp.fastq",
            f"{self.prefix}_{sname}_R2_tmp.fastq",
        ]

        # -----------------------------------------------------
        # create set of all read names in new kmer-matched file
        lines1 = get_line_nos(old_fastqs[0], new_fastqs[0])
        lines2 = get_line_nos(old_fastqs[1], new_fastqs[1])

        # bail out if no data
        if (not lines1) or (not lines2):
            return 0, None

        # -------------------------------------------------------
        # get union of lines1 and line2
        if union:
            lidxs = lines1.union(lines2)
            logger.debug(
                f"{sname} union of PE reads: "
                f"{len(lines1)} U {len(lines2)} = {len(lidxs)}"
            )
        else:
            lidxs = lines1.intersection(lines2)
        del lines1
        del lines2

        # skip the rest if no lidxs exist
        if not lidxs:
            return 0, None

        # write lines to new files
        fname1 = new_fastqs[0].replace("_tmp.fastq", ".fastq.gz")
        fname2 = new_fastqs[1].replace("_tmp.fastq", ".fastq.gz")
        out1 = gzip.open(fname1, 'wt')
        out2 = gzip.open(fname2, 'wt')
        old1 = (
            gzip.open(old_fastqs[0], 'rt') if old_fastqs[0].endswith('.gz') 
            else open(old_fastqs[0], 'rt')
        )
        old2 =  (
            gzip.open(old_fastqs[1], 'rt') if old_fastqs[1].endswith('.gz') 
            else open(old_fastqs[1], 'rt')
        )
        
        # read chunks of 4 lines
        quart1 = zip(old1, old1, old1, old1)
        quart2 = zip(old2, old2, old2, old2)

        # save each 4-line chunk to chunks if it lidxs
        chunks1 = []
        chunks2 = []        
        idx = 0

        # iterate until end of file
        while 1:
            try:
                ch_r1 = next(quart1)
                ch_r2 = next(quart2)
                if idx in lidxs:
                    chunks1.extend(ch_r1)
                    chunks2.extend(ch_r2)
            except StopIteration:
                break
            idx += 1

            # occasionally write to disk and clear
            if len(chunks1) == 10000:
                out1.write("".join(chunks1))
                out2.write("".join(chunks2))
                chunks1 = []
                chunks2 = []

        if chunks1:
            out1.write("".join(chunks1))
            out2.write("".join(chunks2))            
        out1.close()
        out2.close()
        old1.close()
        old2.close()                        

        # return the number of paired reads
        return len(lidxs), (fname1, fname2)



def get_line_nos(fastq_big, fastq_small):
    """
    Get line-nos ...
    """
    # store lines from big that are in small
    line_nos = set()

    # get 4-line iterator
    bi_data = (
        gzip.open(fastq_big, 'rt') if fastq_big.endswith('.gz') 
        else open(fastq_big, 'r')
    )
    quart_big = zip(bi_data, bi_data, bi_data, bi_data)

    # get 4-line iterator
    si_data = (
        gzip.open(fastq_small, 'rt') if fastq_small.endswith('.gz') 
        else open(fastq_small, 'r')
    )
    quart_small = zip(si_data, si_data, si_data, si_data)

    # get first chunk (UNLESS file is empty!)
    try:
        header_big = next(quart_big)[0].strip()
        header_small = next(quart_small)[0].strip()
    except StopIteration:
        return None

    idx = 0
    while 1:
        try:
            # advance both if they match, else only the big one
            if header_big == header_small:
                header_big = next(quart_big)[0].strip()                
                header_small = next(quart_small)[0].strip()
                line_nos.add(idx)
            else:
                header_big = next(quart_big)[0].strip()
            idx += 1
        except StopIteration:
            break
    return line_nos    



if __name__ == "__main__":

    import kmerkit
    kmerkit.set_loglevel("DEBUG")

    # DATA
    # FASTQS = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R*_.fastq.gz"
    JSON = "/tmp/test.json"
    FASTQS = "~/Documents/kmerkit/data/amaranths/hybridus_*.fastq.gz"
    # kmerkit.utils.get_fastq_dict_from_

    # set up filter tool
    kfilt = Kextract(
        json_file=JSON,
        # samples=["hybridus_SLH_AL_1060"],
        # database_samples=[],
        # fastq_dict=None, #FASTQS,
        min_kmers_per_read=10,
    )
    kfilt.run()
    # print(kfilt.statsdf.T)
