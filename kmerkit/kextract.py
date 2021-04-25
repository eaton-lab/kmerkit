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
import gzip
import subprocess
from loguru import logger
from kmerkit.kmctools import KMTBIN
from kmerkit.utils import KmerkitError, get_fastq_dict_from_path
from kmerkit.kschema import Project, KextractData, KextractParams, KextractBase

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
    """
    def __init__(self, json_file, samples=None, min_kmers_per_read=1, keep_paired=True):

        # load project
        self.json_file = json_file
        self.project = Project.parse_file(json_file).dict()

        # type-check params
        self._params = KextractParams(
            min_kmers_per_read=min_kmers_per_read,
            keep_paired=True,
        )
        self.params = self._params.dict()

        # output prefix
        self.prefix = os.path.join(
            self.project['workdir'], f"{self.project['name']}_kextract")

        # get input data from database samples AND newly entered samples.
        self.factq_dict = {}
        self.select_samples(samples)


    def select_samples(self, samples):
        """
        Checks that db_samples are in kcount.data and that fastq_dict
        samples are not repeated in db_samples. ...
        """
        # ...
        init_set = set(self.project['kinit']['data'])
        input_set = set(samples)

        samples = set()
        for iset in input_set:
            # special keyword to select kfilter group 0
            if iset == "0":
                fsamps = set(self.project['kfilter']['params']['trait_0'])
                for fsamp in fsamps:
                    samples.update(
                        [str(i) for i in self.project['kinit']['data'][fsamp]]
                    )

            # special keyword to select kfilter group 1
            elif iset == "1":
                fsamps = set(self.project['kfilter']['params']['trait_1'])
                for fsamp in fsamps:
                    samples.update(
                        [str(i) for i in self.project['kinit']['data'][fsamp]]
                    )

            # a new filepath
            elif os.path.exists(iset):
                samples.update([str(iset)])
            # a sample name from the init_set
            else:
                if iset in init_set:
                    samples.update([str(self.project['kinit']['data'][iset])])
                else:
                    logger.warning(
                        f"Skipping sample {iset} not present in project database")

        self.fastq_dict = get_fastq_dict_from_path(None, samples, "_R")



    def get_reads_with_kmers(self, fastq, sname, readnum):
        """
        Generate a directory full of fastq files for each sample
        where reads are only kept if they contain kmers from the
        specified kmer database.

        CMD: kmc_tools filter database input.fastq -ci10 -cx100 out.fastq
        """
        # Here broken into [kmc_tools filter database <options>]
        cmd = [
            KMTBIN, "-hp", "-t8",
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
        cmd.extend([self.prefix + f"_{sname}_R{readnum}.fastq"])

        # log and call the kmc_tool command
        logger.debug(" ".join(cmd))
        subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            # cwd=self.workdir,
        )



    def match_paired_reads(self, sname):
        """
        Get read pairs for all reads in which EITHER has a kmer match.
        Current approach may be memory crushing... TODO: write chunk 
        files.
        """
        logger.debug(f"pair matching in {sname}")

        # get fastq file paths of original data
        old_fastqs = [str(i) for i in self.fastq_dict[sname]]
        new_fastqs = [
            f"{self.prefix}_{sname}_R1.fastq",
            f"{self.prefix}_{sname}_R2.fastq",
        ]

        # -----------------------------------------------------
        # create set of all read names in new kmer-matched file
        set1 = fastq_to_read_names(new_fastqs[0])
        lines1 = fastq_to_line_nos(old_fastqs[0], set1)
        del set1

        set2 = fastq_to_read_names(new_fastqs[1])
        lines2 = fastq_to_line_nos(old_fastqs[1], set2)
        del set2

        # -------------------------------------------------------
        # get union of lines1 and line2
        lidxs = lines1.union(lines2)
        del lines1
        del lines2

        # skip the rest if no lidxs exist
        if not lidxs:
            return 0, (None, None)

        # overwrite new read files with reads from original files
        # at all line indices in lidxs.
        for new, old in zip(new_fastqs, old_fastqs):

            # open writer 
            with open(new, 'w') as out:

                # load the originals
                readio = (
                    gzip.open(old, 'rt') if old.endswith('.gz') 
                    else open(old, 'r')
                )
                # read in 4 lines at a time
                matched = iter(readio)
                quart = zip(matched, matched, matched, matched)

                # save each 4-line chunk to chunks if it lidxs
                chunk = []
                idx = 0

                # iterate until end of file
                while 1:
                    try:
                        test = next(quart)
                        if idx in lidxs:
                            chunk.extend(test)
                    except StopIteration:
                        break
                    idx += 1

                    # occasionally write to disk and clear
                    if len(chunk) == 2000:
                        out.write("".join(chunk))
                        chunk = []

                if chunk:
                    out.write("".join(chunk))
                readio.close()

        # return the number of paired reads
        return len(lidxs), new_fastqs



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



    def run(self, force=False):
        """
        Iterate over all fastq files to call filter funcs.
        """
        # check for current step
        if not force:
            self.check_overwrite()

        stats = {}
        for sname in self.fastq_dict:

            # get fastq filename
            fastqs = self.fastq_dict[sname]

            # accommodates paired reads:
            for readnum, fastq in enumerate(fastqs):

                # call kmc_tools filter, writes new fastq to prefix
                self.get_reads_with_kmers(fastq, sname, readnum + 1)

            # get read pairs where either contains the kmer
            if self.params['keep_paired']:
                nreads, data_out = self.match_paired_reads(sname)
            # else:
                # nreads = self.concat_reads()

            # store results
            stats[sname] = KextractData(
                data_in=fastqs,
                data_out=data_out,
                kmer_matched_reads=nreads
            )
            logger.debug(stats[sname].json(indent=4))

        # save to project
        self.project['kextract'] = KextractBase(
            params=self._params,
            data=stats,
        )
        with open(self.json_file, 'w') as out:
            out.write(Project(**self.project).json(indent=4))


def fastq_to_read_names(fastq):
    """
    Extracts the fastq read header from all reads as a set.
    Used in match_paired_reads()
    """
    name_set = set()
    with open(fastq, 'r') as readio:
        matched = iter(readio)
        quart = zip(matched, matched, matched, matched)
        while 1:
            try:
                header = next(quart)[0].strip()
                name_set.add(header)
            except StopIteration:
                break
    return name_set


def fastq_to_line_nos(fastq, name_set):
    """
    Extract line nos of headers matching to a target name set.
    Used in match_paired_reads()    
    """
    # find lines with same read names in orig fastq file
    line_nos = set()
    readio = (
        gzip.open(fastq, 'rt') if fastq.endswith('.gz') else open(fastq, 'r')
    )
    matched = iter(readio)
    quart = zip(matched, matched, matched, matched)
    idx = 0
    while 1:
        try:
            header = next(quart)[0].strip()
            if header in name_set:
                line_nos.add(idx)
                idx += 1
        except StopIteration:
            break
    readio.close()
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
