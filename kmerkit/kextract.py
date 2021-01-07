#!/usr/bin/env python

"""
Kcounts -> Kgroup --[ kmatrix -> Kgwas ]--> Kextract

Extract fastq reads containing target kmers to produce new 
fastq files with subset of matching read(pair)s. The fastq
files that this is applied to do not need to have been 
processed earlier by kcount or any other tool. The target kmers
are identified in kgroup and/or kgwas.

TODO: option to not filter invariant kmers in Kgroup? Since excluding
these may limit our ability to construct contigs? Not sure about this...

"""

import os
import sys
import glob
import gzip
import subprocess
import pandas as pd
from loguru import logger
from kmerkit.utils import KmerkitError


# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


class Kextract:
    """
    Extract fastq reads containing target kmers and write to new files.

    Files 

    Parameters
    ==========
    name (str):
        Prefix name for the fastq files that will be written. 
        Example: <workdir>/k_extract_<name>_<sample_name>.fastq
    workdir (str):
        Working directory where new filtered fastq files will be written.
        Examples: '/tmp' or '/tmp/newfastqs'.
    fastq_path (str):
        A wildcard selector to match one or more fastq files (can be .gz).
        Examples: './data/*.fastq.gz'
    group_kmers (str):
        The prefix path for a kmc binary database with kmers.
        Examples: '/tmp/test_group'
    name_split (str):
        String on which to split fastq file names to extract sample 
        names as the first item. Common values are "_" or ".fastq". 
    mindepth (int):


    Attributes
    ===========
    statsdf (pandas.DataFrame):
        Statistics on the number of reads per sample. The complete CSV
        is written to <workdir>/<name>.csv.

    """
    def __init__(self, name, workdir, fastq_path, group_kmers, name_split="_", mindepth=1):

        # store parameters
        self.name = name
        self.group_kmers = group_kmers
        self.fastq_path = os.path.realpath(os.path.expanduser(fastq_path))
        self.name_split = name_split
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.mindepth = mindepth
        
        # output prefix
        os.makedirs(self.workdir, exist_ok=True)
        self.prefix = os.path.join(self.workdir, f"kextract_{self.name}")

        # get the kmctools binary for kmer set comparisons
        self.kmctools_binary = os.path.join(sys.prefix, "bin", "kmc_tools")
        logger.debug("using KMC binary: {}".format(self.kmctools_binary))

        # attributes to be filled
        self.files = []
        self.sample_names = []
        self.names_to_infiles = {}

        # check files
        self.check_database()
        self.expand_filenames()
        self.check_samples()

        # dataframe for results
        self.statsdf = pd.DataFrame(
            index=self.sample_names,
            columns=[
                "total_reads",
                "kmer_matched_reads", 
                "new_fastq_path", 
                "orig_fastq_path",
            ],
            dtype=int,
            data=0,
        )



    def check_database(self):
        """
        Check that a kmers database exists for this group name.
        """
        assert os.path.exists(self.group_kmers + ".kmc_suf"), (
            f"group_kmers database {self.group_kmers} not found"
        )
        assert os.path.exists(self.group_kmers + ".kmc_pre"), (
            f"group_kmers database {self.group_kmers} not found"
        )



    def expand_filenames(self):
        """
        Allows for selecting multiple input files using wildcard
        operators like "./fastqs/*.fastq.gz" to select all fastq.gz
        files in the folder fastqs.
        """
        if isinstance(self.fastq_path, (str, bytes)):
            self.files = glob.glob(self.fastq_path)

            # raise an exception if no files were found
            if not any(self.files):
                msg = f"no fastq files found at: {self.files}"
                logger.error(msg)
                raise KmerkitError(msg)

            # sort the input files
            self.files = sorted(self.files)

        # report on found files
        logger.debug(f"found {len(self.files)} input files")



    def check_samples(self):
        """
        Gets sample names from the input files and checks that all 
        have the same style of suffix (e.g., .fastq.gz).
        """
        # split file names to keep what comes before 'name_split'
        sample_names = [
            os.path.basename(i.split(self.name_split)[0]) for i in self.files
        ]

        # check that all sample_names are unique
        if len(set(sample_names)) != len(sample_names):
            
            # if not, then check each occurs 2X (PE reads)
            if not all([sample_names.count(i) == 2 for i in sample_names]):
                raise KmerkitError(
                    "Sample names are not unique, or in sets of 2 (PE). "
                    "You may need to try a different name_split setting."
                )                
            logger.debug("detected PE data")

        # store dict mapping names and files (or file pairs for PE)
        for sname, file in zip(sample_names, self.files):

            # names to input fastqs
            if sname in self.names_to_infiles:
                self.names_to_infiles[sname].append(file)
            else:
                self.names_to_infiles[sname] = [file]
 


    def get_reads_with_kmers(self, fastq, sname, readnum):
        """
        Generate a directory full of fastq files for each sample
        where reads are only kept if they contain kmers from the
        specified kmer database.

        CMD: kmc_tools filter database input.fastq -ci10 -cx100 out.fastq
        """
        # Here broken into [kmc_tools filter database <options>]
        cmd = [self.kmctools_binary, "filter", self.group_kmers]

        # insert options to filter on kmers:
        # -ci<val> : exclude kmers occurring < val times.
        # -cx<val> : exclude kmers occurring > val times.
        cmd.extend([f'-ci{self.mindepth}'])

        # add the input.fastq
        cmd.extend([fastq])

        # insert options to filter on reads:
        # -ci<val> : exclude reads containing < val kmers.
        # -cx<val> : exclude reads containing > val kmers
        cmd.extend(['-ci1'])

        # add the output fastq path
        cmd.extend([self.prefix + f"_{sname}_R{readnum}.fastq"])

        # log and call the kmc_tool command
        logger.debug(" ".join(cmd))
        subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )



    def get_paired_reads(self):
        """
        Get read pairs for all reads in which EITHER has a kmer match.
        Current approach may be memory crushing...
        """
        for sname in self.statsdf.index:
            logger.debug(f"pair fix {sname}")
            # get fastq file paths of original data
            orig_fastqs = self.statsdf.at[sname, 'orig_fastq_path']
            old1, old2 = orig_fastqs.split(",")

            # get fastq file paths of new filtered data
            new_fastqs = self.statsdf.at[sname, 'new_fastq_path']
            new1, new2 = new_fastqs.split(",")

            # -----------------------------------------------------
            # create set of all read names in new kmer-matched file
            set1 = set()
            readio = open(new1, 'r')
            matched = iter(readio)
            quart = zip(matched, matched, matched, matched)
            while 1:
                try:
                    header = next(quart)[0].strip()
                    set1.add(header)
                except StopIteration:
                    break
            readio.close()

            # find lines with same read names in orig fastq file
            lines1 = set()
            readio = (
                gzip.open(old1, 'rt') if old1.endswith('.gz') 
                else open(old1, 'r')
            )
            matched = iter(readio)
            quart = zip(matched, matched, matched, matched)
            idx = 0
            while 1:
                try:
                    header = next(quart)[0].strip()
                    if header in set1:
                        lines1.add(idx)
                    idx += 1
                except StopIteration:
                    break
            readio.close()
            del set1

            # ---------------------------------------------------
            # create set of all read names in read2s
            set2 = set()
            readio = open(new2, 'r')
            matched = iter(readio)
            quart = zip(matched, matched, matched, matched)
            while 1:
                try:
                    header = next(quart)[0].strip()
                    set2.add(header)
                except StopIteration:
                    break
            readio.close()

            # find lines containing read names in matched file
            lines2 = set()
            readio = (
                gzip.open(old2, 'rt') if old2.endswith('.gz') 
                else open(old2, 'r')
            )
            matched = iter(readio)
            quart = zip(matched, matched, matched, matched)
            idx = 0
            while 1:
                try:
                    header = next(quart)[0].strip()
                    if header in set2:
                        lines2.add(idx)
                    idx += 1
                except StopIteration:
                    break
            readio.close()
            del set2

            # -------------------------------------------------------
            # get union of lines1 and line2
            lidxs = lines1.union(lines2)
            del lines1
            del lines2

            # skip the rest if no lidxs exist
            if not lidxs:
                continue

            # overwrite new read files with reads from original files
            # at all line indices in lidxs.
            for new, old in [(new1, old1), (new2, old2)]:

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



    def run(self):
        """
        Iterate over all fastq files to call filter funcs.
        """
        for sname in self.names_to_infiles:        

            # get fastq filename
            fastqs = self.names_to_infiles[sname]

            # accommodates paired reads:
            for readnum, fastq in enumerate(fastqs):

                # call kmc_tools filter, writes new fastq to prefix
                self.get_reads_with_kmers(fastq, sname, readnum + 1)

            # store file paths
            self.statsdf.loc[sname, "orig_fastq_path"] = ",".join([
                i for i in fastqs if i])
            self.statsdf.loc[sname, "new_fastq_path"] = ",".join([
                self.prefix + f"_{sname}_R{readnum + 1}.fastq"
                for readnum, fastq in enumerate(fastqs)
            ])

        # get read pairs where either contains the kmer
        self.get_paired_reads()

        # get final stats
        #self.count_kmer_reads()



    def count_kmer_reads(self):
        """
        Test file is not fetched from statsdf... 
        
        """
        for sname in self.names_to_infiles:

            # get path to the test fastq file
            # ...

            # count number of matched kmers
            with open(self.statsdf.at[sname, "new_fastq_path"], 'r') as indat:
                nmatched = sum(1 for i in indat)
                self.statsdf.loc[sname, "kmer_matched_reads"] = nmatched

            # get nreads in the kmer-matched file
            with open(self.prefix + f"_{sname}.fastq", 'r') as indat:
                nreads = sum(1 for i in indat) / 4
                self.statsdf.loc[sname, "kmer_matched_reads"] = nreads

            # get nreads in the test fastq file
            if fastq.endswith('.gz'):               
                with gzip.open(fastq, 'r') as indat:
                    nreads = sum(1 for i in indat) / 4
            else:
                with gzip.open(fastq, 'r') as indat:
                    nreads = sum(1 for i in indat) / 4                    
            self.statsdf.loc[sname, "total_reads"] = nreads


            # logger report
            logger.debug(f"found {nmatched} matching reads in sample {sname}")



if __name__ == "__main__":

    # first run: python3 kcount.py; python3 kgroup.py
    # DATA
    # FASTQS = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R*_.fastq.gz"
    FASTQS = "~/Documents/kmerkit/data/hybridus_*.fastq.gz"

    import kmerkit
    kmerkit.set_loglevel("DEBUG")

    # set up filter tool
    kfilt = Kextract(
        name="test",
        workdir="/tmp",
        fastq_path=FASTQS,
        group_kmers="/tmp/kgroup_test",
        name_split="_R",
        mindepth=10,
    )
    kfilt.run()
    print(kfilt.statsdf.T)
