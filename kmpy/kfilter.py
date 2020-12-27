#!/usr/bin/env python

"""
Class with functions for filtering fastq files based on a set of target kmers.

TODO: 
    - implement filters on number of matching kmers, or trimming, etc.
"""


import os
import sys
import glob
import gzip
import subprocess
import pandas as pd
from loguru import logger
from kmpy.utils import KmpyError


# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


class Kfilter:
    """
    Filter fastq reads based on a set of target kmers and write to new files.

    Finds existing database using the CSV in 'workdir/name'.
    Group samples in database based on phenos info (dict or csv)
    and calls union function to get shared kmer counts.

    Parameters
    ==========
    name (str):
        Prefix name for the fastq files that will be written. 
        Example: 'kfiltered' --> <workdir>/<name>_<sample_name>.fastq
    workdir (str):
        Working directory where new filtered fastq files will be written.
        Examples: '/tmp' or '/tmp/newfastqs'.
    fastq_path (str):
        A wildcard selector to match one or more fastq files (can be .gz).
        Examples: './data/*.fastq.gz'
    group_kmers (str):
        The prefix path for a kmc binary database with kmers.
        Examples: '/tmp/test_group'

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
        self.prefix = os.path.join(self.workdir, f"kfilter_{self.name}")

        # get the kmctools binary for kmer set comparisons
        self.kmctools_binary = os.path.join(sys.prefix, "bin", "kmc_tools")
        logger.debug("using KMC binary: {}".format(self.kmctools_binary))

        # attributes to be filled
        self.files = []
        self.sample_names = []
        self.names_to_files = {}
        self.files_to_names = {}        

        # check files
        self.check_database()
        self.expand_filenames()
        self.check_samples()

        # dataframe for results
        self.statsdf = pd.DataFrame(
            index=self.sample_names,
            columns=[
                "total_reads", "kmer_matched_reads", 
                "new_fastq_path", "orig_fastq_path",
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
                msg = f"no fastq files found at: {self.fastq_path}"
                logger.error(msg)
                raise KmpyError(msg)

            # sort the input files
            self.files = sorted(self.files)

        # report on found files
        logger.debug("found {} input files".format(len(self.files)))



    def check_samples(self):
        """
        Gets sample names from the input files and checks that all 
        have the same style of suffix (e.g., .fastq.gz).
        """
        # split file names to keep what comes before 'name_split'
        sample_names = [
            os.path.basename(i.split(self.name_split)[0]) for i in self.files
        ]

        # store dictionaries for mapping names and files
        self.names_to_files = dict(zip(sample_names, self.files))
        self.files_to_names = dict(zip(self.files, sample_names))

        # check that all sample_names are unique
        assert len(set(sample_names)) == len(sample_names), "sample names are not unique."

        # get alphanumeric sorted names
        self.sample_names = sorted(self.names_to_files.keys())        



    def get_fastqs_containing_kmers(self, fastq, sname):
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
        cmd.extend([self.prefix + f"_{sname}.fastq"])

        # log and call the kmc_tool command
        logger.debug(" ".join(cmd))
        subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )


    def run(self):
        """
        Iterate over all fastq files to call filter funcs.
        """
        for sname in self.sample_names:

            # get fastq filename
            fastq = self.names_to_files[sname]

            # call kmc_tools filter on fastq to match kmers
            self.get_fastqs_containing_kmers(fastq, sname)

            # store file paths
            self.statsdf.loc[sname, "orig_fastq_path"] = fastq
            self.statsdf.loc[sname, "new_fastq_path"] = (
                self.prefix + f"_{sname}.fastq"
            )

            # count numer of matched kmers
            with open(self.statsdf.at[sname, "new_fastq_path"], 'r') as indat:
                nmatched = sum(1 for i in indat)
                self.statsdf.loc[sname, "kmer_matched_reads"] = nmatched

            # count stats and report to logger
            if fastq.endswith('.gz'):               
                with gzip.open(fastq, 'r') as indat:
                    nreads = sum(1 for i in indat) / 4
            else:
                with gzip.open(fastq, 'r') as indat:
                    nreads = sum(1 for i in indat) / 4                    
            self.statsdf.loc[sname, "total_reads"] = nreads
            with open(self.prefix + f"_{sname}.fastq", 'r') as indat:
                nreads = sum(1 for i in indat) / 4
                self.statsdf.loc[sname, "kmer_matched_reads"] = nreads

            # logger report
            logger.debug(f"found {nmatched} matching reads in sample {sname}")



if __name__ == "__main__":

    # first run: python3 kcount.py; python3 kgroup.py
    # DATA
    FASTQS = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R1_.fastq.gz"

    # set up filter tool
    kfilt = Kfilter(
        name="test",
        workdir="/tmp",
        fastq_path=FASTQS,
        group_kmers="/tmp/kgroup_test",
        name_split="_R",
        mindepth=5,
    )
    kfilt.run()
    print(kfilt.statsdf.T)
