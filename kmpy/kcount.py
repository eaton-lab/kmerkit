#!/usr/bin/env python

"""
The main class object for executing KMC kmer funcs
KMC GitHub: https://github.com/refresh-bio/KMC
KMC Paper: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399
KMC Suppl: https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/33/17/10.1093_bioinformatics_btx304/3/btx304_supplementary.pdf?Expires=1611337056&Signature=4iJG6giNcZdiDuHwljf-SbELTll74FtIj3YIFvfESeZC~m39EZPJdSXfqAJStvCr5SmH9lHGRCdJGHBLseX~ZunAgFZBFFHikmODBI14Kq84ctkQMihTvBzU1rme~S6MpXcC1Erxavl~ckAEnE7jfwIRJbtm4bSkTk-sEcZKIHqR3H0SwhdN0zMmhqMFkwn~jvNo5Rd~yPFwq8aXtE2CBrMORgVUsu~ACFKnl7sWB2FLtsZ2zp~ENuVz28mYwZWkyFXnUlRq2sKHenjWdw4BChI~QDf5EULM2oXgx4dSDlLTaIaJjsZBGl8tKQGw5Ohz48YEbqO82pn14AyeJkaRsQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA
KMC Docs: https://github.com/refresh-bio/KMC/blob/master/API.pdf
KMC-tools Docs: https://github.com/refresh-bio/KMC/blob/master/kmc_tools.pdf

-hp = hide progress
-t x = option to limit threads
-cs  maximal value of a counter... //???

"""

import os
import sys
import glob
import subprocess
import pandas as pd
from loguru import logger
from kmpy.utils import KmpyError


# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


class Kcount:
    """
    Calls the kmc counting functions to create a database.

    Parameters:
    ===========
    name (str):
        ...
    workdir (str):
        ...
    fastq_path (str):
        ...
    kmersize (int):
        ...
    name_split (str):
        ...
    """
    def __init__(self, name, workdir, fastq_path, kmersize, name_split="_"):

        # store input params
        self.name = name
        self.fastq_path = os.path.realpath(os.path.expanduser(fastq_path))
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.kmersize = kmersize
        self.name_split = name_split
        self.prefix = os.path.join(self.workdir, f"kcount_{self.name}")

        # attributes to be filled
        self.files = []
        self.names_to_infiles = {}
        self.infiles_to_names = {}
        self.names_to_outfiles = {}
        self.statsdf = None

        # get kmc binary from conda 
        self.kmc_binary = os.path.join(sys.prefix, "bin", "kmc")
        logger.info("using KMC binary: {}".format(self.kmc_binary))

        # expands wildcard operators to get filenames from string
        self.expand_filenames()

        # constructs statsdf and names_to_files
        self.check_samples()

        # dataframe with sample names and column names
        self.statsdf = pd.DataFrame(
            index=sorted(self.names_to_infiles),
            columns=[
                'kmers_total', 
                'kmers_unique',
                'kmers_unique_counted',
                'reads_total', 
                'kmers_below_thresh', 
                'kmers_above_thresh',
                'database_path',
                'fastq_path',
            ],
            data=0,
            dtype=int,
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
                msg = "no fastq files found at: {self.files}"
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
        self.names_to_infiles = dict(zip(sample_names, self.files))
        self.infiles_to_names = dict(zip(self.files, sample_names))
        self.names_to_outfiles = dict(
            zip(sample_names, [f"{self.prefix}_{i}" for i in sample_names])
        )

        # check that all sample_names are unique
        assert len(set(sample_names)) == len(sample_names), "sample names are not unique."



    def run(self):
        """
        Calls kmc count on all files
        """
        for filename in self.files:
            self.call_kmc_count(filename)



    def call_kmc_count(self, filename):
        """
        takes a fastq file and calls KMC count on it.
        """
        # get sample name
        sname = self.infiles_to_names[filename]

        # create command: 'kmc -k17 inputfastq outname workdir'
        cmd = [
            self.kmc_binary, "-k{}".format(self.kmersize), 
            filename,
            self.names_to_outfiles[sname],
            self.workdir,
        ]

        # call subprocess on the command
        out = subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )     

        # parse STDOUT to get the kmer stats
        stats = out.stdout.decode().split("Stats:")[-1]

        # store results to the dataframe
        for line in stats.strip().split("\n"):
            key, val = line.split(" : ")
            
            if key.strip() == "Total no. of k-mers":
                self.statsdf.loc[sname, "kmers_total"] = int(val.strip())

            if key.strip() == "No. of unique k-mers":
                self.statsdf.loc[sname, "kmers_unique"] = int(val.strip())

            if key.strip() == "No. of unique counted k-mers":
                self.statsdf.loc[sname, "kmers_unique_counted"] = int(val.strip())

            if key.strip() == "Total no. of reads":
                self.statsdf.loc[sname, "reads_total"] = int(val.strip())

            if key.strip() == "No. of unique k-mers above max. threshold":
                self.statsdf.loc[sname, "kmers_above_thresh"] = int(val.strip())

            if key.strip() == "No. of unique k-mers below min. threshold":
                self.statsdf.loc[sname, "kmers_below_thresh"] = int(val.strip())

        # save fastq filepath
        self.statsdf.loc[sname, "database_path"] = self.names_to_outfiles[sname]
        self.statsdf.loc[sname, "fastq_path"] = filename

        # save database file to the workdir
        path = os.path.join(self.prefix + ".csv")
        self.statsdf.to_csv(path)
        logger.info(f"database updated [{sname}] ({path})")



if __name__ == "__main__":


    # test dataset
    FILES = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R1_.fastq.gz"

    # example
    counter = Kcount(
        name="test", 
        workdir="/tmp/", 
        fastq_path=FILES, 
        kmersize=35, 
        name_split="_R",
    )
    counter.run()

    # statdf is saved to the workdir as a CSV
    print(counter.statsdf.T)
