#!/usr/bin/env python

"""
The main class object for executing KMC kmer funcs
KMC GitHub: https://github.com/refresh-bio/KMC
KMC Paper: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399
KMC Suppl: https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/33/17/10.1093_bioinformatics_btx304/3/btx304_supplementary.pdf?Expires=1611337056&Signature=4iJG6giNcZdiDuHwljf-SbELTll74FtIj3YIFvfESeZC~m39EZPJdSXfqAJStvCr5SmH9lHGRCdJGHBLseX~ZunAgFZBFFHikmODBI14Kq84ctkQMihTvBzU1rme~S6MpXcC1Erxavl~ckAEnE7jfwIRJbtm4bSkTk-sEcZKIHqR3H0SwhdN0zMmhqMFkwn~jvNo5Rd~yPFwq8aXtE2CBrMORgVUsu~ACFKnl7sWB2FLtsZ2zp~ENuVz28mYwZWkyFXnUlRq2sKHenjWdw4BChI~QDf5EULM2oXgx4dSDlLTaIaJjsZBGl8tKQGw5Ohz48YEbqO82pn14AyeJkaRsQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA
KMC Docs: https://github.com/refresh-bio/KMC/blob/master/API.pdf
KMC-tools Docs: https://github.com/refresh-bio/KMC/blob/master/kmc_tools.pdf
"""

import os
import subprocess
import pandas as pd



class Kcount:
    """
    Calls the kmc counting functions to create a database.
    """
    def __init__(self, files, outdir, kmersize, name_split="_"):

        # store input params
        self.files = files
        self.outdir = outdir
        self.kmersize = kmersize
        self.name_split = name_split

        # attributes to be filled
        self.names_to_files = {}
        self.files_to_names = {}
        self.statsdf = None

        # constructs statsdf and names_to_files
        self.check_samples()


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

        # fill dataframe with sample names and column names
        self.statsdf = pd.DataFrame(
            columns=sample_names,
            index=[
                'kmers_total', 
                'kmers_unique',
                'kmers_unique_counted',
                'reads_total', 
                'kmers_below_thresh', 
                'kmers_above_thresh',
            ],
            data=0,
            dtype=int,
        )



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
        # create command: 'kmc -k17 inputfastq outname outdir'
        cmd = [
            "kmc", "-k{}".format(self.kmersize), 
            filename,
            self.files_to_names[filename],
            self.outdir,
        ]

        # call subprocess on the command
        out = subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
        )     

        # parse STDOUT to get the kmer stats
        stats = out.stdout.decode().split("Stats:")[-1]

        # store results to the dataframe
        for line in stats.strip().split("\n"):
            key, val = line.split(" : ")
            
            if key.strip() == "Total no. of k-mers":
                self.statsdf.loc["kmers_total"] = int(val.strip())

            if key.strip() == "No. of unique k-mers":
                self.statsdf.loc["kmers_unique"] = int(val.strip())

            if key.strip() == "No. of unique counted k-mers":
                self.statsdf.loc["kmers_unique_counted"] = int(val.strip())

            if key.strip() == "Total no. of reads":
                self.statsdf.loc["reads_total"] = int(val.strip())

            if key.strip() == "No. of unique k-mers above max. threshold":
                self.statsdf.loc["kmers_above_thresh"] = int(val.strip())

            if key.strip() == "No. of unique k-mers below min. threshold":
                self.statsdf.loc["kmers_below_thresh"] = int(val.strip())



if __name__ == "__main__":


    # test dataset
    fileslist = [
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/1A_0_R1_.fastq.gz",
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/2E_0_R1_.fastq.gz",
    ]

    # example
    counter = Kcount(fileslist, "./", 17, name_split="_R")
    counter.run()

    # print stats
    print(counter.statsdf)

    # show output database.
    # ...
