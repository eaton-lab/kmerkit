#!/usr/bin/env python


"""
Load a phenotypes CSV file, or dictionary in the API, 
to assign samples to groups, and then apply kmer_tools operations
like merge or intersect on kmers.
"""


import os
import pandas as pd





class Kgroup:
    """
    Finds existing database using the CSV in 'workdir/name'.
    Group samples in database based on phenos info (dict or csv)
    and calls union function to get shared kmer counts.
    """
    def __init__(self, name, workdir, phenos):

        # store parameters
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.statsdf = os.path.join(self.workdir, self.name + "_kcounts.csv")
        self.phenos = phenos
        self.imap = {}


    def parse_phenos_to_imap(self):
        """
        Parse phenotype input that is either a dict or CSV
        """
        
        # load the kcounts database to get all sample names
        # ...


    def get_union_complex(self):
        """
        To do a complex operation requires making an input file 
        like the one below.
        
        INPUT:
        name1 = /tmp/name1 -ci0 -cx500 
        name2 = /tmp/name2 -ci0 -cx500
        name3 = /tmp/name3 -ci0 -cx500     

        OUTPUT:
        /tmp/out = name1 + name2 + name3 ...

        OUTPUT_PARAMS:
        -ci0 -cx10 -cs255
        """

        cmd_input = ["INPUT"]


        # TODO: is workdir a global param?
        # cmd: 'kmer_tools [global_params] complex <operations file>'
        cmd = ["kmer_tools", "complex", complex_file]





if __name__ == "__main__":


    import kcount

    # test dataset
    FILES = [
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/1A_0_R1_.fastq.gz",
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/1B_0_R1_.fastq.gz",
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/1C_0_R1_.fastq.gz",
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/2E_0_R1_.fastq.gz",
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/2F_0_R1_.fastq.gz",
        "/home/deren/Documents/ipyrad/isolation/reftest_fastqs/2G_0_R1_.fastq.gz",
    ]

    PHENOS = {
        1: ["1A_0", "1B_0", "1C_0"],
        2: ["2E_0", "2F_0", "2G_0"],
    }

    # build database
    counter = kcount.Kcount(FILES, kmersize=17, name_split="_R", workdir="/tmp")
    counter.run()

    # group database
    grouper = Kgroup(PHENOS, workdir="/tmp")

