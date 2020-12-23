#!/usr/bin/env python


"""
Load a phenotypes CSV file, or dictionary in the API, 
to assign samples to groups, and then apply kmer_tools operations
like merge or intersect on kmers.
"""


import os
import pandas as pd
from loguru import logger





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
        self.kcpath = os.path.join(self.workdir, self.name + "_kcounts.csv")
        self.phenos = phenos

        # attributes to be filled
        self.statsdf = None
        self.imap = {}

        # fills statsdf and imap
        self.parse_phenos_to_imap()



    def parse_phenos_to_imap(self):
        """
        Parse phenotype input that is either a dict or CSV, 
        check all names for match with database files in workdir.
        """
        
        # load the kcounts database to get all sample names
        self.statdf = pd.read_csv(self.kcpath, index_col=0)

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
    FILES = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R1_.fastq.gz"

    PHENOS = {
        1: ["1A_0", "1B_0", "1C_0", "1D_0"],
        2: ["2E_0", "2F_0", "2G_0", "2H_0"],
    }

    # build database
    counter = kcount.Kcount(
        name="test", 
        workdir="/tmp",
        files=FILES, 
        kmersize=17, 
        name_split="_R",
    )
    counter.run()


    # group database
    grouper = Kgroup(
        name="test",
        workdir="/tmp",
        phenos=PHENOS,
    )
    print(grouper.statdf)

