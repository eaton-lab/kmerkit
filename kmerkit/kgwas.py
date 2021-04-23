#!/usr/bin/env python

"""
Kcounts -> [Kfilter] -> Kmatrix -> Kgwas

Genome-wide association study using kmers implemented by taking
kmers converted to genotype calls in Kmatrix and here converting
to a plink formatted input file for gemma analysis.
"""


import os
import sys
import numpy as np
from loguru import logger



class Kgwas:
    """

    """
    def __init__(self, name, workdir, phenos, trait):

        # store params
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.phenos = phenos
        self.trait = trait

        # output prefix for .complex, subtract, countsubtract, intersect
        self.prefix = os.path.join(self.workdir, "kgwas")

        # get gemma binary path
        self.gemma_binary = os.path.join(sys.prefix, "bin", "gemma")
        self.plink_binary = os.path.join(sys.prefix, "bin", "plink")
        logger.debug("gemma bin: {}".format(self.gemma_binary))
        logger.debug("plink bin: {}".format(self.gemma_binary))        

        # attributes to be filled
        self.nsamples
        self.samples

        # load datafiles
        self.load_geno_matrix()
        self.load_phenos()

        # the first six useless columns of plink
        sadcols = pd.DataFrame({
            "fam": 0,
            "ind": range(self.nsamples),
            "father": 0,
            "mother": 0,
            "sex": 0,      # although...
            "phenotype": self.phenos.loc[self.samples, self.trait],
        })



    def load_phenos(self):
        """

        Check that samples in phenos are also in ...
        """



    def load_geno_matrix(self):
        """
        """



    def write_genos_to_plink(self):
        """
        
        """



    def get_kinship_matrix(self):
        """

        """


    def run(self):
        """

        """



if __name__ == "__main__":

    kgw = Kgwas(
        name="hyb",
        workdir="/tmp",

    )

    kgw#.run()
