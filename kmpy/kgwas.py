#!/usr/bin/env python

"""
Genome-wide association study using kmers implemented by taking
kmers converted to genotype calls in Kmatrix and here converting
to a plink formatted input file for gemma analysis.
"""


import os
import sys
import numpy as np




class Kgwas:
	"""

	"""
	def __init__(self, name, workdir, ):

		# store params
		self.name = name
		self.workdir = os.path.realpath(os.path.expanduser(workdir))

        # output prefix for .complex, subtract, countsubtract, intersect
        self.prefix = os.path.join(self.workdir, "kgwas")

        # get gemma binary path
		self.gemma_binary = os.path.join(sys.prefix, "bin", "gemma")
        logger.debug("using gemma binary: {}".format(self.gemma_binary))

        # ...


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
