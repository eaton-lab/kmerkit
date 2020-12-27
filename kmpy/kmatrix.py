#!/usr/bin/env python

"""
Takes kmer databases for N samples and constructs a matrix
of (nsamples, nkmers) of either count or binary data.
"""

import os
import numpy as np


class Kmatrix:
    """
    Constructs a matrix of (nsamples x nkmers) as an np.array.

    Finds databases for counted samples from the path to the 
    counts CSV file. In many cases the matrix is unlikely to fit into 
    memory, so we will use memory-mapping to write to the file 
    efficiently.

    Memory-mapped files are used for accessing small segments of large
    files on disk, without reading the entire file into memory. NumPy’s 
    memmap’s are array-like objects. 
    """
    def __init__(self, name, workdir):

        # store params
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.prefix = os.path.join(self.workdir, self.name + "_matrix")

        # can the matrix fit into memory?
        self.matrix = np.memmap(filename, dtype=np.bool_, mode="w+", shape=(3, 4))

