#!/usr/bin/env python

"""
Primary classes of Kmpy (Kmer-python)

TODO: A tool to grab the largest file and kcount as cs255 and cs655
to test for whether we should use uint8 versus uint16. 
Why not just handle this completely for the user to avoid problems.
"""

__version__ = "0.0.5"

from loguru import logger
from kmpy.kcount import Kcount
from kmpy.kgroup import Kgroup
from kmpy.kextract import Kextract
from kmpy.kmatrix import Kmatrix
from kmpy.kmctools import KMCBIN, KMTBIN
from kmpy.utils import set_loglevel

# start the logger in INFO
set_loglevel("DEBUG")
logger.info(f"KMC bin: {KMCBIN}")
