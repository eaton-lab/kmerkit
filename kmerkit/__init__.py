#!/usr/bin/env python

"""
Primary classes of kmerkit

TODO: A tool to grab the largest file and kcount as cs255 and cs655
to test for whether we should use uint8 versus uint16. 
Why not just handle this completely for the user to avoid problems.
"""

__version__ = "0.0.16"

from loguru import logger

from kmerkit.kinit import init_project
from kmerkit.kcount import Kcount
from kmerkit.kextract import Kextract
from kmerkit.kmatrix import Kmatrix
from kmerkit.kmctools import KMCBIN, KMTBIN
from kmerkit.utils import set_loglevel

# start the logger in INFO
set_loglevel("WARNING")
