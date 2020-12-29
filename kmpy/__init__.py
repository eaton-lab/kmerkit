#!/usr/bin/env python

"""
Primary classes of Kmpy (Kmer-python)
"""

__version__ = "0.0.3"


from kmpy.kcount import Kcount
from kmpy.kgroup import Kgroup
from kmpy.kfilter import Kfilter
from kmpy.kmatrix import Kmatrix
from kmpy.utils import set_logger

# start the logger in INFO
set_logger("DEBUG")

