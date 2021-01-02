#!/usr/bin/env python

"""
Call convenient KMC_tools wrappers
"""

import os
import sys
import subprocess
from loguru import logger
from kmpy.utils import KmpyError


KMCBIN = os.path.join(sys.prefix, "bin", "kmc")
KMTBIN = os.path.join(sys.prefix, "bin", "kmc_tools")
assert os.path.exists(KMCBIN)



def info(database, mindepth=1):
    """
    Returns the number of kmers in the database. Option mindepth
    will call 'transform reduce -ci{} beforehand to filter.'
    """
    if mindepth:
        cmd = [
            KMTBIN,
            "transform",
            database,
            "reduce",
            database + "_tmpinfo",
            "-ci{}".format(mindepth)            
        ]
        logger.debug(" ".join(cmd))
        subprocess.run(cmd, check=True)  #, stdout=subprocess.PIPE)

        # rename database 
        database += "_tmpinfo"

    # get info 
    cmd = [
        KMTBIN,
        "info",
        database,
    ]
    out = subprocess.run(
        cmd, 
        check=True, 
        stderr=subprocess.STDOUT, 
        stdout=subprocess.PIPE,
    )

    # parse info from out
    lines = out.stdout.decode().strip().split("\n")
    for line in lines:
        if line.startswith("total"):
            nkmers = int(line.split()[-1])

    # cleanup
    if mindepth:
        os.remove(database)
    return nkmers




def dump(database, write_kmers=True, write_counts=True):
    """
    Dump kmers to text for a specified database. 
    calls: 'kmc_tools transform /tmp/name dump /tmp/name_kmers.txt'

    Parameters
    ----------
    ...
    write_kmers (bool):
        If False then the counts are excluded from dump txt file.
    write_counts (bool):
        If False then the kmers are excluded from dump txt file.        

    Returns
    -------
    None
    """
    # if counts is False then write 'compact' database
    cmd = [
        KMTBIN,
        "transform", 
        database,
        "-ci1",
        "-cx1000000000",
        "dump",
        "-s",
        database + "_kmers.txt",
        "-ci1",
        "-cx1000000000",
        "-cs65535",
    ]

    # call subprocess on the command
    logger.debug(" ".join(cmd))
    subprocess.run(
        cmd, 
        stderr=subprocess.STDOUT, 
        stdout=subprocess.PIPE,
        check=True,
    )

    # only if not both
    if not (write_kmers and write_counts):

        # overwrite with file cutting first column only
        if write_kmers:
            cmd = ["cut", "-f1", database + "_kmers.txt"]
            with open(database + "_kmers.tmp", "w") as ofile: 
                subprocess.run(cmd, stdout=ofile, check=True)
            os.rename(database + "_kmers.tmp", database + "_kmers.txt")

        # overwrite with file cutting first column only
        if write_counts:
            cmd = ["cut", "-f2", database + "_kmers.txt"]
            with open(database + "_kmers.tmp", "w") as ofile: 
                subprocess.run(cmd, stdout=ofile, check=True)
            os.rename(database + "_kmers.tmp", database + "_kmers.txt")
