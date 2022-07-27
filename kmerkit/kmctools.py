#!/usr/bin/env python

"""
Call convenient KMC_tools wrappers
"""

import os
import sys
import shutil
import subprocess
from loguru import logger
from kmerkit.utils import KmerkitError


# FASTPBIN = os.path.join(sys.prefix, "bin", "fastp")
# KMCBIN = os.path.join(sys.prefix, "bin", "kmc")
# KMTBIN = os.path.join(sys.prefix, "bin", "kmc_tools")

# Using shutil alternative installation can be used, not only conda installations.
# FASTPBIN = shutil.which("fastp") 
# KMCBIN = shutil.which("kmc")
# KMTBIN = shutil.which("kmc_tools")


# assert FASTPBIN, (
#     "fastp binary missing; call 'conda install -c bioconda fastp'")

# assert KMCBIN, (
#     "kmc binary missing; call 'conda install kmc -c bioconda'")

# assert KMTBIN, (
#     "kmc_tools binary missing; verify your kmc installation 'conda install kmc -c bioconda'")


# Define info for all dependencies
dependencies = {"fastp": {"var": "FASTPBIN", 
                          "error_message": "- fastp binary missing; call 'conda install -c bioconda fastp'"},
               "kmc": {"var": "KMCBIN", 
                       "error_message": "- kmc binary missing; call 'conda install kmc -c bioconda'"},
               "kmc_tools": {"var": "KMTBIN", 
                             "error_message": "- kmc_tools binary missing; verify your kmc installation 'conda install kmc -c bioconda'"}
               }

# Check binaries and compose error message if needed
error_message = ""
for binary in dependencies:
    binary_path = shutil.which(binary) #check path in any part of the system
    globals()[dependencies[binary]["var"]] = binary_path #set constant variable 
    if not binary_path: #if not found add error message to final display
        error_message += "\n" + dependencies[binary]["error_message"]

# Asset if some binary is missed
assert not error_message, (error_message)
    

    
def info(database, mindepth=0):
    """
    Returns the number of kmers in the database. Option mindepth
    will call 'transform reduce -ci{} beforehand to filter.'
    """
    if mindepth > 1:
        cmd = [
            KMTBIN,
            "transform",
            database,
            "reduce",
            database + "_tmpinfo",
            "-ci{}".format(mindepth)            
        ]
        logger.debug(" ".join(cmd))
        subprocess.run(
            cmd, 
            check=True, 
            stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

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
    if mindepth > 1:
        os.remove(database + ".kmc_suf")
        os.remove(database + ".kmc_pre")        
    return nkmers




def dump(database, write_kmers=True, write_counts=True, min_depth=1, max_depth=100000000):
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
        KMTBIN, "-hp", "-t8",
        "transform", 
        database,
        "-ci{}".format(min_depth),
        "-cx{}".format(max_depth),
        "dump",
        database + "_kmers.txt",
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


if __name__ == "__main__":

    import kmerkit
    kmerkit.set_loglevel("DEBUG")
    from kmerkit.kschema import Project

    PROJ = Project.parse_file("/tmp/test.json")
    prefix = PROJ.dict()['kcount']['data']['hybridus_SLH_AL_1060']['database']

    # print(info(prefix, 1))
    dump(prefix)
