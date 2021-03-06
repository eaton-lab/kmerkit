#!/usr/bin/env python

"""
General utilities
"""

import os
import sys
import glob
from copy import copy
from loguru import logger


COMPLEX = """\n
INPUT:
{input_string}

OUTPUT:
{output_string}

OUTPUT_PARAMS:
{output_params}
"""



class KmerkitError(Exception):
    """
    Exception handler that does clean exit for CLI, but also prints
    the traceback and cleaner message for API.
    """
    def __init__(self, *args, **kwargs):
        # raise the exception with this string message and a traceback
        Exception.__init__(self, *args, **kwargs)



class Group:
    """
    Stores subset of samples with a given trait value, and which are 
    present in both databases. Can return sample list or a string for
    computing the union of kmers in kmer_tools format: (A + B + C + D)   

    cmode: see kmer_tools docs. 

    HACK:
    double adds an extra copy of the first sample. This is useful for 
    performing count_subtraction on:
        var = (A + B + C) - (A * B * C)`
        countA = var ~ A
    In the above case if the kmer is only in A it will be removed. If instead
    we do:
        var = (nullA + nullB + nullC + A + B + C) - (nullA * nullB * nullC * A * B * C)`
        countA = var ~ A
    In this case the total counts are inflated by the counts in A, but all
    kmers in var will always be present in each countA, countB, etc. file.
    """
    def __init__(self, samples, double=False, cmode=None):
        # store inputs
        self.samples = copy(samples)
        self.cmode = (cmode if cmode else "")

        # optional, add double counter: see doc string
        if double:
            self.samples += [f"null_{i}" for i in self.samples]

        # to be filled
        self.ustring = ""
        self.istring = ""        

        # filler funcs
        self.get_union_string()
        self.get_intersect_string()        


    def get_union_string(self):
        """
        Builds the union string of kmers in any samples
        The counter shows the sum count of the kmer across all samples.
        """
        self.ustring = f" + {self.cmode}".join(self.samples)


    def get_intersect_string(self):
        """
        Builds the intersect string of kmers in all samples. 
        The counter shows the min count of the kmer is any one sample.
        """
        self.istring = f" * {self.cmode}".join(self.samples)


    def get_string(self, arg):
        """
        Returns either the union or intersect string depending on arg
        """
        if arg == "union":
            return self.ustring
        return self.istring




def get_fastq_dict_from_path(fastq_path, name_split="_R"):
    """
    Select multiple files using a wildcard selector
    (e.g., "./data/*.fastq.gz") and parse names from files by 
    splitting on a string separator (such as "_R" to split before
    the _R1 or _R2 designator) and selecting name before the split.

    Parameters
    ----------
    fastq_path (str):
        A wildcard string to select multiple fastq files. Examples: 
        "./files/*.fastq" or "./data/samples-[0-9]*.fastq.gz".
    name_split (str):
        Split names on this character to extract sample names from 
        fastq file names. For example, "_R" is frequently used to 
        split names prior to the "_R1" or "_R2" read specifier.
    """
    # the dictionary to be filled
    fastq_dict = {}

    # expand fastq_path
    fastq_path = os.path.realpath(os.path.expanduser(fastq_path))

    # raise an exception if no files were found
    files = glob.glob(fastq_path)
    if not any(files):
        msg = f"no fastq files found at: {fastq_path}"
        logger.error(msg)
        raise KmerkitError(msg)

    # sort the input files
    files = sorted(files)

    # report on found files
    logger.info("found {} input files".format(len(files)))

    # split file names to keep what comes before 'name_split'
    sample_names = [
        os.path.basename(i.split(name_split)[0]) for i in files
    ]

    # do not allow .fastq to still be present in names
    if any(['.fastq' in i for i in sample_names]):
        raise KmerkitError(
            "Failed extracting sample names from fastq filenames. "
            "Try modifying the name_split argument."
        )

    # check that all sample_names are unique
    if len(set(sample_names)) != len(sample_names):
            
        # if not, then check each occurs 2X (PE reads)
        if not all([sample_names.count(i) == 2 for i in sample_names]):
            raise KmerkitError(
                "Sample names are not unique, or in sets of 2 (PE). "
                "You may need to try a different name_split setting."
            )                
        logger.info("detected PE data")

    # store dict mapping names and files (or file pairs for PE)
    for sname, file in zip(sample_names, files):

        # names to input fastqs
        if sname in fastq_dict:
            fastq_dict[sname].append(file)
        else:
            fastq_dict[sname] = [file]

    # pretty print to logger debug
    pretty = "FASTQ_DICT:\n"
    for key in sorted(fastq_dict):
        pretty += f"{key}\n"
        for val in fastq_dict[key]:
            pretty += f"    {val}\n"
    logger.debug(pretty.strip())

    return fastq_dict



def set_loglevel(loglevel="DEBUG"):#, logfile=None):
    """
    Config and start the logger
    """
    config = {
        "handlers": [
            {
                "sink": sys.stderr, 
                "format": (
                    "{time:hh:mm} <level>{level: <7}</level> <white>|</white> "
                    "<magenta>{file: <11}</magenta> <white>|</white> "
                    "<cyan>{function: <16}</cyan> <white>|</white> "
                    "<level>{message}</level>"
                ),
                "level": loglevel,
                },
            # {
                # "sink": logfile,                   
                # "format": "{time:YYYY-MM-DD} | {function} | {message}",
                # "level": "INFO",
                # }
        ]
    }
    logger.configure(**config)
    logger.enable("kmerkit")
