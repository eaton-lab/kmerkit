#!/usr/bin/env python

"""
General utilities
"""

import os
import sys
import json
import subprocess
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



class KmpyError(Exception):
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





class ReadTrimming:
    """
    Simple read trimming with fastp 
    https://github.com/OpenGene/fastp
    """
    def __init__(self, read1=None, read2=None, workdir="/tmp", subsample=None):

        # input args
        self.read1 = os.path.realpath(os.path.expanduser(read1))
        self.read2 = (
            os.path.realpath(os.path.expanduser(read2)) if read2
            else None
        )
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.paired = read2 is not None
        self.subsample = subsample
        self.fastp_binary = os.path.join(sys.prefix, "bin", "fastp")

        # output file paths (do not bother gzipping since these are tmp files)
        basename = os.path.basename(read1).rsplit(".fastq", 1)[0]
        self.tmp1 = os.path.join(
            workdir, 
            f'trimmed_{os.path.basename(read1)}'
        ).strip(".gz")
        self.tmp2 = (
            os.path.join(
                workdir, 
                f'trimmed_{os.path.basename(read2)}') if read2
            else None
        ).strip(".gz")

        # paths to stats files
        self.json = os.path.join(workdir, f'{basename}.json')
        self.html = os.path.join(workdir, f'{basename}.html')


    def run(self):
        """
        Calls fastp in subprocess and writes tmpfiles to workdir.
        """
        if self.paired:
            cmd = [
                self.fastp_binary, 
                "-i", self.read1,
                "-I", self.read2,
                "-o", self.tmp1,
                "-O", self.tmp2,
            ]

        else:
            cmd = [
                self.fastp_binary,
                "-i", self.read1,
                "-o", self.tmp1,
            ]

        # force stats files to workdir (temp)
        cmd.extend(["-j", self.json, "-h", self.html])

        # subsampling may be useful for RAD-seq data
        if self.subsample:
            cmd.extend(["--reads_to_process", str(self.subsample)])

        # logger record
        logger.debug(" ".join(cmd))

        # run the command
        proc = subprocess.Popen(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            cwd=self.workdir,
        )
        out = proc.communicate()
        if proc.returncode:
            logger.error("FASTP ERROR")
            raise KmpyError(out[0].decode())


    def parse_stats_from_json(self):
        """
        Get stats from JSON file.
        """
        with open(self.json, 'r') as indata:
            jdata = json.loads(indata.read())

        orig_reads = int(jdata['summary']['before_filtering']['total_reads'])
        new_reads = int(jdata['summary']['after_filtering']['total_reads'])
        return orig_reads, new_reads



    def cleanup(self):
        """
        Called to remove the tmp cleaned files and filtered read statsfiles.
        """
        for tmpfile in [self.json, self.html, self.tmp1, self.tmp2]:
            if tmpfile:
                if os.path.exists(tmpfile):
                    os.remove(tmpfile)






def set_logger(loglevel="DEBUG"):#, logfile=None):
    """
    Config and start the logger
    """
    config = {
        "handlers": [
            {
                "sink": sys.stdout, 
                "format": (
                    "{time:hh:mm} <white>|</white> "
                    "<magenta>{file: <10}</magenta> <white>|</white> "
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
    logger.enable("kmpy")



if __name__ == "__main__":

    # test read trimming 
    FASTQ = "~/Documents/ipyrad/isolation/reftest_fastqs/1A_0_R1_.fastq.gz"
    tool = ReadTrimming(read1=FASTQ, read2=None, workdir="/tmp", subsample=10000)
    tool.run()

    # check results and cleanup
    print(tool.parse_stats_from_json())
    tool.cleanup()
