#!/usr/bin/env python

"""
A Class object for read trimming and subsampling reads with fastp.

tool = ReadTrimming(
	read1='samp_R1.fastq', 
	read2='samp_R2.fastq',
	workdir='/tmp',
	subsample=1000,
)
tool.trim_reads()
print(tool.parse_stats_from_json())
tool.cleanup()
"""

import os
import sys
import json
import subprocess
from loguru import logger
from kmerkit.utils import KmerkitError


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


    def trim_reads(self):
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
            raise KmerkitError(out[0].decode())


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



if __name__ == '__main__':

    # test read trimming 
    FASTQ = "~/Documents/ipyrad/isolation/reftest_fastqs/1A_0_R1_.fastq.gz"
    tool = ReadTrimming(read1=FASTQ, read2=None, workdir="/tmp", subsample=10000)
    tool.trim_reads()

    # check results and cleanup
    print(tool.parse_stats_from_json())
    tool.cleanup()
