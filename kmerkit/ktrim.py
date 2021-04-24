#!/usr/bin/env python

"""
A Class object for read trimming and subsampling reads with fastp.

kmerkit trim -j test.json --subsample 1000 --cores 20 
"""

import os
import json
import subprocess
import concurrent.futures
import numpy as np
from loguru import logger
from kmerkit.kmctools import FASTPBIN
from kmerkit.utils import KmerkitError, num_cpus
from kmerkit.kschema import Project, KtrimData, KtrimBase, KtrimParams


class Ktrim:
    """
    ...
    """
    def __init__(self, json_file, subsample=None):

        # load project 
        self.json_file = json_file
        self.project = Project.parse_file(json_file).dict()

        # load samples from init
        self.fastq_dict = self.project['kinit']['data']

        # store params
        self._params = KtrimParams(subsample=subsample)
        self.params = self._params.dict()


    def check_overwrite(self):
        """
        Warn user of overwriting.
        """
        if self.project['ktrim']:
            logger.error(
                "\nKtrim results exist, use force to overwrite, or consider "
                "using branching to produce new results on a separate named "
                "branch without overwriting previous results."
            )
            raise KmerkitError("Preventing data overwrite")        


    def run(self, force=False, threads=None, workers=None):
        """
        Run fastp on multiple samples in parallel. Each job uses 
        approximately 2 threads, so we submit cores / 2 jobs.
        """
        # check for current step
        if not force:
            self.check_overwrite()

        # store sample results for KtrimData
        data = {}

        # set cores values to limit njobs to ncores / 6
        if workers in [0, None]:
            workers = max(1, int(np.ceil(num_cpus() / 6)))
        threads = (threads if threads else 4)
        logger.debug(f"workers={workers}; threads={threads}")        

        # start jobs and store futures
        future_to_sname = {}
        with concurrent.futures.ProcessPoolExecutor(workers) as lbview:
            for sname in self.fastq_dict:
                read1s, read2s = self.fastq_dict[sname]
                args = (
                    read1s, read2s, 
                    self.project['workdir'], 
                    self.params['subsample'],
                    threads,
                )
                future = lbview.submit(trim_reads, *args)
                future_to_sname[future] = sname

        # track futures and collect results
        for future in concurrent.futures.as_completed(future_to_sname):

            # get sample
            sname = future_to_sname[future]

            # get finished results
            logger.info(f"finished trimming: {sname}")
            result = future.result()
            data[sname] = KtrimData(
                data_in=list(result[0]),
                data_out=list(result[1]),
                fastp_stats=result[2],
            )

        # save to project
        self.project['ktrim'] = KtrimBase(data=data, params=self._params)
        with open(self.json_file, 'w') as out:
            out.write(Project(**self.project).json(indent=4))



def trim_reads(read1, read2, workdir, subsample, threads):
    "function to run on executor"
    tool = ReadTrimming(read1, read2, workdir, subsample, threads)
    tool.run()
    jdata = tool.parse_stats_from_json()
    return (tool.read1, tool.read2), (tool.tmp1, tool.tmp2), jdata



class ReadTrimming:
    """
    Simple read trimming with fastp 
    https://github.com/OpenGene/fastp
    """
    def __init__(
        self, 
        read1=None, 
        read2=None, 
        workdir="/tmp", 
        subsample=None,
        threads=None,
        ):

        # input args
        self.read1 = os.path.realpath(os.path.expanduser(read1))
        self.read2 = (
            os.path.realpath(os.path.expanduser(read2)) if read2
            else None
        )
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.paired = read2 is not None
        self.subsample = subsample
        self.threads = threads

        # output file paths (do not bother gzipping since these are tmp files)
        basename = os.path.basename(read1).rsplit(".fastq", 1)[0]
        self.tmp1 = os.path.join(
            workdir, 
            f'trimmed_{os.path.basename(read1)}'
        ).strip(".gz") + ".gz"
        self.tmp2 = (
            os.path.join(
                workdir, 
                f'trimmed_{os.path.basename(read2)}') if read2
            else None
        ).strip(".gz") + ".gz"

        # paths to stats files
        self.json = os.path.join(workdir, f'{basename}.json')
        self.html = os.path.join(workdir, f'{basename}.html')


    def run(self):
        """
        Calls fastp in subprocess and writes tmpfiles to workdir.
        """
        if self.paired:
            cmd = [
                FASTPBIN,
                "-i", self.read1,
                "-I", self.read2,
                "-o", self.tmp1,
                "-O", self.tmp2,
            ]

        else:
            cmd = [
                FASTPBIN,
                "-i", self.read1,
                "-o", self.tmp1,
            ]

        # force stats files to workdir (temp)
        cmd.extend(["-j", self.json, "-h", self.html])

        # subsampling may be useful for RAD-seq data
        if self.subsample:
            cmd.extend(["--reads_to_process", str(self.subsample)])

        # specify number of threads
        if self.threads:
            cmd.extend(["--thread", str(self.threads)])

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
        return {i: jdata[i] for i in ("summary", "filtering_result")}


    def cleanup(self):
        """
        Called to remove the tmp cleaned files and filtered read statsfiles.
        """
        for tmpfile in [self.json, self.html, self.tmp1, self.tmp2]:
            if tmpfile:
                if os.path.exists(tmpfile):
                    os.remove(tmpfile)



if __name__ == '__main__':

    import kmerkit
    kmerkit.set_loglevel("DEBUG")

    # test read trimming 
    FASTQ = "~/Documents/kmerkit/data/amaranths/tricolor*.gz"
    fastq_dict = kmerkit.utils.get_fastq_dict_from_path(FASTQ)

    proj = kmerkit.init_project(
        name='trimtest', 
        workdir='/tmp', 
        fastq_dict=fastq_dict, 
        force=True,
    )

    ktr = Ktrim('/tmp/trimtest.json')
    ktr.run()
