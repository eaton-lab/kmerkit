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
import json
import subprocess
import concurrent.futures
from loguru import logger
from kmerkit.kmctools import FASTPBIN
from kmerkit.utils import KmerkitError
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


    def run(self, force=False, cores=None):
        """
        Run fastp on multiple samples in parallel. Each job uses 
        approximately 2 threads, so we submit cores / 2 jobs.
        """
        # check for current step
        if not force:
            self.check_overwrite()

        # store sample results for KtrimData
        data = {}

        # set cores values
        if cores == 0:
            cores = None
        if cores is not None:
            cores = int(cores / 2)

        future_to_sname = {}
        with concurrent.futures.ProcessPoolExecutor(cores) as lbview:
            for sname in self.fastq_dict:
                read1s, read2s = self.fastq_dict[sname]
                args = (read1s, read2s, self.project['workdir'], self.params['subsample'])
                future = lbview.submit(trim_reads, *args)
                future_to_sname[future] = sname

        # track progress and collect results
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



def trim_reads(read1, read2, workdir, subsample):
    "function to run on executor"
    tool = ReadTrimming(read1, read2, workdir, subsample)
    tool.run()
    jdata = tool.parse_stats_from_json()
    return (tool.read1, tool.read2), (tool.tmp1, tool.tmp2), jdata



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
        # orig_reads = int(jdata['summary']['before_filtering']['total_reads'])
        # new_reads = int(jdata['summary']['after_filtering']['total_reads'])
        # return orig_reads, new_reads



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
    # print(proj.json(indent=2))

    # tool = ReadTrimming(read1=FASTQ, read2=None, workdir="/tmp", subsample=10000)
    # tool.trim_reads()

    # check results and cleanup
    # print(tool.parse_stats_from_json())
    # tool.cleanup()
