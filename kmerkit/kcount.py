#!/usr/bin/env python

"""
Counts kmers in fastq data files using KMC

Links:
KMC GitHub: https://github.com/refresh-bio/KMC
KMC Paper: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399
KMC-tools Docs: https://github.com/refresh-bio/KMC/blob/master/kmc_tools.pdf
KMC Docs: doesn't exist apparently...

Notes: 
    - Canon filter uses count diffs to find presence/absence of forms, 
      this will fail if both forms occur at maxcount (e.g., 255), safer
      to just enforce larger maxcount (65535) by default, and eat the 
      larger disk usage. Test on a superlarge dset before deciding...

    - maxcount could also matter for count=1 with regard to nsamples in 
      kfilter, another reason to just set it higher...
"""

import os
import json
import subprocess
from loguru import logger

from kmerkit.kmctools import KMCBIN
from kmerkit.utils import get_fastq_dict_from_path, KmerkitError
from kmerkit.kschema import KcountBase, KcountData, KcountParams, Project


class Kcount:
    """
    Calls the kmc counting functions to create a database.

    Parameters:
    ===========
    json (str):
        Path to a kmerkit JSON project file.
    kmer_size (int):
        The size of kmers to be counted from reads.
    min_depth (int):
        Minimum depth below which kmers will be excluded. Default=1.
    max_depth (int):
        Maximum depth above which kmers will be excluded. Default=1e9.
    max_count (int):
        Maximum value that will be recorded for a kmer depth. Default=255.
    canonical (bool):
        Count kmers in canonical form (both ways). Default=True.
    """
    def __init__(self, json_file, kmer_size, min_depth, max_depth, max_count, canonical):

        # load project json file to get sample paths
        self.json_file = json_file
        self.project = Project.parse_file(json_file).dict()
        self.fastq_dict = self.project['kinit']['data']
        self.samples = {}

        # use Serializable schema to perform type checking
        self.params = KcountParams(
            kmer_size=kmer_size, 
            min_depth=min_depth,
            max_depth=max_depth,
            max_count=max_count,
            canonical=canonical,
        ).dict()

        # prefix for output files
        self.prefix = os.path.join(
            self.project['workdir'], f"{self.project['name']}_kcount")

        # report which KMC will be used. Same will be used for all.
        logger.info(f"KMC bin: {KMCBIN}")



    def call_kmc_count(self, files, sname, threads=4):
        """
        Calls kmc to count kmers.
          files: original or trimmed fastq files
          sname: sample name
        """
        # to support SE or PE data we list input files in a file
        input_file = f"{self.prefix}_{os.getpid()}.tmp"
        with open(input_file, 'w') as out:
            out.write("\n".join(str(i) for i in files if i))

        # output file path
        output_file = f"{self.prefix}_{sname}"

        # create command: 'kmc -k17 @filelist outname workdir'
        cmd = [
            KMCBIN,
            "-ci{}".format(self.params['min_depth']),
            "-cx{}".format(self.params['max_depth']),
            "-cs{}".format(self.params['max_count']),
            "-k{}".format(self.params['kmer_size']),
            "-j{}".format(output_file + "-stats.json"),
            "-t{}".format(threads),
            "@" + input_file,
            output_file,
            self.project['workdir'],
        ]

        # insert canonical option (TURNS OFF)
        if not self.params['canonical']:
            cmd.insert(4, "-b")

        # call subprocess on the command
        logger.debug(" ".join(cmd))
        out = subprocess.run(
            cmd,
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.project['workdir'],
        )

        # remove input file list
        os.remove(input_file)

        # parse stats from json file
        with open(output_file + "-stats.json", 'r') as indata:
            kdata = json.loads(indata.read())
        os.remove(output_file + "-stats.json")
        return kdata["Stats"]


    def check_overwrite(self):
        """
        Prevent overwriting future files.
        """
        next_steps = [
            self.project.get("kfilter"),
            self.project.get("ktree"),
            self.project.get("kmatrix"),
        ]
        if any(next_steps):
            logger.warning(
                "\nRunning kcount will overwrite previous kmer-counting "
                "results and remove references to any existing downstream "
                "analyses on these files. You must use the force argument "
                "to confirm this action. An alternative recommended workflow "
                "is to use the 'branch' option to create a separate new named "
                "project (new JSON file) from which to start this analysis "
                "without overwriting your previous results"
            )
            raise KmerkitError("Preventing data overwrite")


    def run(self, threads=4, force=False):
        """
        Calls kmc count on all files
        """
        # check for current step
        if not force:
            self.check_overwrite()

        # iterate over samples one at a time.
        for sname in self.fastq_dict:

            # count kmers -------------------------------
            indata = self.fastq_dict[sname]
            kmcstats = self.call_kmc_count(indata, sname, threads)

            # convert kmc stats to KcountData() object to store result
            self.samples[sname] = KcountData(**{
                'reads_total': kmcstats["#Total_reads"],
                'kmers_total': kmcstats["#Total no. of k-mers"],
                'kmers_unique': kmcstats["#Unique_k-mers"],
                'kmers_unique_counted': kmcstats["#Unique_counted_k-mers"],
                'kmers_below_threshold': kmcstats["#k-mers_below_min_threshold"],
                'kmers_above_threshold': kmcstats["#k-mers_above_max_threshold"],
                'database': f"{self.prefix}_{sname}"
            })

        # save to JSON
        self.project['kcount'] = KcountBase(
            params=KcountParams(**self.params),
            data=self.samples,
        )
        project = Project(**self.project)
        with open(self.json_file, 'w') as out:
            out.write(project.json(indent=4))






if __name__ == "__main__":


    import kmerkit
    kmerkit.set_loglevel("DEBUG")

    # FILES = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R*_.fastq.gz"
    FILES = "~/Documents/kmerkit/data/amaranths/hybridus_*.fastq.gz"
    FASTQ_DICT = get_fastq_dict_from_path(FILES, "_R")

    # init a project
    kmerkit.init_project('test', '/tmp', FASTQ_DICT, force=True)

    # load project json and count kmers
    kco = Kcount(
        "/tmp/test.json",        
        kmer_size=35,
        min_depth=1,
        max_depth=1e9,
        max_count=255,
        canonical=True,
    )
    kco.run(threads=4)

    # statdf is saved to the workdir as a CSV
    # print(counter.statsdf.T)
