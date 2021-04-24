#!/usr/bin/env python

"""
Counts kmers in fastq data files using KMC

Links:
KMC GitHub: https://github.com/refresh-bio/KMC
KMC Paper: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399
KMC-tools Docs: https://github.com/refresh-bio/KMC/blob/master/kmc_tools.pdf
KMC Docs: doesn't exist apparently...

KMC params to explore for optimization:
  -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12
  -sm - use strict memory mode (memory limit from -m<n> switch will not be exceeded)
  -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9
  -f<a/q/m/bam> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam); default: FASTQ
  -ci<value> - exclude k-mers occurring less than <value> times (default: 2)
  -cs<value> - maximal value of a counter (default: 255)
  -cx<value> - exclude k-mers occurring more of than <value> times (default: 1e9)
  -b - turn off transformation of k-mers into canonical form
  -r - turn on RAM-only mode 
  -n<value> - number of bins 
  -t<value> - total number of threads (default: no. of CPU cores)
  -sf<value> - number of FASTQ reading threads
  -sp<value> - number of splitting threads
  -sr<value> - number of threads for 2nd stage

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
import tempfile
import concurrent.futures
import numpy as np
from loguru import logger

from kmerkit.kmctools import KMCBIN
from kmerkit.utils import get_fastq_dict_from_path, KmerkitError, num_cpus
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

        # select ktrim data if exists else get kinit data
        if self.project.get('ktrim'):
            self.fastq_dict = {
                sname: self.project['ktrim']['data'][sname]['data_out']
                for sname in self.project['ktrim']['data']
            }
        else:
            self.fastq_dict = self.project['kinit']['data']

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
        logger.info(f"Counting kmers for {len(self.fastq_dict)} samples with KMC")



    def call_kmc_count(self, files, sname, threads, max_ram):
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

        # run from a tempfile
        with tempfile.TemporaryDirectory(
            suffix=str(os.getpid()), 
            dir=self.project['workdir'],
            ) as tmpdir:

            # create command: 'kmc -k17 @filelist outname workdir'
            cmd = [
                KMCBIN,
                "-ci{}".format(self.params['min_depth']),
                "-cx{}".format(self.params['max_depth']),
                "-cs{}".format(self.params['max_count']),
                "-k{}".format(self.params['kmer_size']),
                "-j{}".format(output_file + "-stats.json"),
                "-t{}".format(threads),
                "-m{}".format(max_ram),
                "@" + input_file,
                output_file,
                tmpdir,
            ]

            # insert canonical option (TURNS OFF)
            if not self.params['canonical']:
                cmd.insert(4, "-b")

            # log command
            logger.debug(" ".join(cmd))

            # call subprocess from a tempdir to avoid conflict in kmc tmps
            out = subprocess.run(
                cmd,
                stderr=subprocess.STDOUT, 
                stdout=subprocess.PIPE,
                check=True,
                cwd=tmpdir,
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
        if self.project['kcount']:
            logger.error(
                "\nRunning kcount will overwrite previous kmer-counting "
                "results and remove references to any existing downstream "
                "analyses on these files. You must use the force argument "
                "to confirm this action. An alternative recommended workflow "
                "is to use the 'branch' option to create a separate new named "
                "project (new JSON file) from which to start this analysis "
                "without overwriting your previous results."
            )
            raise KmerkitError("Preventing data overwrite")


    def run(self, workers=None, threads=4, max_ram=12, force=False):
        """
        Calls kmc count on all files
        """
        # check for current step
        if not force:
            self.check_overwrite()

        # set cores values to limit njobs to ncores / 4
        if workers in [0, None]:
            workers = max(1, int(np.ceil(num_cpus() / 4)))
            threads = 4
        # if user set workers, then scale threads to match
        else:
            workers = int(workers)
            threads = (threads if threads else int(num_cpus() / workers))
            # (threads if threads else 4)
        logger.debug(f"workers={workers}; threads={threads}; max_ram={max_ram}")        

        # start jobs and store futures
        future_to_sname = {}
        with concurrent.futures.ProcessPoolExecutor(workers) as lbview:
            for sname in self.fastq_dict:
                args = (self.fastq_dict[sname], sname, threads, max_ram)
                future = lbview.submit(self.call_kmc_count, *args)
                future_to_sname[future] = sname

        # track futures and collect results
        data = {}
        for future in concurrent.futures.as_completed(future_to_sname):

            # get sample
            sname = future_to_sname[future]

            # get finished results
            kmcstats = future.result()

            # log completion
            logger.info(
                "new database: '{}' reads={} uniq-kmers={}"
                .format(
                    sname, kmcstats["#Total_reads"], kmcstats["#Unique_k-mers"]
            ))

            # convert kmc stats to KcountData() object to store result
            data[sname] = KcountData(**{
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
            data=data,
        )
        project = Project(**self.project)
        with open(self.json_file, 'w') as out:
            out.write(project.json(indent=4))


    # def oldrun(self):
    #     # iterate over samples one at a time.
    #     samples = {}
    #     for sname in self.fastq_dict:

    #         # count kmers -------------------------------
    #         indata = self.fastq_dict[sname]
    #         kmcstats = self.call_kmc_count(indata, sname, threads)

    #         # convert kmc stats to KcountData() object to store result
    #         samples[sname] = KcountData(**{
    #             'reads_total': kmcstats["#Total_reads"],
    #             'kmers_total': kmcstats["#Total no. of k-mers"],
    #             'kmers_unique': kmcstats["#Unique_k-mers"],
    #             'kmers_unique_counted': kmcstats["#Unique_counted_k-mers"],
    #             'kmers_below_threshold': kmcstats["#k-mers_below_min_threshold"],
    #             'kmers_above_threshold': kmcstats["#k-mers_above_max_threshold"],
    #             'database': f"{self.prefix}_{sname}"
    #         })
    #         logger.info("new database: '{}' reads={} uniq-kmers={}"
    #             .format(
    #                 sname, kmcstats["#Total_reads"], kmcstats["#Unique_k-mers"]
    #             )
    #         )

    #     # save to JSON
    #     self.project['kcount'] = KcountBase(
    #         params=KcountParams(**self.params),
    #         data=samples,
    #     )
    #     project = Project(**self.project)
    #     with open(self.json_file, 'w') as out:
    #         out.write(project.json(indent=4))






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
        max_count=65535,
        canonical=True,
    )
    kco.run(threads=4)

    # statdf is saved to the workdir as a CSV
    # print(counter.statsdf.T)
