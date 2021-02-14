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

TODO:
    - t x = option to limit threads
"""

import os
import json
import subprocess
import pandas as pd
from loguru import logger

from kmerkit.kmctools import KMCBIN
from kmerkit.read_trimming import ReadTrimming
from kmerkit.utils import get_fastq_dict_from_path


class Kcount:
    """
    Calls the kmc counting functions to create a database.

    Parameters:
    ===========
    name (str):
        A name prefix to be used for output files.
    workdir (str):
        A directory for storing output and temporary files.
    kmersize (int):
        The size of kmers to be counted from reads.
    fastq_dict: dict:
        A dictionary mapping sample names to a list of fastq files.
    trim_reads (bool):
        Reads are trimmed for adapters and low quality bases by 'fastp'
        prior to counting kmers.
    subsample_reads (int):
        Subsample reads from each sample to include at most N reads
        (or PE readpairs). This can be useful for samples with variable 
        coverage.
    mindepth (int):
        Minimum depth below which kmers will be excluded. Default=1.
    maxdepth (int):
        Maximum depth above which kmers will be excluded. Default=1e9.
    maxcount (int):
        Maximum value that will be recorded for a kmer depth. Default=255.
    canonical (bool):
        Count kmers in canonical form (both ways). Default=True.
    """
    def __init__(
        self, 
        name, 
        workdir, 
        kmersize, 
        fastq_dict,
        trim_reads=False, 
        subsample_reads=None,
        mindepth=1,
        maxdepth=1e9,
        maxcount=65530,
        canonical=True,
        ):

        # store input params
        self.name = name
        self.fastq_dict = fastq_dict
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.kmersize = kmersize
        self.trim_reads = trim_reads
        self.subsample_reads = (int(subsample_reads) if subsample_reads else 0)
        self.mindepth = int(mindepth)
        self.maxdepth = int(maxdepth)
        self.maxcount = int(maxcount)
        self.canonical = canonical

        # output file prefix
        self.prefix = os.path.join(self.workdir, f"kcount_{self.name}")

        # report which KMC will be used. Same will be used for all.
        logger.info(f"KMC bin: {KMCBIN}")

        # constructs statsdf and names_to_files
        self.check_fastq_dict()

        # dataframe with sample names and column names
        self.statsdf = pd.DataFrame(
            index=sorted(self.fastq_dict),
            columns=[
                'reads_total', 
                'reads_passed_trimming',
                'kmers_total', 
                'kmers_unique',
                'kmers_unique_counted',
                'kmers_below_thresh', 
                'kmers_above_thresh',
                'database',
                'fastqs'
            ],
            data=0,
            dtype=int,
        )

        # add additional columns of uniform values w/ param settings
        self.statsdf['kmersize'] = int(self.kmersize)
        self.statsdf['trimmed'] = int(self.trim_reads)
        self.statsdf['subsample_reads'] = self.subsample_reads
        self.statsdf['mindepth'] = self.mindepth
        self.statsdf['maxdepth'] = self.maxdepth
        self.statsdf['maxcount'] = self.maxcount
        self.statsdf['canonical'] = int(self.canonical)



    def check_fastq_dict(self):
        """
        Expand file paths and check that they exist. Also,
        TODO: check sample names and do not allow strange characters.
        """
        for oldkey, val in self.fastq_dict.items():

            # get new key without any strange characters
            newkey = (oldkey
                .replace("-", "_")
                .replace("@", "_")
                .replace(" ", "_")
            )

            # check type of val
            assert isinstance(val, list), (
                "Filepaths in fastq_dict must be stored as list objects.\n"
                "Example: fdict = {'a': ['a_R1.fastq', 'a_R2.fastq'], 'b'...}"
            )

            # expand each path in val and check exists
            file_list = []
            for path in val:
                fullpath = os.path.realpath(os.path.expanduser(path))
                assert os.path.exists(fullpath), (
                    f"file {fullpath} in fastq_dict cannot be found.")
                file_list.append(fullpath)

            # remove old key (in case of strange characters) and store new.
            del self.fastq_dict[oldkey]
            self.fastq_dict[newkey] = sorted(file_list)



    def call_kmc_count(self, files, sname):
        """
        Calls kmc to count kmers.
          files: original or trimmed fastq files
          sname: sample name
        """
        # to support SE or PE data we list input files in a file
        input_file = self.prefix + ".tmp"
        with open(input_file, 'w') as out:
            out.write("\n".join(i for i in files if i))

        # output file path
        output_file = f"{self.prefix}_{sname}"

        # create command: 'kmc -k17 @filelist outname workdir'
        cmd = [
            KMCBIN,
            "-ci{}".format(self.mindepth),
            "-cx{}".format(self.maxdepth),
            "-cs{}".format(self.maxcount),
            "-k{}".format(self.kmersize),
            "-j{}".format(output_file + "-stats.json"),
            "@" + input_file,
            output_file,
            self.workdir,
        ]

        # insert canonical option (TURNS OFF)
        if not self.canonical:
            cmd.insert(4, "-b")

        # call subprocess on the command
        logger.debug(" ".join(cmd))
        out = subprocess.run(
            cmd,
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )

        # remove input file list
        os.remove(input_file)

        # parse stats from json file
        with open(output_file + "-stats.json", 'r') as indata:
            kdata = json.loads(indata.read())
        return kdata["Stats"]



    def kmer_stats(self, ofiles, sname, kmerstats, trimstats):
        """
        Store stats to dataframe. KMC has irregular use of -, _, ' ' 
        separators in keys here, so there is worry that this could
        break if they fix this in the future, so I replace all with "-"
        """
        # replace irregular separators and set values to ints
        kdict = {
            key.replace("_", "-").replace(" ", "-").strip("#"): int(val)
            for key, val in kmerstats.items()
        }

        # enter to dataframe        
        self.statsdf.loc[sname, "kmers_total"] = kdict["Total-no.-of-k-mers"]
        self.statsdf.loc[sname, "kmers_unique"] = kdict["Unique-k-mers"]
        self.statsdf.loc[sname, "kmers_unique_counted"] = kdict["Unique-counted-k-mers"]
        self.statsdf.loc[sname, "reads_total"] = kdict["Total-reads"]
        self.statsdf.loc[sname, "kmers_above_thresh"] = kdict["k-mers-above-max-threshold"]
        self.statsdf.loc[sname, "kmers_below_thresh"] = kdict["k-mers-below-min-threshold"]

        # save kmc db prefix path and original fastq names
        self.statsdf.loc[sname, "database"] = f"{self.prefix}_{sname}"
        self.statsdf.loc[sname, "fastqs"] = ",".join(ofiles)
        # os.path.basename(i) for i in ofiles if i]

        # store trimming info
        if not trimstats:
            self.statsdf.loc[sname, "reads_passed_trimming"] = self.statsdf.loc[sname, "reads_total"]
        else:
            self.statsdf.loc[sname, "reads_total"] = trimstats[0]
            self.statsdf.loc[sname, "reads_passed_trimming"] = trimstats[1]

        # correct read counts for PE double-counting
        if len(self.fastq_dict[sname]) > 1:
            self.statsdf.loc[sname, "reads_total"] = int(
                self.statsdf.loc[sname, "reads_total"] / 2)
            self.statsdf.loc[sname, "reads_passed_trimming"] = int(
                self.statsdf.loc[sname, "reads_passed_trimming"] / 2)

        # save database csv to the workdir
        path = os.path.join(self.prefix + ".csv")
        self.statsdf.to_csv(path)
        logger.info(f"new database: {self.name}_{sname}")

        # remove stats file
        os.remove(f"{self.prefix}_{sname}-stats.json")



    def run(self):
        """
        Calls kmc count on all files
        """
        for sname in self.fastq_dict:

            # get files
            files = ofiles = self.fastq_dict[sname]
            read1 = files[0]
            read2 = (None if len(files) == 1 else files[1])

            # trim reads --------------------------------      
            trimstats = None
            if self.trim_reads:
                tool = ReadTrimming(
                    read1=read1, 
                    read2=read2,
                    workdir=self.workdir,
                    subsample=self.subsample_reads,
                )
                tool.trim_reads()
                trimstats = tool.parse_stats_from_json()
                files = [tool.tmp1, tool.tmp2]

            # count kmers -------------------------------
            kmerstats = self.call_kmc_count(files, sname)

            # store results to dataframe
            self.kmer_stats(ofiles, sname, kmerstats, trimstats)

            # cleanup tmp files from read trimming
            if self.trim_reads:
                tool.cleanup()






if __name__ == "__main__":

    # test dataset w/ R1 and R2 files
    # FILES = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R*_.fastq.gz"
    FILES = "~/Documents/kmerkit/data/amaranths/hybridus_*.fastq.gz"
    FASTQ_DICT = get_fastq_dict_from_path(FILES, "_R")

    import kmerkit
    kmerkit.set_loglevel("INFO")

    # example
    counter = Kcount(
        name="hybridus", 
        workdir="/tmp/", 
        fastq_dict=FASTQ_DICT,
        kmersize=31,
        trim_reads=True,
        subsample_reads=5e5,
        mindepth=1,
        maxdepth=1e9,
        maxcount=255,
        canonical=True,
    )
    # print(counter.statsdf.T)
    counter.run()

    # statdf is saved to the workdir as a CSV
    # print(counter.statsdf.T)
