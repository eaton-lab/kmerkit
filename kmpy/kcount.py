#!/usr/bin/env python

"""
The main class object for executing KMC kmer funcs
KMC GitHub: https://github.com/refresh-bio/KMC
KMC Paper: https://academic.oup.com/bioinformatics/article/33/17/2759/3796399
KMC Suppl: https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/33/17/10.1093_bioinformatics_btx304/3/btx304_supplementary.pdf?Expires=1611337056&Signature=4iJG6giNcZdiDuHwljf-SbELTll74FtIj3YIFvfESeZC~m39EZPJdSXfqAJStvCr5SmH9lHGRCdJGHBLseX~ZunAgFZBFFHikmODBI14Kq84ctkQMihTvBzU1rme~S6MpXcC1Erxavl~ckAEnE7jfwIRJbtm4bSkTk-sEcZKIHqR3H0SwhdN0zMmhqMFkwn~jvNo5Rd~yPFwq8aXtE2CBrMORgVUsu~ACFKnl7sWB2FLtsZ2zp~ENuVz28mYwZWkyFXnUlRq2sKHenjWdw4BChI~QDf5EULM2oXgx4dSDlLTaIaJjsZBGl8tKQGw5Ohz48YEbqO82pn14AyeJkaRsQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA
KMC Docs: https://github.com/refresh-bio/KMC/blob/master/API.pdf
KMC-tools Docs: https://github.com/refresh-bio/KMC/blob/master/kmc_tools.pdf


TODO: 
    - KMC cannot handle dashes in names! need to fix names while keeping 
      paths to original data files correct.
    -hp = hide progress
    -t x = option to limit threads
    -cs  maximal value of a counter... //???

"""

import os
import sys
import glob
import subprocess
import pandas as pd
from loguru import logger
from kmpy.utils import KmpyError, ReadTrimming




# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


class Kcount:
    """
    Calls the kmc counting functions to create a database.

    Parameters:
    ===========
    name (str):
        A name prefix to be used for output files.
    workdir (str):
        A directory for storing output and temporary files.
    fastq_path (str):
        A wildcard string to select multiple fastq files. Examples: 
        "./files/*.fastq" or "./data/samples-[0-9]*.fastq.gz".
    kmersize (int):
        The size of kmers to be counted from reads.
    name_split (str):
        Split names on this character to extract sample names from 
        fastq file names. For example, "_R" is frequently used to 
        split names prior to the "_R1" or "_R2" read specifier.
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
    """
    def __init__(
        self, 
        name, 
        workdir, 
        fastq_path, 
        kmersize, 
        name_split="_", 
        trim_reads=False, 
        subsample_reads=None,
        mindepth=1,
        maxdepth=1e9,
        maxcount=255,
        ):

        # store input params
        self.name = name
        self.fastq_path = os.path.realpath(os.path.expanduser(fastq_path))
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.kmersize = kmersize
        self.name_split = name_split
        self.trim_reads = trim_reads
        self.subsample_reads = (int(subsample_reads) if subsample_reads else None)
        self.prefix = os.path.join(self.workdir, f"kcount_{self.name}")
        self.mindepth = int(mindepth)
        self.maxdepth = int(maxdepth)
        self.maxcount = int(maxcount)

        # attributes to be filled
        self.files = []
        self.names_to_infiles = {}
        self.infiles_to_names = {}
        self.names_to_outfiles = {}
        self.statsdf = None

        # get kmc binary from conda 
        self.kmc_binary = os.path.join(sys.prefix, "bin", "kmc")
        logger.info(f"bin: {self.kmc_binary}")

        # expands wildcard operators to get filenames from string
        self.expand_filenames()

        # constructs statsdf and names_to_files
        self.check_samples()

        # dataframe with sample names and column names
        self.statsdf = pd.DataFrame(
            index=sorted(self.names_to_infiles),
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



    def expand_filenames(self):
        """
        Allows for selecting multiple input files using wildcard
        operators like "./fastqs/*.fastq.gz" to select all fastq.gz
        files in the folder fastqs.
        """
        if isinstance(self.fastq_path, (str, bytes)):
            self.files = glob.glob(self.fastq_path)

            # raise an exception if no files were found
            if not any(self.files):
                msg = f"no fastq files found at: {self.fastq_path}"
                logger.error(msg)
                raise KmpyError(msg)

            # sort the input files
            self.files = sorted(self.files)

        # report on found files
        logger.info("found {} input files".format(len(self.files)))



    def check_samples(self):
        """
        Gets sample names from the input files and checks that all 
        have the same style of suffix (e.g., .fastq.gz).
        """
        # split file names to keep what comes before 'name_split'
        sample_names = [
            os.path.basename(i.split(self.name_split)[0]) for i in self.files
        ]

        # do not allow .fastq to still be present in names
        if any(['.fastq' in i for i in sample_names]):
            raise KmpyError(
                "Failed extracting sample names from fastq filenames. "
                "Try modifying the name_split argument."
            )

        # check that all sample_names are unique
        if len(set(sample_names)) != len(sample_names):
            
            # if not, then check each occurs 2X (PE reads)
            if not all([sample_names.count(i) == 2 for i in sample_names]):
                raise KmpyError(
                    "Sample names are not unique, or in sets of 2 (PE). "
                    "You may need to try a different name_split setting."
                )                
            logger.info("detected PE data")

        # store dict mapping names and files (or file pairs for PE)
        for sname, file in zip(sample_names, self.files):

            # names to input fastqs
            if sname in self.names_to_infiles:
                self.names_to_infiles[sname].append(file)
            else:
                self.names_to_infiles[sname] = [file]

            # names for output as 'kcount_{name}'
            self.names_to_outfiles[sname] = f"{self.prefix}_{sname}"
        # logger.debug(self.names_to_infiles)



    def run(self):
        """
        Calls kmc count on all files
        """
        for sname in self.names_to_infiles:

            # get files
            files = ofiles = self.names_to_infiles[sname]
            read1 = files[0]
            read2 = (None if len(files) == 1 else files[1])

            # trim reads
            trimstats = None
            if self.trim_reads:
                tool = ReadTrimming(
                    read1=read1, 
                    read2=read2,
                    workdir=self.workdir,
                    subsample=self.subsample_reads,
                )
                tool.run()
                trimstats = tool.parse_stats_from_json()
                files = [tool.tmp1, tool.tmp2]

            # run kmer counting
            kmcstats = self.call_kmc_count(files, sname)

            # store results to dataframe
            self.store_stats(ofiles, sname, kmcstats, trimstats)

            # cleanup tmp files from read trimming
            if self.trim_reads:
                tool.cleanup()



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

        # create command: 'kmc -k17 @filelist outname workdir'
        cmd = [
            self.kmc_binary, 
            "-ci{}".format(self.mindepth),
            "-cx{}".format(self.maxdepth),
            "-cs{}".format(self.maxcount),
            "-k{}".format(self.kmersize), 
            "@" + input_file,
            self.names_to_outfiles[sname],
            self.workdir,
        ]

        # call subprocess on the command
        out = subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )     

        # parse STDOUT to get the kmer stats
        os.remove(input_file)
        logger.debug(out.stdout.decode())
        return out.stdout.decode().split("Stats:")[-1]



    def store_stats(self, ofiles, sname, kmcstats, trimstats):
        """
        Store stats to dataframe.
        """
        # store results to the dataframe
        for line in kmcstats.strip().split("\n"):
            key, val = line.split(" : ")
            
            if key.strip() == "Total no. of k-mers":
                self.statsdf.loc[sname, "kmers_total"] = int(val.strip())

            if key.strip() == "No. of unique k-mers":
                self.statsdf.loc[sname, "kmers_unique"] = int(val.strip())

            if key.strip() == "No. of unique counted k-mers":
                self.statsdf.loc[sname, "kmers_unique_counted"] = int(val.strip())

            if key.strip() == "Total no. of reads":
                self.statsdf.loc[sname, "reads_total"] = int(val.strip())

            if key.strip() == "No. of k-mers above max. threshold":
                self.statsdf.loc[sname, "kmers_above_thresh"] = int(val.strip())

            if key.strip() == "No. of k-mers below min. threshold":
                self.statsdf.loc[sname, "kmers_below_thresh"] = int(val.strip())

        # save kmc db prefix path and original fastq names
        self.statsdf.loc[sname, "database"] = self.names_to_outfiles[sname]
        self.statsdf.loc[sname, "fastqs"] = ",".join([
            os.path.basename(i) for i in ofiles if i])

        # store trimming info
        if not trimstats:
            self.statsdf.loc[sname, "reads_passed_trimming"] = self.statsdf.loc[sname, "reads_total"]
        else:
            self.statsdf.loc[sname, "reads_total"] = trimstats[0]
            self.statsdf.loc[sname, "reads_passed_trimming"] = trimstats[1]

        # correct read counts for PE double-counting
        if len(self.names_to_infiles[sname]) > 1:
            self.statsdf.loc[sname, "reads_total"] = int(
                self.statsdf.loc[sname, "reads_total"] / 2)
            self.statsdf.loc[sname, "reads_passed_trimming"] = int(
                self.statsdf.loc[sname, "reads_passed_trimming"] / 2)

        # save database csv to the workdir
        path = os.path.join(self.prefix + ".csv")
        self.statsdf.to_csv(path)
        logger.info(f"new database: {sname}")



if __name__ == "__main__":

    # test dataset w/ R1 and R2 files
    # FILES = "~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R*_.fastq.gz"
    FILES = "../data/hybridus_*.fastq.gz"

    # example
    counter = Kcount(
        name="hyb", 
        workdir="/tmp/", 
        fastq_path=FILES, 
        kmersize=31, 
        name_split="_R",
        trim_reads=True,
        subsample_reads=1e6,
        mindepth=1,
        maxdepth=1e9,
        maxcount=255,
    )
    counter.run()

    # statdf is saved to the workdir as a CSV
    print(counter.statsdf.T)
