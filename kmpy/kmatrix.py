#!/usr/bin/env python

"""
Takes kmer databases for N samples and constructs a matrix
of (nsamples, nkmers) of either count or binary data.


TODO: 
    - support count (non-binary) matrix
    - decide if there is need for different min,max settings beyond global?
        - maybe minmap setting could apply to intersect groups, and 
          this could be implemented by setting depths to one and then
          -ci to minmap[i].

"""

import os
import sys
import subprocess
import numpy as np
import pandas as pd
from loguru import logger
from kmpy.utils import KmpyError, Group, COMPLEX



class Kmatrix:
    """
    Constructs a matrix of (nsamples x nkmers) as an np.array.

    Finds databases for counted samples from the path to the 
    counts CSV file. In many cases the matrix is unlikely to fit into 
    memory, so we will use memory-mapping to write to the file 
    efficiently.

    Parameters
    ==========
    name (str):
        Name prefix for output files, should be same as used in kcount.
    workdir (str):
        Working directory containing kmer database files from kcount,
        and the location where new outputs will be saved.
    mindepth (int):
        ...
    maxdepth (int):
        ...
    maxcount (int):
        ...
    counts (bool):
        Fill the genotype matrix with kmer counts, w/ upper limit of 
        maxcount. If False then data are stored binary (presence|absence).
    normalize (bool):
        Normalize count data by number of reads in each sample. Only relevant
        to counts=True data.
    subsample (list):
        An optional list of sample names to subsample which taxa are included
        in the genotype matrix. Default is None which means to use all samples
        in the <workdir>/kcount_<name>.csv dataframe.
    """
    def __init__(
        self, 
        name, 
        workdir, 
        mindepth=1, 
        maxdepth=1e9, 
        maxcount=1024, 
        counts=False, 
        normalize=True,
        subsample=None,
        ):

        # store params
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.prefix = os.path.join(self.workdir, f"kmatrix_{self.name}")
        self.kcpath = os.path.join(self.workdir, f"kcount_{name}.csv")
        self.mindepth = int(mindepth)
        self.maxdepth = int(maxdepth)
        self.maxcount = int(maxcount)
        self.counts = counts
        self.normalize = normalize
        self.subsample = subsample

        # binary path
        self.kmctools_binary = os.path.join(sys.prefix, "bin", "kmc_tools")

        # paths to files 
        self.statsdf = ""
        self.load_count_csv()

        # which samples to use, maps {names: dbprefix}
        self.subsample = sorted(
            self.subsample if self.subsample is not None 
            else self.statsdf.index.tolist()
        )
        self.names_to_db = dict(zip(self.statsdf.index, self.statsdf.database))
        self.nsamples = len(self.subsample)

        # attributes to be filled.
        self.nkmers = 0
        self.matrix = None  # np.zeros((nsamples, nkmers))



    def load_count_csv(self):
        """
        load CSV results from kcount to get sample names and stats
        """
        self.statsdf = pd.read_csv(self.kcpath, index_col=0)


    def get_complex_input(self, double=False):
        """
        Builds part 1/3 of the complex file for a set of samples.
        """
        input_list = []

        # iterate over all samples
        for sname in self.subsample:

            # get fastq file path
            fpath = self.statsdf.loc[sname, "database"]

            # build input string
            cmd = [
                sname, "=", fpath,
                f"-ci{self.mindepth}",
                f"-cx{self.maxdepth}",
                # f"-cs{self.maxcount}",
            ]
            input_list.append(" ".join(cmd))

            # add double hack trick
            if double:
                cmd[0] = f"null_{cmd[0]}"
                input_list.append(" ".join(cmd))

        # final spacer and return
        return "\n".join(input_list)


    def run_complex(self, oper="union"):
        """
        Build a complex file and call it from kmer_tools w/ subprocess.
        If not union then intersect is called.
        """
        # union or intersect
        oper = ("union" if oper=="union" else "intersect")

        # a = /tmp/path -ci1 -cx1000000
        input_str = self.get_complex_input()#double=True)

        # name = (a + b + c + d), or name = (a * b * c * d)
        allsamples = Group(self.subsample, cmode=None)#, double=True)
        operstr = allsamples.get_string(oper)
        outname = self.prefix + "_" + oper
        output_str = f"{outname} = ({operstr})"

        # outparams
        outparams_str = " ".join([
            f"-ci{self.mindepth}",
            f"-cx{self.maxdepth}",
            # f"-cs{self.maxcount}",
        ])

        # format for complex
        complex_string = COMPLEX.format(**{
            'input_string': input_str,
            'output_string': output_str,
            'output_params': outparams_str,
        })

        # write to a file
        complex_file = self.prefix + "_complex.txt"
        with open(complex_file, 'w') as out:
            out.write(complex_string)

        # cmd: 'kmer_tools [global_params] complex <operations file>'
        # logger.debug(complex_string)
        cmd = [self.kmctools_binary, "complex", complex_file]

        # call subprocess on the command
        out = subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )
        logger.info(f"new database: {os.path.split(self.prefix)[-1]}_{oper}")
        # logger.warning(out.stdout.decode())



    def fill_count_mat(self):
        """
        For each sample check all kmers against var_kmers and fill 1 
        in the matrix if kmer is present in the sample. Code here does
        not use kmc_tools b/c it was a pain to hack a working solution
        using subtract methods when it removes kmers with 0 counts, or 
        runs into maxcounts.

        TODO: only fills binary for now. To hack for storing counts we will
        need to handle comparisons that include maxcount edge limits.
        """
        for sidx, sname in enumerate(self.subsample):
            pass



    def fill_binary_mat(self):
        """
        For each sample check all kmers against var_kmers and fill 1 
        in the matrix if kmer is present in the sample. Code here does
        not use kmc_tools b/c it was a pain to hack a working solution
        using subtract methods when it removes kmers with 0 counts, or 
        runs into maxcounts.
        """
        for sidx, sname in enumerate(self.subsample):

            # name of new diff database to be generated
            sampledb = self.names_to_db[sname]

            # dump to kmers file
            self.dump(sampledb, count_kmers=True, write_counts=False)            

            # open file handles
            sampio = open(sampledb + "_kmers.txt", 'r')
            vario = open(self.prefix + "_var_kmers.txt", 'r')

            # get first kmer in each file
            kmi = next(sampio)
            kmu = next(vario)
            # logger.info(f"{kmi.strip()} | {kmu.strip()} | {kmi == kmu}")

            # loop until both iters are exhausted
            idx = 0
            while 1:

                try:
                    # lowest is in sample: advance kmi
                    if kmi < kmu:
                        kmi = next(sampio)

                    # lowest is in var: advance kmu
                    elif kmu < kmi:
                        kmu = next(vario)
                        idx += 1

                    # they match, record it!
                    else: 
                        self.matrix[sidx, idx] = 1
                        kmi = next(sampio)
                        kmu = next(vario)
                        idx += 1
                        # logger.info(f"{idx} {kmu.strip()} {kmi.strip()}")

                except StopIteration:
                    break

            # close file handles
            logger.debug(f"{sname} var kmers: {self.matrix[sidx, :].sum()}")
            sampio.close()
            vario.close()

            # cleanup tmp sample kmer file
            os.remove(sampledb + "_kmers.txt")



    def get_var_kmers(self):
        """
        Compares union to intersect of all kmers to get variable kmers, 
        and removes the union and intersect databases.
        """  
        # build command to subtract union from intersect
        cmd = [
            self.kmctools_binary, 
            "simple",
            self.prefix + "_union",      # all kmers (any sample)
            self.prefix + "_intersect",  # some kmers (shared by all)
            "kmers_subtract",            # present in first, absent in second
            self.prefix + "_var",
        ]

        # call subprocess on the command
        subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )



    def run(self):
        """
        Infer and save a kmer genotypes matrix file containing 
        presence/absence of each kmer in each sample (excluding kmers
        that are not variable among samples).
        """
        # get union of kmer counts across all subsamples
        self.run_complex(oper="union")
        self.dump(self.prefix + "_union", count_kmers=True)

        # get intersect of kmer counts across all subsamples
        self.run_complex(oper="intersect")
        self.dump(self.prefix + "_intersect", count_kmers=True)

        # get kmers that are variable (union - intersect)
        self.get_var_kmers()
        self.nkmers = self.dump(
            self.prefix + "_var", 
            count_kmers=True, 
            write_counts=False,
        )

        # build empty array
        self.matrix = np.memmap(
            filename=self.prefix + "_var_genos.npy",
            mode="w+", 
            shape=(self.nsamples, self.nkmers),
            dtype=np.bool_, 
        )
        logger.info(f"matrix shape: {self.matrix.shape}")

        # fill each sample into matrix
        self.fill_binary_mat()

        # remove temp files
        self.cleanup()



    def cleanup(self):
        """
        remove temp files created during analysis.
        """
        prefixs = [
            self.prefix,
            self.prefix + "_union",
            self.prefix + "_intersect",
        ]
        suffixes = [
            ".kmc_pre",
            ".kmc_suf",
            "_complex.txt",
            "_kmers.txt"
        ]
        for pre in prefixs:
            for post in suffixes:
                if os.path.exists(pre + post):
                    os.remove(pre + post)



    def dump(self, database, count_kmers=False, write_counts=True):
        """
        Dump kmers to text for a specified database. 
        calls: 'kmc_tools transform /tmp/name dump /tmp/name_kmers.txt'

        Parameters
        ----------
        ...
        count_kmers (bool):
            If True then returns a int with n unique kmers after applying 
            the minmax filters.
        write_counts (bool):
            If False then the counts are excluded from dump txt file.

        """
        # if counts is False then write 'compact' database
        cmd = [
            self.kmctools_binary,
            "transform", 
            database,
            "-ci{}".format(self.mindepth),
            "-cx{}".format(self.maxdepth),
            # "-cs{}".format(maxcount),
            "dump",
            "-s",
            database + "_kmers.txt",
            "-ci{}".format(self.mindepth),
            "-cx{}".format(self.maxdepth),
            "-cs{}".format(self.maxcount),           
        ]

        # call subprocess on the command
        subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )

        # overwrite with file cutting first column only
        if not write_counts:
            cmd = ["cut", "-f1", database + "_kmers.txt"]
            with open(database + "_kmers.tmp", "w") as ofile: 
                subprocess.run(cmd, stdout=ofile, check=True)
            os.rename(database + "_kmers.tmp", database + "_kmers.txt")

        # count and report number of kmers actually written (i.e., w/ filters)
        if count_kmers:
            kfile = database + "_kmers.txt"
            with open(kfile, 'r') as indat:
                nkmers = sum(1 for i in indat)
            logger.info(f"{os.path.split(database)[-1]} dump: {nkmers} kmers")
            return nkmers
        return 0





if __name__ == "__main__":

    kma = Kmatrix(
        name="hyb", 
        workdir="/tmp", 
        mindepth=1, 
        maxdepth=1e9, 
        maxcount=255, 
        counts=False, 
        normalize=False,
        subsample=None,
    )

    kma.run()