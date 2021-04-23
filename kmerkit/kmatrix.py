#!/usr/bin/env python

"""
Kcounts -> Kfilter -> Kmatrix

Takes a kmer count database from Kfilter and generates either a binary
matrix, or a count matrix, as a sample x kmer matrix.

Kmatrix removes any columns of the matrix that are invariant or empty 
(although these are usually already removed in Kfilter).

TODO: Think more about the intended goal/endpoint before developing 
further..., how often will kfilter,ktree, or kcount come before this
step?


Binary matrix (for gemma, for sklearn, look up proper input)
-------------
Apple   0000000
Banana  0000000
Cat     0011100
Dog     1111111...

Count matrix (for logistic regression or sklearn, should be normalized)
-------------
Apple,Banana,Cat,Dog
10,4,3,10
5,3,2,1
0,0,10,10
0,0,10,10
...

TODO: 
    - support count (non-binary) matrix
"""

import os
import subprocess
import numpy as np
from loguru import logger
from kmerkit.utils import KmerkitError, Group, COMPLEX
from kmerkit.kmctools import KMTBIN, info, dump
from kmerkit.kschema import Project, KmatrixParams, KmatrixData, KmatrixBase

# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


class Kmatrix:
    """
    Constructs a matrix of (nsamples x nkmers) as an np.array.

    Finds databases for counted samples from the path to the 
    counts CSV file. In many cases the matrix is unlikely to fit into 
    memory, so we will use memory-mapping to write to the file 
    efficiently.

    Parameters
    ==========
    json_file (str):
        Path to a kmerkit project JSON file.
    min_depth (int):
        ...
    max_depth (int):
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
        json_file,
        counts=False, 
        normalize=False,
        condense=False,
        ):

        self.project = Project.parse_file(json_file).dict()

        self.prefix = os.path.join(
            self.project['workdir'], f"{self.project['name']}_kmatrix")

        self.params = KmatrixParams(
            min_depth=min_depth,
            max_depth=max_depth,
            counts=counts,
            normalize=normalize,
            subsample=subsample,
        )

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
        cmd = [KMTBIN, "complex", complex_file]

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
            dump(sampledb, write_counts=False)

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
            KMTBIN, 
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
        #nkmers = info(self.prefix + "_union")
        dump(self.prefix + "_union")

        # get intersect of kmer counts across all subsamples
        self.run_complex(oper="intersect")
        #nkmers = info(self.prefix + "_intersect")
        dump(self.prefix + "_intersect")

        # get kmers that are variable (union - intersect)
        self.get_var_kmers()
        self.nkmers = info(self.prefix + "_var")
        dump(self.prefix + "_var", write_counts=False)

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