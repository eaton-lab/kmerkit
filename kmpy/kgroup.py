#!/usr/bin/env python


"""
Load a phenotypes CSV file, or dictionary in the API, 
to assign samples to groups, and then apply kmer_tools operations
like merge or intersect on kmers.
"""


import os
import sys
import subprocess
import pandas as pd
from loguru import logger
from kmpy.utils import KmpyError


# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


COMPLEX = """\n
INPUT:
{input_string}

OUTPUT:
{output_string}

OUTPUT_PARAMS:
{output_params}
"""


class Kgroup:
    """
    Perform operation on two groups of samples grouped by a trait (column)
    in the phenotype data. Kmer counts within each group are combined
    using operators 'union' or 'intersect' for groups g0 and g1. Then the 
    operator subtract ... is used to compare the two groups. Example:
    
    # kmers in all samples A,B,C (intersect) but not any of D,E,F (union)
    (A * B * C) - (D + E + F) 

    Parameters
    ----------
    name (str):
        A name prefix for output files.
    workdir (str):
        ...
    trait (str):
        A column name from the phenos file containing binary data (0/1)
        that will be used assign samples to groups 0 and 1 depending on
        their value at this trait.
    operation_g0 (str):
        Options are "union" or "intersect". 
    operation_g1 (str):
        Options are "union" or "intersect". 
    operation_g0g1 (str):
        Options are "subtract" or "counters_subtract". See the reverse
        option to reverse which group is subtracted from which.
    mindepth_g0 (int):
        ...
    mindepth_g1 (int):
        ...
    mindepth_g0g1 (int):
        ...
    ...

    reverse (bool):
        The default behavior is to subtract kmers in group 1 from group 0,
        such that you are left with kmers in group 0 but not group 1.
        If alternatively you want to find kmers in group 1 but not group 0
        then you can set reverse to True.
    force (bool):
        Ignore existing files and overwrite using <name>.


    Returns
    ----------
    None. This function produces output files in the working directory
    which was specified in Kcount and loaded by Kgroup.
    """
    def __init__(
        self, 
        name, 
        workdir, 
        phenos,
        trait, 
        operation_g0="union", 
        operation_g1="intersect", 
        operation_g0g1="subtract", 
        mindepth_g0=5,
        mindepth_g1=5,
        mindepth_g0g1=None,
        maxdepth_g0=1000,
        maxdepth_g1=1000,
        maxdepth_g0g1=None,
        # maxcount_g0=1000,
        # maxcount_g1=1000,
        # maxcount_g0g1=None,
        reverse=False,
        force=False,
        ):

        # store parameters
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.phenos = os.path.realpath(os.path.expanduser(phenos))
        self.kcpath = os.path.join(self.workdir, f"kcount_{name}.csv")
        
        # output prefix for .complex, subtract, countsubtract, intersect
        self.prefix = os.path.join(self.workdir, "kgroup")

        # attributes to be filled
        self.samples = []
        self.statsdf = None

        # get the kmctools binary for kmer set comparisons
        self.kmctools_binary = os.path.join(sys.prefix, "bin", "kmc_tools")
        logger.debug("using KMC binary: {}".format(self.kmctools_binary))

        # loads statsdf and phenos
        self.load_count_csv()
        self.load_phenos()

        # default parameters for calling .run_complex()
        self.params = {
            'name': name, 
            'trait': trait,
            'operation_g0': operation_g0,
            'operation_g1': operation_g1,
            'operation_g0g1': operation_g0g1,
            'mindepth_g0': mindepth_g0,
            'mindepth_g1': mindepth_g1,
            'mindepth_g0g1': mindepth_g0g1,
            'maxdepth_g0': maxdepth_g0,
            'maxdepth_g1': maxdepth_g1,
            'maxdepth_g0g1': maxdepth_g0g1,
            # 'maxcount_g0': 1000,
            # 'maxcount_g1': 1000,
            # 'maxcount_g0g1': 1000,
            'reverse': reverse,
            'force': force,
        }



    def load_count_csv(self):
        """
        load CSV results from kcount to get sample names and stats
        """
        self.statsdf = pd.read_csv(self.kcpath, index_col=0)



    def load_phenos(self):
        """
        Parse phenotype input that is either a dict or CSV, 
        check all names for match with database files in workdir.
        """
        # load the kcounts database to get all sample names
        self.phenodf = pd.read_csv(self.phenos, index_col=0)

        # check that names in statsdf match those in phenodf
        setp = set(self.phenodf.index)
        sets = set(self.statsdf.index)

        # raise an error if no names overlap
        if setp.isdisjoint(sets):
            dbnames = ", ".join(self.statsdf.index[:4].tolist())
            phnames = ", ".join(self.phenodf.index[:4].tolist())
            msg = (
                f"Sample names in pheno do not match names in database:\n"
                f"  DATABASE: {dbnames}...\n"
                f"  PHENOS:   {phnames}...\n"               
            )
            logger.error(msg)
            raise KmpyError(msg)

        # warning for samples only in pheno that are not in database
        onlypheno = setp.difference(sets)
        if onlypheno:
            logger.warning(
                f"Skipping samples in pheno not in database: {onlypheno}")

        # store overlapping samples as the test samples
        self.samples = sorted(setp.intersection(sets))



    def get_complex_input(self, group0, group1):
        """
        Builds part 1/3 of the complex file for a set of samples.
        """
        # build input section -------------------------------------------
        input_list = []

        # iterate over all samples
        for sname in group0.samples:

            # get fastq file path
            fpath = self.statsdf.loc[sname, "database_path"]

            # build input string
            cmd = [sname, "=", fpath]
            if self.params['mindepth_g0']:
                cmd.extend([f"-ci{self.params['mindepth_g0']}"])
            if self.params['maxdepth_g0']:
                cmd.extend([f"-cx{self.params['maxdepth_g0']}"])
            # if self.params['maxcount_g0']:
                # cmd.extend([f"-cs{self.params['maxcount_g0']}"])
            input_list.append(" ".join(cmd))

        for sname in group1.samples:

            # get fastq file path
            fpath = self.statsdf.loc[sname, "database_path"]

            # build input string
            cmd = [sname, "=", fpath]
            if self.params['mindepth_g1']:
                cmd.extend([f"-ci{self.params['mindepth_g1']}"])
            if self.params['maxdepth_g1']:
                cmd.extend([f"-cx{self.params['maxdepth_g1']}"])
            # if self.params['maxcount_g1']:
                # cmd.extend([f"-cs{self.params['maxcount_g1']}"])
            input_list.append(" ".join(cmd))

        # final spacer and return
        return "\n".join(input_list)



    def get_complex(self):
        """
        Builds the 'complex' input string (see COMPLEX global above)
        to call complex operations using kmer_tools using user params
        entered in self.params. This is called within self.run_complex()       
        """
        # get groups
        group0 = Group(self.samples, self.phenodf, trait=0)
        group1 = Group(self.samples, self.phenodf, trait=1)

        # INPUT: get sample names and filters
        input_str = self.get_complex_input(group0, group1)

        # OUTPUT: (group1) - (group0) 
        string0 = group0.get_string(self.params['operation_g0'])
        string1 = group1.get_string(self.params['operation_g1'])
        outname = self.prefix + "_" + self.params['name']
        oper = ("-" if self.params["operation_g0g1"] == "subtract" else "~")
        if self.params['reverse']:
            output_str = f"{outname} = ({string1}) {oper} ({string0})"
        else:
            output_str = f"{outname} = ({string0}) {oper} ({string1})"

        # OUTPUT_PARAMS: ...
        opt = []
        if self.params['mindepth_g0g1']:
            opt.extend([f"-ci{self.params['mindepth_g0g1']}"])
        if self.params['maxdepth_g0g1']:
            opt.extend([f"-cx{self.params['maxdepth_g0g1']}"])
        # if self.params['maxcount_g0g1']:
            # opt.extend([f"-cx{self.params['maxcount_g0g1']}"])
        oparams_str = " ".join(opt)

        # Build complex string and print to logger
        complex_string = COMPLEX.format(**{
            'input_string': input_str,
            'output_string': output_str,
            'output_params': oparams_str,
        })
        logger.debug(complex_string)

        return complex_string



    def run(self):
        """
        Create a 'complex' input file to run `kmc_tools complex ...`
        """
        # check that output with 'name' does not already exist.
        complex_file = self.prefix + "_" + self.params['name']
        if not self.params['force']:
            assert not os.path.exists(complex_file), (
            "output already exists for {name}, use force=True to overwrite")

        # check that 'trait' is in pheno.
        assert self.params['trait'] in self.phenodf.columns, (
            f"trait {self.params['trait']} not in phenos")

        # get complex string
        complex_string = self.get_complex()

        # write to a file
        with open(complex_file + "_complex.txt", 'w') as out:
            out.write(complex_string)

        # cmd: 'kmer_tools [global_params] complex <operations file>'
        cmd = [self.kmctools_binary, "complex", complex_file + "_complex.txt"]

        # call subprocess on the command
        out = subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )




    def dump(self, mindepth=1, maxdepth=None, maxcount=None):
        """
        Dump kmers to text for a specified database. This is meant to be 
        run after run_complex and to use the same params['name'].

        kmc_tools transform /tmp/name dump /tmp/name_kmers.txt
        """

        # if counts is False then write 'compact' database
        cmd = [
            self.kmctools_binary, "transform", 
            self.prefix + "_" + self.params['name'],
            "dump",
            self.prefix + "_" + self.params['name'] + "_kmers.txt",
        ]

        if mindepth:
            cmd.extend(["-ci{}".format(mindepth)])
        if maxdepth:
            cmd.extend(["-cx{}".format(maxdepth)])
        if maxcount:
            cmd.extend(["-cs{}".format(maxcount)])

        # call subprocess on the command
        logger.debug(" ".join(cmd))
        subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )

        # count stats and report to logger
        kfile = self.prefix + "_" + self.params['name'] + "_kmers.txt"
        with open(kfile, 'r') as indat:
            nkmers = sum(1 for i in indat)
        logger.info("{} kmers dumped to {}".format(nkmers, kfile))


class Group:
    """
    Stores subset of samples with a given trait value, and which are 
    present in both databases. Can return sample list or a string for
    computing the union of kmers in kmer_tools format: (A + B + C + D)
    
    samples: first filter, only samples we are considering.
    phenodf: used to check trait values.
    trait: used to select trait value for filtering.
    """
    def __init__(self, samples, phenodf, trait):
        # store inputs
        self.dbsamples = samples
        self.phenodf = phenodf
        self.trait = trait

        # to be filled
        self.samples = []
        self.ustring = ""
        self.istring = ""        

        # filler funcs
        self.subselect()       
        self.get_union_string()
        self.get_intersect_string()        


    def subselect(self):
        """
        Returns a list of samples
        """
        # get sample names in phenodf that have trait=1
        self.samples = self.phenodf[self.phenodf.trait == self.trait]

        # only keep if samples are in the database too
        self.samples = [i for i in self.samples.index if i in self.dbsamples]


    def get_union_string(self):
        """
        Builds the union string of kmers in any samples
        The counter shows the sum count of the kmer across all samples.
        """
        self.ustring = " + ".join(self.samples)
        logger.debug(self.ustring)


    def get_intersect_string(self):
        """
        Builds the intersect string of kmers in all samples. 
        The counter shows the min count of the kmer is any one sample.
        """
        self.istring = " * ".join(self.samples)
        logger.debug(self.istring)


    def get_string(self, arg):
        """
        Returns either the union or intersect string depending on arg
        """
        if arg == "union":
            return self.ustring
        return self.istring





if __name__ == "__main__":

    # first run: python3 kcount.py 

    # fake data
    ACC = ["1A_0", "1B_0", "1C_0", "1D_0", "2E_0", "2F_0", "2G_0", "2H_0"]
    TRAITS = [0, 0, 0, 0, 1, 1, 1, 1] 

    # build dataframe
    PHENOS = pd.DataFrame(
        index=ACC,
        columns=["trait"],
        data=TRAITS,
    )
    PHENOS.to_csv("/tmp/phenos.csv")


    # load database with phenotypes data
    kgp = Kgroup(
        name="test",
        workdir="/tmp",
        phenos="/tmp/phenos.csv",
        trait='trait',
        operation_g0="union",
        operation_g1="intersect", 
        operation_g0g1="subtract",
        mindepth_g0=5,
        mindepth_g1=5,
        mindepth_g0g1=None,
        maxdepth_g0=100,
        maxdepth_g1=100,
        maxdepth_g0g1=None,
        # maxcount_g0=1000,
        # maxcount_g1=1000,
        # maxcount_g0g1=None,
        reverse=True,
        force=True,
    )

    # dump the kmers to a file
    kgp.run()
    kgp.dump(mindepth=5)
