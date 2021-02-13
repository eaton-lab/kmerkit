#!/usr/bin/env python

"""
Kcount -> Kfilter

Apply filters to kmer sets to find a target set of kmers.

Loads a phenotypes CSV file and trait column (or imap dictionary in the API)
to assign samples to groups; then apply kmer_tools operations
to filter kmers for inclusion in the dataset based on:
    - mincov: minimum frequency across all samples
    - minmap: dict of min freq in each group (assigned by phenodf.trait)
    - mincov_canon: minimum frequency across all samples of a kmer
        being in present in both directional forms.

TODO:
    - canon shouldn't apply to all samples, It should apply
      to only those samples for which the kmer is present...

    - b/c canon is a bit time and disk consuming, we could apply it
      after the other filters, so we only have to focus on a reduced
      set of kmers. Would this be useful? ... Still need to start by 
      counting all non-con kmers in each sample...

Note: coverage in these filters refers to the number or frequency
of samples in which the kmer is present. It does not refer to the 
count (depth) of the kmer in one or more samples. To filter kmers
based on counts you can use the 'mincount' args in kmerkit.Kcount().
"""

import os
import shutil
import subprocess
import numpy as np
import pandas as pd
from loguru import logger
from kmerkit.kcount import Kcount
from kmerkit.kmctools import KMTBIN, dump, info
from kmerkit.utils import Group, COMPLEX, KmerkitError


# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes


class Kfilter:
    """
    Apply filters to find a target set of kmers from a set of samples
    assigned to groups in a phenos file. Filters apply to coverage
    (presence) across samples, and to their occurrence in canonized 
    form. Presence/absence of kmers here depends on the cutoffs 
    (e.g., mindepth) used in kcount.

    Parameters
    ----------
    name (str):
        A name prefix for output files.
    workdir (str):
        Directory for output files, and where kmer database files are
        currently located (and the .csv from kcount).
    phenos (str):
        Path to a CSV formatted phenotypes file with sample names as
        index and trait names as columns. Values of 'trait' should be 
        binary, 0 and 1 are used to assign samples to groups g0 or g1.
    trait (str):
        A column name from the phenos file to use for this analysis. 
        Values should be binary (0/1).    
    mincov (int, float):
        Global minimum coverage of a kmer across all samples as an
        integer or applied as a proportion (float).
    mincov_canon (int, float):
        Minimum coverage of a kmer across all samples where it must
        occur in both directions (--> and comp(<--)), as an integer 
        or applied as a proportion (float).
    minmap (dict):
        Minimum coverage of a kmer across samples in each group (trait
        category). Keys of the dict should be trait values/categories
        and values of the dict should be floats representing the 
        minimum proportion of samples in the group for which the kmer
        must be present. Examples:
        minmap = {0: 0.0, 1: 0.9} or minmap={"red": 0.9, "blue": 0.0}.
    maxmap (dict):
        Maximum coverage of a kmer across samples in each group (trait
        category). Keys of the dict should be trait values/categories
        and values of the dict should be floats representing the 
        maximum proportion of samples in the group for which the kmer
        can be present. Examples:
        maxmap = {0: 0.0, 1: 1.0} or maxmap={"red": 1.0, "blue": 0.0}.


    Returns
    ----------
    None. Writes kmc binary database to prefix <workdir>/kfilter_{name}
    """
    def __init__(
        self, 
        name, 
        workdir, 
        phenos,
        trait, 
        mincov=0.0,
        mincov_canon=0.25,
        minmap=None,
        maxmap=None,
        ):

        # store file parameters
        self.name = name
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.kcpath = os.path.join(self.workdir, f"kcount_{name}.csv")
        self.phenos = os.path.realpath(os.path.expanduser(phenos))
        self.trait = trait

        # store filter params
        self.mincov = mincov
        self.mincov_canon = mincov_canon
        self.minmap = (minmap if minmap is not None else {})
        self.maxmap = (maxmap if maxmap is not None else {})
        self._minmax_empty = self.minmap == self.maxmap == {}
        
        # output prefix for .complex, subtract, countsubtract, intersect
        self.prefix = os.path.join(self.workdir, f"kfilter_{self.name}")

        # attributes to be filled. Samples is only those in kcounts & phenos.
        self.samples = []
        self.phenodf = None
        self.kcountdf = None

        # fills .kcountdf, .phenodf and .samples
        self.load_dataframes()

        # checks minmap against .phenodf.trait and .samples
        # and converts mincov and minmap to ints
        self.check_filters()



    def load_dataframes(self):
        """
        load kcount and phenotype dataframes to get sample names, 
        group info, and stats. Check all names for matching.
        """
        self.kcountdf = pd.read_csv(self.kcpath, index_col=0)
        self.phenodf = pd.read_csv(self.phenos, index_col=0)

        # check that 'trait' is in pheno.
        assert self.trait in self.phenodf.columns, (
            f"trait {self.trait} not in phenos")

        # check that names in kcountdf match those in phenodf
        setp = set(self.phenodf.index)
        sets = set(self.kcountdf.index)

        # raise an error if no names overlap
        if setp.isdisjoint(sets):
            dbnames = ", ".join(self.kcountdf.index[:4].tolist())
            phnames = ", ".join(self.phenodf.index[:4].tolist())
            msg = (
                "Sample names in pheno do not match names in database:\n"
                f"  DATABASE: {dbnames}...\n"
                f"  PHENOS:   {phnames}...\n"               
            )
            logger.error(msg)
            raise KmerkitError(msg)

        # warning for samples only in pheno that are not in database
        onlypheno = setp.difference(sets)
        if onlypheno:
            logger.warning(
                f"Skipping samples in pheno but not in database: {onlypheno}")

        # store overlapping samples as the test samples
        self.samples = sorted(setp.intersection(sets))

        # keep only rows in pheno that are in sa
        self.phenodf = self.phenodf.loc[self.samples]



    def check_filters(self):
        """
        Check that minmap keys are the same as values in .phenodf.trait, 
        and convert values to ints from floats (proportions).

        Default minmap is 0, and maxmap is 1.0.
        """
        # int encode mincov
        if isinstance(self.mincov, float):
            self.mincov = int(np.floor(self.mincov * len(self.samples)))

        # get the values of 'trait' in phenos
        values = set(self.phenodf[self.trait])

        # iterate over value categories assigning min and max
        for val in values:

            # set value if user did not provide one
            if val not in self.minmap:
                self.minmap[val] = 0.0
            if val not in self.maxmap:
                self.maxmap[val] = 1.0

            # check if there are ANY samples with this trait
            mask = self.phenodf[self.trait] == val
            nsamps = self.phenodf.loc[mask].shape[0]

            # check that minmap key is in 'trait' values
            assert nsamps, (
                "Keys in minmap and maxmap should match 'trait' values. "
                "No samples present in both kcount database and the "
                f"phenos database have {self.trait} == {val}:\n"
                f"{self.phenodf[self.trait]}"
            )

            # int encode minmax
            if isinstance(self.minmap[val], float):
                asint = int(np.floor(self.minmap[val] * nsamps))
                self.minmap[val] = asint
                asint = int(np.floor(self.maxmap[val] * nsamps))                
                self.maxmap[val] = asint

        logger.debug(f"int encoded minmap: {self.minmap}")
        logger.debug(f"int encoded maxmap: {self.maxmap}")



    def call_complex(self, dbdict, mindepth, maxdepth, oper, outname):
        """
        Builds the 'complex' input string (see COMPLEX global)
        to call complex operations using kmer_tools. Samples are 
        selected from dbdict, a dict of {idx: dbprefix-path, ...}. 
        See kgroup for complex operations using substract, etc.

        NB: KMC complex cannot handle dashes in linked names, 
        or keywords like min,max,diff,etc. Much safer to just use ints.
        """
        # INPUT: get sample names and filters
        input_list = []
        for sname, database in dbdict.items():

            # build input string
            cmdstr = f"{sname} = {database} -ci1 -cx1000000000"
            input_list.append(cmdstr)
        input_str = "\n".join(input_list)

        # OUTPUT: (a + b + c)
        outname = self.prefix + "_" + outname
        group = Group([str(i) for i in dbdict.keys()])
        output_str = f"{outname} = {group.get_string(oper)}"

        # OUTPUT_PARAMS: -ci0
        oparams_str = f"-ci{mindepth} -cx{maxdepth}" 

        # Build complex string and print to logger
        complex_string = COMPLEX.format(**{
            'input_string': input_str,
            'output_string': output_str,
            'output_params': oparams_str,
        })
        logger.debug(complex_string)

        # write to a tmp file
        complex_file = self.prefix + "_complex.txt"
        with open(complex_file, 'w') as out:
            out.write(complex_string)

        # cmd: 'kmer_tools [global_params] complex <operations file>'
        cmd = [KMTBIN, "complex", complex_file]

        # call subprocess on the command
        out = subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            cwd=self.workdir,
        )
        logger.info(f"new database: {os.path.basename(outname)}")
        os.remove(complex_file)



    def filter_canon(self):
        """
        Filter kmers that tend to only occur in one form or the
        other, likely due to adapters. 

        Returns:
        None. 
        Writes kmc database prefix=<workdir>/kfilter_{name}_canon-filtered
        """

        # database file for storing counts of kmers being present in both forms
        bothdb = self.prefix + "_canonsums"

        # iterate over samples
        for sidx, sname in enumerate(self.samples):

            # get fastq dict w/ fastqs used in kcounts
            fastq_dict = {sname: self.kcountdf.at[sname, "fastqs"].split(",")}

            # get options used in kcounts, but not mindepth (must be 1)
            kcount_params = {
                'workdir': self.workdir,
                'name': 'tmp-canon',
                'fastq_dict': fastq_dict,
                'kmersize': self.kcountdf.at[sname, "kmersize"],
                'subsample_reads': self.kcountdf.at[sname, "subsample_reads"],
                'trim_reads': self.kcountdf.at[sname, "trimmed"],
                'mindepth': 1, 
                'maxdepth': self.kcountdf.at[sname, "maxdepth"],
                'maxcount': self.kcountdf.at[sname, "maxcount"],
                'canonical': 0,
            }

            # count non-canon kmers in sample with same kcount args except -ci1
            logger.info(f"counting non-canonical kmers [{sname}]")
            kco = Kcount(**kcount_params)
            kco.run()
            nonc_db = os.path.join(self.workdir, f"kcount_tmp-canon_{sname}")

            # subtract counts of non-canon from canon, only keep if -ci1, 
            # these are the kmers observed occurring both ways in this sample.
            canon_db = self.kcountdf.at[sname, "database"]
            cmd = [
                KMTBIN,
                "-hp",                
                "simple",
                canon_db, 
                nonc_db, 
                "counters_subtract",
                self.prefix + "_tmp1",  # tmp1=kmers occurring both ways
                "-ci1"
            ]
            logger.debug(" ".join(cmd))            
            subprocess.run(cmd, check=True)

            # set count =1 on a tmp copy of bothways kmers db
            cmd = [
                KMTBIN,
                "-hp",                
                "transform",
                self.prefix + "_tmp1",
                "set_counts", 
                "1",
                self.prefix + "_tmp2",  # tmp2=kmers w/ count=1
            ]
            logger.debug(" ".join(cmd))
            subprocess.run(cmd, check=True)

            # if the first sample then simply copy tmp2 to BOTHDB
            if not sidx:
                for suff in [".kmc_suf", ".kmc_pre"]:                
                    shutil.copyfile(
                        self.prefix + "_tmp2" + suff,
                        bothdb + suff,
                    )

            # if not, then make copy of BOTHDB, and fill BOTHDB using kmctools
            else:
                for suff in [".kmc_suf", ".kmc_pre"]:                                
                    shutil.copyfile(
                        bothdb + suff,
                        self.prefix + "_tmp3" + suff,
                    )

                # fill bothdb as the sum union of tmp2 and tmp3
                cmd = [
                    KMTBIN,
                    "-hp",
                    "simple",
                    self.prefix + "_tmp2",   # kmers w/ count=1 from this samp
                    self.prefix + "_tmp3",   # kmers w/ count=sum so far
                    "union",
                    bothdb,
                    "-ci1",
                    "-cs{}".format(len(self.samples) + 1),
                ]
                logger.debug(" ".join(cmd))
                subprocess.run(cmd, check=True)              
                
            # cleanup
            logger.debug(f"removing database tmp-canon_{sname}")

        # get stats on full canonized set
        sumkmers = info(self.prefix + "_canonsums")

        # dump and calculate stats on canonsums counts
        dump(self.prefix + "_canonsums", write_kmers=False)
        with open(self.prefix + "_canonsums_kmers.txt", 'r') as indat:
            counts = np.loadtxt(indat, dtype=np.uint16)
            counts = counts / len(self.samples)
            mcanon = counts.mean()
            scanon = counts.std()
            del counts

        # report statistics on canonization
        logger.info(
            "kmers (proportion) in canonized form: "
            f"mean={mcanon:.2f}; std={scanon:.2f}"
        )

        # report a warning if too few are canonized:
        if mcanon < 0.1:
            logger.warning(
                "Very few kmers are canonized in this dataset. This may "
                "indicate your data is strand-specific (e.g., ddRAD) in "
                "which case you should set mincov_canon=0.0. Alternatively, "
                "if WGS reads, then your data may be very low coverage."
            )            

        # filter BOTHDB to get kmers in mincanon prop. of samples (CANON)
        mincanon_asint = int(np.floor(self.mincov_canon * len(self.samples)))
        cmd = [
            KMTBIN, 
            "transform",
            bothdb,
            "reduce",
            self.prefix + "_mincanon-filter", 
            "-ci{}".format(mincanon_asint),
            "-cx1000000000",
            # "-cs255",        # 65K
        ]
        logger.debug(" ".join(cmd))
        subprocess.run(cmd, check=True)              

        # calculate and report statistics on filtered canonized set.
        sumfiltkmers = info(self.prefix + "_mincanon-filter")
        logger.info(f"kmers filtered by mincov_canon: {sumkmers - sumfiltkmers}")

        # clean up tmp files
        countpre = os.path.join(self.workdir, "kcount")
        os.remove(countpre + "_tmp-canon.csv")
        for sname in self.samples:
            os.remove(countpre + "_tmp-canon_" + sname + ".kmc_pre")
            os.remove(countpre + "_tmp-canon_" + sname + ".kmc_suf")            
        os.remove(self.prefix + "_canonsums_kmers.txt")
        for suffix in ["tmp1", "tmp2", "tmp3", "canonsums"]:
            os.remove(self.prefix + "_" + suffix + ".kmc_pre")
            os.remove(self.prefix + "_" + suffix + ".kmc_suf")            



    def run(self):
        """
        Create a 'complex' input file to run `kmc_tools complex ...`
        """
        # MINCANON FILTER ---------------------------------------------
        # get kmers passing the mincov_canon filter (canon-filtered)
        self.filter_canon()

        # PREP for MINCOV and MAPS ------------------------------------
        # create count=1 databases for each sample
        # TODO: we could support a lowdisk option that overwrites the
        # original rather than copy, as long as original is not needed.
        for sname in self.samples:
            cmd = [
                KMTBIN,
                "-hp",                
                "transform",
                self.kcountdf.at[sname, "database"],
                "set_counts", 
                "1",
                self.prefix + f"_count1_{sname}",
            ]
            logger.debug(" ".join(cmd))
            subprocess.run(cmd, check=True)

        # MINCOV FILTER ----------------------------------------------
        # get count=1 kmers for ALL samples in self.samples
        dbdict = {
            idx: self.prefix + f"_count1_{sname}" 
            for idx, sname in enumerate(self.samples)
        }

        # get kmers passing the mincov filter (mincov-filtered)
        self.call_complex(
            dbdict=dbdict,
            oper="union",
            mindepth=self.mincov,
            maxdepth=1000000000,
            outname='mincov-filter',
        )

        # MAP FILTERS -------------------------------------------------
        # get kmers passing the map filters (same keys in both)
        for group in self.minmap:

            # get dict of {sname: count-1-db} for group samples
            mask = self.phenodf[self.trait] == group
            dbdict = {
                idx: self.prefix + f"_count1_{sname}"
                for idx, sname in enumerate(self.phenodf[mask].index)
            }

            # get kmers occurring >= minmap[group] in this group
            self.call_complex(
                dbdict=dbdict, 
                oper="union",
                mindepth=self.minmap[group],
                maxdepth=self.maxmap[group],
                outname=f'map-{group}-filter'
            )

        # INTERSECTION OF FILTERED KMERS -----------------------------
        # get intersection of kmers passing all filters (kmers-filtered)
        dbdict = {
            0 : self.prefix + "_mincanon-filter",
            1 : self.prefix + "_mincov-filter"
        }
        idx = 2
        for group in self.minmap:
            dbdict[idx] = f"{self.prefix}_map-{group}-filter"
            idx += 1
        logger.debug(dbdict)
        self.call_complex(
            dbdict=dbdict,
            oper="intersection",
            mindepth=1,
            maxdepth=1000000000,
            outname="filtered",
        )

        # TODO: log a summary of the filtered kmers
        # ...

        # cleanup tmp files
        os.remove(self.prefix + "_mincov-filter" + ".kmc_pre")
        os.remove(self.prefix + "_mincov-filter" + ".kmc_suf")
        os.remove(self.prefix + "_mincanon-filter" + ".kmc_pre")
        os.remove(self.prefix + "_mincanon-filter" + ".kmc_suf")
        for sname in self.samples:
            os.remove(self.prefix + f"_count1_{sname}" + ".kmc_pre")
            os.remove(self.prefix + f"_count1_{sname}" + ".kmc_suf")
        for group in self.minmap:
            os.remove(self.prefix + f"_map-{group}-filter" + ".kmc_pre")
            os.remove(self.prefix + f"_map-{group}-filter" + ".kmc_suf")



if __name__ == "__main__":

    # first run: python3 kcount.py 
    import kmerkit
    kmerkit.set_loglevel("DEBUG")

    # fake data
    PHENOS = "~/Documents/kmerkit/data/amaranths-phenos.csv"

    # load database with phenotypes data
    kgp = Kfilter(
        name="hybridus",
        workdir="/tmp",
        phenos=PHENOS,
        trait="fake",
        mincov=0.25,
        mincov_canon=1.0, #0.25,
        minmap={
            0: 0.0,
            1: 0.5,
        },
        maxmap={
            0: 0.0,
            1: 1.0,
        },
        #mapcov={0: (0.0, 0.0), 1: (0.5, 1.0)},
    )

    # dump the kmers to a file
    kgp.run()
