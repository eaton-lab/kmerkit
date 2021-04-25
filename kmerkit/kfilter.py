#!/usr/bin/env python

"""
Kcount -> Kfilter

Apply filters to kmer databases to create two sets of kmers, the 
FILTERED set and the KEEP set.

Note: coverage in these filters refers to the number or frequency
of samples in which the kmer is present. It does not refer to the 
count (depth) of the kmer in one or more samples. To filter kmers
based on counts you should use the 'min_depth' arg in the kcount step.

CLI Usage: 
------------------
kmerkit filter \
    --json /tmp/test.json \
    --traits traits.CSV \
    --min_map 0.0 1.0 \
    --max_map 0.1 1.0 \
    --min_map_canon 0.0 0.5 \
    --min_cov 5
"""

# TODO CANON FILTER:
#     - canon shouldn't apply to all samples, It should apply
#       to only those samples for which the kmer is present...

#     - b/c canon is a bit time and disk consuming, we could apply it
#       after the other filters, so we only have to focus on a reduced
#       set of kmers. Would this be useful? ... Still need to start by 
#       counting all non-con kmers in each sample...


import os
import itertools
import subprocess
import numpy as np

from loguru import logger
from kmerkit.kmctools import KMTBIN, info
from kmerkit.utils import Group, COMPLEX, KmerkitError
from kmerkit.kschema import KfilterParams, KfilterData, KfilterBase, Project

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
    json_file (str):
        A kmerkit project JSON file.
    traits (dict[str, int]):
        Map of binary trait values to samples names. Examples:
        traits = {0: ['a', 'b', 'c'], 1: ['d', 'e', 'f']}
    min_cov (int, float):
        Global minimum coverage of a kmer across all samples as an
        integer or applied as a proportion (float).
    min_cov_canon (int, float):
        Minimum coverage of a kmer across all samples where it must
        occur in both directions (--> and comp(<--)), as an integer 
        or applied as a proportion (float).
    min_map (dict):
        Minimum coverage of a kmer across samples in each group (trait
        category). Keys of the dict should be trait values/categories
        and values of the dict should be floats representing the 
        minimum proportion of samples in the group for which the kmer
        must be present. Examples:
            minmap = {0: 0.0, 1: 0.9}
    max_map (dict):
        Maximum coverage of a kmer across samples in each group (trait
        category). Keys of the dict should be trait values/categories
        and values of the dict should be floats representing the 
        maximum proportion of samples in the group for which the kmer
        can be present. Examples:
            maxmap = {0: 0.0, 1: 1.0}

    Returns
    ----------
    None. Writes kmc binary database to prefix <workdir>/kfilter_{name}
    """
    def __init__(self, json_file, traits, min_cov, min_map, max_map, min_map_canon):

        # load user inputs
        self.json_file = json_file
        self.project = Project.parse_file(json_file).dict()
        self.traits_to_samples = traits

        # will be subsampled down to matching taxa
        self.database = self.project['kcount']['data']
        self.samples_to_traits = {}

        # to be filled with results
        self.samples = {}

        # prefix for output files
        self.prefix = os.path.join(
            self.project['workdir'], f"{self.project['name']}_kfilter")

        # subsample traits to get union of kcount db and traits set
        self.select_samples()

        # use Serializable schema to perform type checking
        self._params = KfilterParams(
            min_cov=min_cov,
            min_map=min_map,
            max_map=max_map,
            min_map_canon=min_map_canon,
            trait_0=self.traits_to_samples[0],
            trait_1=self.traits_to_samples[1],
        )
        self.params = self._params.dict()        

        # converts filter params to integers and checks
        self.filters_to_ints()



    def select_samples(self):
        """
        Filter samples to those present in both database and traits_dict
        """
        # check that names in kcountdf match those in traits
        setp = set(self.database)
        sets = set(itertools.chain(*self.traits_to_samples.values()))
        revtraits = {}
        for key in self.traits_to_samples:
            for sample in self.traits_to_samples[key]:
                revtraits[sample] = key

        # raise an error if no names overlap
        if setp.isdisjoint(sets):
            db_names = ", ".join(setp)
            tr_names = ", ".join(sets)
            msg = (
                "Sample names in pheno do not match names in database:\n"
                f"  KMER_DATABASE: {db_names}...\n"
                f"  TRAITS:   {tr_names}...\n"               
            )
            logger.error(msg)
            raise KmerkitError(msg)

        # warning for samples only in traits that are not in database
        only_traits = setp.difference(sets)
        if only_traits:
            logger.warning(
                "Skipping samples in traits but not in kmer_database: "
                f"{only_traits}"
            )

        # store overlapping samples as the test samples
        shared_set = sorted(setp.intersection(sets))
        self.database = {i: self.database[i] for i in shared_set}
        self.samples_to_traits = {i: revtraits[i] for i in shared_set}
        self.traits_to_samples[0] = [
            i for (i, j) in self.samples_to_traits.items() if j == 0
        ]
        self.traits_to_samples[1] = [
            i for (i, j) in self.samples_to_traits.items() if j == 1
        ]



    def filters_to_ints(self):
        """
        Converts filter parameter float values to integers because that
        is what KMC requires.
        """
        # int encode min_cov
        if isinstance(self.params['min_cov'], float):
            self.params['min_cov'] = int(
                np.floor(self.params['min_cov'] * len(self.database)))
            logger.debug(f"int encoded min_cov: {self.params['min_cov']}")

        # int encoded map values
        for key in [0, 1]:
            nsamples = len(self.traits_to_samples[key])
            assert nsamples, f"No samples are set to state={key}"
            for param in ['min_map', 'max_map', 'min_map_canon']:
                value = int(np.floor(self.params[param][key] * nsamples))
                self.params[param][key] = value
                logger.debug(f"int encoded {param}[{key}]: {value}")



    def call_complex(self, dbdict, min_depth, max_depth, oper, out_name):
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
        outname = self.prefix + "_" + out_name
        group = Group([str(i) for i in dbdict.keys()])
        output_str = f"{outname} = {group.get_string(oper)}"

        # OUTPUT_PARAMS: -ci0
        oparams_str = f"-ci{min_depth} -cx{max_depth}"

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
        cmd = [KMTBIN, "-hp", "-t8", "complex", complex_file]

        # call subprocess on the command
        out = subprocess.run(
            cmd, 
            stderr=subprocess.STDOUT, 
            stdout=subprocess.PIPE,
            check=True,
            # cwd=self.workdir,
        )
        logger.info(f"new database: {os.path.basename(out_name)}")
        os.remove(complex_file)



    def get_all_single_counts(self):
        """
        PARALLEL.
        Prepare single count database of every sample which will be used
        for presence/absence set arithmetic in filters.
        """
        for sname in self.database:
            cmd = [
                KMTBIN, "-hp", "-t8",
                "transform",
                self.database[sname]['database'],
                "set_counts", "1",
                f"{self.prefix}_{sname}_count1",
            ]
            logger.debug(" ".join(cmd))
            subprocess.run(cmd, check=True)    


    def get_union_single_counts(self):
        """
        Prepare summed count database from all single count databases.
        This is the full set of observed kmers.
        """
        dbdict = {
            idx: f"{self.prefix}_{sname}_count1"
            for idx, sname in enumerate(self.database)
        }
        # get ALL kmers
        self.call_complex(
            dbdict=dbdict,
            oper="union",
            min_depth=0,
            max_depth=1000000000,
            out_name='union',
        )


    def get_min_cov_filtered_set(self):
        """
        Prepare database of kmers that did NOT occur across enough
        samples in total dataset.
        """
        dbdict = {
            idx: f"{self.prefix}_{sname}_count1"
            for idx, sname in enumerate(self.database)
        }        
        # get kmers NOT passing the mincov filter
        self.call_complex(
            dbdict=dbdict,
            oper="union",
            min_depth=0,
            max_depth=max(1, self.params['min_cov'] - 1),
            out_name='min_cov-filtered',
        )


    def get_group0_filtered_set(self):
        """
        Prepare database of kmers that did NOT pass the filters
        for group 0. This usually has a low min_map and a low max_map.
        """
        # get the group0 samples as a numbered dict
        dbdict = {
            idx: f"{self.prefix}_{sname}_count1"
            for idx, sname in enumerate(self.traits_to_samples[0])
        }

        # if EVERY kmer in this set should be filtered then use the union
        # this is just a faster way to accomplish the same as below.
        if self.params['max_map'][0] in [0, 1]:
            self.call_complex(
                dbdict=dbdict,
                oper="union",
                min_depth=0,
                max_depth=100000000,
                out_name='map-0-filtered'
            )
        else:
            # get kmers occurring from 0-minmap times
            self.call_complex(
                dbdict=dbdict, 
                oper="union",
                min_depth=0,
                max_depth=self.params['min_map'][0],
                out_name='map-0-filter-min'
            )
            # get kmers occurring from maxmap->infinity times
            self.call_complex(
                dbdict=dbdict, 
                oper="union",
                min_depth=self.params['max_map'][0],
                max_depth=100000000,
                out_name='map-0-filter-max'
            )
            # get kmers in the union of the min and max filters
            dbdict = {
                0: f"{self.prefix}_map-0-filter-min",
                1: f"{self.prefix}_map-0-filter-max",
            }
            self.call_complex(
                dbdict=dbdict,
                oper="union",
                min_depth=0,
                max_depth=2,
                out_name='map-0-filtered'
            )
            # os.remove(f"{self.prefix}_map-0-filter-min")
            # os.remove(f"{self.prefix}_map-0-filter-max")


    def get_group1_filtered_set(self):
        """
        Prepare database of kmers that did NOT pass the filters
        for group 1. This usually has a mid min_map and a high max_map.
        """
        # get the group1 samples as a numbered dict
        dbdict = {
            idx: f"{self.prefix}_{sname}_count1"
            for idx, sname in enumerate(self.traits_to_samples[1])
        }

        # if EVERY kmer not at 100% in this set should be filtered 
        # then use the union sum, this is faster than below.
        if self.params['max_map'][1] == 1.0:
            self.call_complex(
                dbdict=dbdict,
                oper="union",
                min_depth=0,
                max_depth=max(1, self.params['min_map'][1] - 1),
                out_name='map-1-filtered'
            )
        else:
            # get kmers occurring from 0-minmap times
            self.call_complex(
                dbdict=dbdict, 
                oper="union",
                min_depth=0,
                max_depth=self.params['min_map'][1],
                out_name='map-1-filter-min'
            )
            # get kmers occurring from maxmap->infinity times
            self.call_complex(
                dbdict=dbdict, 
                oper="union",
                min_depth=self.params['max_map'][1],
                max_depth=100000000,
                out_name='map-1-filter-max'
            )
            # get kmers in the union of the min and max filters
            dbdict = {
                0: f"{self.prefix}_map-1-filter-min",
                1: f"{self.prefix}_map-1-filter-max",
            }
            self.call_complex(
                dbdict=dbdict,
                oper="union",
                min_depth=0,
                max_depth=2,
                out_name='map-1-filtered'
            )
            # os.remove(f"{self.prefix}_map-0-filter-min")
            # os.remove(f"{self.prefix}_map-0-filter-max")


    def get_filtered_union(self):
        """
        Get the union of FILTERED kmers sets
        """
        dbdict = {
            0: f"{self.prefix}_map-0-filtered",
            1: f"{self.prefix}_map-1-filtered",
            2: f"{self.prefix}_min_cov-filtered",            
        }        
        self.call_complex(
            dbdict=dbdict,
            oper="union",
            min_depth=1,
            max_depth=3,
            out_name='filtered'
        )


    def get_kmers_passed_filters(self):
        """
        Get the set of kmers that are in total union of kmers but 
        not in the union of filtered kmers.
        """
        cmd = [
            KMTBIN, "-hp", "-t8",
            "simple",
            f"{self.prefix}_union",
            f"{self.prefix}_filtered",
            "kmers_subtract",
            f"{self.prefix}_passed",
            "-ci1"
        ]
        logger.debug(" ".join(cmd))            
        subprocess.run(cmd, check=True)


    def check_overwrite(self):
        """
        Warn user of overwriting.
        """
        if self.project['kfilter']:
            logger.error(
                "\nKfilter results exist, use force to overwrite, or consider "
                "using branching to produce new results on a separate named "
                "branch without overwriting previous results."
            )
            raise KmerkitError("Preventing data overwrite")        


    def run(self, force=False):
        """
        Create a 'complex' input file to run `kmc_tools complex ...`
        """
        # check for current step
        if not force:
            self.check_overwrite()

        # MINCANON FILTER ---------------------------------------------
        # get kmers passing the mincov_canon filter (canon-filtered)
        # self.filter_canon()
        
        self.get_all_single_counts()
        self.get_union_single_counts()
        self.get_min_cov_filtered_set()
        self.get_group0_filtered_set()
        self.get_group1_filtered_set()
        self.get_filtered_union()
        self.get_kmers_passed_filters()

        pre = f"{self.prefix}_"
        stats = KfilterData(
            kmers_total=info(pre + "union"),
            kmers_passed_total=info(pre + "passed"),
            kmers_filtered_total=info(pre + "filtered"),
            kmers_filtered_by_min_cov=info(pre + "min_cov-filtered"),
            kmers_filtered_by_group0_map=info(pre + "map-0-filtered"),
            kmers_filtered_by_group1_map=info(pre + "map-1-filtered"),
            # kmers_filtered_by_canonized=info(),
            database_filtered=pre + "filtered",
            database_passed=pre + "passed",
        )
        logger.info(stats.json(indent=4))
        data = stats.dict()
        check = data['kmers_total'] - data['kmers_filtered_total']
        if not check == data['kmers_passed_total']:
            raise KmerkitError(
                "error in kmer comparisons: total - filtered != passed")

        # cleanup tmp files
        fnames = [
            "map-0-filtered", "map-0-filter-min", "map-0-filter-max",
            "map-1-filtered", "map-1-filter-min", "map-1-filter-max",
            "union", "min_cov",
        ]
        for fname in fnames:
            fname = self.prefix + fname
            if os.path.exists(fname):
                os.remove(fname)

        # save to JSON
        # limit = ["name", "workdir", "versions", "kinit", "kcount"]
        self.project['kfilter'] = KfilterBase(
            params=self._params,
            data=stats,
        )
        # logger.info(Project(**self.project).json(indent=4))
        with open(self.json_file, 'w') as out:
            out.write(Project(**self.project).json(indent=4))



# def filter_canon(self):
#     """
#     Filter kmers that tend to only occur in one form or the
#     other, likely due to adapters. 

#     Returns:
#     None. 
#     Writes kmc database prefix=<workdir>/{name}_kfilter_canon-filtered
#     """
#     # database file for storing counts of kmers being present in both forms
#     bothdb = self.prefix + "_canonsums"

#     # iterate over samples
#     for sidx, sname in enumerate(self.samples):

#         # get fastq dict w/ fastqs used in kcounts
#         fastq_dict = {sname: self.kcountdf.at[sname, "fastqs"].split(",")}

#         # get options used in kcounts, but not mindepth (must be 1)
#         kcount_params = {
#             'workdir': self.workdir,
#             'name': 'tmp-canon',
#             'fastq_dict': fastq_dict,
#             'kmersize': self.kcountdf.at[sname, "kmersize"],
#             'subsample_reads': self.kcountdf.at[sname, "subsample_reads"],
#             'trim_reads': self.kcountdf.at[sname, "trimmed"],
#             'mindepth': 1, 
#             'maxdepth': self.kcountdf.at[sname, "maxdepth"],
#             'maxcount': self.kcountdf.at[sname, "maxcount"],
#             'canonical': 0,
#         }

#         # count non-canon kmers in sample with same kcount args except -ci1
#         logger.info(f"counting non-canonical kmers [{sname}]")
#         kco = Kcount(**kcount_params)
#         kco.run()
#         nonc_db = os.path.join(self.workdir, f"kcount_tmp-canon_{sname}")

#         # subtract counts of non-canon from canon, only keep if -ci1, 
#         # these are the kmers observed occurring both ways in this sample.
#         canon_db = self.kcountdf.at[sname, "database"]
#         cmd = [
#             KMTBIN,
#             "-hp",                
#             "simple",
#             canon_db, 
#             nonc_db, 
#             "counters_subtract",
#             self.prefix + "_tmp1",  # tmp1=kmers occurring both ways
#             "-ci1"
#         ]
#         logger.debug(" ".join(cmd))            
#         subprocess.run(cmd, check=True)

#         # set count =1 on a tmp copy of bothways kmers db
#         cmd = [
#             KMTBIN,
#             "-hp",                
#             "transform",
#             self.prefix + "_tmp1",
#             "set_counts", 
#             "1",
#             self.prefix + "_tmp2",  # tmp2=kmers w/ count=1
#         ]
#         logger.debug(" ".join(cmd))
#         subprocess.run(cmd, check=True)

#         # if the first sample then simply copy tmp2 to BOTHDB
#         if not sidx:
#             for suff in [".kmc_suf", ".kmc_pre"]:                
#                 shutil.copyfile(
#                     self.prefix + "_tmp2" + suff,
#                     bothdb + suff,
#                 )

#         # if not, then make copy of BOTHDB, and fill BOTHDB using kmctools
#         else:
#             for suff in [".kmc_suf", ".kmc_pre"]:                                
#                 shutil.copyfile(
#                     bothdb + suff,
#                     self.prefix + "_tmp3" + suff,
#                 )

#             # fill bothdb as the sum union of tmp2 and tmp3
#             cmd = [
#                 KMTBIN,
#                 "-hp",
#                 "simple",
#                 self.prefix + "_tmp2",   # kmers w/ count=1 from this samp
#                 self.prefix + "_tmp3",   # kmers w/ count=sum so far
#                 "union",
#                 bothdb,
#                 "-ci1",
#                 "-cs{}".format(len(self.samples) + 1),
#             ]
#             logger.debug(" ".join(cmd))
#             subprocess.run(cmd, check=True)              
            
#         # cleanup
#         logger.debug(f"removing database tmp-canon_{sname}")

#     # get stats on full canonized set
#     sumkmers = info(self.prefix + "_canonsums")

#     # dump and calculate stats on canonsums counts
#     dump(self.prefix + "_canonsums", write_kmers=False)
#     with open(self.prefix + "_canonsums_kmers.txt", 'r') as indat:
#         counts = np.loadtxt(indat, dtype=np.uint16)
#         counts = counts / len(self.samples)
#         mcanon = counts.mean()
#         scanon = counts.std()
#         del counts

#     # report statistics on canonization
#     logger.info(
#         "kmers (proportion) in canonized form: "
#         f"mean={mcanon:.2f}; std={scanon:.2f}"
#     )

#     # report a warning if too few are canonized:
#     if mcanon < 0.1:
#         logger.warning(
#             "Very few kmers are canonized in this dataset. This may "
#             "indicate your data is strand-specific (e.g., ddRAD) in "
#             "which case you should set mincov_canon=0.0. Alternatively, "
#             "if WGS reads, then your data may be very low coverage."
#         )            

#     # filter BOTHDB to get kmers in mincanon prop. of samples (CANON)
#     mincanon_asint = int(np.floor(self.mincov_canon * len(self.samples)))
#     cmd = [
#         KMTBIN, 
#         "transform",
#         bothdb,
#         "reduce",
#         self.prefix + "_mincanon-filter", 
#         "-ci{}".format(mincanon_asint),
#         "-cx1000000000",
#         # "-cs255",        # 65K
#     ]
#     logger.debug(" ".join(cmd))
#     subprocess.run(cmd, check=True)              

#     # calculate and report statistics on filtered canonized set.
#     sumfiltkmers = info(self.prefix + "_mincanon-filter")
#     logger.info(f"kmers filtered by mincov_canon: {sumkmers - sumfiltkmers}")

#     # clean up tmp files
#     countpre = os.path.join(self.workdir, "kcount")
#     os.remove(countpre + "_tmp-canon.csv")
#     for sname in self.samples:
#         os.remove(countpre + "_tmp-canon_" + sname + ".kmc_pre")
#         os.remove(countpre + "_tmp-canon_" + sname + ".kmc_suf")            
#     os.remove(self.prefix + "_canonsums_kmers.txt")
#     for suffix in ["tmp1", "tmp2", "tmp3", "canonsums"]:
#         os.remove(self.prefix + "_" + suffix + ".kmc_pre")
#         os.remove(self.prefix + "_" + suffix + ".kmc_suf")            



if __name__ == "__main__":

    # first run: python3 kcount.py 
    import kmerkit
    kmerkit.set_loglevel("DEBUG")

    # fake data
    PHENOS = "~/Documents/kmerkit/data/amaranths-phenos.csv"
    traits_dict = kmerkit.utils.get_traits_dict_from_csv(PHENOS)

    # load database with phenotypes data
    kgp = Kfilter(
        json_file="/tmp/test.json",
        traits=traits_dict,
        min_cov=0.5,
        min_map={
            0: 0.0,
            1: 0.5,
        },
        max_map={
            0: 0.0,
            1: 1.0,
        },
        min_map_canon={
            0: 0.0,
            1: 0.5,
        },
    )

    # dump the kmers to a file
    kgp.run()
