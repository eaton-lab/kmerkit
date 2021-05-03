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
    --min_canon 0.5 \
    --min_cov 5
"""

# TODO CANON FILTER:
# - applies at end only to the kmers that passed other filters
# - count other-canon kmers in group 1 
# - get counter_subtract: canon-kmers - non-canon-kmers
# - other-canon kmer sets to count=1
# - get canon-kmers with freq > X in group 1
# - get intersection of passed and min-canon-passed


import os
import re
import itertools
import subprocess
import concurrent.futures
from math import floor

from loguru import logger
from kmerkit.kmctools import KMTBIN, info
from kmerkit.utils import Group, COMPLEX, KmerkitError
from kmerkit.kschema import KfilterParams, KfilterData, KfilterBase, Project
# from kmerkit.parallel import Cluster

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
    def __init__(self, json_file, traits_dict, min_cov, min_map, max_map, min_map_canon):

        # load user inputs
        self.json_file = json_file
        self.project = Project.parse_file(json_file).dict()
        self.traits_to_samples = traits_dict

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

        # keep track of database
        self.dbs = {}


    def select_samples(self):
        """
        Filter samples to those present in both database and traits_dict
        """
        # if any traits_dict samples can be regex expanded
        # self.traits_to_samples = None

        # names in count db: {'aaa', 'bbb', 'ccc', 'ccd', 'cce'}
        setp = set(self.database)

        # names in traits: {'aaa', 'bbb', 'cc[c,d,e]'}
        sets = set(itertools.chain(*self.traits_to_samples.values()))

        # try to regex expand non-matching string names
        for key in self.traits_to_samples:
            snames = self.traits_to_samples[key]
            for sname in snames:
                if sname not in setp:
                    search = (re.match(sname, i) for i in setp)
                    matches = [i.string for i in search if i]
                    if matches:
                        logger.debug(f"{sname} expanded to {matches}")
                        self.traits_to_samples[key].extend(matches)

        # names in traits: {'aaa', 'bbb', 'ccc', 'ccd', 'ccf'}
        sets = set(itertools.chain(*self.traits_to_samples.values()))

        # build reverse dict
        revtraits = {}
        for key in self.traits_to_samples:
            for sample in self.traits_to_samples[key]:
                revtraits[sample] = key

        # raise an error if no names overlap
        if setp.isdisjoint(sets):
            msg = (
                "Sample names entered do not match names in database:\n"
                f"  KMER_DATABASE: {', '.join(setp)}\n"
                f"  TRAITS:   {', '.join(sets)}"
            )
            logger.error(msg)
            raise KmerkitError(msg)

        # warning for samples only in traits that are not in database
        only_traits = setp.difference(sets)
        if only_traits:
            logger.warning(
                "Skipping samples in groups or kcount database but not both:"
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
            self.params['min_cov'] = int(max(1, 
                floor(self.params['min_cov'] * len(self.database)))
            )
        logger.info(
            f"Total: nsamples={len(self.database):2d}, "
            f"target kmer occurrences at min-cov>={self.params['min_cov']}"
        )            

        # int encoded map values
        for key in [0, 1]:

            # get n samples and require samples
            nsamples = len(self.traits_to_samples[key])
            assert nsamples, f"No samples are set to state={key}"

            # set maxmap to int
            self.params['max_map'][key] = int(
                max(self.params['min_map'][key],
                    min(nsamples,
                        int(floor(self.params['max_map'][key] * nsamples)))
            ))

            # in key0 the maxmap is actually inverse (allowed)
            if not key:
                self.params['max_map'][key] = nsamples - self.params['max_map'][key]

            # check lower limit
            if not self.params['max_map'][key]:
                self.params['max_map'][key] = 1

            # set minmap last and enforce min of 1
            self.params['min_map'][key] = int(max(1, 
                (floor(self.params['min_map'][key] * nsamples))
            ))

            # not yet implemented
            self.params['min_map_canon'] = 100            

            # report
            logger.info(
                f"Group {key}: nsamples={nsamples:2d}, "
                f"{'target ' if bool(key) else 'filter '}"
                "kmer occurrences at "
                f"min-cov>={self.params['min_map'][key]:2d}, "
                f"max-cov<={self.params['max_map'][key]:2d}"
            )


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


    def get_all_single_counts(self, workers=4, threads=2):
        """
        PARALLELIZE THIS -- seems to max out at 200%
        Prepare single count database of every sample which will be used
        for presence/absence set arithmetic in filters.
        """
        with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:           
            for sname in self.database:
                cmd = [
                    KMTBIN, "-hp", "-t{}".format(threads),
                    "transform",
                    self.database[sname]['database'],
                    "set_counts", "1",
                    f"{self.prefix}_{sname}_count1",
                ]
                logger.debug(" ".join(cmd))
                pool.submit(subprocess.run, cmd)


    def get_union_with_counts(self):
        """
        Prepare summed count database from all single count databases.
        This is the full set of observed kmers.
        """
        dbdict = {
            idx: f"{self.database[sname]['database']}"
            for idx, sname in enumerate(self.database)
        }
        # get ALL kmers
        self.call_complex(
            dbdict=dbdict,
            oper="union",
            min_depth=1,
            max_depth=1000000000,
            out_name='union_counts',
        )


    def get_min_cov_passed_set(self):
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
            min_depth=self.params['min_cov'],
            max_depth=100000000,
            out_name='min_cov-passed',
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

        # simple process if max-map group0 is ALL
        cond0 = self.params['min_map'][0] == 1
        cond1 = self.params['max_map'][0] == len(self.traits_to_samples[0])
        if cond0 and cond1:
            outname = 'map-0-filtered'
        else:
            outname = 'map-0-union'

        # EVERY kmer in group0
        self.call_complex(
            dbdict=dbdict,
            oper="union",
            min_depth=1,
            max_depth=len(self.traits_to_samples[0]),
            out_name=outname,
        )

        # end early if all are filtered.
        if outname == 'map-0-filtered':
            return

        # kmers ALLOWED in final set
        cmd = [
            KMTBIN, "-hp", "-t8", "transform",
            f"{self.prefix}_map-0-union",
            "reduce",
            f"{self.prefix}_map-0-passed",
            f"-ci{self.params['min_map'][0]}",
            f"-cx{self.params['max_map'][0]}",
        ]
        logger.debug(" ".join(cmd))            
        subprocess.run(cmd, check=True)

        # kmers FILTERED from in final set
        cmd = [
            KMTBIN, "-hp", "-t8", "simple",
            f"{self.prefix}_map-0-union",
            f"{self.prefix}_map-0-passed",
            "kmers_subtract",
            f"{self.prefix}_map-0-filtered",
            "-ci1"
        ]
        logger.debug(" ".join(cmd))            
        subprocess.run(cmd, check=True)


    def get_group1_passed_set(self):
        """
        Prepare database of kmers that did NOT pass the filters
        for group 1. This usually has a mid min_map and a high max_map.
        """
        # get the group1 samples as a numbered dict
        dbdict = {
            idx: f"{self.prefix}_{sname}_count1"
            for idx, sname in enumerate(self.traits_to_samples[1])
        }

        self.call_complex(
            dbdict=dbdict,
            oper="union",
            min_depth=self.params['min_map'][1],
            max_depth=self.params['max_map'][1],
            out_name='map-1-passed'
        )


    def get_passed_intersect(self):
        """
        Get the union of FILTERED kmers sets
        """
        dbdict = {
            0: f"{self.prefix}_min_cov-passed",
            1: f"{self.prefix}_map-1-passed",
        }
        self.call_complex(
            dbdict=dbdict,
            oper="intersect",
            min_depth=1,
            max_depth=10000000,
            out_name='passed_intersect'
        )


    def get_kmers_passed(self):
        """
        Get the set of kmers that are in total union of kmers but 
        not in the union of filtered kmers.
        """
        cmd = [
            KMTBIN, "-hp", "-t8", "simple",
            f"{self.prefix}_passed_intersect",
            f"{self.prefix}_map-0-filtered",
            "kmers_subtract",
            f"{self.prefix}_passed",
            "-ci1"
        ]
        logger.debug(" ".join(cmd))            
        subprocess.run(cmd, check=True)


    def get_kmers_passed_counts(self):
        """
        Get the set of kmers that are in total union of kmers but 
        not in the union of filtered kmers.
        """
        cmd = [
            KMTBIN, "-hp", "-t8", "simple",
            f"{self.prefix}_union_counts",
            f"{self.prefix}_passed",
            "intersect",
            f"{self.prefix}_passed_counts",
            "-ci1", "-ocleft"
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

        # distributed on four workers
        self.get_all_single_counts()

        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
            pool.submit(self.get_union_with_counts)
            pool.submit(self.get_min_cov_passed_set)

        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:            
            pool.submit(self.get_group0_filtered_set)
            pool.submit(self.get_group1_passed_set)

        self.get_passed_intersect()
        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
            pool.submit(self.get_kmers_passed)
            pool.submit(self.get_kmers_passed_counts)

        pre = f"{self.prefix}_"
        kmers_total = info(pre + "union_counts")
        kmers_passed = info(pre + "passed")
        kmers_min_cov_passed = info(pre + "min_cov-passed")
        stats = KfilterData(
            kmers_total=kmers_total,
            kmers_passed_total=kmers_passed,
            kmers_filtered_total=kmers_total - kmers_passed,
            kmers_filtered_by_min_cov=kmers_total - kmers_min_cov_passed,
            kmers_filtered_by_group0_map=info(pre + "map-0-filtered"),
            kmers_filtered_by_group1_map=kmers_total - info(pre + "map-1-passed"),
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
    TRAITS = kmerkit.utils.get_traits_dict_from_csv(PHENOS)

    # load database with phenotypes data
    kgp = Kfilter(
        json_file="/tmp/test.json",
        traits=TRAITS,
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
