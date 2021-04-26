#!/usr/bin/env python

"""
kdump
"""

import os
import re
from loguru import logger
from kmerkit.kmctools import dump
from kmerkit.kschema import Project


class Kdump:
    """
    ...
    """
    def __init__(
        self, 
        json_file, 
        samples, 
        min_depth, 
        max_depth, 
        write_kmers, 
        write_counts,
        ):

        # load user inputs
        self.json_file = json_file
        self.project = Project.parse_file(json_file).dict()
        self.kc_database = self.project['kcount']['data']
        self.kf_database = self.project['kfilter']['data']

        # store params
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.write_kmers = write_kmers
        self.write_counts = write_counts

        # updates samples from a list to a dict mapping {name: db}
        self.samples = samples
        self.select_samples()


    def select_samples(self):
        """
        Select samples from kcount name database or from a filepath.
        """
        # names in count db: {'aaa', 'bbb', 'ccc', 'ccd', 'cce'}
        setkc = set(self.kc_database)

        # names in traits: {'aaa', 'cc[c,d,e]', '/tmp/name_kfilter_passed'}
        sets = set(self.samples)

        # try to regex expand non-matching string names
        set_keys = sorted(sets)
        for sname in set_keys:
            if sname not in setkc:
                search = (re.match(sname, i) for i in setkc)
                matches = [i.string for i in search if i]
                if matches:
                    logger.debug(f"{sname} expanded to {matches}")
                    sets.update(set(matches))
                    sets.remove(sname)
                else:
                    if not os.path.isfile(sname + ".kmc_suf"):
                        sets.remove(sname)

        self.samples = {}
        for i in sets:
            if i in self.kc_database:
                self.samples[i] = self.kc_database[i]['database']
            else:
                sname = os.path.basename(i).rsplit(".")[0]
                self.samples[sname] = i


    def run(self):
        """
        ...
        """
        logger.debug(self.samples)
        for sname in self.samples:
            database = self.samples[sname]
            args = (
                database, 
                self.write_kmers, self.write_counts, 
                self.min_depth, self.max_depth, 
            )
            dump(*args)
