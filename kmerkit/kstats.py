#!/usr/bin/env python

"""
Returns statistics and status for a project JSON file.
"""

import itertools
import typer
import pandas as pd
from loguru import logger
from kmerkit.kschema import Project, Kinit

COLUMNS = {
    "kinit": ["kinit_reads"],
    "ktrim": ["ktrim_reads"],
    "kcount": ["kcount_uniq_kmers"], 
    "kfilter": ["kfilter_group", "kfilter_db_kmers"],
    "kextract": ["kextract_reads"],
}

class Kstats:
    """
    Functions to return stats for different modules from a
    project JSON file.
    """
    def __init__(self, json_file):
        self.project = Project.parse_file(json_file)
        self.proj = self.project.dict()


    def json(self):
        """
        If no module is selected then we simply return the JSON
        """
        typer.secho(
            f"Project JSON data:",
            fg=typer.colors.CYAN,
        )
        print(self.project.json(indent=4, exclude_none=True))


    def summary(self):
        """
        Returns a text description of stats to STDOUT
        """
        typer.secho(
            f"Project: {self.proj['name']}", 
            fg=typer.colors.CYAN,
        )
        modules = [i for i in self.proj if (self.proj[i] and i[0] == 'k')]

        # if only init has been run then show init files and exit
        if modules == ['kinit']:
            typer.secho(
                "No modules run yet.\nJSON:", 
                fg=typer.colors.CYAN,
            )
            self.json()
            return


        print(f"Modules: {', '.join(modules)}")
        # select columns depending on modules
        columns = list(
            itertools.chain(*[COLUMNS[module] for module in modules]))

        # otherwise, data exists, so we build a TSV dataframe
        data = pd.DataFrame(
            columns=columns,            
            index=self.proj['kinit']['data']
        )


        if 'ktrim' in modules:
            data.loc[:, 'ktrim_reads'] = [
                self.proj['ktrim']['data'][sname]['reads_total']
                for sname in data.index
            ]

        data.loc[:, 'kinit_reads'] = [
            self.proj['kcount']['data'][sname]['reads_total']
            for sname in data.index
        ]
        print(data.to_string())
        return


    def run(self, module=None):
        """
        Selects a function to run based on module user input
        """
        if module is None:
            self.json()



if __name__ == "__main__":

    status = Kstats("/tmp/null.json")
    status.run()
