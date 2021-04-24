#!/usr/bin/env python

"""
Status object keeps track of sample stats and project status
in JSON format.
"""

import json


class Status:
    """

    """
    def __init__(self):
        self.samples = {}
        self.status = {}

    def __repr__(self):
        return f"<Status {self.samples}>"

    def flowchart(self):
        r"""       
          kcount ----> kfilter -----> kmatrix --> kgwas
                 \             \ /
                  -->   ktree   -
                                 \
                                  --> kextract --> kassemble
        """


if __name__ == "__main__":

    status = Status()
    print(status)
