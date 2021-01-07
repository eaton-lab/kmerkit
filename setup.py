#!/usr/bin/env python

"""
Setup file for installing kmerkit. 
This is used by conda for users and by pip installation for developers.

cd kmerkit
pip install -e .
"""

import re
from setuptools import setup, find_packages


# get version from __init__.py
INITFILE = "kmerkit/__init__.py"
CUR_VERSION = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                    open(INITFILE, "r").read(),
                    re.M).group(1)


# run setup
setup(
    name="kmerkit",
    version=CUR_VERSION,
    url="https://github.com/eaton-lab/kmerkit",
    author="Deren Eaton and Jasmina Dzurlic",
    author_email="de2356@columbia.edu",
    description="Kmer-related operations (kmc) wrapped in Python",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=[
        "toyplot",
        "numpy",
        "requests",
        "future",
        "loguru", 
        "pandas",
    ],
    entry_points={},
    license='GPL',
    classifiers=[
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',                
    ],
)
