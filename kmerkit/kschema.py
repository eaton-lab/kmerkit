#!/usr/bin/env pytyhon

"""
Serialized JSON storage class to keep a record of a project.

PLANNED JSON SCHEMA:
---------------------
{
    name: test,
    workdir: /tmp/test,
    versions: {
        'kmerkit': 0.0.1,
        'kmc': ...,
    },
    kinit: {
        data: {
            A: [a_R1.fastq, a_R2.fastq],
            B: [b_R1.fastq, b_R2.fastq],
            ...,            
        },
        commands: null
    },
    kcount: {
        data: {

        },
        params: {
            ...,
        },
    },
    kfilter: {
        samples: {
            A: [a.pre, a.suf],
        },
        params: {
            ...,
        }
    },
    kmatrix: null,
    klearn: null,
    ...
}
"""

# pylint: disable=no-name-in-module
from enum import IntEnum
from typing import List, Dict, Optional
from pydantic import BaseModel, Field, FilePath, confloat
import sklearn

# ---------------------------------------------------------------------

class Kinit(BaseModel):
    data: Dict[str, List[FilePath]] = Field(default_factory=dict)
    commands: Dict[str, str] = Field(default_factory=dict)

# ---------------------------------------------------------------------
class KtrimParams(BaseModel):
    subsample: int = Field(None)

class KtrimData(BaseModel):
    data_in: List[FilePath] = Field(None)
    data_out: List[FilePath] = Field(None)
    fastp_stats: dict = Field(None)

class KtrimBase(BaseModel):
    data: Dict[str, KtrimData] = Field(None, type=Dict)
    params: KtrimParams = Field(None)

# ---------------------------------------------------------------------
class KcountData(BaseModel):
    reads_total: int = 0
    kmers_total: int = 0
    kmers_unique: int = 0
    kmers_unique_counted: int = 0
    kmers_below_threshold: int = 0
    kmers_above_threshold: int = 0
    database: str = ""


class MaxCountSizes(IntEnum):
    "unsigned integer sizes used for maxcount limit in kmc"
    small = 255
    medium = 65535
    large = 4294967295
    ultra = 18446744073709551615


class KcountParams(BaseModel):
    kmer_size: int = Field(default=35, ge=5)
    min_depth: int = Field(default=1, gt=0)
    max_depth: int = Field(default=int(1e9), ge=1)
    max_count: MaxCountSizes = Field(default=65535)
    canonical: bool = Field(default=True, description="collapse kmers by revcomp")


class KcountBase(BaseModel):
    params: KcountParams = Field(...)
    data: Dict[str, KcountData] = Field(None, type=Dict)

# ---------------------------------------------------------------------

class KfilterData(BaseModel):
    kmers_total: int = 0
    kmers_passed_total: int = 0
    kmers_filtered_total: int = 0
    kmers_filtered_by_min_cov: int = 0
    kmers_filtered_by_group0_map: int = 0
    kmers_filtered_by_group1_map: int = 0
    # kmers_filtered_by_min_map_canon: int = 0
    database_filtered: str = ""
    database_passed: str = ""


class Pool(IntEnum):
    control: 0
    case: 1


class KfilterParams(BaseModel):
    min_cov: confloat(ge=0, le=1)
    min_map: Dict[int, confloat(ge=0, le=1)]
    max_map: Dict[int, confloat(ge=0, le=1)]
    min_map_canon: Dict[int, confloat(ge=0, le=1)]
    trait_0: List[str] = Field(None)
    trait_1: List[str] = Field(None)


class KfilterBase(BaseModel):
    params: KfilterParams = Field(...)
    data: KfilterData = Field(...)

# ---------------------------------------------------------------------

class KmatrixParams(BaseModel):
    counts: bool
    normalize: bool
    condense: bool
    samples: Optional[List[str]]

class KmatrixData(BaseModel):
    matrix: FilePath
    counts: Optional[FilePath]

class KmatrixBase(BaseModel):
    params: KmatrixParams = Field(None)
    data: KmatrixData = Field(None)

# ---------------------------------------------------------------------

class KextractParams(BaseModel):
    min_kmers_per_read: int = 0
    min_depth_kmer: int = 0
    keep_paired: bool

class KextractData(BaseModel):
    data_in: List[FilePath]
    data_out: List[FilePath]
    kmer_matched_reads: int = 0

class KextractBase(BaseModel):
    params: KextractParams = Field(...)
    data: Dict[str, KextractData] = Field(...)

# ---------------------------------------------------------------------

class Versions(dict):
    kmerkit: str = "0.0.0"
    kmc: str = "0.0.0"
    gemma: str = "0.0.0"
    sklearn: str = sklearn.__version__


class Project(BaseModel):
    name: str = "test"
    workdir: str = "/tmp"
    versions: Versions = Field(default_factory=Versions)
    kinit: Kinit = Field(default_factory=Kinit)
    ktrim: KtrimBase = Field(None)
    kcount: KcountBase = Field(None)
    kfilter: KfilterBase = Field(None)
    kmatrix: KmatrixBase = Field(None)
    kextract: KextractBase = Field(None)
    # kgwas:




if __name__ == "__main__":

    # Step 1: user enters a dictionary of fastq data mapping
    # Write JSON top-level to name, workdir
    # Load Kcount from json:
    # 

    # example
    samples = ['a', 'b', 'c']
    # params = KcountParams(max_depth=1, max_count=65535)
    # data = {i: KcountData() for i in samples}
    # kc = KcountBase(params=params, data=data)

    # proj = Project(kcount=kc)

    # params = KfilterParams(
    #     min_cov=0, min_map={0:0, 1:1}, max_map={0:0, 1:1}, min_map_canon={0:0, 1:0}
    # )
    # data = KfilterData()
    # kf = KfilterBase(params=params, data=data)

    # proj = Project(kcount=kc, kfilter=kf, **proj.dict(exclude={'kcount', 'kfilter'}))
    # print(proj.json(indent=2, exclude_none=False))

    # # write to JSON
    # with open("/tmp/test.json", 'w') as out:
    #     out.write(proj.json(indent=2))
    # sample = KcountData().dict()
    # sample['kmers_total'] = 10
    # print(sample)



