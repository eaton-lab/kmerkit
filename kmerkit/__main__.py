#!/usr/bin/env python

"""
TODO:
We could develop command line tools that combine multiple class calls together.

kcount --name hybridus \
       --workdir /tmp \
       --fastqs ~/Documents/ipyrad/isolation/reftest_fastqs/[1-2]*_0_R*_.fastq.gz \
       --trim \
       --canonical \
       --kmersize 31 \
	   --mindepth 1 \


ktree --name hybridus \
	  --workdir /tmp \
      --tree hybridus-tree.nwk \
      --model ... \
      --threshold ... \
      --target-options ... \
"""
import os
from enum import Enum
import tempfile
from pathlib import Path
from typing import List, Tuple, Optional
import typer
import pandas as pd
from kmerkit import __version__
from kmerkit.kschema import Project
from kmerkit.utils import set_loglevel, KmerkitError
from kmerkit.utils import get_fastq_dict_from_path, get_traits_dict_from_csv
from kmerkit.kinit import init_project
from kmerkit.ktrim import Ktrim
from kmerkit.kcount import Kcount
from kmerkit.kfilter import Kfilter
from kmerkit.kextract import Kextract

# add the -h option for showing help
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
 

class LogLevel(str, Enum):
    "categorical options for loglevel to CLI"
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"

# creates the top-level kmerkit app
app = typer.Typer(add_completion=True, context_settings=CONTEXT_SETTINGS)


def version_callback(value: bool):
    "Part 1 of overly complicated typer --version option"
    if value:
        typer.echo(f"kmerkit {__version__}")
        raise typer.Exit()

def docs_callback(value: bool):
    "function to open docs"
    if value:
        typer.echo("Opening https://eaton-lab.org/kmerkit in default browser")
        typer.launch("https://eaton-lab.org/kmerkit")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(None, "--version", "-v", callback=version_callback, is_eager=True, help="print version and exit."),
    docs: bool = typer.Option(None, "--docs", callback=docs_callback, is_eager=True, help="Open documentation in browser."),
    ):
    """
    Call kmerkit commands to access tools in the kmerkit toolkit,
    and kmerkit COMMAND -h to see help options for each tool
    (e.g., kmerkit kcount -h)
    """
    typer.secho(
        f"kmerkit (v.{__version__}): the kmer operations toolkit",
        fg=typer.colors.MAGENTA, bold=True,
    )


@app.command()
def stats(
    json_file: str = typer.Option(..., "--json", "-j", help="kmerkit project JSON file"),
    module: str = typer.Argument(None, help="kmerkit module") 
    # flow: bool = typer.Option(False, help="show status as flow diagram"),
    ):
    """
    Return tabular stats output for a project module to stdout.
    """
    module = module.lower()

    # load the project
    project = Project.parse_file(json_file).dict()
    if module not in project:
        typer.Abort(f"no results exist for {module}")

    typer.secho(
        "{} stats for project {}".format(
            module,
            os.path.basename(json_file).rsplit(".")[0],
        ), 
        fg=typer.colors.CYAN,
    )
    # TODO, make a class to normalize each type to CSV.
    print(pd.json_normalize(project[module]['data']))



@app.command(context_settings=CONTEXT_SETTINGS)
def init(
    name: str = typer.Option("test", "-n", "--name", help="Project name prefix"),
    workdir: str = typer.Option(tempfile.gettempdir(), "-w", "--workdir", help="Project directory"),
    delim: str = typer.Option("_R", help="sample name delimiter"),
    loglevel: LogLevel = typer.Option(LogLevel.INFO, help="logging level"),
    force: bool = typer.Option(False, help="overwrite existing"),
    data: List[Path] = typer.Argument(...,
        show_default=False,
        exists=True,          # <- prob should be false for moving json files.
        dir_okay=True,
        file_okay=True,
        resolve_path=False,
        allow_dash=True,
        help=("File path(s) to input fastq/a data files")
        ),
    # trim_reads: bool = typer.Option(True),
    # subsample_reads: int = typer.Option(None),
    ):
    """
    Initialize a kmerkit project from fastq/a input files.

    Creates a JSON project file in <workdir>/<name>.json. Sample
    names are parsed from input filenames by splitting on the last
    occurrence of the optional 'delim' character (default is '_R').
    Paired reads are automatically detected from _R1 and _R2 in names.
    Examples:

    kmerkit init -n test -w /tmp ./data/fastqs/*.gz\n
    kmerkit init -n test -w /tmp ./data-1/A.fastq ./data-2/B.fastq
    """
    # parse the fastq_dict from string
    set_loglevel(loglevel)
    fastq_dict = get_fastq_dict_from_path(None, data, delim)

    try:
        init_project(name=name, workdir=workdir, fastq_dict=fastq_dict, force=force)
    except KmerkitError:
        typer.Abort()



@app.command(context_settings=CONTEXT_SETTINGS)
def count(
    json_file: Path = typer.Option(..., "-j", "--json"),
    kmer_size: int = typer.Option(17, min=2),
    min_depth: int = typer.Option(1, min=1),
    max_depth: int = typer.Option(int(1e9), min=1),
    max_count: int = typer.Option(65535, min=255),
    canonical: bool = typer.Option(True),
    workers: int = typer.Option(1, help="N worker processes"),
    threads: int = typer.Option(None, help="N threads per worker"),
    force: bool = typer.Option(False, help="overwrite existing"),
    max_ram_per_worker: int = typer.Option(12, help="max RAM in Gb"),
    loglevel: LogLevel = LogLevel.INFO,
    ):
    """
    Count kmers in fastq/a files using KMC. 

    kcount will write kmer database files for each sample to 
    <workdir>/<name>_kcount_.kmc_[suf,pre]. Example:

    kmerkit count -j test.json --kmer-size 35 --min-depth 5
    """
    # report the module
    typer.secho(
        "count: counting kmers from fastq/a files using KMC",
        fg=typer.colors.MAGENTA,
        bold=False,
    )

    # set the loglevel
    set_loglevel(loglevel)
    typer.secho(
        f"loglevel: {loglevel}, logfile: STDERR",
        fg=typer.colors.MAGENTA,
        bold=False,
    )

    # run the command
    counter = Kcount(
        str(json_file),
        kmer_size=kmer_size,
        min_depth=min_depth,
        max_depth=max_depth,
        max_count=max_count,
        canonical=canonical,
    )
    # print(counter.statsdf.T)
    try:
        counter.run(
            threads=threads, 
            workers=workers, 
            force=force, 
            max_ram=max_ram_per_worker,
        )
    except KmerkitError as exc:
        typer.Abort(exc)


# @app.command()
# def group(
#     json_file: Path = typer.Option(..., '-j', '--json'),
#     group: List[List[str]] = typer.Option(None),
#     ):
#     print(group)



@app.command()
def filter(
    json_file: Path = typer.Option(..., '-j', '--json'),
    group0: Optional[List[str]] = typer.Option(None, '--group0', '-0'),
    group1: Optional[List[str]] = typer.Option(None, '--group1', '-1'),
    traits: Optional[Path] = typer.Option(None, '--traits-file'),
    min_cov: float = typer.Option(0.0),
    min_map: Tuple[float,float] = typer.Option((0.0, 0.1)),
    max_map: Tuple[float,float] = typer.Option((0.1, 1.0)),
    loglevel: LogLevel = LogLevel.INFO,
    force: bool = typer.Option(False, help="overwrite existing"),
    # min_map_canon
    ):
    """
    Filter kmers based on frequencies among grouped samples.
    """
    # report the module
    typer.secho(
        "filter: filter kmers based on frequency in case/control groups",
        fg=typer.colors.MAGENTA,
        bold=False,
    )
    # set the loglevel
    set_loglevel(loglevel)
    typer.secho(
        f"loglevel: {loglevel}, logfile: STDERR",
        fg=typer.colors.MAGENTA,
        bold=False,
    )

    # fake data
    if traits:
        traits_dict = get_traits_dict_from_csv(traits)
    else:
        traits_dict = {0: [], 1: []}
    traits_dict[0].extend(group0)
    traits_dict[1].extend(group1)

    # load database with phenotypes data
    kgp = Kfilter(
        json_file=json_file,
        traits=traits_dict,
        min_cov=min_cov,
        min_map={0: min_map[0], 1: min_map[1]},
        max_map={0: max_map[0], 1: max_map[1]},        
        min_map_canon={0: 0.0, 1: 0.5},
    )
    kgp.run(force=force)


@app.command()
def extract(
    json_file: Path = typer.Option(..., '-j', '--json'),
    min_kmers_per_read: int = typer.Option(1),
    keep_paired: bool = typer.Option(True),
    loglevel: LogLevel = LogLevel.INFO,
    force: bool = typer.Option(False, help="overwrite existing"),  
    samples: Optional[List[str]] = typer.Argument(None),
    ):
    """
    Extract reads from fastq/a files containing target kmers.

    Reads must contain at least 'min-kmers-per-read' kmers in them.
    If 'keep-paired' then reads are returns as paired-end. Samples can 
    be entered as arguments in three possible ways: (1) enter sample 
    names that are in the init database; (2) enter an integer for 
    group0 or group1 from the kfilter database; (3) enter a file path 
    to one or more fastq files.

    kmerkit extract -j test.json A B C D      # select from init\n
    kmerkit extract -j test.json 1            # select from filter group\n
    kmerkit extract -j test.json ./data/*.gz  # select new files\n
    """
    typer.secho(
        "extract: extract reads containing target kmers",
        fg=typer.colors.MAGENTA,
        bold=False,
    )
    set_loglevel(loglevel)  

    try:
        kex = Kextract(
            json_file=json_file,
            samples=samples,
            min_kmers_per_read=min_kmers_per_read,
            keep_paired=keep_paired,
        )
        kex.run(force=force)
    except KmerkitError:
        typer.Abort()


@app.command()
def trim(
    json_file: Path = typer.Option(..., "-j", "--json"),
    subsample: int = typer.Option(None, help="subsample to N reads"),
    workers: int = typer.Option(None, help="N worker processes"),
    force: bool = typer.Option(False, help="overwrite existing"),
    loglevel: LogLevel = LogLevel.INFO,    
    ):
    """
    Trim, filter, or subsample reads using fastp.

    kmerkit trim -j test.json --subsample 1000000 --cores 20
    """
    typer.secho(
        "trim: trim and filter reads using fastp (default settings)",
        fg=typer.colors.MAGENTA,
        bold=False,
    )
    set_loglevel(loglevel)

    try:
        ktr = Ktrim(json_file=json_file, subsample=subsample)
        ktr.run(force=force, workers=workers)
    except KmerkitError:
        typer.Abort()
