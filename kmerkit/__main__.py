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
from loguru import logger
from kmerkit import __version__
from kmerkit.kschema import Project
from kmerkit.utils import set_loglevel, KmerkitError
from kmerkit.utils import get_fastq_dict_from_path, get_traits_dict_from_csv
from kmerkit.kinit import init_project
from kmerkit.ktrim import Ktrim
from kmerkit.kcount import Kcount
from kmerkit.kfilter import Kfilter
from kmerkit.kextract import Kextract
from kmerkit.kdump import Kdump
from kmerkit.kstats import Kstats

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
    (e.g., kmerkit count -h)
    """
    typer.secho(
        f"kmerkit (v.{__version__}): the kmer operations toolkit",
        fg=typer.colors.MAGENTA, bold=True,
    )


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



@app.command()
def stats(
    json_file: str = typer.Option(..., "--json", "-j", help="kmerkit project JSON file"),
    module: str = typer.Argument(None, help="Show stats for a specific module") 
    ):
    """
    Return summarized results from kmerkit modules to STDOUT.

    If no argument is entered for 'module' then the project JSON
    file is printed to STDOUT. If a module name is entered (e.g., 
    'trim' or 'ktrim') then the results of this module will be 
    summarized and formatted into a tabular (TSV) format and returned.
    Module names can be entered with or without the 'k' prefix.
    """
    try:
        kst = Kstats(json_file)
        if not module:
            kst.run()
        else:
            for mod in module:
                kst.run(mod.lower().lstrip('k'))
    except KmerkitError:
        typer.Abort()


@app.command(context_settings=CONTEXT_SETTINGS)
def count(
    json_file: Path = typer.Option(..., "-j", "--json"),
    kmer_size: int = typer.Option(17, "-k", "--kmer-size", min=2),
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



@app.command(name='filter')
def kfilter(
    json_file: Path = typer.Option(..., '-j', '--json'),
    group0: Optional[List[str]] = typer.Option(None, '--group0', '-0'),
    group1: Optional[List[str]] = typer.Option(None, '--group1', '-1'),
    traits_file: Optional[Path] = typer.Option(None, '--traits-file'),
    min_cov: float = typer.Option(0.0),
    min_map: Tuple[float,float] = typer.Option((0.0, 1.0)),
    max_map: Tuple[float,float] = typer.Option((0.0, 1.0)),
    loglevel: LogLevel = LogLevel.INFO,
    force: bool = typer.Option(False, help="overwrite existing"),
    # min_map_canon
    ):
    """
    Filter kmers based on frequencies among grouped samples.

    The filter kmers (group 0) are subtracted from the target kmers
    (group 1) to find the final target kmer set. You must enter sample
    names (or regex patterns) to --group0 and --group1 to assign samples
    to each group. Kmers in each group are filtered by min-map and 
    max-map ranges...
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
    if traits_file:
        traits_dict = get_traits_dict_from_csv(traits_file)
    else:
        traits_dict = {0: [], 1: []}
    traits_dict[0].extend(group0)
    traits_dict[1].extend(group1)

    # load database with phenotypes data
    try:
        kgp = Kfilter(
            json_file=json_file,
            traits_dict=traits_dict,
            min_cov=min_cov,
            min_map={0: min_map[0], 1: min_map[1]},
            max_map={0: max_map[0], 1: max_map[1]},        
            min_map_canon={0: 0.0, 1: 0.5},
        )
        kgp.run(force=force)
    except KmerkitError:
        typer.Abort()


@app.command()
def extract(
    json_file: Path = typer.Option(..., '-j', '--json'),
    min_kmers_per_read: int = typer.Option(1),
    keep_paired: bool = typer.Option(True),
    loglevel: LogLevel = LogLevel.INFO,
    force: bool = typer.Option(False, help="overwrite existing"),  
    workers: int = typer.Option(1, help="N worker processes"),
    threads: int = typer.Option(None, help="N threads per worker"),
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
        kex.run(force=force, workers=workers, threads=threads)
    except KmerkitError:
        typer.Abort()
    except KeyboardInterrupt:
        typer.Abort("interrupted")


@app.command()
def trim(
    json_file: Path = typer.Option(..., "-j", "--json"),
    subsample: float = typer.Option(None, help="subsample to N reads"),
    workers: int = typer.Option(None, help="N worker processes"),
    # threads: int = typer.Option(None, help="N threads per worker"),
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
        ktr.run(force=force, workers=workers) #, threads=threads)
    except KmerkitError:
        typer.Abort()


@app.command(name='dump')
def kdump(
    json_file: Path = typer.Option(..., "-j", "--json"),
    min_depth: int = typer.Option(1, help="filter to >= min-depth"),
    max_depth: int = typer.Option(100000, help="filter to <= max-depth"),
    write_kmers: bool = typer.Option(True),
    write_counts: bool = typer.Option(True),
    loglevel: LogLevel = LogLevel.INFO,
    samples: List[str] = typer.Argument(..., help="One or more sample names in kcount database"),
    ):
    """
    Write kmers and/or counts to a file from a KMC database.

    kmerkit dump -j test.json --min-depth 5 sample1
    """   
    typer.secho(
        "dump: write kmers and/or counts to a file.",
        fg=typer.colors.MAGENTA,
        bold=False,
    )
    set_loglevel(loglevel)

    try:
        Kdump(
            json_file, 
            samples, 
            min_depth, max_depth, 
            write_kmers, write_counts,
        ).run()

    except KmerkitError:
        typer.Abort()


@app.command()
def branch(
    json_file: Path = typer.Option(..., "-j", "--json"),
    force: bool = typer.Option(False, help="force overwrite."),
    loglevel: LogLevel = LogLevel.INFO,
    name: str = typer.Argument(..., help="New name of branched project."),
    ):
    """
    Branch to create a new named project JSON file.

    kmerkit branch -j test.json test2
    """   
    typer.secho(
        "branch: write kmers and/or counts to a file.",
        fg=typer.colors.MAGENTA,
        bold=False,
    )
    set_loglevel(loglevel)

    try:
        name = name.strip(".json")
        proj = Project.parse_file(json_file).dict()
        proj['name'] = name
        new_json_file = os.path.join(
            os.path.dirname(json_file),
            name + ".json"
        )
        if os.path.exists(new_json_file):
            if not force:
                msg = "JSON file already exists. Use force."
                logger.error(msg)
                raise KmerkitError(msg)
        with open(new_json_file, 'w') as out:
            out.write(Project(**proj).json(indent=4))
        logger.info(f"wrote new branched project to {new_json_file}")
    except KmerkitError:
        typer.Abort()
