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

from enum import Enum
import tempfile
import typer
from kmerkit import __version__
from .utils import set_loglevel
from .kcount import Kcount, get_fastq_dict_from_path

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

@app.callback()
def callback():
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
def docs():
    "opens the kmerkit documentation in a browser"
    typer.echo("Opening kmerkits documentation in browser")
    typer.launch("https://kmerkit.readthedocs.io")

@app.command()
def info(
    name: str = typer.Option(..., "--name", "-n", help="Prefix for output files."),
    workdir: str = typer.Option(..., "--workdir", "-w", help="Dir name for output files."),
    flow: bool = typer.Option(False, help="show status as flow diagram"),
    stats: bool = typer.Option(False, help="show per-sample stats"),
    ):
    """
    Show status or flow diagram of project
    """
    # print the stats (json?) file location
    typer.secho(
        f"project path: {workdir}/{name}", fg=typer.colors.CYAN)

    # TODO: create a class for stats and progress all in JSON

    # shows a flowchart based on progress
    if flow:
        pass

    if stats:
        pass


@app.command(context_settings=CONTEXT_SETTINGS)
def kcount(
    name: str = typer.Option(
        "test",
        "-n", "--name",
        help="Prefix for output files.",
        ),
    workdir: str = typer.Option(
        tempfile.gettempdir(), 
        "-w", "--workdir",
        help="Dir name for output files.",
        ),
    fastqs: str = typer.Option(
        ..., 
        help="Path to fastq file or dir of fastqs files (e.g., ./fastqs/*.gz)",
        ),
    # fastqs: Path = typer.Option(
    #     ..., 
    #     help="Path to fastq file or dir of fastqs files (e.g., ./fastqs/*.gz)",
    #     show_default=False,
    #     exists=True,
    #     dir_okay=True,
    #     file_okay=True,
    #     resolve_path=True,
    #     allow_dash=True,
    #     ),
    trim_reads: bool = typer.Option(True),
    canonical: bool = True,
    kmersize: int = typer.Option(17, min=2),
    mindepth: int = typer.Option(1, min=0),
    maxdepth: int = typer.Option(int(1e9), min=1),
    maxcount: int = typer.Option(255, min=1),
    subsample_reads: int = typer.Option(None),
    loglevel: LogLevel = LogLevel.INFO,
    ):
    """
    Count kmers in fastq/a files using KMC. 

    kcount will create the workdir if it does not yet exist, and will
    write kmer database files for each sample in files named
    <workdir>/kcount_<name>.kmc_[suf,pre]

    Example call:
      kmerkit kcount -n test -w /tmp --fastqs ./data/*.gz  
    """
    # report the module
    typer.secho(
        "kcount: counting kmers from fastq/a files using KMC",
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

    # TODO: print CLI string to the logger
    # ...

    # parse the fastq_dict from string
    fastq_dict = get_fastq_dict_from_path(fastqs, "_R")

    # run the command
    counter = Kcount(
        name=name,
        workdir=workdir,
        fastq_dict=fastq_dict,
        kmersize=kmersize,
        trim_reads=trim_reads,
        subsample_reads=subsample_reads,
        mindepth=mindepth,
        maxdepth=maxdepth,
        maxcount=maxcount,
        canonical=canonical,
    )
    # print(counter.statsdf.T)
    counter.run()



@app.command()
def kfilter(name: str, workdir: str, traits: str):
    """
    filter kmers based on distribution among samples/traits
    """
    typer.echo("kfilter test")
