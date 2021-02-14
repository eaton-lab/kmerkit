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


kfilter --name hybridus \
		--workdir /tmp
		--

"""

import typer

app = typer.Typer()



@app.command()
def kcount(name: str, workdir: str, kmersize: int):
    """
    Set up a kcount subcommand: 
        kmerkit kcount --name test --workdir /tmp --kmersize ...
    """
    pass


@app.command()
def kfilter(name: str, workdir: str, traits: str):
    """
    Set up a kfilter subcommand: 
        kmerkit kfilter --name test --workdir /tmp --traits ...
    """
    pass


if __name__ == '__main__':
    app()
