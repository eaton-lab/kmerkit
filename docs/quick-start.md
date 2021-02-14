


## Tutorial
This tutorial will walk through an example analysis with `kmerkit`, 
demonstrating just one of many possible operations that can be performed
through the combination of its tool. See the [cookbook](#cookcook)
section for many additional examples that demonstrate analyses on real 
data.

Once you have `kmerkit` installed using `conda` you will be able to call 
the top-level `kmerkit` binary. This gives you access to each of the 
underlying tools (commands). You will almost always call one of these 
subcommands using the syntax `kmerkit {command} {-options ...}`. In
addition to the commands listed below, there are two useful options only
available to the `kmerkit` binary: `kmerkit docs`, which will open the docs
in your default browser (locally); and `kmerkit --install-completion` which
will add hooks to your shell settings to allow auto-completing commands 
(can save a few keystrokes).

```bash
kmerkit --help
```

??? example "stdout"

	```bash
	Usage: kmerkit [OPTIONS] COMMAND [ARGS]...

	  Call kmerkit commands to access tools in the kmerkit toolkit,  and kmerkit
	  COMMAND -h to see help options for each tool (e.g., kmerkit kcount -h)

	Options:
	  --install-completion  Install completion for the current shell.
	  --show-completion     Show completion for the current shell, to copy it or
	                        customize the installation.

	  -h, --help            Show this message and exit.

	Commands:
	  docs     opens the kmerkit documentation in a browser
	  info     Show status or flow diagram of project
	  kcount   Count kmers in fastq/a files using KMC.
	  kfilter  filter kmers based on distribution among samples/traits
	```


### Commands

#### kcount
The `kcount` command is used to generate a kmer database by counting all kmers
from one or more fastq/a files provided as input, and that meet the designated
thresholds (e.g., mindepth). Additional filtering can be applied at later steps,
but this establishes the data that you will use in downstream analyses. For 
each sample a database is composed of two files named 
`<workdir>/<name>-<sample-name>.pre` and `<workdir>/<name>-<sample-name>.suf`

```bash
kmerkit kcount --help
```

#### kfilter
The `kfilter` command is one of two options for filtering kmers based on a
comparative analysis among groups of samples. You tell it how to group 
samples, and define filters to generate a new kmer database composed only of
kmers that pass all filters. An arbitrary number of groups/filters can be
defined.

```bash
kmerkit kfilter --help
```
...
