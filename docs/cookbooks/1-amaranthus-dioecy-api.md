

### Workflow diagram

```mermaid
%%{init: {'theme': 'dark', "flowchart" : { "curve" : "basis" } } }%%
graph LR
    A(kcount)
    B(kfilter)
    C(kextract)
    D(kassemble)
    A --> B --> C --> D
```

### Study description

Here we implement the study by Neves et al. to detect a male linked 
genomic region involved in sex determination in the dioecious plant species *Amaranthus palmeri*. This study uses pool-seq to sequence four populations composing male and female plants from two geographically distinct populations.

<cite>Cátia José Neves, Maor Matzrafi, Meik Thiele, Anne Lorant, Mohsen B Mesgaran, Markus G Stetter, Male Linked Genomic Region Determines Sex in Dioecious Amaranthus palmeri, Journal of Heredity, Volume 111, Issue 7, October 2020, Pages 606–612, <a href=https://doi.org/10.1093/jhered/esaa047>https://doi.org/10.1093/jhered/esaa047</a></cite>


### Get the fastq data
If you wish to follow along you can dowload the data with these instructions.

??? abstract "download fastq data using wget"
    ```bash
    # make a directory to store the raw fastq data files
    mkdir -p ./fastq-data

    URLS=(
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR416/001/ERR4161581/ERR4161581_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR416/001/ERR4161582/ERR4161582_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR416/001/ERR4161583/ERR4161583_1.fastq.gz
        ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR416/001/ERR4161584/ERR4161584_1.fastq.gz
    )

    # download files to the specified fastq directory
    for url in ERR4161581 ERR4161582 ERR4161583 ERR4161584; do
        wget $url;
    done
    ```

??? abstract "or, download fastq data using sra-tools"

    Download the latest version of the sratools from [https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit](https://github.com/ncbi/sra-tools/wiki/01.-
    Downloading-SRA-Toolkit) by selecting the compiled binaries that are 
    appropriate for your system (e.g., Linux or MacOSX). (I really do 
    recommend that you use the latest version since this software is updated
    frequently and does not maintain compatibility with older versions. Follow
    the instructions to setup gcloud or aws to dramatically improve speed.)
    Then run the command below to download the fastq data for this study 
    into a new directory. The total filesize will be about 140Gb.

    ```bash
    # make a directory to store the raw fastq data files
    mkdir -p ./fastq-data

    # download files to the specified fastq directory
    for run in ERR4161581 ERR4161582 ERR4161583 ERR4161584; do
        fasterq-dump --progress --outdir ./fastq-data --temp /tmp $run;
    done
    ```

### Setup: imports and data files

Set the logging level on kmerkit to specify more or less logged info.
```python
import kmerkit
kmerkit.set_loglevel("INFO")
```

Load the files to a dictionary mapping sample names to a list of 
file names.
```python
# create a dictionary mapping study sample names to list of file paths
FASTQ_DICT = {
    "california-male": ["ERS4574576_R1.fastq.gz", "ERS4574576_R2.fastq.gz"],
    "california-female": ["ERS4574577_R1.fastq.gz", "ERS4574577_R2.fastq.gz"],
    "kansas-male": ["ERS4574578_R1.fastq.gz", "ERS4574578_R2.fastq.gz"],
    "kansas-female": ["ERS4574579_R1.fastq.gz", "ERS4574579_R2.fastq.gz"],
}
```

??? info "FASTQ_DICT grouping of sample names to list of file paths."

    ```python


    ```

### Create kmer databases
```python
# init counter class object
tool = Kcount(
    fastq_dict=FASTQ_DICT,
    name="hybridus", 
    workdir="/tmp/", 
    kmersize=35,
    trim_reads=True,
    mindepth=15,
    maxdepth=2000,
    maxcount=65535,       # max count possible with 16 bit integers
    canonical=True,
)

# run the analysis
tool.run()

# optionally access a stats summary for the run
print(tool.statsdf)
```

??? abstract "kmerkit logged output"
    ```console
    ...
    ```

### Filter kmers: unique to males

```python

tool = Kfilter(
    name="hybridus",
    workdir="/tmp",
    trait_dict={
        0: ["california-female", "kansas-female"],
        1: ["california-male", "kansas-male"],
    },
    mincov=0.0,        # no global min applied.
    minmap={
        0: 0.0,        # group 0 can have as little as 0 coverage
        1: 0.5,        # group 1 must have at least 50% coverage
    },
    maxmap={
        0: 0.0,        # group 0 must have 0 coverage
        1: 1.0,        # group 1 can have up to 100% coverage
    },
    mincov_canon={
        0: 0.0,        # no canonical min in group 0
        1: 0.5,        # kmer must be canonical in 50% of group 1 samples
    }
)
tool.run()
print(tool.statsdf)
```

??? abstract "kmerkit logged output"
    ```console
    ...
    ```
The resulting files are KMC database files written with the name prefix `{workdir}/{name}-{sample-name}-kfilter.[pre,suf]`

??? info "peek at workdir file structure"
    ```console
    .
    |---workdir
        |
        ------- a
        |-------b
        --------c
        ...
    ```


### Extract reads unique to males

```python
# set up filter tool
tool = Kextract(
    name="hybridus",
    workdir="/tmp",
    fastq_path=FASTQ_DICT,
    group_kmers="/tmp/kfilter_hybridus_filtered",
    mindepth=1,
)
tool.run()
print(tool.statsdf.T)
```

??? abstract "kmerkit logged output"
    ```console
    ...
    ```


### Assemble contigs from extracted reads
For every `--sample` provided to `kassemble` 

```console
kmerkit kextract \
    --name amaranth-dioecy \
    --workdir ./cookbook1 \
    --sample A ./data/sample-A.fastq.gz \
    --sample B ./data/sample-B.fastq.gz \
    --sample C ./data/sample-C.fastq.gz \
```
??? abstract "kmerkit logged output"
    ```console
    ...
    ```
