

### Workflow diagram

```mermaid
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


### Initialize a new project

```console
kmerkit init --name dioecy --workdir /tmp ./data/amaranth/*.gz
```

### Count kmers

```console
kmerkit count --json /tmp/dioecy.json --kmer-size 35 --min-depth 5 
```

??? abstract "kmerkit logged output"
	```console
	...
	```


### Filter kmers for those unique to males vs females

??? "Create a phenotype file (phenos.tsv)"
	This should be a tab (or any whitespace) separated values in a table with
	sample names in the first row (the column name for this row is ignored)
	and then trait values in subsequent columns. In the example below we
	will select the column "male" but it is fine for other columns of data
	to be present in the file which can be used in other steps (e.g., GWAS).

	```console

	name        trait
	sample-1	1
	sample-2	1
	sample-3	0
	sample-4	0
	...			
	```

```console
kmerkit filter \
	--json /tmp/dioecy.json \
	--trait trait.csv \
	--minmap 0.0 0.9 \            # kmer-freq>=0.9 of samples where trait=1
	--maxmap 0.1 1.0              # kmer-freq<=0.1 of samples where trait=0
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


### Extract reads containing target kmers from samples

```console
kmerkit extract \
	--json /tmp/dioecy.json \
	--min-kmers-per-read 5 \
	./data/amaranths/*.gz
```

??? abstract "kmerkit logged output"
	```console
	...
	```


### Assemble contigs from extracted reads

```console
kmerkit assemble \
	--json /tmp/dioecy.json \
	--assembler
```

??? abstract "kmerkit logged output"
	```console
	...
	```
