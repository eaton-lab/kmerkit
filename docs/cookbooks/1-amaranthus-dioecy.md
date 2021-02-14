

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

Explain and cite Amaranthus studies here...


### Get the fastq data
First we will download the fastq data files from NCBI. There are several ways
to do this, such as the sra-tools, but because these files are especially large
we will just use `wget` here. 

```bash
# make a directory for fastq data files
mkdir -p fastq_data/
URL="ftp://xxx.yyy.zzz"
wget $URL -O fastq_data/
```
??? abstract "wget stdout"
	```console
	...
	```

### Create kmer databases

```console
kmerkit kcount \
	--name amaranth-dioecy \
	--workdir ./cookbook1 \
	--kmersize 17 \
	--mindepth 5 \
	--sample ... \
	--sample ... \
	...
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

	sample  	male 	species 	
	sample-1	1		palmeri
	sample-2	1		palmeri
	sample-3	0		palmeri
	sample-4	0		palmeri
	...			
	```

```console
kmerkit kfilter \
	--name amaranth-dioecy \
	--workdir ./cookbook1 \
	--phenos_file phenos.tsv \    # file mapping sample names to traits
	--pheno male \                # select trait column "male"
	--minfreq 1 0.9 \             # kmer-freq>=0.9 of samples where male=1
	--maxfreq 1 1.0 \             # kmer-freq<=0.0 of samples where male=1
	--minfreq 0 0.0 \             # kmer-freq>=0.0 of samples where male=0
	--maxfreq 0 0.1	\             # kmer-freq<=0.1 of samples where male=0
	--mincanon 1 0.25			  # kmer-freq-canonical>=0.25 of samples where male=1
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
