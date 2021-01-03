
## Kmerkit: kmer assocation toolkit in Python

#### Under active development!

Kmerkit is a general toolkit for performing reference-free 
genome-wide association analyses from kmers. 
It uses [KMC](https://github.com/refresh-bio/KMC/)
to count and filter kmers, and has options to perform associations
using either [gemma]() or by machine learning methods in [scikit-learn]().
Code and tutorials are currently under development.


#### Installation with conda (coming soon)
```bash
# conda install kmerkit -c conda-forge -c bioconda

# for now, do dev installation with pip
git clone https://github.com/eaton-lab/kmerkit
cd kmerkit
pip install -e .
```


#### Interactive analysis in jupyter
The preferred way to run analyses in kmerkit is to use the API 
interactively in a jupyter notebook. This allows access to statistics,
plotting summaries, and encourages users to create reproducible documentation
of their workflow.


```python
import kmerkit  

# DATA
FASTQS = "/tmp/*.fastq.gz"
PHENOS = "/tmp/phenos.csv"

# get dict mapping {sample_names: [fastq_files]}
fastq_dict = kmerkit.get_fastq_dict_from_path(
    fastq_path=FASTQS, 
    name_split="_R",
)

# count kmers
kmerkit.Kcount(
    name='test', 
    workdir='/tmp', 
    fastq_dict=fastqdict,
    kmersize=31,
).run()

# find kmers unique to one group versus another
kmerkit.Kfilter(
    name='test', 
    workdir='/tmp', 
    phenos=PHENOS,
    trait='trait',
    mincov=0.25,               # must be present in 25% overall
    mincanon=0.25,             # must exist in both forms in 25% of samples where present.
    minmap={1: 0.5, 0: 0.0},   # must be present 50% in group 1
    maxmap={1: 1.0, 0: 0.0},   # must be absent in group 0
).run()

# get fastq reads filtered to only those matching target kmers
kmerkit.Kextract(
    name='test',
    workdir='/tmp',
    fastq_dict=fastq_dict,
    group_kmers="/tmp/kgroup_test",
).run()  

# get matrix of (nsamples x nkmers) as geno data, from Kfiltered kmers.
kmerkit.Kmatrix(
    name="test",
    workdir="/tmp",
    ...
)
```

### Code in development

  - Code for gwas/regression, incorporating covariance matrix (phylogeny or kinship)
```python

# get location of kmer-matched reads mapping on reference scaffolds
kmpy.Kmap(...)

# assemble kmer-matched reads into denovo contigs
kmpy.Kassemble(...)

# run plink and gemma with geno and pheno data
kmpy.Kgwas(...)

# visualize sample/group geno data (toyplot, scikit-learn etc.)
kmpy.Klearn(...)
```


## Applications in Eaton lab...

#### 1. Pedicularis RAD data:
  - count 31-mers in each sample (Kcount)
  - filter to find kmers unique to forked pops (Kfilter)
  - extract reads containing fork-kmers (Kextract)
  - map fork-kmer-reads to reference genome and annotation (Kmap)
  - write fork-kmers as genotype matrix (Kmatrix)
  - test for associations using gemma (Kgwas)

```python

import kmpy

# prepare files
FASTQS = "/pinky/.../Pedicularis_*.fastq.gz"
PHENOS = "/pinky/.../Pedicularis_forked.csv"

fastq_dict = kmpy.get_fastq_dict_from_path(
    fastq_path=FASTQS,
    name_split="_R",
)

# count kmers in each sample and write database files
kmpy.Kcount(
    name="fork", 
    workdir="/pinky/", 
    fastq_path=FASTQS, 
    kmersize=31,
).run()

# filter kmers to find those associated with forked beaks
kmpy.Kfilter(
    name="fork", 
    workdir="/pinky/", 
    phenos=PHENOS,
    trait="forked", 
    mincov=0.2,          # globally exclude rare kmers
    mincanon=0.0,        # RAD is strand-specific, so turn this filter off
    minmap={'forked': 0.75, 'non-forked': 0.0},   # min allows some to be missing randomly
    maxmap={'forked': 1.00, 'non-forked': 0.0},   # max allows none to occur in non-forked
).run()

# extract reads from fastq files that contain target kmers
kmpy.Kextract(
  name="fork", 
  workdir="/pinky/", 
  fastq_dict=fastq_dict,
).run()

```



#### 2. Amaranth palmeri & tuberculatus data
- rerun this analysis using 15-mers, to get 2 kmer files (palmeri and tuberc)
  - count kmers in each sample
  - group by M versus F/H
  - target kmers are in 'intersect' on male plants.
  - extract reads with kmers matching


```python

import kmpy

# prepare files
FASTQS = "/pinky/.../Amaranth_fastqs/*.fastq.gz"
PHENOS = "/pinky/.../Amaranth_sampled_sexed.csv"

# get dict of {sample_name: [fastq_files]}
fastq_dict = kmpy.get_fastqdict_from_path(
    fastq_path=FASTQS,
    name_split="_R"
)

# count kmers in each sample and write database files
km.Kcount(
    name="dioecy", 
    workdir="/pinky/", 
    fastq_dict=fastq_dict, 
    kmersize=15, 
).run()

# intersect is no good for RAD data, instead we should do 'union' and require that 
# a kmer is found in at least xx coverage across yy samples.
km.Kfilter(
    name="dioecy", 
    workdir="/pinky/", 
    phenos=mindepth=2, 
    trait="male", 
    ...
).run()

# get original reads for each sample containing target kmers
km.Kextract(...)
```
  
  
  
#### 3. New Amaranth spp. data

- Analyze known male-specific mers in new data
  - extract all fastq reads (Python) from each sample that contain any 15-mers from analysis 2.
  - assemble contigs from those reads (Platanus)
  - annotate and find genes in each contig and find shared genes among spp (Montgomery-like)
  - Build phylogenies for these genes.
      - align gene sequences (muscle)
      - infer gene trees (raxml)


- Search for kmers associated with male/female&mono in pooled spp. (kmergwas)
  - group samples into male pool or female/mono pool
  - get 15-mers in each pool
  - find kmers uniquely shared among males with FILTERS applied (custom scripts)
  
  
- Search for kmers associated with dioecy/mono in pooled spp. (kmergwas)
  - same as above, but different way of pooling samples
  
  
