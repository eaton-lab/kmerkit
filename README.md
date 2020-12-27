# kmertesting
Testing scripts for kmer research


### Code in testing
```python

import kmpy

# DATA
FASTQS = "/tmp/*.fastq.gz"
PHENOS = "/tmp/phenos.csv"

# count kmers
kmpy.Kcount(
  name='test', 
  workdir='/tmp', 
  fastq_path=FASTQS, 
  kmersize=35, 
  name_split="_R",
).run()

# find kmers unique to one group versus another
kmpy.Kgroup(
  name='test', 
  workdir='/tmp', 
  phenos=PHENOS,
  trait='trait',
  operation_g0="union",
  operation_g1="intersect", 
  operation_g0g1="subtract",
  mindepth_g0=5,
  mindepth_g1=5,
  mindepth_g0g1=None,
  maxdepth_g0=100,
  maxdepth_g1=100,
  maxdepth_g0g1=None,
  reverse=True,
  force=True,
).run()

# get fastq reads filtered to only those matching target kmers
kmpy.Kfilter(
  name='test',
  workdir='/tmp',
  fastq_path=FASTQS,
  group_kmers="/tmp/kgroup_test",
  name_split="_R",
  mindepth=5,
).run()  

```

### Code in development

  - Code for gwas/regression, incorporating covariance matrix (phylogeny or kinship)
```python

# get matrix of (nsamples x nkmers) as geno data.
kmpy.Kmatrix(...)


# run plink and gemma with geno and pheno data
kmpy.Kregression(...)


# visualize sample/group geno data (scikit-learn, toyplot, etc.)
kmpy.Kanalyze(...)



```


### Thoughts on applications in Eaton lab...

#### 1. Pedicularis RAD data:
  - count 35-mers in each sample
  - group 100 samples into forked or non-forked pools
  - identify fork-associated kmers 
  - get reads containing kmers (Python)
  - map reads to reference to identify scaffolds with many kmers (Python)
  - search annotation for genes in these regions.

```python

import kmpy

# prepare files
FASTQS = "/pinky/.../Pedicularis_*.fastq.gz"
PHENOS = "/pinky/.../Pedicularis_forked.csv"

# count kmers in each sample and write database files
km.Kcount(name="fork", workdir="/pinky/", fastq_path=FASTQS, kmersize=35, name_split="_R").run()

# intersect is no good for RAD data, instead we should do 'union' and require that 
# a kmer is found in at least xx coverage across yy samples.
km.Kgroup(
  name="fork", 
  workdir="/pinky/", 
  phenos=PHENOS,
  trait="forked", 
  operation_g0="union",               # all kmers in samples w/o fork
  operation_g1="union",               # all kmers in samples w/  fork
  operation_g0g1="counters_subtract", # operate on counts, not presence/absence
  mindepth_g0=1,                      # subtract kmer if present in ANY non-fork
  mindepth_g1=10,                     # require kmer present in at least ... nsamples
  mindepth_g0g1=50,                   # tweak this: how much more common is kmer in 1 than 0?
  reverse=True,
).run()

# get original reads for each sample containing target kmers
km.Kfilter(name="fork", workdir="/pinky/newfastqs/", fastq_path=FASTQS, name_split="_R").run()

```



#### 2. Amaranth palmeri & tuberculatus data
- rerun this analysis using 15-mers, to get 2 kmer files (palmeri and tuberc)
  - Trim reads (maybe can skip...)
  - Pool into male and female pops
  - Jellfish each pool with k=15
  - Custom scripts - Python - identify unique male kmers given filter threshold


```python

import kmpy

# prepare files
FASTQS = "/pinky/.../Amaranth_fastqs/*.fastq.gz"
PHENOS = "/pinky/.../Amaranth_sampled_sexed.csv"

# count kmers in each sample and write database files
km.Kcount(name="dioecy", workdir="/pinky/", fastq_path=FASTQS, kmersize=35, name_split="_R").run()

# intersect is no good for RAD data, instead we should do 'union' and require that 
# a kmer is found in at least xx coverage across yy samples.
km.Kgroup(
  name="dioecy", 
  workdir="/pinky/", 
  phenos=mindepth=2, 
  trait="male", 
  operation_g0="union",       # all kmers in samples that are female/hermaph.
  operation_g1="intersect",   # shared kmers in samples that are male
  operation_g0g1="subtract",
  mindepth_g0=1,              # subtract kmer if present in ANY non-fork
  mindepth_g1=10,             # require kmer present in at least ... nsamples
  reverse=True,
).run()

# get original reads for each sample containing target kmers
km.Kfilter(name="fork", workdir="/pinky/newfastqs/", fastq_path=FASTQS, name_split="_R").run()

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
  
  
