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

```python
import glob
import kmpy as km
allfastqs = glob.glob("Pedicularis_*.fastq.gz")
km.Kcount(allfastqs, kmersize=35, workdir="/tmp", name_split="_R")
km.Kgroup(forked.csv, mindepth=2, workdir="/tmp")
km.Kcompare(workdir="/tmp")
km.Kmapping(...)
```
  - group 100 samples into forked or non-forked pools (bash scripts)
  - count 35-mers in each pool (jellyfish)
  - identify fork-associated kmers (custom scripts)
  - get reads containing kmers (Python)
  - map reads to reference to identify scaffolds with many kmers (Python)
  - search annotation for genes in these regions.



#### 2. Amaranth palmeri & tuberculatus data
- rerun this analysis using 15-mers, to get 2 kmer files (palmeri and tuberc)
  - Trimmomatic the reads
  - Pool into male and female pops
  - Jellfish each pool with k=15
  - Custom scripts - Python - identify unique male kmers given filter threshold
  
  
  
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
  
  
