# kmertesting
Testing scripts for kmer research

## TODO list:


#### 0. Code development: ideas for classes

  - A class for counting kmers to create kmc databases in workdir
```python
# example is something like this...
Kcount(fastqs, kmersize=17, workdir="/tmp", name_split="_R")
```
  - Code for merging databases into groups based on input pheno file.
```python
Kgroup(phenos.csv, mindepth=2, maxdepth=1000, workdir="/tmp")
```
  - figure out whether we need non-canonical kmers...
  - Code for simple presence-absence search with cutoffs (e.g., dioecy study)
```python
Kcompare(...)
```
  - Code for gwas/regression, incorporating covariance matrix (phylogeny or kinship)
```python
Kregression(...)
```
  - Code for analysis and visualization (scikit-learn, toyplot, etc.)
```python
Kanalyze(...)
```

  - Code to get fastqs containing kmers of interest
```python
Kmapping(...)
```


#### 1. Pedicularis study

```python
import glob
import kmertesting as km
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
  
  
