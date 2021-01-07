



### Outline notes


1. Count kmers in every individual using kmc
2. Encode kmer presence/absence as genotype matrix; shape=(samples x kmers)
3. Convert genos to PLINK binary ped format.
3. Calculate kinship matrix from genotypes using gemma


```bash
# option one using plink files
gemma -bfile plink-binary-ped-file-prefix -gk 1 -o outprefix

# option 2 using bimbam files



# gemma manual section 4.4
# There  will  be  two  output  files,  both  inside  an  output  folder  in  the  current  directory. The  pre-fix.log.txt file contains some detailed information about the running parameters and computationtime, while the prefix.cXX.txt or prefix.sXX.txt contains an×nmatrix of estimated relatednessmatrix.
```

4. Perform association test w/ univariate linear mixed model.
```bash
gemma \
    -bfile plink-binary-ped-file-prefix \
    -k kinship-filename \
    -lmm 2 \
    -o outprefix

# gemma manual section 4.6
# For binary traits, one can label controls as 0 and cases as 1, and follow our previous approachesto fit the data with a linear mixed model by treating the binary case control labels as quantitativetraits
```


### PLINK format
GEMMA takes PLINK BINARY PED FORMAT for BOTH geno and pheno data. This
is entered by providing a prefix for three files that end with .bed,
.bim and .fam. 

# conda install vcftool -c conda-forge -c bioconda
# conda install plink -c conda-forge -c bioconda


#### PED FILE FORMAT

The long form is:
fam ind pat mat sex pheno genos...
```bash
plink --file example --1
```


But, we don't need most of this, so we can instead do:





----------------------------

### BIMBAM format

- mean genotype file:
No header.
geno-name, allele1, allele2, mean geno of ind1, mean genos of ind2, ...

```
kmer1, 0, 1, 0, 0, 0, 0, 1
kmer2, 0, 1, 0, 0, 0, 1, 1
kmer3, 0, 1, 0, 0, 1, 1, 1
```



-----------------------------

