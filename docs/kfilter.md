








## CLI


### Entering traits as a CSV file

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




## API 

```python

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

```

