


The `kcount` tool in `kmerkit` uses the `KMC` software tool [citation](cite) 
to create a kmer database for one or more samples from the input of 
one or more fastq/a data files. The kmer databases are used to extract 
kmer counts or statistics, and for downstream analysis with other `kmerkit`
tools. 
<!-- As a wrapper around `KMC`, `kmerkit` is intended to make it easier -->


## Command line usage (CLI)
```bash
kmerkit kcount --help
```


### Entering data and sample names
`kmerkit` offers a flexible approach to entering sample names
and the fastq data files associated with them.
This makes it easy to accommodate single or paired-end data 
files, to combine technical replicates from different sequencing
runs, or even to create pooled samples. By supporting 
wildcard/regex selectors it also helps to avoid typos. 

!!! important optional-class "Three options for data entry"
	
	=== "-s name file file..."
		The `-s` sample option takes a `sample_name` argument followed by one or more `filepath` arguments to assign files to sample names. The other two options for entering data are simply shortcuts to achieving the same with less typing.
		```bash
		kmerkit kcount \
			... \
			-s sample-1 ./data/sample-1_R1.fastq.gz ./data/sample-1_R2.fastq.gz \
			-s sample-2 ./data/sample-2_R1.fastq.gz ./data/sample-2_R2.fastq.gz \
			-s sample-3 ./data/sample-3_R1.fastq.gz ./data/sample-3_R2.fastq.gz
		```

	=== "-s name file*"
		The `-s` sample option takes two arguments: sample_name and 
		`filepath`. The latter can be set as a wildcard/regex type
		string to select multiple file names that will each be assigned
		to the designated sample name.
		```bash
		kmerkit kcount \
			... \
			-s sample-1 ./data/sample-1_R*.fastq.gz \
			-s sample-2 ./data/sample-2_R*.fastq.gz \
			-s sample-3 ./data/sample-3_R*.fastq.gz		
		```

	=== "--fastq_path splitter path*"
		This method uses the `--fastq_path` flag and takes two entries
		as arguments (`splitter` and `filepath`). All files matching the filepath regex will be selected and sample names will be extracted from filenames by splitting on the 'splitter' character. Multiple files matching to the same sample name prefix (e.g., paired-end file names if splitting on '_R') 
		will be assigned to the same sample.
		```bash
		kmerkit kcount \
			...
			--fastq_path _R ./data/sample-*_R*.fastq.gz
		```


## API usage
Import the `kmerkit` package to access its modules.
```python
import kmerkit
``` 

### Entering fastq data and sample names
#### approach 1: extract sample names from file names
Many times your filenames will contain the sample names or IDs, and
thus a convenient way to map sample names to their files is to extract
the sample names from the files themselves.

For this approach we will assume that all of your fastq files
are stored together in a single directory and you want to select
all or a subset of them. They may include single or 
paired-end files (where paired filenames are the same except for 
R1 and R2 in their names).


??? example optional-class "Example: wildcard paths"

	=== "/tmp/*.fastq.gz"
	    ```bash
		/tmp/sample-1_R1.fastq.gz
		/tmp/sample-1_R2.fastq.gz
		/tmp/sample-2_R1.fastq.gz
		/tmp/sample-2_R2.fastq.gz			
		/tmp/sample-3_R1.fastq.gz
		/tmp/sample-3_R2.fastq.gz			
		/tmp/sample-4_R1.fastq.gz
		/tmp/sample-4_R2.fastq.gz			
	    ```

	=== "/tmp/sample[1,2]*.fastq.gz"
		```bash
		/tmp/sample-1_R1.fastq.gz
		/tmp/sample-1_R2.fastq.gz
		/tmp/sample-2_R1.fastq.gz
		/tmp/sample-2_R2.fastq.gz			
		```

	=== "/tmp/sample[1,2]_R1.fastq.gz"
		```bash
		/tmp/sample-1_R1.fastq.gz
		/tmp/sample-2_R1.fastq.gz
		```




```python
# A string file path using wildcard selectors (*) to match multilple filenames
FASTQS = "/tmp/*.fastq.gz"
```

You can then parse the file names based on a separator to generate prefix names. In this example we split on the `"_R"` delimiter. This
splits the filename to keep anything before the delimiter. Paired-end filenames will be the same, and thus be assigned to the same sample. Depending on the delimiter selected you may also be able to assign technical replicates to the same sample. 

```python hl_lines="2 3"
# get dict mapping {sample_names: [fastq_files]}
fastq_dict = kmerkit.get_fastq_dict_from_path(
    fastq_path=FASTQS, 
    name_split="_R",
)
```


??? note optional-class "fastq_dict contents"
	```python
	# fastq_dict is a dictionary mapping sample names to lists of files paths
	{
		"sample-1": ["/tmp/sample-1_R1.fastq.gz", "/tmp/sample-1_R2.fastq.gz"],
		"sample-2": ["/tmp/sample-2_R1.fastq.gz", "/tmp/sample-2_R2.fastq.gz"],
		"sample-3": ["/tmp/sample-3_R1.fastq.gz", "/tmp/sample-3_R2.fastq.gz"],
		"sample-4": ["/tmp/sample-4_R1.fastq.gz", "/tmp/sample-4_R2.fastq.gz"],
	}
	```

This method for assigning fastq files to sample names is only a convenience. In many cases your filenames, or the way in which they are
organized on your system may prevent you from using this method. In 
that case you should you use the second method (below), where you 
simply create the fastq dictionary manually.


```

```

The `kcount` module includes a main Kcount class object for running 
kmer counting operations and viewing that stats associated with each 
sample.
```python
kco = kmerkit.kcount(
	name="test",
	workdir="/tmp",
	fastq_dict = 

)


```