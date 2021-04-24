


## Installation Instructions
Recommended installation is with `conda` to ensure that all dependencies
are pulled in correctly.

=== "conda"
    ```bash
    # (NOT YET AVAILABLE)
    conda install kmerkit -c conda-forge -c bioconda    
    ```


=== "git (for developers)"  
    ```bash
    # install dependencies from conda
    conda install kmc pandas toytree loguru typer -c conda-forge -c bioconda

    # install kmerkit locally from github main branch
    git clone https://github.com/eaton-lab/kmerkit
    cd kmerkit/
    pip install -e . --no-deps
    ```


## Dependencies
These are all *installed automatically by conda* when you install 
`kmerkit` with the `conda` instructions above. 
Care was taken to select dependencies that are stable
and minimal as possible.

#### binaries
- kmc: external tool for kmer counting and set operations
- gemma: linear model inference 

#### Python
- numpy: math and array operations
- scipy: statistical distributions
- pandas: tabular data structures
- scikit-learn: ML model inference
- toyplot: minimalist plotting library
- toytree: minimalist tree class and plotting
- loguru: logging 
- typer: type-checking and CLI 
- pydantic: type-checking and JSON serialization
<!-- - requests: API development -->
<!-- - scikit-learn: statistical inference -->