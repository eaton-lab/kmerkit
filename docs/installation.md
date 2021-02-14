


## Installation
Recommended installation is with `conda` to ensure that the dependencies
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
These are all installed automatically when you install `kmerkit` with 
`conda`. 

- kmc: external tool for kmer counting and set operations
- numpy: math operations
- pandas: organize results into tables
- toyplot: plotting library
- toytree: tree-based operations and plotting
- loguru: logging 
- typer: command-line-interface
<!-- - requests: API development -->
<!-- - scikit-learn: statistical inference -->