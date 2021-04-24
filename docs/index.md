


# Welcome to kmerkit

`kmerkit` is a toolkit for performing evolutionary analyses using 
kmer counts, frequencies, and comparisons with or without the 
context of a reference genome. It is both a stand-alone analysis framework 
as well as an extendable toolkit that can be incorporated into other 
software. The `kmerkit` tools are available as both a command-line interface
as well as a Python API.


Core pipelines:
```mermaid
%%{init: {'theme': 'dark', "flowchart" : { "curve" : "basis" } } }%%
graph LR
   A0(kinit)
   A1(ktrim)
   A(kcount)
   B(kfilter)
   C(kextract)
   D(kassemble)
   E(ktree)
   F(kmap)
   G(kannotate)
   H(kmatrix)
   I(kgwas)
   J(klearn: API)
   A0 --> A1
   A1 --> A
   A0 --> A
   A --> B
   B --> C
   C --> D
   A --> E
   B --> H
   E --> C
   E --> H
   C --> F
   F --> G
   D --> G
   H --> I
   H --> J

linkStyle default stroke-width:2px,fill:none,stroke:grey;   
```

Other convenience utilities:
```mermaid
%%{init: {'theme': 'dark', "flowchart" : { "curve" : "basis" } } }%%
graph TB
   A(kdump)
   B(kstats)
   C(kinfo)
```