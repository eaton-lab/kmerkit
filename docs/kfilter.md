












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

