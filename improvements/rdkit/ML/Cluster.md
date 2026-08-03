# ML/Cluster

## Python 3 index errors

### Lower-triangular matrix indices use true division

`ClusterUtils.py` and `Resemblance.py` compute packed-matrix indices with `/`.
Python 3 returns a float, which cannot index a list or array. Use integer
division for the triangular-number calculation. See
`patches/rdkit/ML/Cluster/001-use-integer-packed-matrix-indices.patch`.
