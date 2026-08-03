# Chem/Suppliers

## Python 3 iterator compatibility

### `DbMolSupplier.py`: random-access supplier calls removed iterator `.next()`

Python 3 iterators expose `__next__()` through the built-in `next()`, not a
`.next()` method. The current path fails for ordinary result-set iterators. See
`patches/rdkit/Chem/Suppliers/001-use-next-builtin.patch`.
