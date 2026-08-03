# Chem/DSViewer

## Python 3 compatibility

### `DSViewer.py`: selection lookup calls removed `dict.iteritems()`

The molecule-selection path raises `AttributeError` under Python 3. Iterate
with `items()` instead. See
`patches/rdkit/Chem/DSViewer/001-use-dict-items.patch`.
