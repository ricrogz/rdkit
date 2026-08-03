# Chem/fmcs

## Python 3 runtime errors

### `fmcs.py`: progress reporting still uses Python 2 `print >>` syntax

Five diagnostic paths parse under Python 3 as shift expressions involving the
`print` function and fail when invoked. Convert them to `print(...,
file=sys.stderr)`. See
`patches/rdkit/Chem/fmcs/001-modernize-stderr-reporting.patch`.
