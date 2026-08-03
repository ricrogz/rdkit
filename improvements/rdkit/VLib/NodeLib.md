# VLib/NodeLib

## Undefined names on validation failures

### `DbPickleSupplier.py`: malformed query results mask their error with `NameError`

Both lazy result initialization paths format `res`, which is undefined, rather
than `self.res`. A result with too few fields therefore loses the intended
`ValueError`. See
`patches/rdkit/VLib/NodeLib/001-report-current-query-result.patch`.
