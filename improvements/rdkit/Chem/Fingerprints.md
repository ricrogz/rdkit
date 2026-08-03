# Chem/Fingerprints

## Undefined variables

### `MolSimilarity.py`: database fallback passes nonexistent `data`

When the optional PostgreSQL lazy sequence is unavailable, `GetFingerprints()`
passes an undefined variable to `ForwardDbFpSupplier`. Wrap the already-built
cursor and SQL command in the standard forward result-set class instead. See
`patches/rdkit/Chem/Fingerprints/001-build-forward-db-result.patch`.

### `FingerprintMols.py`: parse failure reports nonexistent `smi`

The `(ID, mol)` loop logs `smi` when `mol` is false, masking the actual parse
failure with `NameError`. Report the available record ID. See
`patches/rdkit/Chem/Fingerprints/002-report-failed-record-id.patch`.

## Python 3 iterator compatibility

### `FingerprintMols.py`: SD supplier loop calls legacy `.next()`

Use the built-in `next()` so the loop works with Python 3's iterator protocol.
See `patches/rdkit/Chem/Fingerprints/003-use-next-builtin.patch`.
