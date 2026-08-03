# Chem/MolDb

## Undefined variables

### Loader implementations compute coordinates for `m` instead of `mol`

Both `Loader_orig.py` and `Loader_sa.py` call `Compute2DCoords(m)` in their
redraw branch, but their molecule parameter is named `mol`. Any load using
`redraw=True` raises `NameError`. See
`patches/rdkit/Chem/MolDb/001-redraw-current-molecule.patch`.
