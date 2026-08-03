# Chem/Pharm3D

## Uninitialized local variables

### `EmbedLib.py`: verbose matching reads `bm` before assigning it

The verbose `pre update` diagnostic iterates `bm` before the current match has
created that bounds matrix, raising `UnboundLocalError`. It is intended to show
the input bounds, so iterate `bounds` there. See
`patches/rdkit/Chem/Pharm3D/001-fix-pre-update-bounds-diagnostic.patch`.
