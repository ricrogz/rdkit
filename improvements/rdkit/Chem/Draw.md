# Chem/Draw

## Performance improvements

### `SimilarityMaps.py`: fingerprint map computes an unused distance matrix

`GetSimilarityMapForFingerprint()` builds the probe molecule's full distance
matrix but never reads it. This is quadratic work and memory for every map;
remove the call. See
`patches/rdkit/Chem/Draw/001-remove-unused-distance-matrix.patch`.

### `__init__.py`: SVG grid allocates an unused cell list

`_MolsToGridSVG()` allocates `nRows * molsPerRow` empty strings in `blocks` but
never reads the list. Remove the allocation. See
`patches/rdkit/Chem/Draw/002-remove-unused-svg-blocks.patch`.
