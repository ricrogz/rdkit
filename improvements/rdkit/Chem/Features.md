# Chem/Features

## Undefined names

### `ShowFeats.py`: arrowhead construction uses the wrong geometry module name

The module imports `Geometry` but `_buildCanonArrowhead()` refers to
`RDGeometry`, so the first feature-direction arrow raises `NameError`. Use the
imported module consistently. See
`patches/rdkit/Chem/Features/001-use-imported-geometry-module.patch`.

## Geometry errors

### `ShowFeats.py`: every generated arrowhead vertex uses a full-turn angle

The loop evaluates `sin(i * 2π)` and `cos(i * 2π)`, placing every base vertex
at the same point. Divide the turn by `nSteps` so the vertices span the circle.
See `patches/rdkit/Chem/Features/002-distribute-arrowhead-vertices.patch`.
