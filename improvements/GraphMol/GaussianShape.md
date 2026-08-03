# GaussianShape

## Uninitialized state

### `ShapeInput.h`: serialization constructor starts with indeterminate indices

The string constructor relies on archive extraction to populate
`d_activeShape`, `d_numAtoms`, and `d_numFeats`, but they have no defaults. A
partial/older archive or an exception during extraction can leave scalar state
indeterminate. Neutral member initializers also make future constructors safer.
See `patches/GraphMol/GaussianShape/001-initialize-shape-indices.patch`.
