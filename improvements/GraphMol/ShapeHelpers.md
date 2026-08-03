# ShapeHelpers

## Redundant includes

`ShapeEncoder.cpp` repeats `Transform3D.h` and `UniformGrid3D.h`. See
`patches/GraphMol/ShapeHelpers/001-remove-duplicate-geometry-includes.patch`.
