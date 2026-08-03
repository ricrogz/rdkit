# Chem/Subshape

## Incorrect loop bounds

### `SubshapeAligner.py`: `orderedTraversal=False` ignores its computed start index

`_getAllTriangles()` computes `kStart` differently for ordered and unordered
traversal, then always starts at `j + 1`. Use `kStart` so the unordered branch
actually explores the intended permutations. See
`patches/rdkit/Chem/Subshape/001-use-selected-triangle-start.patch`.
