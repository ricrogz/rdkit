# RGroupDecomposition

## Uninitialized reads

### `RGroupGa.cpp`: chromosome fitness has no initial value

`info()` and `getFitness()` can read `fitness` before the first call to
`score()`. Initialize it to zero at construction so diagnostic or selection
code cannot consume an indeterminate `double`. See
`patches/GraphMol/RGroupDecomposition/001-initialize-chromosome-fitness.patch`.
