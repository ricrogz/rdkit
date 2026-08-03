# SimDivPickers

## Memory leaks on allocation failures

### `HierarchicalClusterPicker.cpp`: partial scratch allocation leaks

The clustering routine performs three `calloc` calls and checks them only
after all allocations. If the second or third allocation fails, the invariant
throws while earlier buffers remain allocated. Use vectors for the Fortran
scratch arrays so every exit path releases acquired memory automatically. See
`patches/SimDivPickers/001-use-raii-clustering-buffers.patch`.
