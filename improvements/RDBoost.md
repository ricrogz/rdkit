# RDBoost

## Memory leaks on exceptions

### `Wrap.cpp`: atom-map conversion owns its result with a raw pointer

`translateAtomMap()` allocates the result before converting every Python row.
If sequence extraction or indexing throws, the raw allocation leaks. Manage
the vector with `unique_ptr` during construction and release it only on the
successful raw-pointer return path. See
`patches/RDBoost/001-own-atom-map-during-conversion.patch`.
