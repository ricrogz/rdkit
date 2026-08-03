# SmilesParse

## Dangling pointers

### `SmartsWrite.cpp`: error path throws a pointer into a destroyed string

The unsupported-query-bond path throws `msg.str().c_str()`. The temporary
`std::string` is destroyed at the end of the throw expression, so the caught
`const char *` is immediately dangling. Throw the module's owning
`ValueErrorException` with the constructed string, consistent with other
SMARTS-writing failures. See
`patches/GraphMol/SmilesParse/001-throw-owning-smarts-error.patch`.

## Performance improvements

### `SmartsWrite.cpp`: recursive SMARTS combination copies three strings per node

`_combineChildSmarts()` takes both child strings and the query description by
value even though it never modifies them. This helper is called recursively
across query trees, so large queries repeatedly copy growing strings. Passing
them by const reference removes those allocations without changing behavior.
See `patches/GraphMol/SmilesParse/002-pass-smart-parts-by-reference.patch`.

## Redundant includes

`SmilesParse.cpp` and `CXSmilesOps.cpp` each include `Chirality.h` twice. See
`patches/GraphMol/SmilesParse/003-remove-duplicate-chirality-includes.patch`.
