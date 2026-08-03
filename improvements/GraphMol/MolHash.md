# MolHash

## Integer and formatting correctness

### `hashfunctions.cpp`: hash suffix passes unsigned values to `%d`

The heteroatom-tautomer hash formats an unsigned hydrogen count with `%d`; the
prototype expression also performs subtraction in unsigned arithmetic when a
charge is present. Both can produce undefined or unintended results. Build the
suffix with `std::to_string()` and perform the subtraction in a signed
wide type. See `patches/GraphMol/MolHash/001-build-hash-suffix-type-safely.patch`.

## Stack exhaustion

### `hashfunctions.cpp`: molecule-sized scratch arrays use `alloca`

The extended Murcko scaffold routines allocate one byte per atom on the stack,
including a fresh allocation for every depth-first search. Large or adversarial
molecules can exhaust the stack, and `alloca` is non-standard. Use zero-
initialized vectors and pass their storage to the existing helpers. See
`patches/GraphMol/MolHash/002-replace-molecule-sized-alloca.patch`.
