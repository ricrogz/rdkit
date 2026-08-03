# Descriptors

## Uninitialized reads

### `GETAWAY.cpp`: short molecular diameters leave descriptor bins uninitialized

The code only assigns `HATSk`, `Hk`, `Rk`, and `Rp` bins reached by the
molecule's topological diameter, but later aggregates every bin. Molecules
whose diameter is below eight therefore read indeterminate doubles and can
produce nondeterministic descriptors. Value-initialize all four arrays so
unreached distance bins contribute zero. See
`patches/GraphMol/Descriptors/001-zero-getaway-distance-bins.patch`.
