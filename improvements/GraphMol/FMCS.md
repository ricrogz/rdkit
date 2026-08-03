# FMCS

## Uninitialized state

### `MaximumCommonSubgraph.h`: reset MCS objects contain an indeterminate pointer

The nested `MCS` aggregate does not initialize `QueryMolecule`.
`MaximumCommonSubgraph::find()` resets the current result with `McsIdx = MCS()`
and later tests `McsIdx.QueryMolecule`; when no seed has populated the result,
that read is undefined behavior and can suppress the single-atom fallback or
lead later code to use a garbage molecule pointer. The outer query pointer and
threshold are likewise unsafe before `find()` initializes them. Give all three
scalar members deterministic defaults. See
`patches/GraphMol/FMCS/001-initialize-mcs-state.patch`.

## Undefined behavior

### `Seed.cpp`: 64 external bonds pass a guard and trigger a width-sized shift

`Composition2N` needs to compute `1ULL << NewBonds.size()`. The caller rejects
sizes greater than 64, but a size equal to 64 also shifts by the full width of
the type and is undefined; moreover, `2^64` cannot be represented in its
counter. Reject 64 as well as larger sizes. See
`patches/GraphMol/FMCS/002-reject-width-sized-composition.patch`.
