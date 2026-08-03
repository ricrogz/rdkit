# ChemTransforms

## Correctness defects

### `MolFragmenter.cpp`: subset masks use a platform-dependent signed literal

`fragmentOnSomeBonds()` stores masks in `uint64_t` but creates every bit with
`0x1L`. On LLP64 platforms such as Windows, `long` has only 32 bits, making
shifts for bond indices 32 through 62 undefined even though the function
explicitly permits up to 63 bonds. On LP64, shifting the signed value into bit
63 is also undefined when `maxToCut` is out of range. Use a `uint64_t` one-bit
value and return early when more cuts are requested than candidate bonds. See
`patches/GraphMol/ChemTransforms/006-use-uint64-fragment-masks.patch`.

### `ChemTransforms.cpp`: partial manual matches can dereference `end()`

`replaceCore()` looks up the molecule-side match for a multiply mapped core
atom's neighbor and immediately accesses `->second`. The public overload accepts
a caller-supplied `MatchVectType` and only validates indices, so a partial match
can omit that neighbor and make the lookup return `end()`. Reject the malformed
match with `ValueErrorException` before dereferencing it. See
`patches/GraphMol/ChemTransforms/005-validate-multiply-mapped-neighbor.patch`.

### `MolZip.cpp`: linker bond validation assigns instead of compares

The invariant intended to reject inconsistent duplicate linker bond types uses
`bondType = bond.linkerBondType`. It overwrites the observed type and treats the
assigned enum value as a Boolean, so most mismatches pass undetected (and an
`UNSPECIFIED` stored type fails for the wrong reason). Use equality comparison.
See `patches/GraphMol/ChemTransforms/001-compare-linker-bond-types.patch`.

### `MolZip.cpp`: failed zipping is dereferenced during coordinate generation

The lower-level `molzip()` explicitly returns a null pointer when fragment-on-
bond labels cannot be resolved. The decomposition overload unconditionally
calls `zippedMol->getNumAtoms()` when coordinate generation is enabled, turning
that recoverable failure into a null-pointer dereference. Include `zippedMol` in
the guard. See
`patches/GraphMol/ChemTransforms/002-guard-failed-coordinate-generation.patch`.

### `MolZip.cpp`: coordinate remapping dereferences a possibly missing atom

The coordinate-copy loop dereferences the result of `std::find_if()` without
checking it against `end()`. A missing or inconsistent internal zip-index
property therefore causes undefined behavior. The raw `Conformer` allocated
before that lookup also leaks if this lookup or any property access throws.
Check the lookup invariant and retain the conformer in a `unique_ptr` until
ownership is transferred. See
`patches/GraphMol/ChemTransforms/003-validate-coordinate-atom-lookup.patch`.

## Outdated or incorrect comments

### `MolZip.cpp`: `ZipBond::bond()` comment still claims bonds are always single

The implementation preserves a non-single bond type from either attachment,
so the preceding “only use single bonds” statement and its FIXME no longer
describe the code. Replace them with the actual selection rule. See
`patches/GraphMol/ChemTransforms/004-correct-zip-bond-type-comment.patch`.

## Redundant includes

### `MolFragmenter.cpp`: `trim.hpp` is already supplied by the string umbrella

`boost/algorithm/string.hpp` includes the trim algorithms used by this file, so
the following `boost/algorithm/string/trim.hpp` include is redundant. See
`patches/GraphMol/ChemTransforms/007-remove-redundant-trim-include.patch`.

## Simplifications

### `MolFragmenterJSONParser.cpp`: avoid scanning a string just to test emptiness

The parser calls `strlen(details_json)` only to decide whether the first byte is
the null terminator, and then default-constructs a stream before assigning its
contents. Test `*details_json` directly and construct the stream from the input;
this avoids an unnecessary linear scan and removes reliance on a transitive
declaration of `strlen`. See
`patches/GraphMol/ChemTransforms/008-simplify-empty-json-check.patch`.
