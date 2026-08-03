# GraphMol

## Undefined behavior

### `Kekulize.cpp`: dummy permutation counter shifts past its width

`QuestionEnumerator` represents every subset with an `unsigned int` and
computes `0x1u << d_questions.size()`. A fused aromatic system with at least as
many questionable dummy atoms as value bits in `unsigned int` causes an
out-of-range shift; the per-bit shifts have the same problem. On common
32-bit-`unsigned` targets this can make the fallback enumerator immediately
return no configurations once the list reaches 32 elements. Represent the
binary counter with `boost::dynamic_bitset`, which is already a dependency of
this translation unit and has no fixed-width shift boundary. See
`patches/GraphMol/Kekulize/001-width-independent-question-enumerator.patch`.

## Outdated or incorrect comments

### `Kekulize.cpp`: ring-ordering comment names the wrong wedge endpoint

The ring preordering code tests `wedgedAtoms`, which is populated with each
wedged bond's begin atom, but the adjacent comment says it favors the end atom.
Correct the comment so it describes the code. See
`patches/GraphMol/Kekulize/002-correct-wedge-endpoint-comment.patch`.

## Invalid container access

### `Renumber.cpp`: atom index equal to `numAtoms` passes validation

`renumberAtoms()` rejects an input index only when it is greater than the atom
count. Valid indices end at `numAtoms - 1`, so an index exactly equal to the
count writes past `revOrder`. The routine also accepts duplicate indices,
leaving some reverse-map entries at their default value and silently corrupting
the result. Initialize the reverse map with a sentinel and validate that every
index is both in range and unique. See
`patches/GraphMol/Renumber/001-validate-renumber-permutation.patch`.

## Null dereferences and invalid identity comparisons

### `Canon.cpp`: optional neighboring bonds are converted to references

`updateDoubleBondNeighbors()` passes both optional second-neighbor pointers as
dereferenced references even when one or both are null. This is undefined
behavior before the conflict helper can select a branch that does not use the
missing bond. In the both-conflicting branch, `std::make_pair` also copies each
bond, so address comparisons against the original first bond can never match
and the wrong direction may be removed. Preserve the optional bonds as
pointers and iterate over pointers when bond identity matters. See
`patches/GraphMol/Canon/001-preserve-optional-bond-pointers.patch`.

## Moved-from state

### `Canon.cpp`: stereogroup metadata is set after its atom list is moved

The canonical stereogroup loop moves `sgAtoms` into the new group and only then
tests the moved-from vector and accesses its first atom. In normal library
implementations the vector is empty, so `_stereoGroup` is never written. Set
the property before transferring the atom list. See
`patches/GraphMol/Canon/002-set-stereogroup-before-move.patch`.

## Signed/unsigned arithmetic

### `RWMol.cpp`: negative ring-stereo indices are offset in unsigned arithmetic

Combining a negative encoded atom index with the unsigned atom count promotes
the expression to unsigned; its later conversion back to `int` is
implementation-defined. Subtract the checked integer atom-count offset
directly from the already-negative value. See
`patches/GraphMol/RWMol/001-offset-negative-ring-index-as-int.patch`.

## Performance improvements

### `Atropisomers.cpp` and `Chirality.cpp`: stereogroups are copied during cleanup

Both cleanup passes iterate the molecule's stereogroups by value even though
the input group is read-only and a replacement group is built explicitly.
Each copy duplicates the group's atom and bond vectors. Iterate by const
reference to remove those allocations, which matters most for large enhanced-
stereo ensembles. See
`patches/GraphMol/003-avoid-stereogroup-cleanup-copies.patch`.

## Uninitialized pointers

### `AtomIterators.h`: default heteroatom iterator deletes an indeterminate pointer

The default constructor initializes only `_mol`, while the destructor always
deletes `_qA`. Destroying a default-constructed iterator therefore reads and
may free an indeterminate address. Give the query pointer the same null member
initializer used by `QueryAtomIterator_`. See
`patches/GraphMol/AtomIterators/001-initialize-heteroatom-query.patch`.

## Exception specifications

### `ROMol.h`: move operations promise `noexcept` while repairing owners can throw

Both move operations call ownership setters and graph/container operations that
are not declared non-throwing. If an invariant, allocation, or property move
throws, the unconditional `noexcept` specification terminates the process.
Remove the overly strong specification so callers can handle the failure. See
`patches/GraphMol/ROMol/001-remove-unsafe-move-noexcept.patch`.

## Redundant includes

`QueryOps.h` repeats `<functional>` and `new_canon.h` repeats `<cstring>`.
See `patches/GraphMol/004-remove-duplicate-standard-includes.patch`.
