# Depictor

## Off-by-one errors

### `EmbeddedFrag.cpp`: neighbor-angle loop includes each neighbor paired with itself

The inner iterator is initialized with `auto nbi2 = nbi3++`; post-increment
assigns the old iterator, so each outer iteration starts with `(nbi1, nbi1)`.
Those zero/self angles are not valid neighbor pairs and add needless geometry
work. Start at the next neighbor instead. See
`patches/GraphMol/Depictor/001-skip-self-neighbor-angle.patch`.

## Incomplete copy operations

### `EmbeddedFrag.h`: assignment preserves the destination atom ID

`EmbeddedAtom::operator=` copies every geometric and depiction field except
`aid`, unlike the defaulted copy constructor. Assigning between atoms therefore
combines the source geometry with the destination identity. Copy `aid` as part
of the assignment. See
`patches/GraphMol/Depictor/002-copy-embedded-atom-id.patch`.

## Redundant includes

`RDDepictor.cpp` repeats both `Chirality.h` and `EmbeddedFrag.h`. See
`patches/GraphMol/Depictor/003-remove-duplicate-includes.patch`.
