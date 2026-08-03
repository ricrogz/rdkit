# MarvinParse

## Uninitialized members

### `MarvinDefs.cpp`: public default constructors leave serialized fields indeterminate

The default `MarvinAtom`, `MarvinDataSgroup`, and `MarvinMulticenterSgroup`
constructors leave `isotope`, `x`/`y`, and `center` uninitialized respectively,
while their serialization and conversion methods read those fields. Current
parser paths usually fill them later, but a default-constructed object or an
exceptional/partial construction path can read indeterminate values or an
invalid pointer. Give each member a neutral default. See
`patches/GraphMol/MarvinParse/001-initialize-serialized-members.patch`.

## Redundant includes

`MarvinParser.cpp` repeats `MarvinParser.h`, and `MarvinWriter.cpp` repeats
`LocaleSwitcher.h`. See
`patches/GraphMol/MarvinParse/002-remove-duplicate-includes.patch`.
