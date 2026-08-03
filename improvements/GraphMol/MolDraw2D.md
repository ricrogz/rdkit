# MolDraw2D

## Uninitialized state

### `DrawMol.h`: constructor-specific members have no safe defaults

The molecule-drawing constructor does not initialize `molHeight_` before its
setup helpers run, and the transform-only constructor does not initialize
`includeAnnotations_`. Current call paths generally assign or avoid these
members before reading them, but the class itself does not uphold that invariant
and static analysis detects the indeterminate state. Give both members benign
in-class defaults. See
`patches/GraphMol/MolDraw2D/001-initialize-drawmol-members.patch`.

## Invalid container access

### `AtomSymbol.cpp`: colon adjustment uses the unstripped string for bounds

`colonPos` is calculated after markup is stripped, but the right-neighbor test
uses the original `symbol_.size()`. A colon at the end of the visible symbol can
therefore pass the test merely because markup made the original string longer,
then read `rects_[colonPos + 1]` out of bounds. The invariant also incorrectly
allows `colonPos == rects_.size()` before indexing `rects_[colonPos]`. Validate
and test against `rects_.size()` itself. See
`patches/GraphMol/MolDraw2D/002-use-rectangle-bounds-for-colons.patch`.

## Infinite loops

### `AtomSymbol.cpp`: malformed markup without `>` never advances

When a custom atom label contains `<` with no closing `>`, `find('>')` returns
`npos`; adding one wraps to zero and reconstructs the same string, so the
`while (true)` loop never terminates. Stop colon adjustment when markup is
malformed. See
`patches/GraphMol/MolDraw2D/003-handle-unterminated-label-markup.patch`.

## Performance improvements

### `DrawShape.h`: wedge, wavy-line, and arc constructors copy point vectors twice

Four frequently used shape constructors accept their point vector by value and
then pass it to a base constructor that stores its own copy. Taking the input by
const reference removes the extra allocation and copy without changing
ownership or lifetime. See
`patches/GraphMol/MolDraw2D/004-pass-shape-points-by-reference.patch`.

## Redundant includes

`MolDraw2D.cpp` repeats `DrawMol.h`, and `DrawTextSVG.cpp` repeats
`MolDraw2DDetails.h`. See
`patches/GraphMol/MolDraw2D/005-remove-duplicate-includes.patch`.

## Redundant Python imports

`side_by_side_images.py` imports `glob` but uses `Path.glob()`, while
`update_hash_codes.py` imports `json` without using it. Remove both imports.
See `patches/GraphMol/MolDraw2D/006-remove-unused-python-imports.patch`.
