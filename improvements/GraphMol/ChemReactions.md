# ChemReactions

## Performance improvements and simplifications

### `DaylightParser.cpp`: reaction-prefix truncation allocates a new string

After locating the first whitespace, the parser assigns a prefix produced by
`substr()` back to the same string. `resize()` expresses the truncation
directly and avoids allocating and copying a temporary string. See
`patches/GraphMol/ChemReactions/001-resize-reaction-prefix.patch`.
