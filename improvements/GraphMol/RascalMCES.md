# RascalMCES

## Invalid container access

### `RascalMCES.cpp`: timeout handling can call `front()` after `clear()`

When a timeout occurs with a clique larger than the saved best clique,
`checkTimeout()` clears `maxCliques` and then immediately evaluates
`maxCliques.front()` in the next condition. This is undefined behavior on an
empty vector precisely on the path meant to preserve the better partial result.
Accept an empty vector explicitly before accessing its first element. See
`patches/GraphMol/RascalMCES/001-avoid-front-after-clear.patch`.

## Incorrect copy semantics

### `RascalResult.cpp`: copy operations omit result options and scores

The copy constructor omits `d_ringMatchesRingOnly`, `d_maxFragSep`, and
`d_exactConnectionsMatch`. The copy-assignment operator additionally omits both
tier similarities and `d_ignoreBondOrders`, and retains old molecule pointers
when the source pointer is null. Copies can consequently produce different
SMARTS/results from their source. Copy every constructor field and implement
assignment through a complete temporary copy followed by move assignment. See
`patches/GraphMol/RascalMCES/002-copy-complete-result-state.patch`.
