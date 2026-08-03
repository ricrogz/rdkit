# RDGeneral

## Invalid iterator access

### `utils.h`: permutation search dereferences `end()` before testing it

`countSwapsToInterconvert()` tests `*probeIt2` before checking whether the
iterator has reached `probe.end()`. A mismatched sequence therefore reads past
the container instead of reaching the intended invariant failure. Reverse the
two short-circuit operands. See
`patches/RDGeneral/001-check-end-before-dereference.patch`.

## Redundant includes

`RDValue-taggedunion.h` includes `<cstdint>` twice. See
`patches/RDGeneral/002-remove-duplicate-cstdint-include.patch`.
