# ForceFieldHelpers

## Portability and undefined behavior

### `MMFF/AtomTyper.cpp`: aromatic marker shifts a signed value into its sign bit

`1 << 31` shifts a signed `int` into a non-representable value, which is
undefined before it is assigned to `uint32_t`. Make the literal unsigned so the
bit mask is constructed with defined unsigned arithmetic. See
`patches/GraphMol/ForceFieldHelpers/001-use-unsigned-aromatic-bit.patch`.
