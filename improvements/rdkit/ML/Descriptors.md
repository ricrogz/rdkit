# ML/Descriptors

## Malformed-expression handling

### `Parser.py`: an unmatched opening parenthesis is treated as a complete call

`ParseCompoundDescriptor()` advances `p` to `len(res)` when a known method call
has no matching closing parenthesis. Its `p <= len(res)` check nevertheless
rewrites that unterminated call, dropping its final character and producing a
misleading derived expression. Only rewrite a call after its parentheses have
balanced. See
`patches/rdkit/ML/Descriptors/001-require-balanced-parentheses.patch`.
