# Chem

## Undefined variables

### `BuildFragmentCatalog.py`: cached-list validation refers to removed input text

The catalog-loading path tests `inD.count()` even though `inD` is never
defined. Validate the cached list against the configured molecule limit, which
is the same threshold used immediately before rescoring. See
`patches/rdkit/Chem/001-validate-cached-list-size.patch`.
