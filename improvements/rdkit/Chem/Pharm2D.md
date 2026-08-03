# Chem/Pharm2D

## Undefined variables and discarded results

### `LazyGenerator.py`: feature matches are overwritten and an undefined name is stored

The constructor correctly groups feature atom IDs by family, then replaces the
mapping with a list of invalid zero-argument `GetMolFeature()` calls and assigns
the undefined `pattMatches`. Preserve the grouped results in signature-family
order and store that list. See
`patches/rdkit/Chem/Pharm2D/001-build-lazy-feature-match-list.patch`.
