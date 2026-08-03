# MolChemicalFeatures

## Redundant includes

`FeatureParser.cpp` includes `boost/tokenizer.hpp` twice. See
`patches/GraphMol/MolChemicalFeatures/001-remove-duplicate-tokenizer-include.patch`.
