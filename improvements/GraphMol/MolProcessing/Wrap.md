# MolProcessing/Wrap

## Redundant includes

### `rdMolProcessing.cpp`: duplicate `GeneralFileReader.h` dependency

`MolProcessing.h` already includes `GeneralFileReader.h` as part of exposing
`GeneralMolSupplier::SupplierOptions`, so the wrapper's second direct include is
redundant. Removing it reduces dependency noise without changing compilation.
See `patches/GraphMol/MolProcessing/Wrap/001-remove-redundant-general-file-reader-include.patch`.
