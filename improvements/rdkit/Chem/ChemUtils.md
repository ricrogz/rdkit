# Chem/ChemUtils

## Python 3 compatibility

### `TemplateExpand.py`: command-line paths use the removed `file()` builtin

Both the SMILES input and output paths call Python 2's `file()` constructor.
Under Python 3 these branches raise `NameError`. Use `open()` instead. See
`patches/rdkit/Chem/ChemUtils/001-replace-file-builtin.patch`.
