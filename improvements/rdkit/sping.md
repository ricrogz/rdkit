# sping

## Python 3 compatibility

### `stringformat.py`: the module cannot be imported on supported Python versions

The formatter imports `xmllib`, which was removed from Python 3, and later
imports `sping.PDF` as though `sping` were a top-level package. Consequently,
none of its public formatting functions can be used through `rdkit.sping`.
Use `html.parser.HTMLParser` with small tag/entity dispatch adapters and make
the PDF import package-relative. The patch also imports the colors used by the
module's callable demonstration routine. See
`patches/rdkit/sping/001-port-string-formatter-to-python3.patch`.
