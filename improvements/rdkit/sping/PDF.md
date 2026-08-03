# sping/PDF

## Ineffective assignments

### `pdfmetrics.py`: missing glyph widths are never filled

The AFM parser uses comparison (`==`) where it intends assignment, so zero
widths remain zero instead of inheriting the space width. Replace it with `=`.
See `patches/rdkit/sping/PDF/001-fill-missing-glyph-widths.patch`.
