# MolProcessing

## Correctness defects

### `MolProcessing.cpp`: trailing parse failures disappear in multithreaded mode

`mtWorker()` sizes its result from the greatest record ID present in `accum`,
but the supplier invokes the write callback only for successfully parsed
molecules. If the last record (or every record) fails to parse, its ID never
enters `accum`, so the returned vector is too short and no longer has one entry
per input record. Track `getLastRecordId()` while draining the supplier and use
the greatest consumed record ID for the result size. See
`patches/GraphMol/MolProcessing/001-preserve-trailing-parse-failures.patch`.

## Outdated or incorrect comments

### `MolProcessing.cpp`: API documentation describes an obsolete return type

The documentation says the function returns an `ExplicitBitVect,bitset` pair,
while the implementation and declaration return a vector of owning pointers,
using null entries for records that could not be read. The same block also has
the typo “whegn”. Update the contract description. See
`patches/GraphMol/MolProcessing/002-correct-fingerprint-api-documentation.patch`.
