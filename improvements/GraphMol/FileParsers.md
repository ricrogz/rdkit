# FileParsers

## Undefined and implementation-defined behavior

### `MolFileWriter.cpp`: variadic format does not match `parityFlag`

`parityFlag` is an `unsigned int`, but both V2000 atom-line format strings pass
it to `%d`. A variadic call with a mismatched argument type has undefined
behavior. Use `%u` for this field in both the portable and MSVC branches. See
`patches/GraphMol/FileParsers/001-fix-parity-format-type.patch`.

### `MolFileStereochem.cpp`: unsigned negation is converted back to `int`

The neighbor count is a `size_t`, so multiplying it by `-1` is unsigned
arithmetic. Converting the resulting large unsigned value to the tuple's `int`
element is implementation-defined and can break the intended descending-count
sort. Convert the small bond-neighbor count to `int` before negating it. See
`patches/GraphMol/FileParsers/002-negate-neighbor-count-safely.patch`.

## Invalid object state

### `MaeWriter.cpp`: shared stream is checked after being moved from

The shared-pointer constructor moves `outStream` into `dp_ostream`, then checks
and dereferences the moved-from parameter. In ordinary implementations the
precondition therefore fails even for a valid stream (and a disabled
precondition would allow a null dereference). Validate the stored pointer
instead. See `patches/GraphMol/FileParsers/003-check-stored-mae-stream.patch`.

## Exception safety

### `TDTWriter.cpp`: destructor can terminate while closing the record stream

`TDTWriter::close()` writes the final `|` record delimiter. An ostream with an
exception mask can throw during that insertion; because destructors are
implicitly `noexcept`, allowing it to escape calls `std::terminate`. Keep
explicit `close()` error-reporting behavior, but suppress close failures during
destruction. See
`patches/GraphMol/FileParsers/004-make-tdt-destructor-nothrow.patch`.

## Uninitialized supplier state

### `MolSupplier.h`: default MAE supplier has indeterminate position and length

The public default constructor does not call `init()`, leaving `d_position` and
`d_length` indeterminate for methods such as `length()`, `moveTo()`, or
`getEOFHitOnRead()`. Initialize both counters in their declarations. See
`patches/GraphMol/FileParsers/005-initialize-mae-supplier-counters.patch`.
