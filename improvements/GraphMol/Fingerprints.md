# Fingerprints

## Redundant includes

`RDKitFPGenerator.cpp` repeats `hash.hpp`; `FingerprintUtil.cpp` repeats
`RDKitBase.h` and `dynamic_bitset.hpp`. See
`patches/GraphMol/Fingerprints/001-remove-duplicate-includes.patch`.
