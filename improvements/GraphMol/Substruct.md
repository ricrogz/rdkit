# Substruct

## Correctness defects

### `SubstructUtils.cpp`: substitution scoring can divide by zero and overflow

`ScoreMatchesByDegreeOfCoreSubstitution` computes the sum of atom indices as
`na * (na + 1) / 2` in the unsigned type returned by `getNumAtoms()` and only
then casts to `double`. Large molecules can overflow before conversion, while
an empty molecule produces a zero denominator and a NaN tie-break score for an
empty match. Perform the multiplication in `double` and clamp the denominator
to at least one. See
`patches/GraphMol/Substruct/001-safe-substitution-score-denominator.patch`.
