# ML/Data

## Undefined variables and off-by-one random labels

### `DataUtils.py`: non-shuffle randomization iterates an undefined collection

The alternative randomization branch acknowledges that `examples` is
undefined, then attempts `for _ in len(examples)`, which also passes an integer
where an iterable is required. Generate exactly `nPts` labels, and use
`randrange(nPossible)` so a class count of N yields labels `0..N-1`. See
`patches/rdkit/ML/Data/001-fix-randomized-result-generation.patch`.

## Outdated comments

### `DataUtils.py`: a fixed branch is still labeled as non-working

After correcting the non-shuffle branch, its inline comment would incorrectly
claim that it remains broken because `examples` is undefined. Remove the stale
comment. See `patches/rdkit/ML/Data/002-remove-stale-randomization-comment.patch`.
