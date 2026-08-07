# CIPLabeler correctness patch series

This directory contains 14 ordered `git format-patch` files for the C++
implementation under `Code/GraphMol/CIPLabeler`. They implement the
correctness and robustness fixes identified as C++-1 through C++-19 in
[`CIP_comparison.md`](../CIP_comparison.md). CIPLabeler calls themselves must
remain serialized. The series preserves master's `RDK_TEST_MULTITHREADED`
Ctrl-C regression: one thread runs the sole CIPLabeler call while a second
thread delivers SIGINT; it does not run overlapping labeler calls.

## Baseline and ordering

- RDKit base commit: `e3eb92687579`
- Expected tree after this series: `4768dee389b71b45fc0f98e8d6ce4214bf5c2f0c`
- Reference endpoint used to generate the files:
  `5324307d91ec6c33dc2883b1e1bddc1d1c14d3a0`
- Apply the files in the order recorded in [`series`](series).
- Apply this complete series before the performance series in
  `../cip_labeler_performance_patches`.

From a clean RDKit checkout at the base commit:

```bash
git am ../cip_labeler_bug_patches/*.patch
```

The locally created commit IDs may differ because `git am` supplies new
committer metadata; the resulting tree should match the tree ID above.

## Patch contents

| Patch | Report finding(s) | Change |
|---|---|---|
| `0001` | C++-2 | Eagerly constructs the composite `Rules` sorter instead of mutating a global `const` object on first use and verifies its lifetime sequentially. Existing master tests and their supporting includes are preserved unchanged. |
| `0002` | C++-4 | Returns Rule 6 decisions with magnitude 2 so downstream code recognizes pseudoasymmetric decisions. |
| `0003` | C++-13 | Makes long `PairList` inputs safe and avoids invalid full-width shifts. |
| `0004` | C++-10 | Validates every Rule 4b priority group instead of repeatedly inspecting group zero. |
| `0005` | C++-12 | Makes `maxRecursiveIterations` an exact shared budget across the complete serialized labeling call, including the preliminary pass, restores scoped state on return or exception, and documents the single-thread API contract. |
| `0006` | C++-3, C++-7, C++-17, C++-18 | Keeps unselected configurations available as auxiliary dependencies, writes only selected outputs, fixes `_CIPComputed` and stale-property handling, validates selections, and restores stable last-wins auxiliary-label behavior. |
| `0007` | C++-7, C++-8, C++-9, C++-11 | Uses XOR for atrop pseudo-descriptors, falls back to an unknown isotope's mass number, clears neighbor-order state, and rejects malformed configurations safely. |
| `0008` | Multiple | Adds regression coverage for isotope fallback, selective dependencies and masks, property lifecycle, and malformed stereo markers, including a structurally valid fixture for a CIS marker on a single bond. |
| `0009` | C++-1, C++-6, C++-15 | Computes one final value per Mancude resonance component, preserves the unreduced denominator as resonance metadata, and fully relaxes atoms outside the ring two-core. |
| `0010` | C++-16 | Falls back to the untouched input molecule after failed kekulization instead of caching a partially modified temporary. |
| `0011` | C++-5, C++-14, C++-19 | Widens visit distances, enforces the node cap at each insertion, and prevents `Digraph` from retaining a temporary `CIPMol`. |
| `0012` | C++-2 | Extends eager sorter construction to every standalone `SequenceRule`, removing the remaining lazy `const_cast`. |
| `0013` | C++-11 | Validates double-bond carrier indexes, identity, distinctness, and adjacency before dereferencing them. |
| `0014` | C++-17 | Rejects only set indexes outside the molecule while retaining compatibility with short masks and padded clear bits. |

Several fixes are split across patches so each reviewable concern remains
small. In particular, `0001` plus `0012` are the complete sorter-lifetime
fix, while `0006` plus `0014` are the complete selection-mask fix.

## Verification status

The user built and ran the original series. That run exposed an iteration
budget leaking into later standalone rule comparisons and a malformed test
fixture that failed inside `Bond::setStereo()` before reaching CIPLabeler.
Both corrections are folded into this refreshed series. The refreshed files
have not yet been built or tested; the user will perform that verification.

Suggested commands are:

```bash
cmake --build <build-dir> --target testCIPLabeler
ctest --test-dir <build-dir> -R '^testCIPLabeler$' --output-on-failure
```

The refreshed series was replayed onto a clean worktree. The resulting source
diff was checked for whitespace errors and compared with the corrected working
tree during generation. No build or test command was run while refreshing it.

Focused regression fixtures are still desirable for the both-ends-pseudo
atrop XOR path, overlapping auxiliary-label collisions, and the exact
`N`-versus-`N-1` iteration-budget boundary. The corresponding source changes
were reviewed statically but are not each exercised by a dedicated fixture in
this series.

The later V157894 pathological-runtime follow-up adds no correctness-series
patch and does not change this directory's ordering or endpoint. Its exact
symmetry gates retain the conservative corrected behavior whenever a proof is
ambiguous or exhausted; the implementation and focused coverage are recorded
as performance patch `0014` in `../cip_labeler_performance_patches` and in the
post-series appendix to `../CIP_comparison.md`.
