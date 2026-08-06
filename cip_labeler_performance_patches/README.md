# CIPLabeler performance patch series

This directory contains thirteen ordered `git format-patch` files implementing
the low-risk performance improvements identified during the Java/C++
comparison. The changes preserve the corrected behavior established by the
bug series.

## Baseline and ordering

- Required parent tree: complete bug series,
  `105c21fb52168d1fd276623cc3d057b7991f5e74`
- Expected tree after this series: `50fcf284fa05b053f49a3172cf1d75297ff4099e`
- Reference endpoint used to generate the files:
  `1fe9e6e560bf8f21c6fff68e38d147a30d6d0128`
- Apply the files in the order recorded in [`series`](series).

After applying `../cip_labeler_bug_patches`:

```bash
git am ../cip_labeler_performance_patches/*.patch
```

The locally created commit IDs may differ because `git am` supplies new
committer metadata; the resulting tree should match the tree ID above.

## Patch contents

| Patch | Improvement |
|---|---|
| `0001` | Records a completed constitutional pass so unresolved configurations are not immediately ranked a second time with the same rules. |
| `0002` | Caches each source atom's mass in `CIPMol`, avoiding repeated periodic-table/isotope lookup for duplicate graph nodes. |
| `0003` | Avoids cloning and kekulizing molecules with no aromatic bonds while preserving the untouched-input fallback for failed aromatic kekulization. |
| `0004` | Replaces allocation-heavy linked-list traversal queues with reserved contiguous vectors. |
| `0005` | Reuses queue and edge scratch buffers in recursive rule traversals and removes intermediate node-list allocations. |
| `0006` | Replaces the per-node array of 32-bit visit distances with a compact path bitset plus parent ancestry for the uncommon distance lookup. |
| `0007` | Uses stable `deque` storage for graph nodes and edges, fixes edge reservation, indexes seen molecule atoms with a bitset, and evaluates the negative-resonance predicate lazily. |
| `0008` | Memoizes completed comparisons of exact ordered edge pairs while one sequence-rule recursion is active, eliminating repeated traversal of symmetric ring chains without retaining results across mutable CIP state. |
| `0009` | Shares bounded exact-comparison results across one configuration-label session and cuts equivalent symmetric continuations at focus-free molecular bridges. Auxiliary discovery still traverses every component containing a valid configuration focus, while Rule 4b/5 skip exhaustive traversal when no effective auxiliary descriptor exists. |
| `0010` | Accelerates the highly cyclic stereocage regression by batching target-guided auxiliary occurrence discovery, proving symmetric cyclic constitutional ligands equal through constrained molecular self-isomorphism, pruning descriptor-free Rule 3-6 branches with reroot-aware counts, and retaining bounded exact comparison and sort results across each immutable auxiliary-label distance sphere. |
| `0011` | Reuses direction-local comparisons and sorts across auxiliary reroots, keeps Rule 1a/1b/2 results across distance spheres, memoizes exact path-constrained automorphism queries, keeps Rule 4b/5 reference sorters stable, stores the common one-word visit path inline, and replaces the stable auxiliary sort with distance buckets. |
| `0012` | Generalizes the symmetry acceleration to arbitrary component sizes and ring topologies: persistent large-path visit state, component-local exact self-isomorphism, bounded explicit-field color refinement, exact linear-time distance-signature filtering, sparse dominance evidence, cumulative per-key/component search limits, logarithmically bounded prefilter use, and O(1) directed cyclic-core eligibility after one linear pass. It removes the earlier size, aromaticity, and bounded-automorphism-sample gates without adding molecule- or atom-number-specific cases. |
| `0013` | Accelerates auxiliary labeling in large symmetric polycyclic graphs by retaining the support of each exact constitutional-automorphism witness and avoiding unseen auxiliary configurations only when every other stereo annotation is fixed pointwise. It also batches residual target reachability, reroots the occurrence tree along its parent/LCA paths, stores terminal-leaf edges inline, stops comparisons at tied terminal nodes, and proves identical constitutional continuations across molecular bridges. These gates depend only on graph topology and CIP invariants, not fixture identities or atom numbering. |

The following larger or workload-dependent opportunities from the comparison
are intentionally not included: canonicalizing equivalent directed occurrence
states, subtree-result dynamic programming beyond path-only rerooting, hybrid
dense/sparse automorphism witnesses, arena allocation for large visit-state
records, and a broader lazy-digraph rewrite. Those need dedicated benchmarking
and a larger design review.

## Verification status

The user built and tested the original combined series. This performance
series has been regenerated on top of the corrected bug-series endpoint. The
refreshed combined series has not yet been built, tested, or benchmarked; the
user will perform that verification.

Suggested commands are:

```bash
cmake --build <build-dir> --target testCIPLabeler
ctest --test-dir <build-dir> -R '^testCIPLabeler$' --output-on-failure
```

Performance should be measured separately on representative ordinary,
aromatic, highly symmetric, and deep/ring-rich molecules. The refreshed patch
files were checked for whitespace errors, replayed onto a clean worktree, and
compared with the applied source tree (the local slow regression fixtures and
the user's node-cap adjustment are intentionally outside this patch series).
No build, test, or benchmark command was run while refreshing them.
