# Comparison of the Java `centres` implementation and RDKit's C++ `CIPLabeler`

## Executive summary

The C++ implementation remains recognizably and substantially a port of the Java implementation: both construct a lazily expanded ligand digraph, create duplicate nodes for rings and multiple bonds, calculate fractional atomic numbers for delocalized systems, apply the same production sequence-rule stack (1a, 1b, 2, 3, 4a, 4b, 4c, 5New, and 6), calculate auxiliary descriptors from the outside inward, and convert the resulting ligand order into stereodescriptors.

They are no longer functionally equivalent, however. The two projects have evolved independently since the original port:

- Java is a toolkit-independent library whose caller supplies explicit configuration objects. It now supports ordinary tetrahedral and double-bond configurations, atropisomerism, extended cis/trans systems, extended tetrahedral/cumulene systems, square-planar, trigonal-bipyramidal, octahedral, and square-pyramidal naming, including recent prime-locant and clockwise/anticlockwise suffix work.
- C++ is an RDKit-specific service. It discovers configurations from an `ROMol`, supports selective atom/bond labeling, exposes Python bindings, handles Ctrl-C and a recursive-comparison budget, records ranked neighbors, and has a preliminary constitutional-rule pass. It labels only tetrahedral atoms, ordinary stereogenic double bonds, and atropisomeric bonds, even though RDKit itself can represent allene and non-tetrahedral chiral tags.
- Java treats atomic-number-zero atoms as unresolved wildcards. C++ deliberately treats them as lowest-priority dummy atoms. That makes the implementations observably incompatible for dummy/query molecules.
- Most individual sequence rules are faithful ports, but current Java contains a Rule 6 correction that C++ is missing.

The highest-priority C++ findings are:

1. The charged-Mancude translation assigns a running prefix fraction to each atom instead of one final component fraction. This was reproduced and makes priorities atom-order dependent.
2. `const` global rule objects are modified through `const_cast`; that is undefined behavior even in one thread and is also a first-use data race between threads.
3. Selective labeling omits unselected configurations needed as auxiliary stereochemical dependencies, yet marks the whole molecule `_CIPComputed`. A targeted example loses a valid pseudoasymmetric label, and downstream RDKit fingerprints can then skip the required full calculation.
4. C++ Rule 6 returns `+/-1` instead of Java's corrected `+/-2`, so Rule-6 decisions are not reported as pseudoasymmetric decisions.
5. The C++ visit-distance array uses ordinary `char`. On the current platform it corrupts Rule 1b distances after depth 127 and the visited set at depth 256. A 260-atom acyclic-chain probe generated false ring and primary nodes.
6. `boost::rational` normalization discards a denominator that the graph-building code also uses as a boolean “is delocalized” marker. Even after fixing the prefix-fraction defect, resonant fractions that reduce to an integer will be expanded incorrectly.

These are not detected by the current C++ tests: the checked C++ executable passed all 6,623 assertions in 22 test cases. Specific missing regression tests are proposed at the end of this report.

## Scope and method

This review compares these source snapshots:

- C++: `rdkit/Code/GraphMol/CIPLabeler` at RDKit commit `e3eb92687579` (`master`).
- Java: `centres/core/src/main/java/com/simolecule/centres` at centres commit `d4b3cf07557b` (`develop`).

The comparison covered the public entry points, molecule adapters, configuration classes, digraph/node/edge implementation, Mancude processing, all sequence rules, sorting and pairing, state written back to the molecule, exception/limit behavior, and relevant RDKit consumers of `_CIPComputed`.

Findings described as **confirmed** follow directly from control/data flow or were reproduced against the built C++ library. **Likely** findings are direct semantic mismatches with a later Java correction but do not have a new C++ end-to-end regression in this workspace. **Latent** findings require unusually deep, malformed, or otherwise uncommon input but have a concrete unsafe code path.

## 1. Architecture and public API

### 1.1 Molecule abstraction

Java's `BaseMol<A,B>` is an abstract toolkit adapter. It supplies atom/bond lookup, topology, ring membership, atomic number, isotope mass number, charge, hydrogen count, integer bond order, property access, and diagnostic graph dumping (`BaseMol.java:40-125`). The core library neither owns nor assumes a concrete molecule implementation.

C++'s `CIPMol` is a concrete, non-owning wrapper around an RDKit `ROMol` (`CIPMol.h:66-101`). It directly uses `Atom`, `Bond`, RDKit ring information, RDKit's periodic table, and RDKit stereochemical flags. It additionally:

- caches bond pointers by index (`CIPMol.cpp:19-22,42`);
- lazily clones and kekulizes the entire molecule to obtain integer bond orders (`CIPMol.cpp:65-80`);
- maps zero, hydrogen, and dative bond types to order zero and supports integer orders through six (`CIPMol.cpp:82-109`); and
- lazily caches fractional atomic numbers for the entire molecule (`CIPMol.cpp:24-29`).

The Java adapter is responsible for supplying chemically appropriate integer bond orders. The C++ core takes responsibility for kekulization itself. This makes C++ easier to call in RDKit but less reusable and gives it additional state, performance, and failed-kekulization edge cases.

### 1.2 How configurations enter the algorithm

Java's entry point is:

```java
label(BaseMol<A,B> mol, List<Configuration<A,B>> configs)
```

The caller/toolkit adapter decides which stereogenic units exist and supplies their foci, carriers, and configuration parity (`Labeller.java:42`). This is naturally extensible and also means correctness depends on the caller supplying a complete configuration list.

C++'s entry points accept an `ROMol`, optionally with atom and bond bitsets (`CIPLabeler.h:58-86`). `findConfigs()` examines selected atoms for `CHI_TETRAHEDRAL_CW/CCW` and selected bonds for E/Z, cis/trans, or atropisomer flags (`CIPLabeler.cpp:49-100`). The C++ configuration constructors derive carriers from RDKit neighbor order, stereo atoms, and atropisomer helper data.

Consequences include:

- Java core can label a configuration representation that a toolkit adapter knows how to construct, without changing `Labeller`.
- C++ silently ignores every RDKit configuration type not enumerated in `findConfigs()`.
- C++ selection is presented as controlling which centers are *written*, but currently also removes those centers from the auxiliary-descriptor universe. This is a correctness issue discussed in finding C++-3.

### 1.3 Rule lifetime and call scheduling

Java constructs constitutional and full rule stacks on every call because its rules retain the `BaseMol` adapter (`Labeller.java:42-60`). It then processes each configuration in list order: constitutional label, auxiliary labeling if needed, then full label. Every configuration receives a fresh digraph (`Labeller.java:62-75`).

C++ rules are molecule-independent because node objects cache molecule-specific values. Two rule stacks are global (`CIPLabeler.cpp:42-47`). Every configuration owns a digraph rooted at construction (`Configuration.cpp:73-78`). C++ runs all configurations through a 2,000-comparison constitutional pass first, then revisits unresolved configurations with the requested/shared budget and auxiliary rules (`CIPLabeler.cpp:172-227`).

The preliminary pass often avoids expensive auxiliary work and ensures both ends of unresolved bond configurations have been reached. It can also repeat comparisons for difficult centers and, as implemented, is not charged to the caller's advertised maximum iteration count.

This scheduling also fixes a subtle Java limitation for bond centers. Java's primary `Sp2Bond.label()` returns immediately if the first end is nonunique (`Sp2Bond.java:117-119`), so configurations reachable only through the second end may not be present when auxiliary descriptors are collected. C++ deliberately permits the constitutional pass to expand/sort the second end before returning unresolved (`Sp2Bond.cpp:91,113-129`); atrop handling follows the same pattern. This makes the C++ auxiliary search more complete.

### 1.4 Errors, interruption, and partial results

Java catches every `Exception` around each configuration, prints its message to standard error, and continues (`Labeller.java:62-78`). A failed center is therefore normally a missing label rather than a failed API call. Java has the same 100,000-node digraph guard as C++ but has no recursive-comparison budget or interrupt handling.

C++ normally propagates `TooManyNodesException`, `MaxIterationsExceeded`, invariants, and ordinary runtime errors. Only `MaxIterationsExceeded` in the preliminary pass is swallowed. Ctrl-C is caught, logged, and returned as an interrupted call; the Python wrapper converts it to `KeyboardInterrupt` (`CIPLabeler.cpp:241-249`; `Wrap/rdCIPLabeler.cpp:44-47`). A failure or interrupt can leave a partially updated molecule, though `_CIPComputed` is deliberately cleared before work starts.

## 2. Supported stereochemical configurations

| Configuration | Java core | C++ `CIPLabeler` | Important difference |
|---|---:|---:|---|
| Tetrahedral | Yes | Yes | Java receives four carriers and parity; C++ derives them from RDKit neighbor/chiral-tag order and has special handling for implicit-H and trigonal-pyramidal cases. |
| Ordinary double bond | Yes | Yes | Both produce `E/Z` or pseudo `seqCis/seqTrans`; C++ writes pseudo forms as strings `z/e`, records ranked anchors, and normalizes RDKit stereo atoms/flags. |
| Atropisomeric axis | Yes | Yes | Java receives two foci and four carriers and uses parity. C++ extracts RDKit atrop markers/anchors and uses the atrop stereo enum, toggling a local copy to account for ranked-anchor swaps. Their both-ends-pseudo logic differs; see C++-8. |
| Extended cis/trans | Yes | No | Java traverses the intervening cumulene chain and stores the result on a supplied bond (`ExtendedCisTrans.java`). |
| Extended tetrahedral / cumulene axial | Yes | No | Java supports its documented three-focus representation and longer cumulene traversal, returning `M/P/m/p` on the middle focus (`ExtendedTetrahedral.java`). RDKit `CHI_ALLENE` is not discovered by C++. |
| Square planar | Yes | No | Java emits `SP_4` plus a `conf.index` such as `SP-4-...`; C++ has an enum value but no configuration/discovery implementation. |
| Trigonal bipyramidal | Yes | No | Java emits `TBPY_5`, prime locants, and optional `-C/-A` winding. RDKit's corresponding chiral tag is ignored by C++. |
| Octahedral | Yes | No | Java emits `OC_6`, prime locants, and optional winding. RDKit's corresponding chiral tag is ignored by C++. |
| Square pyramidal (degenerate OC-6) | Yes | No | Current Java's octahedral class emits `SPY-5-...`, including winding; C++ has no equivalent. |

C++ `Descriptor` still declares `SP_4`, `TBPY_5`, and `OC_6` (`Descriptor.h:49-51`), but the atom/bond setters reject these for the implemented configurations and no code can discover or calculate them. Their presence is therefore not functional support.

Java's current atrop implementation is less complete than the table's simple “Yes” suggests. Its primary method constructs two temporary bond-aware digraphs (`Atropisomeric.java:90-112`), rather than populating the unrooted configuration digraph installed by `Labeller`, and it does not override `Configuration.label(node,digraph,rules)`. It therefore cannot contribute an atrop descriptor as an auxiliary configuration. An atrop primary center that needs auxiliary descriptors also cannot resolve: with another configuration present, `labelAux()` calls `getNodes()` on the unrooted graph and attempts to enqueue a null root (`Labeller.java:87-104`; `Digraph.java:113-128`), which throws; without another configuration, the full-rule retry merely creates fresh graphs with no auxiliary annotations. C++ implements both primary and node/digraph labeling in one rerootable graph and can use atrop centers in auxiliary calculations (`AtropisomerBond.cpp:82-186`). This is a material C++ functional extension.

The C++ API header's note that only chiral atoms and double bonds are labeled (`CIPLabeler.h:60-63`) is now stale because atropisomeric bonds are also supported.

## 3. Labels and molecule state

Java delegates storage through `BaseMol.setAtomProp/setBondProp`. Configurations store a `Descriptor` under `"cip.label"`. The coordination-geometry configurations (square planar, trigonal bipyramidal, and octahedral, including its square-pyramidal path) additionally store a detailed configuration index under `"conf.index"`, for example `SP-4-...`, `TBPY-5-...`, `OC-6-...`, or `SPY-5-...` (`BaseMol.java:44-45`; the corresponding `setPrimaryLabel()` methods).

C++ writes RDKit string properties:

- `_CIPCode` on the relevant atom or bond;
- `_CIPNeighborOrder` containing ranked neighboring atom indices; and
- molecule-level `_CIPComputed` on successful return.

`Descriptor::seqTrans` and `Descriptor::seqCis` stringify to lower-case `"e"` and `"z"`, while Java stores the enum values `seqTrans` and `seqCis` (`Descriptor.h:69-73`; `Descriptor.java:53-58`). Ordinary `E/Z`, tetrahedral `R/S/r/s`, and axial `M/P/m/p` spellings otherwise agree.

For double bonds C++ also calls `setStereoAtoms()` and `setStereo()` when committing a label (`Sp2Bond.cpp:37-49`). Because discovery maps `STEREOE/Z` to `STEREOTRANS/CIS`, labeling may normalize the RDKit stereo representation as a side effect. Java core changes only the property requested through its adapter.

Java does not clear an old primary property before recalculation. C++ clears `_CIPCode` on every *currently discovered and selected* configuration during its preliminary pass (`CIPLabeler.cpp:185-187`). That is better for a center that remains configured but becomes non-stereogenic; it does not clear `_CIPNeighborOrder`, and it cannot clear a center that is no longer discovered because its chiral/stereo marker was removed. Those state defects are covered in C++-7.

Java has `Stats`/rule-index instrumentation. C++ omits it but exposes Python bindings and records ranked-neighbor metadata.

## 4. Digraph and resonance-model comparison

### 4.1 Common model

Both implementations build the same basic lazy directed acyclic representation:

- a graph is initially rooted at a stereochemical focus;
- expanding an ordinary node creates explicit neighbors, virtual duplicate nodes for multiple bonds, duplicate nodes at ring closures, and implicit-hydrogen nodes;
- terminal duplicate/H nodes do not expand;
- changing the root reverses incoming edges along the required path rather than rebuilding the graph;
- each primary child receives a copy of an atom-indexed path/visited array; and
- graph generation is intended to stop around 100,000 nodes.

The relevant implementations track closely in `Digraph.java:71-247` / `Node.java:75-223` and `Digraph.cpp:37-204` / `Node.cpp:21-165`.

### 4.2 Representation differences

- Java uses garbage-collected node/edge objects and an unsigned 16-bit Java `char[]` path array. C++ stores nodes and edges in `std::list` for pointer stability and uses `std::vector<char>`.
- Java permits an uninitialized digraph followed by `init()`. C++ requires a root in the constructor and explicitly retains original/current roots.
- Atrop mode is represented differently. Java creates two bond-aware graphs so root multiple bonds are duplicated. C++ has one graph with an `atropisomerMode` flag and reroots across the axis (`Digraph.cpp:56-70,160-180`).
- C++ adds `seenAtom()` and skips an auxiliary configuration when none of its foci was reached during prior expansion (`CIPLabeler.cpp:117-121`). Java directly calls `getNodes()`, which can expand while searching. This is a useful C++ lazy-work optimization.
- Java sorting mutates the actual edge list returned by `Node.getEdges()`. C++ generally copies an edge vector before sorting, preserving graph insertion order at the cost of repeated allocations/copies.
- C++ caches atomic mass in each node; Java resolves isotope data during Rule 2 comparisons.

### 4.3 Mancude/fractional atomic numbers

The type seeding, two-core-like relaxation, resonant-component discovery, and ordinary fractional-number calculation are close ports (`Mancude.java` and `Mancude.cpp`). Java stores mutable `Fraction` objects, later casts their numerator/denominator into node `short` fields, and compares fractions through `double`. C++ uses exact, normalized `boost::rational<int>` values.

Exact rational comparison and full-width `int` storage are C++ improvements. However, translating Java reference semantics into C++ value semantics introduced C++-1, and normalization introduced C++-6 because the denominator is also used as nonnumeric metadata.

## 5. Sequence-rule comparison

The production compositions agree:

- constitutional stack: Rule 1a, Rule 1b, Rule 2;
- full stack: Rules 1a, 1b, 2, 3, 4a, 4b, 4c, 5New, and 6.

Each composite rule is applied exhaustively over the ligand graph before the next rule is tried.

Java's top-level `Rules.getSorter()` returns a fresh flat sorter containing the individual rules (`Rules.java:80-82`). C++ lazily caches a sorter containing the composite `Rules` object as one comparator (`Rules.h:48-53`); that comparator then iterates its children exhaustively. Production deep ranking is normally equivalent, but the exposed `getRules()` shape differs, and this cache design is the source of C++-2.

| Rule/component | Functional comparison |
|---|---|
| Rule 1a | Both compare atomic number, including fractional numbers for relevant duplicate nodes. C++ compares exact rationals. Java treats numerator zero as an unresolved wildcard and propagates a sentinel; C++ treats zero as an ordinary lowest value. |
| Rule 1b | Equivalent. Both retain `IUPAC_2013=false`, distinguish ring duplicates from ordinary nodes, and compare duplicate distances under the same conditions. C++'s smaller path-distance type makes its inputs unsafe at depth, however. |
| Rule 2 | Same intent, different source. Java uses its BODR isotope table on demand and falls back to the supplied mass number. C++ caches RDKit average/exact mass in nodes. C++ correctly gives duplicate nodes mass zero; Java can accidentally map a duplicate's mass number zero back to an element's average weight when compared to a labeled isotope. Conversely, C++ assigns unknown isotopes mass zero instead of Java's mass-number fallback (C++-9). |
| Rule 3 | Equivalent: `Z` outranks `E`, comparing only the end-node auxiliary descriptor (which can semantically describe a bond-centered configuration). |
| Rule 4a | Equivalent exact grouping: `{R,S,M,P,seqCis,seqTrans}` rank 2; `{r,s,m,p,E,Z}` rank 1; unset/unknown/`ns` rank 0. Bond labels are compared before node labels. Java `null` corresponds to C++ `Descriptor::NONE`. |
| Rule 4b | Close port of the first-descriptor-sphere, like/unlike-pair-list algorithm. C++ retains an unused descriptor traversal and has slightly stricter invariant checks. Both contain the same group-count indexing typo (C++-10). |
| Rule 4c | Equivalent: `m/r` rank above `p/s`, with bond labels before node labels. |
| Legacy Rule 5 | Direct comparison behavior agrees and it is not used by the production stack. Java deliberately excludes legacy Rule 5 from a branch sorter; C++ includes every accumulated rule. That only affects custom/legacy compositions. |
| Rule5New | Production behavior is equivalent: non-root comparison uses a temporary R/S reference, root comparison constructs R- and S-reference pair lists, and `+/-2` marks a pseudoasymmetric decision. |
| Rule 6 | Priority choice agrees, but comparison magnitude differs: Java's corrected implementation returns `+/-2`; C++ returns `+/-1`. This changes pseudoasymmetric classification (C++-4). |
| Sort/Priority | Both use stable insertion sorting, increment a counter for each decision with magnitude greater than one, and set the final pseudoasymmetric flag only when that count is exactly one. Java additionally carries wildcard state and the highest rule index for statistics. C++ adds recursive-budget and Ctrl-C checks. |
| PairList | Descriptor normalization and lexicographic like/unlike comparison are equivalent. C++ uses a 64-bit pairing field instead of Java's 32-bit field and omits unused Java append/clear APIs. C++ has a long-list shift UB described in C++-13. |

### Wildcard/dummy policy is deliberately incompatible

Current Java Rule 1a returns `COMP_TO_WILDCARD` when either atomic-number numerator is zero (`Rule1a.java:44-54`). `SequenceRule` propagates the sentinel; `Sort` records it and forces nonuniqueness. Configurations then return `Unknown` either by explicitly checking `wasWildcardFound()` or through the nonunique priority result. This was added to Java in commit `31ba6c4` so a center is undefined if an unspecified atom is reached.

C++ has no wildcard channel. It ranks atomic number zero below real elements (`Rule1a.cpp:20-24`). RDKit tests intentionally require one dummy ligand to participate as a lowest-priority substituent and two dummy ligands to tie (`catch_tests.cpp:554-606`). This behavior is reasonable for a concrete RDKit dummy atom, but query wildcards also commonly have atomic number zero. On a query molecule C++ can therefore emit a concrete CIP label based on an element identity that the query did not specify. A robust RDKit-specific policy would distinguish a concrete dummy attachment from a query wildcard rather than using atomic number alone.

## 6. Potential and confirmed C++ defects

### C++-1 — Charged resonant components receive order-dependent prefix fractions

**Severity:** high. **Confidence:** confirmed by source and runtime probe.

Java creates one mutable `Fraction` for each negatively charged resonant component, assigns the same object reference to every member, and then accumulates into it. All members consequently observe the final component fraction (`Mancude.java:222-240`).

C++ initializes a running numerator and denominator, assigns their *current* value by value to an atom, and only then incorporates that atom (`Mancude.cpp:239-267`). Different atoms in one delocalized component receive different prefixes. For `[CH-]1C=CC=C1`, the current C++ library produced:

```text
atom 0: 0/1
atom 1: 0/1
atom 2: 3/1
atom 3: 4/1
atom 4: 9/2
```

Java semantics give every atom `24/5`. The result is atom-index dependent and violates the algorithm's purpose of giving resonance-equivalent atoms equal priority. It also changes graph topology: `Digraph` checks whether the fractional denominator is greater than one to decide how negative-charge duplicate nodes are made (`Digraph.cpp:163-171,185-193`). In the example, only a subset is recognized as fractional.

**Recommended fix:** collect component members and calculate the final numerator/denominator first; assign the same final value to every member in a second pass. Test equality and CIP-label invariance before/after atom renumbering.

### C++-2 — Global `const Rules` objects are mutated and race during first use

**Severity:** high. **Confidence:** confirmed C++ language violation; concurrent failure is timing dependent.

`constitutional_rules` and `all_rules` are objects originally declared `const` (`CIPLabeler.cpp:42-47`). Their parent `dp_sorter` is initially null. `Rules::getSorter()` casts away constness and installs a `unique_ptr` lazily (`Rules.h:48-52`; `SequenceRule.cpp:45-49,134`).

Modifying an object that was originally defined `const` through `const_cast` is undefined behavior even in a single thread. Two threads making the first calls on separate molecules additionally race on the same `unique_ptr`; one can replace/delete a sorter while the other is returning or dereferencing it. The thread-local comparison budget does not protect this global cache.

Java creates independent rule/sorter objects per call and has no analogous shared cache.

**Recommended fix:** eagerly construct the parent sorter while each `Rules` object is still being constructed, then expose a truly immutable rule graph. An alternative is a correctly declared mutable cache guarded by `std::call_once`, but eager construction is simpler and removes the first-call branch as well.

### C++-3 — Selective labeling drops required auxiliary configurations and falsely marks the whole molecule computed

**Severity:** high. **Confidence:** confirmed by source and runtime probe.

The bitset overload filters configuration discovery itself (`CIPLabeler.cpp:49-100`). `labelAux()` can only derive auxiliary descriptors from the resulting filtered vector (`CIPLabeler.cpp:103-169`). A selected pseudoasymmetric center that depends on an unselected stereocenter therefore cannot be resolved, even though the API describes the bitsets as selecting which atoms/bonds “will be labeled.” Existing `_CIPCode` properties do not help because the sequence-rule calculation uses digraph auxiliary labels, not arbitrary primary properties already on the molecule.

This was reproduced with the existing para-stereochemical example at `catch_tests.cpp:460-476`:

```text
SMILES: C\C=C/[C@@H](\C=C\O)[C@H](C)[C@H](\C=C/C)\C=C\O

full assignment:    atom 3=R, atom 7=r, atom 9=S, _CIPComputed=1
only atom 7 chosen: atom 3=-, atom 7=-, atom 9=-, _CIPComputed=1
```

Atom 7 loses its valid `r` descriptor because atoms 3 and 9 were not built as auxiliary configurations.

Regardless of whether the call selected all, one, or even zero centers, successful return unconditionally sets molecule-wide `_CIPComputed` (`CIPLabeler.cpp:231-252`). RDKit's Morgan bond and atom invariants and fingerprint utility test only for this property's presence before deciding to skip a full call (`MorganGenerator.cpp:163-165,305-309`; `FingerprintUtil.cpp:78-82`). A partial call can consequently cause later chirality-sensitive fingerprints to omit uncomputed labels.

**Recommended fix:** discover all configurations for dependency calculation, carry a separate output mask controlling which primary labels are committed, and set `_CIPComputed` only after a full-molecule calculation. A partial-computation marker would need distinct semantics if callers need to cache subsets.

### C++-4 — Rule 6 does not mark a pseudoasymmetric decision

**Severity:** high for affected stereodescriptors. **Confidence:** high; direct mismatch with an upstream correction.

Current Java returns `+2/-2` when the selected Rule-6 reference breaks a tie and explains that reaching Rule 6 after Rule 5 must be treated as a pseudoasymmetric decision (`Rule6.java:50-65`). Java commit `def7380a` specifically changed this from `+1/-1` and added regression data.

C++ still returns `+1/-1` (`Rule6.cpp:21-34`). C++ `Sort` counts a pseudoasymmetric decision only when `abs(cmp) > 1` (`Sort.cpp:27-51`), and configurations use that flag to choose upper- versus lower-case descriptors. A decisive Rule-6 comparison therefore never contributes to lower-case/pseudoasymmetric classification in C++.

**Recommended fix:** port the Java `+2/-2` correction and its regression structures. Confirm all tetrahedral, double-bond, and axial descriptor casing because they share `Priority`.

### C++-5 — `char` path distances corrupt deep acyclic graphs and Rule 1b

**Severity:** high impact, low frequency. **Confidence:** confirmed by runtime probe on the target platform.

Every C++ node stores its path depths in `std::vector<char>` (`Node.h:131`). A child writes `d_dist + 1` through a `char` cast, the visited check treats zero as unvisited, and ring duplicates read that `char` back as their distance (`Node.cpp:21-23,110,112-116`). Rule 1b then compares those distances (`Rule1b.cpp:24-27`).

On this platform ordinary `char` is signed. Depths above 127 become negative; depth 256 becomes zero and makes an already visited atom appear new. A synthetic acyclic chain of 260 carbons with implicit hydrogens disabled produced 262 graph nodes instead of 260, a false ring duplicate for atom 254 at distance `-1`, and two nonduplicate occurrences of atom 255.

Java uses an unsigned 16-bit Java `char`, so its analogous wrap occurs only at 65,536. Both are below/around the 100,000-node guard, but the C++ threshold is practical for polymers and generated graphs.

**Recommended fix:** use at least `std::uint32_t`, or separate an atom-visited bitset from a full-width stored path depth. Add a graph-only long-chain regression.

### C++-6 — Normalizing the fraction erases the “resonant average” marker

**Severity:** medium to high for affected delocalized systems. **Confidence:** confirmed data-flow defect.

Java's `Fraction` retains an unreduced denominator. C++ uses normalized `boost::rational<int>` (`Mancude.h:51`). Both implementations then use `denominator > 1` not only for arithmetic but also as an “averaged Mancude value” signal. `Node::newTerminalChild()` checks it for every multiple-bond duplicate to choose the source fractional atomic number instead of the target atom's ordinary atomic number (`Node.cpp:25-35`). The negative-charge branches additionally use it to select the number/type of duplicate nodes (`Digraph.cpp:163-171,185-193`).

For example, Java's final raw fraction for `[CH-]1C=C1` is `12/3`. C++ normalizes the same numeric value to `4/1`; the denominator test then incorrectly says it is not a resonant average. Fixing C++-1 without addressing this issue would still build the wrong duplicate-node model whenever the average reduces to an integer. Neutral Mancude systems can receive the wrong duplicate-node priority value as well; charged systems can additionally receive the wrong duplicate topology/count.

**Recommended fix:** keep an explicit `isMancudeAveraged`/resonant-component flag (or raw component denominator) separate from the normalized rational used for numerical Rule 1a comparison.

### C++-7 — Recalculation can leave stale `_CIPCode` and `_CIPNeighborOrder` state

**Severity:** medium. **Confidence:** confirmed by source; neighbor-order case reproduced.

There are two related paths:

1. `resetPrimaryLabel()` clears only `_CIPCode`, while `setPrimaryLabel()` writes both `_CIPCode` and `_CIPNeighborOrder` in all three C++ configuration classes (`Tetrahedral.cpp:46-82`; `Sp2Bond.cpp:37-75`; `AtropisomerBond.cpp:45-80`). Relabeling `C[C@H](F)Cl`, changing F to Cl, and relabeling correctly removes the now-nonunique `_CIPCode` but leaves the old `_CIPNeighborOrder`.
2. Full assignment resets only configurations returned by current discovery. If a caller removes an atom chiral tag or bond stereo flag and calls `assignCIPLabels()` again, that former center is absent from `findConfigs()` and retains its old `_CIPCode`. The call then sets `_CIPComputed`, certifying stale data.

Java core also lacks general stale-property cleanup, but C++ explicitly promises and consumes molecule-level computed state, making internal consistency more important.

**Recommended fix:** clear all properties owned by this calculation together. On a full call, clear `_CIPCode` and `_CIPNeighborOrder` from all relevant atoms/bonds before discovery; on a selective call, clear both for selected output targets without claiming global completion.

### C++-8 — Atropisomer pseudo-descriptor uses OR instead of XOR

**Severity:** medium/high for the affected axial symmetry class. **Confidence:** likely correctness bug with direct Java history.

Java returns lower-case `m/p` only when exactly one end of the axis was distinguished pseudoasymmetrically (`priority1 != priority2`, `Atropisomeric.java:132-143`). Sp2 and extended-tetrahedral Java code use the same exclusive-or rule. Java commit `ecc22be8` is explicitly titled “Only m, p, seqCis, or, seqTrans if either and not both ends were sorted with Rule 5+” and changed atrop handling from OR to XOR.

C++ uses logical OR (`AtropisomerBond.cpp:172-183`), while its own `Sp2Bond` correctly uses inequality/XOR (`Sp2Bond.cpp:161-171`). When both ends are pseudoasymmetrically distinguished, Java/design semantics produce upper-case `M/P`, but C++ produces lower-case `m/p`.

**Recommended fix:** use XOR/inequality and add a both-axis-ends-pseudo regression.

### C++-9 — Unknown isotope mass is cached as zero

**Severity:** medium. **Confidence:** confirmed source behavior.

Java Rule 2 uses the exact isotope table when available and otherwise falls back to the integer mass number (`Rule2.java:61-74`). C++ `Node` calls `PeriodicTable::getMassForIsotope()` directly (`Node.cpp:50-57`), whose documented behavior is to return `0.0` for an unknown isotope (`PeriodicTable.h:239-250`).

An explicitly supplied but unlisted high-mass isotope can consequently rank below known isotopes or tie a duplicate in C++, instead of approximately ranking by its mass number. RDKit's general `Atom::getMass()` already has fallback behavior that can guide the fix.

**Recommended fix:** when the isotope lookup is unavailable/zero, use the supplied isotope mass number as Java does. Add known/unknown isotope ordering tests.

### C++-10 — Rule4b validates `tmp[0]` repeatedly and can index out of bounds

**Severity:** medium but likely rare. **Confidence:** confirmed typo, inherited from Java.

In `Rule4b::getNextLevel()`, the loop intended to verify that every equivalent node has the same number of priority groups reads `tmp[0].size()` on every iteration (`Rule4b.cpp:178-187`). It should read `tmp[i].size()`. If the sizes differ, the guard never detects it; the later `aTmp[i]` can be out of range or extra groups can be ignored (`Rule4b.cpp:189-196`).

Current Java has the same indexing typo (`Rule4b.java:133-147`), so this is an inherited defect rather than a port divergence.

### C++-11 — Malformed configured centers can throw or invoke undefined behavior

**Severity:** medium robustness issue. **Confidence:** concrete paths; requires malformed/programmatically mutated input.

- `AtropisomerBond` returns early when RDKit's helper rejects a bond marker, leaving its carrier vector empty (`AtropisomerBond.cpp:30-34`). `findConfigs()` nevertheless retains the object, and later labeling indexes carriers and sorted edges without first validating sizes (`AtropisomerBond.cpp:115-128,141,162-169`). A bogus programmatically set atrop flag can therefore cause out-of-range access/undefined behavior. Discovery should skip invalid configurations or construction should fail safely.
- Discovery does not verify that a cis/trans/E/Z marker is on a valid stereogenic double bond before constructing `Sp2Bond`. Its constructor requires `findStereoAtoms()` to return exactly two anchors, and labeling later indexes the first edge on each side (`Sp2Bond.cpp:21-34,140-158`). A programmatically marked single or otherwise malformed bond can abort on an invariant or unsafe assumption. Validate bond type, carrier count, and side-edge count before construction/indexing.
- `Tetrahedral::label()` returns `Descriptor::ns` when fewer than three edges exist (`Tetrahedral.cpp:100-109`). The driver treats every value other than `UNKNOWN` as assignable, but `Tetrahedral::setPrimaryLabel()` rejects `ns` and throws (`Tetrahedral.cpp:46-73`; `CIPLabeler.cpp:191-194`). Java stores `ns` through its unrestricted adapter setter. C++ should return `UNKNOWN` for a malformed tetrahedral candidate or handle `ns` deliberately.

### C++-12 — `maxRecursiveIterations` is neither exact nor a full-call maximum

**Severity:** low/medium API correctness. **Confidence:** confirmed control flow.

The helper decrements before checking, so a budget of `N` permits only `N-1` `recursiveCompare()` entries; `N=1` permits none (`CIPLabeler.cpp:267-268`). More importantly, every configuration first receives an independent 2,000-call constitutional pass before the caller's budget is installed (`CIPLabeler.cpp:180-205`). A request such as `maxRecursiveIterations=1` can therefore execute roughly 2,000 comparisons per center before enforcing the alleged maximum.

The later budget is shared across all unresolved configurations, whereas the preliminary allowance is reset per configuration. This is useful defensive scheduling, but it does not match the header's “maximum number of iterations” description.

**Recommended fix:** define whether the option is a per-center or per-call budget, charge the preliminary pass to it, and use post-decrement/explicit accounting so `N` means `N`.

### C++-13 — `PairList` can perform an invalid shift on long descriptor paths

**Severity:** low likelihood, potentially serious under UB sanitizers/optimization. **Confidence:** confirmed latent UB.

`PairList::addAndPair()` stores pairing bits in a 64-bit integer and shifts by `63 - d_descriptors.size()` without a length guard (`Pairlist.h:170-177`). Adding a matching 65th recognized descriptor yields an invalid shift count and C++ undefined behavior. The current Rule4b/Rule5New comparison uses normalized descriptor vectors rather than `d_pairing`, but the unsafe pairing calculation still runs.

The public default constructor also permits an empty list while `getRefDescriptor()` and `compareTo()` index element zero without checking (`Pairlist.h:52,71,119-124`). Internal callers normally populate lists first.

There is also a shared logical error in the otherwise unused pairing accumulator: the first stored descriptor has already been normalized to R/S, but `addAndPair()` compares it with the next *raw* descriptor. Equivalent classes such as initial `R` followed by `M` can therefore fail to set a “like” bit (`Pairlist.h:37-49,170-177`). Current ranking uses normalized-vector `compareTo()`, not `getPairing()`, so this does not currently change CIP output.

**Recommended fix:** remove the unused pairing accumulator, or bound the bit operation/use a dynamic bitset, and assert nonempty preconditions before indexing.

### C++-14 — The node cap can be exceeded by one expansion

**Severity:** low. **Confidence:** confirmed and inherited from Java.

Both versions test the 100,000-node limit only at the beginning of expanding a node (`Digraph.cpp:139-147`; `Digraph.java:177-180`). One expansion can then add several explicit, duplicate, ring, and H nodes. If those additions are terminal, no subsequent expansion necessarily rechecks the count, so the final graph can exceed the advertised cap.

**Recommended fix:** enforce the cap centrally in `addNode()` or preflight the number of nodes an expansion will add.

### C++-15 — Mancude relaxation retains typed atoms with zero typed neighbors

**Severity:** low in current internal use. **Confidence:** confirmed and inherited algorithm issue.

The comment says a resonant atom needs more than one typed neighbor, but the pruning queue initially contains only nodes with exactly one, not zero (`Mancude.cpp:99-132`; Java has the same condition). A zero-degree seeded type survives. For `[N-]1CCC1`, C++ calculates `0/1` for nitrogen instead of retaining atomic number 7; Java creates an invalid `0/0` fraction. Current graph expansion usually does not use that value because its denominator is not greater than one, but the cached result is chemically wrong.

There is a related inherited inconsistency: `RelaxTypes()` counts typed neighbors across all bonds (`Mancude.cpp:105-110,123-129`), while `VisitPart()` connects component members only across ring bonds (`Mancude.cpp:140-149`). An exocyclic typed neighbor can keep a seed alive even though it cannot join the component subsequently used for averaging.

**Recommended fix:** enqueue all typed atoms with count less than two, and calculate/decrement those counts over the same ring-edge relation used by component discovery (unless a different relation is chemically intended and documented). Test isolated charged ring types and exocyclic typed neighbors.

### C++-16 — Failed kekulization may cache a partially modified temporary

**Severity:** low/uncertain robustness risk. **Confidence:** plausible; no end-to-end mislabel reproduced here.

`CIPMol::getBondOrder()` catches `MolSanitizeException` from `MolOps::Kekulize(tmp)` and then unconditionally caches bond types from `tmp` (`CIPMol.cpp:65-78`). If kekulization mutated some bonds before throwing, the cache can contain a mixture of kekulized and aromatic types. Remaining aromatic bonds are treated as order one with a warning (`CIPMol.cpp:93-97`).

At minimum this fallback is approximate for non-kekulizable inputs. A safer design would either guarantee rollback, retry from an untouched copy with a defined fallback, or cache the original bond types consistently when kekulization fails.

### C++-17 — Public selection bitsets are not bounds-validated

**Severity:** low/medium robustness issue. **Confidence:** confirmed unchecked path; Python creates correctly sized masks.

`findConfigs()` trusts every set bit in the exported C++ overload (`CIPLabeler.cpp:49-67`). An out-of-range atom bit reaches RDKit atom lookup; an out-of-range bond bit reaches `CIPMol::getBond()`, which indexes `d_bonds` with unchecked `operator[]` (`CIPMol.cpp:42`). The header does not state an exact-size precondition.

**Recommended fix:** require and validate bitset sizes against the molecule, or validate every discovered index and use checked access. The Python conversion already sizes its masks, so this principally hardens C++ callers.

### C++-18 — Auxiliary-label collisions changed from deterministic overwrite to first-wins

**Severity:** low/latent. **Confidence:** confirmed semantic difference; chemically overlapping trigger not reproduced here.

At one distance Java stores calculated auxiliary labels with `HashMap.put()`, so a later configuration mapped to the same node overwrites an earlier one (`Labeller.java:119-133`). C++ uses `unordered_map.emplace()`, which retains the first value (`CIPLabeler.cpp:146-163`). Its preceding `std::sort` compares distance only and is not stable, so the winner among equal-distance overlapping configurations is not a defined configuration order.

Normally a node corresponds to one relevant stereogenic unit, but overlapping/adjacent axes can violate that assumption. The representation has room for only one node auxiliary descriptor, so the intended collision semantics should be documented and asserted; silently selecting first or last can make such input order-dependent.

### C++-19 — `Digraph` can retain a reference to a temporary `CIPMol`

**Severity:** low/internal lifetime hazard. **Confidence:** confirmed API shape; production configuration lifetimes are safe.

`Digraph(const CIPMol&,...)` stores the adapter by reference (`Digraph.h:54,101`) and does not reject rvalues. Code such as `Digraph g(CIPMol(mol), atom);` leaves `d_mol` dangling at the end of the construction expression, before later lazy expansion uses it. Internal `Configuration` objects retain a longer-lived `CIPMol` owned by the enclosing call, so the normal `assignCIPLabels()` path does not trigger this.

**Recommended fix:** delete a `CIPMol&&` constructor overload or otherwise make the lifetime contract explicit in the type/API.

## 7. C++ performance opportunities

The items below are ordered roughly by likely impact on difficult molecules.

### 7.1 Replace the full visit-vector copy per primary node

Every nonterminal child clones an atom-count-sized `std::vector<char>` (`Node.cpp:112-116`). With `V` molecular atoms and `G` generated digraph nodes, each configuration can spend `O(V * G)` bytes copied/stored solely on path state. Multiple stereocenters each own a separate digraph, magnifying the cost.

Use a persistent parent/path representation, copy-on-write blocks, or a compact bitset plus separate distance storage. DFS backtracking could remove per-node state only as part of a larger rewrite away from the current lazy, persistent, rerootable graph. This change also provides the natural place to fix C++-5.

### 7.2 Replace per-object linked-list allocation and fix edge reservation

`Digraph` uses `std::list<Node>` and `std::list<Edge>` solely to preserve pointer stability (`Digraph.h:116-119`). That causes one allocation per object and poor cache locality. A block arena, object pool, or a carefully chosen stable-address `deque` can preserve pointers while improving locality.

There is also a reversed reservation: C++ reserves four edge slots only for duplicate nodes (`Node.cpp:46-49`), although duplicate terminals normally need one incoming edge and ordinary nodes commonly need about four. Java preallocates four for ordinary nodes (`Node.java:90-92`). Reserve for ordinary nodes or use an inline small vector such as `SmallVector<Edge*,4>`.

### 7.3 Index graph nodes by atom

`seenAtom()` linearly scans all generated nodes (`Digraph.cpp:44-47`). `getNodes(atom)` traverses the current expanded graph for every candidate auxiliary configuration (`Digraph.cpp:80-97`; `CIPLabeler.cpp:109-139`). Maintain at least an atom-index bitset for `seenAtom()`, and preferably atom-to-node occurrence lists with enough reachability/orientation information to accelerate `getNodes()` safely.

### 7.4 Cache expensive comparisons and avoid sorting twice

Insertion sorting is appropriate for the small number of immediate ligands, but a single comparator can recursively traverse a large subtree. `Sort::getGroups()` deep-compares adjacent edges again after sorting (`Sort.cpp:73-90`). Rule stacks repeatedly revisit the same edge pairs, and Rule5New traverses both ligands under both R and S references.

Preserve equality boundaries during prioritization or add a call-scoped memoization cache keyed by rule, edge pair, root/orientation, auxiliary-label generation, and Rule-6/reference state. The state must be part of the key; a global edge-pair cache would be incorrect after rerooting or auxiliary updates.

### 7.5 Reduce recursive-comparison scratch allocation

`SequenceRule::recursiveCompare()` copies two edge vectors at each visited pair and retains two parallel vectors of every queued edge (`SequenceRule.cpp:65-76,105-127`). Small inline edge buffers and one queue of edge pairs would reduce allocation and prevent queue misalignment. A `deque` uses only frontier memory; the current vector/cursor approach has good locality but retains the complete traversal until comparison ends.

Rule4b, Rule5New, `Digraph::changeRoot()`, and `Mancude::RelaxTypes()` use `std::list` as append-while-iterating queues. None requires node-stable list allocation; `deque` or vector-plus-cursor queues would be cheaper.

### 7.6 Make Mancude processing component-oriented

The charged-component pass scans every atom once for every charged resonant part and deduplicates part IDs with a linear search (`Mancude.cpp:207-219,243-267`). Build component-member lists during `VisitParts`, track charged components in an indexed boolean/vector, calculate one final fraction per component, and assign it to members. This changes the pass from approximately `O(parts * atoms)` to `O(atoms + bonds)` while fixing C++-1 cleanly.

Cache the source atom's charge and fractional/resonant state once per `Digraph::expand()` rather than querying inside each bond branch.

### 7.7 Avoid unconditional whole-molecule clone/kekulization on ordinary inputs

The first `getBondOrder()` clones the complete `ROMol` and attempts kekulization even when the molecule contains no aromatic bonds (`CIPMol.cpp:65-78`). Fill the cache directly for a non-aromatic molecule; clone/kekulize only when an aromatic bond is present. Define a consistent failed-kekulization path as part of this work.

Atomic average/exact mass can similarly be cached per original atom/isotope in `CIPMol` instead of querying RDKit's periodic table for every ordinary graph-node occurrence.

### 7.8 Make the preliminary pass budget-aware

The 2,000-call constitutional pass is valuable for easy centers and retains any graph expansion for the second pass, but unresolved centers redo their comparisons and the pass ignores a smaller caller limit. Skip or shrink it when the caller's budget is smaller, and consider an adaptive threshold based on graph growth. Correct budget accounting is prerequisite.

Also distinguish a preliminary pass that timed out from one that completed and returned `UNKNOWN`. A clean constitutional `UNKNOWN` result cannot change when the same constitutional rules are immediately rerun, so that center can proceed directly to auxiliary labeling; only timed-out centers require a constitutional retry.

Because each configuration embeds its own persistent `Digraph`, the all-center preliminary pass can keep many expanded graphs resident until the complete call ends. A configuration already resolved constitutionally no longer needs its own graph: when it later serves as an auxiliary configuration, it is evaluated against the *target center's* graph. Lazy/optional graphs or releasing resolved target graphs can therefore reduce peak memory on molecules with many centers.

### 7.9 Preserve the useful C++ optimizations safely

C++ already has meaningful advantages that should be retained:

- one `CIPMol` shares kekulized bond orders and fractional atomic numbers across all configurations in a call;
- `seenAtom()` avoids auxiliary traversals for completely unreached configurations;
- cursor-indexed vectors avoid Java `LinkedList` allocation during recursive comparison; and
- reusable stateless rule objects can avoid per-call rule construction.

The last optimization only needs immutable eager sorter construction to remove C++-2.

## 8. Recommended fix order

1. Fix charged-component fraction calculation and introduce explicit resonance metadata (C++-1 and C++-6 together). Add atom-renumbering invariance tests before changing other graph logic.
2. Make global rule stacks truly immutable/eagerly initialized (C++-2).
3. Separate the full auxiliary-configuration universe from the selective output mask, and correct `_CIPComputed` semantics (C++-3).
4. Port the Rule 6 magnitude correction and atrop XOR correction (C++-4 and C++-8).
5. Replace `char` visit distances and add long-chain graph tests (C++-5).
6. Make state cleanup atomic across `_CIPCode`, `_CIPNeighborOrder`, and `_CIPComputed` (C++-7).
7. Add unknown-isotope fallback and harden invalid configurations (C++-9 and C++-11).
8. Correct Rule4b indexing, iteration accounting, PairList bounds, and the node-cap check.
9. Validate public bitsets, harden internal adapter lifetimes, and define/assert auxiliary-collision behavior.
10. Profile after correctness fixes, then prioritize visit-state storage, graph allocation, atom indexing, and comparison caching.

## 9. Regression tests to add

- **Charged Mancude:** assert every atom in `[CH-]1C=CC=C1` has `24/5`; repeat after atom renumbering and compare final CIP labels.
- **Reduced fraction:** use a charged resonant component whose fraction reduces to an integer and assert that resonance duplicate semantics are retained.
- **Deep graph:** expand a 260-atom acyclic chain with H generation disabled; require exactly 260 nodes, one primary occurrence per atom, and no ring duplicate.
- **Selective auxiliary labels:** use the para example above; selecting only atom 7 must still compute `r` while committing no primary labels to atoms 3 and 9. A partial call must not set the full `_CIPComputed` marker.
- **Rule 6:** port the structures added with Java commit `def7380a` and assert descriptor case as well as priority.
- **Atrop both-end pseudo:** construct an axis for which both ends use pseudoasymmetric decisions and require upper-case `M/P`.
- **State invalidation:** relabel after (a) making two ligands identical, (b) removing an atom chiral tag, and (c) removing bond stereo; assert both `_CIPCode` and `_CIPNeighborOrder` are absent where appropriate.
- **Isotopes:** compare known and deliberately unlisted isotope mass numbers.
- **Concurrency:** make simultaneous first-ever calls on separate molecules under ThreadSanitizer.
- **Limits/robustness:** exact node-cap enforcement, more than 64 pair descriptors under UBSan, mismatched Rule4b groups, malformed tetrahedral tags, and rejected atrop markers.
- **Selection robustness:** oversized bitsets and deliberately overlapping auxiliary configurations.
- **Mancude pruning:** exocyclic typed neighbors must not preserve atoms excluded from the eventual ring component.

## 10. Verification performed for this report

The current optimized C++ `testCIPLabeler` executable was run with the required RDKit runtime/data environment and passed:

```text
All tests passed (6623 assertions in 22 test cases)
```

Targeted local C++ probes separately confirmed the charged-Mancude prefix fractions, the 260-node path-distance overflow behavior, stale `_CIPNeighborOrder`, and the selective pseudoasymmetric-label failure described above.

The Java Maven tests were not run because `mvn` is not installed in this workspace. Java findings are based on complete source comparison and the checked-in Java test/validation data and commit history.

## 11. Post-series optimization for the V157894 pathological case

This follow-up was performed on 2026-08-06 after the correctness series and
performance patches `0001` through `0013`. It studies the local test named
`Long long running calculation`, whose V3000 input has 92 atoms, 106 bonds,
and nine tetrahedral configurations. The molecule contains a 64-atom,
76-bond biconnected cyclic block. Consequently, asking the unfolded digraph
for every route from one configuration to every other configuration can
materialize thousands of distinct simple-path occurrences even though the
chemically relevant symmetry is local.

### 11.1 Remaining source of broad auxiliary work

Patch `0013` retains the support of an exact constitutional-automorphism
witness. If that witness moves any other stereo annotation, auxiliary
discovery conservatively targets the complete connected component. The exact
matcher is free to return a product of independent automorphisms, however. In
V157894 the local symmetry at atom 81 swaps the two directions of its
six-membered ring, but an arbitrary valid witness can additionally swap a
remote symmetric fragment around atom 20. The remote movement is irrelevant
to the requested ligand equality, yet it activates the component-wide
fallback.

Static topology analysis also identified a genuine local dependency chain
`68 <-> 71 <-> 91`. Those automorphisms necessarily move a neighboring
configuration's carriers, so the conservative fallback cannot always be
removed. Approximate nonduplicate occurrence counts explain the distinction:
the false-positive component-wide fallback at center 81 can request 10,396
configuration occurrences, while the genuine centers 68, 71, and 91 request
about 15,310 occurrences together. A minimal dependency closure for the
latter three would contain about 5,776 occurrences, but safely deriving that
closure for every rerooted occurrence requires provenance that the current
global `seenAtom()` state does not retain.

### 11.2 Implemented improvements

Performance patch `0014` adds four exact, molecule-independent
optimizations:

1. `CIPMol` now records each configuration's foci and complete local
   annotation (foci, carriers, and focus neighbors). When a constitutional
   witness appears to move another configuration, the matcher retries while
   fixing every unrelated annotation pointwise. A successful retry replaces
   only the witness support; a failed or bounded search retains the original
   witness and the existing broad fallback. The fixed atoms are materialized
   as a dense word mask so candidate checks remain constant time inside VF2.
   Ambiguous shared foci or unavailable detailed metadata disable this
   optimization conservatively, malformed metadata is rejected, and the
   earlier `setConfigurationFoci()` interface remains available.
2. A configuration-root ligand tie proved by an annotation-preserving exact
   automorphism is recorded separately. Bond and axial configurations cannot
   break such a tie. Tetrahedral configurations skip auxiliary discovery only
   when their constitutional partition has more than two groups (or fewer
   than four root edges), preserving the existing Rule-6 reference retry for
   one- and two-group cases. This avoids enumerating return paths to already
   reached configurations when no auxiliary descriptor can make the primary
   priorities unique.
3. The exact constitutional sibling shortcut is extended from the current
   digraph root to original-child ring directions below it. The molecular
   automorphism must fix the occurrence's immutable original-root path
   pointwise, preserving visited-state and Rule-1b ring-duplicate distances.
   Below-root use is restricted to ring bonds; ordinary pendant trees retain
   the cheaper bridge proof and existing traversal.
4. Multi-target auxiliary discovery builds a lazy multi-source molecular path
   forest. An unblocked forest path is a positive reachability certificate and
   avoids a fresh molecular BFS. If the occurrence's visited path blocks that
   certificate, the previous residual-connectivity BFS still runs, so the
   optimization cannot prune a reachable target.

The focused Rule-1a regressions exercise an unexpanded exact ring symmetry
below the digraph root and verify that a configuration-preserving root
symmetry records an auxiliary-invariant tie. The large V157894 fixture and the
local node-cap adjustment remain outside the formal patch series, consistent
with the policy used for patch `0013`.

### 11.3 Correctness boundaries and remaining work

No atom numbers, molecule names, component-size thresholds, or aromaticity
special cases are used by these gates. All positive pruning decisions are
backed by an exact automorphism, a pointwise-fixed occurrence path, and (where
required) pointwise-fixed stereo annotations. Failure, unsupported topology,
overlapping configuration ownership, or search-budget exhaustion always
falls back to the pre-existing algorithm.

The genuine `68 <-> 71 <-> 91` dependency still offers a larger future
optimization. It should use occurrence-scoped certificates recording which
specific sibling comparison hid which configurations, rather than rescanning
global seen atoms after target collection. Target collection itself can
materialize unrelated configuration foci along detour paths, so treating all
newly seen atoms as dependencies recreates the component-wide blow-up. A
safe implementation would need a monotone worklist (and rollback/restart of
uncommitted auxiliary spheres) or an equivalent provenance-preserving design.

Per the task constraint, no build, unit test, benchmark, or labeling command
was run for this follow-up. Review was limited to source inspection, static
graph/path analysis, whitespace checks, and clean patch-series replay.
