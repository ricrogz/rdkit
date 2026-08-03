# Consolidated improvement findings

This report consolidates the Markdown audits under `improvements/`. Findings are
ordered by likely impact: memory safety, undefined behavior, hangs, and crashes
come first; user-visible correctness and compatibility defects follow; bounded
edge cases and efficiency issues come next; maintenance-only observations are
last. Line numbers are approximate and refer to the code revision against which
the patches were produced.

## Critical impact

| Finding | Code reference (approx.) | Related patch |
|---|---|---|
| A default heteroatom iterator can delete an indeterminate query pointer. | `rdkit/Code/GraphMol/AtomIterators.h:~100` | [`001-initialize-heteroatom-query.patch`](patches/GraphMol/AtomIterators/001-initialize-heteroatom-query.patch) |
| Optional neighboring bonds are dereferenced even when null, and copied bonds break identity comparisons during canonicalization. | `rdkit/Code/GraphMol/Canon.cpp:~184-275, ~546` | [`001-preserve-optional-bond-pointers.patch`](patches/GraphMol/Canon/001-preserve-optional-bond-pointers.patch) |
| A caller-supplied partial core match can make `replaceCore()` dereference `end()`. | `rdkit/Code/GraphMol/ChemTransforms/ChemTransforms.cpp:~495-512` | [`005-validate-multiply-mapped-neighbor.patch`](patches/GraphMol/ChemTransforms/005-validate-multiply-mapped-neighbor.patch) |
| Coordinate generation dereferences a null result after a recoverable molzip failure. | `rdkit/Code/GraphMol/ChemTransforms/MolZip.cpp:~710` | [`002-guard-failed-coordinate-generation.patch`](patches/GraphMol/ChemTransforms/002-guard-failed-coordinate-generation.patch) |
| Coordinate remapping dereferences a missing atom and can leak its raw conformer on exceptions. | `rdkit/Code/GraphMol/ChemTransforms/MolZip.cpp:~715-744` | [`003-validate-coordinate-atom-lookup.patch`](patches/GraphMol/ChemTransforms/003-validate-coordinate-atom-lookup.patch) |
| An atom index equal to the atom count writes beyond `revOrder`; duplicates also corrupt the reverse map. | `rdkit/Code/GraphMol/Renumber.cpp:~21-35` | [`001-validate-renumber-permutation.patch`](patches/GraphMol/Renumber/001-validate-renumber-permutation.patch) |
| Colon layout uses string bounds instead of rectangle bounds and can index `rects_` out of range. | `rdkit/Code/GraphMol/MolDraw2D/AtomSymbol.cpp:~135-145` | [`002-use-rectangle-bounds-for-colons.patch`](patches/GraphMol/MolDraw2D/002-use-rectangle-bounds-for-colons.patch) |
| An unterminated markup tag causes the atom-label adjustment loop to run forever. | `rdkit/Code/GraphMol/MolDraw2D/AtomSymbol.cpp:~129-137` | [`003-handle-unterminated-label-markup.patch`](patches/GraphMol/MolDraw2D/003-handle-unterminated-label-markup.patch) |
| Timeout handling clears the saved clique vector and then calls `front()` on it. | `rdkit/Code/GraphMol/RascalMCES/RascalMCES.cpp:~822-830` | [`001-avoid-front-after-clear.patch`](patches/GraphMol/RascalMCES/001-avoid-front-after-clear.patch) |
| A permutation search dereferences an iterator before checking it against `end()`. | `rdkit/Code/RDGeneral/utils.h:~69` | [`001-check-end-before-dereference.patch`](patches/RDGeneral/001-check-end-before-dereference.patch) |
| The default metric-matrix calculator can call an indeterminate function pointer. | `rdkit/Code/DataManip/MetricMatrixCalc/MetricMatrixCalc.h:~79-105` | [`001-validate-metric-function.patch`](patches/DataManip/MetricMatrixCalc/001-validate-metric-function.patch) |
| Short molecular diameters leave GETAWAY descriptor bins uninitialized before aggregation. | `rdkit/Code/GraphMol/Descriptors/GETAWAY.cpp:~347-357` | [`001-zero-getaway-distance-bins.patch`](patches/GraphMol/Descriptors/001-zero-getaway-distance-bins.patch) |
| Reset MCS objects contain indeterminate query pointers and threshold state that later code reads. | `rdkit/Code/GraphMol/FMCS/MaximumCommonSubgraph.h:~41-71` | [`001-initialize-mcs-state.patch`](patches/GraphMol/FMCS/001-initialize-mcs-state.patch) |
| SMARTS error handling throws a pointer into a temporary string that is already destroyed when caught. | `rdkit/Code/GraphMol/SmilesParse/SmartsWrite.cpp:~463` | [`001-throw-owning-smarts-error.patch`](patches/GraphMol/SmilesParse/001-throw-owning-smarts-error.patch) |

## High impact

| Finding | Code reference (approx.) | Related patch |
|---|---|---|
| A fixed-width subset counter shifts past the width of `unsigned int` for large fused aromatic systems. | `rdkit/Code/GraphMol/Kekulize.cpp:~456-492` | [`001-width-independent-question-enumerator.patch`](patches/GraphMol/Kekulize/001-width-independent-question-enumerator.patch) |
| Stereogroup metadata is read from a moved-from atom vector and consequently is not written. | `rdkit/Code/GraphMol/Canon.cpp:~1726-1737` | [`002-set-stereogroup-before-move.patch`](patches/GraphMol/Canon/002-set-stereogroup-before-move.patch) |
| ROMol move operations are declared `noexcept` even though owner repair and container/property moves can throw, forcing termination. | `rdkit/Code/GraphMol/ROMol.h:~515, ~542` | [`001-remove-unsafe-move-noexcept.patch`](patches/GraphMol/ROMol/001-remove-unsafe-move-noexcept.patch) |
| Fragment subset masks use a signed, platform-dependent literal, causing invalid shifts on Windows and at high bit positions. | `rdkit/Code/GraphMol/ChemTransforms/MolFragmenter.cpp:~325-352` | [`006-use-uint64-fragment-masks.patch`](patches/GraphMol/ChemTransforms/006-use-uint64-fragment-masks.patch) |
| Duplicate linker bond validation assigns instead of compares, allowing inconsistent bond types. | `rdkit/Code/GraphMol/ChemTransforms/MolZip.cpp:~514` | [`001-compare-linker-bond-types.patch`](patches/GraphMol/ChemTransforms/001-compare-linker-bond-types.patch) |
| Embedded-atom assignment omits the atom ID, combining source geometry with the destination identity. | `rdkit/Code/GraphMol/Depictor/EmbeddedFrag.h:~57-61` | [`002-copy-embedded-atom-id.patch`](patches/GraphMol/Depictor/002-copy-embedded-atom-id.patch) |
| A composition with exactly 64 external bonds passes its guard and performs a 64-bit-width shift. | `rdkit/Code/GraphMol/FMCS/Seed.cpp:~271` | [`002-reject-width-sized-composition.patch`](patches/GraphMol/FMCS/002-reject-width-sized-composition.patch) |
| Variadic V2000 formatting passes an unsigned parity flag to `%d`, which is undefined behavior. | `rdkit/Code/GraphMol/FileParsers/MolFileWriter.cpp:~634-650` | [`001-fix-parity-format-type.patch`](patches/GraphMol/FileParsers/001-fix-parity-format-type.patch) |
| The MAE writer validates and dereferences its stream argument after moving from it. | `rdkit/Code/GraphMol/FileParsers/MaeWriter.cpp:~542-549` | [`003-check-stored-mae-stream.patch`](patches/GraphMol/FileParsers/003-check-stored-mae-stream.patch) |
| The TDT writer destructor can let an ostream exception escape and invoke `std::terminate`. | `rdkit/Code/GraphMol/FileParsers/TDTWriter.cpp:~66-76` | [`004-make-tdt-destructor-nothrow.patch`](patches/GraphMol/FileParsers/004-make-tdt-destructor-nothrow.patch) |
| A default MAE supplier exposes indeterminate position and length counters. | `rdkit/Code/GraphMol/FileParsers/MolSupplier.h:~685-690` | [`005-initialize-mae-supplier-counters.patch`](patches/GraphMol/FileParsers/005-initialize-mae-supplier-counters.patch) |
| The MMFF aromatic marker shifts a signed integer into its sign bit. | `rdkit/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp:~135` | [`001-use-unsigned-aromatic-bit.patch`](patches/GraphMol/ForceFieldHelpers/001-use-unsigned-aromatic-bit.patch) |
| Molecule-sized `alloca` buffers, including one per DFS, can exhaust the stack on large input. | `rdkit/Code/GraphMol/MolHash/hashfunctions.cpp:~1077-1140` | [`002-replace-molecule-sized-alloca.patch`](patches/GraphMol/MolHash/002-replace-molecule-sized-alloca.patch) |
| Multithreaded processing drops trailing parse failures and can return too few result entries. | `rdkit/Code/GraphMol/MolProcessing/MolProcessing.cpp:~51-66` | [`001-preserve-trailing-parse-failures.patch`](patches/GraphMol/MolProcessing/001-preserve-trailing-parse-failures.patch) |
| Rascal result copies omit options, scores, and null pointer state, so copies behave differently from their source. | `rdkit/Code/GraphMol/RascalMCES/RascalResult.cpp:~81-111` | [`002-copy-complete-result-state.patch`](patches/GraphMol/RascalMCES/002-copy-complete-result-state.patch) |
| Python catalog loading references undefined `inD`, breaking cached-list validation. | `rdkit/rdkit/Chem/BuildFragmentCatalog.py:~580` | [`001-validate-cached-list-size.patch`](patches/rdkit/Chem/001-validate-cached-list-size.patch) |
| The feature-arrow path refers to unimported `RDGeometry` and fails with `NameError`. | `rdkit/rdkit/Chem/Features/ShowFeats.py:~23-26` | [`001-use-imported-geometry-module.patch`](patches/rdkit/Chem/Features/001-use-imported-geometry-module.patch) |
| The non-PostgreSQL fingerprint database fallback passes undefined `data`. | `rdkit/rdkit/Chem/Fingerprints/MolSimilarity.py:~121-127` | [`001-build-forward-db-result.patch`](patches/rdkit/Chem/Fingerprints/001-build-forward-db-result.patch) |
| Both molecule database loaders compute redraw coordinates for undefined `m` rather than the current molecule. | `rdkit/rdkit/Chem/MolDb/Loader_orig.py:~70`; `Loader_sa.py:~65` | [`001-redraw-current-molecule.patch`](patches/rdkit/Chem/MolDb/001-redraw-current-molecule.patch) |
| The lazy pharmacophore generator discards valid feature matches, calls an API incorrectly, and stores undefined `pattMatches`. | `rdkit/rdkit/Chem/Pharm2D/LazyGenerator.py:~77-87` | [`001-build-lazy-feature-match-list.patch`](patches/rdkit/Chem/Pharm2D/001-build-lazy-feature-match-list.patch) |
| Packed lower-triangular matrix indices use true division and become invalid float indices under Python 3. | `rdkit/rdkit/ML/Cluster/ClusterUtils.py:~89, ~103`; `Resemblance.py:~128` | [`001-use-integer-packed-matrix-indices.patch`](patches/rdkit/ML/Cluster/001-use-integer-packed-matrix-indices.patch) |
| Non-shuffle activity randomization iterates an undefined collection and generates off-by-one class labels. | `rdkit/rdkit/ML/Data/DataUtils.py:~652-663` | [`001-fix-randomized-result-generation.patch`](patches/rdkit/ML/Data/001-fix-randomized-result-generation.patch) |
| The sping formatter cannot import on supported Python versions because it depends on removed `xmllib` and uses an invalid package import. | `rdkit/rdkit/sping/stringformat.py:~36-46, ~164-287, ~393` | [`001-port-string-formatter-to-python3.patch`](patches/rdkit/sping/001-port-string-formatter-to-python3.patch) |

## Medium impact

| Finding | Code reference (approx.) | Related patch |
|---|---|---|
| Negative ring-stereo indices are offset using unsigned arithmetic and converted back to `int` implementation-dependently. | `rdkit/Code/GraphMol/RWMol.cpp:~154` | [`001-offset-negative-ring-index-as-int.patch`](patches/GraphMol/RWMol/001-offset-negative-ring-index-as-int.patch) |
| The depiction neighbor-angle loop pairs each neighbor with itself, adding invalid zero-angle work. | `rdkit/Code/GraphMol/Depictor/EmbeddedFrag.cpp:~111-118` | [`001-skip-self-neighbor-angle.patch`](patches/GraphMol/Depictor/001-skip-self-neighbor-angle.patch) |
| Negating an unsigned neighbor count and converting it to `int` is implementation-defined and can alter stereochemical sorting. | `rdkit/Code/GraphMol/FileParsers/MolFileStereochem.cpp:~136-143` | [`002-negate-neighbor-count-safely.patch`](patches/GraphMol/FileParsers/002-negate-neighbor-count-safely.patch) |
| Shape deserialization begins with indeterminate scalar indices if archive extraction is partial or fails. | `rdkit/Code/GraphMol/GaussianShape/ShapeInput.h:~302-320` | [`001-initialize-shape-indices.patch`](patches/GraphMol/GaussianShape/001-initialize-shape-indices.patch) |
| Default Marvin objects leave serialized isotope, coordinates, and center pointer fields indeterminate. | `rdkit/Code/GraphMol/MarvinParse/MarvinDefs.cpp:~763, ~1444, ~1746` | [`001-initialize-serialized-members.patch`](patches/GraphMol/MarvinParse/001-initialize-serialized-members.patch) |
| DrawMol constructors do not establish safe defaults for height and annotation state. | `rdkit/Code/GraphMol/MolDraw2D/DrawMol.h:~268, ~309` | [`001-initialize-drawmol-members.patch`](patches/GraphMol/MolDraw2D/001-initialize-drawmol-members.patch) |
| Tautomer hash suffix formatting mismatches unsigned values with `%d` and subtracts charge in unsigned arithmetic. | `rdkit/Code/GraphMol/MolHash/hashfunctions.cpp:~989-1001` | [`001-build-hash-suffix-type-safely.patch`](patches/GraphMol/MolHash/001-build-hash-suffix-type-safely.patch) |
| Chromosome diagnostics and selection can read fitness before the first score. | `rdkit/Code/GraphMol/RGroupDecomposition/RGroupGa.h:~79` | [`001-initialize-chromosome-fitness.patch`](patches/GraphMol/RGroupDecomposition/001-initialize-chromosome-fitness.patch) |
| Substitution scoring can divide by zero for an empty molecule or overflow before conversion for a large molecule. | `rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp:~40-49` | [`001-safe-substitution-score-denominator.patch`](patches/GraphMol/Substruct/001-safe-substitution-score-denominator.patch) |
| Python sequence conversion leaks an allocated atom map if row extraction or indexing throws. | `rdkit/Code/RDBoost/Wrap.cpp:~94-112` | [`001-own-atom-map-during-conversion.patch`](patches/RDBoost/001-own-atom-map-during-conversion.patch) |
| Partial clustering scratch allocation leaks earlier buffers when a later allocation fails. | `rdkit/Code/SimDivPickers/HierarchicalClusterPicker.cpp:~27-45, ~73-79` | [`001-use-raii-clustering-buffers.patch`](patches/SimDivPickers/001-use-raii-clustering-buffers.patch) |
| Python 2's removed `file()` builtin breaks command-line input and output paths. | `rdkit/rdkit/Chem/ChemUtils/TemplateExpand.py:~413, ~446` | [`001-replace-file-builtin.patch`](patches/rdkit/Chem/ChemUtils/001-replace-file-builtin.patch) |
| Selection lookup calls removed `dict.iteritems()` under Python 3. | `rdkit/rdkit/Chem/DSViewer.py:~179` | [`001-use-dict-items.patch`](patches/rdkit/Chem/DSViewer/001-use-dict-items.patch) |
| Arrowhead vertices all use whole-turn angles and collapse onto one point. | `rdkit/rdkit/Chem/Features/ShowFeats.py:~63-69` | [`002-distribute-arrowhead-vertices.patch`](patches/rdkit/Chem/Features/002-distribute-arrowhead-vertices.patch) |
| Fingerprint parse-failure logging references undefined `smi`, masking the original error. | `rdkit/rdkit/Chem/Fingerprints/FingerprintMols.py:~91` | [`002-report-failed-record-id.patch`](patches/rdkit/Chem/Fingerprints/002-report-failed-record-id.patch) |
| The SD supplier loop uses the removed iterator `.next()` method. | `rdkit/rdkit/Chem/Fingerprints/FingerprintMols.py:~190` | [`003-use-next-builtin.patch`](patches/rdkit/Chem/Fingerprints/003-use-next-builtin.patch) |
| Verbose pharmacophore matching reads local `bm` before it is assigned. | `rdkit/rdkit/Chem/Pharm3D/EmbedLib.py:~1180` | [`001-fix-pre-update-bounds-diagnostic.patch`](patches/rdkit/Chem/Pharm3D/001-fix-pre-update-bounds-diagnostic.patch) |
| Unordered subshape traversal calculates but ignores its alternate triangle start index. | `rdkit/rdkit/Chem/Subshape/SubshapeAligner.py:~35` | [`001-use-selected-triangle-start.patch`](patches/rdkit/Chem/Subshape/001-use-selected-triangle-start.patch) |
| The random-access database molecule supplier calls removed iterator `.next()`. | `rdkit/rdkit/Chem/Suppliers/DbMolSupplier.py:~127` | [`001-use-next-builtin.patch`](patches/rdkit/Chem/Suppliers/001-use-next-builtin.patch) |
| Five progress-report paths use Python 2 `print >>` semantics and fail when invoked on Python 3. | `rdkit/rdkit/Chem/fmcs/fmcs.py:~1629-1637, ~2085-2092, ~2722-2816` | [`001-modernize-stderr-reporting.patch`](patches/rdkit/Chem/fmcs/001-modernize-stderr-reporting.patch) |
| An unmatched opening parenthesis is accepted as a complete descriptor call and rewritten misleadingly. | `rdkit/rdkit/ML/Descriptors/Parser.py:~279-284` | [`001-require-balanced-parentheses.patch`](patches/rdkit/ML/Descriptors/001-require-balanced-parentheses.patch) |
| Invalid database query results reference undefined `res`, masking the intended validation error. | `rdkit/rdkit/VLib/NodeLib/DbPickleSupplier.py:~112, ~137` | [`001-report-current-query-result.patch`](patches/rdkit/VLib/NodeLib/001-report-current-query-result.patch) |
| Missing AFM glyph widths are compared rather than assigned, so they remain zero. | `rdkit/rdkit/sping/PDF/pdfmetrics.py:~288` | [`001-fill-missing-glyph-widths.patch`](patches/rdkit/sping/PDF/001-fill-missing-glyph-widths.patch) |

## Low impact

| Finding | Code reference (approx.) | Related patch |
|---|---|---|
| Stereogroup cleanup copies each group and its atom/bond vectors unnecessarily. | `rdkit/Code/GraphMol/Atropisomers.cpp:~485`; `Chirality.cpp:~2163` | [`003-avoid-stereogroup-cleanup-copies.patch`](patches/GraphMol/003-avoid-stereogroup-cleanup-copies.patch) |
| `QueryOps.h` and `new_canon.h` contain duplicate standard-library includes. | `rdkit/Code/GraphMol/QueryOps.h:~24`; `new_canon.h:~23` | [`004-remove-duplicate-standard-includes.patch`](patches/GraphMol/004-remove-duplicate-standard-includes.patch) |
| A Kekulize ring-ordering comment names the wrong wedged-bond endpoint. | `rdkit/Code/GraphMol/Kekulize.cpp:~666` | [`002-correct-wedge-endpoint-comment.patch`](patches/GraphMol/Kekulize/002-correct-wedge-endpoint-comment.patch) |
| Reaction-prefix truncation allocates and copies an avoidable substring. | `rdkit/Code/GraphMol/ChemReactions/DaylightParser.cpp:~182` | [`001-resize-reaction-prefix.patch`](patches/GraphMol/ChemReactions/001-resize-reaction-prefix.patch) |
| A molzip bond-selection comment describes obsolete single-bond-only behavior. | `rdkit/Code/GraphMol/ChemTransforms/MolZip.cpp:~128-135` | [`004-correct-zip-bond-type-comment.patch`](patches/GraphMol/ChemTransforms/004-correct-zip-bond-type-comment.patch) |
| `MolFragmenter.cpp` redundantly includes the specific trim header after the string umbrella. | `rdkit/Code/GraphMol/ChemTransforms/MolFragmenter.cpp:~22` | [`007-remove-redundant-trim-include.patch`](patches/GraphMol/ChemTransforms/007-remove-redundant-trim-include.patch) |
| JSON parsing scans an entire C string just to test whether it is empty. | `rdkit/Code/GraphMol/ChemTransforms/MolFragmenterJSONParser.cpp:~19-30` | [`008-simplify-empty-json-check.patch`](patches/GraphMol/ChemTransforms/008-simplify-empty-json-check.patch) |
| `RDDepictor.cpp` repeats two project includes. | `rdkit/Code/GraphMol/Depictor/RDDepictor.cpp:~28` | [`003-remove-duplicate-includes.patch`](patches/GraphMol/Depictor/003-remove-duplicate-includes.patch) |
| Fingerprint implementation files repeat project and Boost includes. | `rdkit/Code/GraphMol/Fingerprints/FingerprintUtil.cpp:~14-25`; `RDKitFPGenerator.cpp:~26` | [`001-remove-duplicate-includes.patch`](patches/GraphMol/Fingerprints/001-remove-duplicate-includes.patch) |
| Marvin parser and writer files repeat their own project includes. | `rdkit/Code/GraphMol/MarvinParse/MarvinParser.cpp:~32`; `MarvinWriter.cpp:~37` | [`002-remove-duplicate-includes.patch`](patches/GraphMol/MarvinParse/002-remove-duplicate-includes.patch) |
| `FeatureParser.cpp` includes `boost/tokenizer.hpp` twice. | `rdkit/Code/GraphMol/MolChemicalFeatures/FeatureParser.cpp:~17` | [`001-remove-duplicate-tokenizer-include.patch`](patches/GraphMol/MolChemicalFeatures/001-remove-duplicate-tokenizer-include.patch) |
| Four common drawing shape constructors make an avoidable second point-vector copy. | `rdkit/Code/GraphMol/MolDraw2D/DrawShape.h:~142-240`; `DrawShape.cpp:~314-717` | [`004-pass-shape-points-by-reference.patch`](patches/GraphMol/MolDraw2D/004-pass-shape-points-by-reference.patch) |
| MolDraw2D implementation files repeat project includes. | `rdkit/Code/GraphMol/MolDraw2D/DrawTextSVG.cpp:~21`; `MolDraw2D.cpp:~17` | [`005-remove-duplicate-includes.patch`](patches/GraphMol/MolDraw2D/005-remove-duplicate-includes.patch) |
| Two MolDraw2D utility scripts contain unused Python imports. | `rdkit/Code/GraphMol/MolDraw2D/side_by_side_images.py:~1`; `update_hash_codes.py:~8` | [`006-remove-unused-python-imports.patch`](patches/GraphMol/MolDraw2D/006-remove-unused-python-imports.patch) |
| Fingerprint API documentation describes an obsolete return type and contains a typo. | `rdkit/Code/GraphMol/MolProcessing/MolProcessing.cpp:~101-114` | [`002-correct-fingerprint-api-documentation.patch`](patches/GraphMol/MolProcessing/002-correct-fingerprint-api-documentation.patch) |
| The MolProcessing wrapper repeats an include already exposed by `MolProcessing.h`. | `rdkit/Code/GraphMol/MolProcessing/Wrap/rdMolProcessing.cpp:~11` | [`001-remove-redundant-general-file-reader-include.patch`](patches/GraphMol/MolProcessing/Wrap/001-remove-redundant-general-file-reader-include.patch) |
| `Metal.cpp` includes `SubstructMatch.h` twice. | `rdkit/Code/GraphMol/MolStandardize/Metal.cpp:~15` | [`001-remove-duplicate-substruct-include.patch`](patches/GraphMol/MolStandardize/001-remove-duplicate-substruct-include.patch) |
| Recursive SMARTS combination copies three read-only strings at every query-tree node. | `rdkit/Code/GraphMol/SmilesParse/SmartsWrite.cpp:~42-53` | [`002-pass-smart-parts-by-reference.patch`](patches/GraphMol/SmilesParse/002-pass-smart-parts-by-reference.patch) |
| Two SMILES parser files each include `Chirality.h` twice. | `rdkit/Code/GraphMol/SmilesParse/CXSmilesOps.cpp:~21`; `SmilesParse.cpp:~34` | [`003-remove-duplicate-chirality-includes.patch`](patches/GraphMol/SmilesParse/003-remove-duplicate-chirality-includes.patch) |
| `ShapeEncoder.cpp` repeats two geometry includes. | `rdkit/Code/GraphMol/ShapeHelpers/ShapeEncoder.cpp:~11-19` | [`001-remove-duplicate-geometry-includes.patch`](patches/GraphMol/ShapeHelpers/001-remove-duplicate-geometry-includes.patch) |
| `Synthon.cpp` includes `MolDescriptors.h` twice. | `rdkit/Code/GraphMol/SynthonSpaceSearch/Synthon.cpp:~13` | [`001-remove-duplicate-descriptor-include.patch`](patches/GraphMol/SynthonSpaceSearch/001-remove-duplicate-descriptor-include.patch) |
| `MolOps.cpp` includes `Chirality.h` twice. | `rdkit/Code/GraphMol/Wrap/MolOps.cpp:~35` | [`001-remove-duplicate-chirality-include.patch`](patches/GraphMol/Wrap/001-remove-duplicate-chirality-include.patch) |
| `RDValue-taggedunion.h` includes `<cstdint>` twice. | `rdkit/Code/RDGeneral/RDValue-taggedunion.h:~39` | [`002-remove-duplicate-cstdint-include.patch`](patches/RDGeneral/002-remove-duplicate-cstdint-include.patch) |
| Similarity-map generation computes an unused quadratic distance matrix. | `rdkit/rdkit/Chem/Draw/SimilarityMaps.py:~254-258` | [`001-remove-unused-distance-matrix.patch`](patches/rdkit/Chem/Draw/001-remove-unused-distance-matrix.patch) |
| SVG grid creation allocates an unused list of cell strings. | `rdkit/rdkit/Chem/Draw/__init__.py:~562-566` | [`002-remove-unused-svg-blocks.patch`](patches/rdkit/Chem/Draw/002-remove-unused-svg-blocks.patch) |
| A corrected activity-randomization branch still carries a stale “doesn't work” comment. | `rdkit/rdkit/ML/Data/DataUtils.py:~655` | [`002-remove-stale-randomization-comment.patch`](patches/rdkit/ML/Data/002-remove-stale-randomization-comment.patch) |

## Informational audit observations

| Finding | Code reference | Related patch |
|---|---|---|
| `FeatTree.cpp`, `FeatTreeUtils.cpp`, and their test appear to be an orphaned subsystem absent from all CMake targets; deletion requires a separate ownership decision. | `rdkit/Code/GraphMol/Basement/FeatTrees/` | None provided |
| Twenty-two empty Python files are package markers, not demonstrably dead modules, so they should not be deleted solely for being empty. | `rdkit/rdkit/{Chem,ML,utils,sping,Avalon}/**/__init__.py` (enumerated below) | None provided |

The empty package markers are:

- `rdkit/rdkit/Chem/Pharm3D/__init__.py`
- `rdkit/rdkit/Chem/AtomPairs/__init__.py`
- `rdkit/rdkit/Chem/SimpleEnum/__init__.py`
- `rdkit/rdkit/Chem/MolKey/__init__.py`
- `rdkit/rdkit/Chem/Features/__init__.py`
- `rdkit/rdkit/Chem/ChemUtils/__init__.py`
- `rdkit/rdkit/Chem/FeatMaps/__init__.py`
- `rdkit/rdkit/Chem/Suppliers/__init__.py`
- `rdkit/rdkit/Chem/MolDb/__init__.py`
- `rdkit/rdkit/Chem/Fraggle/__init__.py`
- `rdkit/rdkit/Chem/Scaffolds/__init__.py`
- `rdkit/rdkit/Chem/Subshape/__init__.py`
- `rdkit/rdkit/utils/__init__.py`
- `rdkit/rdkit/ML/Cluster/__init__.py`
- `rdkit/rdkit/ML/Data/__init__.py`
- `rdkit/rdkit/ML/MLUtils/__init__.py`
- `rdkit/rdkit/ML/SLT/__init__.py`
- `rdkit/rdkit/ML/Descriptors/__init__.py`
- `rdkit/rdkit/ML/Scoring/__init__.py`
- `rdkit/rdkit/sping/ReportLab/__init__.py`
- `rdkit/rdkit/sping/Qt/__init__.py`
- `rdkit/rdkit/Avalon/__init__.py`
