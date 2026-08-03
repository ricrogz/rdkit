# Basement/FeatTrees

## Unused files

`FeatTree.cpp` and `FeatTreeUtils.cpp` are not referenced by any CMake target,
and the feature-tree directory has no `CMakeLists.txt`. Searches outside this
directory find no production callers; only its orphaned `testFeatTrees.cpp`
uses the API. These files therefore do not participate in the build and appear
to be an abandoned subsystem. Per the audit requirements, no deletion patch is
provided.
