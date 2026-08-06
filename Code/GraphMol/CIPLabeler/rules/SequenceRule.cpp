//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDGeneral/ControlCHandler.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <utility>

#include "SequenceRule.h"

#include "../CIPMol.h"
#include "../Digraph.h"

namespace RDKit {
namespace CIPLabeler {

namespace {

struct ComparisonKey {
  std::uint64_t ruleId;
  const Edge *first;
  const Edge *second;
  const Node *firstEnd;
  const Node *secondEnd;
  bool firstAtRoot;
  bool secondAtRoot;
  const Atom *firstRule6Ref;
  const Atom *secondRule6Ref;

  bool operator==(const ComparisonKey &other) const {
    return ruleId == other.ruleId && first == other.first &&
           second == other.second && firstEnd == other.firstEnd &&
           secondEnd == other.secondEnd && firstAtRoot == other.firstAtRoot &&
           secondAtRoot == other.secondAtRoot &&
           firstRule6Ref == other.firstRule6Ref &&
           secondRule6Ref == other.secondRule6Ref;
  }
};

template <typename T>
void hashCombine(std::size_t &seed, const T &value) {
  const auto hash = std::hash<T>{}(value);
  seed ^= hash + 0x9e3779b9u + (seed << 6) + (seed >> 2);
}

struct ComparisonKeyHash {
  std::size_t operator()(const ComparisonKey &key) const {
    std::size_t result = 0;
    hashCombine(result, key.ruleId);
    hashCombine(result, key.first);
    hashCombine(result, key.second);
    hashCombine(result, key.firstEnd);
    hashCombine(result, key.secondEnd);
    hashCombine(result, key.firstAtRoot);
    hashCombine(result, key.secondAtRoot);
    hashCombine(result, key.firstRule6Ref);
    hashCombine(result, key.secondRule6Ref);
    return result;
  }
};

struct ComparisonSessionState {
  std::size_t depth = 0;
  using ComparisonMap =
      std::unordered_map<ComparisonKey, int, ComparisonKeyHash>;
  ComparisonMap currentResults;
  ComparisonMap previousResults;
  ComparisonMap currentConstitutionalResults;
  ComparisonMap previousConstitutionalResults;

  struct SortKey {
    std::uint64_t sortId;
    const Node *node;
    bool nodeAtRoot;
    std::uint8_t outgoingMask;
    const Atom *rule6Ref;
    bool deep;
    std::uint8_t inputCount;
    std::array<const Edge *, SequenceRule::MAX_CACHED_SORT_EDGES> input;

    bool operator==(const SortKey &other) const {
      return sortId == other.sortId && node == other.node &&
             nodeAtRoot == other.nodeAtRoot &&
             outgoingMask == other.outgoingMask && rule6Ref == other.rule6Ref &&
             deep == other.deep && inputCount == other.inputCount &&
             std::equal(input.begin(), input.begin() + inputCount,
                        other.input.begin());
    }
  };

  struct SortKeyHash {
    std::size_t operator()(const SortKey &key) const {
      std::size_t result = 0;
      hashCombine(result, key.sortId);
      hashCombine(result, key.node);
      hashCombine(result, key.nodeAtRoot);
      hashCombine(result, key.outgoingMask);
      hashCombine(result, key.rule6Ref);
      hashCombine(result, key.deep);
      for (std::size_t i = 0; i < key.inputCount; ++i) {
        hashCombine(result, key.input[i]);
      }
      return result;
    }
  };

  struct SortValue {
    std::array<Edge *, SequenceRule::MAX_CACHED_SORT_EDGES> sorted;
    std::uint8_t sortedCount;
    Priority priority;
  };

  using SortMap = std::unordered_map<SortKey, SortValue, SortKeyHash>;
  SortMap currentSorts;
  SortMap previousSorts;
  SortMap currentConstitutionalSorts;
  SortMap previousConstitutionalSorts;
};

// The caches are accelerators, not part of the result. Two bounded generations
// retain recent hot entries instead of permanently freezing the first entries
// encountered by a large cage traversal. The limits apply independently to
// constitutional and descriptor-dependent entries; the latter are normally
// short-lived because each auxiliary distance sphere invalidates them.
constexpr std::size_t MAX_CACHED_COMPARISONS_PER_FAMILY = 250000;
constexpr std::size_t COMPARISONS_PER_GENERATION =
    MAX_CACHED_COMPARISONS_PER_FAMILY / 2;
constexpr std::size_t MAX_CACHED_SORTS_PER_FAMILY = 50000;
constexpr std::size_t SORTS_PER_GENERATION = MAX_CACHED_SORTS_PER_FAMILY / 2;

thread_local ComparisonSessionState comparisonSessionState;

std::uint64_t nextCacheId() {
  static std::atomic<std::uint64_t> nextId{0};
  return ++nextId;
}

std::uint8_t getOutgoingMask(const Node *node, std::span<Edge *const> edges) {
  std::uint8_t result = 0;
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (edges[i]->isBeg(node)) {
      result |= static_cast<std::uint8_t>(1u << i);
    }
  }
  return result;
}

}  // namespace

SequenceRule::ComparisonSession::ComparisonSession() {
  if (comparisonSessionState.depth++ == 0) {
    comparisonSessionState.currentResults.clear();
    comparisonSessionState.previousResults.clear();
    comparisonSessionState.currentConstitutionalResults.clear();
    comparisonSessionState.previousConstitutionalResults.clear();
    comparisonSessionState.currentSorts.clear();
    comparisonSessionState.previousSorts.clear();
    comparisonSessionState.currentConstitutionalSorts.clear();
    comparisonSessionState.previousConstitutionalSorts.clear();
  }
}

SequenceRule::ComparisonSession::~ComparisonSession() {
  if (--comparisonSessionState.depth == 0) {
    comparisonSessionState.currentResults.clear();
    comparisonSessionState.previousResults.clear();
    comparisonSessionState.currentConstitutionalResults.clear();
    comparisonSessionState.previousConstitutionalResults.clear();
    comparisonSessionState.currentSorts.clear();
    comparisonSessionState.previousSorts.clear();
    comparisonSessionState.currentConstitutionalSorts.clear();
    comparisonSessionState.previousConstitutionalSorts.clear();
  }
}

void SequenceRule::ComparisonSession::invalidateAuxiliaryCaches() {
  comparisonSessionState.currentResults.clear();
  comparisonSessionState.previousResults.clear();
  comparisonSessionState.currentSorts.clear();
  comparisonSessionState.previousSorts.clear();
}

SequenceRule::SequenceRule(bool useConstitutionalRootEquivalence,
                           bool auxiliaryIndependent)
    : d_cacheId{nextCacheId()},
      d_useConstitutionalRootEquivalence{useConstitutionalRootEquivalence},
      d_auxiliaryIndependent{auxiliaryIndependent} {
  dp_sorter.reset(new Sort(this));
}

SequenceRule::~SequenceRule() = default;

Descriptor SequenceRule::getBondLabel(const Edge *edge) const {
  Bond *bond = edge->getBond();
  if (bond == nullptr) {
    return Descriptor::NONE;
  }
  return edge->getAux();
}

int SequenceRule::getComparision(const Edge *a, const Edge *b) const {
  return getComparision(a, b, true);
}

int SequenceRule::getComparision(const Edge *a, const Edge *b,
                                 bool deep) const {
  return deep ? recursiveCompare(a, b) : compare(a, b);
}

const Sort *SequenceRule::getSorter() const { return dp_sorter.get(); }

bool SequenceRule::isAuxiliaryIndependent() const {
  return d_auxiliaryIndependent;
}

int SequenceRule::recursiveCompare(const Edge *a, const Edge *b) const {
  const ComparisonSession comparisonSession;

  if (!CIPLabeler_detail::decrementRemainingCallCountAndCheck()) {
    throw MaxIterationsExceeded();
  }
  if (ControlCHandler::getGotSignal()) {
    throw ControlCCaught();
  }

  // Rules 1a, 1b, and 2 usually distinguish ligands using an inexpensive
  // comparison at the current depth. Their nonzero direct results are not
  // cached, so avoid constructing and probing a cache key unless they tie.
  if (d_auxiliaryIndependent) {
    const auto directResult = compare(a, b);
    if (directResult != 0) {
      return directResult;
    }
    if (a->getEnd()->isTerminal() && b->getEnd()->isTerminal()) {
      return 0;
    }
  }

  const auto firstGraph = a->getBeg()->getDigraph();
  const auto secondGraph = b->getBeg()->getDigraph();
  const bool firstAtRoot = firstGraph->getCurrentRoot() == a->getBeg();
  const bool secondAtRoot = secondGraph->getCurrentRoot() == b->getBeg();
  const ComparisonKey key{
      d_cacheId,
      a,
      b,
      a->getEnd(),
      b->getEnd(),
      firstAtRoot,
      secondAtRoot,
      d_auxiliaryIndependent ? nullptr : firstGraph->getRule6Ref(),
      d_auxiliaryIndependent ? nullptr : secondGraph->getRule6Ref()};
  auto &currentResults =
      d_auxiliaryIndependent
          ? comparisonSessionState.currentConstitutionalResults
          : comparisonSessionState.currentResults;
  auto &previousResults =
      d_auxiliaryIndependent
          ? comparisonSessionState.previousConstitutionalResults
          : comparisonSessionState.previousResults;
  const auto cacheResult = [&](int result) {
    if (currentResults.size() >= COMPARISONS_PER_GENERATION) {
      previousResults = std::move(currentResults);
      currentResults.clear();
    }
    currentResults.emplace(key, result);
    return result;
  };
  const auto findResult = [](const auto &results,
                             const ComparisonKey &candidate, int &value) {
    const auto found = results.find(candidate);
    if (found == results.end()) {
      return false;
    }
    value = found->second;
    return true;
  };
  int cachedResult = 0;
  if (findResult(currentResults, key, cachedResult) ||
      findResult(previousResults, key, cachedResult)) {
    return cachedResult;
  }
  const ComparisonKey reverseKey{
      d_cacheId,
      b,
      a,
      b->getEnd(),
      a->getEnd(),
      secondAtRoot,
      firstAtRoot,
      d_auxiliaryIndependent ? nullptr : secondGraph->getRule6Ref(),
      d_auxiliaryIndependent ? nullptr : firstGraph->getRule6Ref()};
  if (firstGraph == secondGraph) {
    if (findResult(currentResults, reverseKey, cachedResult) ||
        findResult(previousResults, reverseKey, cachedResult)) {
      return -cachedResult;
    }
  }

  if (!d_auxiliaryIndependent) {
    const auto directResult = compare(a, b);
    if (directResult != 0) {
      // Direct Rule 4b/5 comparisons can include complete pair-list
      // traversals.
      return cacheResult(directResult);
    }
    if (a->getEnd()->isTerminal() && b->getEnd()->isTerminal()) {
      return cacheResult(0);
    }
  }

  if (!isRecursiveComparisonNeeded(a, b) ||
      hasEquivalentAcyclicContinuation(a, b) ||
      hasEquivalentConstitutionalSiblingContinuation(a, b)) {
    return cacheResult(0);
  }

  return cacheResult(recursiveCompareEqual(a, b));
}

bool SequenceRule::getCachedSort(std::uint64_t sortId, const Node *node,
                                 bool deep, bool auxiliaryIndependent,
                                 std::vector<Edge *> &edges, bool &unique,
                                 bool &pseudoAsymmetric) {
  if (node == nullptr || edges.size() < 2u ||
      edges.size() > SequenceRule::MAX_CACHED_SORT_EDGES) {
    return false;
  }
  const auto graph = node->getDigraph();
  const auto edgeSpan = std::span<Edge *const>{edges.data(), edges.size()};
  ComparisonSessionState::SortKey key{
      sortId,
      node,
      graph->getCurrentRoot() == node,
      getOutgoingMask(node, edgeSpan),
      auxiliaryIndependent ? nullptr : graph->getRule6Ref(),
      deep,
      static_cast<std::uint8_t>(edges.size()),
      {}};
  std::copy(edges.begin(), edges.end(), key.input.begin());
  const auto findSort = [&](const auto &sorts) {
    const auto found = sorts.find(key);
    if (found == sorts.end()) {
      return false;
    }
    edges.assign(found->second.sorted.begin(),
                 found->second.sorted.begin() + found->second.sortedCount);
    unique = found->second.priority.isUnique();
    pseudoAsymmetric = found->second.priority.isPseudoAsymetric();
    return true;
  };
  const auto &currentSorts =
      auxiliaryIndependent ? comparisonSessionState.currentConstitutionalSorts
                           : comparisonSessionState.currentSorts;
  const auto &previousSorts =
      auxiliaryIndependent ? comparisonSessionState.previousConstitutionalSorts
                           : comparisonSessionState.previousSorts;
  return findSort(currentSorts) || findSort(previousSorts);
}

void SequenceRule::cacheSort(std::uint64_t sortId, const Node *node, bool deep,
                             bool auxiliaryIndependent,
                             std::span<Edge *const> input,
                             const std::vector<Edge *> &sorted,
                             const Priority &priority) {
  if (node == nullptr || input.size() < 2u || sorted.size() != input.size() ||
      input.size() > SequenceRule::MAX_CACHED_SORT_EDGES) {
    return;
  }
  auto &currentSorts = auxiliaryIndependent
                           ? comparisonSessionState.currentConstitutionalSorts
                           : comparisonSessionState.currentSorts;
  auto &previousSorts = auxiliaryIndependent
                            ? comparisonSessionState.previousConstitutionalSorts
                            : comparisonSessionState.previousSorts;
  if (currentSorts.size() >= SORTS_PER_GENERATION) {
    previousSorts = std::move(currentSorts);
    currentSorts.clear();
  }
  const auto graph = node->getDigraph();
  ComparisonSessionState::SortKey key{
      sortId,
      node,
      graph->getCurrentRoot() == node,
      getOutgoingMask(node, input),
      auxiliaryIndependent ? nullptr : graph->getRule6Ref(),
      deep,
      static_cast<std::uint8_t>(input.size()),
      {}};
  std::copy(input.begin(), input.end(), key.input.begin());
  ComparisonSessionState::SortValue value{
      {}, static_cast<std::uint8_t>(sorted.size()), priority};
  std::copy(sorted.begin(), sorted.end(), value.sorted.begin());
  currentSorts.emplace(std::move(key), std::move(value));
}

int SequenceRule::recursiveCompareEqual(const Edge *a, const Edge *b) const {
  std::vector<std::pair<const Edge *, const Edge *>> queue{{a, b}};
  std::vector<Edge *> as;
  std::vector<Edge *> bs;
  as.reserve(4);
  bs.reserve(4);

  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    const auto [aParent, bParent] = queue[pos];
    if (hasEquivalentAcyclicContinuation(aParent, bParent)) {
      continue;
    }
    auto aNode = aParent->getEnd();
    auto bNode = bParent->getEnd();
    if (aNode->isTerminal() && bNode->isTerminal()) {
      continue;
    }
    const auto &aNodeEdges = aNode->getEdges();
    const auto &bNodeEdges = bNode->getEdges();
    as.assign(aNodeEdges.begin(), aNodeEdges.end());
    bs.assign(bNodeEdges.begin(), bNodeEdges.end());

    // shallow sort first of all
    sort(aNode, as, false);
    sort(bNode, bs, false);

    int sizediff = three_way_comparison(as.size(), bs.size());

    {
      auto aIt = as.begin();
      auto bIt = bs.begin();
      for (; aIt != as.end() && bIt != bs.end(); ++aIt, ++bIt) {
        Edge *aEdge = *aIt;
        Edge *bEdge = *bIt;

        if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
          continue;
        }

        const auto cmp = compare(aEdge, bEdge);
        if (cmp != 0) {
          return cmp;
        }
      }
    }

    if (sizediff != 0) {
      return sizediff;
    }

    sort(aNode, as);
    sort(bNode, bs);

    {
      auto aIt = as.begin();
      auto bIt = bs.begin();
      for (; aIt != as.end() && bIt != bs.end(); ++aIt, ++bIt) {
        Edge *aEdge = *aIt;
        Edge *bEdge = *bIt;

        if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
          continue;
        }

        const auto cmp = compare(aEdge, bEdge);
        if (cmp != 0) {
          return cmp;
        }

        queue.emplace_back(aEdge, bEdge);
      }
    }
  }
  return 0;
}

void SequenceRule::setSorter(const Sort *sorter) {
  dp_sorter.reset(sorter);
  d_cacheId = nextCacheId();
}

bool SequenceRule::isRecursiveComparisonNeeded(const Edge *a,
                                               const Edge *b) const {
  (void)a;
  (void)b;
  return true;
}

Priority SequenceRule::sort(const Node *node, std::vector<Edge *> &edges,
                            bool deep) const {
  return getSorter()->prioritize(node, edges, deep);
}

Priority SequenceRule::sort(const Node *node,
                            std::vector<Edge *> &edges) const {
  return sort(node, edges, true);
}

bool SequenceRule::areUpEdges(Node *aNode, Node *bNode, Edge *aEdge,
                              Edge *bEdge) const {
  // step over 'up' edges
  if (aEdge->isEnd(aNode)) {
    // if b is 'down' something's not right!
    if (!bEdge->isEnd(bNode)) {
      throw std::runtime_error("Something unexpected!");
    }
    return true;
  }
  return false;
}

bool SequenceRule::hasEquivalentAcyclicContinuation(const Edge *a,
                                                    const Edge *b) const {
  const auto aBeg = a->getBeg();
  const auto bBeg = b->getBeg();
  const auto aEnd = a->getEnd();
  const auto bEnd = b->getEnd();
  const auto graph = aBeg->getDigraph();

  if (graph != bBeg->getDigraph() || a->getBond() == nullptr ||
      a->getBond() != b->getBond() || aBeg->getAtom() != bBeg->getAtom() ||
      aEnd->getAtom() == nullptr || aEnd->getAtom() != bEnd->getAtom() ||
      aEnd->isDuplicateOrH() || bEnd->isDuplicateOrH() ||
      aEnd->getDistance() != bEnd->getDistance() ||
      !aEnd->isOriginalChildOf(aBeg) || !bEnd->isOriginalChildOf(bBeg)) {
    return false;
  }

  // Constitutional rules do not inspect auxiliary descriptors. Once two
  // occurrences enter the same side of the same bridge at the same distance,
  // their remaining unfolded trees are identical even when that side contains
  // stereochemical configurations. Descriptor-dependent rules retain the
  // stricter configuration-free proof below.
  if (d_auxiliaryIndependent && !graph->getMol().isInRing(a->getBond())) {
    return true;
  }
  return graph->isAcyclicBranchWithoutConfiguration(a);
}

bool SequenceRule::hasEquivalentConstitutionalSiblingContinuation(
    const Edge *a, const Edge *b) const {
  if (!d_useConstitutionalRootEquivalence) {
    return false;
  }

  const auto aBeg = a->getBeg();
  const auto bBeg = b->getBeg();
  const auto aEnd = a->getEnd();
  const auto bEnd = b->getEnd();
  const auto graph = aBeg->getDigraph();
  const auto currentRoot = graph->getCurrentRoot();
  if (graph != bBeg->getDigraph() || aBeg != bBeg ||
      a->getBond() == nullptr || b->getBond() == nullptr ||
      aEnd->isDuplicateOrH() || bEnd->isDuplicateOrH() ||
      !aEnd->isOriginalChildOf(aBeg) || !bEnd->isOriginalChildOf(aBeg)) {
    return false;
  }

  // Keep the existing current-root behavior. Below that root, restrict exact
  // symmetry searches to two ring directions: acyclic continuations already
  // have a cheaper bridge proof, and this gate leaves ordinary pendant-tree
  // comparisons unchanged.
  if (aBeg != currentRoot &&
      (!graph->getMol().isInRing(a->getBond()) ||
       !graph->getMol().isInRing(b->getBond()))) {
    return false;
  }
  bool equivalent = false;
  std::vector<unsigned int> movedAtoms;
  if (aBeg == graph->getOriginalRoot()) {
    equivalent = graph->getMol().hasConstitutionalAutomorphism(
        aBeg->getAtom(), aEnd->getAtom(), bEnd->getAtom(), movedAtoms);
  } else {
    // Every non-original occurrence carries an immutable path from the
    // original root. Requiring that path to be fixed pointwise preserves its
    // visited set and every distance used for ring duplicates under Rule 1b,
    // whether or not this local sibling pair is below a temporary reroot.
    equivalent = graph->getMol().hasConstitutionalAutomorphism(
        aBeg->getAtom(), aEnd->getAtom(), bEnd->getAtom(),
        aBeg->getVisitedAtomCheckpoint(), aBeg->getVisitedAtomDeltas(),
        movedAtoms);
  }
  if (equivalent) {
    // VF2 may return a valid witness that unnecessarily composes the local
    // ligand swap with an independent symmetry elsewhere in the component.
    // Such a witness makes labelAux() conservatively discover remote stereo
    // configurations. When that happened, ask for another exact witness
    // while fixing every other configuration neighborhood pointwise. Failure
    // only keeps the original witness and the existing conservative fallback.
    const auto configurationOwner = graph->getOriginalRoot()->getAtom();
    bool preservesConfigurations =
        !graph->getMol().constitutionalAutomorphismMovesConfiguration(
            movedAtoms, configurationOwner);
    if (!preservesConfigurations) {
      std::vector<unsigned int> configurationPreservingMovedAtoms;
      if (graph->getMol()
              .hasConfigurationPreservingConstitutionalAutomorphism(
                  aBeg->getAtom(), aEnd->getAtom(), bEnd->getAtom(),
                  aBeg->getVisitedAtomCheckpoint(),
                  aBeg->getVisitedAtomDeltas(), configurationOwner,
                  configurationPreservingMovedAtoms)) {
        movedAtoms = std::move(configurationPreservingMovedAtoms);
        preservesConfigurations = true;
      }
    }
    if (preservesConfigurations && aBeg == currentRoot) {
      // The compared edges are ligands of the configuration currently being
      // ranked. Their exact rooted isomorphism also preserves every possible
      // auxiliary annotation. The configuration-specific labeling code can
      // therefore avoid auxiliary discovery when Rule 6 cannot break the
      // remaining tie.
      graph->noteAuxiliaryInvariantRootTie();
    }
    graph->noteConstitutionalRootEquivalence(movedAtoms);
  }
  return equivalent;
}

}  // namespace CIPLabeler
}  // namespace RDKit
