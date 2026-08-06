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
  const Node *firstRoot;
  const Node *secondRoot;
  const Atom *firstRule6Ref;
  const Atom *secondRule6Ref;

  bool operator==(const ComparisonKey &other) const {
    return ruleId == other.ruleId && first == other.first &&
           second == other.second && firstRoot == other.firstRoot &&
           secondRoot == other.secondRoot &&
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
    hashCombine(result, key.firstRoot);
    hashCombine(result, key.secondRoot);
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

  static constexpr std::size_t MAX_INLINE_SORT_EDGES = 8;

  struct SortKey {
    std::uint64_t sortId;
    const Node *node;
    const Node *root;
    const Atom *rule6Ref;
    bool deep;
    std::uint8_t inputCount;
    std::array<const Edge *, MAX_INLINE_SORT_EDGES> input;

    bool operator==(const SortKey &other) const {
      return sortId == other.sortId && node == other.node &&
             root == other.root && rule6Ref == other.rule6Ref &&
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
      hashCombine(result, key.root);
      hashCombine(result, key.rule6Ref);
      hashCombine(result, key.deep);
      for (std::size_t i = 0; i < key.inputCount; ++i) {
        hashCombine(result, key.input[i]);
      }
      return result;
    }
  };

  struct SortValue {
    std::array<Edge *, MAX_INLINE_SORT_EDGES> sorted;
    std::uint8_t sortedCount;
    Priority priority;
  };

  using SortMap = std::unordered_map<SortKey, SortValue, SortKeyHash>;
  SortMap currentSorts;
  SortMap previousSorts;
};

// The caches are accelerators, not part of the result. Two bounded generations
// retain recent hot entries instead of permanently freezing the first entries
// encountered by a large cage traversal.
constexpr std::size_t MAX_CACHED_COMPARISONS = 250000;
constexpr std::size_t COMPARISONS_PER_GENERATION = MAX_CACHED_COMPARISONS / 2;
constexpr std::size_t MAX_CACHED_SORTS = 50000;
constexpr std::size_t SORTS_PER_GENERATION = MAX_CACHED_SORTS / 2;

thread_local ComparisonSessionState comparisonSessionState;

std::uint64_t nextCacheId() {
  static std::atomic<std::uint64_t> nextId{0};
  return ++nextId;
}

}  // namespace

SequenceRule::ComparisonSession::ComparisonSession() {
  if (comparisonSessionState.depth++ == 0) {
    comparisonSessionState.currentResults.clear();
    comparisonSessionState.previousResults.clear();
    comparisonSessionState.currentSorts.clear();
    comparisonSessionState.previousSorts.clear();
  }
}

SequenceRule::ComparisonSession::~ComparisonSession() {
  if (--comparisonSessionState.depth == 0) {
    comparisonSessionState.currentResults.clear();
    comparisonSessionState.previousResults.clear();
    comparisonSessionState.currentSorts.clear();
    comparisonSessionState.previousSorts.clear();
  }
}

SequenceRule::SequenceRule(bool useConstitutionalRootEquivalence)
    : dp_sorter{new Sort(this)},
      d_cacheId{nextCacheId()},
      d_useConstitutionalRootEquivalence{useConstitutionalRootEquivalence} {}

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

int SequenceRule::recursiveCompare(const Edge *a, const Edge *b) const {
  const ComparisonSession comparisonSession;

  if (!CIPLabeler_detail::decrementRemainingCallCountAndCheck()) {
    throw MaxIterationsExceeded();
  }
  if (ControlCHandler::getGotSignal()) {
    throw ControlCCaught();
  }

  const auto firstGraph = a->getBeg()->getDigraph();
  const auto secondGraph = b->getBeg()->getDigraph();
  const ComparisonKey key{d_cacheId,
                          a,
                          b,
                          firstGraph->getCurrentRoot(),
                          secondGraph->getCurrentRoot(),
                          firstGraph->getRule6Ref(),
                          secondGraph->getRule6Ref()};
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
  if (findResult(comparisonSessionState.currentResults, key, cachedResult) ||
      findResult(comparisonSessionState.previousResults, key, cachedResult)) {
    return cachedResult;
  }
  const ComparisonKey reverseKey{d_cacheId,
                                 b,
                                 a,
                                 secondGraph->getCurrentRoot(),
                                 firstGraph->getCurrentRoot(),
                                 secondGraph->getRule6Ref(),
                                 firstGraph->getRule6Ref()};
  if (firstGraph == secondGraph) {
    if (findResult(comparisonSessionState.currentResults, reverseKey,
                   cachedResult) ||
        findResult(comparisonSessionState.previousResults, reverseKey,
                   cachedResult)) {
      return -cachedResult;
    }
  }

  const auto directResult = compare(a, b);
  if (directResult != 0) {
    return directResult;
  }

  if (!isRecursiveComparisonNeeded(a, b) ||
      hasEquivalentAcyclicContinuation(a, b) ||
      hasEquivalentConstitutionalRootContinuation(a, b)) {
    return 0;
  }

  const auto result = recursiveCompareEqual(a, b);
  if (comparisonSessionState.currentResults.size() >=
      COMPARISONS_PER_GENERATION) {
    comparisonSessionState.previousResults =
        std::move(comparisonSessionState.currentResults);
    comparisonSessionState.currentResults.clear();
  }
  comparisonSessionState.currentResults.emplace(key, result);
  return result;
}

bool SequenceRule::getCachedSort(std::uint64_t sortId, const Node *node,
                                 bool deep, std::vector<Edge *> &edges,
                                 bool &unique, bool &pseudoAsymmetric) {
  if (node == nullptr || edges.size() < 2u ||
      edges.size() > ComparisonSessionState::MAX_INLINE_SORT_EDGES) {
    return false;
  }
  const auto graph = node->getDigraph();
  ComparisonSessionState::SortKey key{sortId,
                                      node,
                                      graph->getCurrentRoot(),
                                      graph->getRule6Ref(),
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
  return findSort(comparisonSessionState.currentSorts) ||
         findSort(comparisonSessionState.previousSorts);
}

void SequenceRule::cacheSort(std::uint64_t sortId, const Node *node, bool deep,
                             const std::vector<Edge *> &input,
                             const std::vector<Edge *> &sorted,
                             const Priority &priority) {
  if (node == nullptr || input.size() < 2u || sorted.size() != input.size() ||
      input.size() > ComparisonSessionState::MAX_INLINE_SORT_EDGES) {
    return;
  }
  if (comparisonSessionState.currentSorts.size() >= SORTS_PER_GENERATION) {
    comparisonSessionState.previousSorts =
        std::move(comparisonSessionState.currentSorts);
    comparisonSessionState.currentSorts.clear();
  }
  const auto graph = node->getDigraph();
  ComparisonSessionState::SortKey key{sortId,
                                      node,
                                      graph->getCurrentRoot(),
                                      graph->getRule6Ref(),
                                      deep,
                                      static_cast<std::uint8_t>(input.size()),
                                      {}};
  std::copy(input.begin(), input.end(), key.input.begin());
  ComparisonSessionState::SortValue value{
      {}, static_cast<std::uint8_t>(sorted.size()), priority};
  std::copy(sorted.begin(), sorted.end(), value.sorted.begin());
  comparisonSessionState.currentSorts.emplace(std::move(key), std::move(value));
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

  return graph->isAcyclicBranchWithoutConfiguration(a);
}

bool SequenceRule::hasEquivalentConstitutionalRootContinuation(
    const Edge *a, const Edge *b) const {
  if (!d_useConstitutionalRootEquivalence) {
    return false;
  }

  const auto aBeg = a->getBeg();
  const auto bBeg = b->getBeg();
  const auto aEnd = a->getEnd();
  const auto bEnd = b->getEnd();
  const auto graph = aBeg->getDigraph();
  const auto root = graph->getOriginalRoot();
  if (graph != bBeg->getDigraph() || graph->getCurrentRoot() != root ||
      aBeg != root || bBeg != root || a->getBond() == nullptr ||
      b->getBond() == nullptr || aEnd->isDuplicateOrH() ||
      bEnd->isDuplicateOrH() || !aEnd->isOriginalChildOf(root) ||
      !bEnd->isOriginalChildOf(root)) {
    return false;
  }

  // Keep the exact-isomorphism cost focused on cyclic ligand systems, where
  // path unfolding is combinatorial and the cheaper bridge shortcut cannot
  // apply. Returning false only disables this optimization.
  if (!graph->getMol().isInRing(a->getBond()) ||
      !graph->getMol().isInRing(b->getBond())) {
    return false;
  }

  const bool equivalent = graph->getMol().hasConstitutionalAutomorphism(
      root->getAtom(), aEnd->getAtom(), bEnd->getAtom());
  if (equivalent) {
    graph->noteConstitutionalRootEquivalence();
  }
  return equivalent;
}

}  // namespace CIPLabeler
}  // namespace RDKit
