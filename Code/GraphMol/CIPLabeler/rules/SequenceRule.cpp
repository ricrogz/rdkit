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
  std::unordered_map<ComparisonKey, int, ComparisonKeyHash> results;
};

// The cache is an accelerator, not part of the result. Stop adding entries on
// exceptionally large digraphs so a labeling session has bounded cache memory.
constexpr std::size_t MAX_CACHED_COMPARISONS = 250000;

thread_local ComparisonSessionState comparisonSessionState;

std::uint64_t nextCacheId() {
  static std::atomic<std::uint64_t> nextId{0};
  return ++nextId;
}

}  // namespace

SequenceRule::ComparisonSession::ComparisonSession() {
  if (comparisonSessionState.depth++ == 0) {
    comparisonSessionState.results.clear();
  }
}

SequenceRule::ComparisonSession::~ComparisonSession() {
  if (--comparisonSessionState.depth == 0) {
    comparisonSessionState.results.clear();
  }
}

SequenceRule::SequenceRule()
    : dp_sorter{new Sort(this)}, d_cacheId{nextCacheId()} {}

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
  const auto found = comparisonSessionState.results.find(key);
  if (found != comparisonSessionState.results.end()) {
    return found->second;
  }

  const auto directResult = compare(a, b);
  if (directResult != 0) {
    return directResult;
  }

  if (hasEquivalentAcyclicContinuation(a, b)) {
    return 0;
  }

  const auto result = recursiveCompareEqual(a, b);
  if (comparisonSessionState.results.size() < MAX_CACHED_COMPARISONS) {
    comparisonSessionState.results.emplace(key, result);
  }
  return result;
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

}  // namespace CIPLabeler
}  // namespace RDKit
