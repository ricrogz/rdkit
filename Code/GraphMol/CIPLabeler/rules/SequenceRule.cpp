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

#include <cstddef>
#include <functional>
#include <memory>
#include <unordered_map>
#include <utility>

#include "SequenceRule.h"

#include "../CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

namespace {

struct EdgePair {
  const Edge *first;
  const Edge *second;

  bool operator==(const EdgePair &other) const {
    return first == other.first && second == other.second;
  }
};

struct EdgePairHash {
  std::size_t operator()(const EdgePair &pair) const {
    const auto firstHash = std::hash<const Edge *>{}(pair.first);
    const auto secondHash = std::hash<const Edge *>{}(pair.second);
    return firstHash ^
           (secondHash + 0x9e3779b9u + (firstHash << 6) + (firstHash >> 2));
  }
};

struct RecursiveComparisonCache {
  std::size_t depth = 0;
  std::unordered_map<EdgePair, int, EdgePairHash> results;
};

// Digraph orientation, auxiliary labels, and the Rule 6 reference can change
// between top-level comparisons. Keep results only while one exact rule object
// has an active recursive comparison. This also prevents temporary Rule4b and
// Rule5New replacement rules from inheriting entries when stack addresses are
// reused with a different reference descriptor.
thread_local std::unordered_map<const SequenceRule *,
                                std::unique_ptr<RecursiveComparisonCache>>
    recursiveComparisonCaches;

class ScopedRecursiveComparisonCache {
 public:
  explicit ScopedRecursiveComparisonCache(const SequenceRule *rule)
      : d_rule{rule} {
    auto &cache = recursiveComparisonCaches[rule];
    if (!cache) {
      cache = std::make_unique<RecursiveComparisonCache>();
    }
    d_cache = cache.get();
    if (d_cache->depth++ == 0) {
      d_cache->results.clear();
    }
  }

  ScopedRecursiveComparisonCache(const ScopedRecursiveComparisonCache &) =
      delete;
  ScopedRecursiveComparisonCache &operator=(
      const ScopedRecursiveComparisonCache &) = delete;

  ~ScopedRecursiveComparisonCache() {
    if (--d_cache->depth == 0) {
      recursiveComparisonCaches.erase(d_rule);
    }
  }

  bool find(const Edge *first, const Edge *second, int &result) const {
    const auto found = d_cache->results.find({first, second});
    if (found == d_cache->results.end()) {
      return false;
    }
    result = found->second;
    return true;
  }

  void store(const Edge *first, const Edge *second, int result) {
    d_cache->results.emplace(EdgePair{first, second}, result);
  }

 private:
  const SequenceRule *d_rule;
  RecursiveComparisonCache *d_cache;
};

}  // namespace

SequenceRule::SequenceRule() : dp_sorter{new Sort(this)} {}

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
  if (!CIPLabeler_detail::decrementRemainingCallCountAndCheck()) {
    throw MaxIterationsExceeded();
  }
  if (ControlCHandler::getGotSignal()) {
    throw ControlCCaught();
  }

  const auto directResult = compare(a, b);
  if (directResult != 0) {
    return directResult;
  }

  ScopedRecursiveComparisonCache cache(this);
  int cachedResult;
  if (cache.find(a, b, cachedResult)) {
    return cachedResult;
  }

  const auto result = recursiveCompareEqual(a, b);
  cache.store(a, b, result);
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

void SequenceRule::setSorter(const Sort *sorter) { dp_sorter.reset(sorter); }

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

}  // namespace CIPLabeler
}  // namespace RDKit
