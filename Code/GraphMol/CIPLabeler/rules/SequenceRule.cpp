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

#include <utility>

#include "SequenceRule.h"

#include "../CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

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

  int cmp = compare(a, b);
  if (cmp != 0) {
    return cmp;
  }

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

        cmp = compare(aEdge, bEdge);
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

        cmp = compare(aEdge, bEdge);
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
