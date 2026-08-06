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

#include <algorithm>
#include <utility>

#include <RDGeneral/Invariant.h>

#include "Rule4b.h"

#include "../Digraph.h"
#include "Pairlist.h"

namespace RDKit {
namespace CIPLabeler {

Rule4b::Rule4b() = default;

Rule4b::Rule4b(Descriptor ref) : d_ref{ref} {}

std::vector<Descriptor> Rule4b::getReferenceDescriptors(
    const Node *node) const {
  std::vector<Descriptor> result;
  result.reserve(2);
  auto prev = initialLevel(node);
  while (!prev.empty()) {
    for (const auto &nodes : prev) {
      if (getReference(nodes, result)) {
        return result;
      }
    }
    prev = getNextLevel(prev);
  }
  return {};
}

int Rule4b::compare(const Edge *a, const Edge *b) const {
  const auto &aBeg = a->getBeg();
  const auto &aEnd = a->getEnd();
  const auto &bBeg = b->getBeg();
  const auto &bEnd = b->getEnd();
  if (aBeg->getDigraph()->getCurrentRoot() != aBeg ||
      bBeg->getDigraph()->getCurrentRoot() != bBeg) {
    if (d_ref == Descriptor::NONE) {
      return 0;
    }
    Descriptor aDesc = aEnd->getAux();
    Descriptor bDesc = bEnd->getAux();
    if (aDesc != Descriptor::NONE && bDesc != Descriptor::NONE &&
        aDesc != Descriptor::ns && bDesc != Descriptor::ns) {
      bool alike = PairList::ref(d_ref) == PairList::ref(aDesc);
      bool blike = PairList::ref(d_ref) == PairList::ref(bDesc);
      if (alike && !blike) {
        return +1;
      }
      if (blike && !alike) {
        return -1;
      }
    }
    return 0;
  } else {
    auto list1 = newPairLists(getReferenceDescriptors(aEnd));

    auto list2 = newPairLists(getReferenceDescriptors(bEnd));

    if (list1.empty() != list2.empty()) {
      throw std::runtime_error(
          "Substituents should be topologically equivalent!");
    }
    if (list1.size() == 1) {
      return comparePairs(aEnd, bEnd, list1[0].getRefDescriptor(),
                          list2[0].getRefDescriptor());
    } else if (list1.size() > 1) {
      std::vector<const Node *> queue;
      std::vector<Edge *> edges;
      edges.reserve(4);
      for (auto &plist : list1) {
        fillPairs(aEnd, plist, queue, edges);
      }
      for (auto &plist : list2) {
        fillPairs(bEnd, plist, queue, edges);
      }

      std::sort(list1.rbegin(), list1.rend());
      std::sort(list2.rbegin(), list2.rend());

      for (auto i = 0u; i < list1.size(); ++i) {
        int cmp = list1[i].compareTo(list2[i]);
        if (cmp != 0) {
          return cmp;
        }
      }
    }
    return 0;
  }
}

bool Rule4b::getReference(const std::vector<const Node *> &nodes,
                          std::vector<Descriptor> &result) const {
  int right = 0;
  int left = 0;
  for (const auto &node : nodes) {
    auto desc = node->getAux();
    switch (desc) {
      case Descriptor::NONE:
        continue;
      case Descriptor::R:
      case Descriptor::M:
      case Descriptor::seqCis:
        ++right;
        break;
      case Descriptor::S:
      case Descriptor::P:
      case Descriptor::seqTrans:
        ++left;
        break;
      default:
        break;
    }
  }
  if (right + left == 0) {
    return false;
  } else if (right > left) {
    result.push_back(Descriptor::R);
    return true;
  } else if (right < left) {
    result.push_back(Descriptor::S);
    return true;
  } else {
    result.push_back(Descriptor::R);
    result.push_back(Descriptor::S);
    return true;
  }
}

std::vector<std::vector<const Node *>> Rule4b::initialLevel(
    const Node *node) const {
  return {{node}};
}

std::vector<std::vector<const Node *>> Rule4b::getNextLevel(
    const std::vector<std::vector<const Node *>> &prevLevel) const {
  std::vector<std::vector<const Node *>> nextLevel;
  nextLevel.reserve(4 * prevLevel.size());

  for (const auto &prev : prevLevel) {
    std::vector<std::vector<std::vector<Edge *>>> tmp;
    tmp.reserve(prev.size());
    for (const auto &node : prev) {
      auto edges = node->getNonTerminalOutEdges();
      sort(node, edges);
      tmp.push_back(getSorter()->getGroups(edges));
    }

    // Equivalent nodes must produce the same number of priority groups.
    // Check each node's groups, not just the first node repeatedly.
    std::size_t size = 0;
    if (!tmp.empty()) {
      size = tmp.front().size();
      for (const auto &groups : tmp) {
        if (size != groups.size()) {
          throw std::runtime_error("Something unexpected!");
        }
      }
    }

    for (std::size_t i = 0; i < size; ++i) {
      std::vector<const Node *> eq;
      std::size_t eqSize = 0;
      for (const auto &groups : tmp) {
        eqSize += groups[i].size();
      }
      eq.reserve(eqSize);
      for (const auto &groups : tmp) {
        for (const auto edge : groups[i]) {
          eq.push_back(edge->getEnd());
        }
      }
      if (!eq.empty()) {
        nextLevel.push_back(eq);
      }
    }
  }
  return nextLevel;
}

std::vector<PairList> Rule4b::newPairLists(
    const std::vector<Descriptor> &descriptors) const {
  std::vector<PairList> pairs;
  pairs.reserve(descriptors.size());
  for (Descriptor descriptor : descriptors) {
    pairs.emplace_back(descriptor);
  }
  return pairs;
}

void Rule4b::fillPairs(const Node *beg, PairList &plist,
                       std::vector<const Node *> &queue,
                       std::vector<Edge *> &edges) const {
  const Rule4b replacement_rule(plist.getRefDescriptor());
  const auto &sorter = getRefSorter(&replacement_rule);
  queue.clear();
  queue.push_back(beg);

  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    const auto node = queue[pos];
    plist.add(node->getAux());
    const auto &nodeEdges = node->getEdges();
    edges.assign(nodeEdges.begin(), nodeEdges.end());
    sorter.prioritize(node, edges);
    for (const auto &edge : edges) {
      if (edge->isBeg(node) && !edge->getEnd()->isTerminal()) {
        queue.push_back(edge->getEnd());
      }
    }
  }
}

int Rule4b::comparePairs(const Node *a, const Node *b, Descriptor refA,
                         Descriptor refB) const {
  const Rule4b replacementA(refA);
  const Rule4b replacementB(refB);
  const auto &aSorter = getRefSorter(&replacementA);
  const auto &bSorter = getRefSorter(&replacementB);
  std::vector<std::pair<const Node *, const Node *>> queue{{a, b}};
  std::vector<Edge *> aEdges;
  std::vector<Edge *> bEdges;
  aEdges.reserve(4);
  bEdges.reserve(4);

  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    const auto [aNode, bNode] = queue[pos];

    const auto &desA = PairList::ref(aNode->getAux());
    const auto &desB = PairList::ref(bNode->getAux());

    if (desA == refA && desB != refB) {
      return +1;
    } else if (desA != refA && desB == refB) {
      return -1;
    }

    const auto &aNodeEdges = aNode->getEdges();
    const auto &bNodeEdges = bNode->getEdges();
    aEdges.assign(aNodeEdges.begin(), aNodeEdges.end());
    bEdges.assign(bNodeEdges.begin(), bNodeEdges.end());
    aSorter.prioritize(aNode, aEdges);
    bSorter.prioritize(bNode, bEdges);

    auto aIt = aEdges.begin();
    auto bIt = bEdges.begin();
    while (aIt != aEdges.end() && bIt != bEdges.end()) {
      while (aIt != aEdges.end() &&
             (!(*aIt)->isBeg(aNode) || (*aIt)->getEnd()->isTerminal())) {
        ++aIt;
      }
      while (bIt != bEdges.end() &&
             (!(*bIt)->isBeg(bNode) || (*bIt)->getEnd()->isTerminal())) {
        ++bIt;
      }
      if (aIt == aEdges.end() || bIt == bEdges.end()) {
        break;
      }
      queue.emplace_back((*aIt)->getEnd(), (*bIt)->getEnd());
      ++aIt;
      ++bIt;
    }
  }
  return 0;
}

Sort Rule4b::getRefSorter(const SequenceRule *replacement_rule) const {
  const auto &rules = getSorter()->getRules();

  CHECK_INVARIANT(std::find(rules.begin(), rules.end(), this) != rules.end(),
                  "Rule4b instance not in rule set");

  std::vector<const SequenceRule *> new_rules;
  new_rules.reserve(rules.size());
  for (const auto &rule : rules) {
    if (this != rule) {
      new_rules.push_back(rule);
    }
  }
  new_rules.push_back(replacement_rule);
  return {new_rules};
}

}  // namespace CIPLabeler
}  // namespace RDKit
