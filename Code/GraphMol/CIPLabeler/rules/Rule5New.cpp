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
#include <RDGeneral/Invariant.h>

#include "Rule5New.h"

#include "../Digraph.h"

namespace RDKit {
namespace CIPLabeler {

Rule5New::Rule5New() = default;

Rule5New::Rule5New(Descriptor ref) : d_ref{ref} {}

int Rule5New::compare(const Edge *a, const Edge *b) const {
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
    auto listRA = PairList(Descriptor::R);
    auto listRB = PairList(Descriptor::R);
    auto listSA = PairList(Descriptor::S);
    auto listSB = PairList(Descriptor::S);
    std::vector<const Node *> queue;
    std::vector<Edge *> edges;
    edges.reserve(4);
    fillPairs(aEnd, listRA, queue, edges);
    fillPairs(aEnd, listSA, queue, edges);
    fillPairs(bEnd, listRB, queue, edges);
    fillPairs(bEnd, listSB, queue, edges);
    int cmpR = listRA.compareTo(listRB);
    int cmpS = listSA.compareTo(listSB);
    // -2/+2 for pseudo-asymetric
    // -1/+1 if not (e.g. the R > R and S > S lists)
    if (cmpR < 0) {
      return cmpS < 0 ? -1 : -2;
    } else if (cmpR > 0) {
      return cmpS > 0 ? +1 : +2;
    } else {
      return 0;
    }
  }
}

void Rule5New::fillPairs(const Node *beg, PairList &plist,
                         std::vector<const Node *> &queue,
                         std::vector<Edge *> &edges) const {
  const Rule5New replacement_rule(plist.getRefDescriptor());
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

Sort Rule5New::getRefSorter(const SequenceRule *replacement_rule) const {
  const auto &rules = getSorter()->getRules();

  CHECK_INVARIANT(std::find(rules.begin(), rules.end(), this) != rules.end(),
                  "Rule5New instance not in rule set");

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
