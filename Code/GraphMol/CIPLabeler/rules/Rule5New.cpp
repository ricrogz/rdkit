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
#include <memory>
#include <vector>

#include <RDGeneral/Invariant.h>

#include "Rule5New.h"

#include "../Digraph.h"

namespace RDKit {
namespace CIPLabeler {

static Rule5New referenceR(Descriptor::R);
static Rule5New referenceS(Descriptor::S);
std::unique_ptr<const Sort> referenceSorterR;
std::unique_ptr<const Sort> referenceSorterS;

Rule5New::Rule5New() = default;

Rule5New::Rule5New(Descriptor ref) : d_ref{ref} {}

int8_t Rule5New::compare(const Edge *a, const Edge *b) const {
  if (!a->getBeg()->getDigraph()->hasEffectiveAuxDescriptors() &&
      !b->getBeg()->getDigraph()->hasEffectiveAuxDescriptors()) {
    return 0;
  }
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
    EdgeVector edges;
    fillPairs(aEnd, listRA, queue, edges);
    fillPairs(aEnd, listSA, queue, edges);
    fillPairs(bEnd, listRB, queue, edges);
    fillPairs(bEnd, listSB, queue, edges);
    auto cmpR = listRA.compareTo(listRB);
    auto cmpS = listSA.compareTo(listSB);
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
                         EdgeVector &edges) const {
  const auto &sorter = getRefSorter(plist.getRefDescriptor());
  queue.clear();
  queue.push_back(beg);

  for (unsigned int pos = 0; pos < queue.size(); ++pos) {
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

const Sort &Rule5New::getRefSorter(Descriptor ref) const {
  if (ref == Descriptor::R) {
    if (!referenceSorterR) {
      referenceSorterR = makeRefSorter(&referenceR);
    }
    return *referenceSorterR;
  }
  if (ref == Descriptor::S) {
    if (!referenceSorterS) {
      referenceSorterS = makeRefSorter(&referenceS);
    }
    return *referenceSorterS;
  }
  throw std::logic_error("Invalid Rule 5 reference descriptor");
}

std::unique_ptr<const Sort> Rule5New::makeRefSorter(
    const SequenceRule *replacementRule) const {
  const auto &rules = getSorter()->getRules();
  std::vector<const SequenceRule *> new_rules;
  new_rules.reserve(rules.size());
  for (const auto &rule : rules) {
    if (this != rule) {
      new_rules.push_back(rule);
    }
  }
  new_rules.push_back(replacementRule);
  return std::make_unique<const Sort>(std::move(new_rules));
}

}  // namespace CIPLabeler
}  // namespace RDKit
