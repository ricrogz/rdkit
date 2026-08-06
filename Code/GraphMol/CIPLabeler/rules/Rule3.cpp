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

#include "Rule3.h"

#include "../Digraph.h"

namespace RDKit {
namespace CIPLabeler {

namespace {
int ord(Descriptor lab) {
  switch (lab) {
    case Descriptor::E:
      return 1;
    case Descriptor::Z:
      return 2;
    default:
      return 0;
  }
}
}  // namespace

Rule3::Rule3() = default;

int Rule3::compare(const Edge *a, const Edge *b) const {
  return three_way_comparison(ord(a->getEnd()->getAux()),
                              ord(b->getEnd()->getAux()));
}

bool Rule3::isRecursiveComparisonNeeded(const Edge *a, const Edge *b) const {
  return a->getBeg()->getDigraph()->hasAuxDescriptorOnSide(
             a, AUX_DESCRIPTOR_GEOMETRIC) ||
         b->getBeg()->getDigraph()->hasAuxDescriptorOnSide(
             b, AUX_DESCRIPTOR_GEOMETRIC);
}

}  // namespace CIPLabeler
}  // namespace RDKit
