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
#include <stdexcept>

#include <RDGeneral/Invariant.h>

#include "Edge.h"
#include "Node.h"

namespace RDKit {
namespace CIPLabeler {

Edge::Edge(Node *beg, Node *end, Bond *bond)
    : dp_beg{beg}, dp_end{end}, dp_bond{bond} {}

Node *Edge::getOther(const Node *node) const {
  PRECONDITION(node, "bad node")

  if (isBeg(node)) {
    return getEnd();
  } else if (isEnd(node)) {
    return getBeg();
  } else {
    throw std::runtime_error("Not an end-point of this edge!");
  }
}

Node *Edge::getBeg() const { return dp_beg; }

Node *Edge::getEnd() const { return dp_end; }

Bond *Edge::getBond() const { return dp_bond; }

Descriptor Edge::getAux() const { return d_aux; }

bool Edge::isBeg(const Node *node) const { return node == dp_beg; }

bool Edge::isEnd(const Node *node) const { return node == dp_end; }

void Edge::setAux(Descriptor aux) {
  const auto oldClass = getAuxDescriptorClass(d_aux);
  const auto newClass = getAuxDescriptorClass(aux);
  if (oldClass != newClass) {
    Node *originalChild = nullptr;
    if (dp_end->isOriginalChildOf(dp_beg)) {
      originalChild = dp_end;
    } else if (dp_beg->isOriginalChildOf(dp_end)) {
      originalChild = dp_beg;
    }
    if (originalChild != nullptr) {
      if (oldClass != 0u) {
        originalChild->adjustAuxDescriptorCount(oldClass, -1);
      }
      if (newClass != 0u) {
        originalChild->adjustAuxDescriptorCount(newClass, 1);
      }
    }
  }
  d_aux = aux;
}

void Edge::flip() { std::swap(dp_beg, dp_end); }

}  // namespace CIPLabeler
}  // namespace RDKit
