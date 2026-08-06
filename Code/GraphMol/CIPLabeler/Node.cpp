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
#include <bit>
#include <vector>

#include "Digraph.h"
#include "Edge.h"
#include "Node.h"
#include "CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

Node *Node::newTerminalChild(int idx, Atom *atom, int flags) const {
  int new_dist = flags & DUPLICATE ? getVisitedDistance(idx) : d_dist + 1;
  std::vector<std::uint64_t> new_visit;

  if (flags & BOND_DUPLICATE) {
    const auto &frac = dp_g->getMol().getFractionalAtomicNum(dp_atom);
    if (frac.isAveraged()) {
      return &dp_g->addNode(std::move(new_visit), atom, frac.value(), new_dist,
                            flags, this);
    }
  }

  auto atomic_num = atom ? atom->getAtomicNum() : 1;
  return &dp_g->addNode(std::move(new_visit), atom, atomic_num, new_dist, flags,
                        this);
}

Node::Node(Digraph *g, std::vector<std::uint64_t> &&visit, Atom *atom,
           boost::rational<int> &&frac, int dist, int flags, const Node *parent)
    : dp_g{g},
      dp_atom{atom},
      dp_parent{parent},
      d_dist{dist},
      d_atomic_num{std::move(frac)},
      d_flags{flags},
      d_visit{std::move(visit)} {
  d_edges.reserve(d_visit.empty() ? 1u : 4u);
  if (d_flags & DUPLICATE) {
    d_atomic_mass = 0.;
  } else if (dp_atom != nullptr) {
    d_atomic_mass = dp_g->getMol().getAtomicMass(dp_atom);
  } else {
    const auto &table = RDKit::PeriodicTable::getTable();
    d_atomic_mass = table->getAtomicWeight(1);
  }
  if (d_visit.empty() || d_flags & DUPLICATE) {
    d_flags |= EXPANDED;
  }
}

Digraph *Node::getDigraph() const { return dp_g; }

Atom *Node::getAtom() const { return dp_atom; }

unsigned int Node::getAtomIdx() const {
  if (isSet(IMPL_HYDROGEN)) {
    return Atom::NOATOM;
  }

  return dp_atom->getIdx();
}

int Node::getDistance() const { return d_dist; }

boost::rational<int> Node::getAtomicNumFraction() const { return d_atomic_num; }

int Node::getAtomicNum() const {
  if (dp_atom == nullptr) {
    return 1;
  }
  return dp_atom->getAtomicNum();
};

unsigned Node::getMassNum() const {
  if (dp_atom == nullptr || isDuplicate()) {
    return 0u;
  }
  return dp_atom->getIsotope();
}

double Node::getAtomicMass() const { return d_atomic_mass; }

Descriptor Node::getAux() const { return d_aux; }

bool Node::isSet(int mask) const { return mask & d_flags; }

bool Node::isDuplicate() const { return d_flags & DUPLICATE; }

bool Node::isDuplicateOrH() const { return d_flags & DUPLICATE_OR_H; }

bool Node::isTerminal() const {
  return d_visit.empty() || (isExpanded() && d_edges.size() == 1);
}

bool Node::isExpanded() const { return d_flags & EXPANDED; }

bool Node::isVisited(int idx) const {
  const auto atom_idx = static_cast<unsigned int>(idx);
  return d_visit[atom_idx / 64u] & (std::uint64_t{1} << (atom_idx % 64u));
}

bool Node::isOriginalChildOf(const Node *parent) const {
  return dp_parent == parent;
}

int Node::getVisitedDistance(int idx) const {
  if (!isVisited(idx)) {
    // A forward bond duplicate refers to the new child, which is not yet in
    // this node's path. The previous visit array also returned zero here.
    return 0;
  }
  for (const auto *node = this; node != nullptr; node = node->dp_parent) {
    if (node->dp_atom != nullptr &&
        node->dp_atom->getIdx() == static_cast<unsigned int>(idx)) {
      return node->d_dist;
    }
  }
  return 0;
}

Node *Node::newChild(int idx, Atom *atom) const {
  auto new_visit = d_visit;
  const auto atom_idx = static_cast<unsigned int>(idx);
  new_visit[atom_idx / 64u] |= std::uint64_t{1} << (atom_idx % 64u);
  auto atomic_num = atom ? atom->getAtomicNum() : 1;
  return &dp_g->addNode(std::move(new_visit), atom, atomic_num, d_dist + 1, 0,
                        this);
}

Node *Node::newBondDuplicateChild(int idx, Atom *atom) const {
  return newTerminalChild(idx, atom, BOND_DUPLICATE);
}

Node *Node::newRingDuplicateChild(int idx, Atom *atom) const {
  return newTerminalChild(idx, atom, RING_DUPLICATE);
}

Node *Node::newImplicitHydrogenChild() const {
  return newTerminalChild(-1, nullptr, IMPL_HYDROGEN);
}

void Node::add(Edge *e) { d_edges.push_back(e); }

void Node::setAux(Descriptor desc) {
  const auto oldClass = getAuxDescriptorClass(d_aux);
  const auto newClass = getAuxDescriptorClass(desc);
  if (oldClass != newClass) {
    if (oldClass != 0u) {
      adjustAuxDescriptorCount(oldClass, -1);
    }
    if (newClass != 0u) {
      adjustAuxDescriptorCount(newClass, 1);
    }
  }
  d_aux = desc;
}

std::size_t Node::getAuxDescriptorCount(unsigned mask) const {
  std::size_t result = 0;
  for (unsigned i = 0; i < d_aux_descriptor_counts.size(); ++i) {
    if ((mask & (1u << i)) != 0u) {
      result += d_aux_descriptor_counts[i];
    }
  }
  return result;
}

void Node::adjustAuxDescriptorCount(unsigned descriptorClass, int delta) {
  PRECONDITION(
      descriptorClass != 0u && (descriptorClass & (descriptorClass - 1u)) == 0u,
      "descriptor class must contain exactly one bit")
  const auto index = static_cast<unsigned>(std::countr_zero(descriptorClass));
  PRECONDITION(index < d_aux_descriptor_counts.size(),
               "invalid descriptor class")

  for (auto node = this; node != nullptr;
       node = const_cast<Node *>(node->dp_parent)) {
    if (delta > 0) {
      ++node->d_aux_descriptor_counts[index];
    } else {
      PRECONDITION(node->d_aux_descriptor_counts[index] != 0u,
                   "auxiliary descriptor count underflow")
      --node->d_aux_descriptor_counts[index];
    }
  }
}

const std::vector<Edge *> &Node::getEdges() const {
  if (!isExpanded()) {
    auto non_const_this = const_cast<Node *>(this);
    non_const_this->d_flags |= EXPANDED;
    dp_g->expand(non_const_this);
  }
  return d_edges;
}

std::vector<Edge *> Node::getEdges(Atom *end) const {
  std::vector<Edge *> res;
  for (auto &edge : getEdges()) {
    if (edge->getEnd()->isDuplicate()) {
      continue;
    };
    if (end == edge->getBeg()->getAtom() || end == edge->getEnd()->getAtom()) {
      res.push_back(edge);
    }
  }
  return res;
}

std::vector<Edge *> Node::getNonTerminalOutEdges() const {
  std::vector<Edge *> edges;
  for (auto &edge : getEdges()) {
    if (edge->isBeg(this) && !edge->getEnd()->isTerminal()) {
      edges.push_back(edge);
    }
  }
  return edges;
}

}  // namespace CIPLabeler
}  // namespace RDKit
