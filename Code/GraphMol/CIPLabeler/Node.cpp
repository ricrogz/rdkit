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
#include <bit>
#include <utility>
#include <vector>

#include <RDGeneral/Invariant.h>

#include "Digraph.h"
#include "Edge.h"
#include "Node.h"
#include "CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

NodeVisitState::NodeVisitState(std::size_t wordCount) : d_wordCount{wordCount} {
  if (d_wordCount > INLINE_WORD_COUNT) {
    dp_large = std::make_unique<LargeState>();
    dp_large->checkpoint = std::make_shared<const std::vector<std::uint64_t>>(
        d_wordCount, std::uint64_t{0});
  }
}

NodeVisitState::NodeVisitState(std::vector<std::uint64_t> &&words)
    : d_wordCount{words.size()} {
  if (d_wordCount <= INLINE_WORD_COUNT) {
    std::copy(words.begin(), words.end(), d_inlineWords.begin());
  } else {
    dp_large = std::make_unique<LargeState>();
    dp_large->checkpoint =
        std::make_shared<const std::vector<std::uint64_t>>(std::move(words));
  }
}

NodeVisitState::NodeVisitState(const NodeVisitState &other)
    : d_wordCount{other.d_wordCount}, d_inlineWords{other.d_inlineWords} {
  if (other.dp_large != nullptr) {
    dp_large = std::make_unique<LargeState>(*other.dp_large);
  }
}

NodeVisitState &NodeVisitState::operator=(const NodeVisitState &other) {
  if (this != &other) {
    std::unique_ptr<LargeState> large;
    if (other.dp_large != nullptr) {
      large = std::make_unique<LargeState>(*other.dp_large);
    }
    d_wordCount = other.d_wordCount;
    d_inlineWords = other.d_inlineWords;
    dp_large = std::move(large);
  }
  return *this;
}

NodeVisitState::NodeVisitState(NodeVisitState &&other) noexcept
    : d_wordCount{std::exchange(other.d_wordCount, std::size_t{0})},
      d_inlineWords{other.d_inlineWords},
      dp_large{std::move(other.dp_large)} {
  other.d_inlineWords = {};
}

NodeVisitState &NodeVisitState::operator=(NodeVisitState &&other) noexcept {
  if (this != &other) {
    d_wordCount = std::exchange(other.d_wordCount, std::size_t{0});
    d_inlineWords = other.d_inlineWords;
    other.d_inlineWords = {};
    dp_large = std::move(other.dp_large);
  }
  return *this;
}

bool NodeVisitState::empty() const noexcept { return d_wordCount == 0u; }

bool NodeVisitState::test(unsigned int atomIdx) const {
  const auto word = atomIdx / 64u;
  if (word >= d_wordCount) {
    return false;
  }
  const auto mask = std::uint64_t{1} << (atomIdx % 64u);
  if (d_wordCount <= INLINE_WORD_COUNT) {
    return (d_inlineWords[word] & mask) != 0u;
  }
  const auto &large = *dp_large;
  if (((*large.checkpoint)[word] & mask) != 0u) {
    return true;
  }
  for (std::size_t i = 0; i < large.addedAtomCount; ++i) {
    if (large.addedAtoms[i] == atomIdx) {
      return true;
    }
  }
  return false;
}

void NodeVisitState::set(unsigned int atomIdx) {
  const auto word = atomIdx / 64u;
  PRECONDITION(word < d_wordCount, "visit atom index is out of range")
  if (d_wordCount <= INLINE_WORD_COUNT) {
    d_inlineWords[word] |= std::uint64_t{1} << (atomIdx % 64u);
    return;
  }
  if (test(atomIdx)) {
    return;
  }

  auto &large = *dp_large;
  if (large.addedAtomCount == CHECKPOINT_INTERVAL) {
    materialize();
  }
  PRECONDITION(large.addedAtomCount < CHECKPOINT_INTERVAL,
               "visit-state delta count overflow")
  large.addedAtoms[large.addedAtomCount++] = atomIdx;
  if (large.addedAtomCount == CHECKPOINT_INTERVAL) {
    materialize();
  }
}

void NodeVisitState::materialize() const {
  if (dp_large == nullptr || dp_large->addedAtomCount == 0u) {
    return;
  }

  auto &large = *dp_large;
  auto checkpoint =
      std::make_shared<std::vector<std::uint64_t>>(*large.checkpoint);
  for (std::size_t i = 0; i < large.addedAtomCount; ++i) {
    const auto atomIdx = large.addedAtoms[i];
    (*checkpoint)[atomIdx / 64u] |= std::uint64_t{1} << (atomIdx % 64u);
  }
  large.checkpoint = std::move(checkpoint);
  large.addedAtomCount = 0u;
}

std::span<const std::uint64_t> NodeVisitState::words() const {
  if (d_wordCount == 0u) {
    return {};
  }
  if (d_wordCount <= INLINE_WORD_COUNT) {
    return {d_inlineWords.data(), d_wordCount};
  }
  materialize();
  return {dp_large->checkpoint->data(), dp_large->checkpoint->size()};
}

std::span<const std::uint64_t> NodeVisitState::checkpointWords()
    const noexcept {
  if (d_wordCount == 0u) {
    return {};
  }
  if (d_wordCount <= INLINE_WORD_COUNT) {
    return {d_inlineWords.data(), d_wordCount};
  }
  return {dp_large->checkpoint->data(), dp_large->checkpoint->size()};
}

std::span<const unsigned int> NodeVisitState::addedAtoms() const noexcept {
  if (dp_large == nullptr) {
    return {};
  }
  return {dp_large->addedAtoms.data(), dp_large->addedAtomCount};
}

Node *Node::newTerminalChild(int idx, Atom *atom, int flags) const {
  int new_dist = flags & DUPLICATE ? getVisitedDistance(idx) : d_dist + 1;
  NodeVisitState new_visit;

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

Node::Node(Digraph *g, NodeVisitState &&visit, Atom *atom,
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
  return d_visit.test(atom_idx);
}

bool Node::isOriginalChildOf(const Node *parent) const {
  return dp_parent == parent;
}

std::span<const std::uint64_t> Node::getVisitedAtoms() const {
  return d_visit.words();
}

std::span<const std::uint64_t> Node::getVisitedAtomCheckpoint() const {
  return d_visit.checkpointWords();
}

std::span<const unsigned int> Node::getVisitedAtomDeltas() const {
  return d_visit.addedAtoms();
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
  new_visit.set(atom_idx);
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
