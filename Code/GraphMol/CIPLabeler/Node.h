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
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/container/small_vector.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "Descriptor.h"
#include "Mancude.h"
#include "Edge.h"

namespace RDKit {

class Atom;

namespace CIPLabeler {

class Digraph;

// Most CIP molecules need only one 64-bit path word. Keeping it inline avoids
// one heap allocation for every nonterminal occurrence without making the
// many terminal Node objects unnecessarily large. The converting constructor
// preserves source compatibility with the former std::vector-based API.
class NodeVisitState : public boost::container::small_vector<std::uint64_t, 1> {
  using Base = boost::container::small_vector<std::uint64_t, 1>;

 public:
  using Base::Base;
  NodeVisitState() = default;
  NodeVisitState(std::vector<std::uint64_t> &&words)
      : Base(words.begin(), words.end()) {}
};

class Node {
 public:
  /**
   * Flag indicates whether the node has been expanded.
   */
  static const int EXPANDED = 0x1;

  /**
   * Flag indicates whether the node was duplicated
   * at a ring closure.
   */
  static const int RING_DUPLICATE = 0x2;

  /**
   * Flag indicates whether the node was duplicated
   * at a bond with order &gt; 1.
   */
  static const int BOND_DUPLICATE = 0x4;

  /**
   * Mask to check if a node is duplicated.
   */

  static const int DUPLICATE = RING_DUPLICATE | BOND_DUPLICATE;

  /**
   * Node was created for an implicit hydrogen,
   * the 'atom' value will be null.
   */
  static const int IMPL_HYDROGEN = 0x8;

  /**
   * Mask to check if a node is duplicated or created for an implicit H (not a
   * primary node).
   */
  static const int DUPLICATE_OR_H =
      RING_DUPLICATE | BOND_DUPLICATE | IMPL_HYDROGEN;

  Node() = delete;
  Node(const Node &) = delete;
  Node &operator=(const Node &) = delete;

  Node(Digraph *g, NodeVisitState &&visit, Atom *atom,
       boost::rational<int> &&frac, int dist, int flags, const Node *parent);

  Digraph *getDigraph() const;

  Atom *getAtom() const;

  unsigned int getAtomIdx() const;

  int getDistance() const;

  boost::rational<int> getAtomicNumFraction() const;

  int getAtomicNum() const;

  unsigned getMassNum() const;

  double getAtomicMass() const;

  Descriptor getAux() const;

  bool isSet(int mask) const;

  bool isDuplicate() const;

  bool isDuplicateOrH() const;

  bool isTerminal() const;

  bool isExpanded() const;

  bool isVisited(int idx) const;

  // True when this node was created as a child of parent before any temporary
  // digraph rerooting changed edge directions.
  bool isOriginalChildOf(const Node *parent) const;

  // Atoms in the immutable occurrence path built from the original root.
  // Temporary rerooting changes edge directions but not this history.
  std::span<const std::uint64_t> getVisitedAtoms() const;

  Node *newChild(int idx, Atom *atom) const;

  Node *newBondDuplicateChild(int idx, Atom *atom) const;

  Node *newRingDuplicateChild(int idx, Atom *atom) const;

  Node *newImplicitHydrogenChild() const;

  void add(Edge *e);

  void setAux(Descriptor desc);

  // Number of effective auxiliary descriptors in this node's immutable
  // original-forward occurrence subtree. The mask is a combination of the
  // AUX_DESCRIPTOR_* constants from Descriptor.h.
  std::size_t getAuxDescriptorCount(unsigned mask) const;

  // Adjust descriptor counts for this node and each of its original
  // ancestors. Used by Node and Edge descriptor assignment.
  void adjustAuxDescriptorCount(unsigned descriptorClass, int delta);

  const std::vector<Edge *> &getEdges() const;

  std::vector<Edge *> getEdges(Atom *end) const;

  std::vector<Edge *> getNonTerminalOutEdges() const;

 private:
  Digraph *dp_g;
  Atom *dp_atom;
  const Node *dp_parent;
  int d_dist;
  boost::rational<int> d_atomic_num;
  double d_atomic_mass;
  Descriptor d_aux = Descriptor::NONE;
  std::array<std::size_t, 3> d_aux_descriptor_counts{};
  int d_flags = 0x0;

  std::vector<Edge *> d_edges;

  NodeVisitState d_visit;

  Node *newTerminalChild(int idx, Atom *atom, int flags) const;
  int getVisitedDistance(int idx) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
