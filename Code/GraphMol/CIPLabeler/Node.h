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
#include <memory>
#include <span>
#include <vector>

#include "Descriptor.h"
#include "Edge.h"
#include "Mancude.h"

namespace RDKit {

class Atom;

namespace CIPLabeler {

class Digraph;

// Immutable-path visit state used while unfolding the molecular graph. States
// requiring at most two 64-bit words remain entirely inline. Larger states
// share a dense checkpoint and retain at most seven added atoms locally; the
// eighth addition creates a new checkpoint. Copies therefore have bounded cost
// and never copy the molecule-sized bitset.
class NodeVisitState {
 public:
  NodeVisitState() = default;
  explicit NodeVisitState(std::size_t wordCount);
  NodeVisitState(std::vector<std::uint64_t> &&words);
  NodeVisitState(const NodeVisitState &other);
  NodeVisitState &operator=(const NodeVisitState &other);
  NodeVisitState(NodeVisitState &&other) noexcept;
  NodeVisitState &operator=(NodeVisitState &&other) noexcept;

  bool empty() const noexcept;
  bool test(unsigned int atomIdx) const;
  void set(unsigned int atomIdx);

  // Materializes pending additions into a private immutable checkpoint when
  // necessary. The returned span remains valid until this state is modified.
  std::span<const std::uint64_t> words() const;

  // A non-materializing view used by exact symmetry queries. Together these
  // two spans represent the same set as words().
  std::span<const std::uint64_t> checkpointWords() const noexcept;
  std::span<const unsigned int> addedAtoms() const noexcept;

 private:
  static constexpr std::size_t INLINE_WORD_COUNT = 2;
  static constexpr std::size_t CHECKPOINT_INTERVAL = 8;

  struct LargeState {
    std::shared_ptr<const std::vector<std::uint64_t>> checkpoint;
    std::array<unsigned int, CHECKPOINT_INTERVAL> addedAtoms{};
    std::uint8_t addedAtomCount = 0;
  };

  std::size_t d_wordCount = 0;
  std::array<std::uint64_t, INLINE_WORD_COUNT> d_inlineWords{};
  std::unique_ptr<LargeState> dp_large;

  void materialize() const;
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
  std::span<const std::uint64_t> getVisitedAtomCheckpoint() const;
  std::span<const unsigned int> getVisitedAtomDeltas() const;

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
  friend class Digraph;

  Digraph *dp_g;
  Atom *dp_atom;
  const Node *dp_parent;
  Edge *dp_parent_edge = nullptr;
  unsigned int d_tree_depth = 0;
  // Permanent membership in the immutable original-parent tree. This is
  // independent of the edge directions around the current temporary root.
  bool d_attached_to_origin = false;
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
