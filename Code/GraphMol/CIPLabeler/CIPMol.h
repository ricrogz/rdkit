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

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "Descriptor.h"
#include "Mancude.h"

namespace RDKit {

namespace CIPLabeler {

template <typename T, typename U>
class CIPMolSpan {
 public:
  class CIPMolIter {
   public:
    CIPMolIter() = delete;
    CIPMolIter(ROMol &mol, U pos) : d_mol{mol}, d_pos{std::move(pos)} {}

    T &operator*() {
      d_current = d_mol[*d_pos];
      return d_current;
    }

    CIPMolIter &operator++() {
      ++d_pos;
      return *this;
    }

    bool operator!=(const CIPMolIter &it) const { return d_pos != it.d_pos; }

   private:
    ROMol &d_mol;
    U d_pos;
    T d_current = nullptr;
  };

 public:
  CIPMolSpan() = delete;
  CIPMolSpan(ROMol &mol, std::pair<U, U> &&itr)
      : d_mol{mol},
        d_istart{std::move(itr.first)},
        d_iend{std::move(itr.second)} {}

  CIPMolIter begin() { return {d_mol, d_istart}; }
  CIPMolIter end() { return {d_mol, d_iend}; }

 private:
  ROMol &d_mol;
  const U d_istart;
  const U d_iend;
};

class CIPMol {
 public:
  CIPMol() = delete;

  explicit CIPMol(ROMol &mol);

  // Average atomic number with other atoms that are in an
  // aromatic ring with this one.
  const FractionalAtomicNum &getFractionalAtomicNum(Atom *atom) const;

  unsigned getNumAtoms() const;

  unsigned getNumBonds() const;

  Atom *getAtom(int idx) const;

  CXXAtomIterator<MolGraph, Atom *> atoms() const;

  Bond *getBond(int idx) const;

  CIPMolSpan<Bond *, ROMol::OEDGE_ITER> getBonds(Atom *atom) const;

  CIPMolSpan<Atom *, ROMol::ADJ_ITER> getNeighbors(Atom *atom) const;

  bool isInRing(Bond *bond) const;

  // Record every atom that is the focus of a configuration considered by the
  // current labeling operation. This lets the lazy digraph prove that an
  // acyclic branch cannot contain auxiliary stereochemical information.
  void setConfigurationFoci(boost::dynamic_bitset<> foci);

  // Returns true only when bond is acyclic and the molecular component on the
  // endAtom side contains no registered configuration focus.
  bool isAcyclicBranchWithoutConfiguration(Bond *bond, Atom *endAtom) const;

  // Exact equality proof for two constitutional ligands at a fixed molecular
  // root. This is only an equality shortcut; automorphism indices are never
  // used as CIP priorities.
  bool hasConstitutionalAutomorphism(Atom *root, Atom *from, Atom *to) const;

  // As above, but require every atom set in the path bitset to remain fixed.
  // This proves equality for a rerooted occurrence without losing the
  // path-history state that controls ring-duplicate construction.
  bool hasConstitutionalAutomorphism(
      Atom *root, Atom *from, Atom *to,
      std::span<const std::uint64_t> fixedAtoms,
      std::span<const unsigned int> addedFixedAtoms = {}) const;

  // Component membership is used to keep symmetry and auxiliary-label work
  // local to the connected structure containing the center of interest.
  bool isInSameComponent(Atom *first, Atom *second) const;

  // Integer bond order of a kekulized molecule
  // Dative bonds get bond order 0.
  int getBondOrder(Bond *bond) const;

  double getAtomicMass(Atom *atom) const;

 private:
  ROMol &d_mol;
  mutable std::vector<RDKit::Bond::BondType> d_kekulized_bonds;
  std::vector<FractionalAtomicNum> d_atomnums;
  std::vector<RDKit::Bond *> d_bonds;
  mutable std::vector<double> d_atomic_masses;
  mutable std::vector<unsigned char> d_atomic_mass_cached;
  boost::dynamic_bitset<> d_configuration_foci;
  // Two directed sides per bond: 0=unknown, 1=contains a focus, 2=no focus.
  mutable std::vector<unsigned char> d_configuration_branch_cache;
  struct ConstitutionalAutomorphismKey {
    unsigned int root;
    unsigned int first;
    unsigned int second;

    bool operator==(const ConstitutionalAutomorphismKey &other) const {
      return root == other.root && first == other.first &&
             second == other.second;
    }
  };
  struct ConstitutionalAutomorphismKeyHash {
    std::size_t operator()(const ConstitutionalAutomorphismKey &key) const;
  };
  struct ConstitutionalAutomorphismEvidence {
    // A witness mapping is represented by the sorted atoms it moves. It
    // remains a proof for any path disjoint from that set. Sparse sets avoid
    // making cache entries proportional to the size of a large molecule.
    std::vector<std::vector<unsigned int>> movedAtomSets;
    // An exhaustive failure while fixing F remains a proof for every superset
    // of F. Both collections are maintained as small dominance antichains.
    std::vector<std::vector<unsigned int>> failedFixedAtomSets;
    std::size_t storedAtomIndexCount = 0;
    // A bounded search exhausted its compatibility-callback allowance.
    // Existing witnesses can still be reused, but new searches for this
    // endpoint pair are disabled.
    bool searchDisabled = false;
  };
  using ConstitutionalAutomorphismEvidenceMap =
      std::unordered_map<ConstitutionalAutomorphismKey,
                         ConstitutionalAutomorphismEvidence,
                         ConstitutionalAutomorphismKeyHash>;
  mutable ConstitutionalAutomorphismEvidenceMap
      d_current_constitutional_automorphism_evidence;
  mutable ConstitutionalAutomorphismEvidenceMap
      d_previous_constitutional_automorphism_evidence;
  mutable std::size_t d_current_automorphism_evidence_atom_indices = 0;
  mutable std::size_t d_previous_automorphism_evidence_atom_indices = 0;
  mutable bool d_constitutional_automorphism_data_initialized = false;
  mutable std::vector<unsigned int> d_component_ids;
  mutable std::vector<unsigned int> d_component_local_indices;
  mutable std::vector<std::vector<unsigned int>> d_component_atoms;
  // The exact self-isomorphism search is performed only on the connected
  // component containing the compared ligands. Components are copied lazily,
  // so disconnected spectators neither consume VF2 work nor add memory unless
  // they are themselves queried.
  mutable std::vector<std::shared_ptr<RWMol>> d_component_mols;
  mutable std::vector<std::vector<unsigned int>> d_component_bond_indices;
  mutable std::vector<unsigned int> d_constitutional_symmetry_classes;
  mutable std::vector<unsigned char>
      d_constitutional_symmetry_classes_initialized;
  mutable std::vector<unsigned char> d_component_has_cycle;
  mutable std::vector<unsigned char> d_component_is_supported;
  // Candidate-compatibility callbacks spent by constrained self-matches in
  // each component. A cumulative bound prevents many individually bounded
  // endpoint queries from becoming an unbounded aggregate cost.
  mutable std::vector<std::size_t> d_component_automorphism_search_callbacks;
  // Linear distance-prefilter passes are independently bounded so negative
  // endpoint pairs cannot accumulate quadratic preprocessing.
  mutable std::vector<std::size_t> d_component_automorphism_prefilter_queries;
  mutable std::vector<unsigned char> d_atom_in_cyclic_core;
  mutable std::vector<unsigned int> d_cyclic_core_distances;

  bool hasUniqueBond(Atom *begin, Atom *end) const;
  void initKekulizedBonds() const;
  void initConstitutionalAutomorphismData() const;
  void initConstitutionalSymmetryClasses(unsigned int component) const;
  const RWMol &getConstitutionalComponentMol(unsigned int component) const;
  bool atomsConstitutionallyEquivalent(const Atom &first,
                                       const Atom &second) const;
  bool bondsConstitutionallyEquivalent(const Bond &first,
                                       const Bond &second) const;
  bool branchReachesCyclicCore(Atom *root, Atom *end) const;
  bool hasConstitutionalAutomorphism(
      unsigned int rootIdx, unsigned int fromIdx, unsigned int toIdx,
      std::span<const std::uint64_t> fixedAtoms,
      std::span<const unsigned int> addedFixedAtoms) const;
  ConstitutionalAutomorphismEvidence &getAutomorphismEvidence(
      const ConstitutionalAutomorphismKey &key) const;
  static std::optional<bool> findAutomorphismEvidence(
      const ConstitutionalAutomorphismEvidence &evidence,
      std::span<const std::uint64_t> fixedAtoms,
      std::span<const unsigned int> addedFixedAtoms);
  void addAutomorphismWitness(ConstitutionalAutomorphismEvidence &evidence,
                              std::vector<unsigned int> movedAtoms) const;
  void addAutomorphismFailure(ConstitutionalAutomorphismEvidence &evidence,
                              std::vector<unsigned int> fixedAtoms) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
