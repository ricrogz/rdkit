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

#include <memory>
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

  // Integer bond order of a kekulized molecule
  // Dative bonds get bond order 0.
  int getBondOrder(Bond *bond) const;

  double getAtomicMass(Atom *atom) const;

 private:
  ROMol &d_mol;
  std::vector<RDKit::Bond::BondType> d_kekulized_bonds;
  std::vector<FractionalAtomicNum> d_atomnums;
  std::vector<RDKit::Bond *> d_bonds;
  mutable std::vector<double> d_atomic_masses;
  mutable std::vector<unsigned char> d_atomic_mass_cached;
  boost::dynamic_bitset<> d_configuration_foci;
  // Two directed sides per bond: 0=unknown, 1=contains a focus, 2=no focus.
  mutable std::vector<unsigned char> d_configuration_branch_cache;
};

}  // namespace CIPLabeler
}  // namespace RDKit
