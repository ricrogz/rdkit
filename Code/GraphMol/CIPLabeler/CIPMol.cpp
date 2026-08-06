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
#include <cstddef>
#include <cstdint>
#include <utility>

#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

CIPMol::CIPMol(ROMol &mol) : d_mol{mol} {
  d_bonds.reserve(mol.getNumBonds());
  std::ranges::copy(mol.bonds(), std::back_inserter(d_bonds));
  d_atomic_masses.resize(mol.getNumAtoms());
  d_atomic_mass_cached.resize(mol.getNumAtoms());
}

const FractionalAtomicNum &CIPMol::getFractionalAtomicNum(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  if (d_atomnums.empty()) {
    const_cast<CIPMol *>(this)->d_atomnums = calcFracAtomNums(*this);
  }
  return d_atomnums[atom->getIdx()];
}

unsigned CIPMol::getNumAtoms() const { return d_mol.getNumAtoms(); }

unsigned CIPMol::getNumBonds() const { return d_mol.getNumBonds(); };

Atom *CIPMol::getAtom(int idx) const { return d_mol.getAtomWithIdx(idx); };

CXXAtomIterator<MolGraph, Atom *> CIPMol::atoms() const {
  return d_mol.atoms();
}

Bond *CIPMol::getBond(int idx) const { return d_bonds[idx]; };

CIPMolSpan<Bond *, ROMol::OEDGE_ITER> CIPMol::getBonds(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  return {d_mol, d_mol.getAtomBonds(atom)};
}

CIPMolSpan<Atom *, ROMol::ADJ_ITER> CIPMol::getNeighbors(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  return {d_mol, d_mol.getAtomNeighbors(atom)};
}

bool CIPMol::isInRing(Bond *bond) const {
  PRECONDITION(bond, "bad bond")
  const auto rings = d_mol.getRingInfo();

  if (!rings->isFindFastOrBetter()) {
    MolOps::fastFindRings(d_mol);
  }

  return rings->numBondRings(bond->getIdx()) != 0u;
};

void CIPMol::setConfigurationFoci(boost::dynamic_bitset<> foci) {
  PRECONDITION(foci.size() == getNumAtoms(),
               "configuration focus bitset has the wrong size")
  d_configuration_foci = std::move(foci);
  d_configuration_branch_cache.assign(2u * getNumBonds(), 0u);
}

bool CIPMol::isAcyclicBranchWithoutConfiguration(Bond *bond,
                                                 Atom *endAtom) const {
  if (bond == nullptr || endAtom == nullptr ||
      d_configuration_foci.size() != getNumAtoms() || isInRing(bond)) {
    return false;
  }

  const auto beginAtom = bond->getBeginAtom();
  const auto otherAtom = bond->getEndAtom();
  std::size_t side = 0;
  if (endAtom == beginAtom) {
    side = 0;
  } else if (endAtom == otherAtom) {
    side = 1;
  } else {
    return false;
  }

  const auto cacheIdx = 2u * bond->getIdx() + side;
  const auto cached = d_configuration_branch_cache[cacheIdx];
  if (cached != 0u) {
    return cached == 2u;
  }

  boost::dynamic_bitset<> seen(getNumAtoms());
  std::vector<Atom *> queue{endAtom};
  seen.set(endAtom->getIdx());
  bool containsFocus = false;

  for (std::size_t pos = 0; pos < queue.size() && !containsFocus; ++pos) {
    auto atom = queue[pos];
    if (d_configuration_foci.test(atom->getIdx())) {
      containsFocus = true;
      break;
    }
    for (auto candidate : getBonds(atom)) {
      if (candidate == bond) {
        continue;
      }
      auto neighbor = candidate->getOtherAtom(atom);
      if (!seen.test(neighbor->getIdx())) {
        seen.set(neighbor->getIdx());
        queue.push_back(neighbor);
      }
    }
  }

  d_configuration_branch_cache[cacheIdx] = containsFocus ? 1u : 2u;
  return !containsFocus;
}

bool CIPMol::hasConstitutionalAutomorphism(Atom *root, Atom *from,
                                           Atom *to) const {
  if (root == nullptr || from == nullptr || to == nullptr) {
    return false;
  }
  if (from == to) {
    return true;
  }

  // Exact self-isomorphism is intended as a shortcut for compact, highly
  // cyclic stereochemical units. Avoid introducing potentially expensive VF2
  // work on large molecules; a false result only disables the optimization.
  constexpr std::uint64_t MAX_AUTOMORPHISM_ATOMS = 128;
  const auto numAtoms = static_cast<std::uint64_t>(getNumAtoms());
  if (numAtoms == 0u || numAtoms > MAX_AUTOMORPHISM_ATOMS) {
    return false;
  }
  if (std::ranges::any_of(d_bonds, [](const Bond *bond) {
        return bond->getBondType() == Bond::AROMATIC;
      })) {
    return false;
  }

  const auto rootIdx = static_cast<std::uint64_t>(root->getIdx());
  auto fromIdx = static_cast<std::uint64_t>(from->getIdx());
  auto toIdx = static_cast<std::uint64_t>(to->getIdx());
  if (fromIdx > toIdx) {
    std::swap(fromIdx, toIdx);
  }
  const auto cacheKey = (rootIdx * numAtoms + fromIdx) * numAtoms + toIdx;
  const auto cached = d_constitutional_automorphism_cache.find(cacheKey);
  if (cached != d_constitutional_automorphism_cache.end()) {
    return cached->second;
  }

  SubstructMatchParameters params;
  params.recursionPossible = false;
  params.uniquify = false;
  params.maxMatches = 1;
  params.numThreads = 1;
  params.useChirality = false;
  params.extraAtomCheckOverridesDefaultCheck = true;
  params.extraAtomCheck = [this, rootIdx, fromIdx, toIdx](const Atom &queryAtom,
                                                          const Atom &molAtom) {
    const auto queryIdx = static_cast<std::uint64_t>(queryAtom.getIdx());
    const auto molIdx = static_cast<std::uint64_t>(molAtom.getIdx());
    if ((queryIdx == rootIdx) != (molIdx == rootIdx) ||
        (queryIdx == fromIdx) != (molIdx == toIdx)) {
      return false;
    }

    return queryAtom.getAtomicNum() == molAtom.getAtomicNum() &&
           queryAtom.getIsotope() == molAtom.getIsotope() &&
           queryAtom.getFormalCharge() == molAtom.getFormalCharge() &&
           queryAtom.getNumRadicalElectrons() ==
               molAtom.getNumRadicalElectrons() &&
           queryAtom.getTotalNumHs() == molAtom.getTotalNumHs();
  };
  params.extraBondCheckOverridesDefaultCheck = true;
  params.extraBondCheck = [](const Bond &queryBond, const Bond &molBond) {
    return queryBond.getBondType() == molBond.getBondType() &&
           queryBond.getIsAromatic() == molBond.getIsAromatic();
  };

  const bool equivalent = !SubstructMatch(d_mol, d_mol, params).empty();
  d_constitutional_automorphism_cache.emplace(cacheKey, equivalent);
  return equivalent;
}

int CIPMol::getBondOrder(Bond *bond) const {
  PRECONDITION(bond, "bad bond")
  if (d_kekulized_bonds.empty()) {
    auto &bonds =
        const_cast<std::vector<RDKit::Bond::BondType> &>(d_kekulized_bonds);
    bonds.reserve(d_mol.getNumBonds());

    const bool hasAromaticBond =
        std::ranges::any_of(d_bonds, [](const Bond *candidate) {
          return candidate->getBondType() == Bond::AROMATIC;
        });
    if (hasAromaticBond) {
      RWMol tmp{d_mol};
      const ROMol *bond_source = &tmp;
      try {
        MolOps::Kekulize(tmp);
      } catch (const MolSanitizeException &) {
        // Kekulize() may have changed some bonds before discovering that no
        // valid assignment exists. Fall back to the untouched input instead
        // of caching that partial assignment.
        bond_source = &d_mol;
      }
      for (const auto candidate : bond_source->bonds()) {
        bonds.push_back(candidate->getBondType());
      }
    } else {
      for (const auto candidate : d_bonds) {
        bonds.push_back(candidate->getBondType());
      }
    }
  }

  const auto bond_type = d_kekulized_bonds.at(bond->getIdx());

  // Dative bonds might need to be considered with a different bond order
  // for the end atom at the end of the bond.
  switch (bond_type) {
    case Bond::ZERO:
    case Bond::HYDROGEN:
    case Bond::DATIVE:
    case Bond::DATIVEL:
    case Bond::DATIVER:
      return 0;
    case Bond::SINGLE:
      return 1;
    case Bond::AROMATIC:
      BOOST_LOG(rdWarningLog)
          << "non kekulizable aromatic bond being treated as bond order 1"
          << std::endl;
      return 1;
    case Bond::DOUBLE:
      return 2;
    case Bond::TRIPLE:
      return 3;
    case Bond::QUADRUPLE:
      return 4;
    case Bond::QUINTUPLE:
      return 5;
    case Bond::HEXTUPLE:
      return 6;
    default:
      throw std::runtime_error("Non integer-order bonds are not allowed.");
  }
};

double CIPMol::getAtomicMass(Atom *atom) const {
  PRECONDITION(atom, "bad atom")
  const auto index = atom->getIdx();
  if (!d_atomic_mass_cached.at(index)) {
    d_atomic_masses[index] = atom->getMass();
    d_atomic_mass_cached[index] = true;
  }
  return d_atomic_masses[index];
}

}  // namespace CIPLabeler
}  // namespace RDKit
