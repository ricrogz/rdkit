//
//
//  Copyright (C) 2018 Ricardo Rodriguez-Schmidt
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Sgroup.h"
#include "ROMol.h"

namespace RDKit {

void SGroup::setOwningMol(ROMol *mol) {
  PRECONDITION(mol, "owning molecule is nullptr");

  if (dp_mol != nullptr) {
    remap_atoms_bonds_to_new_mol(mol);
  }

  dp_mol = mol;
}

const std::string &SGroup::getStrProp(const std::string &prop) const {
  if (!hasStrProp(prop)) {
    std::ostringstream errout;
    errout << "String Property '" << prop << "' not found in SGroup "
           << getId();
    throw SGroupException(errout.str());
  }
  return d_strProp.at(prop);
}

//! check if SGroup has the given property set
bool SGroup::hasStrProp(const std::string &prop) const {
  return d_strProp.find(prop) != d_strProp.end();
}

void SGroup::setId(unsigned int id) {
  PRECONDITION(dp_mol, "SGroup is not owned by any molecule");

  if (id == 0) {
    d_id = dp_mol->getNextFreeSGroupId();
  } else {
    if (!dp_mol->isSGroupIdFree(id)) {
      std::ostringstream errout;
      errout << "ID " << id
             << " is already assigned to a SGroup on the same molecule";
      throw SGroupException(errout.str());
    }
    d_id = id;
  }
}

unsigned int SGroup::getIndexInMol() const {
  if (!dp_mol) {
    std::ostringstream errout;
    errout << "SGroup has no owning molecule" << std::endl;
    throw SGroupException(errout.str());
  }

  auto match_sgroup = [&](const SGROUP_SPTR &sg) { return this == sg.get(); };

  auto sgroupItr =
      std::find_if(dp_mol->beginSGroups(), dp_mol->endSGroups(), match_sgroup);

  if (sgroupItr == dp_mol->endSGroups()) {
    std::ostringstream errout;
    errout << "Unable to find own index in owning mol SGroup collection"
           << std::endl;
    throw SGroupException(errout.str());
  }

  return sgroupItr - dp_mol->beginSGroups();
}

void SGroup::addAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");

  Atom *atom = dp_mol->getAtomWithIdx(idx);
  if (!atom) {
    std::ostringstream errout;
    errout << "Cannot find Atom " << idx << " in same molecule as SGroup "
           << getId();
    throw SGroupException(errout.str());
  }
  d_atoms.push_back(atom);
}

void SGroup::addPAtomWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");
  PRECONDITION(!d_atoms.empty(), "no atoms");

  Atom *atom = dp_mol->getAtomWithIdx(idx);
  if (!atom) {
    std::ostringstream errout;
    errout << "Cannot find Atom " << idx << " in same molecule as SGroup "
           << getId();
    throw SGroupException(errout.str());
  } else if (std::find(d_atoms.begin(), d_atoms.end(), atom) == d_atoms.end()) {
    std::ostringstream errout;
    errout << "Atom " << idx << " is not a member of SGroup " << getId();
    throw SGroupException(errout.str());
  }

  d_patoms.push_back(atom);
}

void SGroup::addBondWithIdx(unsigned int idx) {
  PRECONDITION(dp_mol, "bad mol");

  Bond *bond = dp_mol->getBondWithIdx(idx);
  if (!bond) {
    std::ostringstream errout;
    errout << "Cannot find Bond " << idx << " in same molecule as SGroup "
           << getId();
    throw SGroupException(errout.str());
  }
  d_bonds.push_back(bond);
}

void SGroup::addBracket(Bracket &b) {
  PRECONDITION(d_brackets.size() < 2, "Too many brackets");
  d_brackets.push_back(b);
}

void SGroup::addCState(unsigned int xbond, RDGeom::Point3D *vectorPtr) {
  PRECONDITION(dp_mol, "bad mol");
  PRECONDITION(!d_bonds.empty(), "no bonds");

  Bond *bond = dp_mol->getBondWithIdx(xbond);
  if (!bond) {
    std::ostringstream errout;
    errout << "Cannot find Bond " << xbond << " in same molecule as SGroup "
           << getId();
    throw SGroupException(errout.str());
  }

  d_cstates.emplace_back(bond, vectorPtr);
}

void SGroup::addDataField(const std::string &data) {
  d_dataFields.push_back(data);
}

void SGroup::remap_atoms_bonds_to_new_mol(ROMol *other_mol) {
  for (auto &&atom : d_atoms) {
    unsigned int idx = atom->getIdx();
    atom = other_mol->getAtomWithIdx(idx);
  }

  for (auto &&patom : d_patoms) {
    unsigned int idx = patom->getIdx();
    patom = other_mol->getAtomWithIdx(idx);
  }

  for (auto &&bond : d_bonds) {
    unsigned int idx = bond->getIdx();
    bond = other_mol->getBondWithIdx(idx);
  }
}

void SGroup::remap_parent_sgroup_to_new_mol(ROMol *other_mol) {
  if (d_parent) {
    auto parent_idx = d_parent->getIndexInMol();
    d_parent = other_mol->getSGroup(parent_idx);
  }
}

void SGroup::addAttachPoint(const AttachPoint &sap) { d_saps.push_back(sap); }

//! check if the bond is SGroup XBOND or CBOND
SGroup::BondType SGroup::getBondType(Bond *bond) const {
  Atom *atom1 = bond->getBeginAtom();
  Atom *atom2 = bond->getEndAtom();

  bool atom1_in_sgroup =
      std::find(d_atoms.begin(), d_atoms.end(), atom1) != d_atoms.end();
  bool atom2_in_sgroup =
      std::find(d_atoms.begin(), d_atoms.end(), atom2) != d_atoms.end();

  if (atom1_in_sgroup && atom2_in_sgroup) {
    return SGroup::BondType::CBOND;
  } else if (atom1_in_sgroup || atom2_in_sgroup) {
    return SGroup::BondType::XBOND;
  } else {
    std::ostringstream errout;
    errout << "Unexpected SGroup Bond Type condition for bond "
           << bond->getIdx();
    throw SGroupException(errout.str());
  }
}

bool SGroupTypeOK(std::string typ) {
  const char *cfailTyps[15] = {
      // polymer sgroups:
      "SRU", "MON", "COP", "CRO", "GRA", "MOD", "MER", "ANY",
      // formulations/mixtures:
      "COM", "MIX", "FOR",
      // other
      "SUP", "MUL", "DAT", "GEN"};
  std::vector<std::string> failTyps(cfailTyps, cfailTyps + 15);
  return std::find(failTyps.begin(), failTyps.end(), typ) != failTyps.end();
}

bool SGroupSubTypeOK(std::string typ) {
  const char *cfailTyps[3] = {"ALT", "RAN", "BLO"};
  std::vector<std::string> failTyps(cfailTyps, cfailTyps + 3);
  return std::find(failTyps.begin(), failTyps.end(), typ) != failTyps.end();
}

bool SGroupConnectTypeOK(std::string typ) {
  const char *cfailTyps[3] = {"HH", "HT", "EU"};
  std::vector<std::string> failTyps(cfailTyps, cfailTyps + 3);
  return std::find(failTyps.begin(), failTyps.end(), typ) != failTyps.end();
}

}  // namespace RDKit

// TO DO: finish this!
std::ostream &operator<<(std::ostream &target, const RDKit::SGroup &sgroup) {
  target << sgroup.getId() << " " << sgroup.getType();
  return target;
}
