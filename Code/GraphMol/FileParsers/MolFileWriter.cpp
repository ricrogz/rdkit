//
//  Copyright (C) 2003-2017 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  23/12/2013:
//     V3000 mol block writer contributed by Jan Holst Jensen
//
#include "FileParsers.h"
#include "MolFileStereochem.h"
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitQueries.h>
#include <RDGeneral/Ranking.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Depictor/RDDepictor.h>

namespace RDKit {

//*************************************
//
// Every effort has been made to adhere to MDL's standard
// for mol files
//
//*************************************

namespace {
int getQueryBondTopology(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(bond->hasQuery(), "no query");
  int res = 0;
  Bond::QUERYBOND_QUERY *qry = bond->getQuery();
  // start by catching combined bond order + bond topology queries

  if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
      qry->endChildren() - qry->beginChildren() == 2) {
    auto child1 = qry->beginChildren();
    auto child2 = child1 + 1;
    if ((*child1)->getDescription() == "BondOr" &&
        (*child2)->getDescription() == "BondInRing") {
      qry = child2->get();
    } else if ((*child1)->getDescription() == "BondInRing" &&
               (*child2)->getDescription() == "BondOr") {
      qry = child1->get();
    }
  }
  if (qry->getDescription() == "BondInRing") {
    if (qry->getNegation())
      res = 2;
    else
      res = 1;
  }
  return res;
}

// returns 0 if there's a basic bond-order query
int getQueryBondSymbol(const Bond *bond) {
  PRECONDITION(bond, "no bond");
  PRECONDITION(bond->hasQuery(), "no query");
  int res = 8;

  Bond::QUERYBOND_QUERY *qry = bond->getQuery();
  if (qry->getDescription() == "BondOrder") {
    // trap the simple bond-order query
    res = 0;
  } else {
    // start by catching combined bond order + bond topology queries
    if (qry->getDescription() == "BondAnd" && !qry->getNegation() &&
        qry->endChildren() - qry->beginChildren() == 2) {
      auto child1 = qry->beginChildren();
      auto child2 = child1 + 1;
      if ((*child1)->getDescription() == "BondOr" &&
          (*child2)->getDescription() == "BondInRing") {
        qry = child1->get();
      } else if ((*child1)->getDescription() == "BondInRing" &&
                 (*child2)->getDescription() == "BondOr") {
        qry = child2->get();
      }
    }
    if (qry->getDescription() == "BondOr" && !qry->getNegation()) {
      if (qry->endChildren() - qry->beginChildren() == 2) {
        auto child1 = qry->beginChildren();
        auto child2 = child1 + 1;
        if ((*child1)->getDescription() == "BondOrder" &&
            !(*child1)->getNegation() &&
            (*child2)->getDescription() == "BondOrder" &&
            !(*child2)->getNegation()) {
          // ok, it's a bond query we have a chance of dealing with
          int t1 = static_cast<BOND_EQUALS_QUERY *>(child1->get())->getVal();
          int t2 = static_cast<BOND_EQUALS_QUERY *>(child2->get())->getVal();
          if (t1 > t2) std::swap(t1, t2);
          if (t1 == Bond::SINGLE && t2 == Bond::DOUBLE) {
            res = 5;
          } else if (t1 == Bond::SINGLE && t2 == Bond::AROMATIC) {
            res = 6;
          } else if (t1 == Bond::DOUBLE && t2 == Bond::AROMATIC) {
            res = 7;
          }
        }
      }
    } else if (qry->getDescription() == "SingleOrAromaticBond" &&
               !qry->getNegation()) {
      res = 6;
    }
  }
  return res;
}
}  // namespace

const std::string GetMolFileChargeInfo(const RWMol &mol) {
  std::stringstream res;
  std::stringstream chgss;
  std::stringstream radss;
  std::stringstream massdiffss;
  unsigned int nChgs = 0;
  unsigned int nRads = 0;
  unsigned int nMassDiffs = 0;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    const Atom *atom = *atomIt;
    if (atom->getFormalCharge() != 0) {
      ++nChgs;
      chgss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                   atom->getFormalCharge();
      if (nChgs == 8) {
        res << boost::format("M  CHG%3d") % nChgs << chgss.str() << std::endl;
        chgss.str("");
        nChgs = 0;
      }
    }
    unsigned int nRadEs = atom->getNumRadicalElectrons();
    if (nRadEs != 0 && atom->getTotalDegree() != 0) {
      ++nRads;
      if (nRadEs % 2) {
        nRadEs = 2;
      } else {
        nRadEs = 3;  // we use triplets, not singlets:
      }
      radss << boost::format(" %3d %3d") % (atom->getIdx() + 1) % nRadEs;
      if (nRads == 8) {
        res << boost::format("M  RAD%3d") % nRads << radss.str() << std::endl;
        radss.str("");
        nRads = 0;
      }
    }
    if (!atom->hasQuery()) {
      int isotope = atom->getIsotope();
      if (isotope != 0) {
        ++nMassDiffs;
        massdiffss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                          isotope;
        if (nMassDiffs == 8) {
          res << boost::format("M  ISO%3d") % nMassDiffs << massdiffss.str()
              << std::endl;
          massdiffss.str("");
          nMassDiffs = 0;
        }
      }
    }
  }
  if (nChgs) {
    res << boost::format("M  CHG%3d") % nChgs << chgss.str() << std::endl;
  }
  if (nRads) {
    res << boost::format("M  RAD%3d") % nRads << radss.str() << std::endl;
  }
  if (nMassDiffs) {
    res << boost::format("M  ISO%3d") % nMassDiffs << massdiffss.str()
        << std::endl;
  }
  return res.str();
}

bool hasComplexQuery(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  bool res = false;
  if (atom->hasQuery()) {
    res = true;
    // counter examples:
    //  1) atomic number
    //  2) the smarts parser inserts AtomAnd queries
    //     for "C" or "c":
    //
    std::string descr = atom->getQuery()->getDescription();
    if (descr == "AtomAtomicNum") {
      res = false;
    } else if (descr == "AtomAnd") {
      if ((*atom->getQuery()->beginChildren())->getDescription() ==
          "AtomAtomicNum") {
        res = false;
      }
    }
  }
  return res;
}

bool isListQuery(const Atom::QUERYATOM_QUERY *q) {
  // list queries are series of nested ors of AtomAtomicNum queries
  PRECONDITION(q, "bad query");
  bool res = false;
  std::string descr = q->getDescription();
  if (descr == "AtomOr") {
    res = true;
    for (auto cIt = q->beginChildren(); cIt != q->endChildren() && res; ++cIt) {
      std::string descr = (*cIt)->getDescription();
      // we don't allow negation of any children of the query:
      if ((*cIt)->getNegation()) {
        res = false;
      } else if (descr == "AtomOr") {
        res = isListQuery((*cIt).get());
      } else if (descr != "AtomAtomicNum") {
        res = false;
      }
    }
  }
  return res;
}

void getListQueryVals(const Atom::QUERYATOM_QUERY *q, INT_VECT &vals) {
  // list queries are series of nested ors of AtomAtomicNum queries
  PRECONDITION(q, "bad query");
  std::string descr = q->getDescription();
  PRECONDITION(descr == "AtomOr", "bad query");
  if (descr == "AtomOr") {
    for (auto cIt = q->beginChildren(); cIt != q->endChildren(); ++cIt) {
      std::string descr = (*cIt)->getDescription();
      CHECK_INVARIANT((descr == "AtomOr" || descr == "AtomAtomicNum"),
                      "bad query");
      // we don't allow negation of any children of the query:
      if (descr == "AtomOr") {
        getListQueryVals((*cIt).get(), vals);
      } else if (descr == "AtomAtomicNum") {
        vals.push_back(
            static_cast<ATOM_EQUALS_QUERY *>((*cIt).get())->getVal());
      }
    }
  }
}

bool hasListQuery(const Atom *atom) {
  PRECONDITION(atom, "bad atom");
  bool res = false;
  if (atom->hasQuery()) {
    res = isListQuery(atom->getQuery());
  }
  return res;
}

const std::string GetMolFileQueryInfo(const RWMol &mol) {
  std::stringstream ss;
  boost::dynamic_bitset<> listQs(mol.getNumAtoms());
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    if (hasListQuery(*atomIt)) listQs.set((*atomIt)->getIdx());
  }
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    bool wrote_query = false;
    if (!listQs[(*atomIt)->getIdx()] && hasComplexQuery(*atomIt)) {
      std::string sma =
          SmartsWrite::GetAtomSmarts(static_cast<const QueryAtom *>(*atomIt));
      ss << "V  " << std::setw(3) << (*atomIt)->getIdx() + 1 << " " << sma
         << std::endl;
      wrote_query = true;
    }
    std::string molFileValue;
    if (!wrote_query && (*atomIt)->getPropIfPresent(
                            common_properties::molFileValue, molFileValue))
      ss << "V  " << std::setw(3) << (*atomIt)->getIdx() + 1 << " "
         << molFileValue << std::endl;
  }
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    if (listQs[(*atomIt)->getIdx()]) {
      INT_VECT vals;
      getListQueryVals((*atomIt)->getQuery(), vals);
      ss << "M  ALS " << std::setw(3) << (*atomIt)->getIdx() + 1 << " ";
      ss << std::setw(2) << vals.size();
      if ((*atomIt)->getQuery()->getNegation()) {
        ss << " T ";
      } else {
        ss << " F ";
      }
      BOOST_FOREACH (int val, vals) {
        ss << std::setw(4) << std::left
           << (PeriodicTable::getTable()->getElementSymbol(val));
      }
      ss << "\n";
    }
  }
  return ss.str();
}

const std::string GetMolFileRGroupInfo(const RWMol &mol) {
  std::stringstream ss;
  unsigned int nEntries = 0;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    unsigned int lbl;
    if ((*atomIt)->getPropIfPresent(common_properties::_MolFileRLabel, lbl)) {
      ss << " " << std::setw(3) << (*atomIt)->getIdx() + 1 << " "
         << std::setw(3) << lbl;
      ++nEntries;
    }
  }
  std::stringstream ss2;
  if (nEntries)
    ss2 << "M  RGP" << std::setw(3) << nEntries << ss.str() << std::endl;
  return ss2.str();
}

const std::string GetMolFileAliasInfo(const RWMol &mol) {
  std::stringstream ss;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    std::string lbl;
    if ((*atomIt)->getPropIfPresent(common_properties::molFileAlias, lbl)) {
      if (!lbl.empty())
        ss << "A  " << std::setw(3) << (*atomIt)->getIdx() + 1 << "\n"
           << lbl << "\n";
    }
  }
  return ss.str();
}

const std::string GetMolFileZBOInfo(const RWMol &mol) {
  std::stringstream res;
  std::stringstream ss;
  unsigned int nEntries = 0;
  boost::dynamic_bitset<> atomsAffected(mol.getNumAtoms(), 0);
  for (ROMol::ConstBondIterator bondIt = mol.beginBonds();
       bondIt != mol.endBonds(); ++bondIt) {
    if ((*bondIt)->getBondType() == Bond::ZERO) {
      ++nEntries;
      ss << " " << std::setw(3) << (*bondIt)->getIdx() + 1 << " "
         << std::setw(3) << 0;
      if (nEntries == 8) {
        res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << std::endl;
        nEntries = 0;
        ss.str("");
      }
      atomsAffected[(*bondIt)->getBeginAtomIdx()] = 1;
      atomsAffected[(*bondIt)->getEndAtomIdx()] = 1;
    }
  }
  if (nEntries) {
    res << "M  ZBO" << std::setw(3) << nEntries << ss.str() << std::endl;
  }
  if (atomsAffected.count()) {
    std::stringstream hydss;
    unsigned int nhyd = 0;
    std::stringstream zchss;
    unsigned int nzch = 0;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
      if (!atomsAffected[i]) continue;
      const Atom *atom = mol.getAtomWithIdx(i);
      nhyd++;
      hydss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                   atom->getTotalNumHs();
      if (nhyd == 8) {
        res << boost::format("M  HYD%3d") % nhyd << hydss.str() << std::endl;
        hydss.str("");
        nhyd = 0;
      }
      if (atom->getFormalCharge()) {
        nzch++;
        zchss << boost::format(" %3d %3d") % (atom->getIdx() + 1) %
                     atom->getFormalCharge();
        if (nzch == 8) {
          res << boost::format("M  ZCH%3d") % nzch << zchss.str() << std::endl;
          zchss.str("");
          nzch = 0;
        }
      }
    }
    if (nhyd) {
      res << boost::format("M  HYD%3d") % nhyd << hydss.str() << std::endl;
    }
    if (nzch) {
      res << boost::format("M  ZCH%3d") % nzch << zchss.str() << std::endl;
    }
  }
  return res.str();
}

const std::string AtomGetMolFileSymbol(const Atom *atom, bool padWithSpaces) {
  PRECONDITION(atom, "");

  std::string res;
  if (atom->hasProp(common_properties::_MolFileRLabel)) {
    res = "R#";
    //    } else if(!atom->hasQuery() && atom->getAtomicNum()){
  } else if (atom->getAtomicNum()) {
    res = atom->getSymbol();
  } else {
    if (!atom->hasProp(common_properties::dummyLabel)) {
      if (atom->hasQuery() && atom->getQuery()->getNegation() &&
          atom->getQuery()->getDescription() == "AtomAtomicNum" &&
          static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery())->getVal() == 1) {
        res = "A";
      } else if (atom->hasQuery() && atom->getQuery()->getNegation() &&
                 atom->getQuery()->getDescription() == "AtomOr" &&
                 atom->getQuery()->endChildren() -
                         atom->getQuery()->beginChildren() ==
                     2 &&
                 (*atom->getQuery()->beginChildren())->getDescription() ==
                     "AtomAtomicNum" &&
                 static_cast<ATOM_EQUALS_QUERY *>(
                     (*atom->getQuery()->beginChildren()).get())
                         ->getVal() == 6 &&
                 (*++(atom->getQuery()->beginChildren()))->getDescription() ==
                     "AtomAtomicNum" &&
                 static_cast<ATOM_EQUALS_QUERY *>(
                     (*++(atom->getQuery()->beginChildren())).get())
                         ->getVal() == 1) {
        res = "Q";
      } else if (hasComplexQuery(atom)) {
        if (hasListQuery(atom)) {
          res = "L";
        } else {
          res = "*";
        }
      } else {
        res = "R";
      }
    } else {
      std::string symb;
      atom->getProp(common_properties::dummyLabel, symb);
      if (symb == "*")
        res = "R";
      else if (symb == "X")
        res = "R";
      else if (symb == "Xa")
        res = "R1";
      else if (symb == "Xb")
        res = "R2";
      else if (symb == "Xc")
        res = "R3";
      else if (symb == "Xd")
        res = "R4";
      else if (symb == "Xf")
        res = "R5";
      else if (symb == "Xg")
        res = "R6";
      else if (symb == "Xh")
        res = "R7";
      else if (symb == "Xi")
        res = "R8";
      else if (symb == "Xj")
        res = "R9";
      else
        res = symb;
    }
  }
  // pad the end with spaces
  if (padWithSpaces) {
    while (res.size() < 3) res += " ";
  }
  return res;
}

namespace {
unsigned int getAtomParityFlag(const Atom *atom, const Conformer *conf) {
  PRECONDITION(atom, "bad atom");
  PRECONDITION(conf, "bad conformer");
  if (!conf->is3D() || !(atom->getDegree() >= 3 && atom->getTotalDegree() == 4))
    return 0;

  const ROMol &mol = atom->getOwningMol();
  RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
  std::vector<std::pair<unsigned int, RDGeom::Point3D>> vs;
  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
  while (nbrIdx != endNbrs) {
    const Atom *at = mol.getAtomWithIdx(*nbrIdx);
    unsigned int idx = at->getIdx();
    RDGeom::Point3D v = conf->getAtomPos(idx);
    v -= pos;
    if (at->getAtomicNum() == 1) {
      idx += mol.getNumAtoms();
    }
    vs.push_back(std::make_pair(idx, v));
    ++nbrIdx;
  }
  std::sort(vs.begin(), vs.end(),
            Rankers::pairLess<unsigned int, RDGeom::Point3D>());
  double vol;
  if (vs.size() == 4) {
    vol = vs[0].second.crossProduct(vs[1].second).dotProduct(vs[3].second);
  } else {
    vol = -vs[0].second.crossProduct(vs[1].second).dotProduct(vs[2].second);
  }
  if (vol < 0) {
    return 2;
  } else if (vol > 0) {
    return 1;
  }
  return 0;
}
}  // namespace

bool hasNonDefaultValence(const Atom *atom) {
  if (atom->getNumRadicalElectrons() != 0) return true;
  if (atom->hasQuery()) return false;
  switch (atom->getAtomicNum()) {
    case 1:   // H
    case 5:   // B
    case 6:   // C
    case 7:   // N
    case 8:   // O
    case 9:   // F
    case 15:  // P
    case 16:  // S
    case 17:  // Cl
    case 35:  // Br
    case 53:  // I
      return false;
  }
  return true;
}

void GetMolFileAtomProperties(const Atom *atom, const Conformer *conf,
                              int &totValence, int &atomMapNumber,
                              unsigned int &parityFlag, double &x, double &y,
                              double &z) {
  PRECONDITION(atom, "");
  totValence = 0;
  atomMapNumber = 0;
  parityFlag = 0;
  x = y = z = 0.0;

  if (!atom->getPropIfPresent(common_properties::molAtomMapNumber,
                              atomMapNumber)) {
    // XXX FIX ME->should we fail here? previously we would not assign
    // the atomMapNumber if it didn't exist which could result in garbage
    //  values.
    atomMapNumber = 0;
  }

  if (conf) {
    const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    x = pos.x;
    y = pos.y;
    z = pos.z;
    if (conf->is3D() && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        atom->getChiralTag() != Atom::CHI_OTHER && atom->getDegree() >= 3 &&
        atom->getTotalDegree() == 4) {
      parityFlag = getAtomParityFlag(atom, conf);
    }
  }
  if (hasNonDefaultValence(atom)) {
    if (atom->getTotalDegree() == 0) {
      // Specify zero valence for elements/metals without neighbors
      // or hydrogens (degree 0) instead of writing them as radicals.
      totValence = 15;
    } else {
      // write the total valence for other atoms
      totValence = atom->getTotalValence() % 15;
    }
  }
}

const std::string GetMolFileAtomLines(const Atom *atom,
                                      const Conformer *conf = nullptr) {
  PRECONDITION(atom, "");
  std::string res;
  int totValence, atomMapNumber;
  unsigned int parityFlag;
  double x, y, z;
  GetMolFileAtomProperties(atom, conf, totValence, atomMapNumber, parityFlag, x,
                           y, z);

  int massDiff, chg, stereoCare, hCount, rxnComponentType, rxnComponentNumber,
      inversionFlag, exactChangeFlag;
  massDiff = 0;
  chg = 0;
  stereoCare = 0;
  hCount = 0;
  rxnComponentType = 0;
  rxnComponentNumber = 0;
  inversionFlag = 0;
  exactChangeFlag = 0;

  atom->getPropIfPresent(common_properties::molRxnRole, rxnComponentType);
  atom->getPropIfPresent(common_properties::molRxnComponent,
                         rxnComponentNumber);

  std::string symbol = AtomGetMolFileSymbol(atom, true);
#if 0
  const boost::format fmter(
      "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d");
  std::stringstream ss;
  ss << boost::format(fmter) % x % y % z % symbol.c_str() % massDiff % chg %
            parityFlag % hCount % stereoCare % totValence % rxnComponentType %
            rxnComponentNumber % atomMapNumber % inversionFlag %
            exactChangeFlag;
  res += ss.str();
#else
  // it feels ugly to use snprintf instead of boost::format, but at least of the
  // time of this writing (with boost 1.55), the snprintf version runs in 20% of
  // the time.
  char dest[128];
#ifndef _MSC_VER
  snprintf(dest, 128,
           "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
           z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
           totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
           inversionFlag, exactChangeFlag);
#else
  // ok, technically we should be being more careful about this, but tiven that
  // the format string makes it impossible for this to overflow, I think we're
  // safe. I just used the snprintf above to prevent linters from complaining
  // about use of sprintf
  sprintf_s(dest, 128,
            "%10.4f%10.4f%10.4f %3s%2d%3d%3d%3d%3d%3d  0%3d%3d%3d%3d%3d", x, y,
            z, symbol.c_str(), massDiff, chg, parityFlag, hCount, stereoCare,
            totValence, rxnComponentType, rxnComponentNumber, atomMapNumber,
            inversionFlag, exactChangeFlag);

#endif
  res += dest;
#endif
  return res;
};

namespace {
/*
  If a molecule contains dative bonds the V2000 format should not
  be used as it doesn't support dative bonds. If a dative bond is
  detected while writing a V2000 molfile the RequiresV3000Exception
  is thrown and the V2000 writer will redo the export in V3000 format.

  This is arguably a rather brute-force way of detecting the proper output
  format, but the only alternatives I (Jan Holst Jensen) had in mind were:

    1) Check all bond types before output. Slow and would affect all
       V2000 exports.
    2) Maintain a reference count of dative bonds in molecule. Complex
       and error-prone.
*/
class RequiresV3000Exception : public std::runtime_error {
 public:
  explicit RequiresV3000Exception()
      : std::runtime_error("RequiresV3000Exception"){};
};
}  // namespace

int BondGetMolFileSymbol(const Bond *bond) {
  PRECONDITION(bond, "");
  // FIX: should eventually recognize queries
  int res = 0;
  if (bond->hasQuery()) {
    res = getQueryBondSymbol(bond);
  }
  if (!res) {
    switch (bond->getBondType()) {
      case Bond::SINGLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 1;
        }
        break;
      case Bond::DOUBLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 2;
        }
        break;
      case Bond::TRIPLE:
        res = 3;
        break;
      case Bond::AROMATIC:
        res = 4;
        break;
      case Bond::ZERO:
        res = 1;
        break;
      case Bond::DATIVE:
        // Dative bonds requires V3000 format. Throw special exception to
        // force output to be re-done in V3000.
        throw RequiresV3000Exception();
      default:
        break;
    }
  }
  return res;
  // return res.c_str();
}

// only valid for single bonds
int BondGetDirCode(const Bond::BondDir dir) {
  int res = 0;
  switch (dir) {
    case Bond::NONE:
      res = 0;
      break;
    case Bond::BEGINWEDGE:
      res = 1;
      break;
    case Bond::BEGINDASH:
      res = 6;
      break;
    case Bond::UNKNOWN:
      res = 4;
      break;
    default:
      break;
  }
  return res;
}

void GetMolFileBondStereoInfo(const Bond *bond, const INT_MAP_INT &wedgeBonds,
                              const Conformer *conf, int &dirCode,
                              bool &reverse) {
  PRECONDITION(bond, "");
  dirCode = 0;
  reverse = false;
  Bond::BondDir dir = Bond::NONE;
  if (bond->getBondType() == Bond::SINGLE) {
    // single bond stereo chemistry
    dir = DetermineBondWedgeState(bond, wedgeBonds, conf);
    dirCode = BondGetDirCode(dir);
    // if this bond needs to be wedged it is possible that this
    // wedging was determined by a chiral atom at the end of the
    // bond (instead of at the beginning). In this case we need to
    // reverse the begin and end atoms for the bond when we write
    // the mol file
    if ((dirCode == 1) || (dirCode == 6)) {
      auto wbi = wedgeBonds.find(bond->getIdx());
      if (wbi != wedgeBonds.end() &&
          static_cast<unsigned int>(wbi->second) != bond->getBeginAtomIdx()) {
        reverse = true;
      }
    }
  } else if (bond->getBondType() == Bond::DOUBLE) {
    // double bond stereochemistry -
    // if the bond isn't specified, then it should go in the mol block
    // as "any", this was sf.net issue 2963522.
    // two caveats to this:
    // 1) if it's a ring bond, we'll only put the "any"
    //    in the mol block if the user specifically asked for it.
    //    Constantly seeing crossed bonds in rings, though maybe
    //    technically correct, is irritating.
    // 2) if it's a terminal bond (where there's no chance of
    //    stereochemistry anyway), we also skip the any.
    //    this was sf.net issue 3009756
    if (bond->getStereo() <= Bond::STEREOANY) {
      if (bond->getStereo() == Bond::STEREOANY) {
        dirCode = 3;
      } else if (!(bond->getOwningMol().getRingInfo()->numBondRings(
                     bond->getIdx())) &&
                 bond->getBeginAtom()->getDegree() > 1 &&
                 bond->getEndAtom()->getDegree() > 1) {
        // we don't know that it's explicitly unspecified (covered above with
        // the ==STEREOANY check)
        // look to see if one of the atoms has a bond with direction set
        if (bond->getBondDir() == Bond::EITHERDOUBLE) {
          dirCode = 3;
        } else {
          bool nbrHasDir = false;

          ROMol::OEDGE_ITER beg, end;
          boost::tie(beg, end) =
              bond->getOwningMol().getAtomBonds(bond->getBeginAtom());
          while (beg != end && !nbrHasDir) {
            const Bond *nbrBond = bond->getOwningMol()[*beg];
            if (nbrBond->getBondType() == Bond::SINGLE &&
                (nbrBond->getBondDir() == Bond::ENDUPRIGHT ||
                 nbrBond->getBondDir() == Bond::ENDDOWNRIGHT)) {
              nbrHasDir = true;
            }
            ++beg;
          }
          boost::tie(beg, end) =
              bond->getOwningMol().getAtomBonds(bond->getEndAtom());
          while (beg != end && !nbrHasDir) {
            const Bond *nbrBond = bond->getOwningMol()[*beg];
            if (nbrBond->getBondType() == Bond::SINGLE &&
                (nbrBond->getBondDir() == Bond::ENDUPRIGHT ||
                 nbrBond->getBondDir() == Bond::ENDDOWNRIGHT)) {
              nbrHasDir = true;
            }
            ++beg;
          }
          if (!nbrHasDir) {
            dirCode = 3;
          }
        }
      }
    }
  }
}

const std::string GetMolFileBondLines(const Bond *bond,
                                      const INT_MAP_INT &wedgeBonds,
                                      const Conformer *conf) {
  PRECONDITION(bond, "");

  int dirCode;
  bool reverse;
  GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode, reverse);
  int symbol = BondGetMolFileSymbol(bond);

  std::stringstream ss;
  if (reverse) {
    // switch the begin and end atoms on the bond line
    ss << std::setw(3) << bond->getEndAtomIdx() + 1;
    ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
  } else {
    ss << std::setw(3) << bond->getBeginAtomIdx() + 1;
    ss << std::setw(3) << bond->getEndAtomIdx() + 1;
  }
  ss << std::setw(3) << symbol;
  ss << " " << std::setw(2) << dirCode;

  if (bond->hasQuery()) {
    int topol = getQueryBondTopology(bond);
    if (topol) {
      ss << " " << std::setw(2) << 0 << " " << std::setw(2) << topol;
    }
  }

  return ss.str();
}

const std::string GetV3000MolFileAtomLines(const Atom *atom,
                                           const Conformer *conf = nullptr) {
  PRECONDITION(atom, "");
  int totValence, atomMapNumber;
  unsigned int parityFlag;
  double x, y, z;
  GetMolFileAtomProperties(atom, conf, totValence, atomMapNumber, parityFlag, x,
                           y, z);

  std::stringstream ss;
  ss << "M  V30 " << atom->getIdx() + 1;

  std::string symbol = AtomGetMolFileSymbol(atom, false);
  if (!hasListQuery(atom)) {
    ss << " " << symbol;
  } else {
    INT_VECT vals;
    getListQueryVals(atom->getQuery(), vals);
    if (atom->getQuery()->getNegation())
      ss << " "
         << "\"NOT";
    ss << " [";
    for (unsigned int i = 0; i < vals.size(); ++i) {
      if (i != 0) ss << ",";
      ss << PeriodicTable::getTable()->getElementSymbol(vals[i]);
    }
    ss << "]";
    if (atom->getQuery()->getNegation()) ss << "\"";
  }

  ss << " " << x << " " << y << " " << z;
  ss << " " << atomMapNumber;

  // Extra atom properties.
  int chg = atom->getFormalCharge();
  int isotope = atom->getIsotope();
  if (parityFlag != 0) {
    ss << " CFG=" << parityFlag;
  }
  if (chg != 0) {
    ss << " CHG=" << chg;
  }
  if (isotope != 0) {
    // the documentation for V3000 CTABs says that this should contain the
    // "absolute atomic weight" (whatever that means).
    // Online examples seem to have integer (isotope) values and Marvin won't
    // even read something that has a float.
    // We'll go with the int.
    int mass = static_cast<int>(round(atom->getMass()));
    // dummies may have an isotope set but they always have a mass of zero:
    if (!mass) mass = isotope;
    ss << " MASS=" << mass;
  }

  unsigned int nRadEs = atom->getNumRadicalElectrons();
  if (nRadEs != 0 && atom->getTotalDegree() != 0) {
    if (nRadEs % 2) {
      nRadEs = 2;
    } else {
      nRadEs = 3;  // we use triplets, not singlets:
    }
    ss << " RAD=" << nRadEs;
  }

  if (totValence != 0) {
    if (totValence == 15) {
      ss << " VAL=-1";
    } else {
      ss << " VAL=" << totValence;
    }
  }
  if (symbol == "R#") {
    unsigned int rLabel = 1;
    atom->getPropIfPresent(common_properties::_MolFileRLabel, rLabel);
    ss << " RGROUPS=(1 " << rLabel << ")";
  }
  // HCOUNT - *query* hydrogen count. Not written by this writer.

  return ss.str();
};

int GetV3000BondCode(const Bond *bond) {
  // JHJ: As far as I can tell, the V3000 bond codes are the same as for V2000.
  //      Except: The dative bond type is only supported in V3000.
  PRECONDITION(bond, "");
  int res = 0;
  // FIX: should eventually recognize queries
  if (bond->hasQuery()) res = getQueryBondSymbol(bond);
  if (!res) {
    switch (bond->getBondType()) {
      case Bond::SINGLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 1;
        }
        break;
      case Bond::DOUBLE:
        if (bond->getIsAromatic()) {
          res = 4;
        } else {
          res = 2;
        }
        break;
      case Bond::TRIPLE:
        res = 3;
        break;
      case Bond::AROMATIC:
        res = 4;
        break;
      case Bond::DATIVE:
        res = 9;
        break;
      default:
        res = 0;
        break;
    }
  }
  return res;
}

int BondStereoCodeV2000ToV3000(int dirCode) {
  // The Any bond configuration (code 4 in v2000 ctabs) seems to be missing
  switch (dirCode) {
    case 0:
      return 0;
    case 1:
      return 1;  // V2000 Up       => Up.
    case 3:
      return 2;  // V2000 Unknown  => Either.
    case 4:
      return 2;  // V2000 Any      => Either.
    case 6:
      return 3;  // V2000 Down     => Down.
    default:
      return 0;
  }
}

const std::string GetV3000MolFileBondLines(const Bond *bond,
                                           const INT_MAP_INT &wedgeBonds,
                                           const Conformer *conf) {
  PRECONDITION(bond, "");

  int dirCode;
  bool reverse;
  GetMolFileBondStereoInfo(bond, wedgeBonds, conf, dirCode, reverse);

  std::stringstream ss;
  ss << "M  V30 " << bond->getIdx() + 1;
  ss << " " << GetV3000BondCode(bond);
  if (reverse) {
    // switch the begin and end atoms on the bond line
    ss << " " << bond->getEndAtomIdx() + 1;
    ss << " " << bond->getBeginAtomIdx() + 1;
  } else {
    ss << " " << bond->getBeginAtomIdx() + 1;
    ss << " " << bond->getEndAtomIdx() + 1;
  }
  if (dirCode != 0) {
    ss << " CFG=" << BondStereoCodeV2000ToV3000(dirCode);
  }
  if (bond->hasQuery()) {
    int topol = getQueryBondTopology(bond);
    if (topol) {
      ss << " TOPO=" << topol;
    }
  }
  return ss.str();
}

inline std::string FormatV2000IntField(int value) {
  char output[5];
  snprintf(output, 5, " %3d", value);
  return std::string(output);
}

inline std::string FormatV2000NumEntriesField(int value) {
  char output[4];
  snprintf(output, 4, " %2d", value);
  return std::string(output);
}

inline std::string FormatV2000DoubleField(double value) {
  char output[11];
  snprintf(output, 11, "%10.4f", value);
  return std::string(output);
}

inline std::string FormatV2000StringField(const std::string &value,
                                          unsigned int fieldSize, bool pad,
                                          bool addSeparator) {
  std::ostringstream os;
  if (addSeparator) {
    os << ' ';
  }
  if (value.size() >= fieldSize) {
    os << value.substr(0, fieldSize);
  } else if (pad) {
    os << std::setw(fieldSize) << std::left << value;
  } else {
    os << value;
  }
  return os.str();
}

std::string BuildV2000STYLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = 1 + (sGroupItr - mol.beginSGroups());
    temp << FormatV2000IntField(idx)
         << FormatV2000StringField((*sGroupItr)->getType(), 3, true, true);
    if (++count == 8) {
      ret << "M  STY" << FormatV2000NumEntriesField(8) << temp.str()
          << std::endl;
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  STY" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000StringPropLines(const unsigned int entriesPerLine,
                                      const ROMol &mol,
                                      const std::string &propName,
                                      const std::string &propCode,
                                      const unsigned int fieldWitdh) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    if ((*sGroupItr)->hasStrProp(propName)) {  // Write field only if defined
      temp << FormatV2000IntField(idx)
           << FormatV2000StringField((*sGroupItr)->getStrProp(propName),
                                     fieldWitdh, true, true);
      if (++count == entriesPerLine) {
        ret << "M  " << propCode << FormatV2000NumEntriesField(entriesPerLine)
            << temp.str() << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  " << propCode << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SLBLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    unsigned int id = (*sGroupItr)->getId();
    if (id > 0) {  // Write field only if specific id was assigned
      temp << FormatV2000IntField(idx) << FormatV2000IntField(id);
      if (++count == 8) {
        ret << "M  SLB" << FormatV2000NumEntriesField(8) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SLB" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SDSLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    // Write field only if defined
    if ((*sGroupItr)->hasStrProp("ESTATE") &&
        (*sGroupItr)->getStrProp("ESTATE") == "E") {
      temp << FormatV2000IntField(idx);
      if (++count == 15) {
        ret << "M  SDS EXP" << FormatV2000NumEntriesField(15) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SDS EXP" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SPLLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    SGROUP_SPTR *parent = (*sGroupItr)->getParent();
    if (parent) {  // Write field only if a parent is defined
      unsigned int parentIdx = 1 + (*parent)->getIndexInMol();
      temp << FormatV2000IntField(idx) << FormatV2000IntField(parentIdx);
      if (++count == 8) {
        ret << "M  SPL" << FormatV2000NumEntriesField(8) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SPL" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SBTLines(const ROMol &mol) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    if ((*sGroupItr)->hasStrProp("BRKTYP")) {
      std::string bracketType = (*sGroupItr)->getStrProp("BRKTYP");
      if (bracketType == "BRACKET") {
        temp << FormatV2000IntField(idx) << FormatV2000IntField(0);
      } else if (bracketType == "PAREN") {
        temp << FormatV2000IntField(idx) << FormatV2000IntField(1);
      } else {
        std::ostringstream errout;
        errout << "Invalid BRKTYP value '" << bracketType << "' for SGroup "
               << idx;
        throw SGroupException(errout.str());
      }
      if (++count == 8) {
        ret << "M  SBT" << FormatV2000NumEntriesField(8) << temp.str()
            << std::endl;
        temp.str("");
        count = 0;
      }
    }
  }
  if (count) {
    ret << "M  SBT" << FormatV2000NumEntriesField(count) << temp.str()
        << std::endl;
  }
  return ret.str();
}

template <class T>
std::string BuildV2000IdxVectorDataLines(const unsigned int entriesPerLine,
                                         const unsigned int sGroupId,
                                         const std::string &code,
                                         const T &dataVector) {
  std::ostringstream ret;
  std::ostringstream temp;

  unsigned int count = 0;
  for (const auto &element : dataVector) {
    temp << FormatV2000IntField(1 + element->getIdx());
    if (++count == entriesPerLine) {
      ret << "M  " << code << FormatV2000IntField(sGroupId)
          << FormatV2000NumEntriesField(entriesPerLine) << temp.str()
          << std::endl;
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  " << code << FormatV2000IntField(sGroupId)
        << FormatV2000NumEntriesField(count) << temp.str() << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SMTLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->getType() == "MUL" && sgroup->hasStrProp("MULT")) {
    ret << "M  SMT" << FormatV2000IntField(idx)
        << FormatV2000StringField(sgroup->getStrProp("MULT"), 69, false, true)
        << std::endl;
  } else if (sgroup->hasStrProp("LABEL")) {
    ret << "M  SMT" << FormatV2000IntField(idx)
        << FormatV2000StringField(sgroup->getStrProp("LABEL"), 69, false, true)
        << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SDILine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  const std::vector<SGroup::Bracket> brackets = sgroup->getBrackets();

  if (brackets.size() == 2) {
    for (const auto &bracket : brackets) {
      ret << "M  SDI" << FormatV2000IntField(idx)
          << FormatV2000NumEntriesField(4);

      for (unsigned int iPoint = 0; iPoint < 2; ++iPoint) {
        ret << FormatV2000DoubleField(bracket.at(iPoint).x);
        ret << FormatV2000DoubleField(bracket.at(iPoint).y);
      }
      ret << std::endl;
    }
  } else if (!brackets.empty()) {
    std::ostringstream errout;
    errout << "SGroup " << sgroup->getId()
           << " has an unsupported number of brackets (" << brackets.size()
           << ")";
    throw SGroupException(errout.str());
  }

  return ret.str();
}

std::string BuildV2000SBVLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  for (const auto &cstate : sgroup->getCStates()) {
    ret << "M  SBV" << FormatV2000IntField(idx) << cstate.bond->getIdx();
    if (cstate.vector) {
      ret << FormatV2000DoubleField(cstate.vector->x);
      ret << FormatV2000DoubleField(cstate.vector->y);
    } else {
      ret << FormatV2000DoubleField(0);
      ret << FormatV2000DoubleField(0);
    }
    ret << std::endl;
  }

  return ret.str();
}

std::string BuildV2000SDTLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasStrProp("FIELDNAME")) {
    ret << "M  SDT" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getStrProp("FIELDNAME"), 30, true,
                                  true);

    if (sgroup->hasStrProp("FIELDTYPE")) {
      ret << FormatV2000StringField(sgroup->getStrProp("FIELDTYPE"), 2, true,
                                    false);
    } else {
      ret << " T";
    }

    if (sgroup->hasStrProp("FIELDINFO")) {
      ret << FormatV2000StringField(sgroup->getStrProp("FIELDINFO"), 20, true,
                                    false);
    }

    if (sgroup->hasStrProp("QUERYTYPE")) {
      ret << FormatV2000StringField(sgroup->getStrProp("QUERYTYPE"), 2, true,
                                    false);
    }
    if (sgroup->hasStrProp("QUERYOP")) {
      ret << FormatV2000StringField(sgroup->getStrProp("QUERYOP"), 15, true,
                                    false);
    }

    ret << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SDDLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasStrProp("FIELDDISP")) {
    ret << "M  SDD" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getStrProp("FIELDDISP"), 69, false,
                                  true);
    ret << std::endl;
  }

  return ret.str();
}

std::string BuildV2000SCDSEDLines(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  for (const auto &data : sgroup->getDataFields()) {
    unsigned int length = data.size();
    if (length > 200) {
      std::ostringstream errout;
      errout << "Data field '" << data << "' in SGroup " << sgroup->getId()
             << " is longed than limit of 200 characters.";
      throw SGroupException(errout.str());
    }
    unsigned int start = 0;
    unsigned int end = 69;
    for (; length > end; start += 69, end += 69) {
      std::string dataChunk = data.substr(start, end);
      ret << "M  SCD" << FormatV2000IntField(idx)
          << FormatV2000StringField(dataChunk, 69, true, true) << std::endl;
    }
    std::string dataChunk = data.substr(start, end);
    ret << "M  SED" << FormatV2000IntField(idx)
        << FormatV2000StringField(dataChunk, 69, false, true) << std::endl;
  }

  return ret.str();
}

std::string BuildV2000PXALine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasStrProp("PXA")) {
    ret << "M  PXA" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getStrProp("PXA"), 69, false, true);
    ret << std::endl;
  }

  return ret.str();
}

std::string BuildV2000SAPLines(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;
  std::ostringstream temp;

  const std::vector<SGroup::AttachPoint> saps = sgroup->getAttachPoints();

  unsigned int count = 0;
  unsigned int entriesPerLine = 6;
  for (const auto &sap : saps) {
    temp << FormatV2000IntField(sap.aAtom->getIdx());
    if (sap.lvAtom != nullptr) {
      temp << FormatV2000IntField(sap.lvAtom->getIdx());
    } else {
      temp << FormatV2000IntField(0);
    }
    temp << FormatV2000StringField(sap.id, 2, false, true);
    if (++count == entriesPerLine) {
      ret << "M  SAP" << FormatV2000IntField(idx)
          << FormatV2000IntField(entriesPerLine) << temp.str() << std::endl;
      temp.str("");
      count = 0;
    }
  }
  if (count) {
    ret << "M  SAP" << FormatV2000IntField(idx) << FormatV2000IntField(count)
        << temp.str() << std::endl;
  }
  return ret.str();
}

std::string BuildV2000SCLLine(const int idx, const SGroup *sgroup) {
  std::ostringstream ret;

  if (sgroup->hasStrProp("CLASS")) {
    ret << "M  SCL" << FormatV2000IntField(idx);
    ret << FormatV2000StringField(sgroup->getStrProp("CLASS"), 69, false, true);
    ret << std::endl;
  }

  return ret.str();
}

const std::string GetMolFileSGroupInfo(const RWMol &mol) {
  std::ostringstream ret;

  // multiple group per line properties
  ret << BuildV2000STYLines(mol);
  ret << BuildV2000StringPropLines(8, mol, "SUBTYPE", "SST", 3);
  ret << BuildV2000SLBLines(mol);
  ret << BuildV2000StringPropLines(8, mol, "CONNECT", "SCN", 3);
  ret << BuildV2000SDSLines(mol);
  ret << BuildV2000SPLLines(mol);
  ret << BuildV2000StringPropLines(8, mol, "COMPNO", "SNC", 3);
  ret << BuildV2000SBTLines(mol);

  // single group per line properties
  for (auto sGroupItr = mol.beginSGroups(); sGroupItr != mol.endSGroups();
       ++sGroupItr) {
    unsigned int idx = (sGroupItr - mol.beginSGroups()) + 1;
    ret << BuildV2000IdxVectorDataLines(15, idx, "SAL",
                                        (*sGroupItr)->getAtoms());
    ret << BuildV2000IdxVectorDataLines(15, idx, "SPA",
                                        (*sGroupItr)->getPAtoms());
    ret << BuildV2000IdxVectorDataLines(15, idx, "SBL",
                                        (*sGroupItr)->getBonds());
    ret << BuildV2000SMTLine(idx, sGroupItr->get());
    // Write CRS line -- CRS still not supported
    ret << BuildV2000SDILine(idx, sGroupItr->get());
    ret << BuildV2000SBVLine(idx, sGroupItr->get());

    ret << BuildV2000SDTLine(idx, sGroupItr->get());
    ret << BuildV2000SDDLine(idx, sGroupItr->get());
    // SCD/SED must come after SDT
    ret << BuildV2000SCDSEDLines(idx, sGroupItr->get());

    ret << BuildV2000PXALine(idx, sGroupItr->get());
    ret << BuildV2000SAPLines(idx, sGroupItr->get());
    ret << BuildV2000SCLLine(idx, sGroupItr->get());
  }

  return ret.str();
}

template <class T>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const std::vector<T *> &dataVector) {
  return BuildV3000IdxVectorDataBlock(key, dataVector.begin(),
                                      dataVector.end());
}

template <class Iterator>
std::string BuildV3000IdxVectorDataBlock(const std::string &key,
                                         const Iterator &dataVectorBegin,
                                         const Iterator &dataVectorEnd) {
  using T = typename std::iterator_traits<Iterator>::value_type;

  std::ostringstream ret;

  unsigned int size = dataVectorEnd - dataVectorBegin;

  if (size) {
    auto getIdx = [](T element) -> std::string {
      return std::to_string(1 + element->getIdx());
    };

    std::vector<std::string> tempStr(1 + size);
    tempStr[0] = std::to_string(size);

    std::transform(dataVectorBegin, dataVectorEnd, tempStr.begin() + 1, getIdx);

    ret << ' ' << key << "=(" << boost::algorithm::join(tempStr, " ") << ')';
  }

  return ret.str();
}

/* Classify bonds between XBONDS and CBOfindP work on a copy of
 * bonds vector to prevent reordering of original vector */
std::string BuildV3000BondsBlock(std::vector<Bond *> bonds) {
  std::ostringstream ret;

  auto isXBond = [](Bond *bond) {
    return Bond::SGroupBondType::XBOND == bond->getSGroupBondType();
  };

  auto first_cbond = std::partition(bonds.begin(), bonds.end(), isXBond);

  ret << BuildV3000IdxVectorDataBlock("XBONDS", bonds.begin(), first_cbond);
  ret << BuildV3000IdxVectorDataBlock("CBONDS", first_cbond, bonds.end());

  return ret.str();
}

std::string FormatV3000StringPropertyBlock(const std::string &prop,
                                           const SGROUP_SPTR &sgroup) {
  std::ostringstream ret;

  if (sgroup->hasStrProp(prop)) {
    ret << ' ' << prop << '=';
    std::string value = sgroup->getStrProp(prop);
    bool hasSpaces = (value.end() != find(value.begin(), value.end(), ' '));

    if (hasSpaces || value.empty()) {
      ret << "\"";
    }

    ret << value;

    if (hasSpaces || value.empty()) {
      ret << "\"";
    }
  }

  return ret.str();
}

std::string FormatV3000ParentBlock(const SGROUP_SPTR &sgroup) {
  std::ostringstream ret;

  SGROUP_SPTR *parent = sgroup->getParent();

  if (parent) {
    auto mol = sgroup->getOwningMol();
    unsigned int parentIdx = 1 + (*parent)->getIndexInMol();
    ret << " PARENT=" << parentIdx;
  }

  return ret.str();
}

inline std::string FormatV3000DoubleField(double value) {
  return boost::trim_copy(FormatV2000DoubleField(value));
}

std::string FormatV3000BracketBlock(
    const std::vector<SGroup::Bracket> brackets) {
  std::ostringstream ret;

  if (!brackets.empty() && brackets.size() != 2) {
    std::ostringstream errout;
    errout << "SGroup has an unsupported number of brackets ("
           << brackets.size() << ")";
    throw SGroupException(errout.str());
  }

  for (const auto &bracket : brackets) {
    ret << " BRKXYZ=(9";
    for (unsigned int iPoint = 0; iPoint < 2; ++iPoint) {
      ret << ' ' << FormatV3000DoubleField(bracket.at(iPoint).x);
      ret << ' ' << FormatV3000DoubleField(bracket.at(iPoint).y);
      ret << " 0";  // z coordinate is 0 by format specification
    }
    ret << " 0 0 0";  // 3rd point is 0 by format specification
    ret << ")";
  }

  return ret.str();
}

std::string FormatV3000CStateBlock(const std::vector<SGroup::CState> &cstates) {
  std::ostringstream ret;

  for (const auto &cstate : cstates) {
    unsigned int xbondIdx = 1 + cstate.bond->getIdx();
    ret << " CSTATE=(";
    if (cstate.vector) {
      ret << "4 " << xbondIdx;
      ret << ' ' << FormatV3000DoubleField(cstate.vector->x);
      ret << ' ' << FormatV3000DoubleField(cstate.vector->y);
      ret << " 0";
    } else {
      ret << "1 " << xbondIdx;
    }
    ret << ")";
  }

  return ret.str();
}

std::string FormatV3000AttachPointBlock(
    const std::vector<SGroup::AttachPoint> &saps) {
  std::ostringstream ret;

  for (const auto &sap : saps) {
    ret << " SAP=(3 " << (1 + sap.aAtom->getIdx());

    if (sap.aAtom == sap.lvAtom) {
      ret << " aidx";
    } else {
      ret << ' ' << (1 + sap.lvAtom->getIdx());
    }

    ret << ' ' << sap.id << ")";
  }

  return ret.str();
}

const std::string GetV3000MolFileSGroupLines(const unsigned int idx,
                                             const SGROUP_SPTR &sgroup) {
  std::ostringstream os;

  os << ' ' << idx;
  os << ' ' << sgroup->getType();
  os << ' ' << sgroup->getId();

  os << BuildV3000IdxVectorDataBlock("ATOMS", sgroup->getAtoms());
  os << BuildV3000BondsBlock(sgroup->getBonds());
  os << BuildV3000IdxVectorDataBlock("PATOMS", sgroup->getPAtoms());
  os << FormatV3000StringPropertyBlock("SUBTYPE", sgroup);
  os << FormatV3000StringPropertyBlock("MULT", sgroup);
  os << FormatV3000StringPropertyBlock("CONNECT", sgroup);
  os << FormatV3000ParentBlock(sgroup);
  os << FormatV3000StringPropertyBlock("COMPNO", sgroup);
  // XBHEAD -> part of V2000 CRS, not supported yet
  // XBCORR -> part of V2000 CRS, not supported yet
  os << FormatV3000StringPropertyBlock("LABEL", sgroup);
  os << FormatV3000BracketBlock(sgroup->getBrackets());
  os << FormatV3000StringPropertyBlock("ESTATE", sgroup);
  os << FormatV3000CStateBlock(sgroup->getCStates());
  os << FormatV3000StringPropertyBlock("FIELDNAME", sgroup);
  os << FormatV3000StringPropertyBlock("FIELDINFO", sgroup);
  os << FormatV3000StringPropertyBlock("FIELDDISP", sgroup);
  os << FormatV3000StringPropertyBlock("ESTATE", sgroup);
  os << FormatV3000StringPropertyBlock("QUERYTYPE", sgroup);
  os << FormatV3000StringPropertyBlock("QUERYOP", sgroup);
  os << FormatV3000StringPropertyBlock("CLASS", sgroup);
  os << FormatV3000AttachPointBlock(sgroup->getAttachPoints());
  os << FormatV3000StringPropertyBlock("BRKTYP", sgroup);
  os << FormatV3000StringPropertyBlock("SEQID", sgroup);

  std::string sGroupBlock = os.str();
  unsigned int length = sGroupBlock.size();
  os.str("");

  unsigned int start = 0;
  while (length - start > 73) {
    os << "M  V30 " << sGroupBlock.substr(start, start + 72) << '-'
       << std::endl;
    start += 72;
  }
  os << "M  V30 " << sGroupBlock.substr(start, start + 73) << std::endl;

  return os.str();
}

//------------------------------------------------
//
//  gets a mol block as a string
//
//------------------------------------------------
std::string outputMolToMolBlock(const RWMol &tmol, int confId,
                                bool forceV3000) {
  std::string res;

  bool isV3000;
  unsigned int nAtoms, nBonds, nLists, nSGroups, chiralFlag, nsText,
      nRxnComponents;
  unsigned int nReactants, nProducts, nIntermediates;
  nAtoms = tmol.getNumAtoms();
  nBonds = tmol.getNumBonds();
  nLists = 0;
  nSGroups = tmol.getNumSGroups();

  chiralFlag = 0;
  nsText = 0;
  nRxnComponents = 0;
  nReactants = 0;
  nProducts = 0;
  nIntermediates = 0;

  tmol.getPropIfPresent(common_properties::_MolFileChiralFlag, chiralFlag);

  const Conformer *conf;
  if (confId < 0 && tmol.getNumConformers() == 0) {
    conf = nullptr;
  } else {
    conf = &(tmol.getConformer(confId));
  }

  std::string text;
  if (tmol.getPropIfPresent(common_properties::_Name, text)) {
    res += text;
  }
  res += "\n";

  // info
  if (tmol.getPropIfPresent(common_properties::MolFileInfo, text)) {
    res += text;
  } else {
    std::stringstream ss;
    ss << "  " << std::setw(8) << "RDKit";
    ss << std::setw(10) << "";
    if (conf) {
      if (conf->is3D()) {
        ss << "3D";
      } else {
        ss << common_properties::TWOD;
      }
    }
    res += ss.str();
  }
  res += "\n";
  // comments
  if (tmol.getPropIfPresent(common_properties::MolFileComments, text)) {
    res += text;
  }
  res += "\n";

  if (forceV3000)
    isV3000 = true;
  else
    isV3000 = (nAtoms > 999) || (nBonds > 999);

  // the counts line:
  std::stringstream ss;
  if (isV3000) {
    // All counts in the V3000 info line should be 0
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << std::setw(3) << 0;
    ss << "999 V3000\n";
  } else {
    ss << std::setw(3) << nAtoms;
    ss << std::setw(3) << nBonds;
    ss << std::setw(3) << nLists;
    ss << std::setw(3) << nSGroups;
    ss << std::setw(3) << chiralFlag;
    ss << std::setw(3) << nsText;
    ss << std::setw(3) << nRxnComponents;
    ss << std::setw(3) << nReactants;
    ss << std::setw(3) << nProducts;
    ss << std::setw(3) << nIntermediates;
    ss << "999 V2000\n";
  }
  res += ss.str();

  if (!isV3000) {
    // V2000 output.
    for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
         atomIt != tmol.endAtoms(); ++atomIt) {
      res += GetMolFileAtomLines(*atomIt, conf);
      res += "\n";
    }

    INT_MAP_INT wedgeBonds = pickBondsToWedge(tmol);
    for (ROMol::ConstBondIterator bondIt = tmol.beginBonds();
         bondIt != tmol.endBonds(); ++bondIt) {
      res += GetMolFileBondLines(*bondIt, wedgeBonds, conf);
      res += "\n";
    }

    res += GetMolFileChargeInfo(tmol);
    res += GetMolFileRGroupInfo(tmol);
    res += GetMolFileQueryInfo(tmol);
    res += GetMolFileAliasInfo(tmol);
    res += GetMolFileZBOInfo(tmol);

    res += GetMolFileSGroupInfo(tmol);

    // FIX: R-group logic, SGroups and 3D features etc.
  } else {
    // V3000 output.
    res += "M  V30 BEGIN CTAB\n";
    std::stringstream ss;
    ss << "M  V30 COUNTS " << nAtoms << " " << nBonds << " " << nSGroups
       << " 0 " << chiralFlag << "\n";
    //      |
    //      num3DConstraints (not implemented)

    res += ss.str();

    res += "M  V30 BEGIN ATOM\n";
    for (ROMol::ConstAtomIterator atomIt = tmol.beginAtoms();
         atomIt != tmol.endAtoms(); ++atomIt) {
      res += GetV3000MolFileAtomLines(*atomIt, conf);
      res += "\n";
    }
    res += "M  V30 END ATOM\n";

    if (tmol.getNumBonds()) {
      res += "M  V30 BEGIN BOND\n";
      INT_MAP_INT wedgeBonds = pickBondsToWedge(tmol);
      for (ROMol::ConstBondIterator bondIt = tmol.beginBonds();
           bondIt != tmol.endBonds(); ++bondIt) {
        res += GetV3000MolFileBondLines(*bondIt, wedgeBonds, conf);
        res += "\n";
      }
      res += "M  V30 END BOND\n";
    }

    if (tmol.getNumSGroups()) {
      res += "M  V30 BEGIN SGROUP\n";
      for (ROMol::ConstSGroupIterator sGroupIt = tmol.beginSGroups();
           sGroupIt != tmol.endSGroups(); ++sGroupIt) {
        unsigned int idx = 1 + (sGroupIt - tmol.beginSGroups());
        res += GetV3000MolFileSGroupLines(idx, *sGroupIt);
      }
      res += "M  V30 END SGROUP\n";
    }

    res += "M  V30 END CTAB\n";
  }
  res += "M  END\n";
  return res;
}

std::string MolToMolBlock(const ROMol &mol, bool includeStereo, int confId,
                          bool kekulize, bool forceV3000) {
  RDUNUSED_PARAM(includeStereo);
  RDKit::Utils::LocaleSwitcher switcher;
  ROMol tromol(mol);
  RWMol &trwmol = static_cast<RWMol &>(tromol);
  // NOTE: kekulize the molecule before writing it out
  // because of the way mol files handle aromaticity
  if (trwmol.needsUpdatePropertyCache()) {
    trwmol.updatePropertyCache(false);
  }
  if (kekulize) MolOps::Kekulize(trwmol);

  if (includeStereo && !trwmol.getNumConformers()) {
    // generate coordinates so that the stereo we generate makes sense
    RDDepict::compute2DCoords(trwmol);
  }
#if 0
    if(includeStereo){
      // assign "any" status to any stereo bonds that are not
      // marked with "E" or "Z" code - these bonds need to be explictly written
      // out to the mol file
      MolOps::findPotentialStereoBonds(trwmol);
      // now assign stereo code if any have been specified by the directions on
      // single bonds
      MolOps::assignStereochemistry(trwmol);
    }
#endif
  const RWMol &tmol = const_cast<RWMol &>(trwmol);

  try {
    return outputMolToMolBlock(tmol, confId, forceV3000);
  } catch (RequiresV3000Exception) {
    return outputMolToMolBlock(tmol, confId, true);
  }
}

//------------------------------------------------
//
//  Dump a molecule to a file
//
//------------------------------------------------
void MolToMolFile(const ROMol &mol, const std::string &fName,
                  bool includeStereo, int confId, bool kekulize,
                  bool forceV3000) {
  auto *outStream = new std::ofstream(fName.c_str());
  if (!outStream || !(*outStream) || outStream->bad()) {
    std::ostringstream errout;
    errout << "Bad output file " << fName;
    throw BadFileException(errout.str());
  }
  std::string outString =
      MolToMolBlock(mol, includeStereo, confId, kekulize, forceV3000);
  *outStream << outString;
  delete outStream;
}
}  // namespace RDKit
