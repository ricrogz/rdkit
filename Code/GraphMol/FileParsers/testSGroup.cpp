//
//  Copyright (C) 2002-2018 Greg Landrum and Rational Discovery LLC
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/test.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include "FileParsers.h"
#include "MolSupplier.h"
#include "MolWriters.h"

#include <memory>
#include <cstdlib>

using namespace RDKit;

template <class T>
void testSet(const std::vector<T> &groupVector,
             const std::vector<unsigned int> &reference) {
  TEST_ASSERT(groupVector.size() == reference.size());

  for (auto sgItr = groupVector.begin(), refItr = reference.begin();
       refItr != reference.end(); ++sgItr, ++refItr) {
    TEST_ASSERT(1 + (*sgItr)->getIdx() == *refItr);
  }
}

void testBrackets(
    const std::vector<SGroup::Bracket> &brackets,
    const std::vector<std::array<std::array<double, 3>, 3>> &reference) {
  TEST_ASSERT(brackets.size() == 2);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        TEST_ASSERT(abs(brackets[i][j][k] - reference[i][j][k]) < 1.e-6);
      }
    }
  }
}

std::string outputFileNameBuilder(const std::string &rdbase,
                                  const std::string &suffix, bool forceV3000) {
  std::ostringstream ret;
  ret << rdbase << "/Code/GraphMol/FileParsers/test_data/output_";

  if (forceV3000) {
    ret << "V3k_";
  }

  ret << suffix;

  return ret.str();
}

void testSGroup(const std::string &rdbase, bool forceV3000) {
  BOOST_LOG(rdInfoLog) << " ----------> Testing SGroups " << std::endl;

  BOOST_LOG(rdInfoLog) << " ----------> V2000 -- 1 " << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.mol";

    auto m = std::shared_ptr<RWMol>(MolFileToMol(fName));

    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 10);

    SGROUP_SPTR sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "MON");

    std::vector<unsigned int> atoms_reference = {2, 3, 4, 1, 5};

    testSet(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference =
        {};  // No bonds defined in this mol
    testSet(sgroup->getBonds(), bonds_reference);

    std::vector<std::array<std::array<double, 3>, 3>> brackets_reference = {{
        {{{{-3.9679, -0.1670, 0.}}, {{-3.9679, 2.1705, 0.}}, {{0., 0., 0.}}}},
        {{{{-0.7244, 2.1705, 0.}}, {{-0.7244, -0.1670, 0.}}, {{0., 0., 0.}}}},
    }};
    testBrackets(sgroup->getBrackets(), brackets_reference);

    std::string ofile =
        outputFileNameBuilder(rdbase, "Issue3432136_1.mol", true);
    auto writer = SDWriter(ofile);
    writer.setForceV3000(forceV3000);

    writer.write(*m);
    writer.close();
  }

  BOOST_LOG(rdInfoLog) << " ----------> V3000 -- 1 " << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_1.v3k.mol";
    auto m = std::shared_ptr<RWMol>(MolFileToMol(fName));
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 1);

    SGROUP_SPTR sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "MON");

    std::vector<unsigned int> atoms_reference = {2, 3, 4, 1, 5};
    testSet(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference =
        {};  // No bonds defined in this mol
    testSet(sgroup->getBonds(), bonds_reference);

    std::string ofile =
        outputFileNameBuilder(rdbase, "Issue3432136_1.v3k.mol", true);
    auto writer = SDWriter(ofile);
    writer.setForceV3000(forceV3000);

    writer.write(*m);
    writer.close();
  }

  BOOST_LOG(rdInfoLog) << " ----------> V3000 -- 2 " << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.v3k.mol";
    auto m = std::shared_ptr<RWMol>(MolFileToMol(fName));
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 2);

    SGROUP_SPTR sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "SUP");

    std::vector<unsigned int> atoms_reference = {6, 7, 8, 9, 11, 12};
    testSet(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference = {5};
    testSet(sgroup->getBonds(), bonds_reference);

    auto *bond = sgroup->getBonds()[0];
    TEST_ASSERT(sgroup->getBondType(bond) == SGroup::BondType::XBOND);

    std::string ofile =
        outputFileNameBuilder(rdbase, "Issue3432136_2.v3k.mol", true);
    auto writer = SDWriter(ofile);
    writer.setForceV3000(forceV3000);

    writer.write(*m);
    writer.close();
  }

  BOOST_LOG(rdInfoLog) << " ----------> V2000 -- 2 " << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/Issue3432136_2.mol";
    auto m = std::shared_ptr<RWMol>(MolFileToMol(fName));
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 10);

    SGROUP_SPTR sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "SUP");

    std::vector<unsigned int> atoms_reference = {6, 7, 8, 9, 11, 12};
    testSet(sgroup->getAtoms(), atoms_reference);

    std::vector<unsigned int> bonds_reference = {5, 11};
    testSet(sgroup->getBonds(), bonds_reference);

    auto *bond = sgroup->getBonds()[0];
    TEST_ASSERT(sgroup->getBondType(bond) == SGroup::BondType::XBOND);

    std::string ofile =
        outputFileNameBuilder(rdbase, "Issue3432136_2.mol", true);
    auto writer = SDWriter(ofile);
    writer.setForceV3000(forceV3000);

    writer.write(*m);
    writer.close();
  }

  BOOST_LOG(rdInfoLog) << " ----------> sgroup.sdf " << std::endl;
  {
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/sgroup.sdf";

    SDMolSupplier sdsup(fName);
    std::shared_ptr<ROMol> m;
    unsigned int i = 0;

    std::string ofile = outputFileNameBuilder(rdbase, "sgroup.sdf", true);
    auto writer = SDWriter(ofile);
    writer.setForceV3000(forceV3000);

    BOOST_LOG(rdInfoLog) << " ----> mol " << ++i << std::endl;
    m.reset(sdsup.next());
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 1);

    SGROUP_SPTR sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "DAT");

    std::vector<unsigned int> atoms_reference_1 = {1};
    testSet(sgroup->getAtoms(), atoms_reference_1);

    std::vector<unsigned int> bonds_reference_1 = {};
    testSet(sgroup->getBonds(), bonds_reference_1);

    writer.write(*m);

    BOOST_LOG(rdInfoLog) << " ----> mol " << ++i << std::endl;
    m.reset(sdsup.next());
    TEST_ASSERT(m);
    TEST_ASSERT(m->getNumSGroups() == 1);

    sgroup = m->getSGroup(0);

    TEST_ASSERT(sgroup->getType() == "DAT");

    std::vector<unsigned int> atoms_reference_2 = {1};
    testSet(sgroup->getAtoms(), atoms_reference_2);

    std::vector<unsigned int> bonds_reference_2 = {};
    testSet(sgroup->getBonds(), bonds_reference_2);

    writer.write(*m);

    while (!sdsup.atEnd()) {
      BOOST_LOG(rdInfoLog) << " ----> mol " << ++i << std::endl;
      m.reset(sdsup.next());
      writer.write(*m);
    }
    writer.close();
  }
  BOOST_LOG(rdInfoLog) << " Finished <---------- " << std::endl;
}

int main() {
  const char *rdbase = getenv("RDBASE");
  if (rdbase == nullptr) {
    std::cerr << "\n\n RDBASE has not been set, aborting.\n\n";
    return 1;
  }

  RDLog::InitLogs();
  testSGroup(std::string(rdbase), false);  // V2000
  testSGroup(std::string(rdbase), true);   // V3000

  return 0;
}
