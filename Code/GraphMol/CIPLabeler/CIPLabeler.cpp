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
#include <memory>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "GraphMol/Chirality.h"
#include "GraphMol/RDKitBase.h"
#include <RDGeneral/ControlCHandler.h>
#include <RDGeneral/Exceptions.h>

#include "CIPLabeler.h"
#include "CIPMol.h"
#include "configs/Sp2Bond.h"
#include "configs/Tetrahedral.h"
#include "configs/AtropisomerBond.h"

#include "rules/Rules.h"
#include "rules/Rule1a.h"
#include "rules/Rule1b.h"
#include "rules/Rule2.h"
#include "rules/Rule3.h"
#include "rules/Rule4a.h"
#include "rules/Rule4b.h"
#include "rules/Rule4c.h"
#include "rules/Rule5New.h"
#include "rules/Rule6.h"

namespace RDKit {
namespace CIPLabeler {

namespace {

// constitutional rules
const Rules constitutional_rules({new Rule1a, new Rule1b, new Rule2});

// all rules (require aux calc)
const Rules all_rules({new Rule1a, new Rule1b, new Rule2, new Rule3, new Rule4a,
                       new Rule4b, new Rule4c, new Rule5New, new Rule6});

struct ConfigEntry {
  std::unique_ptr<Configuration> config;
  bool selected = false;
  bool constitutionalPassComplete = false;
};

using ConfigList = std::vector<ConfigEntry>;

bool isSelected(const boost::dynamic_bitset<> &selection, size_t index) {
  return index < selection.size() && selection.test(index);
}

ConfigList findConfigs(CIPMol &mol, const boost::dynamic_bitset<> &atoms,
                       const boost::dynamic_bitset<> &bonds) {
  ConfigList configs;

  // All configurations are required here, including unselected ones: they may
  // provide auxiliary descriptors needed to label a selected configuration.
  for (size_t index = 0; index < mol.getNumAtoms(); ++index) {
    auto atom = mol.getAtom(index);
    auto chiraltag = atom->getChiralTag();
    if (chiraltag == Atom::CHI_TETRAHEDRAL_CW ||
        chiraltag == Atom::CHI_TETRAHEDRAL_CCW) {
      auto cfg = std::make_unique<Tetrahedral>(mol, atom);
      if (cfg->getCarriers().size() == 4) {
        configs.push_back({std::move(cfg), isSelected(atoms, index)});
      }
    }
  }

  for (size_t index = 0; index < mol.getNumBonds(); ++index) {
    auto bond = mol.getBond(index);

    auto bond_cfg = bond->getStereo();
    switch (bond_cfg) {
      case Bond::STEREOE:
        bond_cfg = Bond::STEREOTRANS;
        break;
      case Bond::STEREOZ:
        bond_cfg = Bond::STEREOCIS;
        break;
      default:
        break;
    }
    switch (bond_cfg) {
      case Bond::STEREOTRANS:
      case Bond::STEREOCIS: {
        if (bond->getBondType() != Bond::DOUBLE) {
          break;
        }
        auto cfg = std::make_unique<Sp2Bond>(mol, bond, bond->getBeginAtom(),
                                             bond->getEndAtom(), bond_cfg);
        if (cfg->getCarriers().size() == 2) {
          configs.push_back({std::move(cfg), isSelected(bonds, index)});
        }
      } break;

      case Bond::STEREOATROPCCW:
      case Bond::STEREOATROPCW: {
        auto cfg = std::make_unique<AtropisomerBond>(
            mol, bond, bond->getBeginAtom(), bond->getEndAtom(), bond_cfg);
        if (cfg->getCarriers().size() == 2) {
          configs.push_back({std::move(cfg), isSelected(bonds, index)});
        }
      } break;

      default:
        break;
    }
  }

  boost::dynamic_bitset<> configurationFoci(mol.getNumAtoms());
  for (const auto &entry : configs) {
    for (const auto focus : entry.config->getFoci()) {
      configurationFoci.set(focus->getIdx());
    }
  }
  mol.setConfigurationFoci(std::move(configurationFoci));

  return configs;
}

bool labelAux(ConfigList &configs, const Rules &rules, ConfigEntry &center) {
  using Node_Cfg_Pair = std::pair<Node *, Configuration *>;
  std::vector<Node_Cfg_Pair> aux;

  auto &digraph = center.config->getDigraph();
  if (digraph.getCurrentRoot() != digraph.getOriginalRoot()) {
    digraph.changeRoot(digraph.getOriginalRoot());
  }

  // Ordinarily retain the existing reached-focus optimization. If an exact
  // constitutional automorphism ended a comparison without expanding its
  // remote paths, include every configuration: those unseen descriptors may
  // be precisely what Rules 3-6 need to break the constitutional symmetry.
  const bool includeUnseen = digraph.usedConstitutionalRootEquivalence();
  std::vector<unsigned char> eligible(configs.size());
  boost::dynamic_bitset<> targetFoci(digraph.getMol().getNumAtoms());
  for (std::size_t i = 0; i < configs.size(); ++i) {
    const auto &entry = configs[i];
    if (entry.config.get() == center.config.get()) {
      continue;
    }
    const auto &config = entry.config;
    // FIXME: specific to each descriptor
    const auto &foci = config->getFoci();

    if (!includeUnseen && std::ranges::none_of(foci, [&](auto focus) {
          return digraph.seenAtom(focus);
        })) {
      continue;
    }
    eligible[i] = 1u;
    targetFoci.set(foci[0]->getIdx());
  }

  // Collect every eligible focus in one target-guided traversal. A separate
  // getNodes() walk per configuration repeatedly scanned the same dense cage
  // occurrence tree.
  const auto occurrences = digraph.getNodesForAuxiliaryLabeling(targetFoci);
  for (std::size_t i = 0; i < configs.size(); ++i) {
    if (!eligible[i]) {
      continue;
    }
    const auto &config = configs[i].config;
    const auto &foci = config->getFoci();
    for (const auto &node : occurrences[foci[0]->getIdx()]) {
      if (node->isDuplicate()) {
        continue;
      }
      auto low = node;
      if (foci.size() == 2) {
        for (const auto &edge : node->getEdges(foci[1])) {
          const auto &other_node = edge->getOther(node);
          if (other_node->getDistance() < node->getDistance()) {
            low = other_node;
          }
        }
      }
      if (!low->isDuplicate()) {
        aux.emplace_back(low, config.get());
      }
    }
  }

  auto farthest = [](const Node_Cfg_Pair &a, const Node_Cfg_Pair &b) {
    return a.first->getDistance() > b.first->getDistance();
  };
  // Java's Collections.sort() is stable, so preserve configuration discovery
  // order for configurations at the same distance.
  std::stable_sort(aux.begin(), aux.end(), farthest);

  // Using a boost::unordered_map because it is more performant
  // than the STL version.
  boost::unordered_map<Node *, Descriptor> queue;
  for (std::size_t begin = 0; begin < aux.size();) {
    std::size_t end = begin + 1;
    while (end < aux.size() &&
           aux[end].first->getDistance() == aux[begin].first->getDistance()) {
      ++end;
    }

    // Descriptors are immutable within one distance sphere. Keep exact
    // comparison and sort results across all occurrence labels in that sphere,
    // then clear the session before committing labels for the next sphere.
    {
      const SequenceRule::ComparisonSession comparisonSession;
      for (auto pos = begin; pos < end; ++pos) {
        const auto &node = aux[pos].first;
        const auto &config = aux[pos].second;
        auto label = config->label(node, digraph, rules);
        // Match Java HashMap.put(): a later configuration at this distance wins
        // if multiple configurations map to the same digraph node.
        queue[node] = label;
      }
    }

    for (const auto &entry : queue) {
      entry.first->setAux(entry.second);
    }
    queue.clear();
    begin = end;
  }

  return true;
}

// The chiral centers in current rdkit examples that can be resolved using only
// the constitutional rules average about 8 iterations (the highest count is
// 1039, in one of the examples in the CIP validation suite). We use 2000 as
// threshold to allow some margin.
constexpr unsigned int constitutionalRuleTimeout = 2000;

struct IterationBudget {
  void reset(unsigned int maxRecursiveIterations) {
    hasGlobalLimit = maxRecursiveIterations != 0;
    remainingGlobal = maxRecursiveIterations;
    inPreliminaryPass = false;
    remainingPreliminary = 0;
  }

  void beginPreliminaryPass() {
    inPreliminaryPass = true;
    remainingPreliminary = constitutionalRuleTimeout;
  }

  void endPreliminaryPass() {
    inPreliminaryPass = false;
    remainingPreliminary = 0;
  }

  bool consume() {
    if ((inPreliminaryPass && remainingPreliminary == 0) ||
        (hasGlobalLimit && remainingGlobal == 0)) {
      return false;
    }
    if (inPreliminaryPass) {
      --remainingPreliminary;
    }
    if (hasGlobalLimit) {
      --remainingGlobal;
    }
    return true;
  }

  bool hasGlobalLimit = false;
  unsigned int remainingGlobal = 0;
  bool inPreliminaryPass = false;
  unsigned int remainingPreliminary = 0;
};

thread_local IterationBudget iterationBudget;

class ScopedIterationBudget {
 public:
  explicit ScopedIterationBudget(unsigned int maxRecursiveIterations)
      : d_previous{iterationBudget} {
    iterationBudget.reset(maxRecursiveIterations);
  }

  ~ScopedIterationBudget() { iterationBudget = d_previous; }

 private:
  IterationBudget d_previous;
};

class ScopedPreliminaryBudget {
 public:
  ScopedPreliminaryBudget() { iterationBudget.beginPreliminaryPass(); }
  ~ScopedPreliminaryBudget() { iterationBudget.endPreliminaryPass(); }
};

void label(ConfigList &configs, unsigned int maxRecursiveIterations) {
  const ScopedIterationBudget callBudget(maxRecursiveIterations);

  // First, if the specified number of iterations allows it, run all centers
  // through a fast pass with the constitutional rules allow easy stuff to be
  // resolved.
  for (auto &entry : configs) {
    if (!entry.selected) {
      continue;
    }
    auto &conf = entry.config;
    // Make sure this stereo center has no label
    conf->resetPrimaryLabel();

    try {
      const ScopedPreliminaryBudget preliminaryBudget;
      auto desc = conf->label(constitutional_rules);
      entry.constitutionalPassComplete = true;
      if (desc != Descriptor::UNKNOWN) {
        conf->setPrimaryLabel(desc);
      }
    } catch (const MaxIterationsExceeded &) {
    }
  }

  // try again on everything that hasn't been resolved yet
  for (auto &entry : configs) {
    if (!entry.selected) {
      continue;
    }
    auto &conf = entry.config;
    if (conf->hasPrimaryLabel()) {
      // already resolved!
      continue;
    }

    auto desc = Descriptor::UNKNOWN;
    // A completed constitutional pass is deterministic. Repeating it before
    // adding auxiliary descriptors cannot change an UNKNOWN result.
    if (!entry.constitutionalPassComplete) {
      desc = conf->label(constitutional_rules);
    }
    if (desc != Descriptor::UNKNOWN) {
      conf->setPrimaryLabel(desc);
    } else {
      if (labelAux(configs, all_rules, entry)) {
        desc = conf->label(all_rules);

        if (desc != Descriptor::UNKNOWN) {
          conf->setPrimaryLabel(desc);
        }
      }
    }
  }
}

void validateSelection(const boost::dynamic_bitset<> &selection,
                       size_t expectedSize, const char *kind) {
  for (auto index = selection.find_first();
       index != boost::dynamic_bitset<>::npos;
       index = selection.find_next(index)) {
    if (index >= expectedSize) {
      std::ostringstream msg;
      msg << "CIP " << kind << " selection contains out-of-range index "
          << index << "; molecule has " << expectedSize << ' ' << kind
          << " entries";
      throw ValueErrorException(msg.str());
    }
  }
}

bool isFullSelection(const boost::dynamic_bitset<> &selection,
                     size_t expectedSize) {
  // validateSelection() guarantees that every set bit denotes a molecule
  // object, so selecting expectedSize objects means that all are selected.
  return selection.count() == expectedSize;
}

template <typename T>
void clearCIPProperties(T *object) {
  object->clearProp(common_properties::_CIPCode);
  object->clearProp(common_properties::_CIPNeighborOrder);
}

void clearSelectedCIPProperties(ROMol &mol,
                                const boost::dynamic_bitset<> &atoms,
                                const boost::dynamic_bitset<> &bonds,
                                bool fullSelection) {
  if (fullSelection) {
    for (auto atom : mol.atoms()) {
      clearCIPProperties(atom);
    }
    for (auto bond : mol.bonds()) {
      clearCIPProperties(bond);
    }
    return;
  }

  for (auto index = atoms.find_first(); index != boost::dynamic_bitset<>::npos;
       index = atoms.find_next(index)) {
    clearCIPProperties(mol.getAtomWithIdx(index));
  }
  for (auto index = bonds.find_first(); index != boost::dynamic_bitset<>::npos;
       index = bonds.find_next(index)) {
    clearCIPProperties(mol.getBondWithIdx(index));
  }
}

}  // namespace

void assignCIPLabels(ROMol &mol, const boost::dynamic_bitset<> &atoms,
                     const boost::dynamic_bitset<> &bonds,
                     unsigned int maxRecursiveIterations) {
  ControlCHandler hdlr;

  validateSelection(atoms, mol.getNumAtoms(), "atom");
  validateSelection(bonds, mol.getNumBonds(), "bond");

  const bool fullSelection = isFullSelection(atoms, mol.getNumAtoms()) &&
                             isFullSelection(bonds, mol.getNumBonds());

  // reset the mark, for the case that this fails
  mol.clearProp(common_properties::_CIPComputed);
  clearSelectedCIPProperties(mol, atoms, bonds, fullSelection);
  CIPMol cipmol{mol};
  auto configs = findConfigs(cipmol, atoms, bonds);

  try {
    label(configs, maxRecursiveIterations);
  } catch (const ControlCCaught &) {
  }
  if (hdlr.getGotSignal()) {
    BOOST_LOG(rdWarningLog)
        << "Interrupted, cancelling CIP label calculation" << std::endl;
    return;
  }

  if (fullSelection) {
    const bool computed = true;
    mol.setProp(common_properties::_CIPComputed, true, computed);
  }
}

void assignCIPLabels(ROMol &mol, unsigned int maxRecursiveIterations) {
  boost::dynamic_bitset<> atoms(mol.getNumAtoms());
  boost::dynamic_bitset<> bonds(mol.getNumBonds());
  atoms.set();
  bonds.set();
  assignCIPLabels(mol, atoms, bonds, maxRecursiveIterations);
}

}  // namespace CIPLabeler

namespace CIPLabeler_detail {

bool decrementRemainingCallCountAndCheck() {
  return CIPLabeler::iterationBudget.consume();
}

}  // namespace CIPLabeler_detail
}  // namespace RDKit
