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
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <tuple>
#include <unordered_set>
#include <utility>

#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/ControlCHandler.h>

#include "CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

namespace {

constexpr std::size_t MAX_AUTOMORPHISM_EVIDENCE_PER_KIND = 8;
constexpr std::size_t AUTOMORPHISM_EVIDENCE_KEYS_PER_GENERATION = 2048;
constexpr std::size_t AUTOMORPHISM_EVIDENCE_ATOM_INDICES_PER_GENERATION =
    262144;
constexpr std::size_t MIN_AUTOMORPHISM_SEARCH_CALLBACKS = 65536;
constexpr std::size_t AUTOMORPHISM_SEARCH_CALLBACKS_PER_ATOM = 4096;
constexpr std::size_t MAX_AUTOMORPHISM_SEARCH_CALLBACKS = 4194304;
constexpr std::size_t COMPONENT_AUTOMORPHISM_SEARCH_MULTIPLIER = 4;

struct AutomorphismSearchCallbackLimitReached {};
struct ComponentAutomorphismSearchCallbackLimitReached {};

struct DistanceSignature {
  unsigned int symmetryClass;
  unsigned int rootDistance;
  unsigned int branchDistance;

  bool operator==(const DistanceSignature &other) const = default;
};

template <typename T>
void hashCombine(std::size_t &seed, const T &value) {
  const auto hash = std::hash<T>{}(value);
  seed ^= hash + 0x9e3779b9u + (seed << 6) + (seed >> 2);
}

void radixSortDistanceSignatures(std::vector<DistanceSignature> &signatures,
                                 std::size_t keyCount) {
  PRECONDITION(keyCount != 0u, "distance-signature key range is empty")
  std::vector<DistanceSignature> scratch(signatures.size());
  std::vector<std::size_t> offsets(keyCount + 1u);
  const auto sortBy = [&](auto getKey) {
    std::ranges::fill(offsets, 0u);
    for (const auto &signature : signatures) {
      const auto key = static_cast<std::size_t>(getKey(signature));
      PRECONDITION(key < keyCount, "distance-signature key is out of range")
      ++offsets[key + 1u];
    }
    for (std::size_t i = 1; i < offsets.size(); ++i) {
      offsets[i] += offsets[i - 1u];
    }
    for (const auto &signature : signatures) {
      scratch[offsets[getKey(signature)]++] = signature;
    }
    signatures.swap(scratch);
  };

  // Stable least-significant-key passes produce an exact lexicographic order
  // without node allocations or comparison-sort behavior on large components.
  sortBy([](const auto &signature) { return signature.branchDistance; });
  sortBy([](const auto &signature) { return signature.rootDistance; });
  sortBy([](const auto &signature) { return signature.symmetryClass; });
}

using AtomPropertySignature =
    std::tuple<unsigned int, unsigned int, int, unsigned int, unsigned int,
               unsigned int>;
using NeighborColor = std::tuple<int, bool, unsigned int>;

struct RefinementSignature {
  unsigned int oldColor;
  std::vector<NeighborColor> neighbors;

  bool operator==(const RefinementSignature &other) const = default;
  bool operator<(const RefinementSignature &other) const {
    return std::tie(oldColor, neighbors) <
           std::tie(other.oldColor, other.neighbors);
  }
};

AtomPropertySignature getAtomPropertySignature(const Atom &atom) {
  return {atom.getAtomicNum(),    atom.getIsotope(),
          atom.getFormalCharge(), atom.getNumRadicalElectrons(),
          atom.getTotalNumHs(),   atom.getDegree()};
}

bool atomsHaveSameConstitutionalProperties(const Atom &first,
                                           const Atom &second) {
  return first.getAtomicNum() == second.getAtomicNum() &&
         first.getIsotope() == second.getIsotope() &&
         first.getFormalCharge() == second.getFormalCharge() &&
         first.getNumRadicalElectrons() == second.getNumRadicalElectrons() &&
         first.getTotalNumHs() == second.getTotalNumHs() &&
         first.getDegree() == second.getDegree();
}

bool maskContainsAtom(std::span<const std::uint64_t> mask,
                      unsigned int atomIdx) {
  const auto word = atomIdx / 64u;
  return word < mask.size() &&
         (mask[word] & (std::uint64_t{1} << (atomIdx % 64u))) != 0u;
}

bool fixedSetContainsAtom(std::span<const std::uint64_t> mask,
                          std::span<const unsigned int> addedAtoms,
                          unsigned int atomIdx) {
  return maskContainsAtom(mask, atomIdx) ||
         std::ranges::find(addedAtoms, atomIdx) != addedAtoms.end();
}

bool atomSetIsDisjointFromMask(std::span<const unsigned int> atomSet,
                               std::span<const std::uint64_t> mask,
                               std::span<const unsigned int> addedAtoms) {
  for (const auto atomIdx : atomSet) {
    if (fixedSetContainsAtom(mask, addedAtoms, atomIdx)) {
      return false;
    }
  }
  return true;
}

bool atomSetIsSubsetOfMask(std::span<const unsigned int> atomSet,
                           std::span<const std::uint64_t> mask,
                           std::span<const unsigned int> addedAtoms) {
  for (const auto atomIdx : atomSet) {
    if (!fixedSetContainsAtom(mask, addedAtoms, atomIdx)) {
      return false;
    }
  }
  return true;
}

bool atomSetIsSubset(std::span<const unsigned int> subset,
                     std::span<const unsigned int> superset) {
  return std::ranges::includes(superset, subset);
}

bool hasSupportedConstitutionalBondType(const Bond &bond) {
  switch (bond.getBondType()) {
    case Bond::ZERO:
    case Bond::HYDROGEN:
    case Bond::SINGLE:
    case Bond::DOUBLE:
    case Bond::TRIPLE:
    case Bond::QUADRUPLE:
    case Bond::QUINTUPLE:
    case Bond::HEXTUPLE:
    case Bond::AROMATIC:
      return true;
    default:
      return false;
  }
}

}  // namespace

std::size_t CIPMol::ConstitutionalAutomorphismKeyHash::operator()(
    const ConstitutionalAutomorphismKey &key) const {
  std::size_t result = 0;
  hashCombine(result, key.root);
  hashCombine(result, key.first);
  hashCombine(result, key.second);
  return result;
}

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
  d_configuration_atom_sets.clear();
  d_configuration_branch_cache.assign(2u * getNumBonds(), 0u);
}

void CIPMol::setConfigurationData(
    boost::dynamic_bitset<> foci,
    std::vector<ConfigurationAtomSet> configurationAtomSets) {
  PRECONDITION(foci.size() == getNumAtoms(),
               "configuration focus bitset has the wrong size")
  boost::dynamic_bitset<> recordedFoci(getNumAtoms());
  for (auto &atomSet : configurationAtomSets) {
    PRECONDITION(!atomSet.foci.empty(), "configuration has no focus")
    PRECONDITION(
        std::ranges::all_of(atomSet.foci,
                            [&](auto atomIdx) { return atomIdx < getNumAtoms(); }),
        "configuration focus index is out of range")
    PRECONDITION(
        std::ranges::all_of(atomSet.atoms,
                            [&](auto atomIdx) { return atomIdx < getNumAtoms(); }),
        "configuration atom index is out of range")
    std::ranges::sort(atomSet.foci);
    atomSet.foci.erase(std::unique(atomSet.foci.begin(), atomSet.foci.end()),
                       atomSet.foci.end());
    std::ranges::sort(atomSet.atoms);
    atomSet.atoms.erase(std::unique(atomSet.atoms.begin(), atomSet.atoms.end()),
                        atomSet.atoms.end());
    for (const auto focusIdx : atomSet.foci) {
      PRECONDITION(std::ranges::binary_search(atomSet.atoms, focusIdx),
                   "configuration focus is absent from its local atoms")
      recordedFoci.set(focusIdx);
    }
  }
  PRECONDITION(recordedFoci == foci,
               "configuration metadata does not match the focus bitset")
  d_configuration_foci = std::move(foci);
  d_configuration_atom_sets = std::move(configurationAtomSets);
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

bool CIPMol::hasUniqueBond(Atom *begin, Atom *end) const {
  unsigned int count = 0;
  for (const auto bond : getBonds(begin)) {
    if (bond->getOtherAtom(begin) == end && ++count > 1u) {
      return false;
    }
  }
  return count == 1u;
}

void CIPMol::initConstitutionalAutomorphismData() const {
  if (d_constitutional_automorphism_data_initialized) {
    return;
  }

  const auto numAtoms = getNumAtoms();
  constexpr auto unassigned = std::numeric_limits<unsigned int>::max();
  d_component_ids.assign(numAtoms, unassigned);
  d_component_local_indices.assign(numAtoms, unassigned);
  d_constitutional_symmetry_classes.assign(numAtoms, 0u);
  if (numAtoms == 0u) {
    d_constitutional_automorphism_data_initialized = true;
    return;
  }

  d_component_atoms.clear();
  std::vector<Atom *> queue;
  for (unsigned int atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
    if (d_component_ids[atomIdx] != unassigned) {
      continue;
    }
    const auto component = static_cast<unsigned int>(d_component_atoms.size());
    d_component_atoms.emplace_back();
    queue.clear();
    queue.push_back(getAtom(atomIdx));
    d_component_ids[atomIdx] = component;
    for (std::size_t pos = 0; pos < queue.size(); ++pos) {
      const auto atom = queue[pos];
      d_component_local_indices[atom->getIdx()] =
          static_cast<unsigned int>(d_component_atoms.back().size());
      d_component_atoms.back().push_back(atom->getIdx());
      for (const auto neighbor : getNeighbors(atom)) {
        if (d_component_ids[neighbor->getIdx()] == unassigned) {
          d_component_ids[neighbor->getIdx()] = component;
          queue.push_back(neighbor);
        }
      }
    }
  }

  d_component_has_cycle.assign(d_component_atoms.size(), 0u);
  d_component_is_supported.assign(d_component_atoms.size(), 1u);
  d_component_automorphism_search_callbacks.assign(d_component_atoms.size(),
                                                   0u);
  d_component_automorphism_prefilter_queries.assign(d_component_atoms.size(),
                                                    0u);
  d_component_mols.clear();
  d_component_mols.resize(d_component_atoms.size());
  d_component_bond_indices.clear();
  d_component_bond_indices.resize(d_component_atoms.size());
  d_constitutional_symmetry_classes_initialized.assign(d_component_atoms.size(),
                                                       0u);
  std::unordered_set<std::uint64_t> atomPairs;
  atomPairs.reserve(d_bonds.size());
  for (const auto bond : d_bonds) {
    const auto begin = bond->getBeginAtomIdx();
    const auto end = bond->getEndAtomIdx();
    const auto component = d_component_ids[begin];
    d_component_bond_indices[component].push_back(bond->getIdx());

    const auto first = std::min(begin, end);
    const auto second = std::max(begin, end);
    const auto atomPair = (static_cast<std::uint64_t>(first) << 32u) | second;
    if (begin == end || !atomPairs.insert(atomPair).second ||
        !hasSupportedConstitutionalBondType(*bond)) {
      // VF2 treats the molecular topology as an undirected simple graph. Fall
      // back to normal CIP traversal when that graph cannot exactly represent
      // a bond, or when ordinary CIP expansion would reject its bond order.
      d_component_is_supported[component] = 0u;
    }
  }
  for (std::size_t component = 0; component < d_component_atoms.size();
       ++component) {
    d_component_has_cycle[component] =
        d_component_bond_indices[component].size() >=
        d_component_atoms[component].size();
  }

  // Peel leaves to obtain the cyclic 2-core. A symmetry query is useful only
  // when both compared branches lead into this core; this keeps remote rings
  // from making unrelated pendant trees eligible for VF2.
  std::vector<unsigned int> coreDegree(numAtoms, 0u);
  std::vector<unsigned int> peelQueue;
  peelQueue.reserve(numAtoms);
  for (const auto atom : atoms()) {
    coreDegree[atom->getIdx()] = atom->getDegree();
    if (coreDegree[atom->getIdx()] < 2u) {
      peelQueue.push_back(atom->getIdx());
    }
  }
  d_atom_in_cyclic_core.assign(numAtoms, 1u);
  for (std::size_t pos = 0; pos < peelQueue.size(); ++pos) {
    const auto atomIdx = peelQueue[pos];
    if (!d_atom_in_cyclic_core[atomIdx]) {
      continue;
    }
    d_atom_in_cyclic_core[atomIdx] = 0u;
    for (const auto neighbor : getNeighbors(getAtom(atomIdx))) {
      const auto neighborIdx = neighbor->getIdx();
      if (d_atom_in_cyclic_core[neighborIdx] && coreDegree[neighborIdx] > 0u &&
          --coreDegree[neighborIdx] == 1u) {
        peelQueue.push_back(neighborIdx);
      }
    }
  }
  // The vertices outside a graph's 2-core form trees rooted at that core.
  // Multi-source distances therefore answer every directed "does this side
  // reach the core?" query in O(1), instead of running a BFS for each bond in
  // a long pendant tree.
  d_cyclic_core_distances.assign(numAtoms, unassigned);
  std::vector<unsigned int> coreQueue;
  coreQueue.reserve(numAtoms);
  for (unsigned int atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
    if (d_atom_in_cyclic_core[atomIdx]) {
      d_cyclic_core_distances[atomIdx] = 0u;
      coreQueue.push_back(atomIdx);
    }
  }
  for (std::size_t pos = 0; pos < coreQueue.size(); ++pos) {
    if ((pos & 0xffu) == 0u && ControlCHandler::getGotSignal()) {
      throw ControlCCaught();
    }
    const auto atomIdx = coreQueue[pos];
    for (const auto neighbor : getNeighbors(getAtom(atomIdx))) {
      const auto neighborIdx = neighbor->getIdx();
      if (d_cyclic_core_distances[neighborIdx] == unassigned) {
        d_cyclic_core_distances[neighborIdx] =
            d_cyclic_core_distances[atomIdx] + 1u;
        coreQueue.push_back(neighborIdx);
      }
    }
  }

  d_constitutional_automorphism_data_initialized = true;
}

void CIPMol::initConstitutionalSymmetryClasses(unsigned int component) const {
  if (d_constitutional_symmetry_classes_initialized[component]) {
    return;
  }

  // Build an equitable constitutional partition directly from the fields
  // used by the exact matcher. Unlike canonical SMILES ranking, this ignores
  // annotations such as dummy-atom maps and avoids special ring traversals on
  // the highly fused systems for which this shortcut is intended. The result
  // is only a necessary-condition filter; VF2 still proves equality exactly.
  const auto &componentAtoms = d_component_atoms[component];
  std::vector<std::pair<AtomPropertySignature, unsigned int>> initialOrder;
  initialOrder.reserve(componentAtoms.size());
  for (std::size_t localIdx = 0; localIdx < componentAtoms.size(); ++localIdx) {
    initialOrder.emplace_back(
        getAtomPropertySignature(*getAtom(componentAtoms[localIdx])), localIdx);
  }
  std::ranges::sort(initialOrder, [](const auto &first, const auto &second) {
    return first.first < second.first;
  });

  std::vector<unsigned int> colors(componentAtoms.size());
  unsigned int classCount = 0u;
  for (std::size_t pos = 0; pos < initialOrder.size(); ++pos) {
    if (pos != 0u && initialOrder[pos - 1u].first != initialOrder[pos].first) {
      ++classCount;
    }
    colors[initialOrder[pos].second] = classCount;
  }
  ++classCount;

  // Use only logarithmically many synchronous refinement rounds so path-like
  // appendages cannot turn this prefilter into an unbounded propagation pass.
  // Regular symmetric cages normally stabilize after one round. Stopping
  // early is conservative: it can only merge candidate classes and cannot
  // create a false equality.
  const auto maxRounds =
      std::max<std::size_t>(1u, 2u * std::bit_width(componentAtoms.size()));
  std::vector<RefinementSignature> signatures(componentAtoms.size());
  std::vector<unsigned int> order(componentAtoms.size());
  std::vector<unsigned int> nextColors(componentAtoms.size());
  for (std::size_t round = 0; round < maxRounds; ++round) {
    for (std::size_t localIdx = 0; localIdx < componentAtoms.size();
         ++localIdx) {
      if ((localIdx & 0xffu) == 0u && ControlCHandler::getGotSignal()) {
        throw ControlCCaught();
      }
      auto &signature = signatures[localIdx];
      signature.oldColor = colors[localIdx];
      signature.neighbors.clear();
      const auto atom = getAtom(componentAtoms[localIdx]);
      signature.neighbors.reserve(atom->getDegree());
      for (const auto bond : getBonds(atom)) {
        const auto neighborIdx = bond->getOtherAtomIdx(atom->getIdx());
        signature.neighbors.emplace_back(
            static_cast<int>(bond->getBondType()), bond->getIsAromatic(),
            colors[d_component_local_indices[neighborIdx]]);
      }
      std::ranges::sort(signature.neighbors);
      order[localIdx] = localIdx;
    }
    std::ranges::sort(order, [&](const auto first, const auto second) {
      return signatures[first] < signatures[second];
    });

    unsigned int nextClassCount = 0u;
    for (std::size_t pos = 0; pos < order.size(); ++pos) {
      if (pos != 0u && signatures[order[pos - 1u]] != signatures[order[pos]]) {
        ++nextClassCount;
      }
      nextColors[order[pos]] = nextClassCount;
    }
    ++nextClassCount;
    colors.swap(nextColors);
    if (nextClassCount == classCount) {
      break;
    }
    classCount = nextClassCount;
  }

  for (std::size_t localIdx = 0; localIdx < componentAtoms.size(); ++localIdx) {
    d_constitutional_symmetry_classes[componentAtoms[localIdx]] =
        colors[localIdx];
  }
  d_constitutional_symmetry_classes_initialized[component] = 1u;
}

const RWMol &CIPMol::getConstitutionalComponentMol(
    unsigned int component) const {
  if (d_component_mols[component]) {
    return *d_component_mols[component];
  }

  auto componentMol = std::make_shared<RWMol>();
  const auto &componentAtoms = d_component_atoms[component];
  for (const auto atomIdx : componentAtoms) {
    auto atomCopy = std::unique_ptr<Atom>(getAtom(atomIdx)->copy());
    componentMol->addAtom(atomCopy.get(), false, true);
    atomCopy.release();
  }

  const auto &componentBondIndices = d_component_bond_indices[component];
  for (const auto bondIdx : componentBondIndices) {
    const auto bond = getBond(bondIdx);
    const auto beginIdx = bond->getBeginAtomIdx();
    auto bondCopy = std::unique_ptr<Bond>(bond->copy());
    bondCopy->setBeginAtomIdx(d_component_local_indices[beginIdx]);
    bondCopy->setEndAtomIdx(d_component_local_indices[bond->getEndAtomIdx()]);
    // Chirality is deliberately excluded from this constitutional equality
    // proof. Clearing index-bearing stereo metadata also keeps the component
    // copy independent of the atom numbering in the source molecule.
    bondCopy->setStereo(Bond::STEREONONE);
    bondCopy->getStereoAtoms().clear();
    bondCopy->setBondDir(Bond::NONE);
    componentMol->addBond(bondCopy.get(), true);
    bondCopy.release();
  }

  d_component_mols[component] = std::move(componentMol);
  return *d_component_mols[component];
}

bool CIPMol::branchReachesCyclicCore(Atom *root, Atom *end) const {
  const auto bond = d_mol.getBondBetweenAtoms(root->getIdx(), end->getIdx());
  if (bond == nullptr) {
    return false;
  }
  constexpr auto unreachable = std::numeric_limits<unsigned int>::max();
  const auto rootDistance = d_cyclic_core_distances[root->getIdx()];
  const auto endDistance = d_cyclic_core_distances[end->getIdx()];
  return endDistance != unreachable &&
         (endDistance == 0u || endDistance < rootDistance);
}

bool CIPMol::hasConstitutionalAutomorphism(Atom *root, Atom *from,
                                           Atom *to) const {
  std::vector<unsigned int> movedAtoms;
  return hasConstitutionalAutomorphism(root, from, to, movedAtoms);
}

bool CIPMol::hasConstitutionalAutomorphism(
    Atom *root, Atom *from, Atom *to,
    std::vector<unsigned int> &movedAtoms) const {
  movedAtoms.clear();
  if (root == nullptr || from == nullptr || to == nullptr) {
    return false;
  }
  const auto numAtoms = getNumAtoms();
  const auto rootIdx = root->getIdx();
  const auto fromIdx = from->getIdx();
  const auto toIdx = to->getIdx();
  if (rootIdx >= numAtoms || fromIdx >= numAtoms || toIdx >= numAtoms ||
      getAtom(rootIdx) != root || getAtom(fromIdx) != from ||
      getAtom(toIdx) != to) {
    return false;
  }
  if (!hasUniqueBond(root, from) || !hasUniqueBond(root, to)) {
    return false;
  }
  if (from == to) {
    return true;
  }

  return hasConstitutionalAutomorphism(rootIdx, fromIdx, toIdx, {}, {},
                                       &movedAtoms);
}

bool CIPMol::hasConstitutionalAutomorphism(
    Atom *root, Atom *from, Atom *to, std::span<const std::uint64_t> fixedAtoms,
    std::span<const unsigned int> addedFixedAtoms) const {
  std::vector<unsigned int> movedAtoms;
  return hasConstitutionalAutomorphism(root, from, to, fixedAtoms,
                                       addedFixedAtoms, movedAtoms);
}

bool CIPMol::hasConstitutionalAutomorphism(
    Atom *root, Atom *from, Atom *to, std::span<const std::uint64_t> fixedAtoms,
    std::span<const unsigned int> addedFixedAtoms,
    std::vector<unsigned int> &movedAtoms) const {
  movedAtoms.clear();
  if (root == nullptr || from == nullptr || to == nullptr) {
    return false;
  }
  const auto numAtoms = getNumAtoms();
  const auto rootIdx = root->getIdx();
  const auto fromIdx = from->getIdx();
  const auto toIdx = to->getIdx();
  if (rootIdx >= numAtoms || fromIdx >= numAtoms || toIdx >= numAtoms ||
      getAtom(rootIdx) != root || getAtom(fromIdx) != from ||
      getAtom(toIdx) != to || fixedAtoms.size() < (numAtoms + 63u) / 64u) {
    return false;
  }
  if (!hasUniqueBond(root, from) || !hasUniqueBond(root, to)) {
    return false;
  }
  if (from == to) {
    return true;
  }

  return hasConstitutionalAutomorphism(rootIdx, fromIdx, toIdx, fixedAtoms,
                                       addedFixedAtoms, &movedAtoms);
}

bool CIPMol::constitutionalAutomorphismMovesConfiguration(
    std::span<const unsigned int> movedAtoms,
    const Atom *configurationOwner) const {
  if (movedAtoms.empty()) {
    return false;
  }
  const auto ownerConfiguration =
      findConfigurationOwnedBy(configurationOwner);
  if (!ownerConfiguration) {
    return true;
  }

  for (std::size_t i = 0; i < d_configuration_atom_sets.size(); ++i) {
    if (i == *ownerConfiguration) {
      continue;
    }
    for (const auto atomIdx : d_configuration_atom_sets[i].atoms) {
      if (std::ranges::find(movedAtoms, atomIdx) != movedAtoms.end()) {
        return true;
      }
    }
  }
  return false;
}

bool CIPMol::hasConfigurationPreservingConstitutionalAutomorphism(
    Atom *root, Atom *from, Atom *to,
    std::span<const std::uint64_t> fixedAtoms,
    std::span<const unsigned int> addedFixedAtoms,
    const Atom *configurationOwner,
    std::vector<unsigned int> &movedAtoms) const {
  movedAtoms.clear();
  if (root == nullptr || from == nullptr || to == nullptr ||
      configurationOwner == nullptr) {
    return false;
  }
  const auto numAtoms = getNumAtoms();
  const auto rootIdx = root->getIdx();
  const auto fromIdx = from->getIdx();
  const auto toIdx = to->getIdx();
  const auto ownerIdx = configurationOwner->getIdx();
  if (rootIdx >= numAtoms || fromIdx >= numAtoms || toIdx >= numAtoms ||
      ownerIdx >= numAtoms || getAtom(rootIdx) != root ||
      getAtom(fromIdx) != from || getAtom(toIdx) != to ||
      getAtom(ownerIdx) != configurationOwner ||
      (!fixedAtoms.empty() &&
       fixedAtoms.size() < (numAtoms + 63u) / 64u) ||
      !hasUniqueBond(root, from) || !hasUniqueBond(root, to)) {
    return false;
  }

  const auto ownerConfiguration =
      findConfigurationOwnedBy(configurationOwner);
  if (!ownerConfiguration) {
    return false;
  }
  if (from == to) {
    return true;
  }

  const auto wordCount = (numAtoms + 63u) / 64u;
  std::vector<std::uint64_t> protectedAtoms(wordCount, 0u);
  if (!fixedAtoms.empty()) {
    std::copy_n(fixedAtoms.begin(), wordCount, protectedAtoms.begin());
  }
  const auto protectAtom = [&](unsigned int atomIdx) {
    protectedAtoms[atomIdx / 64u] |=
        std::uint64_t{1} << (atomIdx % 64u);
  };
  for (const auto atomIdx : addedFixedAtoms) {
    if (atomIdx < numAtoms) {
      protectAtom(atomIdx);
    }
  }
  for (std::size_t i = 0; i < d_configuration_atom_sets.size(); ++i) {
    if (i != *ownerConfiguration) {
      for (const auto atomIdx : d_configuration_atom_sets[i].atoms) {
        protectAtom(atomIdx);
      }
    }
  }

  return hasConstitutionalAutomorphism(rootIdx, fromIdx, toIdx, protectedAtoms,
                                       {}, &movedAtoms);
}

std::optional<std::size_t> CIPMol::findConfigurationOwnedBy(
    const Atom *configurationOwner) const {
  if (configurationOwner == nullptr ||
      configurationOwner->getIdx() >= getNumAtoms() ||
      getAtom(configurationOwner->getIdx()) != configurationOwner) {
    return std::nullopt;
  }

  const auto ownerIdx = configurationOwner->getIdx();
  std::optional<std::size_t> result;
  for (std::size_t i = 0; i < d_configuration_atom_sets.size(); ++i) {
    if (!std::ranges::binary_search(d_configuration_atom_sets[i].foci,
                                    ownerIdx)) {
      continue;
    }
    if (result) {
      // Shared foci need a stable configuration identity, not an atom-only
      // owner. Until such an identity is available, remain conservative.
      return std::nullopt;
    }
    result = i;
  }
  return result;
}

bool CIPMol::isInSameComponent(Atom *first, Atom *second) const {
  if (first == nullptr || second == nullptr) {
    return false;
  }
  const auto firstIdx = first->getIdx();
  const auto secondIdx = second->getIdx();
  if (firstIdx >= getNumAtoms() || secondIdx >= getNumAtoms() ||
      getAtom(firstIdx) != first || getAtom(secondIdx) != second) {
    return false;
  }
  initConstitutionalAutomorphismData();
  return d_component_ids[firstIdx] == d_component_ids[secondIdx];
}

bool CIPMol::atomsConstitutionallyEquivalent(const Atom &first,
                                             const Atom &second) const {
  // Every refinement class is a subset of one complete atom-property class,
  // so this single comparison covers the explicit constitutional fields as
  // well as the refined neighborhood invariant.
  return d_constitutional_symmetry_classes[first.getIdx()] ==
         d_constitutional_symmetry_classes[second.getIdx()];
}

bool CIPMol::bondsConstitutionallyEquivalent(const Bond &first,
                                             const Bond &second) const {
  return first.getBondType() == second.getBondType() &&
         first.getIsAromatic() == second.getIsAromatic() &&
         d_kekulized_bonds[first.getIdx()] ==
             d_kekulized_bonds[second.getIdx()];
}

CIPMol::ConstitutionalAutomorphismEvidence &CIPMol::getAutomorphismEvidence(
    const ConstitutionalAutomorphismKey &key) const {
  if (const auto current =
          d_current_constitutional_automorphism_evidence.find(key);
      current != d_current_constitutional_automorphism_evidence.end()) {
    return current->second;
  }

  std::optional<ConstitutionalAutomorphismEvidence> previousEvidence;
  if (const auto previous =
          d_previous_constitutional_automorphism_evidence.find(key);
      previous != d_previous_constitutional_automorphism_evidence.end()) {
    previousEvidence.emplace(std::move(previous->second));
    d_previous_automorphism_evidence_atom_indices -=
        previousEvidence->storedAtomIndexCount;
    d_previous_constitutional_automorphism_evidence.erase(previous);
  }

  if (d_current_constitutional_automorphism_evidence.size() >=
          AUTOMORPHISM_EVIDENCE_KEYS_PER_GENERATION ||
      (previousEvidence &&
       d_current_automorphism_evidence_atom_indices +
               previousEvidence->storedAtomIndexCount >
           AUTOMORPHISM_EVIDENCE_ATOM_INDICES_PER_GENERATION)) {
    d_previous_constitutional_automorphism_evidence =
        std::move(d_current_constitutional_automorphism_evidence);
    d_previous_automorphism_evidence_atom_indices =
        d_current_automorphism_evidence_atom_indices;
    d_current_constitutional_automorphism_evidence.clear();
    d_current_automorphism_evidence_atom_indices = 0;
  }

  const auto insertion =
      previousEvidence
          ? d_current_constitutional_automorphism_evidence.emplace(
                key, std::move(*previousEvidence))
          : d_current_constitutional_automorphism_evidence.try_emplace(key);
  if (previousEvidence) {
    d_current_automorphism_evidence_atom_indices +=
        insertion.first->second.storedAtomIndexCount;
  }
  return insertion.first->second;
}

std::optional<bool> CIPMol::findAutomorphismEvidence(
    const ConstitutionalAutomorphismEvidence &evidence,
    std::span<const std::uint64_t> fixedAtoms,
    std::span<const unsigned int> addedFixedAtoms,
    std::vector<unsigned int> *movedAtomsOut) {
  for (const auto &candidate : evidence.movedAtomSets) {
    if (atomSetIsDisjointFromMask(candidate, fixedAtoms, addedFixedAtoms)) {
      if (movedAtomsOut != nullptr) {
        *movedAtomsOut = candidate;
      }
      return true;
    }
  }
  if (evidence.searchDisabled) {
    return false;
  }
  if (std::ranges::any_of(evidence.failedFixedAtomSets,
                          [&](const auto &failedFixedAtoms) {
                            return atomSetIsSubsetOfMask(
                                failedFixedAtoms, fixedAtoms, addedFixedAtoms);
                          })) {
    return false;
  }
  return std::nullopt;
}

void CIPMol::addAutomorphismWitness(
    ConstitutionalAutomorphismEvidence &evidence,
    std::vector<unsigned int> movedAtoms) const {
  std::ranges::sort(movedAtoms);
  if (std::ranges::any_of(evidence.movedAtomSets, [&](const auto &existing) {
        return atomSetIsSubset(existing, movedAtoms);
      })) {
    return;
  }
  std::size_t removedIndices = 0;
  std::size_t removedSets = 0;
  for (const auto &existing : evidence.movedAtomSets) {
    if (atomSetIsSubset(movedAtoms, existing)) {
      removedIndices += existing.size();
      ++removedSets;
    }
  }
  if (evidence.movedAtomSets.size() - removedSets >=
          MAX_AUTOMORPHISM_EVIDENCE_PER_KIND ||
      d_current_automorphism_evidence_atom_indices - removedIndices +
              movedAtoms.size() >
          AUTOMORPHISM_EVIDENCE_ATOM_INDICES_PER_GENERATION) {
    return;
  }
  std::erase_if(evidence.movedAtomSets, [&](const auto &existing) {
    return atomSetIsSubset(movedAtoms, existing);
  });
  evidence.storedAtomIndexCount -= removedIndices;
  d_current_automorphism_evidence_atom_indices -= removedIndices;
  evidence.storedAtomIndexCount += movedAtoms.size();
  d_current_automorphism_evidence_atom_indices += movedAtoms.size();
  evidence.movedAtomSets.push_back(std::move(movedAtoms));
}

void CIPMol::addAutomorphismFailure(
    ConstitutionalAutomorphismEvidence &evidence,
    std::vector<unsigned int> fixedAtoms) const {
  std::ranges::sort(fixedAtoms);
  if (std::ranges::any_of(evidence.failedFixedAtomSets,
                          [&](const auto &existing) {
                            return atomSetIsSubset(existing, fixedAtoms);
                          })) {
    return;
  }
  std::size_t removedIndices = 0;
  std::size_t removedSets = 0;
  for (const auto &existing : evidence.failedFixedAtomSets) {
    if (atomSetIsSubset(fixedAtoms, existing)) {
      removedIndices += existing.size();
      ++removedSets;
    }
  }
  if (evidence.failedFixedAtomSets.size() - removedSets >=
          MAX_AUTOMORPHISM_EVIDENCE_PER_KIND ||
      d_current_automorphism_evidence_atom_indices - removedIndices +
              fixedAtoms.size() >
          AUTOMORPHISM_EVIDENCE_ATOM_INDICES_PER_GENERATION) {
    return;
  }
  std::erase_if(evidence.failedFixedAtomSets, [&](const auto &existing) {
    return atomSetIsSubset(fixedAtoms, existing);
  });
  evidence.storedAtomIndexCount -= removedIndices;
  d_current_automorphism_evidence_atom_indices -= removedIndices;
  evidence.storedAtomIndexCount += fixedAtoms.size();
  d_current_automorphism_evidence_atom_indices += fixedAtoms.size();
  evidence.failedFixedAtomSets.push_back(std::move(fixedAtoms));
}

bool CIPMol::hasConstitutionalAutomorphism(
    unsigned int rootIdx, unsigned int fromIdx, unsigned int toIdx,
    std::span<const std::uint64_t> fixedAtoms,
    std::span<const unsigned int> addedFixedAtoms,
    std::vector<unsigned int> *movedAtomsOut) const {
  initConstitutionalAutomorphismData();
  const auto component = d_component_ids[rootIdx];
  if (d_component_ids[fromIdx] != component ||
      d_component_ids[toIdx] != component ||
      !d_component_has_cycle[component] ||
      !d_component_is_supported[component] ||
      !branchReachesCyclicCore(getAtom(rootIdx), getAtom(fromIdx)) ||
      !branchReachesCyclicCore(getAtom(rootIdx), getAtom(toIdx))) {
    return false;
  }

  if (fromIdx > toIdx) {
    std::swap(fromIdx, toIdx);
  }
  if (!atomsHaveSameConstitutionalProperties(*getAtom(fromIdx),
                                             *getAtom(toIdx))) {
    return false;
  }
  initConstitutionalSymmetryClasses(component);
  if (!atomsConstitutionallyEquivalent(*getAtom(fromIdx), *getAtom(toIdx))) {
    return false;
  }

  const ConstitutionalAutomorphismKey key{rootIdx, fromIdx, toIdx};
  auto &evidence = getAutomorphismEvidence(key);
  if (const auto cached = findAutomorphismEvidence(
          evidence, fixedAtoms, addedFixedAtoms, movedAtomsOut)) {
    return *cached;
  }

  const auto &componentAtoms = d_component_atoms[component];
  const auto cacheFailure = [&](bool pathIndependent = false) {
    std::vector<unsigned int> fixedAtomSet;
    if (!pathIndependent) {
      for (const auto atomIdx : componentAtoms) {
        if (atomIdx != rootIdx &&
            fixedSetContainsAtom(fixedAtoms, addedFixedAtoms, atomIdx)) {
          fixedAtomSet.push_back(atomIdx);
        }
      }
    }
    addAutomorphismFailure(evidence, std::move(fixedAtomSet));
  };
  if (fixedSetContainsAtom(fixedAtoms, addedFixedAtoms, fromIdx) ||
      fixedSetContainsAtom(fixedAtoms, addedFixedAtoms, toIdx)) {
    cacheFailure();
    return false;
  }

  const auto scaledSearchCallbacks =
      componentAtoms.size() > MAX_AUTOMORPHISM_SEARCH_CALLBACKS /
                                  AUTOMORPHISM_SEARCH_CALLBACKS_PER_ATOM
          ? MAX_AUTOMORPHISM_SEARCH_CALLBACKS
          : componentAtoms.size() * AUTOMORPHISM_SEARCH_CALLBACKS_PER_ATOM;
  const auto searchCallbackLimit =
      std::clamp(scaledSearchCallbacks, MIN_AUTOMORPHISM_SEARCH_CALLBACKS,
                 MAX_AUTOMORPHISM_SEARCH_CALLBACKS);
  const auto componentSearchCallbackLimit =
      COMPONENT_AUTOMORPHISM_SEARCH_MULTIPLIER * searchCallbackLimit;
  auto &componentSearchCallbacks =
      d_component_automorphism_search_callbacks[component];
  if (componentSearchCallbacks >= componentSearchCallbackLimit) {
    return false;
  }

  // Exact graph-distance individualization is a strong VF2 prefilter, but it
  // is itself linear in the component. Use it for only logarithmically many
  // uncached endpoint pairs so a succession of cheap negative proofs cannot
  // become quadratic. Later searches remain exact and are bounded by their
  // root/endpoint/fixed constraints plus the callback allowances below.
  auto &prefilterQueries =
      d_component_automorphism_prefilter_queries[component];
  const auto prefilterQueryLimit =
      std::max<std::size_t>(1u, 2u * std::bit_width(componentAtoms.size()));
  const bool useDistanceInvariants = prefilterQueries < prefilterQueryLimit;
  std::vector<unsigned int> rootDistances;
  std::vector<unsigned int> fromDistances;
  std::vector<unsigned int> toDistances;
  if (useDistanceInvariants) {
    ++prefilterQueries;
    std::vector<unsigned int> queue;
    queue.reserve(componentAtoms.size());
    const auto getDistances = [&](unsigned int start) {
      constexpr auto unreachable = std::numeric_limits<unsigned int>::max();
      std::vector<unsigned int> distances(componentAtoms.size(), unreachable);
      queue.clear();
      queue.push_back(start);
      distances[d_component_local_indices[start]] = 0u;
      for (std::size_t pos = 0; pos < queue.size(); ++pos) {
        if ((pos & 0xffu) == 0u && ControlCHandler::getGotSignal()) {
          throw ControlCCaught();
        }
        const auto atomIdx = queue[pos];
        const auto atomDistance = distances[d_component_local_indices[atomIdx]];
        for (const auto neighbor : getNeighbors(getAtom(atomIdx))) {
          const auto neighborIdx = neighbor->getIdx();
          const auto neighborLocal = d_component_local_indices[neighborIdx];
          if (distances[neighborLocal] == unreachable) {
            distances[neighborLocal] = atomDistance + 1u;
            queue.push_back(neighborIdx);
          }
        }
      }
      return distances;
    };
    rootDistances = getDistances(rootIdx);
    fromDistances = getDistances(fromIdx);
    toDistances = getDistances(toIdx);

    // Compare the two exact distance-signature multisets in deterministic
    // linear time using the component size as the key range. A mismatch proves
    // that no root-fixing automorphism can map from to to.
    std::vector<DistanceSignature> querySignatures;
    std::vector<DistanceSignature> targetSignatures;
    querySignatures.reserve(componentAtoms.size());
    targetSignatures.reserve(componentAtoms.size());
    for (const auto atomIdx : componentAtoms) {
      const auto localIdx = d_component_local_indices[atomIdx];
      querySignatures.push_back({d_constitutional_symmetry_classes[atomIdx],
                                 rootDistances[localIdx],
                                 fromDistances[localIdx]});
      targetSignatures.push_back({d_constitutional_symmetry_classes[atomIdx],
                                  rootDistances[localIdx],
                                  toDistances[localIdx]});
    }
    radixSortDistanceSignatures(querySignatures, componentAtoms.size());
    radixSortDistanceSignatures(targetSignatures, componentAtoms.size());
    if (querySignatures != targetSignatures) {
      cacheFailure(true);
      return false;
    }
  }

  std::size_t searchCallbacks = 0;
  const auto consumeSearchCallback = [&]() {
    if (componentSearchCallbacks >= componentSearchCallbackLimit) {
      throw ComponentAutomorphismSearchCallbackLimitReached();
    }
    if (searchCallbacks >= searchCallbackLimit) {
      throw AutomorphismSearchCallbackLimitReached();
    }
    ++componentSearchCallbacks;
    ++searchCallbacks;
  };
  initKekulizedBonds();
  const auto &componentMol = getConstitutionalComponentMol(component);

  // One constrained match answers the exact question directly. Unlike a
  // bounded sample of the automorphism group, this remains effective for very
  // large symmetry groups and for arbitrary molecule sizes.
  SubstructMatchParameters params;
  params.recursionPossible = false;
  params.uniquify = false;
  params.maxMatches = 1;
  params.numThreads = 1;
  params.useChirality = false;
  params.extraAtomCheckOverridesDefaultCheck = true;
  params.extraAtomCheck =
      [this, rootIdx, fromIdx, toIdx, component, fixedAtoms, addedFixedAtoms,
       &rootDistances, &fromDistances, &consumeSearchCallback, &toDistances,
       useDistanceInvariants](const Atom &queryAtom, const Atom &molAtom) {
        if (ControlCHandler::getGotSignal()) {
          throw ControlCCaught();
        }
        const auto queryLocal = queryAtom.getIdx();
        const auto molLocal = molAtom.getIdx();
        const auto &componentAtoms = d_component_atoms[component];
        if (queryLocal >= componentAtoms.size() ||
            molLocal >= componentAtoms.size()) {
          return false;
        }
        consumeSearchCallback();
        const auto queryIdx = componentAtoms[queryLocal];
        const auto molIdx = componentAtoms[molLocal];
        const auto queryFixed =
            queryIdx == rootIdx ||
            fixedSetContainsAtom(fixedAtoms, addedFixedAtoms, queryIdx);
        const auto molFixed =
            molIdx == rootIdx ||
            fixedSetContainsAtom(fixedAtoms, addedFixedAtoms, molIdx);
        if ((queryIdx == fromIdx && molIdx != toIdx) ||
            (useDistanceInvariants &&
             (rootDistances[queryLocal] != rootDistances[molLocal] ||
              fromDistances[queryLocal] != toDistances[molLocal])) ||
            ((queryFixed || molFixed) && molIdx != queryIdx)) {
          return false;
        }
        return atomsConstitutionallyEquivalent(*getAtom(queryIdx),
                                               *getAtom(molIdx));
      };
  params.extraBondCheckOverridesDefaultCheck = true;
  params.extraBondCheck = [this, component, &consumeSearchCallback](
                              const Bond &queryBond, const Bond &molBond) {
    if (ControlCHandler::getGotSignal()) {
      throw ControlCCaught();
    }
    consumeSearchCallback();
    const auto &componentBonds = d_component_bond_indices[component];
    const auto queryIdx = queryBond.getIdx();
    const auto molIdx = molBond.getIdx();
    if (queryIdx >= componentBonds.size() || molIdx >= componentBonds.size()) {
      return false;
    }
    return bondsConstitutionallyEquivalent(*getBond(componentBonds[queryIdx]),
                                           *getBond(componentBonds[molIdx]));
  };

  std::vector<MatchVectType> matches;
  try {
    matches = SubstructMatch(componentMol, componentMol, params);
  } catch (const ComponentAutomorphismSearchCallbackLimitReached &) {
    // Aggregate exhaustion is unknown, never an inequality proof. Existing
    // witnesses remain reusable, but new searches in this component fall back
    // to ordinary CIP traversal.
    evidence.searchDisabled = true;
    return false;
  } catch (const AutomorphismSearchCallbackLimitReached &) {
    // Per-key exhaustion is likewise unknown. Avoid paying the same bounded
    // search again for this endpoint pair without admitting negative evidence.
    evidence.searchDisabled = true;
    return false;
  }
  if (ControlCHandler::getGotSignal()) {
    throw ControlCCaught();
  }
  if (matches.empty()) {
    cacheFailure();
    return false;
  }

  if (matches.front().size() != componentAtoms.size()) {
    return false;
  }
  std::vector<unsigned int> componentMapping(
      componentAtoms.size(), std::numeric_limits<unsigned int>::max());
  for (const auto &[queryIdx, molIdx] : matches.front()) {
    if (queryIdx < 0 || molIdx < 0 ||
        static_cast<unsigned int>(queryIdx) >= componentAtoms.size() ||
        static_cast<unsigned int>(molIdx) >= componentAtoms.size()) {
      return false;
    }
    const auto queryLocal = static_cast<unsigned int>(queryIdx);
    const auto molLocal = static_cast<unsigned int>(molIdx);
    componentMapping[queryLocal] = componentAtoms[molLocal];
  }
  if (std::ranges::any_of(componentMapping, [](const auto atomIdx) {
        return atomIdx == std::numeric_limits<unsigned int>::max();
      })) {
    return false;
  }

  std::vector<unsigned int> witnessMovedAtoms;
  for (const auto atomIdx : componentAtoms) {
    if (componentMapping[d_component_local_indices[atomIdx]] != atomIdx) {
      witnessMovedAtoms.push_back(atomIdx);
    }
  }
  if (movedAtomsOut != nullptr) {
    *movedAtomsOut = witnessMovedAtoms;
  }
  addAutomorphismWitness(evidence, std::move(witnessMovedAtoms));
  return true;
}

void CIPMol::initKekulizedBonds() const {
  if (!d_kekulized_bonds.empty()) {
    return;
  }

  std::vector<RDKit::Bond::BondType> bonds;
  bonds.reserve(d_mol.getNumBonds());
  const bool hasAromaticBond =
      std::ranges::any_of(d_bonds, [](const Bond *candidate) {
        return candidate->getBondType() == Bond::AROMATIC;
      });
  if (hasAromaticBond) {
    RWMol tmp{d_mol};
    const ROMol *bondSource = &tmp;
    try {
      MolOps::Kekulize(tmp);
    } catch (const MolSanitizeException &) {
      // Kekulize() may have changed some bonds before discovering that no
      // valid assignment exists. Fall back to the untouched input instead of
      // caching that partial assignment.
      bondSource = &d_mol;
    }
    for (const auto candidate : bondSource->bonds()) {
      bonds.push_back(candidate->getBondType());
    }
  } else {
    for (const auto candidate : d_bonds) {
      bonds.push_back(candidate->getBondType());
    }
  }
  d_kekulized_bonds = std::move(bonds);
}

int CIPMol::getBondOrder(Bond *bond) const {
  PRECONDITION(bond, "bad bond")
  initKekulizedBonds();

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
