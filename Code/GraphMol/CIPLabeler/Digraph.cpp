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
#include <sstream>

#include "Digraph.h"
#include "CIPMol.h"
#include "Node.h"
#include "Edge.h"

namespace RDKit {
namespace CIPLabeler {

namespace {

/**
 * Upper limit on the size of the digraph, stops out of memory error with a
 * more graceful failure. 0=Infinite
 */
constexpr std::size_t MAX_NODE_COUNT = 100000;

/**
 * Used for debugging only, 0=Infinite
 */
const int MAX_NODE_DIST = 0;
}  // namespace

Node &Digraph::addNode(NodeVisitState &&visit, Atom *atom,
                       boost::rational<int> &&frac, int dist, int flags,
                       const Node *parent) {
  if (MAX_NODE_COUNT > 0 && d_nodes.size() >= MAX_NODE_COUNT) {
    std::stringstream errmsg;
    errmsg << "Digraph generation failed: more than " << MAX_NODE_COUNT
           << "nodes found.";
    throw TooManyNodesException(errmsg.str());
  }
  d_nodes.emplace_back(this, std::move(visit), atom, std::move(frac), dist,
                       flags, parent);
  if (atom == nullptr) {
    d_seen_null = true;
  } else {
    const auto atom_idx = atom->getIdx();
    if (atom_idx < d_seen_atoms.size() && d_mol.getAtom(atom_idx) == atom) {
      d_seen_atoms[atom_idx] = true;
    }
  }
  return d_nodes.back();
}

bool Digraph::seenAtom(Atom *atom) const {
  if (atom == nullptr) {
    return d_seen_null;
  }
  const auto atom_idx = atom->getIdx();
  if (atom_idx < d_seen_atoms.size() && d_mol.getAtom(atom_idx) == atom) {
    return d_seen_atoms[atom_idx];
  }
  // Preserve the old pointer-identity behavior for an unexpected atom that
  // does not belong to this molecule.
  return std::ranges::any_of(
      d_nodes, [&](const auto &n) { return n.getAtom() == atom; });
}

void Digraph::addEdge(Node *beg, Bond *bond, Node *end) {
  d_edges.emplace_back(beg, end, bond);
  auto &e = d_edges.back();
  beg->add(&e);
  end->add(&e);
}

Digraph::Digraph(const CIPMol &mol, Atom *atom, bool atropisomerMode)
    : d_mol{mol}, d_seen_atoms(mol.getNumAtoms()) {
  PRECONDITION(atom, "cannot init digraph on a nullptr")

  auto visit = NodeVisitState((d_mol.getNumAtoms() + 63u) / 64u);
  visit[atom->getIdx() / 64u] |= std::uint64_t{1} << (atom->getIdx() % 64u);

  auto dist = 1;
  auto flags = 0x0;
  auto atomic_num = atom->getAtomicNum();

  dp_root = &addNode(std::move(visit), atom, atomic_num, dist, flags);
  dp_origin = dp_root;
  d_atropisomerMode = atropisomerMode;
}

const CIPMol &Digraph::getMol() const { return d_mol; };

Node *Digraph::getOriginalRoot() const { return dp_origin; };

Node *Digraph::getCurrentRoot() const { return dp_root; }

int Digraph::getNumNodes() const { return d_nodes.size(); }

std::vector<Node *> Digraph::getNodes(Atom *atom) const {
  return collectNodes(atom, false);
}

std::vector<Node *> Digraph::getNodesForAuxiliaryLabeling(Atom *atom) const {
  if (atom == nullptr || getCurrentRoot() != getOriginalRoot()) {
    return collectNodes(atom, true);
  }
  boost::dynamic_bitset<> targets(d_mol.getNumAtoms());
  targets.set(atom->getIdx());
  return getNodesForAuxiliaryLabeling(targets)[atom->getIdx()];
}

std::vector<std::vector<Node *>> Digraph::getNodesForAuxiliaryLabeling(
    const boost::dynamic_bitset<> &targets) const {
  PRECONDITION(targets.size() == d_mol.getNumAtoms(),
               "auxiliary target bitset has the wrong size")

  std::vector<std::vector<Node *>> result(d_mol.getNumAtoms());
  if (targets.none()) {
    return result;
  }

  // Visit bits and original-parent links describe paths from the original
  // root. Fall back conservatively if a caller asks during a temporary reroot.
  if (getCurrentRoot() != getOriginalRoot()) {
    for (auto idx = targets.find_first(); idx != boost::dynamic_bitset<>::npos;
         idx = targets.find_next(idx)) {
      result[idx] = collectNodes(d_mol.getAtom(idx), true);
    }
    return result;
  }

  std::vector<Node *> nodes{getOriginalRoot()};
  std::vector<unsigned int> reachabilityQueue;
  reachabilityQueue.reserve(d_mol.getNumAtoms());
  std::vector<unsigned int> seen(d_mol.getNumAtoms());
  unsigned int generation = 0;

  for (std::size_t pos = 0; pos < nodes.size(); ++pos) {
    auto node = nodes[pos];
    const auto atom = node->getAtom();
    if (atom != nullptr && !node->isDuplicate() &&
        targets.test(atom->getIdx())) {
      result[atom->getIdx()].push_back(node);
    }

    for (const auto edge : node->getEdges()) {
      if (!edge->isBeg(node) || isAcyclicBranchWithoutConfiguration(edge)) {
        continue;
      }
      auto child = edge->getEnd();
      if (!child->isOriginalChildOf(node) || child->isDuplicateOrH() ||
          child->getAtom() == nullptr) {
        continue;
      }
      if (targets.test(child->getAtom()->getIdx()) ||
          canReachUnvisitedTarget(child, targets, reachabilityQueue, seen,
                                  generation)) {
        nodes.push_back(child);
      }
    }
  }
  return result;
}

bool Digraph::canReachUnvisitedTarget(const Node *node,
                                      const boost::dynamic_bitset<> &targets,
                                      std::vector<unsigned int> &queue,
                                      std::vector<unsigned int> &seen,
                                      unsigned int &generation) const {
  const auto atom = node->getAtom();
  if (atom == nullptr || node->isDuplicateOrH()) {
    return false;
  }

  // Use generation marks so each small molecular reachability search only
  // clears the entries it actually visits.
  if (++generation == 0u) {
    std::fill(seen.begin(), seen.end(), 0u);
    generation = 1u;
  }
  queue.clear();
  queue.push_back(atom->getIdx());
  seen[atom->getIdx()] = generation;

  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    auto current = d_mol.getAtom(queue[pos]);
    for (auto bond : d_mol.getBonds(current)) {
      auto neighbor = bond->getOtherAtom(current);
      const auto neighborIdx = neighbor->getIdx();
      if (seen[neighborIdx] == generation || node->isVisited(neighborIdx)) {
        continue;
      }
      if (targets.test(neighborIdx)) {
        return true;
      }
      seen[neighborIdx] = generation;
      queue.push_back(neighborIdx);
    }
  }
  return false;
}

std::vector<Node *> Digraph::collectNodes(
    Atom *atom, bool pruneConfigurationFreeBranches) const {
  std::vector<Node *> result;
  std::vector<Node *> queue = {getCurrentRoot()};

  for (size_t i = 0; i < queue.size(); ++i) {
    auto node = queue[i];
    if (atom == node->getAtom()) {
      result.push_back(node);
    }
    for (const auto &e : node->getEdges()) {
      if (!e->isBeg(node)) {
        continue;
      }
      if (pruneConfigurationFreeBranches &&
          isAcyclicBranchWithoutConfiguration(e)) {
        continue;
      }
      queue.push_back(e->getEnd());
    }
  }
  return result;
}

bool Digraph::isAcyclicBranchWithoutConfiguration(const Edge *edge) const {
  if (edge == nullptr) {
    return false;
  }
  const auto end = edge->getEnd();
  if (end == nullptr || end->isDuplicateOrH() ||
      !end->isOriginalChildOf(edge->getBeg())) {
    return false;
  }
  return d_mol.isAcyclicBranchWithoutConfiguration(edge->getBond(),
                                                   end->getAtom());
}

bool Digraph::hasAuxDescriptorOnSide(const Edge *edge, unsigned mask) const {
  if (edge == nullptr || dp_origin == nullptr || mask == 0u) {
    return false;
  }
  const auto beg = edge->getBeg();
  const auto end = edge->getEnd();
  if (end->isOriginalChildOf(beg)) {
    return end->getAuxDescriptorCount(mask) != 0u;
  }
  if (beg->isOriginalChildOf(end)) {
    const auto total = dp_origin->getAuxDescriptorCount(mask);
    const auto excluded = beg->getAuxDescriptorCount(mask);
    return total > excluded;
  }
  return false;
}

void Digraph::noteConstitutionalRootEquivalence() {
  d_usedConstitutionalRootEquivalence = true;
}

bool Digraph::usedConstitutionalRootEquivalence() const {
  return d_usedConstitutionalRootEquivalence;
}

/**
 * Access the reference atom for Rule 6 (if one is set).
 */
Atom *Digraph::getRule6Ref() const { return dp_rule6Ref; }

/**
 * Used exclusively for Rule 6, we set one atom as the reference.
 * @param ref reference atom
 */
void Digraph::setRule6Ref(Atom *ref) { dp_rule6Ref = ref; }

/**
 * Sets the root node of this digraph by flipping the directions
 * of edges as required.
 *
 * @param newroot the new root
 */
void Digraph::changeRoot(Node *newroot) {
  std::vector<Edge *> toflip;
  std::vector<Node *> queue{newroot};
  for (std::size_t pos = 0; pos < queue.size(); ++pos) {
    const auto node = queue[pos];
    for (const auto &e : node->getEdges()) {
      if (e->isEnd(node)) {
        toflip.push_back(e);
        queue.push_back(e->getBeg());
      }
    }
  }
  for (auto &e : toflip) {
    e->flip();
  }
  dp_root = newroot;
}

void Digraph::expand(Node *beg) {
  const auto &atom = beg->getAtom();
  const auto &edges = beg->getEdges();
  const auto &prev =
      edges.size() > 0 && !edges[0]->isBeg(beg) ? edges[0]->getBond() : nullptr;

  if (MAX_NODE_DIST > 0 && beg->getDistance() > MAX_NODE_DIST) {
    return;
  }

  bool averaged_negative_checked = false;
  bool averaged_negative = false;
  const auto has_averaged_negative_charge = [&]() {
    if (!averaged_negative_checked) {
      averaged_negative_checked = true;
      averaged_negative = atom->getFormalCharge() < 0 &&
                          d_mol.getFractionalAtomicNum(atom).isAveraged();
    }
    return averaged_negative;
  };

  // create 'explicit' nodes
  for (const auto &bond : d_mol.getBonds(atom)) {
    const auto &nbr = bond->getOtherAtom(atom);
    const int nbrIdx = nbr->getIdx();
    const int bord = d_mol.getBondOrder(bond);
    const int virtual_nodes = bord - 1;

    if (!beg->isVisited(nbrIdx)) {
      auto end = beg->newChild(nbrIdx, nbr);
      addEdge(beg, bond, end);

      // duplicate nodes for bond orders (except for root atoms...)
      // for example >S=O
      if (dp_origin != beg || d_atropisomerMode) {
        if (has_averaged_negative_charge()) {
          end = beg->newBondDuplicateChild(nbrIdx, nbr);
          addEdge(beg, bond, end);
        } else {
          for (int i = 0; i < virtual_nodes; ++i) {
            end = beg->newBondDuplicateChild(nbrIdx, nbr);
            addEdge(beg, bond, end);
          }
        }
      }
    } else if (bond == prev) {  // bond order expansion (backwards)
      if (dp_origin->getAtom() != nbr || d_atropisomerMode) {
        for (int i = 0; i < virtual_nodes; ++i) {
          auto end = beg->newBondDuplicateChild(nbrIdx, nbr);
          addEdge(beg, bond, end);
        }
      }
    } else {  // ring closures
      auto end = beg->newRingDuplicateChild(nbrIdx, nbr);
      addEdge(beg, bond, end);

      if (has_averaged_negative_charge()) {
        end = beg->newBondDuplicateChild(nbrIdx, nbr);
        addEdge(beg, bond, end);
      } else {
        for (int i = 0; i < virtual_nodes; ++i) {
          end = beg->newBondDuplicateChild(nbrIdx, nbr);
          addEdge(beg, bond, end);
        }
      }
    }
  }

  // Create implicit hydrogen nodes
  const int hcnt = atom->getTotalNumHs();
  for (int i = 0; i < hcnt; ++i) {
    auto end = beg->newImplicitHydrogenChild();
    addEdge(beg, nullptr, end);
  }
}

}  // namespace CIPLabeler
}  // namespace RDKit
