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
#include <span>
#include <stdexcept>
#include <vector>

#include "../CIPLabeler.h"

#include "../Descriptor.h"
#include "../Edge.h"
#include "../Node.h"
#include "../Sort.h"
#include "Pairlist.h"

namespace RDKit {
namespace CIPLabeler {

class CIPMol;
class AtropisomerBond;
class Sp2Bond;
class Tetrahedral;

namespace {
template <typename T>
inline int three_way_comparison(const T &x, const T &y) {
  return x < y ? -1 : (x == y ? 0 : 1);
}
}  // namespace

class SequenceRule {
 public:
  static constexpr std::size_t MAX_CACHED_SORT_EDGES = 8;

  // Keeps exact comparison results alive across the nested sorts performed by
  // one labeling operation. Nested sessions share the cache; the outermost
  // session owns its lifetime.
  class ComparisonSession {
   public:
    ComparisonSession();
    ~ComparisonSession();

    // Auxiliary descriptors change between distance spheres. Discard only
    // descriptor-dependent entries while retaining constitutional work.
    static void invalidateAuxiliaryCaches();

    ComparisonSession(const ComparisonSession &) = delete;
    ComparisonSession &operator=(const ComparisonSession &) = delete;

   private:
    friend class SequenceRule;
    friend class Sort;
    friend class AtropisomerBond;
    friend class Sp2Bond;
    friend class Tetrahedral;
  };

  explicit SequenceRule(bool useConstitutionalRootEquivalence = false,
                        bool auxiliaryIndependent = false);

  virtual ~SequenceRule();

  Descriptor getBondLabel(const Edge *edge) const;

  int getComparision(const Edge *a, const Edge *b) const;

  virtual int getComparision(const Edge *a, const Edge *b, bool deep) const;

  virtual const Sort *getSorter() const;

  bool isAuxiliaryIndependent() const;

  int recursiveCompare(const Edge *a, const Edge *b) const;

  virtual void setSorter(const Sort *sorter);

  Priority sort(const Node *node, std::vector<Edge *> &edges, bool deep) const;

  Priority sort(const Node *node, std::vector<Edge *> &edges) const;

  virtual int compare(const Edge *a, const Edge *b) const = 0;

 protected:
  std::unique_ptr<const Sort> dp_sorter = nullptr;

  // Rules whose value can only come from a particular descriptor class can
  // prove that recursion is unnecessary from per-side descriptor counts.
  virtual bool isRecursiveComparisonNeeded(const Edge *a, const Edge *b) const;

 private:
  std::uint64_t d_cacheId;
  bool d_useConstitutionalRootEquivalence;
  bool d_auxiliaryIndependent;

  bool areUpEdges(Node *aNode, Node *bNode, Edge *aEdge, Edge *bEdge) const;
  bool hasEquivalentAcyclicContinuation(const Edge *a, const Edge *b) const;
  bool hasEquivalentConstitutionalRootContinuation(const Edge *a,
                                                   const Edge *b) const;
  int recursiveCompareEqual(const Edge *a, const Edge *b) const;

  static bool getCachedSort(std::uint64_t sortId, const Node *node, bool deep,
                            bool auxiliaryIndependent,
                            std::vector<Edge *> &edges, bool &unique,
                            bool &pseudoAsymmetric);
  static void cacheSort(std::uint64_t sortId, const Node *node, bool deep,
                        bool auxiliaryIndependent, std::span<Edge *const> input,
                        const std::vector<Edge *> &sorted,
                        const Priority &priority);

  friend class Sort;
};

}  // namespace CIPLabeler
}  // namespace RDKit
