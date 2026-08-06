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

#include <cstdint>
#include <stdexcept>
#include <memory>
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
  // Keeps exact comparison results alive across the nested sorts performed by
  // one configuration-label calculation. Nested sessions share the cache;
  // the outermost session owns its lifetime.
  class ComparisonSession {
   public:
    ComparisonSession(const ComparisonSession &) = delete;
    ComparisonSession &operator=(const ComparisonSession &) = delete;

   private:
    ComparisonSession();
    ~ComparisonSession();

    friend class SequenceRule;
    friend class Sort;
    friend class AtropisomerBond;
    friend class Sp2Bond;
    friend class Tetrahedral;
  };

  SequenceRule();

  virtual ~SequenceRule();

  Descriptor getBondLabel(const Edge *edge) const;

  int getComparision(const Edge *a, const Edge *b) const;

  virtual int getComparision(const Edge *a, const Edge *b, bool deep) const;

  virtual const Sort *getSorter() const;

  int recursiveCompare(const Edge *a, const Edge *b) const;

  void setSorter(const Sort *sorter);

  Priority sort(const Node *node, std::vector<Edge *> &edges, bool deep) const;

  Priority sort(const Node *node, std::vector<Edge *> &edges) const;

  virtual int compare(const Edge *a, const Edge *b) const = 0;

 protected:
  std::unique_ptr<const Sort> dp_sorter = nullptr;

 private:
  std::uint64_t d_cacheId;

  bool areUpEdges(Node *aNode, Node *bNode, Edge *aEdge, Edge *bEdge) const;
  bool hasEquivalentAcyclicContinuation(const Edge *a, const Edge *b) const;
  int recursiveCompareEqual(const Edge *a, const Edge *b) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
