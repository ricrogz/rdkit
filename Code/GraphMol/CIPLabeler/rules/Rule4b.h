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

#include "SequenceRule.h"

namespace RDKit {
namespace CIPLabeler {

/**
 * A descriptor pair rule. This rule defines that like descriptor pairs have
 * priority over unlike descriptor pairs.
 *
 */
class Rule4b : public SequenceRule {
 public:
  Rule4b();

  Rule4b(Descriptor ref);

  int compare(const Edge *a, const Edge *b) const override;

  void setSorter(const Sort *sorter) override;

 protected:
  bool isRecursiveComparisonNeeded(const Edge *a, const Edge *b) const override;

 private:
  struct ReferenceRuleTag {};

  Rule4b(Descriptor ref, ReferenceRuleTag);

  const Descriptor d_ref = Descriptor::NONE;
  std::unique_ptr<Rule4b> dp_referenceR;
  std::unique_ptr<Rule4b> dp_referenceS;
  std::unique_ptr<const Sort> dp_referenceSorterR;
  std::unique_ptr<const Sort> dp_referenceSorterS;

  std::vector<Descriptor> getReferenceDescriptors(const Node *node) const;

  bool getReference(const std::vector<const Node *> &nodes,
                    std::vector<Descriptor> &result) const;

  std::vector<std::vector<const Node *>> initialLevel(const Node *node) const;

  std::vector<std::vector<const Node *>> getNextLevel(
      const std::vector<std::vector<const Node *>> &prevLevel) const;

  std::vector<PairList> newPairLists(
      const std::vector<Descriptor> &descriptors) const;

  void fillPairs(const Node *beg, PairList &plist,
                 std::vector<const Node *> &queue,
                 std::vector<Edge *> &edges) const;

  int comparePairs(const Node *a, const Node *b, Descriptor refA,
                   Descriptor refB) const;

  const Sort &getRefSorter(Descriptor ref) const;
  std::unique_ptr<const Sort> makeRefSorter(
      const SequenceRule *replacementRule) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
