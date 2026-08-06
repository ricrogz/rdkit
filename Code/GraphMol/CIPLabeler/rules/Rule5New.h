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
class Rule5New : public SequenceRule {
 public:
  Rule5New();

  Rule5New(Descriptor ref);

  int compare(const Edge *a, const Edge *b) const override;

  void setSorter(const Sort *sorter) override;

 protected:
  bool isRecursiveComparisonNeeded(const Edge *a, const Edge *b) const override;

 private:
  struct ReferenceRuleTag {};

  Rule5New(Descriptor ref, ReferenceRuleTag);

  const Descriptor d_ref = Descriptor::NONE;
  std::unique_ptr<Rule5New> dp_referenceR;
  std::unique_ptr<Rule5New> dp_referenceS;
  std::unique_ptr<const Sort> dp_referenceSorterR;
  std::unique_ptr<const Sort> dp_referenceSorterS;

  void fillPairs(const Node *beg, PairList &plist,
                 std::vector<const Node *> &queue,
                 std::vector<Edge *> &edges) const;

  const Sort &getRefSorter(Descriptor ref) const;
  std::unique_ptr<const Sort> makeRefSorter(
      const SequenceRule *replacementRule) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
