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

#include <initializer_list>
#include <vector>

#include "Rule5.h"
#include "SequenceRule.h"

namespace RDKit {
namespace CIPLabeler {

class Rules : public SequenceRule {
 public:
  Rules() = delete;

  Rules(std::initializer_list<SequenceRule *> rules) {
    for (auto &rule : rules) {
      add(rule);
    }
    // SequenceRule's eagerly constructed sorter stores this composite
    // comparator and will therefore continue to see rules added through add().
  }

  ~Rules() override {
    for (auto &rule : d_rules) {
      delete rule;
    }
  }

  void add(SequenceRule *rule) {
    if (rule == nullptr) {
      throw std::runtime_error("No sequence rule provided");
    }
    d_rules.push_back(rule);
    rule->setSorter(new Sort(d_rules));
  }

  int getNumSubRules() const { return d_rules.size(); }

  const Sort *getSorter() const override { return dp_sorter.get(); }

  int compare(const Edge *o1, const Edge *o2) const override {
    // Try using each rules. The rules will expand the search exhaustively
    // to all child substituents
    for (const auto &rule : d_rules) {
      // compare expands exhaustively across the whole graph
      int value = rule->recursiveCompare(o1, o2);
      if (value != 0) {
        return value;
      }
    }
    return 0;
  }

  int getComparision(const Edge *a, const Edge *b, bool deep) const override {
    (void)deep;

    // Try using each rules. The rules will expand the search exhaustively
    // to all child substituents
    for (const auto &rule : d_rules) {
      // compare expands exhaustively across the whole graph
      int value = rule->recursiveCompare(a, b);

      if (value != 0) {
        return value;
      }
    }

    return 0;
  }

 private:
  std::vector<const SequenceRule *> d_rules;
};

}  // namespace CIPLabeler
}  // namespace RDKit
