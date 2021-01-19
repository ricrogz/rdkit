//
//  Copyright (C) 2001-2018 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <exception>
#include <iostream>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

class CustomException : public std::runtime_error {
 public:
  CustomException(const char *msg) : std::runtime_error(msg){};
  CustomException(const std::string &msg) : std::runtime_error(msg){};
  ~CustomException() noexcept {};
};

void trigger_error() {
  // Try to parse a SMILES which is known to fail,
  // so we get an exception to catch
  try {
    auto mol = SmilesToMol("[NH6]");
  } catch (const AtomValenceException &) {
    throw CustomException("This is the exception we expect to see");
  }
}

int main() {
  try {
    trigger_error();
  } catch (const CustomException &) {
    // Just ignore the exception we expect.
  }
  return 0;
}
