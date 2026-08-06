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

#include <string>
#include <stdexcept>

namespace RDKit {
namespace CIPLabeler {

/**
 * Defines a descriptor which can be assigned to an atom to indicate the type of
 * chirality (if there is any). Each descriptor defines its general @{link
 * Type} which can be useful when comparing centres of different geometry.
 *
 */
enum class Descriptor {
  NONE,  // Unspecified
  UNKNOWN,
  ns,  // Other
  /**
   * Tetrahedral
   */
  R,
  S,
  r,
  s,
  /**
   * Cis/Trans
   */
  seqTrans,
  seqCis,
  E,
  Z,
  /* Axial */
  M,
  P,
  m,
  p,

  SP_4,
  TBPY_5,
  OC_6
};

// Descriptor classes used to prove that a ligand branch cannot contribute to
// a particular sequence rule. A descriptor belongs to at most one class.
inline constexpr unsigned AUX_DESCRIPTOR_GEOMETRIC = 0x1;  // E/Z (Rule 3)
inline constexpr unsigned AUX_DESCRIPTOR_PAIR =
    0x2;  // R/S/M/P/seqCis/seqTrans (Rules 4b/5)
inline constexpr unsigned AUX_DESCRIPTOR_PSEUDO = 0x4;  // r/s/m/p (Rule 4c)
inline constexpr unsigned AUX_DESCRIPTOR_ANY =
    AUX_DESCRIPTOR_GEOMETRIC | AUX_DESCRIPTOR_PAIR | AUX_DESCRIPTOR_PSEUDO;

inline unsigned getAuxDescriptorClass(Descriptor descriptor) {
  switch (descriptor) {
    case Descriptor::E:
    case Descriptor::Z:
      return AUX_DESCRIPTOR_GEOMETRIC;
    case Descriptor::R:
    case Descriptor::S:
    case Descriptor::M:
    case Descriptor::P:
    case Descriptor::seqCis:
    case Descriptor::seqTrans:
      return AUX_DESCRIPTOR_PAIR;
    case Descriptor::r:
    case Descriptor::s:
    case Descriptor::m:
    case Descriptor::p:
      return AUX_DESCRIPTOR_PSEUDO;
    default:
      return 0;
  }
}

static std::string to_string(const Descriptor &desc) {
  switch (desc) {
    case Descriptor::NONE:
      return "NONE";
    case Descriptor::UNKNOWN:
      return "UNKNOWN";
    case Descriptor::ns:
      return "ns";
    case Descriptor::R:
      return "R";
    case Descriptor::S:
      return "S";
    case Descriptor::r:
      return "r";
    case Descriptor::s:
      return "s";
    case Descriptor::seqTrans:
      return "e";
    case Descriptor::seqCis:
      return "z";
    case Descriptor::E:
      return "E";
    case Descriptor::Z:
      return "Z";
    case Descriptor::M:
      return "M";
    case Descriptor::P:
      return "P";
    case Descriptor::m:
      return "m";
    case Descriptor::p:
      return "p";
    case Descriptor::SP_4:
      return "SP_4";
    case Descriptor::TBPY_5:
      return "TBPY_5";
    case Descriptor::OC_6:
      return "OC_6";
  }
  throw std::runtime_error("Unknown descriptor");
}

}  // namespace CIPLabeler
}  // namespace RDKit
