//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_BUILD_THREADSAFE_SSS
#ifndef MULTITHREADED_SD_MOL_SUPPLIER
#define MULTITHREADED_SD_MOL_SUPPLIER
#include "MultithreadedMolSupplier.h"
namespace RDKit {
namespace v2 {
namespace FileParsers {

//! This class is still a bit experimental and the public API may change
//! in future releases.
class RDKIT_FILEPARSERS_EXPORT MultithreadedSDMolSupplier
    : public MultithreadedMolSupplier {
 public:
  explicit MultithreadedSDMolSupplier(
      const std::string &fileName, const Parameters &params = Parameters(),
      const MolFileParserParams &parseParams = MolFileParserParams());

  explicit MultithreadedSDMolSupplier(
      std::istream *inStream, bool takeOwnership = true,
      const Parameters &params = Parameters(),
      const MolFileParserParams &parseParams = MolFileParserParams());

  MultithreadedSDMolSupplier();
  virtual ~MultithreadedSDMolSupplier() { close(); }

  void setProcessPropertyLists(bool val) { df_processPropertyLists = val; }
  bool getProcessPropertyLists() const { return df_processPropertyLists; }
  bool getEOFHitOnRead() const { return df_eofHitOnRead; }

  //! reads next record and returns whether or not EOF was hit
  bool extractNextRecord(std::string &record, unsigned int &lineNum) override;
  void readMolProps(RWMol &mol, std::istringstream &inStream);
  //! parses the record and returns the resulting molecule
  RWMol *processMoleculeRecord(const std::string &record,
                               unsigned int lineNum) override;

 private:
  void initFromSettings(bool takeOwnership, const Parameters &params,
                        const MolFileParserParams &parseParams);

  bool df_processPropertyLists = true;
  bool df_eofHitOnRead = false;
  MolFileParserParams d_parseParams;
};
}  // namespace FileParsers
}  // namespace v2

inline namespace v1 {
class RDKIT_FILEPARSERS_EXPORT MultithreadedSDMolSupplier
    : public v1::MultithreadedMolSupplier {
  //! this is an abstract base class to concurrently supply molecules one at a
  //! time
 public:
  using ContainedType = v2::FileParsers::MultithreadedSDMolSupplier;
  MultithreadedSDMolSupplier() {}

  explicit MultithreadedSDMolSupplier(
      const std::string &fileName, bool sanitize = true, bool removeHs = true,
      bool strictParsing = true, unsigned int numWriterThreads = 1,
      size_t sizeInputQueue = 5, size_t sizeOutputQueue = 5) {
    ContainedType::Parameters params;
    params.numWriterThreads = numWriterThreads;
    params.sizeInputQueue = sizeInputQueue;
    params.sizeOutputQueue = sizeOutputQueue;
    v2::FileParsers::MolFileParserParams parseParams;
    parseParams.sanitize = sanitize;
    parseParams.removeHs = removeHs;
    parseParams.strictParsing = strictParsing;

    dp_supplier.reset(new ContainedType(fileName, params, parseParams));
  }

  explicit MultithreadedSDMolSupplier(
      std::istream *inStream, bool takeOwnership = true, bool sanitize = true,
      bool removeHs = true, bool strictParsing = true,
      unsigned int numWriterThreads = 1, size_t sizeInputQueue = 5,
      size_t sizeOutputQueue = 5) {
    ContainedType::Parameters params;
    params.numWriterThreads = numWriterThreads;
    params.sizeInputQueue = sizeInputQueue;
    params.sizeOutputQueue = sizeOutputQueue;
    v2::FileParsers::MolFileParserParams parseParams;
    parseParams.sanitize = sanitize;
    parseParams.removeHs = removeHs;
    parseParams.strictParsing = strictParsing;

    dp_supplier.reset(
        new ContainedType(inStream, takeOwnership, params, parseParams));
  }

  void setProcessPropertyLists(bool val) {
    PRECONDITION(dp_supplier, "no supplier");
    static_cast<ContainedType *>(dp_supplier.get())
        ->setProcessPropertyLists(val);
  }
  bool getProcessPropertyLists() const {
    PRECONDITION(dp_supplier, "no supplier");
    return static_cast<ContainedType *>(dp_supplier.get())
        ->getProcessPropertyLists();
  }
};
}  // namespace v1
}  // namespace RDKit
#endif
#endif
