#ifdef RDK_BUILD_THREADSAFE_SSS
//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MultithreadedSmilesMolSupplier.h"

namespace RDKit {
namespace v2 {
namespace FileParsers {
MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    const std::string &fileName, const Parameters &params,
    const SmilesMolSupplierParams &parseParams) {
  dp_inStream = openAndCheckStream(fileName);
  CHECK_INVARIANT(dp_inStream, "bad instream");
  CHECK_INVARIANT(!(dp_inStream->eof()), "early EOF");
  // set df_takeOwnership = true
  initFromSettings(true, params, parseParams);
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier(
    std::istream *inStream, bool takeOwnership, const Parameters &params,
    const SmilesMolSupplierParams &parseParams) {
  CHECK_INVARIANT(inStream, "bad instream");
  CHECK_INVARIANT(!(inStream->eof()), "early EOF");
  dp_inStream = inStream;
  initFromSettings(takeOwnership, params, parseParams);
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSmilesMolSupplier::MultithreadedSmilesMolSupplier() {
  dp_inStream = nullptr;
  initFromSettings(true, d_params, d_parseParams);
}

void MultithreadedSmilesMolSupplier::closeStreams() {
  if (df_owner && dp_inStream) {
    delete dp_inStream;
    df_owner = false;
    dp_inStream = nullptr;
  }
}

void MultithreadedSmilesMolSupplier::initFromSettings(
    bool takeOwnership, const Parameters &params,
    const SmilesMolSupplierParams &parseParams) {
  df_owner = takeOwnership;
  d_params = params;
  d_parseParams = parseParams;
  d_params.numWriterThreads = getNumThreadsToUse(d_params.numWriterThreads);
  d_inputQueue.reset(
      new ConcurrentQueue<std::tuple<std::string, unsigned int, unsigned int>>(
          d_params.sizeInputQueue));
  d_outputQueue.reset(
      new ConcurrentQueue<std::tuple<RWMol *, std::string, unsigned int>>(
          d_params.sizeOutputQueue));
  d_line = -1;
}

// --------------------------------------------------
//
//  Reads and processes the title line
//
void MultithreadedSmilesMolSupplier::processTitleLine() {
  PRECONDITION(dp_inStream, "bad stream");

  std::string tempStr;
  // loop until we get a valid line
  do {
    tempStr = getLine(dp_inStream);
    ++d_line;
  } while (!dp_inStream->eof() && !dp_inStream->fail() &&
           (tempStr[0] == '#' || strip(tempStr).empty()));

  boost::char_separator<char> sep(d_parseParams.delimiter.c_str(), "",
                                  boost::keep_empty_tokens);
  tokenizer tokens(tempStr, sep);
  for (tokenizer::iterator tokIter = tokens.begin(); tokIter != tokens.end();
       ++tokIter) {
    std::string pname = strip(*tokIter);
    d_props.push_back(pname);
  }
}

bool MultithreadedSmilesMolSupplier::extractNextRecord(std::string &record,
                                                       unsigned int &lineNum) {
  PRECONDITION(dp_inStream, "bad stream");
  if (dp_inStream->eof()) {
    return false;
  }

  // need to process title line
  // if we have not called next yet and the read count = 0
  // then we are seeking the first record
  if (d_readCount == 0 && d_parseParams.titleLine) {
    this->processTitleLine();
  }

  do {
    record = getLine(dp_inStream);
    ++d_line;
  } while (!dp_inStream->eof() && !dp_inStream->fail() &&
           (record[0] == '#' || strip(record).empty()));

  lineNum = d_line;

  return true;
}

RWMol *MultithreadedSmilesMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  // -----------
  // tokenize the input line:
  // -----------
  boost::char_separator<char> sep(d_parseParams.delimiter.c_str(), "",
                                  boost::keep_empty_tokens);
  tokenizer tokens(record, sep);
  STR_VECT recs;
  for (tokenizer::iterator tokIter = tokens.begin(); tokIter != tokens.end();
       ++tokIter) {
    std::string rec = strip(*tokIter);
    recs.push_back(rec);
  }
  if (recs.size() <= static_cast<unsigned int>(d_parseParams.smilesColumn)) {
    std::ostringstream errout;
    errout << "ERROR: line #" << lineNum << "does not contain enough tokens\n";
    throw FileParseException(errout.str());
  }

  // -----------
  // get the smiles and create a molecule
  // -----------
  auto res = MolFromSmiles(recs[d_parseParams.smilesColumn],
                           d_parseParams.parseParameters);
  if (!res) {
    std::stringstream errout;
    errout << "Cannot create molecule from : '"
           << recs[d_parseParams.smilesColumn] << "'";
    throw SmilesParseException(errout.str());
  }

  // -----------
  // get the name (if there's a name column)
  // -----------
  if (d_parseParams.nameColumn == -1) {
    // if no name defaults it to the line number we read it from string
    std::ostringstream tstr;
    tstr << lineNum;
    std::string mname = tstr.str();
    res->setProp(common_properties::_Name, mname);
  } else {
    if (d_parseParams.nameColumn >= static_cast<int>(recs.size())) {
      BOOST_LOG(rdWarningLog)
          << "WARNING: no name column found on line " << lineNum << std::endl;
    } else {
      res->setProp(common_properties::_Name, recs[d_parseParams.nameColumn]);
    }
  }

  // -----------
  // read in the properties
  // -----------
  for (unsigned int col = 0; col < recs.size(); col++) {
    if (static_cast<int>(col) == d_parseParams.smilesColumn ||
        static_cast<int>(col) == d_parseParams.nameColumn) {
      continue;
    }
    std::string pname, pval;
    if (d_props.size() > col) {
      pname = d_props[col];
    } else {
      pname = "Column_";
      std::stringstream ss;
      ss << col;
      pname += ss.str();
    }

    pval = recs[col];
    res->setProp(pname, pval);
  }
  return res.release();
}

}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
#endif
