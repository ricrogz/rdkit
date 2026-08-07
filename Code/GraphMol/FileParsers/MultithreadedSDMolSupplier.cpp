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
#include "MultithreadedSDMolSupplier.h"

#include "FileParserUtils.h"

#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace RDKit {
namespace v2 {
namespace FileParsers {
MultithreadedSDMolSupplier::MultithreadedSDMolSupplier(
    const std::string &fileName, const Parameters &params,
    const MolFileParserParams &parseParams) {
  dp_inStream = openAndCheckStream(fileName);
  initFromSettings(true, params, parseParams);
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSDMolSupplier::MultithreadedSDMolSupplier(
    std::istream *inStream, bool takeOwnership, const Parameters &params,
    const MolFileParserParams &parseParams) {
  PRECONDITION(inStream, "bad stream");
  dp_inStream = inStream;
  initFromSettings(takeOwnership, params, parseParams);
  POSTCONDITION(dp_inStream, "bad instream");
}

MultithreadedSDMolSupplier::MultithreadedSDMolSupplier() {
  dp_inStream = nullptr;
  initFromSettings(false, d_params, d_parseParams);
}

void MultithreadedSDMolSupplier::initFromSettings(
    bool takeOwnership, const Parameters &params,
    const MolFileParserParams &parseParams) {
  MultithreadedMolSupplier::initFromSettings(takeOwnership, params);
  d_parseParams = parseParams;
  df_processPropertyLists.store(true, std::memory_order_relaxed);
  d_line = 0;
}

bool MultithreadedSDMolSupplier::extractNextRecord(std::string &record,
                                                   unsigned int &lineNum,
                                                   unsigned int &index) {
  PRECONDITION(dp_inStream, "no stream");
  if (dp_inStream->eof()) {
    return false;
  }

  bool readAnyLine = false;
  std::string currentStr;
  std::string prevStr;
  record.clear();
  lineNum = d_line;

  // keep reading while we can, and the current line is not a record separator,
  // and the previous one is a valid end of a mol block (blank line or M  END)
  while (!dp_inStream->eof() && !dp_inStream->fail() &&
         ((prevStr.find_first_not_of(" \t\r\n") != std::string::npos &&
           prevStr.find("M  END") != 0) ||
          !currentStr.starts_with("$$$$"))) {
    std::swap(prevStr, currentStr);
    if (std::getline(*dp_inStream, currentStr)) {
      readAnyLine = true;
      record += currentStr;
      record.push_back('\n');
    }
    ++d_line;
  }

  // A truly empty stream is logical EOF. If getline() successfully read one
  // or more blank lines, preserve them as a null record, matching
  // ForwardSDMolSupplier.
  if (!readAnyLine) {
    if (dp_inStream->eof() && d_lastReadRecordId == 0) {
      // Match ForwardSDMolSupplier's empty-input behavior. Do not set this
      // after a real record: the multithreaded supplier has already prefetched
      // EOF at that point, and the Python wrapper would otherwise discard the
      // final molecule.
      df_eofHitOnRead = true;
    }
    return false;
  }

  index = ++d_lastReadRecordId;
  return true;
}

void MultithreadedSDMolSupplier::readMolProps(RWMol &mol,
                                              std::istringstream &inStream) {
  PRECONDITION(inStream, "no stream");
  const bool processPropertyLists =
      df_processPropertyLists.load(std::memory_order_relaxed);
  bool hasProp = false;
  bool warningIssued = false;
  std::string dlabel;
  std::string inputLine;
  std::getline(inStream, inputLine);
  std::string_view tempStr = inputLine;

  // FIX: report files missing the $$$$ marker
  while (!inStream.eof() && !inStream.fail() && !tempStr.starts_with("$$$$")) {
    tempStr = FileParserUtils::strip(tempStr);
    if (!tempStr.empty()) {
      if (tempStr.front() == '>') {  // data header line: start of a data item
        // ignore all other crap and seek for for a data label enclosed
        // by '<' and '>'
        // FIX: "CTfile.pdf" (page 51) says that the data header line does not
        // have to contain a data label (instead can have something line field
        // id into a MACCS db). But we do not currently know what to do in this
        // situation - so ignore such data items for now
        hasProp = true;
        warningIssued = false;
        tempStr.remove_prefix(1);        // remove the first ">" sign
        size_t sl = tempStr.find('<');   // begin datalabel
        size_t se = tempStr.rfind('>');  // end datalabel
        if ((sl == std::string_view::npos) || (se == std::string_view::npos) ||
            (se == (sl + 1))) {
          // we either do not have a data label or the label is empty
          // no data label ignore until next data item
          // i.e. until we hit a blank line
          std::getline(inStream, inputLine);
          tempStr = inputLine;
          auto stmp = FileParserUtils::strip(tempStr);
          while (!stmp.empty()) {
            std::getline(inStream, inputLine);
            tempStr = inputLine;
            if (inStream.eof()) {
              throw FileParseException("End of data field name not found");
            }
            stmp = FileParserUtils::strip(tempStr);
          }
        } else {
          dlabel = std::string(tempStr.substr(sl + 1, se - sl - 1));
          // we know the label - now read in the relevant properties
          // until we hit a blank line
          std::getline(inStream, inputLine);
          tempStr = inputLine;

          std::string prop;
          auto stmp = FileParserUtils::strip(tempStr);
          int nplines = 0;  // number of lines for this property
          while (!stmp.empty() ||
                 (!tempStr.empty() &&
                  (tempStr.front() == ' ' || tempStr.front() == '\t'))) {
            nplines++;
            if (nplines > 1) {
              prop += "\n";
            }
            // take off \r if it's still in the property:
            if (tempStr.back() == '\r') {
              tempStr.remove_suffix(1);
            }
            prop.append(tempStr);
            // erase inputLine in case the file does not end with a carriage
            // return (we will end up in an infinite loop if we don't do
            // this and we do not check for EOF in this while loop body)
            inputLine.clear();
            std::getline(inStream, inputLine);
            tempStr = inputLine;
            stmp = FileParserUtils::strip(tempStr);
          }
          mol.setProp(dlabel, prop);
          if (processPropertyLists) {
            // apply this as an atom property list if that's appropriate
            FileParserUtils::processMolPropertyList(mol, dlabel);
          }
        }
      } else {
        if (d_parseParams.strictParsing) {
          // at this point we should always be at a line starting with '>'
          // following a blank line. If this is not true and df_strictParsing
          // is true, then throw an exception, otherwise truncate the rest of
          // the data field following the blank line until the next '>' or EOF
          // and issue a warning
          // FIX: should we be deleting the molecule (which is probably fine)
          // because we couldn't read the data ???
          throw FileParseException("Problems encountered parsing data fields");
        } else {
          if (!warningIssued) {
            if (hasProp) {
              BOOST_LOG(rdWarningLog)
                  << "Property <" << dlabel
                  << "> will be truncated after the first blank line\n";
            } else {
              BOOST_LOG(rdWarningLog) << "Spurious data before the first "
                                         "property will be ignored\n";
            }
            warningIssued = true;
          }
        }
      }
    }
    std::getline(inStream, inputLine);
    tempStr = inputLine;
  }
}

std::unique_ptr<RWMol> MultithreadedSDMolSupplier::processMoleculeRecord(
    const std::string &record, unsigned int lineNum) {
  PRECONDITION(dp_inStream, "no stream");
  std::istringstream inStream(record);
  auto res =
      v2::FileParsers::MolFromMolDataStream(inStream, lineNum, d_parseParams);
  if (res) {
    this->readMolProps(*res, inStream);
  }
  return res;
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
#endif
