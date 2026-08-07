//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef CONTROLCHANDLER_H
#define CONTROLCHANDLER_H

#include <csignal>
#include <stdexcept>

namespace RDKit {

class ControlCCaught : public std::runtime_error {
 public:
  explicit ControlCCaught()
      : std::runtime_error("The process was interrupted with Ctrl+c"){};
};

//! This class catches a control-C/SIGINT and sets the flag d_gotSignal
//! if one is received.  It is intended to be used inside a long
//! C++ calculation called from Python which intercepts the signal
//! handler.  The C++ code must check the value of d_gotSignal
//! periodically and act accordingly. Nested instances share the handler and
//! flag installed by the first instance; only the last instance resets them.
//! This allows an outer operation to observe an interrupt caught by an inner
//! operation. Resetting the signal handler and flag after the last instance is
//! essential because they are static variables. Handler scopes must all run on
//! the same thread.
//! Example usage, inside a boost::python wrapper:
//!  ResultsObject results;
//!  {
//!   NOGIL gil;
//!    results = someFunction();
//!  }
//!  if (results.getCancelled()) {
//!    throw_runtime_error("someFunction cancelled");
//!  }
//! It's important that the exception is thrown once the GIL has been
//! released, otherwise a crash is inevitable at some future point.
class ControlCHandler {
 public:
  ControlCHandler() {
    if (d_nestingDepth++ == 0) {
      d_prev_handler = std::signal(SIGINT, signalHandler);
    }
  }
  ControlCHandler(const ControlCHandler &) = delete;
  ControlCHandler(ControlCHandler &&) = delete;
  ControlCHandler &operator=(const ControlCHandler &) = delete;
  ControlCHandler &operator=(ControlCHandler &&) = delete;
  ~ControlCHandler() {
    if (--d_nestingDepth == 0) {
      std::signal(SIGINT, d_prev_handler);
      d_gotSignal = 0;
    }
  }
  static bool getGotSignal() { return d_gotSignal != 0; }
  static void signalHandler(int signalNumber) {
    if (signalNumber == SIGINT) {
      d_gotSignal = 1;
      std::signal(SIGINT, d_prev_handler);
    }
  }
  static void reset() {
    d_gotSignal = 0;
    std::signal(SIGINT, signalHandler);
  }

 private:
  // sig_atomic_t is the standard type that may be read and written by an
  // asynchronous signal handler in this otherwise single-threaded lifecycle.
  inline static volatile std::sig_atomic_t d_gotSignal{0};
  inline static unsigned int d_nestingDepth{0};
  inline static void (*d_prev_handler)(int);
};
}  // namespace RDKit
#endif  // CONTROLCHANDLER_H
