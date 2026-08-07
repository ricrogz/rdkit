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
#include "MultithreadedMolSupplier.h"

#include <RDGeneral/RDLog.h>

#include <exception>
#include <ostream>

namespace RDKit {

namespace v2 {
namespace FileParsers {

void MultithreadedMolSupplier::initFromSettings(bool takeOwnership,
                                                const Parameters &params) {
  df_owner = takeOwnership;
  d_params = params;
  d_params.numWriterThreads = getNumThreadsToUse(params.numWriterThreads);
  d_inputQueue.reset(new inputQueue_t(d_params.sizeInputQueue));
  d_outputQueue.reset(new outputQueue_t(d_params.sizeOutputQueue));
}

void MultithreadedMolSupplier::close() {
  const std::lock_guard<std::mutex> lifecycleLock(d_lifecycleMutex);
  df_forceStop = true;

  // Mark both queues done before joining. This wakes consumers as well as
  // producers, and future pushes are rejected instead of being enqueued.
  if (d_inputQueue) {
    d_inputQueue->setDone();
  }
  if (d_outputQueue) {
    d_outputQueue->setDone();
  }

  endThreads();

  // Threads have stopped, so queue storage can now be destroyed safely.
  if (d_inputQueue) {
    d_inputQueue->clear();
  }
  if (d_outputQueue) {
    d_outputQueue->clear();
  }

  // Derived destructors call close() while their parsing state is still alive;
  // only release the input stream after every worker has stopped.
  closeStreams();
  df_started = false;
}

void MultithreadedMolSupplier::closeStreams() {
  if (df_owner && dp_inStream) {
    delete dp_inStream;
    df_owner = false;
    dp_inStream = nullptr;
  }
  df_started = false;
}

void MultithreadedMolSupplier::reader() {
  std::string record;
  unsigned int lineNum, index;
  try {
    while (!df_forceStop && extractNextRecord(record, lineNum, index)) {
      const auto callback =
          std::atomic_load_explicit(&d_readCallback, std::memory_order_acquire);
      if (callback) {
        try {
          record = (*callback)(record, index);
        } catch (const std::exception &e) {
          BOOST_LOG(rdErrorLog)
              << "Read callback exception: " << e.what() << std::endl;
        } catch (...) {
          BOOST_LOG(rdErrorLog)
              << "Unknown read callback exception" << std::endl;
        }
      }
      inputRecord_t inputRecord{std::move(record), lineNum, index};
      if (df_forceStop || !d_inputQueue->push(std::move(inputRecord))) {
        break;
      }
      ++d_readRecordCount;
    }
  } catch (const std::exception &e) {
    BOOST_LOG(rdErrorLog) << "Reader thread exception: " << e.what()
                          << std::endl;
  } catch (...) {
    BOOST_LOG(rdErrorLog) << "Unknown reader thread exception" << std::endl;
  }
  df_readerDone.store(true, std::memory_order_release);
  d_inputQueue->setDone();
}

void MultithreadedMolSupplier::writer() {
  try {
    inputRecord_t inputRecord;
    while (!df_forceStop && d_inputQueue->pop(inputRecord)) {
      if (df_forceStop.load(std::memory_order_relaxed)) {
        break;
      }
      try {
        auto mol = processMoleculeRecord(std::get<0>(inputRecord),
                                         std::get<1>(inputRecord));
        const auto callback = std::atomic_load_explicit(
            &d_writeCallback, std::memory_order_acquire);
        if (!df_forceStop && mol && callback) {
          (*callback)(*mol, std::get<0>(inputRecord), std::get<2>(inputRecord));
        }
        outputRecord_t outputRecord{std::move(mol),
                                    std::move(std::get<0>(inputRecord)),
                                    std::get<2>(inputRecord)};
        if (!d_outputQueue->push(std::move(outputRecord))) {
          break;
        }
      } catch (...) {
        // Preserve one output entry per input record, including parse and
        // write-callback failures.
        outputRecord_t nullRecord{nullptr, std::move(std::get<0>(inputRecord)),
                                  std::get<2>(inputRecord)};
        if (!d_outputQueue->push(std::move(nullRecord))) {
          break;
        }
      }
    }
  } catch (const std::exception &e) {
    df_forceStop.store(true, std::memory_order_release);
    d_inputQueue->setDone();
    d_outputQueue->setDone();
    BOOST_LOG(rdErrorLog) << "Writer thread exception: " << e.what()
                          << std::endl;
  } catch (...) {
    df_forceStop.store(true, std::memory_order_release);
    d_inputQueue->setDone();
    d_outputQueue->setDone();
    BOOST_LOG(rdErrorLog) << "Unknown writer thread exception" << std::endl;
  }

  const auto finishedWriterCount =
      d_finishedWriterCount.fetch_add(1, std::memory_order_relaxed) + 1;
  if (finishedWriterCount == d_params.numWriterThreads) {
    d_outputQueue->setDone();
  }
}

std::unique_ptr<RWMol> MultithreadedMolSupplier::next() {
  const std::lock_guard<std::mutex> nextLock(d_nextMutex);
  {
    const std::lock_guard<std::mutex> lifecycleLock(d_lifecycleMutex);
    if (!df_started && !df_forceStop) {
      df_started = true;
      try {
        startThreads();
      } catch (...) {
        df_forceStop = true;
        d_inputQueue->setDone();
        d_outputQueue->setDone();
        endThreads();
        df_started = false;
        throw;
      }
    }
  }

  outputRecord_t outputRecord;
  if (!df_forceStop && d_outputQueue->pop(outputRecord)) {
    {
      const std::lock_guard<std::mutex> itemTextLock(d_lastItemTextMutex);
      d_lastItemText = std::move(std::get<1>(outputRecord));
    }
    d_lastRecordId.store(std::get<2>(outputRecord), std::memory_order_release);
    auto res = std::move(std::get<0>(outputRecord));
    const auto callback =
        std::atomic_load_explicit(&d_nextCallback, std::memory_order_acquire);
    if (res && callback) {
      try {
        (*callback)(*res, *this);
      } catch (...) {
        // Ignore exception and proceed with mol as is.
      }
    }
    d_returnedRecordCount.fetch_add(1, std::memory_order_release);
    return res;
  }
  if (!df_forceStop && d_returnedRecordCount.load(std::memory_order_acquire)) {
    // A failed pop after natural worker completion is a physical EOF read,
    // not a null molecule record. Forward-iterator wrappers must suppress it.
    df_eofHitOnRead.store(true, std::memory_order_release);
  }
  return nullptr;
}

// this calls joins on the reader and writer threads
//  and waits until completion.  To actually force a stop
//  call close which handles the input and output queues
void MultithreadedMolSupplier::endThreads() {
  if (!df_started) {
    return;
  }

  // stop the writers before stopping the readers
  //  otherwise there might be a deadlock
  for (auto &thread : d_writerThreads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
  if (d_readerThread.joinable()) {
    d_readerThread.join();
  }
}

void MultithreadedMolSupplier::startThreads() {
  // Reserve before any thread starts so vector allocation cannot fail while
  // worker threads are already running.
  d_writerThreads.reserve(d_params.numWriterThreads);
  // run the reader function in a separate thread
  d_readerThread = std::thread(&MultithreadedMolSupplier::reader, this);
  // run the writer function in separate threads
  for (unsigned int i = 0; i < d_params.numWriterThreads; i++) {
    d_writerThreads.emplace_back(&MultithreadedMolSupplier::writer, this);
  }
}

bool MultithreadedMolSupplier::isAtEnd() const {
  if (df_forceStop.load(std::memory_order_acquire)) {
    return true;
  }
  // The reader increments its count after every successful input-queue push;
  // the release store publishes the final value. The returned count avoids
  // depending on writer-exit timing once the final output has been consumed.
  if (!df_readerDone.load(std::memory_order_acquire)) {
    return false;
  }
  return d_returnedRecordCount.load(std::memory_order_acquire) ==
         d_readRecordCount;
}

bool MultithreadedMolSupplier::atEnd() { return isAtEnd(); }

bool MultithreadedMolSupplier::getEOFHitOnRead() const {
  // Delaying the EOF flag until all output is consumed prevents forward
  // iterators from discarding a molecule that the reader prefetched before
  // discovering EOF.
  if (df_forceStop.load(std::memory_order_acquire) || !isAtEnd()) {
    return false;
  }
  return !df_forceStop.load(std::memory_order_acquire) &&
         df_eofHitOnRead.load(std::memory_order_acquire);
}

unsigned int MultithreadedMolSupplier::getLastRecordId() const {
  return d_lastRecordId.load(std::memory_order_acquire);
}

std::string MultithreadedMolSupplier::getLastItemText() const {
  const std::lock_guard<std::mutex> itemTextLock(d_lastItemTextMutex);
  return d_lastItemText;
}

void MultithreadedMolSupplier::setNextCallback(nextCallBackFn_t cb) {
  std::shared_ptr<const nextCallBackFn_t> callback;
  if (cb) {
    callback = std::make_shared<nextCallBackFn_t>(std::move(cb));
  }
  std::atomic_store_explicit(&d_nextCallback, std::move(callback),
                             std::memory_order_release);
}

void MultithreadedMolSupplier::setWriteCallback(writeCallBackFn_t cb) {
  std::shared_ptr<const writeCallBackFn_t> callback;
  if (cb) {
    callback = std::make_shared<writeCallBackFn_t>(std::move(cb));
  }
  std::atomic_store_explicit(&d_writeCallback, std::move(callback),
                             std::memory_order_release);
}

void MultithreadedMolSupplier::setReadCallback(readCallBackFn_t cb) {
  std::shared_ptr<const readCallBackFn_t> callback;
  if (cb) {
    callback = std::make_shared<readCallBackFn_t>(std::move(cb));
  }
  std::atomic_store_explicit(&d_readCallback, std::move(callback),
                             std::memory_order_release);
}

void MultithreadedMolSupplier::reset() {
  UNDER_CONSTRUCTION("reset() not supported for MultithreadedMolSupplier();");
}
}  // namespace FileParsers
}  // namespace v2
}  // namespace RDKit
#endif
