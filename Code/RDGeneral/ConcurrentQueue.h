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
#ifndef CONCURRENT_QUEUE
#define CONCURRENT_QUEUE
#include <cstddef>
#include <condition_variable>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace RDKit {
template <typename E>
class ConcurrentQueue {
 private:
  size_t d_capacity;
  bool d_done;
  std::vector<E> d_elements;
  size_t d_head, d_tail, d_size;
  mutable std::mutex d_lock;
  std::condition_variable d_notEmpty, d_notFull;

 private:
  ConcurrentQueue(const ConcurrentQueue<E> &);
  ConcurrentQueue &operator=(const ConcurrentQueue<E> &);

 public:
  explicit ConcurrentQueue(size_t capacity)
      : d_capacity(capacity),
        d_done(false),
        d_elements(capacity),
        d_head(0),
        d_tail(0),
        d_size(0) {
    if (!capacity) {
      throw std::invalid_argument("ConcurrentQueue capacity must be nonzero");
    }
  }

  //! pushes an element into the queue, blocking while the queue is full
  //! returns false without enqueueing the element if the queue is done
  bool push(E element);

  //! pops an existing element, even after the queue has been marked done
  //! returns false when the queue is both empty and done; blocks while it is
  //! empty but not done
  bool pop(E &element);

  //! checks whether the ConcurrentQueue is empty
  bool isEmpty() const;

  //! returns the value of the variable done
  bool getDone() const;

  //! sets the variable d_done = true
  void setDone();

  //! clears the vector
  void clear();
};

template <typename E>
bool ConcurrentQueue<E>::push(E element) {
  std::unique_lock<std::mutex> lk(d_lock);
  d_notFull.wait(lk, [this]() { return d_done || d_size < d_capacity; });
  if (d_done) {
    return false;
  }
  d_elements.at(d_tail) = std::move(element);
  d_tail = (d_tail + 1) % d_capacity;
  ++d_size;
  lk.unlock();
  d_notEmpty.notify_one();
  return true;
}

template <typename E>
bool ConcurrentQueue<E>::pop(E &element) {
  std::unique_lock<std::mutex> lk(d_lock);
  d_notEmpty.wait(lk, [this]() { return d_done || d_size != 0; });
  if (!d_size) {
    return false;
  }
  element = std::move(d_elements.at(d_head));
  d_head = (d_head + 1) % d_capacity;
  --d_size;
  lk.unlock();
  d_notFull.notify_one();
  return true;
}

template <typename E>
bool ConcurrentQueue<E>::isEmpty() const {
  std::lock_guard<std::mutex> lk(d_lock);
  return !d_size;
}

template <typename E>
bool ConcurrentQueue<E>::getDone() const {
  std::lock_guard<std::mutex> lk(d_lock);
  return d_done;
}

template <typename E>
void ConcurrentQueue<E>::setDone() {
  {
    std::lock_guard<std::mutex> lk(d_lock);
    d_done = true;
  }
  d_notEmpty.notify_all();
  d_notFull.notify_all();
}

template <typename E>
void ConcurrentQueue<E>::clear() {
  std::lock_guard<std::mutex> lk(d_lock);
  d_elements.clear();
  d_head = 0;
  d_tail = 0;
  d_size = 0;
}

}  // namespace RDKit
#endif
#endif
