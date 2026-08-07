# Multithreaded molecule supplier review

## Scope

This review covers the changes on top of `refactor_multithread_suppliers` in
the multithreaded molecule suppliers, together with the queue behavior those
suppliers rely on. The expected record-level behavior is:

- `MultithreadedSDMolSupplier` matches `ForwardSDMolSupplier`.
- `MultithreadedSmilesMolSupplier` matches `SmilesMolSupplier`.
- Multithreading may change output order when more than one writer is used,
  but it must not add, lose, or duplicate records.

The multithreaded suppliers remain lazy: their worker threads start on the
first call to `next()`. Consequently, an empty `SmilesMolSupplier` can report
`atEnd()` during construction while its multithreaded counterpart discovers
the same condition on the first read. Their returned record streams and
Python iteration behavior are kept consistent.

## Termination and EOF handling

The original output-queue test performed two independently locked operations:
it checked whether the queue was empty and then whether it was done. A writer
could enqueue its last result between those operations and then mark the queue
done, causing `atEnd()` to return true while a molecule was still queued.

The replacement termination protocol uses two counts:

1. The reader counts records only after they have been successfully placed on
   the input queue.
2. It publishes the final count with a release store before marking the input
   queue done.
3. `next()` atomically counts every output entry it removes, including null
   entries for malformed records.
4. `atEnd()` uses an acquire load of reader completion before comparing the
   published and returned counts.

This also avoids waiting for all writer threads to exit after the final output
has already been consumed once the reader has published its final count. If a
caller reaches physical queue EOF during the smaller interval before that
publication becomes visible, `next()` marks the failed pop as EOF so the
forward-iterator wrapper does not expose an extra Python `None`. A null entry
actually produced for a malformed record is not suppressed; direct C++ calls
can still receive `nullptr` when they explicitly read physical EOF, as with the
other forward suppliers.

`getEOFHitOnRead()` is virtual and const again, preserving the previous API,
and delays a reader EOF indication until every published record has been
returned. This is necessary because the reader can discover EOF while writers
still have valid molecules in flight. Both concrete suppliers inherit this
common accounting instead of using an SD-specific prefetch flag or a SMILES
override that always reports `false`.

Record IDs read by the reader are now reader-thread-only state, so they do not
pay for unnecessary atomic read/modify/write operations. The last returned
record ID remains atomic, restoring the thread-safety guarantee that existed
before the reviewed branch changed it to a plain integer. Access to the last
record text is protected by a mutex.

## Queue shutdown and ownership

`ConcurrentQueue::setDone()` previously woke only consumers. A producer
blocked on a full output queue could therefore remain blocked forever during
`close()`. The queue now has these semantics:

- `setDone()` wakes both producers and consumers.
- A push waiting when the queue becomes done returns `false` without adding an
  element.
- Existing queued elements can still be drained after `setDone()`.
- A zero-capacity queue is rejected instead of deadlocking on its first push.
- Head, tail, and size are maintained as a bounded ring instead of unbounded
  counters.
- One waiter is notified for each successful push or pop, avoiding the former
  notify-all thundering herd while still waking enough waiters.

The output queue now owns `std::unique_ptr<RWMol>` values rather than raw
pointers. Closing a supplier marks both queues done, joins every joinable
thread, and only then clears the queues. This removes manual deletion paths,
prevents leaks on exceptions, and makes a cancelled push destroy its molecule
safely. It also fixes the early-close deadlock where late writer results could
refill a bounded output queue after the old one-time drain.

Thread-vector capacity is reserved before any thread starts. Thread creation
failures mark both queues done and join whichever threads were successfully
created, rather than leaving a partially started supplier that can block on a
later call. Writers also recheck cancellation after waking from the input
queue, avoiding unnecessary parsing when shutdown races with a pop.

## Callback and public-state thread safety

Callback setters previously assigned to `std::function` objects while reader
or writer threads could read those same objects. The callbacks are now stored
as shared snapshots accessed through the standard atomic `shared_ptr`
load/store operations supported by GCC 11. Each record uses one stable
snapshot, so replacing or clearing a callback cannot race with an in-progress
invocation and the old callable stays alive until that invocation finishes.

Exceptions from all callback types remain isolated from supplier control flow.
The reader now also catches non-standard callback exceptions and prevents any
exception from escaping the reader or writer thread entry point, which would
otherwise call `std::terminate`. An unexpected failure outside the writer's
per-record handling initiates queue shutdown instead of leaving the published
and returned record counts permanently mismatched.

Calls to `next()` are serialized. This makes concurrent consumers safe with
respect to startup, last-record metadata, and next-callback observations.
Lazy startup and `close()` are serialized by a lifecycle mutex, so workers
cannot be started twice or be started after shutdown has begun.
Callback targets themselves can still be invoked concurrently by multiple
writer threads; mutable state captured by a write callback must therefore be
synchronized by the callback author.

`setProcessPropertyLists()` and its getter now use an atomic flag. Each SD
record samples the flag once, avoiding a data race and ensuring internally
consistent handling of all property lists in that record.

## Supplier consistency fixes

### SD supplier

- A truly empty stream is distinguished from a successfully read blank record.
  Empty input is logical EOF; successfully read blank input remains a null
  record, matching `ForwardSDMolSupplier`.
- EOF discovered by read-ahead is not allowed to hide the final real molecule.
- Record assembly appends each line and newline without constructing a
  temporary concatenated string.
- SD property parsing now follows the current `ForwardSDMolSupplier`
  string-view implementation. Trimming no longer allocates a new string for
  every property line.
- Both SD readers now refresh the trimmed line while skipping a property whose
  header has no valid label. Previously that loop retained the first nonblank
  value and ran to EOF instead of stopping at the next blank line.

### SMILES supplier

- Comment and blank-line skipping no longer creates a synthetic null record at
  EOF.
- Title-line and data-line reads update line numbers in the same way as
  `SmilesMolSupplier`, including names derived from line numbers.
- Title columns, unnamed columns, and empty trailing columns retain the same
  property names and values as the single-threaded supplier.
- The completed line is moved into the reader record instead of copied.

## Performance changes

- Queue push accepts values by value and moves them into ring storage; pop
  moves them out.
- Input text moves from the reader to the input queue, then from the writer to
  the output queue, and finally into last-item storage. Large SD records are no
  longer copied at every queue boundary.
- Molecules remain under unique ownership throughout the pipeline.
- The writer-completion mutex and compound counter were replaced by one atomic
  increment per writer thread.
- The writer-thread vector is reserved once and threads are emplaced directly.
- Queue notifications use `notify_one()` for ordinary data flow and
  `notify_all()` only for shutdown.
- SD property trimming uses `std::string_view`.

## Regression coverage

The tests were extended to cover:

- empty and blank SD input compared with `ForwardSDMolSupplier`;
- empty and blank SMILES input compared with `SmilesMolSupplier`;
- SMILES title properties, unnamed-column behavior, and line-derived names;
- invalid SD property headers followed by a valid property in both SD
  suppliers;
- destruction while one-element input and output queues are under load;
- rejection of a zero-capacity queue;
- the existing high-writer-count single-record race;
- both blank-only and comment-only Python SMILES files (the test body was
  corrected so assertions run for each subcase).

Per request, no build or test command was run. Verification in this change is
limited to static code review, formatting, diff inspection, and whitespace
checks.
