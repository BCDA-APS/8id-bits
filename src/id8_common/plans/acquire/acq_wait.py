"""Bounded, fault-aware waiting for detector acquisitions.

One polling loop behind the "wait until the hardware finishes" steps in the
acquisition path -- used by eiger4m_modes, lambda2m_modes and
dual_acq_eiger4m_rigaku3m. It exists because the same loop had been written
five times with five different levels of care, and the differences were real
bugs:

* ``rigaku3m_modes._wait_for_state`` was the good one -- bounded, and it raised
  the moment the detector reached a terminal state. It is the model this module
  generalises, though rigaku3m_modes still runs its own bounded copy rather than
  importing from here.
* ``dual_acq``'s ``cam_busy`` for the same detector tested ``detector_state != 0``,
  which is true for Error/Disconnected/Aborted, so a dead detector counted as
  *busy* until a multi-minute timeout expired -- with the shutter open.
* ``acquire_eiger_external`` and ``acquire_lambda_external`` polled 600 times and
  then simply carried on, force-closing a short HDF file. The caller went on to
  write NeXus metadata and submit a DM job for a dataset it had no reason to
  believe was complete.

Porting is only partial. dual_acq_eiger4m_rigaku3m routes every wait through
here; eiger4m_modes and lambda2m_modes use it only for the HDF drain at the end
of their External modes. Their cam waits, and the Internal modes end to end,
are still hand-rolled ``while ... get() == 1: sleep()`` loops with no timeout
and no fault check.

Design rules, each learned from one of those:

* ``timeout`` is keyword-only with no default, so no call site can go unbounded
  by forgetting it.
* Running out of time RAISES. There is deliberately no "warn and continue" mode:
  the only thing it would buy is keeping the truncated-file path reachable.
* A terminal fault raises immediately rather than waiting out the timeout.
* ``done`` is checked BEFORE ``faults`` -- a detector that finished and then
  reported an error has still finished.

Imports only the standard library, so it stays usable from the Ophyd-only
session and from the detector-mode modules that ``acq_helpers`` may not import.
"""

import time as ttime

POLL_INTERVAL = 0.1

# ADCore ADStatus enum, as reported on <cam>:DetectorState_RBV.
STATE_IDLE = 0
STATE_ACQUIRE = 1
STATE_READOUT = 2

# Terminal: these never become Idle on their own, so a loop that only tests
# "not Idle" spins until its timeout.
STATE_FAILED = {6: "Error", 9: "Disconnected", 10: "Aborted"}


class DetectorWaitError(RuntimeError):
    """Base for every failure raised out of wait_until()."""


class DetectorTimeout(DetectorWaitError):
    """The condition never became true within the timeout."""


class DetectorFault(DetectorWaitError):
    """The hardware reached a terminal state while we were waiting."""


def _text(value):
    """Return ``value``, calling it first if it is a zero-arg callable.

    This is what lets ``what=`` be either a plain string or a function that
    builds the message only when an error is actually about to be raised.
    """
    if callable(value):
        return value()

    return value


def wait_until(done, *, timeout, what, faults=(), poll=POLL_INTERVAL, settle=0.0):
    """Poll ``done()`` until it is true. Raise on a terminal fault or on timeout.

    Args:
        done: zero-arg predicate; the wait ends when it returns true.
        timeout: seconds. Keyword-only and required, on purpose.
        what: description for the error message. May be a string or a zero-arg
            callable, so an expensive message (a frame count, say) is only built
            if it is actually needed.
        faults: zero-arg callables returning a reason string when the hardware
            has failed, or None when it is healthy.
        poll: seconds between polls.
        settle: slept once before the first poll, for hardware that needs a
            moment to report that it started.

    Returns:
        float: seconds actually waited.

    Example -- the real call in ``eiger4m_modes.acquire_eiger_external()``,
    waiting for the HDF plugin to flush the last frames to disk::

        wait_until(
            lambda: hdf_frames_written(eiger4M.hdf1)[0] >= hdf_frames_written(eiger4M.hdf1)[1],
            timeout=HDF_DRAIN_TIMEOUT,
            what=hdf_progress(eiger4M.hdf1, "eiger4M"),
            faults=(cam_fault(eiger4M.cam, "eiger4M"),),
        )

    In words: keep polling until the plugin has written at least as many frames
    as it was told to expect; give up after HDF_DRAIN_TIMEOUT seconds; and if it
    does give up -- or if the camera reaches Error/Disconnected/Aborted first --
    raise, saying how many frames of how many were written.

    Every argument here is a FUNCTION rather than a value, because each one has
    to be re-evaluated on every pass round the loop. ``lambda: ...`` is just a
    small function written inline; the ``done`` argument above is the same thing
    as::

        def frames_all_written():
            captured, expected = hdf_frames_written(eiger4M.hdf1)
            return captured >= expected

        wait_until(frames_all_written, ...)

    Dropping the ``lambda`` and passing the comparison directly would read the
    frame counters ONCE, before the wait even starts, and hand this loop a fixed
    True or False. ``cam_fault(...)`` and ``hdf_progress(...)`` are the same
    idea: they do not check the camera or build a message when you call them,
    they return a function that does so later, once per poll.
    """
    if settle:
        ttime.sleep(settle)

    start = ttime.time()
    deadline = start + timeout

    while True:
        if done():
            return ttime.time() - start

        for fault in faults:
            reason = fault()
            if reason:
                raise DetectorFault(f"{_text(what)}: {reason}")

        if ttime.time() > deadline:
            raise DetectorTimeout(f"{_text(what)}: timed out after {timeout:g}s")

        ttime.sleep(poll)


def cam_fault(cam, label=""):
    """faults= entry: reports when a camera has reached a terminal state."""
    prefix = f"{label} " if label else ""

    def check():
        try:
            state = int(cam.detector_state.get())
        except Exception:
            return None

        if state in STATE_FAILED:
            message = ""
            try:
                message = str(cam.status_message.get() or "").strip()
            except Exception:
                pass
            return f"{prefix}detector went to {STATE_FAILED[state]}" + (f" ({message})" if message else "")

        return None

    return check


def hdf_frames_written(hdf):
    """(captured, expected) for an HDF plugin, as ints."""
    return int(hdf.num_captured.get()), int(hdf.num_capture.get())


def hdf_progress(hdf, label=""):
    """what= callable reporting how far an HDF plugin got. Cheap enough to poll."""
    prefix = f"{label} " if label else ""

    def describe():
        try:
            captured, expected = hdf_frames_written(hdf)
            return f"{prefix}hdf1 wrote {captured} of {expected} frames"
        except Exception:
            return f"{prefix}hdf1"

    return describe


def acq_timeout(acq_period, num_frames, factor=3.0, floor=60.0):
    """A generous but finite allowance for an acquisition of this size."""
    return max(floor, factor * float(acq_period) * int(num_frames) + floor)
