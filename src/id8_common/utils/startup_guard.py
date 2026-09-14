"""
Block EPICS writes attempted while a startup script is running.

Starting a session must never disturb hardware. It used to: with
``ALLOW_AREA_DETECTOR_WARMUP: True``, ``ad_setup()`` called apstools'
``AD_prime_plugin2()``, which fires a real exposure -- and
``AD_plugin_primed()`` never reports the Eiger as primed (its cam and HDF1
``data_type`` differ permanently), so *every* session start fired one. That
crashed measurement G0209 on 2026-09-12. The warmup branch is gone now; this
module is the backstop that keeps it, or anything like it, from coming back
in library code where a source grep would not find it.

Armed at the top of a startup script and reported at the bottom::

    from id8_common.utils.startup_guard import arm_startup_guard
    from id8_common.utils.startup_guard import report_startup_writes

    arm_startup_guard()
    ...
    report_startup_writes()

A blocked write prints a red banner naming the signal and the line
responsible, records it, and returns without writing. It never raises: a
guard that aborted startup could lock everyone out of the beamline over one
false positive.

Two implementation details are the whole point of this module -- a naive
guard is silently useless:

* Both ``put`` and ``set`` are defined on ``EpicsSignal``, *not* on
  ``EpicsSignalBase`` (checked against the ophyd in use). Patching the base
  class intercepts nothing, because the subclass method wins.
* ``.set()`` is not optional. ``BlueskyInterface.stage()`` writes through
  ``.set()``, and ``AD_prime_plugin2`` -- the exact call that crashed G0209 --
  uses ``sig.set(val).wait()``. Patching only ``put`` would miss it entirely.

``EpicsSignalRO`` needs no patch; its ``put``/``set`` already raise.
``EpicsSignalWO`` (``devices/epics_signal_wo.py``) overrides ``set`` but
reaches EPICS through ``self.put``, so the ``put`` patch covers it.
"""

import logging
import os
import sys

from ophyd.signal import EpicsSignal
from ophyd.status import Status

logger = logging.getLogger(__name__)

RED = "\033[91m"
RESET = "\033[0m"

# Exact-match pvnames that may be written during startup.
#
# 8ideSoft:Reg1 is expt.measurement_num -- a soft register, and expt only
# writes it when the counter has gone backwards relative to the on-disk
# mirror. That write is a data-safety correction (a counter that restarts
# silently overwrites data), so it is allowed through rather than suppressed.
DEFAULT_ALLOW = ("8ideSoft:Reg1",)

_state = None


class _GuardState:
    """Everything the patch needs, kept off the module namespace."""

    def __init__(self, allow, startup_file):
        self.allow = frozenset(allow)
        self.startup_file = startup_file
        self.blocked = []  # (pvname, "file:line", method)
        self.orig_put = EpicsSignal.put
        self.orig_set = EpicsSignal.set


def _in_startup_call(startup_file):
    """Is the startup module that armed the guard still on this stack?

    The check exists for the case where startup raises somewhere between
    arming and reporting and leaves the patch installed. Writes typed at the
    interactive prompt afterwards must still reach hardware -- a guard that
    silently ate someone's ``att(20)`` would be worse than the bug it
    prevents. It also means EPICS monitor-callback threads, which have no
    startup frame, are never blocked.
    """
    frame = sys._getframe(1)
    while frame is not None:
        if frame.f_code.co_filename == startup_file:
            return True
        frame = frame.f_back
    return False


def _blame():
    """First stack frame outside this module and outside ophyd itself."""
    frame = sys._getframe(1)
    while frame is not None:
        fname = frame.f_code.co_filename
        if fname != __file__ and f"{os.sep}ophyd{os.sep}" not in fname:
            return f"{fname}:{frame.f_lineno} in {frame.f_code.co_name}"
        frame = frame.f_back
    return "<unknown caller>"


def _block(signal, value, method):
    """Record and announce a refused write. Never raises."""
    pvname = getattr(signal, "pvname", getattr(signal, "name", "<unknown>"))
    where = _blame()
    _state.blocked.append((pvname, where, method))
    print(
        f"{RED}*** startup guard: BLOCKED {pvname}.{method}({value!r}) ***\n"
        f"    from {where}\n"
        f"    Startup must not touch hardware. The write did NOT happen; the\n"
        f"    session is continuing. If this write is legitimate, add the\n"
        f"    pvname to DEFAULT_ALLOW in utils/startup_guard.py.{RESET}"
    )
    logger.warning("startup guard blocked write to %s from %s", pvname, where)


def arm_startup_guard(allow=DEFAULT_ALLOW):
    """Patch ``EpicsSignal.put``/``.set`` to refuse writes during startup.

    Call from a startup script, before any device is built. ``allow`` is an
    exact-match pvname allowlist.
    """
    global _state
    if _state is not None:
        print(f"{RED}*** startup guard already armed; ignoring second arm ***{RESET}")
        return

    # The caller is the startup module; its filename is what the stack-origin
    # check below looks for. Taking it from the caller rather than hardcoding
    # startup.py/startup_ophyd.py keeps this usable from any entry point.
    startup_file = sys._getframe(1).f_code.co_filename
    _state = _GuardState(allow, startup_file)

    orig_put = _state.orig_put
    orig_set = _state.orig_set

    def guarded_put(self, value, *args, **kwargs):
        if _state is not None and getattr(self, "pvname", None) not in _state.allow:
            if _in_startup_call(_state.startup_file):
                _block(self, value, "put")
                return None
        return orig_put(self, value, *args, **kwargs)

    def guarded_set(self, value, *args, **kwargs):
        if _state is not None and getattr(self, "pvname", None) not in _state.allow:
            if _in_startup_call(_state.startup_file):
                _block(self, value, "set")
                # An already-finished Status, so a caller doing
                # `sig.set(v).wait()` -- AD_prime_plugin2 does exactly that --
                # returns instead of hanging until its timeout.
                st = Status(self)
                st.set_finished()
                st.wait()
                return st
        return orig_set(self, value, *args, **kwargs)

    EpicsSignal.put = guarded_put
    EpicsSignal.set = guarded_set
    print("[startup_guard] armed: EPICS writes during startup will be blocked")


def report_startup_writes():
    """Restore ophyd's own put/set and report anything that was blocked."""
    global _state
    if _state is None:
        return []

    EpicsSignal.put = _state.orig_put
    EpicsSignal.set = _state.orig_set
    blocked = list(_state.blocked)
    _state = None

    if not blocked:
        print("[startup_guard] disarmed: no EPICS writes attempted during startup")
        return blocked

    lines = "\n".join(f"    {pv}.{method}()  <-  {where}" for pv, where, method in blocked)
    print(
        f"{RED}\n*** startup guard: {len(blocked)} EPICS write(s) BLOCKED during "
        f"startup ***\n{lines}\n"
        f"    None of them happened. Hardware was not touched, but something in\n"
        f"    the startup path wanted to write -- please report this.\n{RESET}"
    )
    return blocked
