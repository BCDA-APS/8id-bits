"""
This module provides plans for controlling the beam shutter and safety interlocks
at the 8ID-E station.
"""

import epics as pe
import time
from id8_common.registry import get_connected_device
from id8_common.registry import oregistry

# softglue is present in every session, so it stays resolved at import time.
# The shutter, PIND and filter set can legitimately be absent (offline IOC,
# other station), so they are resolved per call by get_connected_device():
# a missing one then raises naming itself, instead of binding None here and
# surfacing later as an AttributeError somewhere unrelated.
softglue = oregistry.get("softglue")

def showbeam():
    """Open the beam shutter to show the beam."""
    get_connected_device("shutter_8ide").operation.put(0)
    time.sleep(0.05)

def blockbeam():
    """Block the beam by closing the shutter."""
    get_connected_device("shutter_8ide").operation.put(1)
    time.sleep(0.05)

def shutteron():
    """Enable the shutter control logic."""
    get_connected_device("shutter_8ide").logic.put(0)
    time.sleep(0.5)

def shutteroff():
    """Disable the shutter control logic."""
    get_connected_device("shutter_8ide").logic.put(1)
    time.sleep(0.5)

def post_align():
    """Configure system for post-alignment state by blocking the beam."""
    pe.caput("8idiSoft:FLIGHT:bo1:8", 1)
    blockbeam()
    time.sleep(0.5)

def pre_align():
    """Configure system for pre-alignment state by disabling shutter control."""
    pe.caput("8idiSoft:FLIGHT:bo1:8", 0)
    softglue.enable_rigaku.put('0')
    shutteroff()
    time.sleep(0.5)

def PIND_status(status=0):
    """Move the PIN diode in or out of the beam path.

    Args:
        status: 0 retracts the diode (the default -- out of the beam),
            1 inserts it.
    """
    get_connected_device("pd").active.put(status)

#: Seconds to allow the filter set to reach a requested attenuation.
FILTER_SETTLE_TIMEOUT = 30.0

#: Fractional disagreement between requested and achieved attenuation that is
#: worth complaining about. The set is discrete, so exact hits are not expected
#: -- att(20) has legitimately landed on 22. A stuck blade looks nothing like
#: that: on 2026-09-11 att(10) achieved 1.78 and ran a measurement at the wrong
#: attenuation without saying so.
FILTER_TOLERANCE = 0.25


def _wait_for_filters(atten, timeout):
    """Block until the filter set has stopped moving, and return the readback.

    Does not use ``PVPositioner.move()``. ``filterBusy`` is an enum whose Done
    state is 0, but the positioner classes take ophyd's ``done_value=1``
    default, so the positioner reads "Changing" as "arrived" and ``move()``
    returns the moment the blades start. Polling the same signal directly
    avoids depending on that, fixed or not.

    A request that changes nothing never raises filterBusy at all, so the wait
    for it to go busy is a short grace period rather than a condition.
    """
    busy = atten.done
    readback = atten.readback
    poll = 0.1

    grace = time.monotonic() + 1.0
    while time.monotonic() < grace:
        if busy.get() == 1:
            break
        time.sleep(poll)

    # Done, and the readback unchanged over three consecutive polls -- the IOC
    # can report Done between blades of a multi-blade change.
    deadline = time.monotonic() + timeout
    last, stable = None, 0
    while time.monotonic() < deadline:
        if busy.get() == 0:
            now = readback.get()
            stable = stable + 1 if now == last else 0
            last = now
            if stable >= 3:
                return now
        else:
            last, stable = None, 0
        time.sleep(poll)

    raise TimeoutError(
        f"filter set still moving {timeout:g}s after the request "
        f"(attenuation_actual={readback.get()})"
    )


def att(att_ratio: float = 1.0, *, strict: bool = False, timeout: float = FILTER_SETTLE_TIMEOUT):
    """Set the attenuation ratio on the 8-ID-E filter set, and check it took.

    Waits for the blades to settle and compares what was achieved against what
    was asked for. The old version wrote the setpoint and slept 0.5 s, so a
    caller could not tell a satisfied request from an impossible one.

    Args:
        att_ratio: Attenuation factor, i.e. how many times the beam is
            weakened. Anything below 1.0 is clamped to 1.0, since a filter
            cannot make the beam brighter than no filter at all.
        strict: Raise instead of warning when the achieved attenuation misses
            the request by more than FILTER_TOLERANCE. Worth setting for an
            attenuation series, where the attenuation is the variable and a
            silent miss makes the dataset meaningless.
        timeout: Seconds to wait for the blades to settle.

    Returns:
        float: the achieved attenuation, read back from the hardware.
    """
    if att_ratio < 1.0:
        att_ratio = 1.0

    atten = get_connected_device("filter_8ide").attenuation
    atten.setpoint.put(att_ratio)
    achieved = _wait_for_filters(atten, timeout)

    if abs(achieved - att_ratio) > FILTER_TOLERANCE * att_ratio:
        msg = (
            f"attenuation: asked for {att_ratio:g}, achieved {achieved:g}. "
            f"Check the filter set -- a stuck blade reads exactly like this."
        )
        if strict:
            raise RuntimeError(msg)
        print(f"WARNING: {msg}")

    return achieved
