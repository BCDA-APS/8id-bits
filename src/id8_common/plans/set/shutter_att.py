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

def att(att_ratio: float=1.0):
    """Set the attenuation ratio on the 8-ID-E filter set.

    Args:
        att_ratio: Attenuation factor, i.e. how many times the beam is
            weakened. Anything below 1.0 is clamped to 1.0, since a filter
            cannot make the beam brighter than no filter at all.
    """
    if att_ratio < 1.0:
        att_ratio = 1.0
    get_connected_device("filter_8ide").attenuation.move(att_ratio)
    time.sleep(0.5)
