"""
Shutter control logic for the 8ID-I beamline.

This module provides plans for controlling the beam shutter and safety interlocks
using the LabJack device.
"""

import epics as pe
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

labjack = oregistry["labjack"]


def showbeam():
    """Open the beam shutter to show the beam."""
    if labjack.operation.get(use_cache=False) != 0:
        yield from bps.mv(labjack.operation, 0)


def blockbeam():
    """Block the beam by closing the shutter."""
    if labjack.operation.get(use_cache=False) != 1:
        yield from bps.mv(labjack.operation, 1)


def shutteron():
    """Enable the shutter control logic."""
    yield from bps.mv(labjack.logic, 0)


def shutteroff():
    """Disable the shutter control logic."""
    yield from bps.mv(labjack.logic, 1)


def post_align():
    """Configure system for post-alignment state.

    Sets flight tube output and blocks the beam.
    """
    pe.caput("8idiSoft:FLIGHT:bo1:8", 1)
    yield from blockbeam()


def pre_align():
    """Configure system for pre-alignment state.

    Clears flight tube output and disables shutter control.
    """
    pe.caput("8idiSoft:FLIGHT:bo1:8", 0)
    yield from shutteroff()
