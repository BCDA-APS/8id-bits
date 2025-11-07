"""
Shutter control logic for the 8ID-E beamline.

This module provides plans for controlling the beam shutter and safety interlocks
at the 8ID-E station.
"""

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

shutter_8ide = oregistry["shutter_8ide"]


def showbeam():
    """Open the beam shutter to show the beam."""
    yield from bps.mv(shutter_8ide.operation, 0)


def blockbeam():
    """Block the beam by closing the shutter."""
    yield from bps.mv(shutter_8ide.operation, 1)


def shutteron():
    """Enable the shutter control logic."""
    yield from bps.mv(shutter_8ide.logic, 0)


def shutteroff():
    """Disable the shutter control logic."""
    yield from bps.mv(shutter_8ide.logic, 1)


def post_align():
    """Configure system for post-alignment state by blocking the beam."""
    yield from blockbeam()


def pre_align():
    """Configure system for pre-alignment state by disabling shutter control."""
    yield from shutteroff()
