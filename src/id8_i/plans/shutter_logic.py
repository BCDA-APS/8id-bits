"""
This module provides plans for controlling the beam shutter and safety interlocks
at the 8ID-E station.
"""

from apsbits.core.instrument_init import oregistry
import epics as pe
import time

shutter_8ide = oregistry["shutter_8ide"]


def showbeam():
    """Open the beam shutter to show the beam."""
    shutter_8ide.operation.put(0)
    time.sleep(0.5)


def blockbeam():
    """Block the beam by closing the shutter."""
    shutter_8ide.operation.put(1)
    time.sleep(0.5)


def shutteron():
    """Enable the shutter control logic."""
    shutter_8ide.logic.put(0)
    time.sleep(0.5)


def shutteroff():
    """Disable the shutter control logic."""
    shutter_8ide.logic.put(1)
    time.sleep(0.5)


def post_align():
    """Configure system for post-alignment state by blocking the beam."""
    pe.caput("8idiSoft:FLIGHT:bo1:8", 1)
    blockbeam()
    time.sleep(0.5)


def pre_align():
    """Configure system for pre-alignment state by disabling shutter control."""
    pe.caput("8idiSoft:FLIGHT:bo1:8", 0)
    shutteroff()
    time.sleep(0.5)
