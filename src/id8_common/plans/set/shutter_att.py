"""
This module provides plans for controlling the beam shutter and safety interlocks
at the 8ID-E station.
"""

from apsbits.core.instrument_init import oregistry
import epics as pe
import time

shutter_8ide = oregistry["shutter_8ide"]
pd = oregistry["pd"]
filter_beam = oregistry["filter_8ide"]

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

def PIND_status(status=0):
    """Remove the PIND from the beam path."""
    pd.active.put(status)

def att(att_ratio: float=1.0):
    """Set the attenuation ratio with multiple attempts.

    Args:
        att_ratio: Attenuation ratio to set (0-15)
    """
    if att_ratio < 1.0:
        att_ratio = 1.0
    else:
        pass
    filter_beam.attenuation.move(att_ratio)
    time.sleep(0.5)
