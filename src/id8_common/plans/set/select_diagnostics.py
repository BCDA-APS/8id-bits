"""
Plan that allows for moving to pre-programmed positions seen as strings.

This module provides plans for moving the upstream diagnostics 
among microscope, wavefront sensor, and pind.
"""
# from typing import Literal
import time
from apsbits.core.instrument_init import oregistry

cam_stage_8idi = oregistry["cam_stage_8idi"]

def select_diagnostics(name: str):

    if name == "microscope":
        cam_stage_8idi.x.move(-0.105)
        cam_stage_8idi.y.move(-0.184)
    if name == "pind":
        cam_stage_8idi.x.move(-205.0)
        cam_stage_8idi.y.move(-0.184)
    if name == "wavefront":
        cam_stage_8idi.x.move(-114.75)
        cam_stage_8idi.y.move(-0.184)
    else:
        print("Diagnostic names must be microscope, pind or wavefront")

    time.sleep(2.0)

