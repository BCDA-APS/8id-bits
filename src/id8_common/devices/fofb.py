"""
Fast orbit feedback loop status 
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO

class FOFB(Device):
    """Horizontal and vertical FOFB loops"""
    h_loop = Component(EpicsSignalRO, ":H:LoopStatusM", string=False)
    v_loop = Component(EpicsSignalRO, ":V:LoopStatusM", string=False)


