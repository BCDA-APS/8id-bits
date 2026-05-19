"""
PIND power supply readout and control
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal

class PIND(Device):
    """PIN diode upstream of 8IDE Huber"""

    state = Component(EpicsSignalRO, "status2.RVAL")
    active = Component(EpicsSignal, "set2")
    
    



