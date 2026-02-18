"""
8ID undulator readout.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal

class Undulator(Device):
    """Device representing a 8ID undulator."""

    upstream = Component(EpicsSignalRO, "USID:EnergyM")
    downstream = Component(EpicsSignalRO, "DSID:EnergyM")
    


