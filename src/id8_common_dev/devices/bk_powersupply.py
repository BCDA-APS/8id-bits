"""
BK 9130 power supply readout and control
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal

class bk_pid(Device):
    """Device representing a B&K 9130 class device.
    Component(EpicsSignal, "")
    """

    VAL = Component(EpicsSignal, "VAL")
    RDBK = Component(EpicsSignalRO, "CVAL")
    



