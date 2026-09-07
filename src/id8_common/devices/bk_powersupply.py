"""
BK 9130 power supply readout and control
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal

class bk_pid(Device):
    """PID loop on a B&K 9130 triple-output DC power supply.

    ``VAL`` is the setpoint written to the loop, ``CVAL`` its current value.
    Its configs/devices.yml entry -- currently commented out -- uses prefix
    "8idiSoft:BPMpid." Note the trailing dot: the prefix stops mid-record-name,
    so the components below complete it into "8idiSoft:BPMpid.VAL" and ".CVAL".
    """

    VAL = Component(EpicsSignal, "VAL")
    RDBK = Component(EpicsSignalRO, "CVAL")
    



