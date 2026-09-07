"""LabJack LJT705 in 8-ID-I."""

from id8_common.devices.epics_signal_wo import EpicsSignalWO
from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal


class LabJack(Device):
    """The two digital outputs of the 8-ID-I LabJack LJT705 used for the shutter.

    Built from configs/devices.yml as `labjack`, prefix "8idiSoft:LJT705:", so
    `operation` is 8idiSoft:LJT705:Bo0 and `logic` is :Bo1. Only those two
    digital-output bits are declared -- the LJT705 has many more analog and
    digital channels; add them here if you need one.

    Nothing in id8_common reads or writes this device today; the 8-ID-E shutter
    goes through fast_shutter.FastShutter instead.

    `operation` is an EpicsSignalWO -- a write-only signal whose .set() reports
    success without verifying the readback. See epics_signal_wo.py for why.
    """

    operation = Component(EpicsSignalWO, "Bo0", kind="omitted")
    logic = Component(EpicsSignal, "Bo1")
