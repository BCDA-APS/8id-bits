"""
Lakeshore 336 (temperature readout only)
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO


class Lakeshore(Device):
    """Device representing a Lakeshore 336 temperature controller.

    This device provides read-only access to temperature measurements from up to four
    input channels of a Lakeshore 336 temperature controller.
    """

    readback_ch1 = Component(EpicsSignalRO, "IN1")
    readback_ch2 = Component(EpicsSignalRO, "IN2")
    readback_ch3 = Component(EpicsSignalRO, "IN3")
    readback_ch4 = Component(EpicsSignalRO, "IN4")
