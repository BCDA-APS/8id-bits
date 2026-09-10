"""Flight path device configuration for 8ID-I beamline.

This module defines the flight path device which controls the length of the flight path
using an EPICS motor.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsMotor


class FlightPath(Device):
    """Device representing the flight path length control.

    This device controls the length of the flight path using an EPICS motor.
    """

    length = Component(EpicsMotor, "m1", name="length")
    swing = Component(EpicsMotor, "m2", name="swing")

    ds_x = Component(EpicsMotor, "m5")
    ds_y = Component(EpicsMotor, "m6")
    us = Component(EpicsMotor, "m7")
