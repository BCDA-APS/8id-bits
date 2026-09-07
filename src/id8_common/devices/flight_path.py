"""Flight path device configuration for 8ID-I beamline.

This module defines the flight path device which controls the length of the flight path
using an EPICS motor.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsMotor


class FlightPath(Device):
    """The 8-ID-I evacuated flight path between sample and detector.

    Built from configs/devices.yml as `flight_path_8idi`, prefix
    "8idiSoft:FLIGHT:", so `length` is 8idiSoft:FLIGHT:m1 and so on.

    `length` drives the flight-path length and `swing` the detector swing
    angle. Treat `swing` with care: it is shared by eiger4M and rigaku3M,
    and plans/set/select_device.py deliberately never drives it -- only an
    explicit move_detector_axes() call does.

    ds_x/ds_y/us are the beam-stop motors, the same three EPICS motors that
    flight_tube.FlightTubeBeamStop also exposes.
    """

    length = Component(EpicsMotor, "m1", name="length")
    swing = Component(EpicsMotor, "m2", name="swing")

    ds_x = Component(EpicsMotor, "m5")
    ds_y = Component(EpicsMotor, "m6")
    us = Component(EpicsMotor, "m7")
