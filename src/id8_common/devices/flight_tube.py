# """
# SAXS Flight Tube
# """

from ophyd import Component
from ophyd import Device
from ophyd import EpicsMotor


class FlightTubeDetector(Device):
    """Device representing the detector position in the flight tube.

    This device controls the vertical (z) and angular (tth) position of the detector
    in the SAXS flight tube.
    """

    z = Component(EpicsMotor, "m1")
    tth = Component(EpicsMotor, "m2")


class FlightTubeBeamStop(Device):
    """Device representing the beam stop position in the flight tube.

    This device controls the position of the beam stop in the SAXS flight tube,
    including downstream (ds) and upstream (us) positions.
    """

    ds_x = Component(EpicsMotor, "m5")
    ds_y = Component(EpicsMotor, "m6")
    us = Component(EpicsMotor, "m7")
