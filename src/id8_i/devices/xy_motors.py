"""
Damm (dynamic aperture) in station 8-ID-D
"""

from ophyd import Device
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt


class XY_Motors(Device):
    def __init__(
        self,
        prefix: str,
        x_motor: str,
        y_motor: str,
        *args,
        **kwargs,
    ):
        # Determine the prefix for the motors
        pieces = prefix.strip(":").split(":")
        self.motor_prefix = ":".join(pieces[:-1])

        self._x_motor = x_motor
        self._y_motor = y_motor

        super().__init__(prefix, *args, **kwargs)

    # Real motors that directly control the slits
    x = FCpt(EpicsMotor, "{motor_prefix}:{_x_motor}", labels={"motors"})
    y = FCpt(EpicsMotor, "{motor_prefix}:{_y_motor}", labels={"motors"})
