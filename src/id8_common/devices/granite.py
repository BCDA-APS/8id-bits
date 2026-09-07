"""
granite 1 DOF Motor
"""

from ophyd import Device
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt


class granite_device(Device):
    """The single-axis granite base that carries the 8-ID-I sample environment.

    Only one motor: `x`, which slides the whole environment along the granite.
    plans/set/select_device.py drives it to park the qnw, rheometer or robot at
    the beam. It opens granite_8idi_valve before moving and closes it after;
    the valve exists only to gate this motion.

    Built from configs/devices.yml as `granite` with prefix
    "8idiSoft:CR8-I2:US" and x_motor "m8". __init__ drops the last colon-
    separated piece of the prefix ("US") to get the motor prefix, so `x` ends
    up on "8idiSoft:CR8-I2:m8".
    """

    def __init__(
        self,
        prefix: str,
        x_motor: str,
        *args,
        **kwargs,
    ):
        """Initialize the Granite stage device.

        Args:
            prefix: The EPICS PV prefix for the device
            x_motor: The name of the x-axis motor PV
            **kwargs: Additional keyword arguments passed to the Device constructor
        """
        # Determine the prefix for the motors
        pieces = prefix.strip(":").split(":")
        self.motor_prefix = ":".join(pieces[:-1])

        self._x_neg_motor = x_motor

        super().__init__(prefix, *args, **kwargs)

    x = FCpt(EpicsMotor, "{motor_prefix}:{_x_neg_motor}", labels={"motors"})
