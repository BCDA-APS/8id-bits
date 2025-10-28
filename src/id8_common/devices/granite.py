"""
granite 1 DOF Motor
"""

from ophyd import Device
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt


class granite_device(Device):
    """A device class for controlling Granite stages in the beamline.

    This class provides control over Granite stages used for precise positioning
    and motion control. It includes functionality for controlling position,
    velocity, and acceleration parameters.
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
