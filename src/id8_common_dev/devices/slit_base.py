"""
Slit Base Devices at 8-ID
"""

from ophyd import Device
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt


class SlitBase(Device):
    """A base device class for controlling slit devices.

    This class provides basic functionality for controlling slit devices,
    including horizontal and vertical motor control. It serves as a base
    class for more specific slit device implementations.
    """

    def __init__(
        self,
        prefix: str,
        h_motor: str,
        v_motor: str,
        *args,
        **kwargs,
    ):
        """Initialize the slit base device.

        Args:
            prefix: The EPICS PV prefix for the device
            h_motor: The name of the horizontal motor PV
            v_motor: The name of the vertical motor PV
            *args: Additional positional arguments passed to the Device constructor
            **kwargs: Additional keyword arguments passed to the Device constructor
        """
        # Determine the prefix for the motors
        pieces = prefix.strip(":").split(":")
        self.motor_prefix = ":".join(pieces[:-1])

        self._h_motor = h_motor
        self._v_motor = v_motor

        super().__init__(prefix, *args, **kwargs)

    # Real motors that directly control the slits
    h = FCpt(EpicsMotor, "{motor_prefix}:{_h_motor}", labels={"motors"})
    v = FCpt(EpicsMotor, "{motor_prefix}:{_v_motor}", labels={"motors"})
