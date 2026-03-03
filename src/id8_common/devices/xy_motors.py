from ophyd import Device
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt


class XY_Motors(Device):
    """A device class for controlling XY motor pairs.

    This class provides control over pairs of motors that operate in X and Y
    directions. It is used for precise two-dimensional positioning of
    beamline components.
    """

    x = FCpt(EpicsMotor, "{_x_pv}", labels={"motors"})
    y = FCpt(EpicsMotor, "{_y_pv}", labels={"motors"})

    def __init__(
        self,
        prefix: str,
        x_motor: str,
        y_motor: str,
        *args,
        **kwargs,
    ):
        pieces = prefix.strip(":").split(":")
        motor_prefix = ":".join(pieces[:-1])

        # Use prefix if available, otherwise use the direct motor string
        self._x_pv = f"{motor_prefix}:{x_motor}" if motor_prefix else x_motor
        self._y_pv = f"{motor_prefix}:{y_motor}" if motor_prefix else y_motor

        super().__init__(prefix, *args, **kwargs)