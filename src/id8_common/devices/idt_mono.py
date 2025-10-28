"""
IDT Mono in station 8-ID-A
"""

from ophyd import Device
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt


class IDTMono(Device):
    """A device class for controlling the IDT monochromator in the beamline.

    This class provides control over the IDT monochromator used for X-ray
    beam conditioning. It includes functionality for controlling Bragg angles,
    crystal gaps, flags, and various alignment motors.
    """

    def __init__(
        self,
        prefix: str,
        bragg_motor: str,
        xtal_gap_motor: str,
        flag_motor: str,
        coarse_pitch_motor: str,
        coarse_roll_motor: str,
        x_pitch_motor: str,
        y_pitch_motor: str,
        pitch_piezo: str,
        roll_piezo: str,
        *args,
        **kwargs,
    ):
        """Initialize the IDT monochromator device.

        Args:
            prefix: The EPICS PV prefix for the device
            bragg_motor: The name of the Bragg angle motor PV
            xtal_gap_motor: The name of the crystal gap motor PV
            flag_motor: The name of the flag motor PV
            coarse_pitch_motor: The name of the coarse pitch motor PV
            coarse_roll_motor: The name of the coarse roll motor PV
            x_pitch_motor: The name of the X pitch motor PV
            y_pitch_motor: The name of the Y pitch motor PV
            **kwargs: Additional keyword arguments passed to the Device constructor
        """
        # Determine the prefix for the motors
        pieces = prefix.strip(":").split(":")
        self.motor_prefix = ":".join(pieces[:-1])

        self._bragg_motor = bragg_motor
        self._xtal_gap_motor = xtal_gap_motor
        self._flag_motor = flag_motor
        self._coarse_pitch_motor = coarse_pitch_motor
        self._coarse_roll_motor = coarse_roll_motor
        self._x_pitch_motor = x_pitch_motor
        self._y_pitch_motor = y_pitch_motor
        self._pitch_piezo = pitch_piezo
        self._roll_piezo = roll_piezo

        super().__init__(prefix, *args, **kwargs)

    bragg = FCpt(EpicsMotor, "{motor_prefix}:{_bragg_motor}", labels={"motors"})
    xtal_gap = FCpt(EpicsMotor, "{motor_prefix}:{_xtal_gap_motor}", labels={"motors"})
    flag = FCpt(EpicsMotor, "{motor_prefix}:{_flag_motor}", labels={"motors"})
    coarse_pitch = FCpt(
        EpicsMotor, "{motor_prefix}:{_coarse_pitch_motor}", labels={"motors"}
    )
    coarse_roll = FCpt(
        EpicsMotor, "{motor_prefix}:{_coarse_roll_motor}", labels={"motors"}
    )
    x_pitch = FCpt(EpicsMotor, "{motor_prefix}:{_x_pitch_motor}", labels={"motors"})
    y_pitch = FCpt(EpicsMotor, "{motor_prefix}:{_y_pitch_motor}", labels={"motors"})

    # TODO: These are Labjack signals instead of motors
    pitch_piezo = FCpt(EpicsMotor, "{motor_prefix}:{_pitch_piezo}", labels={"motors"})
    roll_piezo = FCpt(EpicsMotor, "{motor_prefix}:{_roll_piezo}", labels={"motors"})