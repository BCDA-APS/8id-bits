"""
HHL Slits in station 8-ID-A
"""

# from apstools.devices.hhl_slits import HHLSlits

from ophyd import Device
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt


class HHLSlits(Device):
    """
    High Heat Load SHHLSlitst.

    There are no independent parts to move, so each axis only has center and size.

    Based on the 25-ID-A whitebeam slits.

    The motor parameters listed below specify which motor records
    control which axis. The last piece of the PV prefix will be
    removed, and the motor number added on. For example, if the prefix
    is "255ida:slits:US:", and the pitch motor is "255ida:slits:m3",
    then *pitch_motor* should be "m3".

    Parameters
    ==========
    prefix:
      EPICS prefix required to communicate with HHL Slit IOC, ex: "25ida:slits:US:"
    pitch_motor:
      The motor record suffix controlling the real pitch motor, ex "m3"
    yaw_motor:
      The motor record suffix controlling the real yaw motor, ex "m4"
    horizontal_motor:
      The motor record suffix controlling the real horizontal motor, ex: "m1"
    diagonal_motor:
      The motor record suffix controlling the real diagonal motor, ex: "m2"
    """

    def __init__(
        self,
        prefix: str,
        secondary_prefix: str,
        pitch_motor: str,
        yaw_motor: str,
        horizontal_motor: str,
        diagonal_motor: str,
        hgap_motor: str,
        hcen_motor: str,
        vgap_motor: str,
        vcen_motor: str,
        *args,
        **kwargs,
    ):
        """Initialize the HHL slits device.

        Args:
            prefix: The EPICS PV prefix for the device
            secondary_prefix: The secondary EPICS PV prefix for additional controls
            pitch_motor: The name of the pitch motor PV
            yaw_motor: The name of the yaw motor PV
            horizontal_motor: The name of the horizontal motor PV
            diagonal_motor: The name of the diagonal motor PV
            hgap_motor: The name of the horizontal gap motor PV
            hcen_motor: The name of the horizontal center motor PV
            vgap_motor: The name of the vertical gap motor PV
            vcen_motor: The name of the vertical center motor PV
            **kwargs: Additional keyword arguments passed to the Device constructor
        """
        # Determine the prefix for the motors
        pieces = prefix.strip(":").split(":")
        self.motor_prefix = ":".join(pieces[:-1])

        pieces = secondary_prefix.strip(":").split(":")
        self.secondary_motor_prefix = ":".join(pieces[:-1])

        self._pitch_motor = pitch_motor
        self._yaw_motor = yaw_motor
        self._horizontal_motor = horizontal_motor
        self._diagonal_motor = diagonal_motor

        self._hgap_motor = hgap_motor
        self._hcen_motor = hcen_motor
        self._vgap_motor = vgap_motor
        self._vcen_motor = vcen_motor

        super().__init__(prefix, *args, **kwargs)

    # Real motors that directly control the slits
    pitch = FCpt(EpicsMotor, "{motor_prefix}:{_pitch_motor}", labels={"motors"})
    yaw = FCpt(EpicsMotor, "{motor_prefix}:{_yaw_motor}", labels={"motors"})
    horizontal = FCpt(
        EpicsMotor, "{motor_prefix}:{_horizontal_motor}", labels={"motors"}
    )
    diagonal = FCpt(EpicsMotor, "{motor_prefix}:{_diagonal_motor}", labels={"motors"})

    hgap = FCpt(EpicsMotor, "{secondary_motor_prefix}:{_hgap_motor}", labels={"motors"})
    hcen = FCpt(EpicsMotor, "{secondary_motor_prefix}:{_hcen_motor}", labels={"motors"})
    vgap = FCpt(EpicsMotor, "{secondary_motor_prefix}:{_vgap_motor}", labels={"motors"})
    vcen = FCpt(EpicsMotor, "{secondary_motor_prefix}:{_vcen_motor}", labels={"motors"})
