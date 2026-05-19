"""Module for controlling transfocator devices in the beamline.

This module provides classes and utilities for controlling transfocator
devices, which are used for X-ray beam focusing. It includes functionality
for lens control, positioning, and alignment.
"""

from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt
from ophyd import MotorBundle


class Transfocator(MotorBundle):
    """
    Transfocator (compound refractive lens, CRL).

    Example::

        rl1 = Transfocator(
            "8iddSoft:TRANS:",
            name="rl1",
            pv_y="m1",
            pv_x="m2",
            pv_yaw="m3",
            pv_pitch="m4",
            pv_lens1="m5",
            pv_lens2="m6",
            pv_lens3="m7",
            pv_lens4="m8",
            pv_lens5="m9",
            pv_lens6="m10",
            pv_lens7="m11",
            pv_lens8="m12",
            pv_lens9="m13",
            pv_lens10="m14",
        )

    .. note:: POLAR instrument has some related functions (`transfocator.py`).
        https://github.com/APS-4ID-POLAR/polar_instrument/pull/7/files
    """

    # Motors that position the CRL assembly.
    # By some convention, the names x&y are used rather than h&v.
    x = FCpt(EpicsMotor, "{prefix}{_pv_x}", labels={"motors"})
    y = FCpt(EpicsMotor, "{prefix}{_pv_y}", labels={"motors"})
    pitch = FCpt(EpicsMotor, "{prefix}{_pv_pitch}", labels={"motors"})
    yaw = FCpt(EpicsMotor, "{prefix}{_pv_yaw}", labels={"motors"})

    # Motors that insert or remove each lens.
    # Note: Operated as two-position motors.  Consider refactoring as such.
    lens1 = FCpt(EpicsMotor, "{prefix}{_pv_lens1}", labels={"lens"})
    lens2 = FCpt(EpicsMotor, "{prefix}{_pv_lens2}", labels={"lens"})
    lens3 = FCpt(EpicsMotor, "{prefix}{_pv_lens3}", labels={"lens"})
    lens4 = FCpt(EpicsMotor, "{prefix}{_pv_lens4}", labels={"lens"})
    lens5 = FCpt(EpicsMotor, "{prefix}{_pv_lens5}", labels={"lens"})
    lens6 = FCpt(EpicsMotor, "{prefix}{_pv_lens6}", labels={"lens"})
    lens7 = FCpt(EpicsMotor, "{prefix}{_pv_lens7}", labels={"lens"})
    lens8 = FCpt(EpicsMotor, "{prefix}{_pv_lens8}", labels={"lens"})
    lens9 = FCpt(EpicsMotor, "{prefix}{_pv_lens9}", labels={"lens"})
    lens10 = FCpt(EpicsMotor, "{prefix}{_pv_lens10}", labels={"lens"})

    # TODO: Describe parameters of each lens.
    # TODO: Add method(s) to select & report focal parameters based on lens combination.

    def __init__(
        self,
        prefix: str,
        *args,
        pv_x: str = "",
        pv_y: str = "",
        pv_yaw: str = "",
        pv_pitch: str = "",
        pv_lens1: str = "",
        pv_lens2: str = "",
        pv_lens3: str = "",
        pv_lens4: str = "",
        pv_lens5: str = "",
        pv_lens6: str = "",
        pv_lens7: str = "",
        pv_lens8: str = "",
        pv_lens9: str = "",
        pv_lens10: str = "",
        **kwargs,
    ):
        """Initialize the Transfocator device.

        Args:
            prefix: The EPICS PV prefix for the device
            pv_y: The PV name for Y position motor
            pv_x: The PV name for X position motor
            pv_yaw: The PV name for yaw motor
            pv_pitch: The PV name for pitch motor
            pv_lens1: The PV name for lens 1 motor
            pv_lens2: The PV name for lens 2 motor
            pv_lens3: The PV name for lens 3 motor
            pv_lens4: The PV name for lens 4 motor
            pv_lens5: The PV name for lens 5 motor
            pv_lens6: The PV name for lens 6 motor
            pv_lens7: The PV name for lens 7 motor
            pv_lens8: The PV name for lens 8 motor
            pv_lens9: The PV name for lens 9 motor
            pv_lens10: The PV name for lens 10 motor
            **kwargs: Additional keyword arguments passed to the Device constructor
        """
        self._pv_x = pv_x
        self._pv_y = pv_y
        self._pv_pitch = pv_pitch
        self._pv_yaw = pv_yaw
        self._pv_lens1 = pv_lens1
        self._pv_lens2 = pv_lens2
        self._pv_lens3 = pv_lens3
        self._pv_lens4 = pv_lens4
        self._pv_lens5 = pv_lens5
        self._pv_lens6 = pv_lens6
        self._pv_lens7 = pv_lens7
        self._pv_lens8 = pv_lens8
        self._pv_lens9 = pv_lens9
        self._pv_lens10 = pv_lens10
        super().__init__(prefix, *args, **kwargs)
