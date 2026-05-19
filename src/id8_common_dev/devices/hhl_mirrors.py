"""
High-heat load mirrors in station 8-ID-A.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsMotor


class HHL_Mirror1(Device):
    """
    High-Heat Load Mirror 1 in 8-ID-A.

    The area detector on the 'flag' motor defined separately
    (in ./area_detectors.py) as 'flag1'.
    """

    x = Component(EpicsMotor, "FMBO:m4")
    y = Component(EpicsMotor, "FMBO:m2")
    coarse_pitch = Component(EpicsMotor, "FMBO:m6")
    fine_pitch = Component(EpicsMotor, "FMBO:Piezo:m1")
    flag = Component(EpicsMotor, "CR8-A1:m5")


class HHL_Mirror2(Device):
    """
    High-Heat Load Mirror 2 in 8-ID-A.

    The area detector on the 'flag' motor defined separately
    (in ./area_detectors.py) as 'flag2'.
    """

    x = Component(EpicsMotor, "FMBO:m3")
    y = Component(EpicsMotor, "FMBO:m1")
    coarse_pitch = Component(EpicsMotor, "FMBO:m5")
    fine_pitch = Component(EpicsMotor, "FMBO:Piezo:m2")
    flag = Component(EpicsMotor, "CR8-A1:m6")
    bender1 = Component(EpicsMotor, "FMBO:m7")
    bender2 = Component(EpicsMotor, "FMBO:m8")
