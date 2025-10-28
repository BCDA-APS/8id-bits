"""
Aerotech motorized stages.

Plan these demonstrations:

Motivation: Move sample during measurement to avoid beam damage.

1. step scan Eiger4M v. sample.x
1. fly scan Eiger4M v. sample.x
  1. start motor moving
  2. start Eiger4M with one acquisition of n frames
  3. Eiger4M completes acquisition
  4. motor stops
1. fly scan Eiger4M v. sample.x and run DM workflow
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsMotor


class AerotechSampleStage(Device):
    """Device representing the sample stage with Aerotech motors.

    This device controls the position of the sample stage using Aerotech motors,
    including x, y, z translations and roll, pitch, yaw rotations.
    """

    x = Component(EpicsMotor, "m1", labels=("sample", "motor"))
    y = Component(EpicsMotor, "m3", labels=("sample", "motor"))
    z = Component(EpicsMotor, "m2", labels=("sample", "motor"))
    # different IOC!! "8ide:"
    # roll = Component(EpicsMotor, prefix="8ide:m15", labels=("sample", "motor"))
    # pitch = Component(EpicsMotor, prefix="8ide:m14", labels=("sample", "motor"))
    # yaw = Component(EpicsMotor, prefix="8ide:m16", labels=("sample", "motor"))


class AerotechDetectorStage(Device):
    """Device representing the detector stage with Aerotech motors.

    This device controls the position of the detector stage using Aerotech motors,
    including x and y translations. Z translation and yaw rotation are to be added.
    """

    x = Component(EpicsMotor, "m4", labels=("detector", "motor"))
    y = Component(EpicsMotor, "m5", labels=("detector", "motor"))
    # TBA: z
    # TBA: yaw (rotation about y)


class AerotechRheometerStage(Device):
    """Device representing the rheometer stage with Aerotech motors.

    This device controls the position of the rheometer stage using Aerotech motors,
    including x, y, z translations and roll, pitch, yaw rotations.
    """

    x = Component(EpicsMotor, "m8", labels=("rheometer", "motor"))
    y = Component(EpicsMotor, "m9", labels=("rheometer", "motor"))
    z = Component(EpicsMotor, "m7", labels=("rheometer", "motor"))
    # TODO: Are these motors assigned correctly?
    roll = Component(EpicsMotor, "m10", labels=("rheometer", "motor"))
    pitch = Component(EpicsMotor, "m11", labels=("rheometer", "motor"))
    yaw = Component(EpicsMotor, "m12", labels=("rheometer", "motor"))


