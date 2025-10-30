"""
Define motors in Huber diffractometer
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsMotor


class Huber_Diffractometer(Device):
    """A device class for controlling Huber diffractometers in the beamline.

    This class provides control over Huber diffractometers used for X-ray
    diffraction experiments. It includes functionality for controlling angles,
    positions, and other diffractometer parameters.
    """

    nu = Component(EpicsMotor, "m4", name="nu")
    delta = Component(EpicsMotor, "m5", name="delta")
    mu = Component(EpicsMotor, "m6", name="mu")
    eta = Component(EpicsMotor, "m7", name="eta")
    chi = Component(EpicsMotor, "m8", name="chi")
    phi = Component(EpicsMotor, "m9", name="phi")
    sample_y = Component(EpicsMotor, "m10", name="sample_y")
    sample_z = Component(EpicsMotor, "m11", name="sample_z")
    sample_x = Component(EpicsMotor, "m15", name="sample_x")
