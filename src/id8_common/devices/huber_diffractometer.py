"""
Define motors in Huber diffractometer
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsMotor


class Huber_Diffractometer(Device):
    """The Huber diffractometer at 8-ID-E.

    Built from configs/devices.yml as `huber`, prefix "8ideSoft:CR8-E1:", so
    `nu` is 8ideSoft:CR8-E1:m4 and so on. nu/delta/mu/eta/chi/phi are the six
    circles; x/y/z translate the sample on top of them.

    `nu` and `delta` double as lambda2M's swing angles in
    plans/set/device_position.yaml. plans/set/select_device.py never drives
    them; only an explicit move_detector_axes() call does.
    """

    nu = Component(EpicsMotor, "m4", name="nu")
    delta = Component(EpicsMotor, "m5", name="delta")
    mu = Component(EpicsMotor, "m6", name="mu")
    eta = Component(EpicsMotor, "m7", name="eta")
    chi = Component(EpicsMotor, "m8", name="chi")
    phi = Component(EpicsMotor, "m9", name="phi")
    y = Component(EpicsMotor, "m10", name="sample_y")
    z = Component(EpicsMotor, "m11", name="sample_z")
    x = Component(EpicsMotor, "m15", name="sample_x")
