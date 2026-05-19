"""Monochromator device configuration for 8ID beamline.

This module defines the monochromator device which provides read-only access to the
monochromator energy readback.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO


class Mono(Device):
    """Device representing the monochromator.

    This device provides read-only access to the monochromator energy readback value.
    """

    energy_readback = Component(EpicsSignalRO, "BraggERdbkAO", name="energy_readback")
