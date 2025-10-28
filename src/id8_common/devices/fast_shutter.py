"""LabJack LJT705 in 8-ID-I."""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal


class FastShutter(Device):
    """A device class for controlling fast shutters in the beamline.

    This class provides control over fast shutters used for beam control and safety.
    It includes functionality for opening and closing the shutter, as well as
    monitoring its status.
    """

    operation = Component(EpicsSignal, "State")
    logic = Component(EpicsSignal, "Lock")
