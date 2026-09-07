"""LabJack LJT705 in 8-ID-I."""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO


class FastShutter(Device):
    """A device class for controlling fast shutters in the beamline.

    This class provides control over fast shutters used for beam control and safety.
    It includes functionality for opening and closing the shutter, as well as
    monitoring its status.
    """

    operation = Component(EpicsSignal, "State")
    logic = Component(EpicsSignal, "Lock")

    # Where the blade actually is, as opposed to what `operation` last asked
    # for. The two disagree whenever softglue is driving the shutter (logic on
    # "Override"), which is every external-trigger acquisition -- so this is
    # the only signal that can confirm the shutter really moved. Read-only, and
    # its enum is the honest one: "Open" / "Closed".
    state_rbv = Component(EpicsSignalRO, "State_RBV")
