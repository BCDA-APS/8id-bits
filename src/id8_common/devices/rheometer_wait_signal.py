"""
Softglue counter. Increments when the MCR rheometer finishes current program.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO


class Rheometer_Wait(Device):
    """Device representing the rheometer wait signal.

    This device monitors a softglue counter that increments when the MCR rheometer
    finishes its current program. It provides read-only access to the pulse count.
    """

    pulse_count = Component(EpicsSignalRO, "8idMZ4:SG:UpCntr-1_COUNTS")
    rheo_output = Component(EpicsSignalRO, "8idiSoft:LJT705:Ai4")
