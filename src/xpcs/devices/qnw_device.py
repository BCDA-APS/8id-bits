"""
QNW temperature controller
"""

from apstools.devices import PVPositionerSoftDoneWithStop
from ophyd import Component
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO
from ophyd import Signal


class QnwDevice(PVPositionerSoftDoneWithStop):
    readback = Component(EpicsSignalRO, "SH_RBV", kind="hinted", auto_monitor=True)
    setpoint = Component(EpicsSignal, "TARG", kind="normal", put_complete=True)
    tolerance = Component(Signal, value=0.1, kind="config")
    ramprate = Component(EpicsSignal, "RAMP", kind="normal", put_complete=True)


    def qnw_register(self, temp_zone):
        """
        Return the indexed qnw temp zone.
        """
        return getattr(self, f"{temp_zone}")


# TODO:
# Add read-only temperatures to watch.  (What PVs?)