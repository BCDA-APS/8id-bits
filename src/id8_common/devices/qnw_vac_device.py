"""
QNW temperature controller -- the VACUUM stage.

Careful: qnw_device.py in this directory declares a class with the same name,
QnwDevice, and the same four Components. See the note there. This module is
the one configs/devices.yml uses for qnw_vac1/2/3 on "8idiSoft:QNWvac_N:", and
it is what plans/set/qnw_plans.set_qnw() actually drives.
"""

from apstools.devices import PVPositionerSoftDoneWithStop
from ophyd import Component
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO
from ophyd import Signal


class QnwDevice(PVPositionerSoftDoneWithStop):
    """One vacuum-stage QNW temperature controller.

    A soft positioner: `move(setpoint)` blocks until `readback` (SH_RBV) is
    within `tolerance` of `setpoint` (TARG), since the controller has no
    done-moving flag of its own. `ramprate` (RAMP) limits how fast it gets
    there. Driven from plans/set/qnw_plans.py.
    """

    readback = Component(EpicsSignalRO, "SH_RBV", kind="hinted", auto_monitor=True)
    setpoint = Component(EpicsSignal, "TARG", kind="normal", put_complete=True)
    tolerance = Component(Signal, value=0.1, kind="config")
    ramprate = Component(EpicsSignal, "RAMP", kind="normal", put_complete=True)
