"""
QNW temperature controller -- the AIR stage.

Careful: qnw_vac_device.py in this directory declares a class with the same
name, QnwDevice, and the same four Components. The two are distinguished only
by which one configs/devices.yml points at for which prefix -- this one builds
qnw_env1/2/3 on "8idiSoft:QNWenv_N:", the vacuum module builds qnw_vac1/2/3 on
"8idiSoft:QNWvac_N:". Keep them in step, or merge them deliberately; do not
assume an edit here reaches the vacuum controllers.
"""

from apstools.devices import PVPositionerSoftDoneWithStop
from ophyd import Component
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO
from ophyd import Signal


class QnwDevice(PVPositionerSoftDoneWithStop):
    """One air-stage QNW temperature controller.

    A soft positioner: `move(setpoint)` blocks until `readback` (SH_RBV) is
    within `tolerance` of `setpoint` (TARG), since the controller has no
    done-moving flag of its own. `ramprate` (RAMP) limits how fast it gets
    there. Driven from plans/set/qnw_plans.py.
    """

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
