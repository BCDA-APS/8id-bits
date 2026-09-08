"""Alicat PCD pressure controller (8-ID).

Two units are installed, PCD1 and PCD2, served by the 8idAlicat IOC over a MOXA
terminal server. Start the IOC with ``start_Alicat.sh start``; the same script
takes ``status``, ``stop`` and ``caqtdm``.

The full record set in ``Alicat_PCD.db`` is much larger than this class -- max
pressure, ramp rate, Run/Pause, Vent, valve positions, units, firmware, asyn
status. Only the three signals an XPCS measurement needs are exposed here:
the pressure readback and the setpoint (write and readback). Add the others if
and when a plan needs them; every extra Component is another PV this device
waits on at startup.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO


class AlicatPressureController(Device):
    """One Alicat PCD unit.

    Example::

        pcd1 = AlicatPressureController("8idAlicat:PCD1:", name="pcd1")
        pcd1.pressure.get()          # current pressure
        pcd1.setpoint_rbv.get()      # what the controller thinks it is aiming at
        pcd1.setpoint.put(1000.0)    # ask for a new pressure

    ``setpoint`` is the demand and ``setpoint_rbv`` is the controller's own
    readback of it: after a put, read the _RBV to confirm the unit accepted the
    value rather than trusting the write.

    Units are whatever the controller is configured for -- the caQtDM screen
    shows Pa -- and are NOT converted here. `Units_RBV` on the IOC reports them.
    """

    #: Measured pressure.
    pressure = Component(EpicsSignalRO, "Pressure_RBV", kind="hinted")

    #: Demanded pressure, and the controller's readback of that demand.
    setpoint = Component(EpicsSignal, "Setpoint", kind="normal")
    setpoint_rbv = Component(EpicsSignalRO, "Setpoint_RBV", kind="normal")
