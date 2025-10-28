"""Slit Devices at 8-ID."""

from apstools.devices.positioner_soft_done import PVPositionerSoftDone
from apstools.synApps.db_2slit import Optics2Slit1D
from apstools.synApps.db_2slit import Optics2Slit2D_HV
from ophyd import Component as cpt
from ophyd import EpicsSignal


class ID8Optics2Slit1D(Optics2Slit1D):
    """A device class for controlling 1D slits with enhanced positioning.

    This class extends the base Optics2Slit1D class to provide more precise
    control over slit positioning with tighter tolerances and additional
    synchronization capabilities.
    """

    xn = cpt(
        PVPositionerSoftDone,
        "xn",
        setpoint_pv=".VAL",
        readback_pv=".RBV",
        tolerance=8e-4,
    )
    xp = cpt(
        PVPositionerSoftDone,
        "xp",
        setpoint_pv=".VAL",
        readback_pv=".RBV",
        tolerance=8e-4,
    )
    size = cpt(
        PVPositionerSoftDone,
        "size",
        setpoint_pv=".VAL",
        readback_pv=".RBV",
        tolerance=2e-3,
    )
    center = cpt(
        PVPositionerSoftDone,
        "center",
        setpoint_pv=".VAL",
        readback_pv=".RBV",
        tolerance=8e-4,
    )
    sync = cpt(EpicsSignal, "doSync", put_complete=True, kind="omitted")


class ID8Optics2Slit2D_HV(Optics2Slit2D_HV):
    """A device class for controlling 2D horizontal-vertical slits.

    This class provides control over 2D slits with horizontal and vertical
    motion capabilities. It includes functionality for controlling gap sizes,
    center positions, and individual blade positions.
    """

    h = cpt(ID8Optics2Slit1D, "H")
    v = cpt(ID8Optics2Slit1D, "V")
