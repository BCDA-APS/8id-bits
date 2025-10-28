"""
Scanning plans for the 8ID-I beamline.

This module provides plans for scanning various motors and detectors at the
8ID-I beamline, including sample and rheometer stages, with attenuation control.
"""

from typing import Optional

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from bluesky import plans as bp
from ophyd import Device

from .shutter_logic import blockbeam
from .shutter_logic import pre_align
from .shutter_logic import showbeam

rheometer = oregistry["rheometer"]
sample = oregistry["sample"]
filter_beam = oregistry["filter_8ide"]
tetramm1 = oregistry["tetramm1"]


def att(att_ratio: Optional[float] = None):
    """Set the attenuation ratio with multiple attempts.

    Args:
        att_ratio: Attenuation ratio to set (0-15)
    """
    if att_ratio < 1.0:
        att_ratio = 1.0
    else:
        pass
    yield from bps.mv(filter_beam.attenuation, att_ratio)
    yield from bps.sleep(0.5)


def x_lup(
    rel_begin: float = -3,
    rel_end: float = 3,
    num_pts: int = 60,
    att_ratio: int = 7,
    det: Device = tetramm1,
):
    """Perform a relative scan along the sample X axis.

    Args:
        rel_begin: Start position relative to current position (mm)
        rel_end: End position relative to current position (mm)
        num_pts: Number of points in the scan
        att_level: Attenuation level to use (0-15)
        det: Detector to use for the scan
    """
    yield from pre_align()
    yield from att(att_ratio)

    yield from showbeam()
    yield from bp.rel_scan([det], sample.x, rel_begin, rel_end, num_pts)
    yield from blockbeam()


def y_lup(
    rel_begin: float = -3,
    rel_end: float = 3,
    num_pts: int = 60,
    att_ratio: int = 7,
    det: Device = tetramm1,
):
    """Perform a relative scan along the sample Y axis.

    Args:
        rel_begin: Start position relative to current position (mm)
        rel_end: End position relative to current position (mm)
        num_pts: Number of points in the scan
        att_level: Attenuation level to use (0-15)
        det: Detector to use for the scan
    """
    yield from pre_align()
    yield from att(att_ratio)

    yield from showbeam()
    yield from bp.rel_scan([det], sample.y, rel_begin, rel_end, num_pts)
    yield from blockbeam()


def rheo_x_lup(
    rel_begin: float = -3,
    rel_end: float = 3,
    num_pts: int = 60,
    att_ratio: int = 7,
    det: Device = tetramm1,
):
    """Perform a relative scan along the rheometer X axis.

    Args:
        rel_begin: Start position relative to current position (mm)
        rel_end: End position relative to current position (mm)
        num_pts: Number of points in the scan
        att_level: Attenuation level to use (0-15)
        det: Detector to use for the scan
    """
    yield from pre_align()
    yield from att(att_ratio)

    yield from showbeam()
    yield from bp.rel_scan([det], rheometer.x, rel_begin, rel_end, num_pts)
    yield from blockbeam()


def rheo_y_lup(
    rel_begin: float = -3,
    rel_end: float = 3,
    num_pts: int = 60,
    att_ratio: int = 7,
    det: Device = tetramm1,
):
    """Perform a relative scan along the rheometer Y axis.

    Args:
        rel_begin: Start position relative to current position (mm)
        rel_end: End position relative to current position (mm)
        num_pts: Number of points in the scan
        att_level: Attenuation level to use (0-15)
        det: Detector to use for the scan
    """
    yield from pre_align()
    yield from att(att_ratio)

    yield from showbeam()
    yield from bp.rel_scan([det], rheometer.y, rel_begin, rel_end, num_pts)
    yield from blockbeam()


def rheo_set_x_lup(
    att_ratio: int = 7,
    det: Device = tetramm1,
):
    """Perform a series of scans at specific rheometer X positions.

    This plan moves the rheometer to three specific X positions and performs
    relative scans around each position.

    Args:
        att_level: Attenuation level to use (0-15)
        det: Detector to use for the scan
    """
    yield from pre_align()
    yield from att(att_ratio)

    yield from bps.mv(rheometer.x, -14.0)
    yield from showbeam()
    yield from bp.rel_scan([det], rheometer.x, -0.5, 0.5, 100)
    yield from blockbeam()

    yield from bps.mv(rheometer.x, -2.6)
    yield from showbeam()
    yield from bp.rel_scan([det], rheometer.x, -0.5, 0.5, 100)
    yield from blockbeam()

    # yield from bps.mv(rheometer.x, -3.53)
    # yield from showbeam()
    # yield from bp.rel_scan([det], rheometer.x, -8, 8, 160)
    # yield from blockbeam()
