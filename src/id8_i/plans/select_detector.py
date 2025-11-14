"""
Detector selection and configuration plans for the 8ID-I beamline.

This module provides plans for selecting and configuring different detectors
(Eiger4M and Rigaku3M), including setting up their positions, beam centers,
and workflow parameters.
"""

from typing import Literal

import bluesky.plan_stubs as bps
from apsbits.core.instrument_init import oregistry

pv_registers = oregistry["pv_registers"]
detector = oregistry["detector"]


def select_detector(det: Literal["eiger4M", "rigaku3M"]):
    """Select and configure a detector for data collection.

    This plan sets up the detector-specific parameters including workflow name,
    Q-map file, detector position, and beam center coordinates.

    Args:
        det: Detector to select ("eiger" or "rigaku")

    Yields:
        Generator: Bluesky plan messages
    """
    if det == "eiger4M":
        if pv_registers.det_name.get() == det:
            pass
        else:
            yield from bps.mv(pv_registers.workflow_name, "xpcs8-boost-corr")
            yield from bps.mv(pv_registers.qmap_file, "eiger4m_qmap_default.hdf")

            det_x_position = pv_registers.eiger_det_x0.get()
            det_y_position = pv_registers.eiger_det_y0.get()
            db_x_coord = pv_registers.eiger_db_x0.get()
            db_y_coord = pv_registers.eiger_db_y0.get()

            yield from bps.mv(pv_registers.current_det_x0, det_x_position)
            yield from bps.mv(pv_registers.current_det_y0, det_y_position)
            yield from bps.mv(pv_registers.current_db_x0, db_x_coord)
            yield from bps.mv(pv_registers.current_db_y0, db_y_coord)
            yield from bps.mv(detector.x, det_x_position, detector.y, det_y_position)

            yield from bps.mv(pv_registers.det_name, det)
            yield from bps.mv(pv_registers.det_pixel_size, 75e-6)

    elif det == "rigaku3M":
        if pv_registers.det_name.get() == det:
            pass
        else:
            yield from bps.mv(pv_registers.workflow_name, "xpcs8-boost-corr")
            yield from bps.mv(pv_registers.qmap_file, "rigaku3m_qmap_default.hdf")

            det_x_position = pv_registers.rigaku_det_x0.get()
            det_y_position = pv_registers.rigaku_det_y0.get()
            db_x_coord = pv_registers.rigaku_db_x0.get()
            db_y_coord = pv_registers.rigaku_db_y0.get()

            yield from bps.mv(pv_registers.current_det_x0, det_x_position)
            yield from bps.mv(pv_registers.current_det_y0, det_y_position)
            yield from bps.mv(pv_registers.current_db_x0, db_x_coord)
            yield from bps.mv(pv_registers.current_db_y0, db_y_coord)
            yield from bps.mv(detector.x, det_x_position, detector.y, det_y_position)

            yield from bps.mv(pv_registers.det_name, det)
            yield from bps.mv(pv_registers.det_pixel_size, 76e-6)

    elif det == "tempus":
        if pv_registers.det_name.get() == det:
            pass
        else:
            yield from bps.mv(pv_registers.workflow_name, "xpcs8-boost-corr")
            yield from bps.mv(pv_registers.qmap_file, "tempus_qmap_default.hdf")

            det_x_position = pv_registers.tempus_det_x0.get()
            det_y_position = pv_registers.tempus_det_y0.get()
            db_x_coord = pv_registers.tempus_db_x0.get()
            db_y_coord = pv_registers.tempus_db_y0.get()

            yield from bps.mv(pv_registers.current_det_x0, det_x_position)
            yield from bps.mv(pv_registers.current_det_y0, det_y_position)
            yield from bps.mv(pv_registers.current_db_x0, db_x_coord)
            yield from bps.mv(pv_registers.current_db_y0, db_y_coord)
            yield from bps.mv(detector.x, det_x_position, detector.y, det_y_position)

            yield from bps.mv(pv_registers.det_name, det)
            yield from bps.mv(pv_registers.det_pixel_size, 55e-6)

    else:
        print("Detector name must be eiger4M, rigaku3M, or tempus")
