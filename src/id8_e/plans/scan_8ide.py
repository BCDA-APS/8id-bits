"""
Scanning plans for the 8ID-E beamline.

This module provides plans for scanning various motors and detectors at the
8ID-E beamline, including sample and with attenuation control.
"""

from typing import Optional

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from bluesky import plans as bp
from ophyd import Device

from .shutter_logic import blockbeam
from .shutter_logic import pre_align
from .shutter_logic import showbeam
from .sample_info_unpack import gen_folder_prefix
from .nexus_acq_eiger_int import setup_eiger_int_series

huber = oregistry["huber"]
filter_beam = oregistry["filter_8ide"]
tetramm1 = oregistry["tetramm1"]
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]

pv_registers = oregistry["pv_registers"]

def save_images(save_img):
    """
    Toggle saving (1) or not saving (0) Eiger detector images. Uses path defined in String Register 5.
    """

    if save_img not in (0, 1):
        raise ValueError("save_img must be 1 or 0 (to save or not to save)")

    if save_img == 1:
        # enable writing
        yield from bps.mv(eiger4M.cam.fw_enable, 1)
        yield from bps.mv(eiger4M.cam.save_files, 1)

        cycle_name = pv_registers.cycle_name.get()
        exp_name = pv_registers.experiment_name.get()
        # file_name = pv_registers.file_name.get()
        # file_path = pv_registers.file_path.get()
        # print('file path: ', file_path)
        # file_path = f"/gdata/dm/8ID/8IDE/{cycle_name}/{exp_name}/data/{file_name}/"
        file_path = '/home/beams10/8IDIUSER/sdmarks/'

        yield from bps.mv(eiger4M.cam.trigger_mode, "Internal Enable")
        yield from bps.mv(eiger4M.cam.acquire_time, 1)
        yield from bps.mv(eiger4M.cam.acquire_period, 1)
        yield from bps.mv(eiger4M.cam.num_triggers, 1)
        yield from bps.mv(eiger4M.cam.num_images, 1)
   
        folder_prefix = gen_folder_prefix()

        file_header = f"{folder_prefix}"
        file_name = f"{folder_prefix}"
        
        
        
        yield from bps.mv(eiger4M.hdf1.num_capture, 1)
        yield from bps.mv(pv_registers.file_name, file_name)
        yield from bps.mv(pv_registers.file_path, file_path)
        yield from bps.mv(pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf",)
        yield from bps.mv(eiger4M.hdf1.file_name, file_name)
        yield from bps.mv(eiger4M.hdf1.file_path, file_path)

        yield from setup_eiger_int_series(1, 1, file_header, file_name)

    else:
        # disable writing
        yield from bps.mv(eiger4M.cam.fw_enable, 0)
        yield from bps.mv(eiger4M.cam.save_files, 0)
        # yield from bps.mv(eig4M.hdf1.enable, 0)
        # yield from bps.mv(eig4M.hdf1.capture, 0)
        

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
    yield from bp.rel_scan([det], huber.sample_x, rel_begin, rel_end, num_pts)
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
    yield from bp.rel_scan([det], huber.sample_y, rel_begin, rel_end, num_pts)
    yield from blockbeam()

def z_lup(
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
    yield from bp.rel_scan([det], huber.sample_z, rel_begin, rel_end, num_pts)
    yield from blockbeam()

def eta_align(
    rel_begin: float = -3,
    rel_end: float = 3,
    num_pts: int = 60,
    att_ratio: int = 7,
    det: Device = eiger4M,
    # save_image: int = 1
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
    yield from bp.rel_scan([det], huber.eta, rel_begin, rel_end, num_pts)
    yield from blockbeam()

def chi_align(
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
    yield from bp.rel_scan([det], huber.chi, rel_begin, rel_end, num_pts)
    yield from blockbeam()

def eta_scan(
    rel_begin: float = -3,
    rel_end: float = 3,
    num_pts: int = 60,
    att_ratio: int = 7,
    det: Device = eiger4M,
    save_img: int = 1
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

    yield from save_images(save_img)

    yield from showbeam()
    yield from bp.rel_scan([det, eiger4M.stats1.total], huber.eta, rel_begin, rel_end, num_pts)
    yield from blockbeam()

# def l_scan(
#     rel_begin: float = -3,
#     rel_end: float = 3,
#     num_pts: int = 60,
#     att_ratio: int = 7,
#     det: Device = eiger4M,
# ):
#     """Perform a relative scan along the sample Y axis.

#     Args:
#         rel_begin: Start position relative to current position (mm)
#         rel_end: End position relative to current position (mm)
#         num_pts: Number of points in the scan
#         att_level: Attenuation level to use (0-15)
#         det: Detector to use for the scan
#     """
#     yield from pre_align()
#     yield from att(att_ratio)

#     yield from showbeam()
#     yield from bp.rel_scan([det], sample.eta, rel_begin, rel_end, num_pts)
#     yield from blockbeam()

# def l_scan(det, sixc, h0, k0, l0, dh, dk, dl, npts):
#     ts = linspace(-0.5, 0.5, npts)
#     for t in ts:
#         yield from bps.mv(
#             sixc.h, h0 + t * dh,
#             sixc.k, k0 + t * dk,
#             sixc.l, l0 + t * dl,
#         )
#         yield from bps.trigger_and_read([det])