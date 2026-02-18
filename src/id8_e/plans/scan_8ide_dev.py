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
from bluesky.plans import tune_centroid 
from bluesky.callbacks.fitting import PeakStats
from apstools.plans.alignment import lineup2

from .shutter_logic import blockbeam, pre_align, showbeam, PIND_out
from .shutter_logic import pre_align
from .shutter_logic import showbeam
from .sample_info_unpack import gen_folder_prefix
from .nexus_acq_eiger_int import setup_eiger_int_series
import time


huber = oregistry["huber"]
filter_beam = oregistry["filter_8ide"]
tetramm1 = oregistry["tetramm1"]
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]

pv_registers = oregistry["pv_registers"]

ps = PeakStats('motor', 'det.stats2.total')

def save_images(det, save_img, num_frames=1, file_path=None):
    """
    toggle saving (1) or not saving (0) images for either eiger4m or lambda2m.

    det: detector instance (eiger4M or lambda2M)
    save_img: 1 save, 0 don't save
    num_frames: number of frames to capture per point (default 1)
    file_path: override base path (default uses a safe local path)
    """

    if save_img not in (0, 1):
        raise valueerror("save_img must be 1 or 0 (to save or not to save)")

    def _has(obj, attr):
        return getattr(obj, attr, None) is not None

    if file_path is None:
        file_path = "/gdata/dm/8ID/8IDE/2026-1/kisiel202602/data/bluesky/"

    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

    if save_img == 1:
        if _has(det.cam, "fw_enable"):
            yield from bps.mv(det.cam.fw_enable, 1)
        if _has(det.cam, "save_files"):
            yield from bps.mv(det.cam.save_files, 1)

        folder_prefix = gen_folder_prefix()
        file_header = f"{folder_prefix}"
        file_name = f"{folder_prefix}"

        if _has(det, "hdf1"):
            if _has(det.hdf1, "num_capture"):
                yield from bps.mv(det.hdf1.num_capture, num_frames)

            yield from bps.mv(pv_registers.file_name, file_name)
            yield from bps.mv(pv_registers.file_path, file_path)
            yield from bps.mv(
                pv_registers.metadata_full_path,
                f"{file_path}/{file_name}_metadata.hdf",
            )

            if _has(det.hdf1, "file_name"):
                yield from bps.mv(det.hdf1.file_name, file_name)
            if _has(det.hdf1, "file_path"):
                yield from bps.mv(det.hdf1.file_path, file_path)

        if is_eiger:
            if _has(det.cam, "trigger_mode"):
                yield from bps.mv(det.cam.trigger_mode, "Internal Enable")
            if _has(det.cam, "num_triggers"):
                yield from bps.mv(det.cam.num_triggers, 1)
            if _has(det.cam, "num_images"):
                yield from bps.mv(det.cam.num_images, num_frames)

            # eiger internal series configuration 
            yield from setup_eiger_int_series(
                acq_time=float(det.cam.acquire_time.get()) if _has(det.cam, "acquire_time") else 1.0,
                num_frames=num_frames,
                file_header=file_header,
                file_name=file_name,
            )

        else:
            if _has(det.cam, "num_images"):
                yield from bps.mv(det.cam.num_images, num_frames)

    else:
        # disable writing 
        if _has(det.cam, "fw_enable"):
            yield from bps.mv(det.cam.fw_enable, 0)
        if _has(det.cam, "save_files"):
            yield from bps.mv(det.cam.save_files, 0)
        if _has(det, "hdf1") and _has(det.hdf1, "enable"):
            yield from bps.mv(det.hdf1.enable, 0)

def att(att_ratio: Optional[float] = None):
    """Set the attenuation ratio with multiple attempts.

    Args:
        att_ratio: Attenuation ratio to set (0-15)
    """
    if att_ratio < 1.0:
        att_ratio = 1.0
    else:
        pass
    filter_beam.attenuation.move(att_ratio)
    time.sleep(0.5)


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


def dscan(motor, 
          rel_begin, 
          rel_end, 
          num_pts, 
          count_time, 
          det=eiger4M, 
          att_ratio=7, 
          save_img=1
          ):
    
    """
    usage:
        RE(dscan(motor, rel_begin, rel_end, num_pts, count_time))

    args:
        motor: ophyd positioner (e.g. huber.eta, sample.y)
        rel_begin: relative start (motor units)
        rel_end: relative end (motor units)
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (default: eiger4M)
        att_ratio: attenuation ratio
        save_img: toggle saving images
    """
    pre_align()
    att(att_ratio)
    PIND_out()
    det.cam.num_images.put(1)
    det.cam.acquire_time.put(count_time)
    det.cam.acquire_period.put(count_time)

    yield from save_images(det, save_img)
    showbeam()
    yield from bp.rel_scan([det], motor, rel_begin, rel_end, num_pts)
    blockbeam()

def d2scan(
    motor1, 
    rel_begin1,
    rel_end1,
    motor2,
    rel_begin2,
    rel_end2,
    num_pts,
    count_time,
    det=eiger4M,
    att_ratio=7,
    save_img=1,
):
    """
    usage:
        RE(d2scan(m1, s1, f1, m2, s2, f2, npts, ct))

    args:
        motor1, motor2: ophyd positioners (e.g. huber.eta, huber.chi)
        rel_begin1, rel_end1: relative start/end for motor1
        rel_begin2, rel_end2: relative start/end for motor2
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (default: eiger4M)
        att_ratio: attenuation ratio
        save_img: toggle saving images
    """
    pre_align()
    att(att_ratio)
    yield from bps.mv(det.cam.num_images, 1)
    yield from bps.mv(det.cam.acquire_time, count_time)
    yield from bps.mv(det.cam.acquire_period, count_time)
    yield from save_images(det, save_img)
    showbeam()
    yield from bp.rel_scan([det], motor1, rel_begin1, rel_end1, motor2, rel_begin2, rel_end2, num_pts)
    blockbeam()
    

def ascan(motor, 
          rel_begin, 
          rel_end, 
          num_pts, 
          count_time, 
          det=eiger4M, 
          att_ratio=7, 
          save_img=1):
    
    """
    usage:
        RE(dscan(motor, rel_begin, rel_end, num_pts, count_time))

    args:
        motor: ophyd positioner (e.g. huber.eta, sample.y)
        rel_begin: relative start (motor units)
        rel_end: relative end (motor units)
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (default: eiger4M)
        att_ratio: attenuation ratio
        save_img: toggle saving images
    """
    pre_align()
    att(att_ratio)
    yield from bps.mv(eiger4M.cam.num_images, 1)
    yield from bps.mv(det.cam.acquire_time, count_time)
    yield from bps.mv(det.cam.acquire_period, count_time)
    yield from save_images(det, save_img)
    showbeam()
    yield from bp.scan([det], motor, rel_begin, rel_end, num_pts)
    blockbeam()

def a2scan(
    motor1, 
    rel_begin1,
    rel_end1,
    motor2,
    rel_begin2,
    rel_end2,
    num_pts,
    count_time,
    det=eiger4M,
    att_ratio=7,
    save_img=1,
):
    """
    usage:
        RE(d2scan(m1, s1, f1, m2, s2, f2, npts, ct))

    args:
        motor1, motor2: ophyd positioners (e.g. huber.eta, huber.chi)
        rel_begin1, rel_end1: relative start/end for motor1
        rel_begin2, rel_end2: relative start/end for motor2
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (default: eiger4M)
        att_ratio: attenuation ratio
        save_img: toggle saving images
    """
    pre_align()
    att(att_ratio)
    yield from bps.mv(eiger4M.cam.num_images, 1)
    yield from bps.mv(det.cam.acquire_time, count_time)
    yield from bps.mv(det.cam.acquire_period, count_time)
    yield from save_images(det, save_img)
    showbeam()
    yield from bp.scan([det], motor1, rel_begin1, rel_end1, motor2, rel_begin2, rel_end2, num_pts)
    blockbeam()


def lineup_8id(motor, 
          rel_begin, 
          rel_end, 
          num_pts, 
          count_time=1, 
          det=lambda2M, 
          att_ratio=1, 
          save_img=1
          ):
    
    """
    usage:
        RE(lineup_8id(motor, rel_begin, rel_end, num_pts, count_time))

    args:
        motor: ophyd positioner (e.g. huber.eta, sample.y)
        rel_begin: relative start (motor units)
        rel_end: relative end (motor units)
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (default: eiger4M)
        att_ratio: attenuation ratio
        save_img: toggle saving images
    """
    pre_align()
    att(att_ratio)
    yield from bps.mv(
        det.cam.num_images, 1,
        det.cam.acquire_time, count_time,
        det.cam.acquire_period, 1,
    )
    # det.cam.num_images.put(1)
    # det.cam.acquire_time.put(count_time)
    # det.cam.acquire_period.put(1)

    yield from save_images(det, save_img)
    showbeam()
    yield from lineup2([det.stats2.total, det], motor, rel_begin, rel_end, num_pts)
    blockbeam()