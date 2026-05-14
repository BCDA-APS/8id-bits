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

from .shutter_logic import blockbeam, pre_align, showbeam, PIND_status
from .sample_info_unpack import gen_folder_prefix
from .nexus_acq_eiger_int import setup_eiger_int_series
import bluesky.preprocessors as bpp
import numpy as np
import time

huber = oregistry["huber"]
filter_beam = oregistry["filter_8ide"]
tetramm1 = oregistry["tetramm1"]
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]

pv_registers = oregistry["pv_registers"]

# ps = PeakStats('motor', 'det.stats2.total')

def save_images(det, save_img, num_pts, num_frames=1, file_path=None):
    """
    Toggle saving (1) or not saving (0) images.

    det: detector instance (eiger4M, lambda2M, or tetramm)
    save_img: 1 save, 0 don't save
    num_frames: number of frames to capture per point (default 1)
    file_path: override base path (default uses a safe local path)
    """
    if save_img not in (0, 1):
        raise ValueError("save_img must be 1 or 0 (to save or not to save)")
    
    if file_path is None:
        # file_path = "/gdata/dm/8ID/8IDE/2026-1/marks202604/data/bluesky/"
        # file_path = "/home/beams/8IDIUSER/sdmarks/tetramm_testing"
        file_path = f"{pv_registers.file_path.get()}/bluesky"
    
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_tetramm = "tetramm" in det.name.lower()

    def _has(obj, attr):
        return getattr(obj, attr, None) is not None
    print(file_path)
    print(det.hdf1.file_path)
    if save_img == 1:

        folder_prefix = gen_folder_prefix()
        file_name = folder_prefix
        file_header = folder_prefix

        if is_tetramm:
            yield from bps.mv(
                det.hdf1.file_path, file_path,
                det.hdf1.file_name, file_name,
                det.hdf1.num_capture, num_pts,
                det.hdf1.file_write_mode, 2,  # Stream mode
            )
            yield from bps.mv(
                pv_registers.file_name, file_name,
                pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf",
            )
            return  # capture managed manually in dscan
        
        # area detectors (eiger, lambda)
        if _has(det, "cam"):
            if _has(det.cam, "fw_enable"):
                yield from bps.mv(det.cam.fw_enable, 1)
            if _has(det.cam, "save_files"):
                yield from bps.mv(det.cam.save_files, 1)

        if _has(det, "hdf1"):
            if _has(det.hdf1, "num_capture"):
                yield from bps.mv(det.hdf1.num_capture, num_pts)
            yield from bps.mv(
                pv_registers.file_name, file_name,
                # pv_registers.file_path, file_path,
                pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf",
            )
            print(file_path)
            print(det.hdf1.file_path)

            if _has(det.hdf1, "file_name"):
                yield from bps.mv(det.hdf1.file_name, file_name)
            if _has(det.hdf1, "file_path"):
                yield from bps.mv(det.hdf1.file_path, file_path)
            print(file_path)
            print(det.hdf1.file_path)
        if _has(det, "cam"):
            det.cam.trigger_mode.put(0)
            if is_eiger:
                if _has(det.cam, "trigger_mode"):
                    yield from bps.mv(det.cam.trigger_mode, "Internal Enable")
                if _has(det.cam, "num_triggers"):
                    yield from bps.mv(det.cam.num_triggers, 1)
                if _has(det.cam, "num_images"):
                    yield from bps.mv(det.cam.num_images, 1)
                setup_eiger_int_series(
                    acq_time=float(det.cam.acquire_time.get()) if _has(det.cam, "acquire_time") else 1.0,
                    num_frames=num_frames,
                    file_header=file_header,
                    file_name=file_name,
                )
            else:
                if _has(det.cam, "num_images"):
                    yield from bps.mv(det.cam.num_images, num_frames)
        print(file_path)
        print(det.hdf1.file_path)

def dscan(motor, rel_begin, rel_end, num_pts, count_time,
          det=eiger4M, att_ratio=1e6, save_img=1):
    
    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name

    if det == eiger4M or det == lambda2M:
        yield from bps.mv(
            det.cam.num_images, 1,
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
        )
    else:
        print('tetramm1 or other detector with simple acquire trigger')
    print(det.hdf1.num_capture.get(), ' checkpoint2')
    yield from save_images(det, save_img, num_pts)

    if is_tetramm and save_img == 1:
        det.hdf1.enable.put(1)
        det.hdf1.capture.put(1)
        print(f"TetrAMM HDF capture armed: {det.hdf1.file_name.get()}")

    det.hdf1.num_capture.put(num_pts)
    showbeam()
    yield from bp.rel_scan([det], motor, rel_begin, rel_end, num_pts)
    blockbeam()

    if is_tetramm and save_img == 1:
        det.hdf1.capture.put(0)
        print("TetrAMM HDF capture stopped.")

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
    PIND_out()
    yield from bps.mv(det.cam.num_images, 1)
    yield from bps.mv(det.cam.acquire_time, count_time)
    yield from bps.mv(det.cam.acquire_period, count_time)
    yield from save_images(det, save_img)
    showbeam()
    yield from bp.rel_scan([det], motor1, rel_begin1, rel_end1, motor2, rel_begin2, rel_end2, num_pts)
    blockbeam()
    det.hdf1.stage_sigs["enable"] = 1
    

def ascan(motor, 
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
    PIND_status(0)
    yield from bps.mv(eiger4M.cam.num_images, 1)
    yield from bps.mv(det.cam.acquire_time, count_time)
    yield from bps.mv(det.cam.acquire_period, count_time)
    yield from save_images(det, save_img)
    showbeam()
    yield from bp.scan([det], motor, rel_begin, rel_end, num_pts)
    blockbeam()
    det.hdf1.stage_sigs["enable"] = 1

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
    PIND_status(0)
    yield from bps.mv(eiger4M.cam.num_images, 1)
    yield from bps.mv(det.cam.acquire_time, count_time)
    yield from bps.mv(det.cam.acquire_period, count_time)
    yield from save_images(det, save_img)
    showbeam()
    yield from bp.scan([det], motor1, rel_begin1, rel_end1, motor2, rel_begin2, rel_end2, num_pts)
    blockbeam()
    det.hdf1.stage_sigs["enable"] = 1

def lineup_8id(motor, 
          rel_begin, 
          rel_end, 
          num_pts, 
          count_time=1, 
          det=eiger4M, 
          att_ratio=1, 
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
    PIND_status(0)
    yield from bps.mv(
        det.cam.num_images, 1,
        det.cam.acquire_time, count_time,
        det.cam.acquire_period, 1,
    )

    yield from save_images(det, save_img)
    showbeam()
    yield from lineup2([det.stats2.total, det], motor, rel_begin, rel_end, num_pts)
    blockbeam()
    det.hdf1.stage_sigs["enable"] = 1

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


def dscan_fast(motor, rel_begin, rel_end, num_pts, count_time,
               det=eiger4M, att_ratio=1e6, save_img=1):
    """
    Pre-armed software-trigger scan for fast acquisitions.

    For Eiger4M: Arms detector for all N frames upfront, then sends a software
    trigger at each point via det.cam.special_trigger_button. Eliminates
    per-point wait_for_plugins overhead.

    For TetrAMM: No pre-arming (continuous mode). Steps through positions and
    reads the current value at each point. HDF capture is armed/stopped manually.

    For Lambda2M: simple step scan with cam timing set up per-point. No 
    pre-arming (no software trigger support on those cams).

    args:
        motor: ophyd positioner
        rel_begin, rel_end: relative start/end (motor units)
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (eiger4M, lambda2M, or tetramm1)
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
    """
    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name.lower()
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())

    if is_tetramm:
        yield from save_images(det, save_img, num_pts)
        if save_img == 1:
            det.hdf1.enable.put(1)
            det.hdf1.capture.put(1)
            print(f"TetrAMM HDF capture armed: {det.hdf1.file_name.get()}")

        start_pos = motor.position
        positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
        pos_cache = {motor: None}

        def inner_tetramm():
            try:
                for pos in positions:
                    yield from bps.move_per_step({motor: pos}, pos_cache)
                    yield from bps.trigger_and_read([det, motor])
            finally:
                yield from bps.mv(motor, start_pos)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_tetramm()),
                [det, motor],
            )
        finally:
            if save_img == 1:
                det.hdf1.capture.put(0)
                print("TetrAMM HDF capture stopped.")
        return

    # Lambda2M 
    if not is_eiger:
        yield from bps.mv(
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
            det.cam.num_images, 1,
        )
        yield from save_images(det, save_img, num_pts)

        start_pos = motor.position
        positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
        pos_cache = {motor: None}

        def inner_area_det():
            try:
                for pos in positions:
                    yield from bps.move_per_step({motor: pos}, pos_cache)
                    yield from bps.trigger_and_read([det, motor])
            finally:
                yield from bps.mv(motor, start_pos)

        yield from bpp.stage_wrapper(
            bpp.run_wrapper(inner_area_det()),
            [det, motor],
        )
        return

    # Eiger4M
    # Set timing first (save_images reads acquire_time for setup_eiger_int_series)
    yield from bps.mv(
        det.cam.acquire_time, count_time,
        det.cam.acquire_period, count_time,
    )
    yield from save_images(det, save_img, num_pts)

    # Correct settings after save_images (setup_eiger_int_series resets:
    # trigger_mode="Internal Series", num_triggers=1, manual_trigger="Disable",
    # hdf1.num_capture=1)
    yield from bps.mv(
        det.cam.trigger_mode, "Internal Enable",
        det.cam.manual_trigger, "Enable",
        det.cam.num_images, 1,
        det.cam.num_triggers, num_pts,
        det.hdf1.num_capture, num_pts,
    )

    # Patch stage_sigs for pre-armed / software-trigger mode
    cam_keys = ("trigger_mode", "manual_trigger", "num_triggers", "wait_for_plugins")
    saved = {k: det.cam.stage_sigs[k] for k in cam_keys if k in det.cam.stage_sigs}
    det.cam.stage_sigs["trigger_mode"] = "Internal Enable"
    det.cam.stage_sigs["manual_trigger"] = "Enable"
    det.cam.stage_sigs["num_triggers"] = num_pts
    det.cam.stage_sigs["wait_for_plugins"] = "No"

    start_pos = motor.position
    positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
    pos_cache = {motor: None}

    def step(detectors, step, pos_cache):
        """Move motor, fire software trigger, wait for frame, read motor."""
        yield from bps.move_per_step(step, pos_cache)
        target = det.cam.array_counter.get() + 1
        det.cam.special_trigger_button.put(1, wait=False)
        # Wait for array_counter to confirm frame was acquired
        while det.cam.array_counter.get() < target:
            yield from bps.sleep(0.005)
        yield from bps.create("primary")
        yield from bps.read(motor)
        yield from bps.save()

    def inner():
        # Start pre-armed acquisition (accepts num_pts software triggers)
        det.cam.acquire.put(1)
        try:
            for pos in positions:
                yield from step([det], {motor: pos}, pos_cache)
        finally:
            # Wait for HDF to flush all frames before stage_wrapper unstages
            t0 = time.time()
            timeout = num_pts * count_time * 3 + 10
            while det.hdf1.num_captured.get() < num_pts:
                yield from bps.sleep(0.05)
                if time.time() - t0 > timeout:
                    print("WARNING: HDF write timeout — not all frames saved.")
                    break
            det.cam.acquire.put(0)
            yield from bps.mv(motor, start_pos)  # return to start position

    try:
        yield from bpp.stage_wrapper(
            bpp.run_wrapper(inner()),
            [det, motor],
        )
    finally:
        # Restore cam stage_sigs
        for k in cam_keys:
            if k in saved:
                det.cam.stage_sigs[k] = saved[k]
            else:
                det.cam.stage_sigs.pop(k, None)
        # Return detector to normal state: Internal Enable, manual trigger off
        det.cam.trigger_mode.put("Internal Enable")
        det.cam.manual_trigger.put("Disable")





