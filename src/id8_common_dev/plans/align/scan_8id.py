"""
Scanning plans for the 8ID beamlines.

This module provides plans for scanning various motors and detectors at the
8ID -E and -I beamlines, including sample and with attenuation control.
"""

from typing import Optional
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from bluesky import plans as bp
import bluesky.preprocessors as bpp
from ophyd import Device
from bluesky.plans import tune_centroid 
from bluesky.callbacks.fitting import PeakStats
from apstools.plans.alignment import lineup2
from ..shutter_logic import *
from ..sample_info_unpack import gen_folder_prefix
from ..nexus_acq_eiger_int import setup_eiger_int_series
# from ..devices.area_detector import ad_setup
# from apstools.devices import AD_plugin_primed, AD_prime_plugin2
import time
import numpy as np

huber = oregistry["huber"] # sample stage for 8-id-e
sample = oregistry["sample"] # sample stage for 8-id-i
pv_registers = oregistry["pv_registers"]
filter_beam = oregistry["filter_8ide"]
tetramm1 = oregistry["tetramm1"]
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]
softglue_8idi = oregistry["softglue_8idi"]
softglue_8id_acq = oregistry["softglue_8id_acq"]
filter_beam = oregistry["filter_8ide"]
rheometer = oregistry["rheometer"]


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
        file_path = "/gdata/dm/8ID/8IDE/2026-1/comm202602/data/bluesky/"
        # file_path = "/home/beams/8IDIUSER/sdmarks/tetramm_testing"
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_tetramm = "tetramm" in det.name.lower()

    def has(obj, attr):
        return getattr(obj, attr, None) is not None

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
                pv_registers.file_path, file_path,
                pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf",
            )
            return  # capture managed manually in dscan

        # for eiger4m, lambda2m 
        if has(det, "cam"):
            if has(det.cam, "fw_enable"):
                yield from bps.mv(det.cam.fw_enable, 1)
            if has(det.cam, "save_files"):
                yield from bps.mv(det.cam.save_files, 1)

        if has(det, "hdf1"):
            if has(det.hdf1, "num_capture"):
                yield from bps.mv(det.hdf1.num_capture, num_pts)
            yield from bps.mv(
                pv_registers.file_name, file_name,
                pv_registers.file_path, file_path,
                pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf",
            )
            if has(det.hdf1, "file_name"):
                yield from bps.mv(det.hdf1.file_name, file_name)
            if has(det.hdf1, "file_path"):
                yield from bps.mv(det.hdf1.file_path, file_path)

        if has(det, "cam"):
            det.cam.trigger_mode.put(0)
            if is_eiger:
                if has(det.cam, "trigger_mode"):
                    yield from bps.mv(det.cam.trigger_mode, "Internal Enable")
                if has(det.cam, "num_triggers"):
                    yield from bps.mv(det.cam.num_triggers, 1)
                if has(det.cam, "num_images"):
                    yield from bps.mv(det.cam.num_images, 1)
                setup_eiger_int_series(
                    acq_time=float(det.cam.acquire_time.get()) if has(det.cam, "acquire_time") else 1.0,
                    num_frames=num_frames,
                    file_header=file_header,
                    file_name=file_name,
                )
            else:
                if has(det.cam, "num_images"):
                    yield from bps.mv(det.cam.num_images, num_frames)

def dscan(motor, rel_begin, rel_end, num_pts, count_time,
               det=eiger4M, att_ratio=1e6, save_img=1):
    """
    Pre-armed software-trigger scan for fast acquisitions.

    For Eiger4M: Arms detector for all N frames upfront, then sends a software
    trigger at each point via det.cam.special_trigger_button.

    For TetrAMM: No pre-arming (continuous mode). Steps through positions and
    reads the current value at each point. HDF capture is armed/stopped manually.

    For Lambda2M: under construction 

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
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

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

    # lambda2M 
    # if is_lambda:

    #     yield from bps.mv(
    #         det.cam.acquire_time, count_time,
    #         det.cam.acquire_period, count_time,
    #         det.cam.num_images, 1,
    #     )
    #     yield from save_images(det, save_img, num_pts)

    #     start_pos = motor.position
    #     positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
    #     pos_cache = {motor: None}

    #     def inner_area_det():
    #         try:
    #             for pos in positions:
    #                 yield from bps.move_per_step({motor: pos}, pos_cache)
    #                 yield from bps.trigger_and_read([det, motor])
    #         finally:
    #             yield from bps.mv(motor, start_pos)

    #     yield from bpp.stage_wrapper(
    #         bpp.run_wrapper(inner_area_det()),
    #         [det, motor],
    #     )
    #     return

    if is_lambda:
        
        yield from bps.mv(
            det.cam.operating_mode, 3, # 24-bit dual threshold mode
            det.cam.trigger_mode, "External_ImagePer",
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
            det.cam.num_images, num_pts,
            det.hdf1.num_capture, num_pts,
        )
        
        yield from save_images(det, save_img, num_pts)
        # print(det.cam.num_triggers.get())

        start_pos = motor.position
        positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
        pos_cache = {motor: None}

        def step_lambda(detectors, step, pos_cache, frame_num):
            yield from bps.move_per_step(step, pos_cache)
            target = det.cam.num_images.get() + 1
            softglue_8idi.start_pulses.put("1!") 
            while det.hdf1.num_capture.get() < target:
                yield from bps.sleep(0.005)
            yield from bps.create("primary")
            yield from bps.read(motor)            
            yield from bps.save()

        def inner_lambda():            
        # Start pre-armed acquisition (accepts num_pts softglue triggers)
            shutteron()
            showbeam()
            try:
                for ii, pos in enumerate(positions):
                    yield from step_lambda([det], {motor: pos}, pos_cache, ii + 1)
            finally:                
                softglue_8idi.stop_pulses.put("1!")
                det.cam.acquire.put(0)
                det.hdf1.capture.put(0)
                blockbeam()
                yield from bps.mv(motor, start_pos)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_lambda()),
                [motor], 
            )

        finally:
            det.cam.operating_mode.put(3)
            det.cam.trigger_mode.put(0)
            softglue_8id_acq.preset.put(50)
            blockbeam()
            shutteroff()
            print('# images captured: ', det.cam.num_images.get())
        return

    # eiger4M
    if is_eiger:
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

        # Patch stage_sigs for pre-armed software-trigger mode
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
    pre_align()
    att(att_ratio)
    showbeam()

    yield from bp.rel_scan([det], sample.x, rel_begin, rel_end, num_pts)
    blockbeam()


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
    pre_align()
    att(att_ratio)

    showbeam()
    yield from bp.rel_scan([det], sample.y, rel_begin, rel_end, num_pts)
    blockbeam()


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
    pre_align()
    att(att_ratio)

    showbeam()
    yield from bp.rel_scan([det], rheometer.x, rel_begin, rel_end, num_pts)
    blockbeam()


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
    pre_align()
    att(att_ratio)

    showbeam()
    yield from bp.rel_scan([det], rheometer.y, rel_begin, rel_end, num_pts)
    blockbeam()


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
    pre_align()
    att(att_ratio)

    rheometer.x.put(14.0)
    showbeam()
    yield from bp.rel_scan([det], rheometer.x, -0.5, 0.5, 100)
    blockbeam()

    rheometer.x.put(-2.6)
    showbeam()
    yield from bp.rel_scan([det], rheometer.x, -0.5, 0.5, 100)
    blockbeam()



