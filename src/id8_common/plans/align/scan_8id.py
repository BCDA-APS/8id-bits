"""
Scanning plans for the 8ID beamlines.

This module provides plans for scanning various motors and detectors at the
8ID -E and -I beamlines, including sample and with attenuation control.
"""

from typing import Optional

from matplotlib.pylab import det
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from bluesky import plans as bp
import bluesky.preprocessors as bpp
from ophyd import Device
from bluesky.plans import tune_centroid 
from bluesky.callbacks.fitting import PeakStats
from apstools.plans.alignment import lineup2
# from ..shutter_logic import *
# from id8_common.devices import softglue
from id8_common.plans.set.shutter_att import *
# from id8_common.plans.set.auto_filter import *
from id8_common.plans.acquire.ad_acq import *
# from ..nexus_acq_eiger_int import setup_eiger_int_series
import time
import numpy as np

huber = oregistry["huber"] # sample stage for 8-id-e
sample = oregistry["sample"] # sample stage for 8-id-i
pv_registers = oregistry["pv_registers"]
filter_beam = oregistry["filter_8ide"]
tetramm1 = oregistry["tetramm1"]
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]
softglue = oregistry["softglue"]
softglue_8id_acq = oregistry["softglue_8id_acq"]
filter_beam = oregistry["filter_8ide"]
rheometer = oregistry["rheometer"]


def auto_att(
    det,
    pilot_exptime: float = 0.05,
    rate_limit: float = 1e5,
    filter_factor: float = 5.0,
    retry_max: int = 10,
    grace_factor: float = 0.25,
):
    """Find the optimal attenuation using short pilot exposures.

    Args:
        det:            eiger4M or lambda2M
        pilot_exptime:  duration of each test frame (s)
        rate_limit:     max acceptable count rate (max pixel cts/s)
        filter_factor:  transmission multiplier per step, must be > 1
        retry_max:      max iterations before giving up
        grace_factor:   lower rate bound = rate_limit * grace_factor

    Example:

        auto_attenuate(eiger4M, pilot_exptime=0.05, rate_limit=4e5)
  """
    
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    low_rate = rate_limit * grace_factor

    orig_acq_time = det.cam.acquire_time.get()
    orig_acq_period = det.cam.acquire_period.get()

    det.cam.acquire_time.put(pilot_exptime)
    det.cam.acquire_period.put(pilot_exptime)

    if is_eiger:
        det.cam.trigger_mode.put("Internal Series")
        det.cam.num_images.put(1)
        det.cam.num_triggers.put(1)
        det.cam.manual_trigger.put("Disable")
    else:
        # lambda2M
        det.cam.trigger_mode.put("Internal")
        det.cam.num_images.put(1)

    det.stats1.enable.put(1)
    det.stats1.compute_statistics.put(1)

    # Start from maximum attenuation (minimum transmission) for safety
    filter_beam.transmission.move(1e-10)
    time.sleep(0.5)

    showbeam()
    try:
        for attempt in range(retry_max):
            det.cam.acquire.put(1)
            t0 = time.time()
            timeout = pilot_exptime * 5 + 2
            while det.cam.acquire.get() == 1:
                time.sleep(0.02)
                if time.time() - t0 > timeout:
                    print("  WARNING: pilot frame timed out")
                    break

            max_cts = det.stats1.max_value.get()
            rate = max_cts / pilot_exptime
            current_trans = filter_beam.transmission.readback.get()

            print(
                f"Attempt {attempt + 1}: trans={current_trans:.4f}"
                f"max_cts={max_cts:.0f}  rate={rate:.0f} cts/s"
            )

            if rate > rate_limit:
                new_trans = current_trans / filter_factor
                print(f"Rate too high. Reducing transmission to {new_trans:.6f}")
                filter_beam.transmission.move(new_trans)

            elif rate < low_rate:
                if rate > 0:
                    new_trans = current_trans * (0.75 * rate_limit / rate)
                else:
                    # rate=0: jump to a coarse fraction of max rather than tiny multiplier steps
                    new_trans = min(current_trans * 1000, 0.01)
                new_trans = min(new_trans, 1.0)
                if new_trans >= current_trans * 0.999:
                    if current_trans >= 0.999:
                        print("WARNING: already at max transmission, beam too weak.")
                        break
                print(f"Rate too low. Raising transmission to {new_trans:.4f}")
                filter_beam.transmission.move(new_trans)
            else:
                print(f"    Rate in [{low_rate:.0f}, {rate_limit:.0f}] cts/s -- converged.")
                break
        else:
            print(f"WARNING: auto_attenuate did not converge in {retry_max} attempts")
    finally:
        blockbeam()
        det.cam.acquire_time.put(orig_acq_time)
        det.cam.acquire_period.put(orig_acq_period)

    trans = filter_beam.transmission.readback.get()
    atten = filter_beam.attenuation.readback.get()
    print(f"  Final: transmission={trans:.4f}  attenuation={atten}")

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

# def save_images(det, save_img, num_pts, num_frames=1, file_path=None):
def save_images(det, save_img, num_pts, num_frames=1, file_path=None, folder_prefix=None):

    """
    Toggle saving (1) or not saving (0) images.

    det: detector instance (eiger4M, lambda2M, or tetramm)
    save_img: 1 save, 0 don't save
    num_frames: number of frames to capture per point (default 1)
    file_path: override base path (default uses a safe local path)
    """
    if save_img not in (0, 1):
        raise ValueError("save_img must be 1 or 0 (to save or not to save)")

    if save_img == 1:

        if folder_prefix is None:
            folder_prefix = gen_folder_prefix()

        if file_path is None:
            file_path = "/gdata/dm/8ID/8IDE/2026-1/comm202602/data/bluesky/"
        
        is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
        is_tetramm = "tetramm" in det.name.lower()

        def has(obj, attr):
            return getattr(obj, attr, None) is not None
        
        folder_prefix = gen_folder_prefix()
        file_name = folder_prefix
        file_header = folder_prefix

        print('Scan folder created: ' + folder_prefix)
        print('File path: ', file_path)

        if is_tetramm:
            yield from bps.mv(
                det.hdf1.file_path, file_path,
                det.hdf1.file_name, file_name,
                det.hdf1.num_capture, num_pts,
                det.hdf1.file_write_mode, 2,  # Stream mode
            )
            # yield from bps.mv(
            #     pv_registers.file_name, file_name,
            #     pv_registers.file_path, file_path,
            #     pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf",
            # )
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
            # yield from bps.mv(
                # pv_registers.file_name, file_name,
                # pv_registers.file_path, file_path,
                # pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf",
            # )
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
                setup_eiger_internal(
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
    
    folder_prefix = gen_folder_prefix() if save_img == 1 else ""
    md = {"Image file": folder_prefix, "detectors": [det.name], "motors": [motor.name], "plan_name": "dscan", "num_points": num_pts}
    yield from save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    if is_tetramm:
        # yield from save_images(det, save_img, num_pts)
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
                bpp.run_wrapper(inner_tetramm(), md=md),
                [det, motor],
            )
        finally:
            if save_img == 1:
                det.hdf1.capture.put(0)
                print('# images captured: ', det.cam.num_images.get())
        return

    if is_lambda:
        
        # yield from save_images(det, save_img, num_pts)

        yield from bps.mv(
            det.cam.operating_mode, 3, # 24-bit dual threshold mode
            det.cam.trigger_mode, "External_ImagePer",
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
            det.cam.num_images, num_pts,
            det.hdf1.num_capture, num_pts,
            softglue.num_triggers, 1,
            softglue.acq_time, count_time,
            softglue.acq_period, count_time,
        )

        start_pos = motor.position
        positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
        pos_cache = {motor: None}
        
        shutteron()
        showbeam()

        def step_lambda(step, pos_cache):
            yield from bps.move_per_step(step, pos_cache)
            softglue.start_pulses.put("1!")
            yield from bps.create("primary")
            yield from bps.read(motor)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
            yield from bps.save()

        def inner_lambda():
        # Start pre-armed acquisition (accepts num_pts softglue triggers)
            det.cam.acquire.put(1)
            det.hdf1.capture.put(1)
            try:
                for pos in positions:
                    yield from step_lambda({motor: pos}, pos_cache)

            finally:                
                softglue.stop_pulses.put("1!")
                det.cam.acquire.put(0)
                det.hdf1.capture.put(0)
                blockbeam()
                yield from bps.mv(motor, start_pos)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_lambda(), md=md), [motor],
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
        # yield from save_images(det, save_img, num_pts)

        # Correct settings after save_images (setup_eiger_int_series resets:
        # trigger_mode="Internal Series", num_triggers=1, manual_trigger="Disable",
        # hdf1.num_capture=1)
        yield from bps.mv(
            det.cam.trigger_mode, "Internal Series",
            det.cam.manual_trigger, "Enable",
            det.cam.num_images, 1,
            det.cam.num_triggers, num_pts,
            det.hdf1.num_capture, num_pts,
        )

        # Set stage_sigs for pre-armed software-trigger mode
        cam_keys = ("trigger_mode", "manual_trigger", "num_triggers", "wait_for_plugins")
        saved = {k: det.cam.stage_sigs[k] for k in cam_keys if k in det.cam.stage_sigs}
        det.cam.stage_sigs["trigger_mode"] = "Internal Series"
        det.cam.stage_sigs["manual_trigger"] = "Enable"
        det.cam.stage_sigs["num_triggers"] = num_pts
        det.cam.stage_sigs["wait_for_plugins"] = "No"

        start_pos = motor.position
        positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
        pos_cache = {motor: None}

        def step(detectors, step, pos_cache):
            """Move motor, fire software trigger, wait for frame, read motor."""
            yield from bps.move_per_step(step, pos_cache)
            det.cam.special_trigger_button.put(1, wait=False)
            yield from bps.create("primary")
            yield from bps.read(motor)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
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
                bpp.run_wrapper(inner(), md=md),
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
            print('# images captured: ', det.hdf1.num_captured.get())


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
    Two-motor relative scan with per-step triggering for eiger4M, lambda2M, and tetramm.

    usage:
        RE(d2scan(m1, s1, f1, m2, s2, f2, npts, ct))

    args:
        motor1, motor2: ophyd positioners (e.g. huber.eta, huber.chi)
        rel_begin1, rel_end1: relative start/end for motor1
        rel_begin2, rel_end2: relative start/end for motor2
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

    folder_prefix = gen_folder_prefix() if save_img == 1 else ""
    md = {"Image file": folder_prefix, "detectors": [det.name], "motors": [motor1.name, motor2.name], "plan_name": "d2scan", "num_points": num_pts}
    yield from save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    start1 = motor1.position
    start2 = motor2.position
    positions1 = np.linspace(start1 + rel_begin1, start1 + rel_end1, num_pts)
    positions2 = np.linspace(start2 + rel_begin2, start2 + rel_end2, num_pts)
    pos_cache = {motor1: None, motor2: None}

    if is_tetramm:
        # yield from save_images(det, save_img, num_pts)
        if save_img == 1:
            det.hdf1.enable.put(1)
            det.hdf1.capture.put(1)
            print(f"TetrAMM HDF capture armed: {det.hdf1.file_name.get()}")

        def inner_tetramm():
            try:
                for p1, p2 in zip(positions1, positions2):
                    yield from bps.move_per_step({motor1: p1, motor2: p2}, pos_cache)
                    yield from bps.trigger_and_read([det, motor1, motor2])
            finally:
                yield from bps.mv(motor1, start1, motor2, start2)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_tetramm(), md=md),
                [det, motor1, motor2],
            )
        finally:
            if save_img == 1:
                det.hdf1.capture.put(0)
                print("TetrAMM HDF capture stopped.")
        return

    if is_lambda:
        # yield from save_images(det, save_img, num_pts)
        yield from bps.mv(
            det.cam.operating_mode, 3,
            det.cam.trigger_mode, "External_ImagePer",
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
            det.cam.num_images, num_pts,
            det.hdf1.num_capture, num_pts,
            softglue.num_triggers, 1,
            softglue.acq_time, count_time,
            softglue.acq_period, count_time,
        )

        def step_lambda(step_dict, pc):
            yield from bps.move_per_step(step_dict, pc)
            softglue.start_pulses.put("1!")
            yield from bps.create("primary")
            yield from bps.read(motor1)
            yield from bps.read(motor2)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
            yield from bps.save()

        def inner_lambda():
            det.cam.acquire.put(1)
            det.hdf1.capture.put(1)
            try:
                for p1, p2 in zip(positions1, positions2):
                    yield from step_lambda({motor1: p1, motor2: p2}, pos_cache)
            finally:
                softglue.stop_pulses.put("1!")
                det.cam.acquire.put(0)
                det.hdf1.capture.put(0)
                blockbeam()
                yield from bps.mv(motor1, start1, motor2, start2)

        shutteron()
        showbeam()
        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_lambda(), md=md), [motor1, motor2],
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
        yield from bps.mv(
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
        )
        # yield from save_images(det, save_img, num_pts)
        yield from bps.mv(
            det.cam.trigger_mode, "Internal Series",
            det.cam.manual_trigger, "Enable",
            det.cam.num_images, 1,
            det.cam.num_triggers, num_pts,
            det.hdf1.num_capture, num_pts,
        )

        cam_keys = ("trigger_mode", "manual_trigger", "num_triggers", "wait_for_plugins")
        saved = {k: det.cam.stage_sigs[k] for k in cam_keys if k in det.cam.stage_sigs}
        det.cam.stage_sigs["trigger_mode"] = "Internal Series"
        det.cam.stage_sigs["manual_trigger"] = "Enable"
        det.cam.stage_sigs["num_triggers"] = num_pts
        det.cam.stage_sigs["wait_for_plugins"] = "No"

        def step(detectors, step_dict, pc):
            yield from bps.move_per_step(step_dict, pc)
            det.cam.special_trigger_button.put(1, wait=False)
            yield from bps.create("primary")
            yield from bps.read(motor1)
            yield from bps.read(motor2)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
            yield from bps.save()

        def inner():
            showbeam()
            det.cam.acquire.put(1)
            try:
                for p1, p2 in zip(positions1, positions2):
                    yield from step([det], {motor1: p1, motor2: p2}, pos_cache)
            finally:
                t0 = time.time()
                timeout = num_pts * count_time * 3 + 10
                while det.hdf1.num_captured.get() < num_pts:
                    yield from bps.sleep(0.05)
                    if time.time() - t0 > timeout:
                        print("WARNING: HDF write timeout — not all frames saved.")
                        break
                det.cam.acquire.put(0)
                blockbeam()
                yield from bps.mv(motor1, start1, motor2, start2)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner(), md = md),
                [det, motor1, motor2],
            )
        finally:
            for k in cam_keys:
                if k in saved:
                    det.cam.stage_sigs[k] = saved[k]
                else:
                    det.cam.stage_sigs.pop(k, None)
            det.cam.trigger_mode.put("Internal Enable")
            det.cam.manual_trigger.put("Disable")
            print('# images captured: ', det.cam.num_images.get())


def ascan(
    motor,
    abs_begin,
    abs_end,
    num_pts,
    count_time,
    det=eiger4M,
    att_ratio=7,
    save_img=1,
):
    """
    Absolute single-motor scan with per-step triggering for eiger4M, lambda2M, and tetramm.

    usage:
        RE(ascan(motor, abs_begin, abs_end, num_pts, count_time))

    args:
        motor: ophyd positioner (e.g. huber.eta, sample.y)
        abs_begin: absolute start position (motor units)
        abs_end: absolute end position (motor units)
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

    folder_prefix = gen_folder_prefix() if save_img == 1 else ""
    md = {"Image file": folder_prefix, "detectors": [det.name], "motors": [motor.name], "plan_name": "ascan", "num_points": num_pts}
    yield from save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    start_pos = motor.position
    positions = np.linspace(abs_begin, abs_end, num_pts)
    pos_cache = {motor: None}

    if is_tetramm:
        # yield from save_images(det, save_img, num_pts)
        if save_img == 1:
            det.hdf1.enable.put(1)
            det.hdf1.capture.put(1)
            print(f"TetrAMM HDF capture armed: {det.hdf1.file_name.get()}")

        def inner_tetramm():
            try:
                for pos in positions:
                    yield from bps.move_per_step({motor: pos}, pos_cache)
                    yield from bps.trigger_and_read([det, motor])
            finally:
                yield from bps.mv(motor, start_pos)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_tetramm(), md = md),
                [det, motor],
            )
        finally:
            if save_img == 1:
                det.hdf1.capture.put(0)
                print("TetrAMM HDF capture stopped.")
        return

    if is_lambda:
        # yield from save_images(det, save_img, num_pts)
        yield from bps.mv(
            det.cam.operating_mode, 3,
            det.cam.trigger_mode, "External_ImagePer",
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
            det.cam.num_images, num_pts,
            det.hdf1.num_capture, num_pts,
            softglue.num_triggers, 1,
            softglue.acq_time, count_time,
            softglue.acq_period, count_time,
        )

        def step_lambda(step_dict, pc):
            yield from bps.move_per_step(step_dict, pc)
            softglue.start_pulses.put("1!")
            yield from bps.create("primary")
            yield from bps.read(motor)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
            yield from bps.save()

        def inner_lambda():
            det.cam.acquire.put(1)
            det.hdf1.capture.put(1)
            try:
                for pos in positions:
                    yield from step_lambda({motor: pos}, pos_cache)
            finally:
                softglue.stop_pulses.put("1!")
                det.cam.acquire.put(0)
                det.hdf1.capture.put(0)
                blockbeam()
                yield from bps.mv(motor, start_pos)

        shutteron()
        showbeam()
        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_lambda(), md=md), [motor],
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
        yield from bps.mv(
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
        )
        # yield from save_images(det, save_img, num_pts)
        yield from bps.mv(
            det.cam.trigger_mode, "Internal Series",
            det.cam.manual_trigger, "Enable",
            det.cam.num_images, 1,
            det.cam.num_triggers, num_pts,
            det.hdf1.num_capture, num_pts,
        )

        cam_keys = ("trigger_mode", "manual_trigger", "num_triggers", "wait_for_plugins")
        saved = {k: det.cam.stage_sigs[k] for k in cam_keys if k in det.cam.stage_sigs}
        det.cam.stage_sigs["trigger_mode"] = "Internal Series"
        det.cam.stage_sigs["manual_trigger"] = "Enable"
        det.cam.stage_sigs["num_triggers"] = num_pts
        det.cam.stage_sigs["wait_for_plugins"] = "No"

        def step(detectors, step_dict, pc):
            yield from bps.move_per_step(step_dict, pc)
            det.cam.special_trigger_button.put(1, wait=False)
            yield from bps.create("primary")
            yield from bps.read(motor)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
            yield from bps.save()

        def inner():
            showbeam()
            det.cam.acquire.put(1)
            try:
                for pos in positions:
                    yield from step([det], {motor: pos}, pos_cache)
            finally:
                t0 = time.time()
                timeout = num_pts * count_time * 3 + 10
                while det.hdf1.num_captured.get() < num_pts:
                    yield from bps.sleep(0.05)
                    if time.time() - t0 > timeout:
                        print("WARNING: HDF write timeout — not all frames saved.")
                        break
                det.cam.acquire.put(0)
                blockbeam()
                yield from bps.mv(motor, start_pos)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner(), md = md),
                [det, motor],
            )
        finally:
            for k in cam_keys:
                if k in saved:
                    det.cam.stage_sigs[k] = saved[k]
                else:
                    det.cam.stage_sigs.pop(k, None)
            det.cam.trigger_mode.put("Internal Enable")
            det.cam.manual_trigger.put("Disable")
            print('# images captured: ', det.cam.num_images.get())


def a2scan(
    motor1,
    abs_begin1,
    abs_end1,
    motor2,
    abs_begin2,
    abs_end2,
    num_pts,
    count_time,
    det=eiger4M,
    att_ratio=7,
    save_img=1,
):
    """
    Absolute two-motor scan with per-step triggering for eiger4M, lambda2M, and tetramm.

    usage:
        RE(a2scan(m1, s1, f1, m2, s2, f2, npts, ct))

    args:
        motor1, motor2: ophyd positioners (e.g. huber.eta, huber.chi)
        abs_begin1, abs_end1: absolute start/end for motor1
        abs_begin2, abs_end2: absolute start/end for motor2
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

    folder_prefix = gen_folder_prefix() if save_img == 1 else ""
    md = {"Image file": folder_prefix, "detectors": [det.name], "motors": [motor1.name, motor2.name], "plan_name": "a2scan", "num_points": num_pts}
    yield from save_images(det, save_img, num_pts, folder_prefix=folder_prefix)
    
    start1 = motor1.position
    start2 = motor2.position
    positions1 = np.linspace(abs_begin1, abs_end1, num_pts)
    positions2 = np.linspace(abs_begin2, abs_end2, num_pts)
    pos_cache = {motor1: None, motor2: None}

    if is_tetramm:
        # yield from save_images(det, save_img, num_pts)
        if save_img == 1:
            det.hdf1.enable.put(1)
            det.hdf1.capture.put(1)
            print(f"TetrAMM HDF capture armed: {det.hdf1.file_name.get()}")

        def inner_tetramm():
            try:
                for p1, p2 in zip(positions1, positions2):
                    yield from bps.move_per_step({motor1: p1, motor2: p2}, pos_cache)
                    yield from bps.trigger_and_read([det, motor1, motor2])
            finally:
                yield from bps.mv(motor1, start1, motor2, start2)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_tetramm(), md = md),
                [det, motor1, motor2],
            )
        finally:
            if save_img == 1:
                det.hdf1.capture.put(0)
                print("TetrAMM HDF capture stopped.")
        return

    if is_lambda:
        # yield from save_images(det, save_img, num_pts)
        yield from bps.mv(
            det.cam.operating_mode, 3,
            det.cam.trigger_mode, "External_ImagePer",
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
            det.cam.num_images, num_pts,
            det.hdf1.num_capture, num_pts,
            softglue.num_triggers, 1,
            softglue.acq_time, count_time,
            softglue.acq_period, count_time,
        )

        def step_lambda(step_dict, pc):
            yield from bps.move_per_step(step_dict, pc)
            softglue.start_pulses.put("1!")
            yield from bps.create("primary")
            yield from bps.read(motor1)
            yield from bps.read(motor2)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
            yield from bps.save()

        def inner_lambda():
            det.cam.acquire.put(1)
            det.hdf1.capture.put(1)
            try:
                for p1, p2 in zip(positions1, positions2):
                    yield from step_lambda({motor1: p1, motor2: p2}, pos_cache)
            finally:
                softglue.stop_pulses.put("1!")
                det.cam.acquire.put(0)
                det.hdf1.capture.put(0)
                blockbeam()
                yield from bps.mv(motor1, start1, motor2, start2)

        shutteron()
        showbeam()
        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner_lambda(), md = md), [motor1, motor2],
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
        yield from bps.mv(
            det.cam.acquire_time, count_time,
            det.cam.acquire_period, count_time,
        )
        # yield from save_images(det, save_img, num_pts)
        yield from bps.mv(
            det.cam.trigger_mode, "Internal Series",
            det.cam.manual_trigger, "Enable",
            det.cam.num_images, 1,
            det.cam.num_triggers, num_pts,
            det.hdf1.num_capture, num_pts,
        )

        cam_keys = ("trigger_mode", "manual_trigger", "num_triggers", "wait_for_plugins")
        saved = {k: det.cam.stage_sigs[k] for k in cam_keys if k in det.cam.stage_sigs}
        det.cam.stage_sigs["trigger_mode"] = "Internal Series"
        det.cam.stage_sigs["manual_trigger"] = "Enable"
        det.cam.stage_sigs["num_triggers"] = num_pts
        det.cam.stage_sigs["wait_for_plugins"] = "No"

        def step(detectors, step_dict, pc):
            yield from bps.move_per_step(step_dict, pc)
            det.cam.special_trigger_button.put(1, wait=False)
            yield from bps.create("primary")
            yield from bps.read(motor1)
            yield from bps.read(motor2)
            yield from bps.read(det.stats1)
            yield from bps.read(det.stats2)
            yield from bps.read(det.stats3)
            yield from bps.save()

        def inner():
            showbeam()
            det.cam.acquire.put(1)
            try:
                for p1, p2 in zip(positions1, positions2):
                    yield from step([det], {motor1: p1, motor2: p2}, pos_cache)
            finally:
                t0 = time.time()
                timeout = num_pts * count_time * 3 + 10
                while det.hdf1.num_captured.get() < num_pts:
                    yield from bps.sleep(0.05)
                    if time.time() - t0 > timeout:
                        print("WARNING: HDF write timeout — not all frames saved.")
                        break
                det.cam.acquire.put(0)
                blockbeam()
                yield from bps.mv(motor1, start1, motor2, start2)

        try:
            yield from bpp.stage_wrapper(
                bpp.run_wrapper(inner(), md=md),
                [det, motor1, motor2],
            )
        finally:
            for k in cam_keys:
                if k in saved:
                    det.cam.stage_sigs[k] = saved[k]
                else:
                    det.cam.stage_sigs.pop(k, None)
            det.cam.trigger_mode.put("Internal Enable")
            det.cam.manual_trigger.put("Disable")
            print('# images captured: ', det.cam.num_images.get())

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



