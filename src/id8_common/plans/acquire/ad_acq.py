"""
Consolidated acquisition script for 8-ID detectors.

Supported modes:
    eiger4M       : "Internal Series", "External Enable"
    lambda2M      : "Internal"
    rigaku3M      : "ZDT", "ZDT4bit", "ZDT8bit"
    rigaku3M_epics: "EPICS"
"""

from datetime import datetime
import os
import time as ttime
import numpy as np

from apsbits.core.instrument_init import oregistry

from id8_common.utils.dm_util import dm_run_job
from id8_common.utils.dm_util import dm_setup
from id8_common.utils.nexus_utils import create_nexus_format_metadata

from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import shutteron
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.plans.set.shutter_att import post_align
# from id8_common.plans.set.shutter_att import att


pv_registers = oregistry["pv_registers"]

# =============================================================================
# General helpers
# =============================================================================

def get_connected_device(device_name):
    device = oregistry[device_name]

    if not device.connected:
        raise RuntimeError(f"{device_name} is not connected.")

    return device


def get_ophyd_object(name):
    parts = name.split(".")
    obj = oregistry[parts[0]]

    for part in parts[1:]:
        obj = getattr(obj, part)

    return obj


def get_sample_position_register(sample_index):
    register_name = f"sample{sample_index}_pos"
    return getattr(pv_registers, register_name)


def gen_folder_prefix():
    """
    Generate folder prefix from registers and current attenuation.

    Uses:
        pv_registers.header
        pv_registers.measurement_num
        pv_registers.sample_name
        filter_8ide.attenuation.readback

    Example:
        header = "A"
        measurement_num = 12
        sample_name = "G10"
        attenuation = 7

        returns "A0012_G10_a0007"

    The measurement number increments once per call.
    """
    filter_beam = get_connected_device("filter_8ide")

    header = pv_registers.header.get().strip()
    meas_num = int(pv_registers.measurement_num.get())
    sample_name = pv_registers.sample_name.get().strip()
    att_level = int(filter_beam.attenuation.readback.get())

    folder_prefix = f"{header}{meas_num:04d}_{sample_name}_a{att_level:04d}"

    pv_registers.measurement_num.put(meas_num + 1)

    return folder_prefix


def get_common_file_path(file_header, file_name):
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    mount_point = pv_registers.mount_point.get()
    use_subfolder = pv_registers.use_subfolder.get()

    if use_subfolder == "yes":
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == "no":
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_name}"
    else:
        raise ValueError("use_subfolder must be yes or no")

    return file_path


def get_rigaku_file_path(file_header, file_name):
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    mount_point = pv_registers.mount_point.get()
    use_subfolder = pv_registers.use_subfolder.get()

    if use_subfolder == "yes":
        file_path = f"{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == "no":
        file_path = f"{exp_name}/data/{file_name}"
    else:
        raise ValueError("use_subfolder must be yes or no")

    full_path = f"{mount_point}/{cycle_name}/{file_path}"

    return file_path, full_path


# =============================================================================
# Sample motion
# =============================================================================

def sample_mesh_move():
    """
    Move to the next point in a 2D mesh.

    Values read from registers:
        pv_registers.sample_move
        pv_registers.sample_index
        pv_registers.inner_motor
        pv_registers.outer_motor
        pv_registers.inner_center
        pv_registers.outer_center
        pv_registers.inner_range
        pv_registers.outer_range
        pv_registers.inner_pts
        pv_registers.outer_pts

    sample_move:
        "yes" -> move sample
        "no"  -> no sample motion and no position-register update

    sample_index:
        1 -> use pv_registers.sample1_pos
        2 -> use pv_registers.sample2_pos
        ...

    The position register stores the last used zero-based mesh index.
    The inner axis moves fastest.

    inner_range and outer_range are interpreted as full scan widths.
    """
    sample_move = pv_registers.sample_move.get().strip()

    if sample_move != "yes":
        return

    sample_index = int(pv_registers.sample_index.get())
    sample_position_register = get_sample_position_register(sample_index)
    last_index = int(sample_position_register.get())

    inner_motor_name = pv_registers.inner_motor.get().strip()
    outer_motor_name = pv_registers.outer_motor.get().strip()

    inner_center = float(pv_registers.inner_center.get())
    outer_center = float(pv_registers.outer_center.get())

    inner_range = float(pv_registers.inner_range.get())
    outer_range = float(pv_registers.outer_range.get())

    inner_pts = int(pv_registers.inner_pts.get())
    outer_pts = int(pv_registers.outer_pts.get())

    total_pts = inner_pts * outer_pts

    inner_positions = np.linspace(
        inner_center - inner_range / 2,
        inner_center + inner_range / 2,
        inner_pts,
    )

    outer_positions = np.linspace(
        outer_center - outer_range / 2,
        outer_center + outer_range / 2,
        outer_pts,
    )

    pos_index = (last_index + 1) % total_pts

    inner_index = pos_index % inner_pts
    outer_index = pos_index // inner_pts

    inner_pos = inner_positions[inner_index]
    outer_pos = outer_positions[outer_index]

    inner_motor = get_ophyd_object(inner_motor_name)
    outer_motor = get_ophyd_object(outer_motor_name)

    print(f"Moving {outer_motor.name} to {outer_pos}")
    outer_motor.move(outer_pos)

    print(f"Moving {inner_motor.name} to {inner_pos}")
    inner_motor.move(inner_pos)

    sample_position_register.put(pos_index)


# =============================================================================
# Detector setup functions
# =============================================================================

def setup_eiger_internal(acq_time, num_frames, file_header, file_name):
    eiger4M = get_connected_device("eiger4M")
    file_path = get_common_file_path(file_header, file_name)

    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_time)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)

    eiger4M.cam.trigger_mode.put("Internal Series")
    eiger4M.cam.num_images.put(num_frames)
    eiger4M.cam.num_triggers.put(1)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_eiger_external(acq_time, acq_period, num_frames, file_header, file_name):
    eiger4M = get_connected_device("eiger4M")
    softglue = get_connected_device("softglue")

    file_path = get_common_file_path(file_header, file_name)

    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_period)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)

    eiger4M.cam.num_triggers.put(num_frames)
    eiger4M.cam.trigger_mode.put("External Enable")

    softglue.acq_time.put(acq_time)
    softglue.acq_period.put(acq_period)
    softglue.num_triggers.put(num_frames)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_lambda_internal(acq_time, num_frames, file_header, file_name):
    lambda2M = get_connected_device("lambda2M")
    file_path = get_common_file_path(file_header, file_name)

    lambda2M.hdf1.enable.put(1)
    lambda2M.cam.trigger_mode.put("Internal")

    lambda2M.cam.acquire_time.put(acq_time)
    lambda2M.cam.acquire_period.put(acq_time)

    lambda2M.hdf1.file_name.put(file_name)
    lambda2M.hdf1.file_path.put(file_path)

    lambda2M.cam.num_images.put(num_frames)
    lambda2M.hdf1.num_capture.put(num_frames)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_rigaku_zdt(acq_time, num_frames, file_header, file_name):
    rigaku3M = get_connected_device("rigaku3M")
    softglue = get_connected_device("softglue")

    softglue.enable_rigaku.put('1')

    rigaku3M.cam.trigger_mode.put('Start with Trigger')

    file_path, full_path = get_rigaku_file_path(file_header, file_name)

    rigaku3M.cam.acquire_time.put(acq_time)
    rigaku3M.cam.acquire_period.put(acq_time)

    rigaku3M.cam.fast_file_name.put(f"{file_name}.bin")
    rigaku3M.cam.fast_file_path.put(file_path)

    rigaku3M.cam.num_images.put(num_frames)
    rigaku3M.cam.output_control.put("Sparsified")
    rigaku3M.cam.output_resolution.put("2 Bit")

    os.makedirs(full_path, mode=0o770, exist_ok=True)

    metadata_fname = f"{full_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_rigaku_zdt4bit(acq_time, num_frames, file_header, file_name):
    rigaku3M = get_connected_device("rigaku3M")
    softglue = get_connected_device("softglue")

    softglue.enable_rigaku.put('1')

    rigaku3M.cam.trigger_mode.put('Start with Trigger')

    file_path, full_path = get_rigaku_file_path(file_header, file_name)

    rigaku3M.cam.acquire_time.put(acq_time)
    rigaku3M.cam.acquire_period.put(acq_time)

    rigaku3M.cam.fast_file_name.put(f"{file_name}.bin")
    rigaku3M.cam.fast_file_path.put(file_path)

    rigaku3M.cam.num_images.put(num_frames)
    rigaku3M.cam.output_control.put("Sparsified")
    rigaku3M.cam.output_resolution.put("4 Bit")

    os.makedirs(full_path, mode=0o770, exist_ok=True)

    metadata_fname = f"{full_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_rigaku_zdt8bit(acq_time, num_frames, file_header, file_name):
    rigaku3M = get_connected_device("rigaku3M")
    softglue = get_connected_device("softglue")

    softglue.enable_rigaku.put('1')

    rigaku3M.cam.trigger_mode.put('Start with Trigger')

    file_path, full_path = get_rigaku_file_path(file_header, file_name)

    rigaku3M.cam.acquire_time.put(acq_time)
    rigaku3M.cam.acquire_period.put(acq_time)

    rigaku3M.cam.fast_file_name.put(f"{file_name}.bin")
    rigaku3M.cam.fast_file_path.put(file_path)

    rigaku3M.cam.num_images.put(num_frames)
    rigaku3M.cam.output_control.put("Sparsified")
    rigaku3M.cam.output_resolution.put("8 Bit")

    os.makedirs(full_path, mode=0o770, exist_ok=True)

    metadata_fname = f"{full_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_rigaku_epics(acq_time, num_frames, file_header, file_name):
    rigaku3M = get_connected_device("rigaku3M")
    softglue = get_connected_device("softglue")

    softglue.enable_rigaku.put('1')

    rigaku3M.cam.trigger_mode.put('Start with Trigger')

    file_path, full_path = get_rigaku_file_path(file_header, file_name)

    rigaku3M.cam.acquire_time.put(acq_time)
    rigaku3M.cam.acquire_period.put(acq_time)

    rigaku3M.cam.fast_file_name.put(f"{file_name}.bin")
    rigaku3M.cam.fast_file_path.put(file_path)

    rigaku3M.cam.num_images.put(num_frames)
    rigaku3M.cam.output_control.put("areaDetector")
    rigaku3M.cam.output_resolution.put("16 Bit")

    rigaku3M.hdf1.file_name.put(file_name)
    rigaku3M.hdf1.file_path.put(full_path)
    rigaku3M.hdf1.num_capture.put(num_frames)

    os.makedirs(full_path, mode=0o770, exist_ok=True)

    metadata_fname = f"{full_path}/{file_name}_metadata.hdf"

    return metadata_fname



# =============================================================================
# Detector acquire functions
# =============================================================================

def acquire_eiger_internal():
    eiger4M = get_connected_device("eiger4M")

    shutteroff()
    showbeam()
    ttime.sleep(0.1)

    eiger4M.hdf1.capture.put(1)
    eiger4M.cam.acquire.put(1)
    ttime.sleep(0.5)

    while eiger4M.cam.acquire.get() == 1:
        ttime.sleep(0.1)

    blockbeam()

    while eiger4M.hdf1.capture.get() == 1:
        ttime.sleep(0.1)


def acquire_eiger_external():
    eiger4M = get_connected_device("eiger4M")
    softglue = get_connected_device("softglue")

    shutteron()
    showbeam()
    ttime.sleep(0.1)

    eiger4M.hdf1.capture.put(1)
    eiger4M.cam.acquire.put(1)
    ttime.sleep(1.0)

    softglue.start_pulses.put("1!")

    while True:
        #### QZ on 2026/01/06 ####
        # without the 0.5 s wait time, the repeating acqs go out of sync. 
        # Don't know why and maybe the 0.5 s can be made shorter
        #### QZ on 2026/01/06 ####
        ttime.sleep(0.5)
        det_status = eiger4M.cam.acquire_busy.get()
        if det_status == 1:
            ttime.sleep(0.1)
        if det_status == 0:
            break
    blockbeam()

    frame_num_set = eiger4M.hdf1.queue_size.get()
    count = 0
    while count < 100:  # 100 and 0.1 s are some empirical number for timeout conditions
        frame_num_processed = eiger4M.hdf1.queue_free.get()
        if frame_num_processed == frame_num_set:
            break
        else:
            ttime.sleep(0.1)
            count = +1
        eiger4M.hdf1.capture.put(0)
    
    softglue.pv_clear1.put("1!")
    softglue.pv_clear2.put("1!")
    

def acquire_lambda_internal():
    lambda2M = get_connected_device("lambda2M")

    lambda2M.cam.acquire.put(0)
    ttime.sleep(1.0)

    showbeam()

    lambda2M.hdf1.capture.put(1)
    lambda2M.cam.acquire.put(1)

    while lambda2M.cam.acquire.get() == 1:
        ttime.sleep(0.05)

    blockbeam()

    while lambda2M.hdf1.capture.get() == 1:
        ttime.sleep(0.05)


def acquire_rigaku_zdt():
    rigaku3M = get_connected_device("rigaku3M")

    showbeam()
    ttime.sleep(0.1)

    rigaku3M.cam.acquire.put(1)

    while rigaku3M.cam.detector_state.get() != 1:
        ttime.sleep(0.1)

    while rigaku3M.cam.detector_state.get() != 0:
        ttime.sleep(0.1)

    blockbeam()


def acquire_rigaku_epics():
    rigaku3M = get_connected_device("rigaku3M")

    showbeam()
    ttime.sleep(0.1)

    rigaku3M.hdf1.capture.put(1)
    rigaku3M.cam.acquire.put(1)

    while rigaku3M.cam.detector_state.get() != 1:
        ttime.sleep(0.1)

    while rigaku3M.cam.detector_state.get() != 0:
        ttime.sleep(0.1)

    blockbeam()

    while rigaku3M.hdf1.capture.get() == 1:
        ttime.sleep(0.1)


# =============================================================================
# Mode table
# =============================================================================

ACQ_MODES = {
    "eiger4M": {
        "Internal Series": {
            "setup": setup_eiger_internal,
            "acquire": acquire_eiger_internal,
            "needs_acq_period": False,
            "required_devices": ["eiger4M"],
        },
        "External Enable": {
            "setup": setup_eiger_external,
            "acquire": acquire_eiger_external,
            "needs_acq_period": True,
            "required_devices": ["eiger4M", "softglue"],
        },
    },

    "lambda2M": {
        "Internal": {
            "setup": setup_lambda_internal,
            "acquire": acquire_lambda_internal,
            "needs_acq_period": False,
            "required_devices": ["lambda2M"],
        },
    },

    "rigaku3M": {
        "ZDT": {
            "setup": setup_rigaku_zdt,
            "acquire": acquire_rigaku_zdt,
            "needs_acq_period": False,
            "required_devices": ["rigaku3M"],
            "min_acq_time": 20e-6,
        },
        "ZDT4bit": {
            "setup": setup_rigaku_zdt4bit,
            "acquire": acquire_rigaku_zdt,
            "needs_acq_period": False,
            "required_devices": ["rigaku3M"],
            "min_acq_time": 40e-6,
        },
        "ZDT8bit": {
            "setup": setup_rigaku_zdt8bit,
            "acquire": acquire_rigaku_zdt,
            "needs_acq_period": False,
            "required_devices": ["rigaku3M"],
            "min_acq_time": 80e-6,
        },
    },

    "rigaku3M_epics": {
        "EPICS": {
            "setup": setup_rigaku_epics,
            "acquire": acquire_rigaku_epics,
            "needs_acq_period": False,
            "required_devices": ["rigaku3M"],
            "hardware_device": "rigaku3M",
        },
    },
}


# =============================================================================
# Cleanup helper
# =============================================================================

def cleanup_acquisition(det=None, mode_info=None):
    """Stop softglue (if used by this mode) and abort detector acquisition."""
    if mode_info is not None and "softglue" in mode_info.get("required_devices", []):
        try:
            softglue = get_connected_device("softglue")
            softglue.stop_pulses.put("1!")
        except Exception:
            pass

    if det is not None:
        try:
            det.cam.acquire.put(0)
        except Exception:
            pass
        if hasattr(det, "hdf1"):
            try:
                det.hdf1.capture.put(0)
            except Exception:
                pass


# =============================================================================
# Main user-facing acquisition function
# =============================================================================

def det_acq_series(wait_time=0):
    """
    Run repeated detector acquisitions.

    Values read from registers:
        pv_registers.det_name
        pv_registers.det_mode
        pv_registers.header
        pv_registers.sample_name
        pv_registers.measurement_num
        pv_registers.acq_time
        pv_registers.acq_period
        pv_registers.num_frames
        pv_registers.num_repeats

    Sample motion values read from registers:
        pv_registers.sample_move
        pv_registers.sample_index
        pv_registers.inner_motor
        pv_registers.outer_motor
        pv_registers.inner_center
        pv_registers.outer_center
        pv_registers.inner_range
        pv_registers.outer_range
        pv_registers.inner_pts
        pv_registers.outer_pts

    sample_move:
        "Yes" -> sample moves before each repeat
        "No"  -> sample does not move

    wait_time:
        Wait time before each repeated acquisition.

    File naming:
        gen_folder_prefix() generates:
            A0012_G10_a0007

        det_acq_series() appends:
            _f001000
            _r00001
    """
    det = None
    mode_info = None
    try:
        post_align()
        shutteroff()

        workflowProcApi, dmuser = dm_setup()

        detector = pv_registers.det_name.get().strip()
        mode = pv_registers.det_mode.get().strip()

        acq_time = float(pv_registers.acq_time.get())
        acq_period = float(pv_registers.acq_period.get())
        num_frames = int(pv_registers.num_frames.get())
        num_reps = int(pv_registers.num_repeats.get())

        mode_info = ACQ_MODES[detector][mode]

        for device_name in mode_info["required_devices"]:
            get_connected_device(device_name)

        det = get_connected_device(mode_info.get("hardware_device", detector))
        setup_func = mode_info["setup"]
        acquire_func = mode_info["acquire"]

        folder_prefix = gen_folder_prefix()
        file_header = f"{folder_prefix}_f{num_frames:06d}"

        for rep in range(num_reps):
            ttime.sleep(wait_time)

            sample_mesh_move()

            file_name = f"{file_header}_r{rep + 1:05d}"

            if mode_info["needs_acq_period"]:
                metadata_fname = setup_func(
                    acq_time=acq_time,
                    acq_period=acq_period,
                    num_frames=num_frames,
                    file_header=file_header,
                    file_name=file_name,
                )
            else:
                metadata_fname = setup_func(
                    acq_time=acq_time,
                    num_frames=num_frames,
                    file_header=file_header,
                    file_name=file_name,
                )

            time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"\n{time_now}, Starting measurement {file_name}")

            acquire_func()

            time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{time_now}, Complete measurement {file_name}")

            create_nexus_format_metadata(metadata_fname, det=det)

            dm_run_job(workflowProcApi, dmuser, file_name)

    except KeyboardInterrupt:
        cleanup_acquisition(det, mode_info)
        raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).")
    except Exception as e:
        print(f"Error occurred during measurement: {e}")
    finally:
        pass


# =============================================================================
# Example use cases
# =============================================================================

# Before running, populate the registers. Your current register device defines:
#     det_name, det_mode,
#     header, sample_name, measurement_num,
#     acq_time, acq_period, num_frames, num_repeats,
#     sample_move,
#     sample_index,
#     inner_motor, outer_motor,
#     inner_center, outer_center,
#     inner_range, outer_range,
#     inner_pts, outer_pts.


# -----------------------------------------------------------------------------
# Example 1: Eiger internal, no sample motion
# -----------------------------------------------------------------------------
#
# pv_registers.header.put("A")
# pv_registers.sample_name.put("G10")
# pv_registers.measurement_num.put(12)
#
# pv_registers.det_name.put("eiger4M")
# pv_registers.det_mode.put("Internal Series")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(0.1)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(5)
#
# pv_registers.sample_move.put("No")
#
# det_acq_series(wait_time=0)


# -----------------------------------------------------------------------------
# Example 2: Eiger internal, fresh sample spot before every repeat
# -----------------------------------------------------------------------------
#
# pv_registers.header.put("A")
# pv_registers.sample_name.put("G10")
# pv_registers.measurement_num.put(12)
#
# pv_registers.det_name.put("eiger4M")
# pv_registers.det_mode.put("Internal Series")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(0.1)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(20)
#
# pv_registers.sample_move.put("Yes")
# pv_registers.sample_index.put(1)
# pv_registers.sample1_pos.put(-1)
#
# pv_registers.inner_motor.put("sample.x")
# pv_registers.outer_motor.put("sample.y")
# pv_registers.inner_center.put(0)
# pv_registers.outer_center.put(0)
# pv_registers.inner_range.put(2)
# pv_registers.outer_range.put(2)
# pv_registers.inner_pts.put(21)
# pv_registers.outer_pts.put(21)
#
# det_acq_series(wait_time=0)


# -----------------------------------------------------------------------------
# Example 3: Eiger external enable mode
# -----------------------------------------------------------------------------
#
# pv_registers.header.put("A")
# pv_registers.sample_name.put("G10")
# pv_registers.measurement_num.put(12)
#
# pv_registers.det_name.put("eiger4M")
# pv_registers.det_mode.put("External Enable")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(1.0)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(10)
#
# pv_registers.sample_move.put("Yes")
# pv_registers.sample_index.put(1)
# pv_registers.sample1_pos.put(-1)
#
# pv_registers.inner_motor.put("sample.x")
# pv_registers.outer_motor.put("sample.y")
# pv_registers.inner_center.put(0)
# pv_registers.outer_center.put(0)
# pv_registers.inner_range.put(2)
# pv_registers.outer_range.put(2)
# pv_registers.inner_pts.put(21)
# pv_registers.outer_pts.put(21)
#
# det_acq_series(wait_time=0)


# -----------------------------------------------------------------------------
# Example 4: Lambda internal
# -----------------------------------------------------------------------------
#
# pv_registers.header.put("A")
# pv_registers.sample_name.put("G10")
# pv_registers.measurement_num.put(12)
#
# pv_registers.det_name.put("lambda2M")
# pv_registers.det_mode.put("Internal")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(0.1)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(20)
#
# pv_registers.sample_move.put("Yes")
# pv_registers.sample_index.put(1)
# pv_registers.sample1_pos.put(-1)
#
# pv_registers.inner_motor.put("sample.x")
# pv_registers.outer_motor.put("sample.y")
# pv_registers.inner_center.put(0)
# pv_registers.outer_center.put(0)
# pv_registers.inner_range.put(2)
# pv_registers.outer_range.put(0)
# pv_registers.inner_pts.put(21)
# pv_registers.outer_pts.put(1)
#
# det_acq_series(wait_time=0)


# -----------------------------------------------------------------------------
# Example 5: Rigaku ZDT, no sample motion
# -----------------------------------------------------------------------------
#
# pv_registers.header.put("A")
# pv_registers.sample_name.put("G10")
# pv_registers.measurement_num.put(12)
#
# pv_registers.det_name.put("rigaku3M")
# pv_registers.det_mode.put("ZDT")
# pv_registers.acq_time.put(2e-5)
# pv_registers.acq_period.put(2e-5)
# pv_registers.num_frames.put(100000)
# pv_registers.num_repeats.put(5)
#
# pv_registers.sample_move.put("No")
#
# det_acq_series(wait_time=0)


# -----------------------------------------------------------------------------
# Example 6: Rigaku EPICS (slow, HDF5 output via areaDetector plugin)
# -----------------------------------------------------------------------------
#
# pv_registers.header.put("A")
# pv_registers.sample_name.put("G10")
# pv_registers.measurement_num.put(12)
#
# pv_registers.det_name.put("rigaku3M_epics")
# pv_registers.det_mode.put("EPICS")
# pv_registers.acq_time.put(0.001)
# pv_registers.acq_period.put(0.001)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(5)
#
# pv_registers.sample_move.put("No")
#
# det_acq_series(wait_time=0)