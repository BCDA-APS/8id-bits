"""
Consolidated acquisition script for 8-ID detectors.

Supported modes:
    eiger4M  : "Internal Series", "External Enable"
    lambda2M : "Internal"
    rigaku3M : "ZDT"
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
from id8_common.plans.set.shutter_att import att


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


def get_sample_position_register(move_index_register):
    register_name = f"sample{move_index_register}_pos"
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

    header = pv_registers.header.get()
    meas_num = int(pv_registers.measurement_num.get())
    sample_name = pv_registers.sample_name.get()
    att_level = int(filter_beam.attenuation.readback.get())

    folder_prefix = f"{header}{meas_num:04d}_{sample_name}_a{att_level:04d}"

    pv_registers.measurement_num.put(meas_num + 1)

    return folder_prefix


def get_common_file_path(file_header, file_name):
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    mount_point = pv_registers.mount_point.get()
    use_subfolder = pv_registers.use_subfolder.get()

    if use_subfolder == "Yes":
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == "No":
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_name}"
    else:
        raise ValueError("use_subfolder must be Yes or No")

    return file_path


def get_rigaku_file_path(file_header, file_name):
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    mount_point = pv_registers.mount_point.get()
    use_subfolder = pv_registers.use_subfolder.get()

    if use_subfolder == "Yes":
        file_path = f"{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == "No":
        file_path = f"{exp_name}/data/{file_name}"
    else:
        raise ValueError("use_subfolder must be Yes or No")

    full_path = f"{mount_point}/{cycle_name}/{file_path}"

    return file_path, full_path


# =============================================================================
# Sample motion
# =============================================================================

def sample_mesh_move(axis_inner, axis_outer, move_index_register=0):
    """
    Move to the next point in a 2D mesh.

    move_index_register:
        0 -> no sample motion and no register update
        1 -> use pv_registers.sample1_pos
        2 -> use pv_registers.sample2_pos

    The register stores the last used zero-based mesh index.
    The inner axis moves fastest.
    """
    if move_index_register == 0:
        return

    sample_position_register = get_sample_position_register(move_index_register)
    last_index = int(sample_position_register.get())

    inner_pts = int(axis_inner["pts"])
    outer_pts = int(axis_outer["pts"])
    total_pts = inner_pts * outer_pts

    inner_positions = np.linspace(axis_inner["min"], axis_inner["max"], inner_pts)
    outer_positions = np.linspace(axis_outer["min"], axis_outer["max"], outer_pts)

    pos_index = (last_index + 1) % total_pts

    inner_index = pos_index % inner_pts
    outer_index = pos_index // inner_pts

    inner_pos = inner_positions[inner_index]
    outer_pos = outer_positions[outer_index]

    inner_motor = get_ophyd_object(axis_inner["motor"])
    outer_motor = get_ophyd_object(axis_outer["motor"])

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
    softglue_8idi = get_connected_device("softglue_8idi")

    file_path = get_common_file_path(file_header, file_name)

    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_period)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)

    eiger4M.cam.num_triggers.put(num_frames)
    eiger4M.cam.trigger_mode.put("External Enable")

    softglue_8idi.acq_time.put(acq_time)
    softglue_8idi.acq_period.put(acq_period)
    softglue_8idi.num_triggers.put(num_frames)

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

    while eiger4M.cam.acquire.get() == 1:
        ttime.sleep(0.1)

    blockbeam()

    while eiger4M.hdf1.capture.get() == 1:
        ttime.sleep(0.1)


def acquire_eiger_external():
    eiger4M = get_connected_device("eiger4M")
    softglue_8idi = get_connected_device("softglue_8idi")

    shutteron()
    showbeam()
    ttime.sleep(0.1)

    eiger4M.hdf1.capture.put(1)
    eiger4M.cam.acquire.put(1)

    softglue_8idi.start_pulses.put("1!")

    while eiger4M.cam.acquire_busy.get() == 1:
        ttime.sleep(0.1)

    blockbeam()
    shutteroff()

    while eiger4M.hdf1.capture.get() == 1:
        ttime.sleep(0.1)


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
            "required_devices": ["eiger4M", "softglue_8idi"],
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
        },
    },
}


# =============================================================================
# Main user-facing acquisition function
# =============================================================================

def det_acq_series(
    wait_time=0,
    move_index_register=0,
    axis_inner=None,
    axis_outer=None,
):
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

    wait_time:
        Wait time before each repeated acquisition.

    move_index_register:
        0 means no sample motion.
        Nonzero integer N means use pv_registers.sampleN_pos.

    File naming:
        gen_folder_prefix() generates:
            A0012_G10_a0007

        det_acq_series() appends:
            _f001000
            _r00001
    """
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

    if mode_info["needs_acq_period"] and acq_period is None:
        raise ValueError("This mode requires acq_period.")

    for device_name in mode_info["required_devices"]:
        get_connected_device(device_name)

    det = get_connected_device(detector)
    setup_func = mode_info["setup"]
    acquire_func = mode_info["acquire"]

    folder_prefix = gen_folder_prefix()
    file_header = f"{folder_prefix}_f{num_frames:06d}"

    for rep in range(num_reps):
        ttime.sleep(wait_time)

        sample_mesh_move(
            axis_inner=axis_inner,
            axis_outer=axis_outer,
            move_index_register=move_index_register,
        )

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


# =============================================================================
# Example use cases
# =============================================================================

# Before running, populate the registers. Your current register device defines:
#     det_name, det_mode, header, sample_name,
#     measurement_num, acq_time, acq_period, num_frames, num_repeats.
#
# Example register setup:
#
# pv_registers.header.put("A")
# pv_registers.sample_name.put("G10")
# pv_registers.measurement_num.put(12)
# pv_registers.det_name.put("eiger4M")
# pv_registers.det_mode.put("Internal Series")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(0.1)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(5)
#
# With att_level=7, the file prefix becomes:
#     A0012_G10_a0007_f001000


# -----------------------------------------------------------------------------
# Example 1: Eiger internal, no sample motion
# -----------------------------------------------------------------------------
#
# pv_registers.det_name.put("eiger4M")
# pv_registers.det_mode.put("Internal Series")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(0.1)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(5)
#
# det_acq_series(
#     att_level=7,
#     wait_time=0,
#     move_index_register=0,
#     axis_inner={
#         "motor": "sample.x",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
#     axis_outer={
#         "motor": "sample.y",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
# )


# -----------------------------------------------------------------------------
# Example 2: Eiger internal, fresh sample spot before every repeat
# -----------------------------------------------------------------------------
#
# pv_registers.det_name.put("eiger4M")
# pv_registers.det_mode.put("Internal Series")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(0.1)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(20)
#
# det_acq_series(
#     att_level=7,
#     wait_time=0,
#     move_index_register=1,
#     axis_inner={
#         "motor": "sample.x",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
#     axis_outer={
#         "motor": "sample.y",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
# )


# -----------------------------------------------------------------------------
# Example 3: Eiger external enable mode
# -----------------------------------------------------------------------------
#
# pv_registers.det_name.put("eiger4M")
# pv_registers.det_mode.put("External Enable")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(1.0)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(10)
#
# det_acq_series(
#     att_level=7,
#     wait_time=0,
#     move_index_register=1,
#     axis_inner={
#         "motor": "sample.x",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
#     axis_outer={
#         "motor": "sample.y",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
# )


# -----------------------------------------------------------------------------
# Example 4: Lambda internal
# -----------------------------------------------------------------------------
#
# pv_registers.det_name.put("lambda2M")
# pv_registers.det_mode.put("Internal")
# pv_registers.acq_time.put(0.1)
# pv_registers.acq_period.put(0.1)
# pv_registers.num_frames.put(1000)
# pv_registers.num_repeats.put(20)
#
# det_acq_series(
#     att_level=7,
#     wait_time=0,
#     move_index_register=1,
#     axis_inner={
#         "motor": "sample.x",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
#     axis_outer={
#         "motor": "sample.y",
#         "min": 0,
#         "max": 0,
#         "pts": 1,
#     },
# )


# -----------------------------------------------------------------------------
# Example 5: Rigaku ZDT
# -----------------------------------------------------------------------------
#
# pv_registers.det_name.put("rigaku3M")
# pv_registers.det_mode.put("ZDT")
# pv_registers.acq_time.put(2e-5)
# pv_registers.acq_period.put(2e-5)
# pv_registers.num_frames.put(100000)
# pv_registers.num_repeats.put(5)
#
# det_acq_series(
#     att_level=7,
#     wait_time=0,
#     move_index_register=0,
#     axis_inner={
#         "motor": "sample.x",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
#     axis_outer={
#         "motor": "sample.y",
#         "min": -1,
#         "max": 1,
#         "pts": 21,
#     },
# )