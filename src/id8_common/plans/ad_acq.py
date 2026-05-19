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

from ..utils.dm_util import dm_run_job
from ..utils.dm_util import dm_setup
from ..utils.nexus_utils import create_nexus_format_metadata

from .shutter_logic import blockbeam
from .shutter_logic import post_align
from .shutter_logic import showbeam
from .shutter_logic import shutteron
from .shutter_logic import shutteroff


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
        0  -> no sample motion and no register update
        8  -> use pv_registers.sample8_pos
        10 -> use pv_registers.sample10_pos

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

    pv_registers.file_name.put(file_name)
    pv_registers.metadata_full_path.put(f"{file_path}/{file_name}_metadata.hdf")

    eiger4M.cam.trigger_mode.put("Internal Series")
    eiger4M.cam.num_images.put(num_frames)
    eiger4M.cam.num_triggers.put(1)


def setup_eiger_external(acq_time, acq_period, num_frames, file_header, file_name):
    eiger4M = get_connected_device("eiger4M")
    softglue_8idi = get_connected_device("softglue_8idi")

    file_path = get_common_file_path(file_header, file_name)

    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_period)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)

    pv_registers.file_name.put(file_name)
    pv_registers.metadata_full_path.put(f"{file_path}/{file_name}_metadata.hdf")

    eiger4M.cam.num_triggers.put(num_frames)
    eiger4M.cam.trigger_mode.put("External Enable")

    softglue_8idi.acq_time.put(acq_time)
    softglue_8idi.acq_period.put(acq_period)
    softglue_8idi.num_triggers.put(num_frames)


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

    pv_registers.file_name.put(file_name)
    pv_registers.metadata_full_path.put(f"{file_path}/{file_name}_metadata.hdf")


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

    pv_registers.file_name.put(file_name)
    pv_registers.metadata_full_path.put(f"{full_path}/{file_name}_metadata.hdf")

    os.makedirs(full_path, mode=0o770, exist_ok=True)


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
    detector="eiger4M",
    mode="Internal Series",
    acq_time=1,
    acq_period=None,
    num_frames=10,
    num_reps=1,
    wait_time=0,
    header="A",
    sample_name="sample",
    move_index_register=0,
    axis_inner=None,
    axis_outer=None,
):
    """
    Run repeated detector acquisitions.

    detector:
        "eiger4M"
        "lambda2M"
        "rigaku3M"

    mode:
        eiger4M  : "Internal Series", "External Enable"
        lambda2M : "Internal"
        rigaku3M : "ZDT"

    acq_period:
        Required only for detector="eiger4M", mode="External Enable".

    move_index_register:
        0 means no sample motion.
        Nonzero integer N means use pv_registers.sampleN_pos.

    measurement_num:
        Read from pv_registers.measurement_num.
        If header="A" and measurement_num=12, file_header starts with A0012.
        The register increments once per det_acq_series() call.
    """
    mode_info = ACQ_MODES[detector][mode]

    if mode_info["needs_acq_period"] and acq_period is None:
        raise ValueError("This mode requires acq_period.")

    for device_name in mode_info["required_devices"]:
        get_connected_device(device_name)

    det = get_connected_device(detector)
    setup_func = mode_info["setup"]
    acquire_func = mode_info["acquire"]

    post_align()
    shutteroff()

    workflowProcApi, dmuser = dm_setup()

    meas_num = int(pv_registers.measurement_num.get())
    header_name = f"{header}{meas_num:04d}"

    pv_registers.measurement_num.put(meas_num + 1)
    pv_registers.sample_name.put(sample_name)

    file_header = f"{header_name}_{sample_name}_f{num_frames:06d}"

    for rep in range(num_reps):
        ttime.sleep(wait_time)

        sample_mesh_move(
            axis_inner=axis_inner,
            axis_outer=axis_outer,
            move_index_register=move_index_register,
        )

        file_name = f"{file_header}_r{rep + 1:05d}"

        if mode_info["needs_acq_period"]:
            setup_func(
                acq_time=acq_time,
                acq_period=acq_period,
                num_frames=num_frames,
                file_header=file_header,
                file_name=file_name,
            )
        else:
            setup_func(
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

        metadata_fname = pv_registers.metadata_full_path.get()
        create_nexus_format_metadata(metadata_fname, det=det)

        dm_run_job(workflowProcApi, dmuser)


# =============================================================================
# Example use cases
# =============================================================================

# -----------------------------------------------------------------------------
# Example 1: Eiger internal, no sample motion
# -----------------------------------------------------------------------------
#
# det_acq_series(
#     detector="eiger4M",
#     mode="Internal Series",
#     acq_time=0.1,
#     num_frames=1000,
#     num_reps=5,
#     wait_time=0,
#     header="A",
#     sample_name="G10",
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
# det_acq_series(
#     detector="eiger4M",
#     mode="Internal Series",
#     acq_time=0.1,
#     num_frames=1000,
#     num_reps=20,
#     wait_time=0,
#     header="A",
#     sample_name="G10",
#     move_index_register=8,
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
# det_acq_series(
#     detector="eiger4M",
#     mode="External Enable",
#     acq_time=0.1,
#     acq_period=1.0,
#     num_frames=1000,
#     num_reps=10,
#     wait_time=0,
#     header="A",
#     sample_name="G10",
#     move_index_register=8,
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
# det_acq_series(
#     detector="lambda2M",
#     mode="Internal",
#     acq_time=0.1,
#     num_frames=1000,
#     num_reps=20,
#     wait_time=0,
#     header="A",
#     sample_name="G10",
#     move_index_register=8,
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
# det_acq_series(
#     detector="rigaku3M",
#     mode="ZDT",
#     acq_time=2e-5,
#     num_frames=100000,
#     num_reps=5,
#     wait_time=0,
#     header="A",
#     sample_name="G10",
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