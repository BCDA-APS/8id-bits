"""
Simple, modular Ophyd scripts for users.
"""

import os
import time
import subprocess
import shlex

from apsbits.core.instrument_init import oregistry

from ..utils.dm_util import dm_run_job
from ..utils.dm_util import dm_setup
from ..utils.nexus_utils import create_nexus_format_metadata
from .sample_info_unpack import gen_folder_prefix
from .sample_info_unpack import mesh_grid_move
from .shutter_logic import blockbeam
from .shutter_logic import post_align
from .shutter_logic import showbeam
from .shutter_logic import shutteroff

rigaku3M = oregistry["rigaku3M"]
pv_registers = oregistry["pv_registers"]


def setup_rigaku_ZDT_series(acq_time, num_frames, file_header, file_name):
    """Setup the Rigaku3M for ZDT series acquisition.

    Configure the detector's cam module for internal acquisition mode and
    set up the file paths for data storage.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frames: Number of frames to acquire
        file_name: Base name for the output files
    """
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    use_subfolder = pv_registers.use_subfolder.get()
    mount_point = pv_registers.mount_point.get()

    if use_subfolder == 'Yes':
        file_path = f"{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == 'No':   
        file_path = f"{exp_name}/data/{file_name}"
    else: 
        print("Sub folder options can only be either Yes or No")

    acq_period = acq_time

    rigaku3M.cam.acquire_time.put(acq_time)
    rigaku3M.cam.acquire_period.put(acq_period)
    rigaku3M.cam.fast_file_name.put(f"{file_name}.bin")
    rigaku3M.cam.fast_file_path.put(file_path)
    rigaku3M.cam.num_images.put(num_frames)
    rigaku3M.cam.output_control.put('Sparsified')
    rigaku3M.cam.output_resolution.put('2 Bit')

    pv_registers.file_name.put(file_name)
    pv_registers.metadata_full_path.put(
        f"{mount_point}/{cycle_name}/{file_path}/{file_name}_metadata.hdf"
    )

    os.makedirs(f"{mount_point}/{cycle_name}/{file_path}", mode=0o770, exist_ok=True)

    # full_path = f"/gdata/dm/8ID/8IDI/{cycle_name}/{file_path}"
    # remote_cmd = f"mkdir -m 770 -p {shlex.quote(full_path)}"
    # subprocess.run(
    #     ["ssh", "s8ididm", remote_cmd],
    #     check=True,
    # )


############# Homebrew acquisition plan #############
def rigaku_zdt_acquire():
    """Run the Rigaku ZDT acquisition sequence."""
    showbeam()
    time.sleep(0.1)
    rigaku3M.cam.acquire.put(1)
    # time.sleep(2.0)

    # Wait for detector to start (status == 1)
    while True:
        det_status = rigaku3M.cam.detector_state.get()
        if det_status != 1:
            time.sleep(0.1)
        if det_status == 1:
            break

    # Wait for detector to finish (status == 0)
    while True:
        det_status = rigaku3M.cam.detector_state.get()
        if det_status != 0:
            time.sleep(0.1)
        if det_status == 0:
            break

    blockbeam()


############# Homebrew acquisition plan ends #############


def rigaku_acq_ZDT_series(
    acq_time=2e-5,
    num_frames=100000,
    num_reps=2,
    wait_time=0,
    process=True,
    sample_move=False,
):
    """Run ZDT series acquisition with the Rigaku detector.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frame: Number of frames to acquire
        num_reps: Number of repetitions
        wait_time: Time to wait between repetitions
        process: Whether to process data after acquisition
        sample_move: Whether to move sample between repetitions
    """
    # try:
    post_align()
    shutteroff()
    workflowProcApi, dmuser = dm_setup()
    folder_prefix = gen_folder_prefix()

    for ii in range(num_reps):
        if sample_move:
            mesh_grid_move()

        file_header = f"{folder_prefix}_f{num_frames:06d}"
        file_name = f"{folder_prefix}_f{num_frames:06d}_r{ii+1:05d}"
        print(file_name)
        
        setup_rigaku_ZDT_series(acq_time, num_frames, file_header, file_name)

        print(f"\nStarting Measurement {file_name}")
        rigaku_zdt_acquire()
        print(f"Measurement {file_name} Complete")

        metadata_fname = pv_registers.metadata_full_path.get()
        create_nexus_format_metadata(metadata_fname, det=rigaku3M)

        dm_run_job(workflowProcApi, dmuser)

        time.sleep(wait_time)

    # except Exception as e:
    #     print(f"Error occurred during measurement: {e}")
    # finally:
    #     pass