"""
Simple, modular Bluesky plans for users.
"""

from datetime import datetime
import time

from apsbits.core.instrument_init import oregistry

from ..utils.dm_util import dm_run_job
from ..utils.dm_util import dm_setup
from ..utils.nexus_utils import create_nexus_format_metadata
from .sample_info_unpack import gen_folder_prefix
from .sample_info_unpack import mesh_grid_move
from .shutter_logic import *
from .qnw_plans import find_qnw_index

eiger4M = oregistry["eiger4M"]
pv_registers = oregistry["pv_registers"]
filter = oregistry["filter_8ide"]

qnw_env1 = oregistry["qnw_env1"]
qnw_env2 = oregistry["qnw_env2"]
qnw_env3 = oregistry["qnw_env3"]
qnw_controllers = [qnw_env1, qnw_env2, qnw_env3]

def setup_eiger_int_series(acq_time, num_frames, file_header, file_name):
    """Setup the Eiger4M for internal series acquisition.

    Configure the detector's cam module for internal acquisition mode and
    set up the HDF plugin for data storage.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frames: Number of frames to acquire
        file_name: Base name for the output files
    """
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    mount_point = pv_registers.mount_point.get()
    use_subfolder = pv_registers.use_subfolder.get()

    if use_subfolder == "Yes":
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_header}/{file_name}" 
    elif use_subfolder == "No":   
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_name}"
    else: 
        print("Sub folder options can only be either Yes or No") 

    eiger4M.cam.acquire_time.put(acq_time)
    # moved up by damian 2025-11-15, acq_period was defined after calling it:
    acq_period = acq_time
    eiger4M.cam.acquire_period.put(acq_period)
    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)
    pv_registers.file_name.put(file_name)
    pv_registers.metadata_full_path.put(f"{file_path}/{file_name}_metadata.hdf")

    # Unique for Internal Mode
    # moved up by damian 2025-11-15, to set internal series first:
    eiger4M.cam.trigger_mode.put("Internal Series")  # 0
    eiger4M.cam.num_images.put(num_frames)  # Needs to be set after the trigger mode
    eiger4M.cam.num_triggers.put(1)  # Need to put num_trigger to 1 for internal mode

    
############# Homebrew acquisition plan #############
def eiger_acquire():
    """Homebrew plan to acquire data with Eiger detector in internal mode."""
    shutteroff()
    showbeam()
    time.sleep(0.1)
    eiger4M.hdf1.capture.put(1)
    eiger4M.cam.acquire.put(1)

    while True:
        #### QZ on 2026/01/06 ####
        # without the 0.5 s wait time, the repeating acqs go out of sync. 
        # Don't know why and maybe the 0.5 s can be made shorter
        #### QZ on 2026/01/06 ####
        time.sleep(0.5)
        det_status = eiger4M.cam.acquire.get()
        if det_status == 1:
            time.sleep(0.1)
        if det_status == 0:
            break
    blockbeam()

    frame_num_set = eiger4M.hdf1.queue_size.get()
    count = 0
    while count < 100:
        frame_num_processed = eiger4M.hdf1.queue_free.get()
        if frame_num_processed == frame_num_set:
            break
        else:
           time.sleep(0.1)
        count = +1
        eiger4M.hdf1.capture.put(0)


############# Homebrew acquisition plan ends #############


def eiger_acq_int_series_wei(
    acq_time=1,
    num_frames=10,
    num_reps=3,
    wait_time=0,
    sample_name=None,
    sample_move=False,
):
    """Run internal series acquisition with the Eiger detector.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frames: Number of frames to acquire
        num_reps: Number of repetitions
        wait_time: Time to wait between repetitions
        process: Whether to process data after acquisition
        sample_move: Whether to move sample between repetitions
    """
    try:
        post_align()
        shutteroff()
        workflowProcApi, dmuser = dm_setup()

        ############# QZ 2026/03/26: Adding option to change sample name from inside the acq code #############
        meas_num = int(pv_registers.measurement_num.get())
        pv_registers.measurement_num.put(meas_num+1)

        att_level = int(filter.attenuation.readback.get())

        header_name = f'R{meas_num:04d}'
        folder_prefix = f"{header_name}_{sample_name}_a{att_level:04}"
        ############# QZ 2026/03/26: Adding option to change sample name from inside the acq code #############

        for ii in range(num_reps):
            time.sleep(wait_time)

            if sample_move:
                mesh_grid_move()

            file_header = f"{folder_prefix}_f{num_frames:06d}"
            file_name = f"{folder_prefix}_f{num_frames:06d}_r{ii+1:05d}"

            setup_eiger_int_series(acq_time, num_frames, file_header, file_name)

            _ = datetime.now()
            time_now = _.strftime("%Y-%m-%d %H:%M:%S")
            print(f"\n{time_now}, Starting measurement {file_name}")
            eiger_acquire()
            _ = datetime.now()
            time_now = _.strftime("%Y-%m-%d %H:%M:%S")
            print(f"{time_now}, Complete measurement {file_name}")

            metadata_fname = pv_registers.metadata_full_path.get()
            create_nexus_format_metadata(metadata_fname, det=eiger4M)

            dm_run_job(workflowProcApi, dmuser)
    except KeyboardInterrupt:
        raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).")
    except Exception as e:
        print(f"Error occurred during measurement: {e}")
    finally:
        pass
