"""
Simple, modular Ophyd scripts for users.
"""

from datetime import datetime
import time

from apsbits.core.instrument_init import oregistry

from ..utils.dm_util import dm_run_job
from ..utils.dm_util import dm_setup
from ..utils.nexus_utils import create_nexus_format_metadata
from .sample_info_unpack import gen_folder_prefix
# from .sample_info_unpack import mesh_grid_move
from .shutter_logic import *

lambda2M = oregistry["lambda2M"]
pv_registers = oregistry["pv_registers"]
keysight = oregistry["keysight"]


def setup_lambda_int_series(acq_time, num_frames, file_header, file_name):
    """Setup the Lambda2M for internal series acquisition.

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

    file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_header}/{file_name}"
    
    acq_period = acq_time

    # Use .put() for signals
    lambda2M.hdf1.enable.put(1)
    # lambda2M.cam.save_files.put(0)
    # lambda2M.cam.fw_enable.put(1)
    lambda2M.cam.trigger_mode.put(0)  # 0
    lambda2M.cam.acquire_time.put(acq_time)
    lambda2M.cam.acquire_period.put(acq_period)
    lambda2M.hdf1.file_name.put(file_name)
    lambda2M.hdf1.file_path.put(file_path)
    lambda2M.cam.num_images.put(num_frames)
    # lambda2M.cam.num_triggers.put(1)  # Need to put num_trigger to 1 for internal mode
    lambda2M.hdf1.num_capture.put(num_frames)

    pv_registers.file_name.put(file_name)
    pv_registers.metadata_full_path.put(f"{file_path}/{file_name}_metadata.hdf")


############# Homebrew acquisition plan #############
def lambda_acquire():
    """Homebrew script to acquire data with Lambda2M detector in internal mode."""

    # First stop any TV mode that might be running. 
    # Adding 1 second of sleep time because it can take a while for the detector tp stop
    # Comment this out if running rapid repeats
    lambda2M.cam.acquire.put(0)
    time.sleep(1.0)
        
    showbeam()
    time.sleep(0)
    lambda2M.hdf1.capture.put(1)
    lambda2M.cam.acquire.put(1)

    while True:
        #### QZ on 2026/01/06 ####
        # without the 0.5 s wait time, the repeating acqs go out of sync. 
        # Don't know why and maybe the 0.5 s can be made shorter
        #### QZ on 2026/01/06 ####
        time.sleep(0.05)
        det_status = lambda2M.cam.acquire.get()
        if det_status == 1:
            time.sleep(0.05)
        if det_status == 0:
            break
    blockbeam()


    while True:
        #### QZ on 2026/02/18 ####
        # Suresh suggested to just check for hdf plugin done status
        #### QZ on 2026/02/18 ####
        time.sleep(0.05)
        det_plugin_status = lambda2M.hdf1.capture.get()
        if det_plugin_status == 1:
            time.sleep(0.05)
        if det_plugin_status == 0:
            break

    # frame_num_set = lambda2M.hdf1.queue_size.get()
    # count = 0
    # while count < 100:
    #     frame_num_processed = lambda2M.hdf1.queue_free.get()
    #     if frame_num_processed == frame_num_set:
    #         break
    #     else:
    #         time.sleep(0.1)
    #         count = +1
    #     lambda2M.hdf1.capture.put(0)


############# Homebrew acquisition plan ends #############


def lambda_acq_int_series(
    acq_time=1,
    num_frames=10,
    num_reps=3,
    wait_time=0,
    sample_move=True,
):
    """Run internal series acquisition with the Lambda2M detector.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frames: Number of frames to acquire
        num_rep: Number of repetitions
        wait_time: Time to wait between repetitions
        process: Whether to process data after acquisition
        sample_move: Whether to move sample between repetitions
    """
    try:
        shutteroff()
        workflowProcApi, dmuser = dm_setup()
        folder_prefix = gen_folder_prefix()

        for ii in range(num_reps):
            time.sleep(wait_time)

            # if sample_move:
            #     mesh_grid_move()

            file_header = f"{folder_prefix}_f{num_frames:06d}"
            file_name = f"{folder_prefix}_f{num_frames:06d}_r{ii+1:05d}"
            
            setup_lambda_int_series(acq_time, num_frames, file_header, file_name)

            _ = datetime.now()
            time_now = _.strftime("%Y-%m-%d %H:%M:%S")
            print(f"\n{time_now}, Starting measurement {file_name}")
            
            lambda_acquire()
            
            _ = datetime.now()
            time_now = _.strftime("%Y-%m-%d %H:%M:%S")
            print(f"{time_now}, Complete measurement {file_name}")

            metadata_fname = pv_registers.metadata_full_path.get()
            create_nexus_format_metadata(metadata_fname, det=lambda2M)

            dm_run_job(workflowProcApi, dmuser)
    except KeyboardInterrupt:
        shutteroff()
        blockbeam()
        lambda2M.cam.acquire.put(0)
        lambda2M.hdf1.capture.put(0)
        raise RuntimeError("\n Script stopped by user (Ctrl+C).")
    except Exception as e:
        print(f"Error occurred during measurement: {e}")
    finally:
        pass

def softtrig_lambda_acq_int_series(
    acq_time=1,
    num_frames=10,
    num_reps=3,
    wait_time=0,
    sample_move=False,
):
    """Run internal series acquisition with the Lambda2M detector.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frames: Number of frames to acquire
        num_rep: Number of repetitions
        wait_time: Time to wait between repetitions
        process: Whether to process data after acquisition
        sample_move: Whether to move sample between repetitions
    """
    try:
        PIND_status(0)
        shutteroff()
        workflowProcApi, dmuser = dm_setup()
        folder_prefix = gen_folder_prefix()

        for ii in range(num_reps):
            time.sleep(wait_time)

            # if sample_move:
            #     mesh_grid_move()

            file_header = f"{folder_prefix}_f{num_frames:06d}"
            file_name = f"{folder_prefix}_f{num_frames:06d}_r{ii+1:05d}"
            
            setup_lambda_int_series(acq_time, num_frames, file_header, file_name)

            _ = datetime.now()
            time_now = _.strftime("%Y-%m-%d %H:%M:%S")
            print(f"\n{time_now}, Starting measurement {file_name}")
            
            keysight.output.put(1)
            lambda_acquire()
            keysight.output.put(0)
            
            _ = datetime.now()
            time_now = _.strftime("%Y-%m-%d %H:%M:%S")
            print(f"{time_now}, Complete measurement {file_name}")

            metadata_fname = pv_registers.metadata_full_path.get()
            create_nexus_format_metadata(metadata_fname, det=lambda2M)

            dm_run_job(workflowProcApi, dmuser)
        

    except KeyboardInterrupt:
        shutteroff()
        blockbeam()
        lambda2M.cam.acquire.put(0)
        lambda2M.hdf1.capture.put(0)
        raise RuntimeError("\n Script stopped by user (Ctrl+C).")
    except Exception as e:
        print(f"Error occurred during measurement: {e}")
    finally:
        pass