"""
Simple, modular Bluesky plans for users.
"""

import os
import pathlib
import subprocess
from datetime import datetime

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

from ..utils.nexus_utils import create_nexus_format_metadata
from .sample_info_unpack import gen_folder_prefix
from .sample_info_unpack import mesh_grid_move
from .shutter_logic import blockbeam
from .shutter_logic import post_align
from .shutter_logic import showbeam
from .shutter_logic import shutteroff

# --MC--
# from timepix_aps.xpcs_proc import timepix_converter

pv_registers = oregistry["pv_registers"]

rigaku3M = oregistry["rigaku3M"]


def setup_tempus_int_series(file_name):
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
    mount_point = "/home/8-id-i/"

    file_path = f"{mount_point}{cycle_name}/{exp_name}/data_tempus/{file_name}"

    yield from bps.mv(pv_registers.file_name, file_name)
    yield from bps.mv(pv_registers.metadata_full_path, f"{file_path}/{file_name}_metadata.hdf")


############# Homebrew acquisition plan #############
def tempus_acquire(file_name):
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    mount_point = "/home/8-id-i/"
    # mount_point = pv_registers.mount_point.get()

    file_path = f"{mount_point}{cycle_name}/{exp_name}/data_tempus/{file_name}"
    file_full_path = os.path.join(file_path, file_name)

    yield from showbeam()
    yield from bps.sleep(0.1)
    subprocess.run(
        [
            "/local/tempus/tempus-receiver-split/target/debug/tempus-receiver",
            "-b",
            file_full_path,
            "-t",
            "2",
            "-r",
        ]
    )
    yield from blockbeam()

    # --MC--
    # convert the timepix format to rigaku format
    target_folder = str(file_path)
    output_dir = "/gdata/dm/8ID/8IDI/2025-2/tempus202507d/data/timepix/"  # This line is ignored
    time_total = "1"  # total time in seconds for all frames
    time_bin = "383e-9"  # delta_t in seconds; time bin size\                   # This line is ignored
    assert (
        int(float(time_total) / float(time_bin)) <= 2**24
    ), "It will overflow with this configuration. increase time-bin or decrease acquisition time"
    # subprocess.run(["bash", "/home/beams4/MQICHU/bin/launch_timepix_converter_remote",
    #                 target_folder, output_dir, time_total, time_bin])
    # subprocess.run(["bash", "/home/beams4/MQICHU/bin/launch_timepix_converter_remote_background",
    #                 target_folder, output_dir, time_total, time_bin])

    subprocess.run(
        [
            "bash",
            "/home/beams/8IDIUSER/bin/launch_timepix_converter_remote_fixed",
            target_folder,
            output_dir,
            time_total,
            time_bin,
        ]
    )


############# Homebrew acquisition plan ends #############


def tempus_acq_int_series(
    num_frames=2000000,
    num_rep=3,
    wait_time=0,
    sample_move=False,
):
    # try:home/beams4
    yield from post_align()
    yield from shutteroff()
    # workflowProcApi, dmuser = dm_setup()
    folder_prefix = gen_folder_prefix()

    for ii in range(num_rep):
        yield from bps.sleep(wait_time)

        if sample_move:
            yield from mesh_grid_move()

        file_name = f"{folder_prefix}_f{num_frames:06d}_r{ii+1:05d}"
        yield from setup_tempus_int_series(file_name)

        metadata_fname = pv_registers.metadata_full_path.get()
        metadata_fname_path = pathlib.Path(metadata_fname)
        metadata_fname_path.parent.mkdir(parents=True, exist_ok=True)
        create_nexus_format_metadata(metadata_fname, det=rigaku3M)

        _ = datetime.now()
        time_now = _.strftime("%Y-%m-%d %H:%M:%S")
        print(f"\n{time_now}, Starting measurement {file_name}")
        yield from tempus_acquire(file_name)
        _ = datetime.now()
        time_now = _.strftime("%Y-%m-%d %H:%M:%S")
        print(f"{time_now}, Complete measurement {file_name}")

        # dm_run_job("tempus", workflowProcApi, dmuser)
    # except KeyboardInterrupt:
    #     raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).")
    # except Exception as e:
    #     print(f"Error occurred during measurement: {e}")
    # finally:
    #     pass
