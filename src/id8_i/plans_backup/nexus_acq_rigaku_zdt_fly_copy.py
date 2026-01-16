"""
Simple, modular Bluesky plans for users.
"""

import os

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from bluesky import plans as bp

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
sample = oregistry["sample"]


def setup_rigaku_ZDT_fly(acq_time, num_frames, file_name):
    """Setup the Rigaku3M for ZDT fly scan acquisition.

    Configure the detector's cam module for internal acquisition mode and
    set up the file paths for data storage.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frames: Number of frames to acquire
        file_name: Base name for the output files
    """
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()

    file_path = f"{exp_name}/data/{file_name}"
    acq_period = acq_time

    yield from bps.mv(rigaku3M.cam.acquire_time, acq_time)
    yield from bps.mv(rigaku3M.cam.acquire_period, acq_period)
    yield from bps.mv(rigaku3M.cam.file_name, f"{file_name}.bin")
    yield from bps.mv(rigaku3M.cam.file_path, file_path)
    yield from bps.mv(rigaku3M.cam.num_images, num_frames)

    yield from bps.mv(pv_registers.file_name, file_name)
    yield from bps.mv(pv_registers.file_path, f"/gdata/dm/8IDI/{cycle_name}/{file_path}")
    yield from bps.mv(
        pv_registers.metadata_full_path,
        f"/gdata/dm/8IDI/{cycle_name}/{file_path}/{file_name}_metadata.hdf",
    )

    os.makedirs(f"/gdata/dm/8IDI/{cycle_name}/{file_path}", mode=0o770, exist_ok=True)


def rigaku_acq_ZDT_fly(
    acq_time=2e-5,
    num_frame=100000,
    num_rep=3,
    flyspeed=0.5,
    wait_time=0.0,
    sample_move=False,
    process=True,
):
    """Run ZDT fly scan acquisition with the Rigaku detector.

    This plan performs a continuous motion scan while collecting data,
    moving the sample stage at a constant velocity during acquisition.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frame: Number of frames to acquire
        num_rep: Number of repetitions
        flyspeed: Sample stage velocity during scan (mm/s)
        wait_time: Time to wait between repetitions
        sample_move: Whether to move sample between repetitions
        process: Whether to process data after acquisition
    """
    try:
        yield from post_align()
        yield from shutteroff()
        workflowProcApi, dmuser = dm_setup(process)
        folder_prefix = gen_folder_prefix()

        for ii in range(num_rep):
            yield from bps.sleep(wait_time)

            if sample_move:
                yield from mesh_grid_move()

            file_name = f"{folder_prefix}_f{num_frame:06d}_s{flyspeed*1000:04d}_r{ii+1:05d}"
            yield from setup_rigaku_ZDT_fly(acq_time, num_frame, file_name)

            extra_acq_time = 1.0
            total_travel = flyspeed * (acq_time * num_frame + extra_acq_time)
            yield from bps.mv(sample.y.velocity, flyspeed)
            yield from bps.mvr(sample.y.user_setpoint, total_travel)
            yield from showbeam()
            yield from bp.count([rigaku3M])
            yield from blockbeam()
            yield from bps.mv(sample.y.velocity, 5)
            yield from bps.mvr(sample.y, -total_travel)

            metadata_fname = pv_registers.metadata_full_path.get()
            create_nexus_format_metadata(metadata_fname, det=rigaku3M)

            dm_run_job("rigaku", process, workflowProcApi, dmuser, file_name)

    except Exception as e:
        print(f"Error occurred during measurement: {e}")
    finally:
        pass
