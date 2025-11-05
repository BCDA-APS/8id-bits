"""
Simple, modular Bluesky plans for users.
"""

import os

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

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


def setup_rigaku_ZDT_series(acq_time, num_frames, file_name):
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

    file_path = f"{exp_name}/data/{file_name}"
    acq_period = acq_time

    yield from bps.mv(rigaku3M.cam.acquire_time, acq_time)
    yield from bps.mv(rigaku3M.cam.acquire_period, acq_period)
    yield from bps.mv(rigaku3M.cam.fast_file_name, f"{file_name}.bin")
    yield from bps.mv(rigaku3M.cam.fast_file_path, file_path)
    yield from bps.mv(rigaku3M.cam.num_images, num_frames)

    yield from bps.mv(pv_registers.file_name, file_name)
    yield from bps.mv(pv_registers.metadata_full_path, f"/gdata/dm/8ID/8IDI/{cycle_name}/{file_path}")
    yield from bps.mv(
        pv_registers.metadata_full_path,
        f"/gdata/dm/8ID/8IDI/{cycle_name}/{file_path}/{file_name}_metadata.hdf",
    )

    os.makedirs(f"/gdata/dm/8ID/8IDI/{cycle_name}/{file_path}", mode=0o770, exist_ok=True)


############# Homebrew acquisition plan #############
def rigaku_zdt_acquire():
    """Run the Rigaku ZDT acquisition sequence."""

    yield from showbeam()
    yield from bps.sleep(0.1)
    yield from bps.mv(rigaku3M.cam.acquire, 1)
    # yield from bps.sleep(2.0)

    while True:
        det_status = rigaku3M.cam.detector_state.get()
        if det_status != 1:
            yield from bps.sleep(0.1)
        if det_status == 1:
            break

    while True:
        det_status = rigaku3M.cam.detector_state.get()
        if det_status != 0:
            yield from bps.sleep(0.1)
        if det_status == 0:
            break

    yield from blockbeam()


############# Homebrew acquisition plan ends #############


def rigaku_acq_ZDT_series(
    acq_time=2e-5,
    num_frames=100000,
    num_rep=2,
    wait_time=0,
    process=True,
    sample_move=False,
):
    """Run ZDT series acquisition with the Rigaku detector.

    Args:
        acq_time: Acquisition time per frame in seconds
        num_frame: Number of frames to acquire
        num_rep: Number of repetitions
        wait_time: Time to wait between repetitions
        process: Whether to process data after acquisition
        sample_move: Whether to move sample between repetitions
    """
    # try:
    yield from post_align()
    yield from shutteroff()
    workflowProcApi, dmuser = dm_setup()
    folder_prefix = gen_folder_prefix()

    for ii in range(num_rep):
        if sample_move:
            yield from mesh_grid_move()

        file_name = f"{folder_prefix}_f{num_frames:06d}_r{ii+1:05d}"
        print(file_name)
        yield from setup_rigaku_ZDT_series(acq_time, num_frames, file_name)

        print(f"\nStarting Measurement {file_name}")
        yield from rigaku_zdt_acquire()
        print(f"Measurement {file_name} Complete")

        metadata_fname = pv_registers.metadata_full_path.get()
        create_nexus_format_metadata(metadata_fname, det=rigaku3M)

        dm_run_job(workflowProcApi, dmuser)

        yield from bps.sleep(wait_time)

    # except Exception as e:
    #     print(f"Error occurred during measurement: {e}")
    # finally:
    #     pass
