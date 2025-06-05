"""
External trigger acquisition plans for the Eiger4M detector.

This module provides plans for controlling the Eiger4M detector in external
trigger mode, including setup of the detector, softglue triggers, and data
acquisition workflows.
"""

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
from .shutter_logic import shutteron

eiger4M = oregistry["eiger4M"]
softglue_8idi = oregistry["softglue_8idi"]
pv_registers = oregistry["pv_registers"]


def setup_eiger_ext_trig(
    acq_time: float,
    acq_period: float,
    num_frames: int,
    file_name: str,
):
    """Setup the Eiger4M for external trigger mode.

    Args:
        acq_time: Acquisition time per frame in seconds
        acq_period: Time between frames in seconds
        num_frames: Number of frames to acquire
        file_name: Base name for the output files
    """
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()

    file_path = f"/gdata/dm/8IDI/{cycle_name}/{exp_name}/data/{file_name}/"

    yield from bps.mv(eiger4M.cam.trigger_mode, "External Enable")
    yield from bps.mv(eiger4M.cam.acquire_time, acq_time)
    yield from bps.mv(eiger4M.cam.acquire_period, acq_period)
    yield from bps.mv(eiger4M.hdf1.file_name, file_name)
    yield from bps.mv(eiger4M.hdf1.file_path, file_path)
    yield from bps.mv(eiger4M.cam.num_triggers, num_frames)
    yield from bps.mv(eiger4M.hdf1.num_capture, num_frames)

    yield from bps.mv(pv_registers.file_name, file_name)
    yield from bps.mv(pv_registers.file_path, file_path)
    yield from bps.mv(
        pv_registers.metadata_full_path,
        f"{file_path}/{file_name}_metadata.hdf",
    )


def setup_softglue_ext_trig(
    acq_time: float,
    acq_period: float,
    num_frames: int,
):
    """Setup external triggering parameters.

    Args:
        acq_time: Acquisition time per frame in seconds
        acq_period: Time between frames in seconds
        num_frames: Number of frames to acquire
    """
    yield from bps.mv(softglue_8idi.acq_time, acq_time)
    yield from bps.mv(softglue_8idi.acq_period, acq_period)
    # Generate n+1 triggers, in case softglue triggered before area detector.
    # Add 10 more pulses to be safe since we're also sending stop signal.
    yield from bps.mv(softglue_8idi.num_triggers, num_frames + 10)


def softglue_start_pulses():
    """Start generating trigger pulses."""
    yield from bps.mv(softglue_8idi.start_pulses, "1!")


def softglue_stop_pulses():
    """Stop generating trigger pulses."""
    yield from bps.mv(softglue_8idi.stop_pulses, "1!")


def eiger_acq_ext_trig(
    acq_time: float = 1,
    acq_period: float = 2,
    num_frames: int = 10,
    num_rep: int = 2,
    wait_time: float = 0,
    sample_move: bool = True,
    process: bool = True,
):
    """Run an external trigger acquisition sequence.

    Args:
        acq_time: Acquisition time per frame in seconds
        acq_period: Time between frames in seconds
        num_frames: Number of frames to acquire
        num_rep: Number of repetitions
        wait_time: Time to wait between repetitions
        sample_move: Whether to move sample between repetitions
        process: Whether to process data after acquisition
    """
    try:
        yield from setup_softglue_ext_trig(acq_time, acq_period, num_frames)
        yield from post_align()
        yield from shutteron()

        workflowProcApi, dmuser = dm_setup(process)
        folder_prefix = gen_folder_prefix()

        for ii in range(num_rep):
            yield from bps.sleep(wait_time)

            if sample_move:
                yield from mesh_grid_move()

            file_name = f"{folder_prefix}_f{num_frames:06d}_r{ii+1:05d}"
            yield from setup_eiger_ext_trig(acq_time, acq_period, num_frames, file_name)

            yield from showbeam()
            yield from bps.sleep(0.1)
            yield from softglue_start_pulses()
            yield from bp.count([eiger4M])
            yield from softglue_stop_pulses()
            yield from blockbeam()

            metadata_fname = pv_registers.metadata_full_path.get()
            create_nexus_format_metadata(metadata_fname, det=eiger4M)

            dm_run_job("eiger", process, workflowProcApi, dmuser, file_name)
    except Exception as e:
        print(f"Error occurred during measurement: {e}")
    finally:
        pass
