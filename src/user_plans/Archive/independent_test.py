"""
Data acquisition steps as independent bluesky plans.

.. autosummary::
    ~kickoff_dm_workflow
    ~setup_detector
    ~simple_acquire

This plan is an example that combines the above plans.

.. autosummary::
    ~example_full_acquisition
"""

import pathlib
import epics as pe

from apstools.devices import DM_WorkflowConnector
from apstools.utils import share_bluesky_metadata_with_dm
from bluesky import plan_stubs as bps
from bluesky import plans as bp
from bluesky import preprocessors as bpp

from aps_8id_bs_instrument.callbacks.nexus_data_file_writer import nxwriter
from aps_8id_bs_instrument.callbacks.devices.ad_eiger_4M import nxwriter
from aps_8id_bs_instrument.callbacks.initialize_bs_tools import cat


EMPTY_DICT = {}  # Defined as symbol to pass the style checks.
DATA_FOLDER = '/gdata/dm/8IDI/2024-3/comm202410/data/'

def simple_acquire(det, file_name, md: dict = EMPTY_DICT):
    """Just run the acquisition and save the file, nothing else."""

    nxwriter.warn_on_missing_content = False
    nxwriter.file_path = DATA_FOLDER+file_name+'/'
    nxwriter.file_name = nxwriter.file_path / (file_name + ".hdf")

    md["metadatafile"] = nxwriter.file_name.name

    @bpp.subs_decorator(nxwriter.receiver)
    def acquire():
        yield from bp.count([det], md=md)

    # Start the acquire. Eiger will wait for external trigger pulse sequence
    yield from acquire()

    # Send the external trigger pulse sequence from softglue
    pe.put("8idi:softGlueA:MUX2-1_IN0_Signal","1!\"")

    # Wait for NeXus metadata file content to flush to disk.
    # If next acquisition proceeds without waiting, the
    # metadata file will be spoiled.
    yield from nxwriter.wait_writer_plan_stub()


def setup_det_ext_trig(det, acq_time, acq_period, num_frames, file_name):
    
    """Setup the cam module for external trigger (3) mode and populate the hdf plugin"""
    yield from bps.mv(det.cam.trigger_mode, 3)

    yield from bps.mv(det.cam.acquire_time, acq_time)
    yield from bps.mv(det.cam.acquire_period, acq_period)
    yield from bps.mv(det.hdf1.file_name, file_name)
    yield from bps.mv(det.hdf1.file_path, DATA_FOLDER+file_name+'/')
    yield from bps.mv(det.cam.num_triggers, num_frames)


def setup_softglue_ext_trig(acq_time, acq_period):

    pe.caput("8idi:SGControl1.C", acq_time)
    pe.caput("8idi:SGControl1.A", acq_period)


def kickoff_dm_workflow(
    experiment_name,
    file_name,
    qmap_file,
    run,
    analysisMachine="amazonite",
):
    """Start a DM workflow for this bluesky run."""
    dm_workflow = DM_WorkflowConnector(name="dm_workflow")
    forever = 999_999_999_999  # long time, s, disables periodic reports
    workflow_name = "xpcs8-02-gladier-boost"

    yield from bps.mv(dm_workflow.concise_reporting, True)
    yield from bps.mv(dm_workflow.reporting_period, forever)

    # DM argsDict content
    argsDict = dict(
        filePath=file_name,
        experimentName=experiment_name,
        qmap=qmap_file,
        smooth="sqmap",
        gpuID=-1,
        beginFrame=1,
        endFrame=-1,
        strideFrame=1,
        avgFrame=1,
        type="Multitau",
        dq="all",
        verbose=False,
        saveG2=False,
        overwrite=False,
        analysisMachine=analysisMachine,
    )

    yield from dm_workflow.run_as_plan(
        workflow=workflow_name,
        wait=False,
        timeout=forever,
        **argsDict,
    )

    # Upload bluesky run metadata to APS DM.
    share_bluesky_metadata_with_dm(experiment_name, workflow_name, run)

    # Users requested the DM workflow job ID be printed to the console.
    dm_workflow._update_processing_data()
    job_id = dm_workflow.job_id.get()
    job_stage = dm_workflow.stage_id.get()
    job_status = dm_workflow.status.get()
    print(f"DM workflow id: {job_id!r}  status: {job_status}  stage: {job_stage}")


def example_full_acquisition():
    """These are the data acquisition steps for a user."""
    det = eiger4M
    yield from setup_det_ext_trig(det, 0.1, 1000, "A001_001")

    uids = yield from simple_acquire(det)
    print(f"Bluesky run: {uids=}")
    run = cat[uids[0]]

    try:
        yield from nxwriter.wait_writer_plan_stub()
        image_file_name = pathlib.Path(det.hdf1.full_file_name.get()).name
        print(f"{image_file_name=!r}")

        yield from kickoff_dm_workflow(
            "my_dm_experiment",
            image_file_name,
            "my_qmap_file.h5",
            run,
            analysisMachine="amazonite",
        )
    except Exception as exc:
        print(f"Exception: {exc}")