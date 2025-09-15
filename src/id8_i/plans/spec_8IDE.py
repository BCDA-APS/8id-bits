"""
SPEC integration plans for the 8ID-E Eiger4M detector.

This module provides plans for integrating with SPEC software, particularly
for submitting data processing jobs to the Data Management (DM) system when
triggered by SPEC through EPICS PVs.
"""

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from ..utils.dm_util import dm_run_job
from ..utils.dm_util import dm_setup
from ..utils.nexus_utils import create_nexus_format_metadata

# Add more detectors if needed
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]
# rigaku3M = oregistry["rigaku3M"]
pv_registers = oregistry["pv_registers"]


def submit_Nexus_DM(det=None):
    """Submit data processing jobs to DM when triggered by SPEC.

    This plan monitors a trigger PV from SPEC and, when activated, submits
    a data processing job to the Data Management system. It creates the
    necessary metadata files and configures the workflow based on the
    experiment parameters.

    Yields:
        Generator: Bluesky plan messages
    """
    while True:
        bluesky_start = pv_registers.start_bluesky.get()
        if bluesky_start == "Yes":

            metadata_fname = pv_registers.metadata_full_path.get()
            create_nexus_format_metadata(metadata_fname, det)

            workflowProcApi, dmuser = dm_setup()
            dm_run_job("eiger", workflowProcApi, dmuser)
            pv_registers.start_bluesky.put("No")

        else:
            yield from bps.sleep(0.1)
