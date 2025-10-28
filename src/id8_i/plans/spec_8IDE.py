"""
SPEC integration plans for the 8ID-E Eiger4M detector.

This module provides plans for integrating with SPEC software, particularly
for submitting data processing jobs to the Data Management (DM) system when
triggered by SPEC through EPICS PVs.
"""

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

from ..utils.nexus_utils import create_nexus_format_metadata

# Add more detectors if needed
eiger4M = oregistry["eiger4M"]
# lambda2M = oregistry["lambda2M"]
# lambda750k = oregistry["lambda750k"]
# rigaku3M = oregistry["rigaku3M"]
pv_registers = oregistry["pv_registers"]


def submit_Nexus_DM():
    """Submit data processing jobs to DM when triggered by SPEC.

    This plan monitors a trigger PV from SPEC and, when activated, submits
    a data processing job to the Data Management system. It creates the
    necessary metadata files and configures the workflow based on the
    experiment parameters.

    Yields:
        Generator: Bluesky plan messages
    """

    while True:
        metadata_fname = pv_registers.metadata_full_path.get()

        # Convert detector name to detector Ophyd object
        det_name = pv_registers.det_name.get()
        if det_name == "eiger4M":
            det = eiger4M
        # elif det_name == "lambda2M":
        #     det = lambda2M
        # elif det_name == "lambda750k":
        #     det = lambda750k
        elif det_name == "rigaku3M":
            det = rigaku3M
        else:
            det = None
            print("Detector name not found")
        # Convert detector name to detector Ophyd object

        bluesky_start = pv_registers.start_bluesky.get()
        if bluesky_start == "Yes":
            # Create metadata
            create_nexus_format_metadata(metadata_fname, det)

            # Submit analysis via DM workflow
            # workflowProcApi, dmuser = dm_setup()
            # dm_run_job(workflowProcApi, dmuser)

            # Flip Bluesky PV back to No
            print(f"Metadata written at {metadata_fname} \n")
            pv_registers.start_bluesky.put("No")

        else:
            yield from bps.sleep(0.5)

        # Delay time between every execution in the While loop
        yield from bps.sleep(0.5)
