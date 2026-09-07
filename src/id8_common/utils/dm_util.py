"""
DM code from Hannah Parraga.
Set up DM and submit jobs
"""

import datetime
from pathlib import Path

from dm.common.utility.configurationManager import ConfigurationManager
from dm.proc_web_service.api.workflowProcApi import WorkflowProcApi

from .misc import get_machine_name
from id8_common.expt_config import expt

def dm_setup() -> tuple:
    """Set up the Data Management workflow API.

    Credentials and the service URL come from the beamline's DM configuration, so the
    shell running Bluesky must have sourced the DM setup script first.

    Returns:
        Tuple of (workflowProcApi, dmuser). Pass both straight on to dm_run_job().
    """
    # Object that tracks beamline-specific configuration
    configManager = ConfigurationManager.getInstance()
    dmuser, password = configManager.parseLoginFile()
    serviceUrl = configManager.getProcWebServiceUrl()
    # user/password/url info passed to DM API
    workflowProcApi = WorkflowProcApi(dmuser, password, serviceUrl)

    return workflowProcApi, dmuser


def dm_job_log_path() -> Path:
    """<mount_point>/<cycle>/<experiment>/dm_jobs.log -- one line per submitted job.

    At the experiment root, beside data/ and analysis/, rather than inside either:
    both of those are DM-managed trees and this file is ours.
    """
    return Path(f"{expt.mount_point}{expt.cycle_name}/{expt.experiment_name}/dm_jobs.log")


def log_dm_job(job_id: str, file_name: str, filepath: str, machine_name: str, workflow_name: str):
    """Append one line recording a submitted job, so the uuid outlives the terminal.

    The uuid is the only handle on a DM job -- `dmjob.sh <uuid>` is how you ask
    what happened to it -- and until 2026-09-07 it was printed to the session and
    nothing more, so it was gone as soon as the scrollback was.

    A plain text log rather than a field in the NeXus metadata file, for two
    reasons: the metadata file is written BEFORE the job is submitted, so the uuid
    does not exist yet and adding it would mean reopening a file DM has already
    been pointed at; and one greppable file answers "which job was that?" without
    opening sixty HDFs.

        grep A0101 dm_jobs.log                 # the job for one measurement
        awk '!/^#/{print $3}' dm_jobs.log      # every uuid
        dmjob.sh $(tail -1 dm_jobs.log | awk '{print $3}')   # status of the last one

    Never raises. A full disk or a read-only mount must not take down an
    acquisition that has already written its data -- the line is lost, the run
    continues, and the uuid is still on screen.
    """
    path = dm_job_log_path()
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        new = not path.exists()
        with open(path, "a") as handle:
            if new:
                handle.write(
                    "# Submitted DM analysis jobs, appended by id8_common.utils.dm_util.\n"
                    "# date       time      job_uuid  measurement  machine  workflow  data_file\n"
                )
            handle.write(
                f"{stamp}  {job_id}  {file_name}  {machine_name}  {workflow_name}  {filepath}\n"
            )
    except OSError as exc:
        print(f"[dm_util] could not append to {path}: {exc}")


def dm_run_job(workflowProcApi: WorkflowProcApi, dmuser: str, file_name: str):
    """Submit one analysis job to the Data Management system for a finished measurement.

    file_name is the measurement's base path with no suffix; the suffix that the current
    detector actually wrote (.h5, .bin.000, ...) is appended below. Everything else --
    experiment, qmap, workflow, analysis type -- is read from the experiment config, so
    the caller only has to say which file to analyse.
    """

    analysis_machine = expt.analysis_machine
    det_name = expt.det_name

    if analysis_machine == "none":
        # The user turned analysis off for this experiment: write the data, submit nothing.
        pass
    else:
        exp_name = expt.experiment_name
        qmap_file = expt.qmap_file
        workflow_name = expt.workflow_name
        analysis_type = expt.analysis_type
        use_subfolder = expt.use_subfolder

        if det_name == "rigaku3M":
            filepath = f"{file_name}.bin.000"
        elif det_name == "rigaku3M_ftf":
            # Fast transfer now asks the IOC for a .h5 name (the contents were
            # always HDF5 -- verified by \x89HDF magic 2026-09-03, when it was
            # still written as .bin), so the six per-module files are
            # <file_name>.h5.000 .. .h5.005 and this points at the first.
            #
            # The .000 suffix itself is verified, from runs made under the old
            # .bin name; that the IOC appends it the same way to a .h5 name is
            # inferred, not yet observed. Confirm on the first fast-transfer
            # run after this change and correct here if it differs.
            filepath = f"{file_name}.h5.000"
        elif det_name == "rigaku3M_epics":
            filepath = f"{file_name}.h5"
        elif det_name == "eiger4M":
            filepath = f"{file_name}.h5"
        elif det_name == "lambda2M":
            filepath = f"{file_name}.h5"
        elif det_name == "tempus":
            filepath = f"{file_name}.bin"
        else:
            pass

        if analysis_machine == "polaris":
            gpuID = 0
            machine_name = analysis_machine
        elif analysis_machine == "local":
            # "local" does not name a machine -- get_machine_name() picks one of the
            # beamline analysis boxes for us.
            gpuID = -2
            machine_name = get_machine_name()
        else:
            # Any other value is taken as the hostname to run on, as typed.
            gpuID = -2
            machine_name = analysis_machine

        if use_subfolder == "yes":
            use_subfolder_flag = "True"
        elif use_subfolder == "no":   
            use_subfolder_flag = "False"
        else: 
            print("Sub folder options can only be either Yes or No")

        argsDict = {
            "experimentName": exp_name,
            "filePath": filepath,
            "qmap": f"{qmap_file}",
            "analysisMachine": machine_name,
            "gpuID": gpuID,
            "demand": "True",
            "type": analysis_type,
            "saveG2": "False",
            "download": "False",
            "useSubdir": use_subfolder_flag,
            "normalizeFrame": "False"
            #"suffix": "suffix_added",
            # "downloadDirectory": f"/home/8-id-i/{cycle_name}/{exp_name}/analysis/{analysis_type}/"
        }
        
        job = workflowProcApi.startProcessingJob(dmuser, f"{workflow_name}", argsDict=argsDict)
        job_id = job["id"]
        print(f"Job {job_id}")

        log_dm_job(job_id, file_name, filepath, machine_name, workflow_name)

        return job_id
