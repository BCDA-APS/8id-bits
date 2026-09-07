"""
DM code from Hannah Parraga.
Set up DM and submit jobs
"""

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
        print(f"Job {job['id']}")
