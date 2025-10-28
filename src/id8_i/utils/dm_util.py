"""
DM code from Hannah Parraga.
Set up DM and submit jobs
"""

from apsbits.core.instrument_init import oregistry
from dm.common.utility.configurationManager import ConfigurationManager
from dm.proc_web_service.api.workflowProcApi import WorkflowProcApi

from .util_8idi import get_machine_name

pv_registers = oregistry["pv_registers"]


def dm_setup() -> tuple:
    """Set up the Data Management workflow API.

    Args:
        process: Whether to initialize the workflow API

    Returns:
        Tuple containing (workflowProcApi, dmuser) if process is True,
        otherwise (None, None)
    """
    # Object that tracks beamline-specific configuration
    configManager = ConfigurationManager.getInstance()
    dmuser, password = configManager.parseLoginFile()
    serviceUrl = configManager.getProcWebServiceUrl()
    # user/password/url info passed to DM API
    workflowProcApi = WorkflowProcApi(dmuser, password, serviceUrl)

    return workflowProcApi, dmuser


def dm_run_job(workflowProcApi: WorkflowProcApi, dmuser: str):
    analysis_machine = pv_registers.analysis_machine.get()
    det_name = pv_registers.det_name.get()

    if analysis_machine == "none":
        pass
    else:
        exp_name = pv_registers.experiment_name.get()
        qmap_file = pv_registers.qmap_file.get()
        workflow_name = pv_registers.workflow_name.get()
        analysis_machine = pv_registers.analysis_machine.get()
        analysis_type = pv_registers.analysis_type.get()
        file_name = pv_registers.file_name.get()

        if det_name == "rigaku3M":
            filepath = f"{file_name}.bin.000"
        elif det_name == "eiger4M":
            filepath = f"{file_name}.h5"
        elif det_name == "tempus":
            filepath = f"{file_name}.bin"
        else:
            pass

        if analysis_machine == "polaris":
            gpuID = 0
            machine_name = analysis_machine
        elif analysis_machine == "local":
            gpuID = -2
            machine_name = get_machine_name()
        else:
            gpuID = -2
            machine_name = analysis_machine

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
            # "downloadDirectory": f"/home/8-id-i/{cycle_name}/{exp_name}/analysis/{analysis_type}/"
        }

        job = workflowProcApi.startProcessingJob(dmuser, f"{workflow_name}", argsDict=argsDict)
        print(f"Job {job['id']} processing {file_name}")
