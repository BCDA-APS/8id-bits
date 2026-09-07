
from bluesky import plan_stubs as bps

from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
from aps_8id_bs_instrument.plans.nexus_utils import create_nexus_format_metadata
from dm.proc_web_service.api.workflowProcApi import WorkflowProcApi
from dm.common.utility.configurationManager import ConfigurationManager

import time

def submit_Nexus_DM():

    while True:
        bluesky_start = pv_registers.start_bluesky.get()
        if bluesky_start == 'Yes':

            # DM workflow setup. 
            # configManager is an object that tracks beamline-specific configuration.
            # In WorkflowProcApi, user/password/url info is passed to DM API
            configManager = ConfigurationManager.getInstance()  
            dmuser, password = configManager.parseLoginFile()
            serviceUrl = configManager.getProcWebServiceUrl()
            workflowProcApi = WorkflowProcApi(dmuser, password, serviceUrl) 

            # Spec will need to write these fields in StrReg.
            # exp_name, workflow_name, analysis_machine need to be written only once per user.
            # metadata_fname and filename needs to be written per measurement.
            # Change qmap file when needed.
            exp_name = pv_registers.experiment_name.get()
            workflow_name = pv_registers.workflow_name.get()
            analysis_machine = pv_registers.analysis_machine.get()
            qmap_file = pv_registers.qmap_file.get()
            metadata_fname = pv_registers.metadata_full_path.get()
            filename = pv_registers.file_name.get()
            det_name = pv_registers.det_name.get()

            if det_name == 'eiger4M':
                det = eiger4M
            elif det_name == 'rigaku3M':
                det = rigaku3M
            elif det_name == 'lambda2M':
                det = lambda2M
            else:
                raise ValueError("Detector Name Invalid. Must be eiger4M, rigaku3M, lambda2M.")    

            # Miaoqi's code that writes the metadata
            create_nexus_format_metadata(metadata_fname, det=det)
            # for ii in range(5):
            #     print(ii)
            #     time.sleep(1.0)
            # time.sleep(5.0)

            # Code that starts DM workflow
            argsDict = {"experimentName": exp_name, 
                        "filePath": f"{filename}.h5", 
                        "qmap": f"{qmap_file}",
                        "analysisMachine": f"{analysis_machine}",
                        "gpuID": -2,
                        "type": "Multitau",
                        }
            job = workflowProcApi.startProcessingJob(dmuser, f"{workflow_name}", argsDict=argsDict)
            print(f"Job {job['id']} processing {filename}")
            print(filename)
            pv_registers.start_bluesky.put('No')
            time.sleep(1.0)
        else:
            time.sleep(1.0)


def test_sleep():
    yield from bps.sleep(5.0)
