"""
Simple, modular Ophyd scripts for users.
"""

from datetime import datetime
import time

from apsbits.core.instrument_init import oregistry
from id8_common.plans.set.shutter_att import *

tetramm1 = oregistry["tetramm1"]
tetramm2 = oregistry["tetramm2"]
tetramm3 = oregistry["tetramm3"]
tetramm4 = oregistry["tetramm4"]
pv_registers = oregistry["pv_registers"]
huber = oregistry["huber"]


def tetramm_acq_series(
    det=tetramm3,    
    filename=None,
    num_capture=1,
):
    
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    mount_point = pv_registers.mount_point.get()
    
    full_filename = f"{mount_point}{cycle_name}/{exp_name}/data/{filename}"
    
    det.hdf1.file_name.put(full_filename)
    det.hdf1.num_capture.put(num_capture)

    # showbeam()
    time.sleep(0.1)
    det.hdf1.capture.put(1)

    while True:
        time.sleep(0.1)
        det_status = det.hdf1.capture.get()
        if det_status == 1:
            time.sleep(0.1)
        if det_status == 0:
            break
    # blockbeam()
