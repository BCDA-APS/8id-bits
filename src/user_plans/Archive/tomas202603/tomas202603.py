from apsbits.core.instrument_init import oregistry

from legacy.id8_i.plans.nexus_acq_eiger_int_wei import eiger_acq_int_series_wei
from legacy.id8_i.plans.rheometer_wait import wait_for_mcr
import numpy as np
import time


sample_name = 'CNC_AQ9_01'
repeat_1m = 1

shear_rate = 0.1
shear_time = 30


frame_rate = 0.1
frame_num  = 6000

pv_registers.qmap_file.put('eiger4m_qmap_S90x18_D3x9.hdf')
####pv_registers.qmap_file.put('eiger4m_qmap_S360_D36.hdf')

def rheo_xpcs_measure():

    frame_num = 6000
    frame_rate = 0.1
    att(351)
    wait_for_mcr()

    shear_rate = 0.1
    suffix = f"_shear{round(shear_rate*1e2):03d}_time{shear_time:02}"
    sampe_name_full = sample_name +  suffix

    eiger_acq_int_series_wei(
        acq_time=frame_rate, 
    num_frames=frame_num, 
    num_reps=2, 
    sample_name=sampe_name_full,
    sample_move = True)

    wait_for_mcr()

    shear_rate = 0.5
    suffix = f"_shear{round(shear_rate*1e2):03d}_time{shear_time:02}"
    sampe_name_full = sample_name +  suffix

    eiger_acq_int_series_wei(
        acq_time=frame_rate, 
    num_frames=frame_num, 
    num_reps=2, 
    sample_name=sampe_name_full,
    sample_move = True)


def xpcs_measure():

    frame_rate = 0.01
    frame_num = 1000
    # suffix = f"_shear{round(shear_rate*1e1):02d}_time{shear_time:02}_followup"
    suffix = f"_static_"
    sampe_name_full = sample_name +  suffix

    eiger_acq_int_series_wei(
        acq_time=frame_rate, 
    num_frames=frame_num, 
    num_reps=10, 
    sample_name=sampe_name_full,
    sample_move = True)

def rheo_saxs_measure():

    frame_rate = 0.1
    frame_num = 300
    shear_rate = 0.1

    suffix = f"__steadyshear{round(shear_rate*1e1):02d}"
    sampe_name_full = sample_name +  suffix

    wait_for_mcr()
    eiger_acq_int_series_wei(
        acq_time=frame_rate, 
    num_frames=frame_num, 
    num_reps=1, 
    sample_name=sampe_name_full,
    sample_move = True)


def att_scan():

    
    suffix = f"_static"
    sampe_name_full = sample_name +  suffix  

    att_list = np.array([550, 450]) 
    pos_list = np.arange(att_list.shape[0]) * 0.1 + 356
    frame_rate = 0.1
    frame_num  = 6000

    for idx, att_idx in enumerate(att_list):

        att(att_idx)
        pos = pos_list[idx]
        rheometer.y.move(pos)

        eiger_acq_int_series_wei(
        acq_time=frame_rate, 
        num_frames=frame_num, 
        num_reps=repeat_1m, 
        sample_name=sampe_name_full,
        sample_move = True)



