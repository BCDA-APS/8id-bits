
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np



def framerate_scan():

    # yield from select_sample(3)

    open_shutter()

    sample_name = 'Cys1_40'
    shear_rate = 'steady'
    qmap_file = '/gdata/dm/8IDI/2024-3/comm202411/data/eiger4m_qmap_1107_s360_d36_phi1_full_HHR.h5'
    acquire_time=0.00025

    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_0p00025',
                        header='B020', 
                        analysisMachine='adamite',
                        qmap_file=qmap_file, 
                        acquire_time = acquire_time, 
                        acquire_period=0.00025, 
                        num_images=1000, 
                        wf_gpuID=-2)


    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate +  '_rate_0p001',
                        header='B021', 
                        analysisMachine='adamite',
                        qmap_file=qmap_file, 
                        acquire_time=acquire_time, 
                        acquire_period=0.001, 
                        num_images=1000, 
                        wf_gpuID=-2)    
    
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_0p01',
                    header='B022', 
                    analysisMachine='adamite',
                    qmap_file=qmap_file, 
                    acquire_time=acquire_time, 
                    acquire_period=0.01, 
                    num_images=1000, 
                    wf_gpuID=-2)
    
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_0p1',
                    header='B022', 
                    analysisMachine='adamite',
                    qmap_file=qmap_file, 
                    acquire_time=acquire_time, 
                    acquire_period=0.1, 
                    num_images=1000, 
                    wf_gpuID=-2)
    
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_1',
                    header='B022', 
                    analysisMachine='adamite',
                    qmap_file=qmap_file, 
                    acquire_time=acquire_time, 
                    acquire_period=1, 
                    num_images=1000, 
                    wf_gpuID=-2)


    close_shutter()