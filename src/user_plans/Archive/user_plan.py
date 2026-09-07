
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np

def grid_temp_scan(sample_index = None,
                   sample_name = None,
                   acq_rep = 5,
                   set_temp = None,
                   temp_rate = 1,
                   temp_wait = 300,
                   samx_cen = None,
                   samy_cen = None,
                   acq_time = 0.00025,
                   acq_period = 0.00025,
                   num_images = 4000
                                     
                   ):

    # yield from select_sample(sample_index)

    # temp_zone_index = int(np.floor(sample_index/3))+1
    temp_zone_index = (sample_index-1)//3+1
    
    yield from set_qnw(temp_zone_index, set_temp, True, temp_rate)
    print('Target temperature reached')
    yield from bps.sleep(temp_wait)
    
    samx_num = 20
    samy_num = 40

    samx_list = np.linspace(samx_cen-0.5, samx_cen+0.5, num=samx_num)
    samy_list = np.linspace(samy_cen-0.5, samy_cen+0.5, num=samy_num)
 
    for ii in range(acq_rep): 
        
        pos_index = np.mod(ii,samx_num*samy_num)
        yield from bps.mv(
            sample.x, samx_list[np.mod(pos_index,samx_num)],
            sample.y, samy_list[int(np.floor(pos_index/samx_num))]
        )

        open_shutter()

        yield from xpcs_bdp_demo_plan(sample_name,
                            header='E030', 
                            analysisMachine='adamite',
                            qmap_file='/gdata/dm/8IDI/2024-3/comm202410/data/eiger4m_qmap_1025_s360_d36_linear_KJP.h5', 
                            acquire_time=acq_time, 
                            acquire_period=acq_period, 
                            num_images=num_images, 
                            wf_gpuID=0)

        close_shutter()

        yield from bps.sleep(temp_wait) # sleep for 5 mins


def temp_sweep(sample_index = None,
                   sample_name_head = None,
                   acq_rep = 5,
                   samx_cen = None,
                   samy_cen = None,
                   ):

    temp_zone_index = (sample_index-1)//3+1

    for temp in np.arange(10, 25, 0.5):
        print(f'Waiting for {temp}C...')
        yield from set_qnw(temp_zone_index, temp, True, 0.5)
        print('...Target temperature reached')
        
        temp_wait = 600
        yield from bps.sleep(temp_wait) # sleep and wait for 5 mins

        eq_time = 10

        while eq_time>0:
            print(10 - eq_time)
            # yield from bps.sleep(60)
            samx_num = 20
            samy_num = 20

            samx_list = np.linspace(samx_cen-0.5, samx_cen+0.5, num=samx_num)
            samy_list = np.linspace(samy_cen-0.5, samy_cen+0.5, num=samy_num)

            sample_name = f'{sample_name_head}_{temp}_{10-eq_time}'

            for ii in range(acq_rep): 
                
                pos_index = np.mod(ii,samx_num*samy_num)
                yield from bps.mv(
                    sample.x, samx_list[np.mod(pos_index,samx_num)],
                    sample.y, samy_list[int(np.floor(pos_index/samx_num))]
                )



                open_shutter()

                yield from xpcs_bdp_demo_plan(sample_name,
                                    header='K002', 
                                    analysisMachine='adamite',
                                    qmap_file='/gdata/dm/8IDI/2024-3/comm202410/data/eiger4m_qmap_1025_s360_d36_linear_KJP.h5', 
                                    acquire_time=0.00025, 
                                    acquire_period=0.00025, 
                                    num_images=4000, 
                                    wf_gpuID=0)

                close_shutter()


            
            
            eq_time = eq_time - 1


def framerate_scan():

    # yield from select_sample(3)

    open_shutter()

    sample_name = 'Cys1_40'
    shear_rate = 'steady'
    qmap_file = '/gdata/dm/8IDI/2024-3/comm202410/data/eiger4m_qmap_1107_s360_d36_phi1_full_HHR.h5'
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_0p00025',
                        header='B020', 
                        analysisMachine='adamite',
                        qmap_file=qmap_file, 
                        acquire_time=0.00025, 
                        acquire_period=0.00025, 
                        num_images=1000, 
                        wf_gpuID=-2)
    
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate +  '_rate_0p001',
                        header='B021', 
                        analysisMachine='adamite',
                        qmap_file=qmap_file, 
                        acquire_time=0.001, 
                        acquire_period=0.001, 
                        num_images=1000, 
                        wf_gpuID=-2)    
    
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_0p01',
                    header='B022', 
                    analysisMachine='adamite',
                    qmap_file=qmap_file, 
                    acquire_time=0.01, 
                    acquire_period=0.01, 
                    num_images=1000, 
                    wf_gpuID=-2)
    
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_0p1',
                    header='B022', 
                    analysisMachine='adamite',
                    qmap_file=qmap_file, 
                    acquire_time=0.1, 
                    acquire_period=0.1, 
                    num_images=1000, 
                    wf_gpuID=-2)
    
    yield from xpcs_bdp_demo_plan(sample_name + '_' + shear_rate + '_rate_1',
                    header='B022', 
                    analysisMachine='adamite',
                    qmap_file=qmap_file, 
                    acquire_time=1, 
                    acquire_period=1, 
                    num_images=1000, 
                    wf_gpuID=-2)

    close_shutter()