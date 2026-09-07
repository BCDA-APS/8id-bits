
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time
from datetime import datetime



def kofi_set_time_cycling_with_att_list(cycle1=[10,11,12],att1=[50,100,200],cycle2=[13,14,16,17,18],att2=[50,100,50,100,200], total_loops=3000):

    def run_cycle(sample_list,att_list):
        for j,sample in enumerate(sample_list):
            yield from bps.mv(filter_8idi.attenuation_set, att_list[i])
            yield from bps.sleep(2)
            yield from bps.mv(filter_8idi.attenuation_set, att_list[i])
            yield from bps.sleep(2)
            yield from select_sample(sample)
            yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                    num_frame=100000, 
                                    num_rep=1)
        




    start_cycle_2 = datetime(2025,3,17,7,0,0)
    end_measurement = datetime(2025,3,17,8,0,0)

    yield from run_cycle(cycle1,att1)
    yield from run_cycle(cycle2,att2)


    for i in range(0,total_loops):

        if datetime.now()>end_measurement:
            break
        
        if datetime.now()<start_cycle_2:
            att_list = att1
            yield from run_cycle(cycle1,att_list)

        if datetime.now()>start_cycle_2:
            att_list =att2
            yield from run_cycle(cycle2,att_list)

def kofi_set_time_cycling(cycle1=[10,11,12],cycle2=[13,14,16,17,18], total_loops=3000, att=150):

    def run_cycle(sample_list):
        for j in sample_list:
            yield from select_sample(sample)
            yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                    num_frame=100000, 
                                    num_rep=1)
        
    yield from bps.mv(filter_8idi.attenuation_set, att)
    yield from bps.sleep(2)
    yield from bps.mv(filter_8idi.attenuation_set, att)
    yield from bps.sleep(2)



    start_cycle_2 =  datetime(2025,3,17,0)
    end_measurement = datetime(2025,3,17,8,0,0)

    run_cycle(cycle1)
    run_cycle(cycle2)


    for i in total_loops:

        if datetime.now()>end_measurement:
            break
        
        run_cycle(cycle1)

        if datetime.now()>start_cycle_2:
            run_cycle(cycle2)
        


        

    # for j in range(0,total_loops):
    #     current_time = datetime.now()

    



def kofi_saxs(att=100,num_rep=10):
    # sample_list=[5,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
    sample_list=[10,11,12,13,14,15]

    yield from bps.mv(filter_8idi.attenuation_set, att)
    yield from bps.sleep(2)
    yield from bps.mv(filter_8idi.attenuation_set, att)
    yield from bps.sleep(2)

    for sample in sample_list:
        yield from select_sample(sample)
        yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                num_frame=100000, 
                                num_rep=num_rep)
    

def kofi_screen(sample_list, att_list, rep):

    if len(sample_list) != len(att_list):
        raise ValueError('Sample list and atteunation list must be equal')
    else:
        for ii in range(len(sample_list)):
            print('switching sample position...')
            yield from select_sample(sample_list[ii])
            time.sleep(2)
            print(f'switched to position {sample_list[ii]}')

            yield from bps.mv(filter_8idi.attenuation_set, att_list[ii])
            yield from bps.sleep(2)
            yield from bps.mv(filter_8idi.attenuation_set, att_list[ii])
            yield from bps.sleep(2)
        
            yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                                num_frame=100000, 
                                                num_rep=rep)
            
            

def kofi_temp_ramp(qnw_env, set_temp=145, ramp_rate=5, rep=800, att_level=100):

    if qnw_env == "qnw_env1":
        yield from bps.mv(qnw_env1.ramprate, ramp_rate)
        yield from bps.sleep(2)
        yield from bps.mv(qnw_env1.setpoint, set_temp)
    if qnw_env == "qnw_env2":
        yield from bps.mv(qnw_env2.ramprate, ramp_rate)
        yield from bps.sleep(2)
        yield from bps.mv(qnw_env2.setpoint, set_temp)
    if qnw_env == "qnw_env3":
        yield from bps.mv(qnw_env3.ramprate, ramp_rate)
        yield from bps.sleep(2)
        yield from bps.mv(qnw_env3.setpoint, set_temp)

    yield from bps.mv(filter_8idi.attenuation_set, att_level)
    yield from bps.sleep(2)
    yield from bps.mv(filter_8idi.attenuation_set, att_level)
    yield from bps.sleep(2)

    yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                    num_frame=100000, 
                                    num_rep=rep)
    

def kofi_round_robin(sample_list, att_list, cycle_number):

    if len(sample_list) != len(att_list):
        raise ValueError('Sample list and attenuation list must be equal')
    else:

        for j in range(0,cycle_number):
            for ii,sample in enumerate(sample_list):
                print('switching sample position...')
                yield from select_sample(sample)
                print(f'switched to position {sample}')
                if ii ==0:
                    print(f'Measuring standard for cycle {j+1}')

                yield from bps.mv(filter_8idi.attenuation_set, att_list[ii])
                yield from bps.sleep(2)
                yield from bps.mv(filter_8idi.attenuation_set, att_list[ii])
                yield from bps.sleep(2)
            
                yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                                num_frame=100000,
                                                num_rep=1)
                
                
            

def kofi_rep_with_delay(sample_index, wait_time=10, rep=1440, att_level=100):

    yield from select_sample(sample_index)

    yield from bps.mv(filter_8idi.attenuation_set, att_level)
    yield from bps.sleep(2)
    yield from bps.mv(filter_8idi.attenuation_set, att_level)
    yield from bps.sleep(2)

    yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                    num_frame=100000, 
                                    num_rep=rep, 
                                    wait_time=wait_time)

def kofi_round_robin_one_attn(sample_list, cycle_number,attn=100):
    yield from bps.mv(filter_8idi.attenuation_set, att)
    yield from bps.sleep(2)
        
    for j in range(0,cycle_number):
        for ii,sample in enumerate(sample_list):
            print('switching sample position...')
            yield from select_sample(sample)
            print(f'switched to position {sample}')


            # yield from bps.sleep(10)
            yield from rigaku_acq_ZDT_series(acq_time=2e-5, 
                                            num_frame=100000,
                                            num_rep=5,sample_move=True)
        print(f'finishing cycle {j+1}')