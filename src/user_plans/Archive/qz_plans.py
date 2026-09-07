
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time

def brian_screen():

    sample_list = [10, 11, 12, 13, 14, 15, 16, 17, 18] 
    # att_list = [4, 8, 12, 16, 20]
    att_list = [0, 2, 4, 6, 8, 10, 12]
    rep = 10

    for att_value in att_list:

        for sample_no in sample_list:

            print('switching sample position...')
            yield from select_sample(sample_no)
            time.sleep(2)
            print(f'switched to position {sample_no}')
        
            yield from eiger_acq_int_series(eiger4M, 
                                            acq_period=0.00025, 
                                            num_frame=8000, 
                                            num_rep=rep, 
                                            att_level=att_value, 
                                            sample_move = True)
            
def sanjeeva_screen():
        
    sample_list = [19, 20, 21, 22, 23, 24, 25, 26, 27]
    temp_list = [55, 60, 70, 75]
    rep = 30

    for temp in temp_list:
        
        yield from set_qnw(qnw_number = 3, setpoint = temp, wait = True, ramprate = 5)
        yield from bps.sleep(120)

        for sample_no in sample_list:

            print('switching sample position...')
            yield from select_sample(sample_no)
            time.sleep(2)
            print(f'switched to position {sample_no}')
            
            yield from eiger_acq_int_series(eiger4M, 
                                            acq_period=0.00025, 
                                            num_frame=8000, 
                                            num_rep=rep, 
                                            att_level=4, 
                                            sample_move = True)
            
    yield from set_qnw(qnw_number = 3, setpoint = 25, wait = True, ramprate = 5)
            

def overnight_scan():

    yield from brian_screen()
    yield from sanjeeva_screen()
        


