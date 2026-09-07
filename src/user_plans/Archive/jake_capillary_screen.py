
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time
from datetime import datetime

def capillary_screen():
        
    sample_list = [0] #[10, 11, 12, 13, 15] #[16, 17, 18, 19, 20, 21, 22, 24, 25]
    att_list = [12000] #[100, 500, 2000, 10000]
    rep = 70

    for sample_no in sample_list:

        print('switching sample position...')
        yield from select_sample(sample_no)
        time.sleep(2)
        print(f'switched to posit%`ion {sample_no}')

        for att_value in att_list:
            yield from att(att_value)

            for ii in range(rep):
                #yield from wait_for_mcr()
                
                yield from eiger_acq_int_series(acq_time=0.1, 
                                                num_frames= 10000, 
                                                num_rep=1, 
                                                sample_move = True)

                #yield from eiger_acq_ext_trig(acq_time=0.03, 
                                               # acq_period=0.15,
                                                #num_frames= 1000, 
                                               # num_rep=1, 
                                               # sample_move = True)
        
        #yield from att(100000)
        #yield from eiger_acq_int_series(acq_time=1.0, 
        #                                num_frames=1000, 
        #                                num_rep=1, 
        #                                sample_move = True)
        
def jake_test():
    yield from eiger_acq_int_series(acq_time=0.001, 
                                num_frames=1000, 
                                num_rep=1, 
                                sample_move = True)
                