
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time
from datetime import datetime


def capillary_overnight():
        
#    sample_list = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19] 
    sample_list = [12, 14]

    att_rigaku = [1000, 1000]
    # att_eiger1 = [1000,1000]
    # att_eiger2 = [100000,100000]


    for ii in range(len(sample_list)):

        print('switching sample position...')
        yield from select_sample(sample_list[ii])
        time.sleep(2)
        print(f'switched to position {sample_list[ii]}')

        # print('running Rigaku3M acquisition ...')
        # yield from att(att_rigaku[ii])
        # yield from select_detector('rigaku')
        # yield from rigaku_acq_ZDT_series(acq_time = 2e-5, 
        #                                 num_frames = 100000, 
        #                                 num_rep = 1, 
        #                                 sample_move = True)
        
        print('running Eiger4M acquisition ...')
        # yield from att(att_eiger1[ii])
        yield from att(10000)
        yield from select_detector('eiger')
        yield from eiger_acq_int_series(acq_time = 0.001, 
                                        num_frames = 10000, 
                                        num_rep = 5, 
                                        sample_move = True)
        # yield from att(att_eiger1[ii])
        yield from att(100000)
        yield from select_detector('eiger')
        yield from eiger_acq_int_series(acq_time = 0.1, 
                                         num_frames = 5000, 
                                         num_rep = 3, 
                                         sample_move = True)


        