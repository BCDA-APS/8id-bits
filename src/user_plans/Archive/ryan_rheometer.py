
from apsbits.core.instrument_init import oregistry

from legacy.id8_i.plans.nexus_acq_eiger_int import eiger_acq_int_series
from legacy.id8_i.plans.sample_info_unpack import select_sample
from legacy.id8_i.plans.scan_8idi import att
from legacy.id8_i.plans.rheometer_wait import wait_for_mcr
import numpy as np
import time

rheometer = oregistry["rheometer"]


def rheometer_PP():
        
    sample_list = [0]
    att_list = [500]
    height_list = 349.1 + np.linspace(-0.2,0.2,num = 30)
    rep = 8
    x_list = 11.1 + np.linspace(-0.2,0.2,rep)
    
    sleep_list = np.logspace(0,2.4,rep)
    # rep = 20

    for sample_no in sample_list:

        print('switching sample position...')
        select_sample(sample_no)
        time.sleep(2)
        print(f'switched to position {sample_no}')

        for att_value in att_list:
            att(att_value)

            rheometer.y.move(height_list[0])

            for ii in height_list:
                wait_for_mcr()
                wait_for_mcr()

                rheometer.y.move(ii)

                for jj in range(rep):
                    
                    x_pos = x_list[jj]

                    time.sleep(sleep_list[jj])
                    
                    rheometer.x.move(x_pos)

                    eiger_acq_int_series(acq_time=0.0025, 
                    num_frames= 4000, 
                    num_reps=1, 
                    # wait_time=5,
                    sample_move = False)
                    
                wait_for_mcr()
        
def ryan_test():
    eiger_acq_int_series(acq_time=0.001, 
                        num_frames=1000, 
                        num_reps=1, 
                        sample_move = True)
                