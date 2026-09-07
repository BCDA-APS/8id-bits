
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time

def KJP_scan_code():

    sample_list = [3, 4, 5, 6, 7, 8] #sample position numbers 1 - 9
    
    temps = [5, 10, 15] #temperature range

    temp_positions = [0, 1, 2] #temp positions [0, 1, 2] for [1-3, 4-6, 7-9]

    wait_time = 60*10 #time in seconds to wait for temperature equilibration

    #set each temp holder to intial temperature
    for temp_pos in temp_positions:
        yield from set_qnw(qnw_number = temp_pos+1, setpoint = temps[0], wait = False, ramprate = 1)
    
    
    #wait for initial temperature
    print(f'waiting {wait_time} s for temperature')
    #time.sleep(wait_time)


    for next_temp in temps[1:]:
        for temp_pos in temp_positions: #[0, 1, 2]
            for sample_no in sample_list:
                
                if (sample_no-1)//3 != temp_pos: 
                    #skip if sample not in current temp_pos for the loop
                    continue
                
                #switch to sample position
                print('switching sample position...')
                yield from select_sample(sample_no)
                time.sleep(2)
                print(f'switched to position {sample_no}')
                
                yield from eiger_acq_int_series(eiger4M, acq_period=0.00025, 
                                        num_frame=4000, num_rep=30, att_level=0, 
                                         sample_move = True)
                
                yield from eiger_acq_int_series(eiger4M, acq_period=0.001, 
                                        num_frame=10000, num_rep=10, att_level=4, 
                                         sample_move = True)
            
            #set current temp_pos to next temperature
            print(f'Setting {temp_pos+1} to {next_temp} C')
            yield from set_qnw(qnw_number = temp_pos+1, setpoint = next_temp, wait = False, ramprate = 1)
            
            #if the temp_pos is the first in measurement set, record the time
            if temp_pos == temp_positions[0]:
                temp_time = time.time()
        
        #if after all temp_pos measurements are done, the wait time is not reached, wait remaining time
        elap_time = time.time() - temp_time
        remain_time = wait_time - elap_time
        if remain_time>0:
            print(f'waiting {remain_time} s for temperature')
            time.sleep(remain_time)
        else:
            continue
    
    #measure scans for the final temperature
    for temp_pos in temp_positions:

            for sample_no in sample_list:
                if (sample_no-1)//3 != temp_pos:
                    continue
                print('switching sample position...')
                yield from select_sample(sample_no)
                time.sleep(2)
                print(f'switched to {sample_no}')
                
                yield from eiger_acq_int_series(eiger4M, acq_period=0.00025, 
                                        num_frame=4000, num_rep=30, att_level=0, 
                                         sample_move = True)
                
                yield from eiger_acq_int_series(eiger4M, acq_period=0.001, 
                                        num_frame=10000, num_rep=10, att_level=4, 
                                         sample_move = True)
                


    #reset temperature to 25C
    for temp_pos in temp_positions:
        yield from set_qnw(qnw_number = temp_pos+1, setpoint = 25, wait = False, ramprate = 1)

    print('finished')
    

# def KJP_scan_code2():
#     #time.sleep(10)
#     print('switching sample position...')
#     yield from select_sample(4)
#     print(f'switched to {4}')
#     time.sleep(2)
#     yield from eiger_acq_int_series(eiger4M, acq_period=0.00025, 
#                                         num_frame=4000, num_rep=2, att_level=0, 
#                                          sample_move = True)

#     #time.sleep(10)
#     print('switching sample position...')
#     yield from select_sample(5)
#     print(f'switched to {5}')
#     time.sleep(2)
#     yield from eiger_acq_int_series(eiger4M, acq_period=0.00025, 
#                                         num_frame=4000, num_rep=2, att_level=0, 
#                                          sample_move = True)
    
#     #time.sleep(10)
#     print('switching sample position...')
#     yield from select_sample(6)
#     print(f'switched to {6}')
#     time.sleep(2)
#     yield from eiger_acq_int_series(eiger4M, acq_period=0.00025, 
#                                         num_frame=4000, num_rep=2, att_level=0, 
#                                          sample_move = True)
    
#     #time.sleep(10)
#     print('switching sample position...')
#     yield from select_sample(7)
#     print(f'switched to {7}')
#     time.sleep(2)
#     yield from eiger_acq_int_series(eiger4M, acq_period=0.00025, 
#                                         num_frame=4000, num_rep=2, att_level=0, 
#                                          sample_move = True)
    
#     #time.sleep(10)
#     print('switching sample position...')
#     yield from select_sample(8)
#     print(f'switched to {8}')
#     time.sleep(2)
#     yield from eiger_acq_int_series(eiger4M, acq_period=0.00025, 
#                                         num_frame=4000, num_rep=2, att_level=0, 
#                                          sample_move = True)

    