
from apsbits.core.instrument_init import oregistry

from legacy.id8_i.plans.nexus_acq_eiger_int_wei import eiger_acq_int_series_wei
from legacy.id8_i.plans.rheometer_wait import wait_for_mcr
import numpy as np
import time

rheometer = oregistry["rheometer"]
freq_list = np.array([1, 0.667, 0.5, 0.4, 0.333, 0.286, 0.25, 0.2222, 0.2, 0.182, 0.167, 0.154, 0.143, 0.1333, 0.125, 0.118, 0.1110, 0.105, 0.1, 0.0909, 0.0833, 0.0769, 0.0714, 0.0667, 0.0625, 0.0588, 0.0556, 0.0526, 0.05])

frame_rate_list = np.ones(freq_list.shape[0])
frame_rate_list[:9]   = 0.001
frame_rate_list[9:19] = 0.002
frame_rate_list[19:]  = 0.005

frame_num_list     = np.ones(freq_list.shape[0]) * 10000


time_per_scan_list = frame_rate_list * frame_num_list

time_list = 1 / freq_list * 30
repeat_1ms_list = np.floor((time_list - 25) / time_per_scan_list).astype(np.int16)
repeat_1ms_list[0] = 2
repeat_1ms_list[1] = 3
repeat_1ms_list[2] = 4
repeat_1ms_list[7] = 10
repeat_1ms_list[8] = 11
repeat_1ms_list[12] = 8
repeat_1ms_list[13] = 9
repeat_1ms_list[14] = 9
repeat_1ms_list[15] = 10
repeat_1ms_list[16] = 11
repeat_1ms_list[17] = 12
repeat_1ms_list[18] = 12
repeat_1ms_list[24] = 8

residual_time = time_list - (time_per_scan_list * repeat_1ms_list)

shear_rate = 0.0


sample_name = 'P52'

def rheo_xpcs_measure():


    origin_pos = 352
    #hexapodc_pos = -3
    #rheometer.yaw.move(hexapodc_pos)
    rheometer.y.move(origin_pos)
    wait_for_mcr()
    wait_for_mcr()
    wait_for_mcr()


    for interval_idx, freq in enumerate(freq_list):

        wait_for_mcr()
        frame_rate = frame_rate_list[interval_idx]
        frame_num  = round(frame_num_list[interval_idx])
        repeat_1m  = round(repeat_1ms_list[interval_idx])


        ##UNCOMMENT AND RE-ADD SHEAR RATE FOR SHEAR TESTING!!!##

        suffix = f"Int{(interval_idx+1):02d}_freq{round(freq*1e3):04d}_shear{round(shear_rate*1e1):02d}"

        sampe_name_full = sample_name + '_' + suffix

        eiger_acq_int_series_wei(acq_time=frame_rate, 
        num_frames=frame_num, 
        num_reps=repeat_1m, 
        sample_name=sampe_name_full,
        sample_move = True)

        #hexapodc_pos += 0.2
        #rheometer.yaw.move(hexapodc_pos)
        origin_pos -= 0.01
        rheometer.y.move(origin_pos)



