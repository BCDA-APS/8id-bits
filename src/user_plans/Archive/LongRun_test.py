"""
Functions run multiple samples in a loop.

user defined functions during Su beam time in March 18-24, 2025

"""


import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time
from plan_example_slowmode import eiger_scan_slowmode
from plan_example_fastmode import eiger_scan_fastmode


def loop_eiger_scan(loop_num  = 1, sample_indexes = [1,2]):
    for i in range(loop_num):
        for s_index in sample_indexes:
            yield from eiger_scan_slowmode(sam_index=s_index, 
                    set_temp=25,
                    set_ramp_rate=1,
                    att_level=2000, 
                    acq_time=0.2, 
                    acq_period=1.0,
                    num_frames=1500, 
                    num_rep=1)
            

def loop_eiger_scan_fast(loop_num  = 1, sample_indexes = [1,2]):
    for i in range(loop_num):
        for s_index in sample_indexes:
            yield from eiger_scan_fastmode(sam_index=s_index,
                    set_temp=25, 
                    set_ramp_rate=1,       
                    att_level=1000, 
                    acq_time=0.3,
                    num_frames=1000, 
                    num_rep=6)