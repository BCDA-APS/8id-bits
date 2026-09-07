
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time


def brian_overnight(num_r, att):

    yield from post_align()
    
    yield from select_sample(5)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.00025, 
                                num_frame=8000, 
                                num_rep=num_r, 
                                att_level=att, 
                                sample_move = True)

    yield from select_sample(8)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.00025, 
                                num_frame=8000, 
                                num_rep=num_r, 
                                att_level=att, 
                                sample_move = True)
    
def sanjeeva_overnight(num_r):

    yield from post_align()
    
    yield from select_sample(16)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=0, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)    
    
    yield from select_sample(17)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=0, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)  
    
    yield from select_sample(18)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=0, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)  
    
    yield from select_sample(19)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)
    
    yield from select_sample(20)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)
    
    yield from select_sample(22)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=0, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)  
    
    yield from select_sample(23)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=0, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)  
    
    yield from select_sample(24)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=0, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True)  
    
    yield from select_sample(25)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=33, 
                                sample_move = True)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True) 
    
    yield from select_sample(26)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.0025, 
                                num_frame=4000, 
                                num_rep=num_r, 
                                att_level=88, 
                                sample_move = True) 

def qz_long():

    yield from post_align()
    
    # yield from select_sample(5)
    # yield from eiger_acq_int_series(eiger4M, 
    #                             acq_period=0.001, 
    #                             num_frame=1000, 
    #                             num_rep=25, 
    #                             att_level=0, 
    #                             sample_move = True)

    # yield from select_sample(8)
    # yield from eiger_acq_int_series(eiger4M, 
    #                             acq_period=0.001, 
    #                             num_frame=1000, 
    #                             num_rep=25, 
    #                             att_level=0, 
    #                             sample_move = True)
    
    yield from select_sample(11)
    yield from eiger_acq_int_series(eiger4M, 
                                acq_period=0.001, 
                                num_frame=1000, 
                                num_rep=25, 
                                att_level=0, 
                                sample_move = True)


def overnight_scan():

    # yield from brian_overnight(num_r=200, att=0)
    # yield from sanjeeva_overnight(num_r=50)

    yield from set_qnw(qnw_number = 2, setpoint = 70, wait = False, ramprate = 5)
    yield from set_qnw(qnw_number = 3, setpoint = 70, wait = False, ramprate = 5)
    yield from bps.sleep(300)

    yield from sanjeeva_overnight(num_r=50)

    yield from set_qnw(qnw_number = 2, setpoint = 25, wait = False, ramprate = 5)
    yield from set_qnw(qnw_number = 3, setpoint = 25, wait = False, ramprate = 5)
        


 
    