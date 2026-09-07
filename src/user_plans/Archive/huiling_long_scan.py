
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time

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





 
    