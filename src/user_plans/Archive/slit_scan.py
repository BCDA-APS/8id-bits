
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time


def sl4_vary_eiger():

    yield from bps.mv(sl4.h.size, 0.075)
    yield from bps.mv(sl4.v.size, 0.075)
    yield from eiger_acq_int_series(eiger4M, 
                                    acq_period=0.01, 
                                    num_frame=1000, 
                                    num_rep=3, 
                                    att_level=0, 
                                    sample_move = False)
    
    # yield from bps.mv(sl4.h.size, 0.100)
    # yield from bps.mv(sl4.v.size, 0.100)
    # yield from eiger_acq_int_series(eiger4M, 
    #                                 acq_period=0.01, 
    #                                 num_frame=1000, 
    #                                 num_rep=3, 
    #                                 att_level=0, 
    #                                 sample_move = False)
    
    # yield from bps.mv(sl4.h.size, 0.125)
    # yield from bps.mv(sl4.v.size, 0.125)
    # yield from eiger_acq_int_series(eiger4M, 
    #                                 acq_period=0.01, 
    #                                 num_frame=1000, 
    #                                 num_rep=3, 
    #                                 att_level=0, 
    #                                 sample_move = False)

    # yield from bps.mv(sl4.h.size, 0.150)
    # yield from bps.mv(sl4.v.size, 0.150)
    # yield from eiger_acq_int_series(eiger4M, 
    #                                 acq_period=0.01, 
    #                                 num_frame=1000, 
    #                                 num_rep=3, 
    #                                 att_level=0, 
    #                                 sample_move = False)

    # yield from bps.mv(sl4.h.size, 0.175)
    # yield from bps.mv(sl4.v.size, 0.175)
    # yield from eiger_acq_int_series(eiger4M, 
    #                                 acq_period=0.01, 
    #                                 num_frame=1000, 
    #                                 num_rep=3, 
    #                                 att_level=0, 
    #                                 sample_move = False)