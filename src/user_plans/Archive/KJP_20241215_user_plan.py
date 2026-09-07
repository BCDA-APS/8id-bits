import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *


def KJP_plan():

    yield from post_align()

    for i in range(1, 10):

        yield from select_sample(i)

        yield from eiger_acq_ext_trig(eiger4M, acq_time = 0.001,acq_period = 0.1, 
                                      num_frame = 1000,
                                        num_rep = 10, att_level = 0, sample_move = True)
        
