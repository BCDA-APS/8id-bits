
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time


def doe_sc_demo():
    for ii in range(100):
    #diff frame rate, same total data time
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.05, 
                                        num_frame=1000, 
                                        num_rep=1, 
                                        att_level=0, 
                                        sample_move = False)
        