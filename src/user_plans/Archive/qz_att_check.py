
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time


def qz_long():

    yield from post_align()

    sam_list = [2,3]
    att_list = [0,4,7,9,12]
    num_r = 500

    for sam in sam_list:
        yield from select_sample(sam)

        for att in att_list:

            yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=8000, 
                                        num_rep=num_r, 
                                        att_level=att, 
                                        sample_move = True)

    yield from select_sample(7)
    yield from eiger_acq_int_series(eiger4M, 
                                 acq_period=0.001, 
                                 num_frame=1000, 
                                 num_rep=50, 
                                 att_level=20, 
                                 sample_move=True)
    yield from eiger_acq_flyscan(eiger4M, 
                                 acq_period=0.001, 
                                 num_frame=1000, 
                                 num_rep=num_r*2, 
                                 att_level=20, 
                                 sample_move=False, 
                                 flyspeed=0.05)
