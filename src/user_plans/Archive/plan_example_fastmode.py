"""

user defined functions during Su beam time in March 2025

"""

import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time


def eiger_scan_fastmode(sam_index=None, 
               set_temp=None, 
               set_ramp_rate=None, 
               att_level=10, 
               acq_time=0.2, 
               num_frames=10, 
               num_rep=3):

    yield from bps.mv(filter_8idi.attenuation_set, att_level)
    yield from bps.sleep(5)
    
    yield from select_sample(sam_index)

    yield from temp_ramp(set_ramp_rate)
    yield from te(set_temp, wait=False)

    yield from eiger_acq_int_series(acq_time=acq_time,
                                    num_frames=num_frames, 
                                    num_rep=num_rep, 
                                    wait_time=0, 
                                    process=True, 
                                    sample_move=True)