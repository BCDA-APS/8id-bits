
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time


def swan_scan_repeat():
    #diff frame rate, same total data time
    for i in range(10):
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=4000, 
                                        num_rep=50, 
                                        att_level=0, 
                                        sample_move = False)
        
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.0005, 
                                        num_frame=2000, 
                                        num_rep=50, 
                                        att_level=0, 
                                        sample_move = False)

        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=1000, 
                                        num_rep=50, 
                                        att_level=0, 
                                        sample_move = False)
        # same total frames with diff frame rates
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=4000, 
                                        num_rep=50, 
                                        att_level=0, 
                                        sample_move = False)

        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=4000, 
                                        num_rep=50, 
                                        att_level=0, 
                                        sample_move = False)

        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.0005, 
                                        num_frame=4000, 
                                        num_rep=50, 
                                        att_level=0, 
                                        sample_move = False)
        # att with same num frames
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=4000, 
                                        num_rep=100, 
                                        att_level=10, 
                                        sample_move = False)

        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.0005, 
                                        num_frame=4000, 
                                        num_rep=100, 
                                        att_level=10, 
                                        sample_move = False)
        
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=4000, 
                                        num_rep=100, 
                                        att_level=10, 
                                        sample_move = False)
        # att = 15
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=4000, 
                                        num_rep=100, 
                                        att_level=15, 
                                        sample_move = False)
        
#measurement was stopped at this point at 10:36 am on 11/20/'24

        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=4000, 
                                        num_rep=100, 
                                        att_level=15, 
                                        sample_move = False)
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.0005, 
                                        num_frame=4000, 
                                        num_rep=100, 
                                        att_level=15, 
                                        sample_move = False)
        # att = 20
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=4000, 
                                        num_rep=150, 
                                        att_level=20, 
                                        sample_move = False)
        # att with same total data time
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=1000, 
                                        num_rep=150, 
                                        att_level=10, 
                                        sample_move = False)
        
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.0005, 
                                        num_frame=2000, 
                                        num_rep=150, 
                                        att_level=10, 
                                        sample_move = False)
        
        yield from bps.sleep(5)
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=4000, 
                                        num_rep=150, 
                                        att_level=10, 
                                        sample_move = False)


        