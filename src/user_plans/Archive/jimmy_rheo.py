
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np    

# edited so that num_rep refers to the number of rheology runs, not the 
# number of runs one after the other (ie, there are wait steps between each)
def rheo_xpcs_wait(acq_period, num_frame, num_rep, att_level, wait_steps):
    
    for jj in range(num_rep):

        if jj==0:
            current_wait_steps = wait_steps
        else:
            current_wait_steps = wait_steps+1

        yield from mesh_grid_move(sam_index=0,
                                x_cen=-3.53,
                                x_radius=0.050,
                                x_pts=5,
                                y_cen=359.3,
                                y_radius=1,
                                y_pts=5)
        yield from post_align()
        for ii in range(current_wait_steps):
            yield from wait_for_mcr()

        yield from eiger_acq_int_series(acq_period=acq_period, 
                                        num_frame=num_frame, 
                                        num_rep=1, 
                                        att_level=att_level, 
                                        sample_move = False)
