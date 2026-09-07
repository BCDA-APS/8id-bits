
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time


def brian_overnight(num_r, att):

    # yield from post_align()
    # yield from select_sample(5)


    yield from rigaku_acq_ZDT_series(
                                acq_period=2e-5, 
                                num_frame=100000, 
                                num_rep=num_r, 
                                att_level=att,
                                sample_move=False, 
                                process=True
                                )

def overnight_scan():
    # yield from post_align()
    # yield from select_sample(5)
    x_start = 149.6
    y_start = 21.05
    yield from bps.mv(sample.x, x_start + 0.125)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,17)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=400, att=50)
    yield from brian_overnight(num_r=20, att=10)
    yield from brian_overnight(num_r=1, att=2)

    yield from bps.mv(sample.x, x_start + 0.150)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,16.5)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=400, att=50)
    yield from brian_overnight(num_r=20, att=10)
    yield from brian_overnight(num_r=1, att=2)

    yield from bps.mv(sample.x, x_start + 0.175)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,16)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=300, att=50)
    yield from brian_overnight(num_r=15, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start + 0.200)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,15.5)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=300, att=50)
    yield from brian_overnight(num_r=15, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start - 0.025)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,15)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=300, att=50)
    yield from brian_overnight(num_r=15, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start - 0.05)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,14.5)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=200, att=50)
    yield from brian_overnight(num_r=10, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start - 0.075)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,14)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=200, att=50)
    yield from brian_overnight(num_r=10, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start - 0.100)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,13.5)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=200, att=50)
    yield from brian_overnight(num_r=10, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start - 0.125)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,13)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=200, att=50)
    yield from brian_overnight(num_r=10, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start - 0.150)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,12.5)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=200, att=50)
    yield from brian_overnight(num_r=10, att=10)
    yield from brian_overnight(num_r=1, att=2)
    
    yield from bps.mv(sample.x, x_start - 0.175)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,12)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=200, att=50)
    yield from brian_overnight(num_r=10, att=10)
    yield from brian_overnight(num_r=1, att=2)
    

def single_scan():
    yield from post_align()
    # yield from select_sample(5)
    x_start = 149.6
    y_start = 21.05
    yield from bps.mv(sample.x, x_start + 0.100)
    # yield from bps.mv(sample.y, y_start)
    yield from bps.mv(qnw_env2.ramprate,0.2)
    yield from bps.sleep(10)
    yield from bps.mv(qnw_env2.setpoint,18)
    yield from bps.sleep(900)
    yield from brian_overnight(num_r=500, att=50)

    