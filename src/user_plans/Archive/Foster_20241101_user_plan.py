import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time

def temp_ramp():
    yield from set_qnw_vac(qnw_number=2, setpoint = 150, wait = False, ramprate = 10)
    yield from set_qnw_vac(qnw_number=2, setpoint = 160, wait = False, ramprate = 1)


def dummy():
    yield from set_qnw_vac(qnw_number=2, setpoint = 150, wait = True, ramprate = 10)
    yield from set_qnw_vac(qnw_number=2, setpoint = 160, wait = False, ramprate = 1)
    yield from bps.sleep(5)
    yield from eiger_acq_int_series(eiger4M, acq_period=1, num_frame=2000, num_rep=2, att_level=15, sample_move = True)

    yield from select_sample(5)

def snap():
    yield from eiger_acq_ext_trig(eiger4M, acq_period=.05, num_frame=1, num_rep=1, att_level=10, sample_move = False)
