
import bluesky.plan_stubs as bps

from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *

import numpy as np


def bcdi_rock(th_cen = None, 
              th_radius = None, 
              th_pts = None, 
              acq_time_total = 1.0,
              att_level = None):

    th_list = np.linspace(th_cen-th_radius, th_cen+th_radius, th_pts)

    yield from bps.mv(filter_8ide.attenuation_set, att_level)
    yield from bps.sleep(2)
    yield from bps.mv(filter_8ide.attenuation_set, att_level)
    yield from bps.sleep(2)


    yield from shutteroff()
    yield from showbeam()
    yield from bps.sleep(0.1)

    for ii in range(th_pts):
        num_frame = 10
        yield from bps.mv(huber.mu, th_list[ii])
        yield from eiger_acq_int_rock(acq_period = acq_time_total/num_frame, 
                                      num_frame = num_frame, 
                                      att_level = att_level, 
                                      angle = th_list[ii]) 
        
    yield from blockbeam()

    meas_num = int(pv_registers.measurement_num.get())
    yield from bps.mv(pv_registers.measurement_num, meas_num + 1)

