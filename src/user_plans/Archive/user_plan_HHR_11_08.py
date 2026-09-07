
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np



def framerate_scan():

    # yield from select_sample(3)

    yield from showbeam()

    sample_name = 'Cys1_40'
    shear_rate = 'steady'
    qmap_file = '/gdata/dm/8IDI/2024-3/comm202411/data/eiger4m_qmap_1107_s360_d36_phi1_full_HHR.h5'
    acq_time = 0.00025
    att = 10

    yield from select_sample(3)
    yield from eiger_acq_int_series(acq_period=0.00025, num_frame=1000, num_rep=1, att_level=att, sample_move=False)
    yield from eiger_acq_int_series(acq_period=0.001, num_frame=1000, num_rep=1, att_level=att, sample_move=False)
    yield from eiger_acq_int_series(acq_period=0.01, num_frame=1000, num_rep=1,  att_level=att,sample_move=False)
    yield from eiger_acq_int_series(acq_period=0.1, num_frame=1000, num_rep=1,  att_level=att,sample_move=False)
    yield from eiger_acq_ext_trig(acq_time=0.1, acq_period=1, num_frame=1000, num_rep=1,  att_level=att,sample_move=False)
    yield from blockbeam()


def framerate_scan():

    # yield from select_sample(3)

    yield from showbeam()

    sample_name = 'Cys0_48'
    shear_rate = 'shear0p001'
    qmap_file = '/gdata/dm/8IDI/2024-3/comm202411/data/eiger4m_qmap_1107_s360_d36_phi1_full_HHR.h5'
    acq_time = 0.00025
    att = 10

    yield from select_sample(3)
    yield from eiger_acq_int_series(acq_period=0.00025, num_frame=1000, num_rep=1, att_level=att, sample_move=False)
    yield from eiger_acq_int_series(acq_period=0.001, num_frame=1000, num_rep=1, att_level=att, sample_move=False)
    yield from eiger_acq_int_series(acq_period=0.01, num_frame=1000, num_rep=1,  att_level=att,sample_move=False)
    yield from eiger_acq_int_series(acq_period=0.1, num_frame=1000, num_rep=1,  att_level=att,sample_move=False)
    yield from eiger_acq_ext_trig(acq_time=0.1, acq_period=1, num_frame=1000, num_rep=1,  att_level=att,sample_move=False)
    yield from blockbeam()


def vertical_scan():

    # yield from select_sample(3)
    sample_name = 'Cys0_48'
    shear_rate = 'steady'
    qmap_file = '/gdata/dm/8IDI/2024-3/comm202411/data/eiger4m_qmap_1107_s360_d36_phi1_full_HHR.h5'
    acq_time = 0.00025

    

    vertical_position = np.linspace(352, 358, 20)
    yield from select_sample(3)
    for vpos in vertical_position:
        yield from bps.mv(rheometer.y, vpos)
        yield from showbeam()
        yield from eiger_acq_int_series(acq_period=0.00025, num_frame=100, num_rep=1, att_level=10, sample_move=False)
        yield from blockbeam()
    


def att_scan():

    # yield from select_sample(3)

    

    sample_name = 'Cys0_48'
    shear_rate = 'steady'
    qmap_file = '/gdata/dm/8IDI/2024-3/comm202411/data/eiger4m_qmap_1107_s360_d36_phi1_full_HHR.h5'
    acq_time = 0.00025

    yield from select_sample(3)
    att_lvl_list = [25,24,23,22,21,20,19,18,17,16,15]

    for att_lvl in att_lvl_list:
        yield from showbeam()
        yield from eiger_acq_int_series(acq_period=0.1, num_frame=1000, num_rep=1,  att_level=att_lvl,sample_move=False)
        yield from blockbeam()