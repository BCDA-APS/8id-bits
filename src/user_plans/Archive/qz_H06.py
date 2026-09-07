
import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time

def slow_ramp(sample_list=[], t_set=50, num_rep=3000):

    for sample_no in sample_list:

        print('switching sample position...')
        yield from select_sample(sample_no)
        time.sleep(2)
        print(f'switched to position {sample_no}')

        yield from temp_ramp(1)
        yield from te(6)
        yield from rigaku_acq_ZDT_series(num_rep=20, sample_move=True)

        yield from temp_ramp(1)
        yield from te(20, wait=True)

        yield from temp_ramp(0.1)
        yield from te(t_set)
        yield from rigaku_acq_ZDT_series(num_rep=num_rep, sample_move=True)


        yield from temp_ramp(1)
        yield from te(6)


def heat_cycle_rapid(sample_list=[], t_set=50, num_rep=460, cycle_num=1):

    for ii in range(cycle_num):
        for sample_no in sample_list:

            print('switching sample position...')
            yield from select_sample(sample_no)
            time.sleep(2)
            print(f'switched to position {sample_no}')

            yield from temp_ramp(1)
            yield from te(6)
            yield from rigaku_acq_ZDT_series(num_rep=20, sample_move=True)

            yield from temp_ramp(1)
            yield from te(t_set)
            yield from rigaku_acq_ZDT_series(num_rep=num_rep, sample_move=True)

            yield from temp_ramp(10)
            yield from te(6)


def heat_cycle_rapid_attcheck(sample_list=[], att_list=[], t_set=50, num_rep=460):

    yield from bps.sleep(300)

    for att_value in att_list:
        for sample_no in sample_list:

            print('switching sample position...')
            yield from select_sample(sample_no)
            time.sleep(2)
            print(f'switched to position {sample_no}')

            yield from att(att_value)

            yield from temp_ramp(1)
            yield from te(6)
            yield from rigaku_acq_ZDT_series(num_rep=20, sample_move=True)

            yield from temp_ramp(1)
            yield from te(t_set)
            yield from rigaku_acq_ZDT_series(num_rep=num_rep, sample_move=True)

            yield from temp_ramp(1)
            yield from te(6)


def bg_measure(sample_list=[], num_rep=20):

    for sample_no in sample_list:

        print('switching sample position...')
        yield from select_sample(sample_no)
        time.sleep(2)
        print(f'switched to position {sample_no}')

        yield from rigaku_acq_ZDT_series(num_rep=num_rep, sample_move=True)
        

def isothermal_gel(sample_list=None, t_set_list=None, num_rep_list=None):

    num_condition = len(sample_list)
    for ii in range(num_condition):
        sample_no = sample_list[ii]
        t_set = t_set_list[ii]
        num_rep = num_rep_list[ii]

        print('switching sample position...')
        yield from select_sample(sample_no)
        time.sleep(2)
        print(f'switched to position {sample_no}')

        yield from temp_ramp(5)
        yield from te(6, wait=True)
        yield from rigaku_acq_ZDT_series(num_rep=20, sample_move=True)

        yield from temp_ramp(10)
        yield from te(t_set-1.0, wait=True)
        yield from temp_ramp(1)
        yield from te(t_set, wait=True)
        yield from rigaku_acq_ZDT_series(num_rep=num_rep, sample_move=True)

        yield from temp_ramp(5)
        yield from te(6)


def last_overnight():

    yield from heat_cycle_rapid(sample_list=[2], t_set=50, num_rep=460, cycle_num=1)

    yield from isothermal_gel(sample_list=[8], t_set_list=[26], num_rep_list=[1500])

    yield from heat_cycle_rapid(sample_list=[2], t_set=50, num_rep=460, cycle_num=1)

    yield from isothermal_gel(sample_list=[5], t_set_list=[34], num_rep_list=[1000])

    yield from heat_cycle_rapid(sample_list=[2], t_set=50, num_rep=460, cycle_num=1)

      





