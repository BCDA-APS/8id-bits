"""
Bias-voltage sequences run alongside a lambda2M acquisition.

Both plans below are meant to be started once the detector is already
capturing: they first wait for `preset` frames to land, then cycle the
Keithley source until the HDF plugin has captured everything it was asked for.

The voltage program itself lives in a text file next to this module, one
number per line, so it can be edited between runs without touching Python.

This module is a copy of volt_seq.py in the same directory; the two differ
only in commented-out code. Neither is imported by startup.py.
"""

import time as time 
from datetime import datetime
from apsbits.core.instrument_init import oregistry
import numpy as np


keithley = oregistry["keithley2400"]
lambda2M = oregistry["lambda2M"]

def volt_cycle_single(voltage_file = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_single.txt')):

    # The voltage_file parameter is not actually used: the file is re-read here
    # on every call, so an edit to voltage_program_single.txt takes effect
    # without restarting the session. (The default value is only evaluated once,
    # at import, which is exactly why it cannot be relied on.)
    voltage_file = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_single.txt')

    # The file is five numbers, one per line, in the order unpacked below. Each
    # line carries a trailing '# label' that np.loadtxt strips as a comment.

    preset = voltage_file[0]
    bias_time = voltage_file[1]
    pulse_period = voltage_file[2]
    neg_voltage = voltage_file[3]
    pos_voltage = voltage_file[4]

    while (
        lambda2M.hdf1.capture.get() == 1
        and lambda2M.hdf1.num_captured.get() < preset
    ):
        time.sleep(0.05)

    while (
        lambda2M.hdf1.capture.get() == 1
        and lambda2M.hdf1.num_captured.get() < lambda2M.hdf1.num_capture.get()
    ):
        keithley.set_volt.put(0)
        keithley.output.put(1)
        keithley.set_volt.put(pos_voltage)
        # time.sleep(1)
        time.sleep(bias_time)
        keithley.set_volt.put(0)

        time.sleep(pulse_period)

        keithley.set_volt.put(neg_voltage)
        keithley.output.put(1)
        time.sleep(bias_time)
        keithley.set_volt.put(0)

        time.sleep(pulse_period)

def check_done():
    """True once the HDF plugin has stopped capturing, or has all its frames."""
    return not (
        lambda2M.hdf1.capture.get() == 1
        and lambda2M.hdf1.num_captured.get() < lambda2M.hdf1.num_capture.get()
    )

def volt_cycle_series(voltage_program = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_series.txt')):
    
    # Re-read on every call for the same reason as in volt_cycle_single():
    # the default above is evaluated once at import and is then ignored.
    voltage_program = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_series.txt')

    # Six numbers, one per line, in the order unpacked below -- same trailing
    # '# label' convention as voltage_program_single.txt.
    
    preset = voltage_program[0]
    bias_time = voltage_program[1]
    pulse_period = voltage_program[2]
    v_min = voltage_program[3]
    v_max = voltage_program[4]
    v_step = voltage_program[5]

    volt_list_up = np.arange(v_min, v_max, v_step)
    volt_list_down = volt_list_up*-1

    while (
            lambda2M.hdf1.capture.get() == 1
            and lambda2M.hdf1.num_captured.get() < preset
        ):
        time.sleep(0.05)    
    
    while not check_done():
        for i in np.concatenate([volt_list_up, volt_list_down]):
            if check_done():
                break
            print('voltage is ', i)
            keithley.output.put(1)
            keithley.set_volt.put(i)
            time.sleep(bias_time)
            keithley.set_volt.put(0)
            time.sleep(pulse_period)
        keithley.set_volt.put(0)

            



    


    