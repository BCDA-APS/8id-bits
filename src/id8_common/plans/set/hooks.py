import time as time 
from datetime import datetime
from apsbits.core.instrument_init import oregistry
import numpy as np


keithley = oregistry["keithley2400"]
lambda2M = oregistry["lambda2M"]

def volt_cycle_single(voltage_file = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_single.txt')):

    voltage_file = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_single.txt')

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
    return not (
        lambda2M.hdf1.capture.get() == 1
        and lambda2M.hdf1.num_captured.get() < lambda2M.hdf1.num_capture.get()
    )

def volt_cycle_series(voltage_program = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_series.txt')):
    
    voltage_program = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_series.txt')
    
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

            



    


    