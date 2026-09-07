import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time
import epics as pe

func_pv="dpKeysight:KEY1:1:FUNC"  #0 for sin, 4 for square wave, 2 for triangle wave
freq_pv="dpKeysight:KEY1:1:FREQ"
ncycle_pv="dpKeysight:KEY1:1:BURST:NCYCLES"
amp_pv="dpKeysight:KEY1:1:AMP"
trigger_pv="dpKeysight:KEY1:TRIG.PROC"
pulse_width_pv="dpKeysight:KEY1:1:PULSEW"
offset_pv="dpKeysight:KEY1:1:OFFSET"

# pe.caput("dpKeysight:KEY1:1:FUNC",func)
# epics_put("dpKeysight:KEY1:1:BURST:NCYCLES",ncycle) ##number of waves
# epics_put("dpKeysight:KEY1:1:FREQ", freq) ##freq in Hz
# epics_put("dpKeysight:KEY1:1:AMP", amp) #peak to peak amplitude in Volts
#yield from bps.sleep(3)
# epics_put("dpKeysight:KEY1:TRIG.PROC",send_trigger) #send the waves

def pmnpt_anudeep_att_scan():
    yield from bps.mv(pv_registers.sample_name, 'PMNPT')
    atts=np.array([0,2,5,14,25,40])
    for att in atts:
        yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=20000, 
                                        num_rep=1, 
                                        att_level=att, 
                                        sample_move = False)
        
        yield from bps.sleep(5)

def pmnpt_anudeep_sin_E_field_test(freq,amp,total_time,nrep=1):
    att=14
    func=0
    duty_cycle=0.0
    offset=0.0
    ncycle=np.ceil((total_time+30)*freq*(nrep))
    send_trigger=1
    pe.caput(func_pv,func)
    pe.caput(freq_pv,freq)
    pe.caput(amp_pv,amp)
    pe.caput(ncycle_pv,ncycle)
    pe.caput(pulse_width_pv,duty_cycle/freq)
    pe.caput(offset_pv,offset)
    yield from bps.sleep(3)
    pe.caput(trigger_pv,send_trigger)
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=int(np.floor(total_time/0.001)), 
                                        num_rep=nrep, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.sleep(30)
    
def pmnpt_anudeep_square_E_field_test(freq,amp,total_time,nrep=1):
    att=14
    func=4
    duty_cycle=0.5
    offset=0.0
    ncycle=np.ceil((total_time+30)*freq*(nrep))
    send_trigger=1
    pe.caput(func_pv,func)
    pe.caput(freq_pv,freq)
    pe.caput(amp_pv,amp)
    pe.caput(pulse_width_pv,duty_cycle/freq)
    pe.caput(offset_pv,offset)
    pe.caput(ncycle_pv,ncycle)
    yield from bps.sleep(3)
    pe.caput(trigger_pv,send_trigger)
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.001, 
                                        num_frame=int(np.floor(total_time/0.001)), 
                                        num_rep=nrep, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.sleep(30)
    
def pmnpt_anudeep_freq_sweep():
    freqs=[1,2,5,10,100,200,500,1000,2000,5000,10000,100000]
    # freqs=np.array([1,2])
    amps=[0.5,1,1.5]
    # amps=np.array([0.5,1])
    reps=3
    total_time=100
    for amp_ind in range(len(amps)):
        for freq_ind in range(len(freqs)):
            yield from bps.mv(pv_registers.sample_name, f'PMNPT_sin-freq-{freqs[freq_ind]}Hz_Vpp-{amps[amp_ind]}V')
            yield from pmnpt_anudeep_sin_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,nrep=reps)
            yield from bps.mv(pv_registers.sample_name, f'PMNPT_square-freq-{freqs[freq_ind]}Hz_Vpp-{amps[amp_ind]}V')
            yield from pmnpt_anudeep_square_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,nrep=reps)