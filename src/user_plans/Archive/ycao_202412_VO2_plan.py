import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time
import epics as pe

func_pv="dpKeysight:KEY1:1:FUNC"  #0 for sin, 4 for square wave, 2 for triangle wave, 8 for DC wave
freq_pv="dpKeysight:KEY1:1:FREQ"
ncycle_pv="dpKeysight:KEY1:1:BURST:NCYCLES"
amp_pv="dpKeysight:KEY1:1:AMP"
trigger_pv="dpKeysight:KEY1:TRIG.PROC"
pulse_width_pv="dpKeysight:KEY1:1:PULSEW"
offset_pv="dpKeysight:KEY1:1:OFFSET"

def VO2_sin_E_field_test(freq,amp,total_time,frame_rate,nrep=1):
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

    num_frame = int(total_time / frame_rate)
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=frame_rate, 
                                        num_frame=num_frame, 
                                        num_rep=nrep, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.sleep(30)
    
    

def VO2_trang_E_field_test(freq,amp,total_time,frame_rate, nrep=1):
    att=14
    func=2
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

    num_frame = int(total_time / frame_rate)
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=frame_rate, 
                                        num_frame=num_frame, 
                                        num_rep=nrep, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.sleep(3)
    

def VO2_square_E_field_offset_test(freq,amp,total_time,frame_rate,nrep=1):
    att=14
    func=4
    duty_cycle=0.5
    offset=0.5*amp
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
                                        acq_period=frame_rate, 
                                        num_frame=5000, 
                                        num_rep=nrep, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.sleep(3)
    

def VO2_dc_mode_field_test(freq,amp,total_time,frame_rate, nrep=1):
    att=14
    func=8
    # duty_cycle=0.5
    offset=0.0
    ncycle=np.ceil((total_time+30)*freq*(nrep))
    send_trigger=1
    pe.caput(func_pv,func)
    # pe.caput(freq_pv,freq)
    pe.caput(amp_pv,amp)
    # pe.caput(pulse_width_pv,duty_cycle/freq)
    # pe.caput(offset_pv,offset)
    # pe.caput(ncycle_pv,ncycle)
    yield from bps.sleep(3)
    pe.caput(trigger_pv,send_trigger)

    num_frame = int(total_time / frame_rate)
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=frame_rate, 
                                        num_frame=num_frame, 
                                        num_rep=nrep, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.sleep(3)


def VO2_volt_freq_sweep():
    freqs=[1,5,10,100,500,1000,2000,5000,10000]
    frame_rate_list = [0.01]
    # amps=[1,3,5,7,8,9,10]
    amps=[12,16,18]
    reps=2
    num_frame = 5000
    
    for frame_rate in frame_rate_list:
        total_time= num_frame * frame_rate
        for amp_ind in range(len(amps)):
            for freq_ind in range(len(freqs)):
                # yield from bps.mv(pv_registers.sample_name, f'VO2_1207_sin-freq-{freqs[freq_ind]}Hz_frr-{frame_rate}_Vpp-{amps[amp_ind]}V')
                # yield from VO2_sin_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,frame_rate,nrep=reps)
                yield from bps.mv(pv_registers.sample_name, f'VO2_1207_trang-freq-{freqs[freq_ind]}Hz_frr-{frame_rate}_Vpp-{amps[amp_ind]}V')
                yield from VO2_trang_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,frame_rate,nrep=reps)
                # yield from bps.mv(pv_registers.sample_name, f'VO2_1207_square-freq-{freqs[freq_ind]}Hz_Vpp-{amps[amp_ind]}V-ofs')
                # yield from VO2_square_E_field_offset_test(freqs[freq_ind],amps[amp_ind],total_time,nrep=reps)

 
 
 
def VO2_volt_offset_sweep():
    freqs=[1,10,100,1000,10000]
    frame_rate_list = [0.01]
    amps=[0.1,1,3,7,9] 
    reps=2
    num_frame = 5000
    
    for frame_rate in frame_rate_list:
        total_time= num_frame * frame_rate
        for amp_ind in range(len(amps)):
            for freq_ind in range(len(freqs)):

                # yield from bps.mv(pv_registers.sample_name, f'VO2_1207_sin-freq-{freqs[freq_ind]}Hz_frr-{frame_rate}_Vpp-{amps[amp_ind]}V')
                # yield from VO2_sin_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,frame_rate,nrep=reps)
                # yield from bps.mv(pv_registers.sample_name, f'VO2_1207_trang-freq-{freqs[freq_ind]}Hz_frr-{frame_rate}_Vpp-{amps[amp_ind]}V')
                # yield from VO2_trang_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,frame_rate,nrep=reps)
                yield from bps.mv(pv_registers.sample_name, f'VO2_1207_square-freq-{freqs[freq_ind]}Hz_frr-{frame_rate}_Vpp-{amps[amp_ind]}V-ofs')
                yield from VO2_square_E_field_offset_test(freqs[freq_ind],amps[amp_ind],total_time,frame_rate,nrep=reps)               


def VO2_volt_freqq_sweep():
    freqs=[1,5,10,100,500,1000,2000,5000,10000]
    frame_rate_list = [0.00025]
    amps=[1,3,5,7,8,9,10] 
    reps=2
    num_frame = 12000
    
    for frame_rate in frame_rate_list:
        total_time= num_frame * frame_rate
        for amp_ind in range(len(amps)):
            for freq_ind in range(len(freqs)):
                yield from bps.mv(pv_registers.sample_name, f'VO2_1207_sin-freq-{freqs[freq_ind]}Hz_frr-{frame_rate}_Vpp-{amps[amp_ind]}V')
                yield from VO2_sin_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,frame_rate,nrep=reps)
                yield from bps.mv(pv_registers.sample_name, f'VO2_1207_trang-freq-{freqs[freq_ind]}Hz_frr-{frame_rate}_Vpp-{amps[amp_ind]}V')
                yield from VO2_trang_E_field_test(freqs[freq_ind],amps[amp_ind],total_time,frame_rate,nrep=reps)
                yield from bps.mv(pv_registers.sample_name, f'VO2_1207_square-freq-{freqs[freq_ind]}Hz_Vpp-{amps[amp_ind]}V-ofs')
                yield from VO2_square_E_field_offset_test(freqs[freq_ind],amps[amp_ind],total_time,nrep=reps)


def VO2_volt_dcmode_sweep():
    frame_rate=0.1   # each frame period
    amps=[3,5,6,7,8,9] 
    reps=1
    num_frame = 12000
    total_time = num_frame * frame_rate

    for amp_ind in range(len(amps)):
        yield from bps.mv(pv_registers.sample_name, f'VO2_1207_dc-frr-{frame_rate}_Vpp-{amps[amp_ind]}V')
        yield from VO2_dc_mode_field_test(1,amps[amp_ind],total_time,frame_rate,nrep=reps)
                

def dc_apply_test():
    func=4
    freq=1
    amp=1
    duty_cycle=0.5
    offset=1
    ncycle=10
    send_trigger=1
    pe.caput(func_pv,func)
    pe.caput(freq_pv,freq)
    pe.caput(amp_pv,amp)
    pe.caput(pulse_width_pv,duty_cycle/freq)
    pe.caput(offset_pv,offset)
    pe.caput(ncycle_pv,ncycle)
    yield from bps.sleep(3)
    pe.caput(trigger_pv,send_trigger)


