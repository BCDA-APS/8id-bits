import bluesky.plan_stubs as bps
from aps_8id_bs_instrument.plans import *
from aps_8id_bs_instrument.devices import *
import numpy as np
import time
import epics as pe

def TaS2_take_XPCS_cdw(voltage,att):
    yield from bps.sleep(30)
    yield from bps.mv(pv_registers.sample_name, f'TaS2_Evans_1_1245_0068+0667-{voltage}mV_4kHz-acq1')
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=4000, 
                                        num_rep=10, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.mv(pv_registers.sample_name, f'TaS2_Evans_1_1245_0068+0667-{voltage}mV_125Hz-acq')
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.008, 
                                        num_frame=1250, 
                                        num_rep=20, 
                                        att_level=att, 
                                        sample_move = False)
    yield from bps.mv(pv_registers.sample_name, f'TaS2_Evans_1_1245_0068+0667-{voltage}mV_4kHz-acq2')
    yield from eiger_acq_int_series(eiger4M, 
                                        acq_period=0.00025, 
                                        num_frame=4000, 
                                        num_rep=10, 
                                        att_level=att, 
                                        sample_move = False)
    #ASK YUE ABOUT INCLUDING THE FOLLOWING COMMAND!!!!!
    # yield from bps.mv(pv_registers.sample_name, f'TaS2_Evans_1_1245_0068+0667-{voltage}V_40Hz-acq')
    # yield from eiger_acq_int_series(eiger4M, 
    #                                     acq_period=0.025, 
    #                                     num_frame=4000, 
    #                                     num_rep=1, 
    #                                     att_level=att, 
    #                                     sample_move = False)   

    #DID YOU ASK???
        
    yield from bps.sleep(5)