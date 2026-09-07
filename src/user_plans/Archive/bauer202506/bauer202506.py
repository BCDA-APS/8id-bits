
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from legacy.id8_i.plans.master_plan import run_measurement_info
from legacy.id8_i.devices import *

robotic_pipette = oregistry["robotic_pipette"]

def pipet_eiger():

    robotic_pipette.step_1.put(1)
    yield from bps.sleep(40)
    #yield from run_measurement_info('measurement_info_dry.json') # 1 s
    yield from run_measurement_info('measurement_info_dry-rigaku.json') # 1 s
    
    
    
    
    robotic_pipette.step_2.put(1)
    #yield from bps.sleep(3)
    #yield from run_measurement_info('measurement_info_wet.json') # 30 s, 2 s
    yield from run_measurement_info('measurement_info_wet-rigaku.json') # 30 s, 2 s


