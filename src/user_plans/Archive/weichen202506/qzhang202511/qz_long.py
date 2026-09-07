
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from legacy.id8_i.plans.master_plan import run_measurement_info
from id8_common.devices import *

def qz_long():
    yield from bps.sleep(600)
    yield from run_measurement_info()
    