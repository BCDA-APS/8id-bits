
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from legacy.id8_i.plans.master_plan import run_measurement_info
from id8_common.devices import *


def round_robin(num_loop: int = 10):

    for _ in range(num_loop):
        yield from run_measurement_info()

    #for sample in samples:
    #    # move to the next sample we want to measure
    #    select_sample(sample)

    print("your scan goes here")
