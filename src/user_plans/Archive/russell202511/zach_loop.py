
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from legacy.id8_i.plans.master_plan import run_measurement_info
from id8_common.devices import *


def round_robin(num_loop: int = 10):

    for _ in range(num_loop):
        yield from run_measurement_info("measurment_info_Xuefei1-26.json")
        yield from run_measurement_info("measurment_info_Xuefei_2-23.json")
        yield from run_measurement_info("measurment_info_Preetika_20.json")
        yield from run_measurement_info("measurment_info_Preetika_17.json")
        yield from run_measurement_info("measurment_info_Zach_PDK_acid.json")

