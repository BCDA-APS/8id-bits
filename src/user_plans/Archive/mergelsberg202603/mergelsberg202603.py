from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
import time


temp_list = [4, 8, 12, 16, 20, 24, 28, 32]

def long_scan():

    for t in temp_list: 
        print('temperature = ', t)
        set_qnw(qnw_number=2, setpoint=t, wsetait=False, ramprate=2)
        set_qnw(qnw_number=3, setpoint=t, wait=True, ramprate=2)
        time.sleep(5*60)
        run_measurement_info()

    