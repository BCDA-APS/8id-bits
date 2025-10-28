from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

from .shutter_logic import blockbeam
from .shutter_logic import showbeam

pv_registers = oregistry["pv_registers"]


def break_pv():
    ii = 0
    while True:
        yield from showbeam()
        yield from bps.sleep(2)
        yield from blockbeam()
        yield from bps.sleep(2)
        print(f"Cycling shutter Labjack PV for {ii} times")
        ii += 1
