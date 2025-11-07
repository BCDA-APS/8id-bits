"""
Rheometer wait plans for the 8ID-I beamline.

This module provides plans for waiting on the MCR Rheometer to change states
during measurements. It monitors a pulse count signal to detect state changes.
"""

import datetime

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

mcr_wait_signal = oregistry["mcr_wait_signal"]


def wait_for_mcr(delay_time: float = 0.01):
    """Wait for the MCR Rheometer to change state.

    This plan monitors the pulse count signal from the MCR Rheometer and waits
    until it changes from its current value, indicating a state change. After
    detecting the change, it waits for an additional delay time.

    Args:
        delay_time: Additional time to wait after state change (in seconds)

    Yields:
        Generator: Bluesky plan messages
    """
    print("Waiting for the MCR Rheometer to change FROM THE OLD STATE")
    print(datetime.datetime.now())

    current_value = mcr_wait_signal.pulse_count.get()

    while mcr_wait_signal.pulse_count.get() == current_value:
        yield from bps.sleep(0.1)

    print("MCR Rheometer changed to the NEW MEASUREMENT STATE")
    print(datetime.datetime.now())

    yield from bps.sleep(delay_time)
