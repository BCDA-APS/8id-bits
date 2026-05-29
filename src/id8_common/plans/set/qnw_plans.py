"""
QNW temperature control plans for the 8ID-I beamline.

This module provides plans for controlling QNW temperature controllers, including
setting temperatures, ramping rates, and waiting for temperature stabilization.

Example:
    set_qnw(1, 20)

To watch the EPICS PVs in a simple GUI:
    pvview \
        8idi:QNWenv_1:SH_RBV \
        8idi:QNWenv_1:TARG \
        8idi:QNWenv_2:SH_RBV \
        8idi:QNWenv_2:TARG \
        8idi:QNWenv_3:SH_RBV \
        8idi:QNWenv_3:TARG \
        &
"""

import time
from apsbits.core.instrument_init import oregistry

## For air QNW
qnw_env1 = oregistry["qnw_env1"]
qnw_env2 = oregistry["qnw_env2"]
qnw_env3 = oregistry["qnw_env3"]

# For vacuum QNW
# qnw_env1 = oregistry["qnw_vac1"]
# qnw_env2 = oregistry["qnw_vac2"]
# qnw_env3 = oregistry["qnw_vac3"]

qnw_controllers = [qnw_env1, qnw_env2, qnw_env3]


def set_qnw(qnw_number: int, setpoint: float, wait: bool = True, ramprate: float = 0.3):
    """
    Change temperature on a QNW controller using Ophyd commands.
    """
    # if qnw_number < 1 or qnw_number > len(qnw_controllers):
    #     raise ValueError(f"qnw_number must be between 1 .. {len(qnw_controllers)}," f" received {qnw_number}.")
    qnw = qnw_controllers[qnw_number - 1]
    
    if qnw.ramprate.get() != ramprate:
        qnw.ramprate.put(ramprate)
    time.sleep(1)

    if wait:
        qnw.move(setpoint)
    else:
        qnw.setpoint.put(setpoint)

