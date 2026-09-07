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
from id8_common.registry import get_connected_device

# Which QNW set these plans drive. Swap to the air controllers by changing
# these three names to "qnw_env1"/"qnw_env2"/"qnw_env3".
#
# These used to be bound at module scope as qnw_env1..3 -- aliases pointing at
# the *vacuum* devices. Because startup star-imports this module, those aliases
# also landed in the interactive namespace and shadowed the real air devices of
# the same name that the device loader had put there. Resolving by name per
# call removes the shadowing: at the prompt, qnw_env1 is now the air device it
# is named after, and qnw_vac1 is the vacuum one. set_qnw() still drives the
# vacuum controllers exactly as before.
QNW_NAMES = ("qnw_vac1", "qnw_vac2", "qnw_vac3")


def _qnw(qnw_number: int):
    """Resolve QNW controller 1..3, raising if it is absent or disconnected."""
    if not 1 <= qnw_number <= len(QNW_NAMES):
        raise ValueError(f"qnw_number must be between 1 .. {len(QNW_NAMES)}, received {qnw_number}.")

    return get_connected_device(QNW_NAMES[qnw_number - 1])


def set_qnw(qnw_number: int, setpoint: float, wait: bool = True, ramprate: float = 0.3):
    """
    Change temperature on a QNW controller using Ophyd commands.

    Args:
        qnw_number: Which controller, 1..3 (see QNW_NAMES).
        setpoint: Target temperature, degrees C.
        wait: True blocks until the controller reports it has reached the
            setpoint; False writes the setpoint and returns immediately, so
            the sample keeps ramping while the caller does something else.
        ramprate: How fast the controller is allowed to approach the setpoint,
            in whatever units the RAMP PV takes. Only written when it differs
            from what the controller already has.
    """
    qnw = _qnw(qnw_number)
    
    if qnw.ramprate.get() != ramprate:
        qnw.ramprate.put(ramprate)
    time.sleep(1)

    if wait:
        qnw.move(setpoint)
    else:
        qnw.setpoint.put(setpoint)

