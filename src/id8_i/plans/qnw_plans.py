"""
QNW temperature control plans for the 8ID-I beamline.

This module provides plans for controlling QNW temperature controllers, including
setting temperatures, ramping rates, and waiting for temperature stabilization.

Example:
    RE(set_qnw(1, 20))

To watch the EPICS PVs in a simple GUI:
    pvview \\
        8idi:QNWenv_1:SH_RBV \\
        8idi:QNWenv_1:TARG \\
        8idi:QNWenv_2:SH_RBV \\
        8idi:QNWenv_2:TARG \\
        8idi:QNWenv_3:SH_RBV \\
        8idi:QNWenv_3:TARG \\
        &
"""

from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

from .sample_info_unpack import sort_qnw

qnw_env1 = oregistry["qnw_env1"]
qnw_env2 = oregistry["qnw_env2"]
qnw_env3 = oregistry["qnw_env3"]


def set_qnw(qnw_number: int, setpoint: float, wait: bool = True, ramprate: float = 0.3):
    """
    Change temperature on a QNW controller from a bluesky plan.
    """
    if qnw_number < 1 or qnw_number > len(qnw_controllers):
        raise ValueError(f"qnw_number must be between 1 .. {len(qnw_controllers)}," f" received {qnw_number}.")
    qnw = qnw_controllers[qnw_number - 1]
    if qnw.ramprate.get() != ramprate:
        yield from bps.mv(qnw.ramprate, ramprate)
    yield from bps.sleep(1)

    if wait:
        yield from bps.mv(qnw, setpoint)
    else:
        yield from bps.mv(qnw.setpoint, setpoint)


# Reference these controllers in a list by an index number with 1 offset.

# For air qnw
qnw_controllers = [qnw_env1, qnw_env2, qnw_env3]

# For vac qnw
# qnw_controllers = [qnw_vac1, qnw_vac2, qnw_vac3]


def find_qnw_index() -> int:
    """Find the index of the currently selected QNW environment.

    Returns:
        int: QNW controller index (1-3)

    Raises:
        ValueError: If no QNW environment is selected
    """
    sam_dict = sort_qnw()
    if sam_dict["temp_zone"] == "qnw_env1":
        qnw_number = 1
    elif sam_dict["temp_zone"] == "qnw_env2":
        qnw_number = 2
    elif sam_dict["temp_zone"] == "qnw_env3":
        qnw_number = 3
    else:
        raise ValueError("No QNW environment selected")
    return qnw_number


def set_temperature(setpoint: float, wait: bool = False):
    """Set the temperature of the selected QNW controller.

    Args:
        setpoint: Target temperature in degrees Celsius
        wait: Whether to wait for temperature to stabilize
    """
    qnw_number = find_qnw_index()
    qnw = qnw_controllers[qnw_number - 1]
    if wait:
        yield from bps.mv(qnw, setpoint)
    else:
        yield from bps.mv(qnw.setpoint, setpoint)
    yield from bps.sleep(2)


def set_ramp_rate(ramprate: float = 0.3):
    """Set the temperature ramp rate for the selected QNW controller.

    Args:
        ramprate: Temperature ramp rate in degrees Celsius per minute
    """
    qnw_number = find_qnw_index()
    qnw = qnw_controllers[qnw_number - 1]
    yield from bps.mv(qnw.ramprate, ramprate)
    yield from bps.sleep(2)
    yield from bps.mv(qnw.ramprate, ramprate)
    yield from bps.sleep(2)


def set_temperature_with_ramp(
    setpoint: float,
    ramprate: float,
    wait: bool,
):
    """Set temperature and ramp rate for the selected QNW controller.

    Args:
        setpoint: Target temperature in degrees Celsius
        ramprate: Temperature ramp rate in degrees Celsius per minute
        wait: Whether to wait for temperature to stabilize
    """
    qnw_number = find_qnw_index()
    qnw = qnw_controllers[qnw_number - 1]
    # Run this thrice just to make sure all settings go through
    yield from bps.mv(qnw.ramprate, ramprate)
    yield from bps.sleep(2)
    if wait:
        yield from bps.mv(qnw, setpoint)
    else:
        yield from bps.mv(qnw.setpoint, setpoint)
    yield from bps.sleep(2)

    yield from bps.mv(qnw.ramprate, ramprate)
    yield from bps.sleep(2)
    if wait:
        yield from bps.mv(qnw, setpoint)
    else:
        yield from bps.mv(qnw.setpoint, setpoint)
    yield from bps.sleep(2)

    yield from bps.mv(qnw.ramprate, ramprate)
    yield from bps.sleep(2)
    if wait:
        yield from bps.mv(qnw, setpoint)
    else:
        yield from bps.mv(qnw.setpoint, setpoint)
    yield from bps.sleep(2)


def set_temperature_env(
    qnw_number: int,
    setpoint: float,
    wait: bool = False,
):
    """Set temperature for a specific QNW controller.

    Args:
        qnw_number: QNW controller index (1-3)
        setpoint: Target temperature in degrees Celsius
        wait: Whether to wait for temperature to stabilize
    """
    qnw = qnw_controllers[qnw_number - 1]
    if wait:
        yield from bps.mv(qnw, setpoint)
    else:
        yield from bps.mv(qnw.setpoint, setpoint)


def set_ramp_rate_env(
    qnw_number: int,
    ramprate: float = 0.3,
):
    """Set ramp rate for a specific QNW controller.

    Args:
        qnw_number: QNW controller index (1-3)
        ramprate: Temperature ramp rate in degrees Celsius per minute
    """
    qnw = qnw_controllers[qnw_number - 1]
    yield from bps.mv(qnw.ramprate, ramprate)
