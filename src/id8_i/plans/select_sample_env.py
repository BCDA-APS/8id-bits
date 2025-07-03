"""
Plan that allows for moving to pre-programmed positions seen as strings.

This module provides plans for moving the granite stage and sample stage to
predefined positions for different sample environments (QNW, rheometer, robot).
"""

from typing import Literal

import bluesky.plan_stubs as bps
from apsbits.core.instrument_init import oregistry

sample = oregistry["sample"]
granite = oregistry["granite"]
granite_8idi_valve = oregistry["granite_8idi_valve"]


def select_sample_env(env: Literal["qnw", "rheometer", "robot"]):
    """Move to a predefined sample environment position.

    This plan moves the granite stage to a predefined position for the selected
    sample environment. For the robot environment, it also moves the sample stage.
    The granite valve is enabled during motion and disabled afterward.

    Args:
        env: Sample environment to select ("qnw", "rheometer", or "robot")

    Raises:
        KeyError: If an unknown environment is specified

    Yields:
        Generator: Bluesky plan messages
    """
    choices = {
        "qnw": 923.0,
        "rheometer": 65,
        "robot": 62,
    }
    target = choices.get(env)
    if target is None:
        raise KeyError(f"Unknown environment {env=!r}")

    yield from bps.mv(granite_8idi_valve.enable, 1)
    yield from bps.sleep(2)

    if env == "qnw":
        yield from bps.mv(sample.x, 0.5)
        yield from bps.mv(granite.x, choices["qnw"])
    if env == "rheometer":
        yield from bps.mv(sample.x, 0.5)
        yield from bps.mv(granite.x, choices["rheometer"])
    elif env == "robot":
        yield from bps.mv(granite.x, choices["robot"])
        yield from bps.mv(sample.x, 298)

    yield from bps.mv(granite_8idi_valve.enable, 0)
    yield from bps.sleep(2)
