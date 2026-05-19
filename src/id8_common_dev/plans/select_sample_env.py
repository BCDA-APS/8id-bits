"""
Plan that allows for moving to pre-programmed positions seen as strings.

This module provides plans for moving the granite stage and sample stage to
predefined positions for different sample environments (QNW, rheometer, robot).
"""

from typing import Literal

# import bluesky.plan_stubs as bps
import time
from apsbits.core.instrument_init import oregistry

sample = oregistry["sample"]
granite = oregistry["granite"]
granite_8idi_valve = oregistry["granite_8idi_valve"]


def select_sample_env(env: Literal["qnw", "rheometer", "robot"]):
    """Move to a predefined sample environment position.

    This plan moves the granite stage to a predefined position for the selected
    sample environment. For the robot environment, it also moves the sample stage.
    The granite valve is enabled during motion and disabled afterward.

    """
    choices = {
        "qnw": 923.0,
        "rheometer": 73.7,
        "robot": 62,
    }

    ## QZ changed on 08/01 for Naomi's horizontal capillary setup
    # choices = {
    #     "qnw": 923.0,
    #     "rheometer": 88,
    #     "robot": 62,
    #     "clip": 859
    # }

    target = choices.get(env)
    if target is None:
        raise KeyError(f"Unknown environment {env=!r}")

    granite_8idi_valve.enable.put(1)
    time.sleep(2.0)

    if env == "qnw":
        granite.x.move(choices["qnw"])
    if env == "rheometer":
        granite.x.move(choices["rheometer"])
    elif env == "robot":
        granite.x.move(choices["robot"])

    granite_8idi_valve.enable.put(0)
    time.sleep(2)
