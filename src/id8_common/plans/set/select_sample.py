"""
Select and move to a sample position using the sample_info.yaml configuration.
"""

import yaml
from pathlib import Path

from id8_common.expt_config import expt
from id8_common.registry import oregistry


# There is deliberately no module-level path constant here. The YAML this
# module reads is expt.sample_info_file, resolved from configs/experiment.yml
# on every call, so editing the config takes effect without a restart. It was a
# hardcoded absolute path until 2026-09-06.
#
# Shape of sample_info.yaml::
#
#     defaults:                      # also carries scan ranges other plans read
#       inner_motor: sample.x        # dotted "registry_key.attribute"
#       outer_motor: sample.y
#     samples:
#       sample_1:
#         inner_center: 1.5          # where select_sample(1) drives inner_motor
#         outer_center: -0.2
#       sample_2:
#         inner_center: 2.0
#         outer_center: -0.2
#         inner_motor: rheometer.x   # optional, overrides defaults for this one


def _resolve_motor(motor_str: str):
    """Resolve a dotted motor string (e.g. 'sample.x') to an Ophyd device attribute."""
    parts = motor_str.split(".")
    obj = oregistry[parts[0]]
    for attr in parts[1:]:
        obj = getattr(obj, attr)
    return obj


def select_sample(env: int):
    """Move to a sample position defined in sample_info.yaml.

    Reads inner/outer motor names from the ``defaults`` section and resolves
    them through the ophyd registry, so the same function works for any stage
    configured in the YAML (sample, rheometer, huber, etc.).

    Args:
        env: Sample index matching a ``sample_{env}`` key in the YAML.
    """
    with open(expt.sample_info_file, "r") as f:
        loaded_dict = yaml.safe_load(f)

    defaults = loaded_dict["defaults"]
    sample_key = f"sample_{env}"
    sample_data = loaded_dict["samples"][sample_key]

    inner_center = sample_data["inner_center"]
    outer_center = sample_data["outer_center"]

    inner_motor_str = sample_data.get("inner_motor", defaults["inner_motor"])
    outer_motor_str = sample_data.get("outer_motor", defaults["outer_motor"])

    inner_motor = _resolve_motor(inner_motor_str)
    outer_motor = _resolve_motor(outer_motor_str)

    print(
        f"Moving {sample_key}: {inner_motor_str} -> {inner_center}, "
        f"{outer_motor_str} -> {outer_center}"
    )

    inner_motor.move(inner_center)
    outer_motor.move(outer_center)

    expt.sample_index = env
