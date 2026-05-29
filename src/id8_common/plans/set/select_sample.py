"""
Select and move to a sample position using the sample_info.yaml configuration.
"""

import yaml
from pathlib import Path

from apsbits.core.instrument_init import oregistry

pv_registers = oregistry["pv_registers"]

SAMPLE_INFO_PATH = Path("/home/beams10/8IDIUSER/bluesky/src/user_plans/sample_info.yaml")


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
    with open(SAMPLE_INFO_PATH, "r") as f:
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

    pv_registers.sample_index.put(env)
