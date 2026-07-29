"""
Stage motion plans for 8ID beamlines.

All motor names, positions, and device-specific settings are read from
device_position.yaml so that hardware values stay out of the code.
"""

from pathlib import Path

import yaml
import time as ttime
from apsbits.core.instrument_init import oregistry

pv_registers = oregistry["pv_registers"]

DEVICE_POSITION_PATH = Path(__file__).parent / "device_position.yaml"

# rigaku3M_epics is the same physical detector as rigaku3M, run in a different
# acquisition mode. It has no entry of its own in device_position.yaml; it shares
# rigaku3M's motors/db_x/db_y/distance/registers via this alias.
DETECTOR_ALIASES = {
    "rigaku3M_epics": "rigaku3M",
}


def _load_config():
    with open(DEVICE_POSITION_PATH, "r") as f:
        return yaml.safe_load(f)


def _detector_config(config: dict, det_name: str):
    """Look up a detector's config block in device_position.yaml, resolving DETECTOR_ALIASES."""
    key = DETECTOR_ALIASES.get(det_name, det_name)
    return config["detectors"][key]


def _resolve(dotted: str):
    """Resolve a dotted path (e.g. 'cam_stage_8idi.x') via the ophyd registry."""
    parts = dotted.split(".")
    obj = oregistry[parts[0]]
    for attr in parts[1:]:
        obj = getattr(obj, attr)
    return obj


def _move_motors(motors_cfg: list, timeout: float = 300):
    for m in motors_cfg:
        if m.get("device") is None or m.get("position") is None:
            continue
        _resolve(m["device"]).move(m["position"], wait=True, timeout=timeout)


def _find_motor(motors_cfg: list, name: str):
    for m in motors_cfg:
        if m["name"] == name:
            return m
    raise KeyError(f"No motor named '{name}' in device_position.yaml config.")


def select_device(name: str):
    """Move a beamline device to a named pre-configured position.

    Searches detectors, diagnostics, and sample_envs sections of
    device_position.yaml in order. Section-specific behaviour:

    - detectors: writes registers, updates beam-centre and motor-position
      registers in pv_registers, skips motion if already active.
    - diagnostics: moves motors only.
    - sample_envs: opens a valve, moves motors, then closes the valve.

    Args:
        name: Position key in device_position.yaml
              (e.g. ``"eiger4M"``, ``"microscope"``, ``"qnw"``).
    """
    config = _load_config()
    resolved_name = DETECTOR_ALIASES.get(name, name)

    if resolved_name in config["detectors"]:
        cfg = _detector_config(config, name)

        # if pv_registers.det_name.get() == name:
        #     return

        motors_cfg = cfg["motors"]

        for reg_path, value in cfg.get("registers", {}).items():
            _resolve(reg_path).put(value)

        pv_registers.current_det_x0.put(_find_motor(motors_cfg, "horizontal")["position"])
        pv_registers.current_det_y0.put(_find_motor(motors_cfg, "vertical")["position"])
        pv_registers.current_db_x0.put(cfg["db_x"])
        pv_registers.current_db_y0.put(cfg["db_y"])

        _move_motors(motors_cfg)
        pv_registers.det_name.put(name)  # keep the caller's name (e.g. "rigaku3M_epics")

    elif name in config["diagnostics"]:
        _move_motors(config["diagnostics"][name]["motors"])

    elif name in config["sample_envs"]:
        cfg = config["sample_envs"][name]
        valve_path = cfg.get("valve")

        if valve_path:
            _resolve(valve_path).put(1, wait=True)
        ttime.sleep(5)  # wait for valve to open

        _move_motors(cfg["motors"])
        ttime.sleep(5)  # wait after move

        if valve_path:
            _resolve(valve_path).put(0, wait=True)
        ttime.sleep(1)  # wait for valve to close

    else:
        all_names = (
            list(config["detectors"])
            + list(DETECTOR_ALIASES)
            + list(config["diagnostics"])
            + list(config["sample_envs"])
        )
        print(f"Unknown device '{name}'. Available: {all_names}")
