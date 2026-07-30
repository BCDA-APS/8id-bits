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

# Named axis roles every detector's `motors` list may define. Used by
# move_detector_axes() and by master_plan.py to recognize per-protocol position
# overrides in measurement_info.yaml.
AXIS_NAMES = ["horizontal", "vertical", "swing_angle_horizontal", "swing_angle_vertical"]

# Axes select_device() itself is allowed to move. Swing angles are deliberately excluded:
# the flight-path swing (eiger4M/rigaku3M) and huber diffractometer (lambda2M) are shared/
# finicky motors, so select_device() never drives them, even when device_position.yaml
# configures a device for them. Only move_detector_axes() (called from
# master_plan.run_measurement() for explicit per-protocol overrides) may move swing axes.
TRANSLATION_AXES = ["horizontal", "vertical"]


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


def move_detector_axes(name: str, overrides: dict, timeout: float = 300):
    """Move a subset of a detector's axes to explicit positions.

    Applied on top of whatever select_device() already set, so callers only need to
    pass the axes they want to override (e.g. {"swing_angle_horizontal": 5.0}).

    Args:
        name: Detector name as used by select_device() (resolves DETECTOR_ALIASES).
        overrides: {axis_name: position} for any subset of AXIS_NAMES.

    Raises:
        ValueError: if an axis has no real device (e.g. lambda2M's horizontal/vertical
            translation, which doesn't physically exist).
    """
    config = _load_config()
    cfg = _detector_config(config, name)
    motors_cfg = cfg["motors"]

    for axis_name, position in overrides.items():
        motor = _find_motor(motors_cfg, axis_name)

        if motor.get("device") is None:
            raise ValueError(
                f"'{name}' has no '{axis_name}' axis to move (device is null in device_position.yaml)."
            )

        _resolve(motor["device"]).move(position, wait=True, timeout=timeout)


def select_device(name: str):
    """Move a beamline device to a named pre-configured position.

    Searches detectors, diagnostics, and sample_envs sections of
    device_position.yaml in order. Section-specific behaviour:

    - detectors: writes registers, updates beam-centre and motor-position
      registers in pv_registers, moves horizontal/vertical translation only
      (TRANSLATION_AXES) — never touches swing_angle_horizontal/vertical, even
      if device_position.yaml configures a device for them. Use
      move_detector_axes() to move swing axes explicitly.
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

        pv_registers.current_det_x0.put(_find_motor(motors_cfg, "horizontal").get("position", 0))
        pv_registers.current_det_y0.put(_find_motor(motors_cfg, "vertical").get("position", 0))
        pv_registers.current_db_x0.put(cfg["db_x"])
        pv_registers.current_db_y0.put(cfg["db_y"])

        translation_motors = [m for m in motors_cfg if m.get("name") in TRANSLATION_AXES]
        _move_motors(translation_motors)
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
