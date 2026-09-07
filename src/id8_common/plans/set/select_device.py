"""
Stage motion plans for 8ID beamlines.

All motor names, positions, and device-specific settings are read from
device_position.yaml so that hardware values stay out of the code.

Worked example of the YAML shape, cut down to the fields this module actually
reads (device_position.yaml's own header comment documents every field)::

    detectors:
      rigaku3M:
        motors:                            # detector axes are named
          - name: horizontal
            device: detector.x             # dotted "registry_key.attribute"
            position: -250.0               # where select_device() sends it
          - name: swing_angle_horizontal
            device: flight_path_8idi.swing # never moved by select_device()
        allow_motion: false                # optional, defaults to true
        registers:                         # optional, .put() before any motion
          softglue.enable_rigaku: '1'
    diagnostics:
      microscope:
        motors:                            # no "name:" here -- all are moved
          - device: cam_stage_8idi.x
            position: 155
    sample_envs:
      qnw:
        motors:
          - device: granite.x
            position: 923.0
        valve: granite_8idi_valve.enable   # opened, then closed, around motion

Reading that: ``select_device("qnw")`` opens ``granite_8idi_valve.enable``,
drives ``granite.x`` to 923.0, then closes the valve again.
``select_device("rigaku3M")`` puts '1' into ``softglue.enable_rigaku`` and --
because allow_motion is false for that entry -- stops there; with allow_motion
true it would also drive ``detector.x`` to -250.0, still leaving the swing
angle wherever it is.
"""

from pathlib import Path

import yaml
import time as ttime
from id8_common.expt_config import expt
from id8_common.registry import get_ophyd_object


DEVICE_POSITION_PATH = Path(__file__).parent / "device_position.yaml"

# rigaku3M_epics and rigaku3M_ftf are the same physical detector as rigaku3M,
# run with a different output format (see rigaku3m_modes.py). Neither has an
# entry of its own in device_position.yaml; they share rigaku3M's
# motors/db_x/db_y/distance/registers via these aliases.
DETECTOR_ALIASES = {
    "rigaku3M_epics": "rigaku3M",
    "rigaku3M_ftf": "rigaku3M",
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
    """Resolve a dotted path (e.g. 'cam_stage_8idi.x') from device_position.yaml.

    Thin alias for ``registry.get_ophyd_object``. This used to be a fifth,
    private copy of that walk which did a bare ``oregistry[parts[0]]`` and never
    checked ``connected`` -- so a motor whose IOC was down resolved fine here and
    only failed later, at the ``.move()``. Everything this file resolves is
    something it is about to move or put to, which is exactly the case the
    connectivity check exists for.
    """
    return get_ophyd_object(dotted)


def _motion_allowed(cfg: dict) -> bool:
    """Whether this device_position.yaml entry's motors may be driven.

    Backed by the entry's ``allow_motion`` field, which defaults to True when
    absent so an entry that predates the field keeps moving as it always did.
    Set it False to make an entry selectable but stationary -- the collision
    guard for test/commissioning, when selecting a device must not move it.
    """
    return bool(cfg.get("allow_motion", True))


def _move_motors(motors_cfg: list, timeout: float = 300):
    """Drive every entry in a `motors:` list that has both a device and a position.

    Entries missing either one are skipped silently on purpose: an axis with
    ``device: null`` does not physically exist on that detector, and an axis
    with no ``position`` (every swing angle) has no preset to go to.
    """
    for m in motors_cfg:
        if m.get("device") is None or m.get("position") is None:
            continue
        _resolve(m["device"]).move(m["position"], wait=True, timeout=timeout)


def _find_motor(motors_cfg: list, name: str):
    """Return the entry named `name` in a detector's `motors:` list."""
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

    # Raise rather than skip: this is an explicit, targeted motion request (a
    # measurement_info.yaml axis override), so silently ignoring it would run
    # the measurement with the detector somewhere other than the protocol asked
    # for. select_device() skips instead, because there the move is incidental.
    if not _motion_allowed(cfg):
        raise ValueError(
            f"'{name}' has allow_motion: false in device_position.yaml, so "
            f"move_detector_axes() will not drive {sorted(overrides)}. Set "
            f"allow_motion: true for '{name}' to allow motion again."
        )

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

    - detectors: applies the entry's `registers:` block (real ophyd signals only),
      records det_name on `expt`, moves horizontal/vertical translation only
      (TRANSLATION_AXES) — never touches swing_angle_horizontal/vertical, even
      if device_position.yaml configures a device for them. Use
      move_detector_axes() to move swing axes explicitly.
    - diagnostics: moves motors only.
    - sample_envs: opens a valve, moves motors, then closes the valve.

    An entry with ``allow_motion: false`` in device_position.yaml is selected
    without being moved: a detector still gets its registers and det_name, a
    sample_env leaves its valve shut, and nothing drives a motor. The field
    defaults to true when absent. See _motion_allowed(). (Beam centre is not
    written here at all any more -- nexus_utils reads db_x/db_y straight from
    device_position.yaml; see the note further down this function.)

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

        # The current_det_*/current_db_* EPICS registers were written here until
        # 2026-09-06 and read by nothing -- nexus_utils takes beam centre and
        # detector position straight from device_position.yaml. Removed with the
        # rest of the register retirement.

        if _motion_allowed(cfg):
            translation_motors = [m for m in motors_cfg if m.get("name") in TRANSLATION_AXES]
            _move_motors(translation_motors)
        else:
            print(f"'{name}' has allow_motion: false -- selected without moving translation axes.")

        # Run state is what det_acq_series() reads, so a standalone
        # select_device("eiger4M") at the prompt has to update it too --
        # otherwise a manual acquisition would run against whatever detector
        # the previous measurement left in the run state.
        expt.det_name = name  # keep the caller's name (e.g. "rigaku3M_epics")

    elif name in config["diagnostics"]:
        cfg = config["diagnostics"][name]
        if _motion_allowed(cfg):
            _move_motors(cfg["motors"])
        else:
            print(f"'{name}' has allow_motion: false -- selected without moving.")

    elif name in config["sample_envs"]:
        cfg = config["sample_envs"][name]
        valve_path = cfg.get("valve")

        # The valve exists only to gate this motion, so with motion off there is
        # nothing for it to gate -- leave it closed rather than cycle it.
        if not _motion_allowed(cfg):
            print(f"'{name}' has allow_motion: false -- selected without moving or opening the valve.")
            return

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
