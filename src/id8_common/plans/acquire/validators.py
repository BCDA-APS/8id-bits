"""Field checks shared by the single- and dual-detector YAML front ends.

``master_plan.py`` and ``dual_master_plan_eiger4m_rigaku3m.py`` keep their own
acquire loops and their own top-level ``validate_*`` entry points on purpose --
a dual protocol's shape (a ``detectors:`` list of legs, each with its own
timing) really is different from a serial one, and folding them together would
make both harder to read.

What they should not each own is the *field-level* checks underneath: is this a
real detector/mode, is acq_time above the mode's hardware floor, are the mode's
devices connected, is analysis_type in the allowed set, is the mesh usable.
Those were written twice and drifted apart -- see ``require_mode_devices``
below, and ``registry.get_ophyd_object`` -- so they live here now.

Every check takes plain values plus an optional ``where`` label, so the caller
owns the error prefix: the serial path raises ``acq_time must be > 0.`` and the
dual path raises ``Leg 'eiger4M': acq_time must be > 0.`` from the same call.

Import direction is one-way: this module imports ``ad_acq`` (for the ACQ_MODES
table) and is imported by the two master-plan modules. Nothing in ``ad_acq`` or
the per-detector mode modules may import it back.
"""

from pathlib import Path

import yaml

from id8_common.expt_config import expt
from id8_common.plans.acquire.ad_acq import ACQ_MODES
from id8_common.registry import get_connected_device
from id8_common.registry import get_ophyd_object

#: Analysis types a protocol may ask for. Enforced here and nowhere else --
#: dm_util.py forwards expt.analysis_type verbatim as the DM workflow's
#: argsDict["type"] and would submit any string.
VALID_ANALYSIS_TYPES = ["Multitau", "Twotime", "Both"]

#: sample_info.yaml fields a mesh scan needs. Only checked when
#: ``sample_move: yes`` -- a stationary measurement needs none of them.
MESH_SAMPLE_FIELDS = [
    "inner_motor",
    "outer_motor",
    "inner_center",
    "outer_center",
    "inner_range",
    "outer_range",
    "inner_pts",
    "outer_pts",
]


def _prefix(where):
    """Render the caller's context label, or nothing if it gave none."""
    return f"{where}: " if where else ""


# =============================================================================
# Value coercion
# =============================================================================


def read_yaml(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


def yes_no(value, field_name):
    """Normalise a YAML yes/no field to the string ``"yes"`` or ``"no"``.

    YAML 1.1 turns a bare ``no`` into the bool ``False``, so a hand-edited
    ``sample_move: no`` arrives here as a bool and still has to compare equal
    to ``"no"`` for every downstream ``== "yes"`` test.
    """
    if value is True:
        return "yes"

    if value is False:
        return "no"

    if isinstance(value, str):
        value = value.strip().lower()

        if value == "yes":
            return "yes"

        if value == "no":
            return "no"

    raise ValueError(f"{field_name} must be yes or no.")


def as_bool(value, field_name):
    """Same field shapes as ``yes_no``, returned as a real bool."""
    if isinstance(value, bool):
        return value

    return yes_no(value, field_name) == "yes"


def require_fields(mapping, fields, what, where=""):
    """Raise on the first of ``fields`` missing from ``mapping``."""
    for field in fields:
        if field not in mapping:
            raise ValueError(f"{_prefix(where)}missing {what} field: '{field}'.")


def require_positive_int(value, field_name, where=""):
    """Coerce to int and require >= 1. Returns the int."""
    number = int(value)

    if number < 1:
        raise ValueError(f"{_prefix(where)}{field_name} must be >= 1.")

    return number


def normalize_yes_no(measurement):
    """In place: coerce the two yes/no protocol fields to their string form.

    Every later check compares them with ``== "yes"``, so this has to run
    before validation, not after.
    """
    measurement["sample_move"] = yes_no(measurement["sample_move"], "sample_move")
    measurement["position_reset"] = yes_no(measurement.get("position_reset", "no"), "position_reset")


# =============================================================================
# Protocol checks
# =============================================================================


def validate_analysis_type(analysis_type, where=""):
    if analysis_type not in VALID_ANALYSIS_TYPES:
        raise ValueError(
            f"{_prefix(where)}analysis_type must be one of {VALID_ANALYSIS_TYPES} (got '{analysis_type}')."
        )


def validate_detector_mode(detector, mode, where=""):
    """Check a detector/mode pair against ACQ_MODES. Returns its mode_info."""
    if detector not in ACQ_MODES:
        raise ValueError(f"{_prefix(where)}invalid detector '{detector}'. Known: {sorted(ACQ_MODES)}")

    if mode not in ACQ_MODES[detector]:
        raise ValueError(
            f"{_prefix(where)}invalid mode '{mode}' for '{detector}'. Known: {sorted(ACQ_MODES[detector])}"
        )

    return ACQ_MODES[detector][mode]


def require_mode_devices(detector, mode, where=""):
    """Every device this mode DECLARES must be registered and connected.

    Declared means ``required_devices`` in the mode table, plus the detector
    itself -- not every device the mode's code touches.

    Returns the detector object itself (``hardware_device`` when the ACQ_MODES
    key is an alias, e.g. ``rigaku3M_epics`` -> the ``rigaku3M`` device).

    This is the check that had drifted furthest. Before 2026-09-06 there were
    four spellings of it:

    * ``ad_acq.det_acq_series`` -- required_devices AND hardware_device,
      through ``get_connected_device``.
    * ``dual_acq_eiger4m_rigaku3m.dual_acq_series`` -- same, the good one.
    * ``master_plan.validate_required_devices_connected`` -- required_devices
      only, through a bare ``oregistry[name]``, so a device skipped at startup
      raised a naked ``KeyError('softglue')`` with no hint where to look.
    * ``dual_master_plan_eiger4m_rigaku3m.validate_leg`` -- hardware_device only, so a dual
      protocol with an eiger4M External Series leg never checked softglue at
      validation time and only failed once the run was already underway.

    Both validators now call this, so a dry run rejects exactly what the real
    run would, with the registry's own "skipped at startup" message.

    It is only as good as the mode table, though. ``softglue_8id_mz2`` is
    resolved by the eiger External Series / External Enable setups and by
    ``acquire_lambda_external``, but is in no mode's ``required_devices``, so it
    is still not checked until the run has started.
    """
    mode_info = validate_detector_mode(detector, mode, where=where)

    for device_name in mode_info["required_devices"]:
        try:
            get_connected_device(device_name)
        except (KeyError, RuntimeError) as exc:
            raise RuntimeError(f"{_prefix(where)}{detector} {mode} needs '{device_name}': {exc}") from exc

    return get_connected_device(mode_info.get("hardware_device", detector))


def qmap_dir():
    """Directory a bare qmap name resolves against: <experiment>/data/.

    Not a guess -- this is where DM looks. Its own job log reads
    ``Waiting on file: /gdata/.../<experiment>/data/<qmap>``. We submit only the
    bare file name in argsDict, so this is the one place it can be.
    """
    return Path(f"{expt.mount_point}{expt.cycle_name}/{expt.experiment_name}/data")


def validate_qmap_exists(qmap_file, where=""):
    """Fail before anything moves if the named qmap is not on disk.

    Until 2026-09-08 the only check was that the name was a non-empty string, so
    a typo -- or a qmap that was simply never copied into the experiment -- stayed
    invisible until the analysis stage. The production dual plan named two qmaps
    that did not exist and dry-ran clean; a 7.5 hour acquisition would have
    completed and then produced nothing. This turns that into a dry-run error.
    """
    name = str(qmap_file).strip()

    if not name:
        raise ValueError(f"{_prefix(where)}qmap_file must not be empty.")

    directory = qmap_dir()

    if not directory.is_dir():
        # Almost always /gdata not mounted on this host, or a wrong
        # cycle_name/experiment_name in configs/experiment.yml. Say which,
        # rather than blaming the qmap.
        raise FileNotFoundError(
            f"{_prefix(where)}cannot check qmap {name!r}: {directory} does not exist. "
            f"Is /gdata mounted here, and are cycle_name/experiment_name right in "
            f"configs/experiment.yml?"
        )

    path = directory / name

    if not path.is_file():
        available = sorted(q.name for q in directory.glob("*qmap*"))
        raise FileNotFoundError(
            f"{_prefix(where)}qmap {name!r} not found at {path}. "
            f"Available in that directory: {available or 'none'}"
        )

    return path


def validate_acq_time(acq_time, detector, mode, where=""):
    """acq_time > 0, and above this mode's hardware floor. Returns the float.

    The floor is ``min_acq_time`` in the mode table -- the Rigaku ZDT bit
    depths and the Rigaku EPICS mode each have their own.
    """
    value = float(acq_time)

    if value <= 0:
        raise ValueError(f"{_prefix(where)}acq_time must be > 0.")

    min_acq_time = ACQ_MODES[detector][mode].get("min_acq_time")

    if min_acq_time is not None and value < min_acq_time:
        raise ValueError(
            f"{_prefix(where)}{detector} {mode} requires acq_time >= {min_acq_time:.2e} s (got {value:.2e} s)."
        )

    return value


# =============================================================================
# Sample mesh
# =============================================================================


def validate_sample_motion(measurement, sample, forbidden_motors=()):
    """Prove the mesh is runnable before anything is allowed to move.

    ``forbidden_motors`` is the dual path's lockout: ``huber.delta`` and
    ``huber.nu`` are positioned once by ``setup_huber_for_dual()`` and must not
    also be driven as mesh axes. The serial path passes nothing and allows any
    axis.
    """
    if measurement["sample_move"] != "yes":
        return

    require_fields(sample, MESH_SAMPLE_FIELDS, "sample")

    for role in ("inner_motor", "outer_motor"):
        if sample[role] in forbidden_motors:
            raise ValueError(
                f"sample_info.yaml sets {role} = '{sample[role]}', which this acquisition "
                f"must never move. Use a different sample axis, or set sample_move: no."
            )

        # Raises if the motor is unknown or its IOC is down -- see
        # registry.get_ophyd_object.
        get_ophyd_object(sample[role])

    require_positive_int(sample["inner_pts"], "inner_pts")
    require_positive_int(sample["outer_pts"], "outer_pts")

    expt.sample_position(int(measurement["sample_index"]))  # validates the index


def reset_sample_position(measurement):
    """Rewind this sample's mesh index to -1 when the protocol asks for it."""
    if measurement["sample_move"] != "yes":
        return

    if measurement["position_reset"] != "yes":
        return

    expt.set_sample_position(int(measurement["sample_index"]), -1)
