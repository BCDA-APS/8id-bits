"""
Shared helpers for the per-detector mode-definition modules (eiger4m_modes.py,
lambda2m_modes.py, rigaku3m_modes.py) and for ad_acq.py itself.

This module must not import any of those detector modules -- ad_acq.py
assembles ACQ_MODES from them, so a dependency in the other direction would
create an import cycle.
"""

import importlib.util
import os

import numpy as np

from id8_common.expt_config import expt
from id8_common.registry import get_connected_device
from id8_common.registry import get_ophyd_object

#: What ``from acq_helpers import *`` exports.
#:
#: ad_acq.py star-imports this module, and startup.py / startup_ophyd.py then
#: star-import ad_acq into the interactive namespace -- so anything public here
#: lands at the beamline scientist's prompt. Without this list that silently
#: included ``os``, ``importlib`` and ``np``, and every constant added to this
#: module would join them. ``expt`` and ``get_connected_device`` are re-exported
#: deliberately: they already reach the prompt today and are useful there, and
#: ``get_ophyd_object`` moved to registry.py (single definition) but is still
#: listed here so that it keeps reaching the prompt by the same chain. Nothing
#: imports it *from* this module -- every user of it, this file included, takes
#: it from registry.py by name.
__all__ = [
    "expt",
    "get_connected_device",
    "HOOK_FUNCTION_NAME",
    "active_hooks",
    "load_hooks",
    "run_hooks",
    "get_ophyd_object",
    "get_sample_position",
    "gen_folder_prefix",
    "read_sample_identity",
    "get_common_file_path",
    "get_rigaku_file_path",
    "sample_mesh_move",
]


# =============================================================================
# Hooks
# =============================================================================

HOOK_FUNCTION_NAME = "run"
active_hooks = []


def load_hooks(hooks_spec):
    """Load each hook file a protocol asks for and return its ``run`` function.

    hooks_spec is the optional ``hooks:`` block of a measurement in
    measurement_info.yaml: a list of mappings, each with a ``location`` key
    holding the path to a Python file that defines a function called ``run``
    (HOOK_FUNCTION_NAME). None or an empty list gives back an empty list.

    The files are loaded by path with importlib, not by a normal ``import``,
    because they are user scripts that live outside the package and are named
    in YAML at run time. Loading a file RUNS its top-level code, so a hook file
    should define ``run`` and little else.

    A missing file or a file with no ``run`` raises here, at load time --
    det_acq_series() calls this before it touches the shutter, so a typo in the
    protocol fails before any beam reaches the sample rather than mid-run.
    """
    callables = []

    if not hooks_spec:
        return callables

    for entry in hooks_spec:
        location = entry["location"]

        if not os.path.isfile(location):
            raise FileNotFoundError(f"Hook file not found: {location}")

        base = os.path.splitext(os.path.basename(location))[0]
        spec = importlib.util.spec_from_file_location(f"_hook_{base}", location)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        if not hasattr(module, HOOK_FUNCTION_NAME):
            raise AttributeError(f"Hook file {location} has no function '{HOOK_FUNCTION_NAME}'.")

        callables.append(getattr(module, HOOK_FUNCTION_NAME))

    return callables


def run_hooks(hook_callables):
    """Call each loaded hook function. Safe to call with None or []."""
    if not hook_callables:
        return

    for fn in hook_callables:
        fn()


# =============================================================================
# General helpers
# =============================================================================


def get_sample_position(sample_index):
    """Last mesh index used for this sample; -1 means start from the beginning.

    Persistent session state -- it has to survive between measurements and
    between sessions. Was EPICS Reg16..Reg41 until 2026-09-06; now the
    `persistent.sample_positions` block of state/run_state.yml. See expt_config.py.
    """
    return expt.sample_position(sample_index)


def _name_and_header():
    """``(header, sample_name)`` for the file name, from run state or sample_info.yaml.

    Run state wins when a measurement has set it: ``run_measurement()`` copies
    both out of the sample block it is working on, and during a multi-sample run
    that is the only thing that knows which sample is actually in the beam.

    With no measurement in flight -- a standalone alignment scan, which is the
    normal case for ``plans/align/ophyd_scan.py`` -- fall back to the same
    source ``run_measurement()`` itself reads: the ``sample_{expt.sample_index}``
    block of ``sample_info.yaml``, whose ``sample_name`` and ``header`` keys are
    exactly what the pre-2026 ``sort_qnw()``/``gen_folder_prefix()`` pair used
    (it read the index from ``pv_registers.qnw_index``, now
    ``expt.sample_index``, set by ``select_sample()``).

    Before this fallback existed, an align scan on a fresh session died with
    "'header' has not been set yet" -- and did so AFTER ``pre_align()`` and
    ``att()`` had already changed beamline state. The workaround was to set
    ``expt.header``/``expt.sample_name`` by hand, which invited invented names
    that then appear in file names for ever.
    """
    try:
        return expt.header, expt.sample_name
    except AttributeError:
        pass  # no measurement in flight -- read the sample block instead

    return read_sample_identity()


def read_sample_identity():
    """``(header, sample_name)`` read FRESH from sample_info.yaml, ignoring run state.

    Split out of :func:`_name_and_header` so a caller can ask for what the file
    says *now* rather than what the last measurement happened to leave in run
    state. ``_name_and_header`` still prefers run state -- during a multi-sample
    measurement that is the only thing that knows which sample is in the beam --
    but the align scans in ``plans/align/ophyd_scan.py`` call this directly, so
    editing sample_info.yaml between scans takes effect without a restart.

    Only ``header`` and ``sample_name`` are read. The ``inner_*``/``outer_*``
    keys in the same block describe a mesh and have nothing to do with naming.
    """
    # Imported here, not at module scope: master_plan imports this module
    # (via ad_acq), so a top-level import would be circular.
    import yaml

    from id8_common.plans.acquire.master_plan import get_sample

    path = expt.sample_info_file
    try:
        with open(path, "r") as handle:
            sample_info = yaml.safe_load(handle) or {}
        sample = get_sample(sample_info, expt.sample_index)
    except Exception as exc:
        raise RuntimeError(
            f"Cannot work out the scan name: no measurement is running, so it has to come "
            f"from sample_{expt.sample_index} in {path} -- and that could not be read ({exc}). "
            f"Either run a measurement, pick a sample with select_sample(<n>), or set "
            f"expt.header / expt.sample_name by hand."
        ) from exc

    missing = [k for k in ("header", "sample_name") if not sample.get(k)]
    if missing:
        raise RuntimeError(
            f"sample_{expt.sample_index} in {path} has no {' or '.join(missing)}. "
            f"Add it there (that is where run_measurement() reads it from too), or set "
            f"expt.header / expt.sample_name by hand."
        )
    return str(sample["header"]), str(sample["sample_name"])


def gen_folder_prefix():
    """
    Generate folder prefix from the current experiment state and attenuation.

    Uses:
        expt.header             (run state)
        expt.measurement_num    (8ideSoft:Reg1 -- the shared counter)
        expt.sample_name        (run state)
        filter_8ide.attenuation.readback

    Example:
        header = "A"
        measurement_num = 12
        sample_name = "G10"
        attenuation = 7

        returns "A0012_G10_a0007"

    The measurement number increments once per call. Its store is the EPICS
    register 8ideSoft:Reg1, not state/run_state.yml -- see PV_FIELDS in
    expt_config.py for why, and for the trap that follows from it: ONE counter
    feeds TWO naming streams. det_acq_series() and dual_acq_series() name their
    output <experiment>/data/A####..., while the align scans (ophyd_scan.py,
    scan_8id.py) call this same function and name theirs
    <experiment>/data/bluesky/A####..., so the highest number already on disk
    may be under data/bluesky/ rather than data/.
    """
    filter_beam = get_connected_device("filter_8ide")

    header, sample_name = _name_and_header()
    meas_num = expt.measurement_num
    att_level = int(filter_beam.attenuation.readback.get())

    folder_prefix = f"{header}{meas_num:04d}_{sample_name}_a{att_level:04d}"

    expt.measurement_num = meas_num + 1

    return folder_prefix


def get_common_file_path(file_header, file_name):
    """Absolute output DIRECTORY for one measurement, on the detector's own mount.

    Despite the name this is a directory, not a file: its last component is
    file_name, and the caller then puts file_name into the HDF plugin
    separately (or builds "<dir>/<file_name>_metadata.hdf" beside it).

    use_subfolder "yes" gives every repeat of a measurement its own
    file_header/ directory under data/; "no" puts them all directly in data/.
    """
    cycle_name = expt.cycle_name
    exp_name = expt.experiment_name
    mount_point = expt.mount_point
    use_subfolder = expt.use_subfolder

    if use_subfolder == "yes":
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == "no":
        file_path = f"{mount_point}{cycle_name}/{exp_name}/data/{file_name}"
    else:
        raise ValueError("use_subfolder must be yes or no")

    return file_path


def get_rigaku_file_path(file_header, file_name):
    """Same output directory as get_common_file_path(), returned twice over.

    Returns (file_path, full_path):
        file_path -- relative to <mount_point>/<cycle_name>. This is the form
                     handed to the IOC in cam.fast_file_path, not the absolute
                     one.
        full_path -- the same place seen from this workstation, for os.makedirs
                     and for building the metadata filename.

    Both are directories, as in get_common_file_path().
    """
    cycle_name = expt.cycle_name
    exp_name = expt.experiment_name
    mount_point = expt.mount_point
    use_subfolder = expt.use_subfolder

    if use_subfolder == "yes":
        file_path = f"{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == "no":
        file_path = f"{exp_name}/data/{file_name}"
    else:
        raise ValueError("use_subfolder must be yes or no")

    full_path = f"{mount_point}/{cycle_name}/{file_path}"

    return file_path, full_path


# =============================================================================
# Sample motion
# =============================================================================


def sample_mesh_move():
    """
    Move to the next point in a 2D mesh.

    Values read from the run state (expt_config.py), set by
    master_plan.run_measurement() from measurement_info.yaml:
        expt.sample_move
        expt.inner_motor    expt.outer_motor
        expt.inner_center   expt.outer_center
        expt.inner_range    expt.outer_range
        expt.inner_pts      expt.outer_pts

    and from the persistent session state (the `persistent:` block of
    state/run_state.yml -- was EPICS Reg6 and Reg16..Reg41):
        expt.sample_index
        expt.sample_position(sample_index)

    sample_move:
        "yes" -> move sample
        "no"  -> no sample motion and no stored-position update

    sample_index:
        1-based key into the per-sample position table (a dict now, keyed by
        index; was the 26 named registers sample1_pos..sample26_pos).

    The stored position is the last used zero-based mesh index.
    The inner axis moves fastest.

    inner_range and outer_range are interpreted as full scan widths.
    """
    if expt.sample_move != "yes":
        return

    sample_index = expt.sample_index
    last_index = get_sample_position(sample_index)

    inner_motor_name = expt.inner_motor
    outer_motor_name = expt.outer_motor

    inner_center = expt.inner_center
    outer_center = expt.outer_center

    inner_range = expt.inner_range
    outer_range = expt.outer_range

    inner_pts = expt.inner_pts
    outer_pts = expt.outer_pts

    total_pts = inner_pts * outer_pts

    inner_positions = np.linspace(
        inner_center - inner_range / 2,
        inner_center + inner_range / 2,
        inner_pts,
    )

    outer_positions = np.linspace(
        outer_center - outer_range / 2,
        outer_center + outer_range / 2,
        outer_pts,
    )

    # Step to the next mesh point, wrapping round to the first one after the
    # last. A sample that has never been used has last_index -1, so it starts
    # at point 0.
    pos_index = (last_index + 1) % total_pts

    # Turn that flat point number back into a row and a column. The inner axis
    # moves fastest, so it takes the remainder and the outer axis the quotient.
    inner_index = pos_index % inner_pts
    outer_index = pos_index // inner_pts

    inner_pos = inner_positions[inner_index]
    outer_pos = outer_positions[outer_index]

    inner_motor = get_ophyd_object(inner_motor_name)
    outer_motor = get_ophyd_object(outer_motor_name)

    print(f"Moving {outer_motor.name} to {outer_pos}")
    outer_motor.move(outer_pos)

    print(f"Moving {inner_motor.name} to {inner_pos}")
    inner_motor.move(inner_pos)

    expt.set_sample_position(sample_index, pos_index)
