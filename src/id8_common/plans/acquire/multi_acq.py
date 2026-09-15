"""
Parallel multi-detector acquisition: several detectors in ONE beam window.

Runs them in the same beam window instead of back to back, so a set of long
measurements costs roughly the slowest one rather than the sum of all of them.
The detectors are NOT frame-synced -- they merely overlap. Each writes its own
file into its own folder, its own NeXus metadata, and submits its own DM job.

One repeat, start to finish:

    open the shutter  ->  arm every detector  ->  wait for the Rigaku to
      confirm it started  ->  wait for every cam  ->  close the shutter
      ->  wait for every HDF plugin  ->  metadata + DM, one leg at a time

There are NO THREADS here and nothing runs out of order. Arming a detector is
an EPICS put that tells the hardware to start and returns immediately, so
multi_acq_series() arms one leg, then the next, then the next, and from that
moment every detector is integrating while the Python side sits in a single
polling loop (wait_all). Read it top to bottom like any other plan: the only
thing happening concurrently is the hardware.

Supported (device, mode) pairs are the keys of MULTI_LEGS. Anything else is
rejected.

SAME YAML AS THE SERIAL PATH. A parallel protocol is read from the same
measurement_info.yaml and sample_info.yaml as a single-detector one, validated
by the same master_plan.py, and told apart only by carrying a ``detectors:``
list where a serial protocol has a scalar ``detector:`` (master_plan.
is_multi_protocol). The mode tables and the per-leg setup functions are shared
verbatim; this module replaces only the acquire half of ACQ_MODES, which cannot
be run once per leg because it owns the shutter.

History, because this file has moved before. It was
trio_acq_rigaku3m_eiger4m_lambda2m.py until 2026-09-11, when it folded into
ad_acq.py; what actually needed to go was the trio_master_plan_* front end
beside it, which carried its own copy of master_plan.py's validation and run
expansion and its own YAML. That duplication is gone and none of it comes back
here -- only the acquisition code moved out again, on 2026-09-14.

What made a parallel run possible at all: setup_rigaku_epics() stopped driving
the fast shutter (softglue MUX on + trigger_mode 'Start with Trigger'). While
the Rigaku gated the beam through softglue, the shutter belonged to one
nominated leg, every other leg had to be armed after that leg confirmed it was
running, and none of it resembled the serial path. Now every supported mode is
internally timed and none of them touches softglue.

Imports run one way: this module imports ad_acq (for ACQ_MODES, which ad_acq
assembles from the per-detector mode tables); ad_acq imports nothing from here.
"""

import warnings
from contextlib import contextmanager
from datetime import datetime
import time as ttime

from id8_common.expt_config import expt
from id8_common.plans.acquire.acq_helpers import gen_folder_prefix
from id8_common.plans.acquire.acq_helpers import get_connected_device
from id8_common.plans.acquire.acq_helpers import sample_mesh_move
from id8_common.plans.acquire.acq_wait import DEFAULT_STOP_TIMEOUT
from id8_common.plans.acquire.acq_wait import STATE_IDLE
from id8_common.plans.acquire.acq_wait import acq_timeout
from id8_common.plans.acquire.acq_wait import cam_fault
from id8_common.plans.acquire.acq_wait import stop_acquiring
from id8_common.plans.acquire.acq_wait import wait_until
from id8_common.plans.acquire.ad_acq import ACQ_MODES
from id8_common.plans.acquire.ad_acq import cleanup_acquisition
from id8_common.plans.set.select_device import select_device
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import post_align
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.registry import get_ophyd_object
from id8_common.utils.dm_util import dm_run_job
from id8_common.utils.dm_util import dm_setup
from id8_common.utils.nexus_writer import create_nexus_format_metadata
# The composed 8-ID schema, built from Miaoqi Chu's core factories.
from id8_common.utils.xpcs_schema import xpcs_schema


# Wall-clock allowance for a cam to finish, as a multiple of the expected
# acquisition time plus a fixed floor. Generous on purpose: the point of a
# timeout here is to stop one hung detector from stranding the others' data, not
# to police timing.
CAM_TIMEOUT_FACTOR = 3.0
CAM_TIMEOUT_FLOOR = 60.0

#: Seconds to give a detector to confirm it has started, once armed. The Rigaku
#: is the reason this is not instant: it stages internally for a few seconds
#: after cam.acquire.put(1) before it reports anything.
#:
#: Fixed, not a per-leg YAML option. It was overridable until 2026-09-14 and
#: every protocol ever written set it to exactly this value.
START_TIMEOUT = 30.0

#: Seconds to give the HDF plugins to flush after the cams stop. 600, which is
#: what every protocol asked for while this was a per-leg YAML option (dropped
#: 2026-09-14) -- the code default of 300 was never the value in use.
HDF_TIMEOUT = 600.0

#: Slept once after arming, before anyone is asked whether they are still busy.
#: For an ADCore cam the arm and the busy flag are the same PV, so a poll that
#: beats the driver to it reads 0 and the shutter closes on an acquisition that
#: has only just begun. Matches the sleep the serial acquire functions already do
#: at the same point.
ARM_SETTLE = 0.5

# huber axes that nothing inside a parallel acquisition may drive.
#
# Both are positioned exactly once per measurement, before acquisition starts, by
# master_plan.setup_huber_for_multi(). From that point on nothing may touch them:
# not a leg's `motors:` block, not the sample mesh. Both routes are refused at
# validation time and again here at run time.
FORBIDDEN_MOTORS = ("huber.delta", "huber.nu")


# -----------------------------------------------------------------------------
# NeXus metadata override paths
# -----------------------------------------------------------------------------
# create_runtime_metadata_dict() derives detector geometry from
# device_position.yaml via expt.det_name, which can only describe one detector at
# a time. It applies its additional_metadata argument LAST, so passing the
# per-leg values through that existing hook overrides every geometry field
# without touching nexus_utils.py.
#
# The cost of that is these path strings duplicating ones that live in
# nexus_utils.py. If a path is renamed there, an override would silently stop
# applying and the file would get the wrong detector's geometry. Instead, drift
# is caught twice: a warning at import (below) and _assert_override_paths(),
# which multi_acq_series() calls before it touches anything, so a parallel run
# refuses to start. update_schema_at_runtime() walks xpcs_schema, so xpcs_schema
# is what a path has to exist in.

OVERRIDE_PATHS = {
    "detector_name": "/entry/instrument/detector_1/detector_name",
    "position_x": "/entry/instrument/detector_1/position_x",
    "position_y": "/entry/instrument/detector_1/position_y",
    "beam_center_x": "/entry/instrument/detector_1/beam_center_x",
    "beam_center_y": "/entry/instrument/detector_1/beam_center_y",
    "beam_center_position_x": "/entry/instrument/detector_1/beam_center_position_x",
    "beam_center_position_y": "/entry/instrument/detector_1/beam_center_position_y",
    "qmap_file": "/entry/instrument/detector_1/qmap_file",
    "distance": "/entry/instrument/detector_1/distance",
    "x_pixel_size": "/entry/instrument/detector_1/x_pixel_size",
    "y_pixel_size": "/entry/instrument/detector_1/y_pixel_size",
    "flightpath_swing": "/entry/instrument/detector_1/flightpath_swing",
    "flightpath_swing_vertical": "/entry/instrument/detector_1/flightpath_swing_vertical",
}


def _schema_has(path):
    """True when a slash-separated NeXus path such as /entry/instrument/... exists in xpcs_schema."""
    node = xpcs_schema

    for component in path.lstrip("/").split("/"):
        if not isinstance(node, dict) or component not in node:
            return False
        node = node[component]

    return True


def missing_override_paths():
    """Override paths that no longer exist in xpcs_schema."""
    return sorted(p for p in OVERRIDE_PATHS.values() if not _schema_has(p))


def _assert_override_paths():
    """Raise if the schema has drifted. Called at acquisition time, not at import."""
    missing = missing_override_paths()

    if missing:
        raise RuntimeError(
            "multi_acq_series() writes NeXus fields that no longer exist in xpcs_schema: "
            f"{missing}. The schema changed under this module -- update OVERRIDE_PATHS "
            "to match id8_common/utils/nexus_utils.py before running a parallel acquisition."
        )


# Warn at import, enforce at acquisition. This module is imported by startup.py, so a hard
# failure here would take down the whole session -- including single-detector runs, which do
# not use any of this. Schema drift must not be able to do that.
if missing_override_paths():
    warnings.warn(
        f"multi_acq: NeXus override paths missing from xpcs_schema: {missing_override_paths()}. "
        "Parallel acquisition will refuse to run until OVERRIDE_PATHS is updated. "
        "Single-detector acquisition is unaffected.",
        RuntimeWarning,
        stacklevel=2,
    )


# -----------------------------------------------------------------------------
# Per-leg arm and poll behaviour
# -----------------------------------------------------------------------------
# The serial acquire functions above interleave showbeam()/blockbeam() with a
# blocking wait, so they cannot be used here -- whichever finished first would
# close the shutter on the others. These callables are the same detector handling
# with the shutter and the waiting lifted out, which is what lets the shutter
# contract live in one place.
#
# Each "arm" entry does two EPICS puts inside one lambda by wrapping them in a
# tuple: both run, left to right, and the tuple itself is discarded. Their order
# is deliberate -- the HDF plugin is put into capture first so it is already
# listening when the cam starts producing frames.
#
# needs_start_confirm is the one entry worth reading closely.
#
# For an ADCore cam (Eiger, Lambda) the SAME PV is the arm and the busy flag:
# cam.acquire goes to 1 inside the put and back to 0 when the series ends. There
# is nothing to confirm -- "armed" and "busy" are one state -- so those legs go
# straight into the busy wait, which is what the code this grew out of did.
#
# The Rigaku is different and it matters. detector_state stays at Idle for
# SECONDS after cam.acquire.put(1) while the detector stages, so "not busy yet"
# and "finished" read identically. A wait that only asked "is anyone still busy?"
# would find every leg idle, close the shutter, and truncate the Rigaku's
# acquisition before it took a frame -- silently, since nothing errors. So that
# leg gets a confirmed-start phase of its own, run with the shutter already open,
# before anybody waits for anybody to finish.
#
# Each "arm" entry does two EPICS puts inside one lambda by wrapping them in a
# tuple: both run, left to right, and the tuple itself is discarded. Their order
# is deliberate -- the HDF plugin is put into capture first so it is already
# listening when the cam starts producing frames.

MULTI_LEGS = {
    ("eiger4M", "Internal Series"): {
        "arm": lambda d: (d.hdf1.capture.put(1), d.cam.acquire.put(1)),
        "needs_start_confirm": False,
        "started": lambda d: d.cam.acquire.get() == 1,
        "cam_busy": lambda d: d.cam.acquire.get() == 1,
        "hdf_busy": lambda d: d.hdf1.capture.get() == 1,
    },
    ("lambda2M", "Internal"): {
        # Same shape as the Eiger's. The one thing the Lambda needs that the
        # Eiger does not -- being taken out of live/TV mode first -- is done in
        # setup_lambda_internal() via stop_lambda_live(), so it happens before
        # the shutter opens rather than inside the beam window.
        "arm": lambda d: (d.hdf1.capture.put(1), d.cam.acquire.put(1)),
        "needs_start_confirm": False,
        "started": lambda d: d.cam.acquire.get() == 1,
        "cam_busy": lambda d: d.cam.acquire.get() == 1,
        "hdf_busy": lambda d: d.hdf1.capture.get() == 1,
    },
    ("rigaku3M_epics", "EPICS"): {
        "arm": lambda d: (d.hdf1.capture.put(1), d.cam.acquire.put(1)),
        "needs_start_confirm": True,
        # "not Idle" rather than "== Acquire" so the readout and saving states
        # count as started too.
        "started": lambda d: d.cam.detector_state.get() != STATE_IDLE,
        # Same test for busy, and the terminal states are NOT excluded on
        # purpose -- the fault check in wait_all() catches those and raises,
        # which is what a bare "!= 0" used to miss: a dead detector read as busy
        # forever, with the shutter open, until the cam timeout expired.
        "cam_busy": lambda d: d.cam.detector_state.get() != STATE_IDLE,
        "hdf_busy": lambda d: d.hdf1.capture.get() == 1,
    },
}


# -----------------------------------------------------------------------------
# Is this mode allowed to share a beam window?
# -----------------------------------------------------------------------------
# MULTI_LEGS says which pairs are supported, but only by being a list -- add a
# row to it and the row is supported, whether or not the mode can actually
# survive being run that way. The two checks below test the PROPERTY instead,
# reading the drives_shutter / self_paced flags each mode declares in its own
# table (eiger4m_modes.py, lambda2m_modes.py, rigaku3m_modes.py).
#
# Both are enforced, in three places:
#
#   * master_plan.validate_leg(), so a dry run rejects a bad protocol before
#     anything moves;
#   * prepare_legs() below, so a direct multi_acq_series() call that skipped
#     validation is rejected too;
#   * assert_multi_legs_sane(), over the whole MULTI_LEGS table, so a row added
#     to it that contradicts its own mode table is caught at import (a warning)
#     and again before any parallel acquisition starts (an exception).
#
# The failure this guards against is not a crash, it is silent bad data: a
# shutter-driving leg fights the one showbeam()/blockbeam() around the set, and
# a mode that needs pacing simply never gets it. Neither errors on its own.


def parallel_objection(device, mode):
    """Why this (device, mode) may not share a beam window, or None if it may.

    Returns a sentence naming the reason, ready to be prefixed with the leg's
    label by the caller. Assumes the pair exists in ACQ_MODES -- callers check
    that first, and get a better message from validate_detector_mode() when it
    does not.
    """
    mode_info = ACQ_MODES[device][mode]

    missing = [flag for flag in ("drives_shutter", "self_paced") if flag not in mode_info]

    if missing:
        # Fail closed. A mode added without the flags is refused from parallel
        # runs rather than assumed safe -- the assumption is the expensive way
        # round, and the fix is one line in the mode table.
        return (
            f"{device} {mode} does not declare {' or '.join(missing)} in its mode table, "
            f"so whether it can share a beam window with another detector is unknown. "
            f"Add the flag(s) to that mode's entry -- see the comment above the table."
        )

    if mode_info["drives_shutter"]:
        reason = (
            f"{device} {mode} gates the fast shutter itself, through softglue, so it "
            f"cannot share a beam window: the single showbeam()/blockbeam() around the "
            f"whole set would fight the detector's own trigger path."
        )

        # The case worth naming outright, because the two Rigaku families look
        # interchangeable in a protocol and are not.
        if device.startswith("rigaku3M"):
            reason += (
                " Use device: rigaku3M_epics with mode: EPICS for parallel acquisition. "
                "rigaku3M (sparsified .bin) and rigaku3M_ftf (fast-transfer .h5) both set "
                "softglue.enable_rigaku = '1' and trigger_mode 'Start with Trigger', and "
                "are serial-only -- run them through a scalar detector: protocol."
            )

        return reason

    if not mode_info["self_paced"]:
        return (
            f"{device} {mode} is not self-paced -- it needs a softglue pulse train or a "
            f"software trigger per frame while it runs -- so it cannot be armed once and "
            f"left going while the plan waits on the other detectors. Only internally "
            f"timed modes can share a beam window."
        )

    return None


def assert_parallel_safe(device, mode, where=""):
    """Raise unless this (device, mode) may share a beam window."""
    objection = parallel_objection(device, mode)

    if objection:
        raise ValueError(f"{where}: {objection}" if where else objection)


def inconsistent_multi_legs():
    """MULTI_LEGS rows their own mode table says cannot share a beam window."""
    bad = {}

    for device, mode in MULTI_LEGS:
        if device not in ACQ_MODES or mode not in ACQ_MODES.get(device, {}):
            bad[(device, mode)] = f"{device} {mode} is not in ACQ_MODES at all."
            continue

        objection = parallel_objection(device, mode)

        if objection:
            bad[(device, mode)] = objection

    return bad


def assert_multi_legs_sane():
    """Raise if MULTI_LEGS lists a pair its mode table forbids. Called per acquisition."""
    bad = inconsistent_multi_legs()

    if bad:
        detail = " ".join(f"[{d}/{m}] {why}" for (d, m), why in sorted(bad.items()))
        raise RuntimeError(
            f"MULTI_LEGS lists {len(bad)} pair(s) that cannot share a beam window: {detail}"
        )


# Warn at import, raise at acquisition -- the same split as the OVERRIDE_PATHS
# check above, and for the same reason: this module is imported by startup.py, so
# a hard failure here would take down a session that only ever runs
# single-detector acquisitions.
if inconsistent_multi_legs():
    warnings.warn(
        f"multi_acq: MULTI_LEGS contains pairs their mode tables forbid in a shared beam "
        f"window: {sorted(inconsistent_multi_legs())}. Parallel acquisition will refuse "
        f"to run until this is fixed. Single-detector acquisition is unaffected.",
        RuntimeWarning,
        stacklevel=2,
    )


def multi_leg_behaviour(device, mode):
    """Look up the arm/poll callables for a (device, mode) pair."""
    key = (device, mode)

    if key not in MULTI_LEGS:
        supported = ", ".join(f"{d}/{m}" for d, m in sorted(MULTI_LEGS))
        raise ValueError(
            f"'{device}' in mode '{mode}' cannot run in a parallel measurement. Supported: {supported}"
        )

    return MULTI_LEGS[key]


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def timestamp():
    """Log timestamp, matching the format det_acq_series() prints."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def resolve_value(value):
    """Resolve a geometry entry that is either a literal number or a dotted ophyd path.

    'huber.nu' reads that motor's current position; 5.0 is used as-is. Lets a leg
    report a live motor position in its metadata without needing a
    device_position.yaml entry.
    """
    if isinstance(value, str):
        return float(get_ophyd_object(value).position)

    return float(value)


def wait_all(legs, key, timeout, description, want=False, settle=0.0):
    """Wait until `key` reads `want` for EVERY leg, or raise.

    `key` names one of the per-leg callables in MULTI_LEGS. Two shapes are used:

        want=True   "started"  -- every detector has confirmed it is running
        want=False  "cam_busy" / "hdf_busy" -- nobody is still working

    `settle` is slept once before the first poll. It is what stops a just-armed
    ADCore cam being read as finished: cam.acquire is both the arm and the busy
    flag, so a poll that lands before the driver has latched the put sees 0 and
    concludes the acquisition is over. The serial acquire functions in the
    per-detector modules sleep 0.5 s after arming for exactly this reason.

    Fault-aware: each leg contributes a cam_fault() check, so a detector that
    reaches Error/Disconnected/Aborted raises within one poll instead of reading
    as busy until `timeout` expires. These loops run with the shutter OPEN (see
    multi_acq_series), which is why waiting out a multi-minute timeout on a dead
    detector is the wrong behaviour.
    """
    faults = []

    for leg in legs:
        if leg.get("det") is not None:
            faults.append(cam_fault(leg["det"].cam, leg["label"]))

    def every_leg_agrees():
        for leg in legs:
            if bool(leg["behaviour"][key](leg["det"])) is not want:
                return False

        return True

    return wait_until(
        every_leg_agrees,
        timeout=timeout,
        what=description,
        faults=tuple(faults),
        settle=settle,
    )


def stop_live_mode(legs):
    """Take every leg out of live/TV mode, and confirm it, BEFORE any setup runs.

    This is the guard for a bug that cost real measurements: a detector left
    free-running from its GUI ("TV mode") ignores a write to acquire_time, so
    the setup function returns happily and the acquisition then runs at the live
    view's exposure rather than the protocol's. Nothing downstream notices --
    the frames are written, the metadata records the REQUESTED value, and the
    exposure the data was actually taken at is simply lost.

    Ordering is the whole point. This must run before leg["setup"]() pushes
    acquire_time, and both must run before showbeam(), so a camera that refuses
    to stop costs a raised exception rather than beam on the sample.

    Cheap and idempotent: a camera already idle is a single PV read.

    There used to be a second half to this -- confirm_acq_time(), which read every
    camera's acquire_time back after setup and refused to run on a mismatch. It
    was removed on 2026-09-11: stopping every camera first removes the only known
    reason a write does not stick, and the readback check cost a settling wait per
    repeat for a symptom that can no longer occur.
    """
    for leg in legs:
        if stop_acquiring(leg["det"].cam, leg["label"], timeout=leg["stop_timeout"]):
            print(f"{timestamp()}, {leg['label']} was still acquiring (live/TV mode?) -- stopped it")


def cam_timeout_for(legs):
    """Allowance for the slowest leg's cam to finish."""
    expected = max(leg["acq_time"] * leg["num_frames"] for leg in legs)

    return acq_timeout(expected, 1, factor=CAM_TIMEOUT_FACTOR, floor=CAM_TIMEOUT_FLOOR)


def metadata_overrides(leg):
    """Per-leg NeXus values, passed through create_nexus_format_metadata's additional_metadata.

    Applied after everything create_runtime_metadata_dict() derived from
    device_position.yaml, so whatever is returned here wins.

    Only what the leg actually states is overridden. A leg sitting at its normal
    device_position.yaml mount needs no `geometry` block at all -- the values that
    file already holds are correct, and this returns just the two things that are
    per-leg by definition (which detector wrote the file, and which qmap goes with
    it). Override a geometry field only when this detector is somewhere
    device_position.yaml does not describe, e.g. an Eiger remounted on the huber arm.

    Lengths are millimetres in the YAML and metres in the file, matching what
    nexus_utils.py does with the same fields.
    """
    geometry = leg.get("geometry") or {}

    values = {
        # leg["device"], not leg["label"]: the label is deliberately aliased for
        # FILE PATHS (rigaku3M_epics -> rigaku3M, see master_plan.leg_label()),
        # but detector_name is what tells a downstream reader which OUTPUT FORMAT
        # produced the file -- .h5 for rigaku3M_epics vs sparsified .bin for
        # rigaku3M. The single-detector path records the full device key, and
        # until 2026-09-06 the parallel path silently recorded the aliased one, so
        # the same acquisition was described differently by the two paths.
        "detector_name": leg["device"],
        "qmap_file": leg["qmap_file"],
    }

    # geometry key -> (override name, scale applied to the YAML value)
    geometry_to_override = {
        "db_x": ("beam_center_x", 1.0),
        "db_y": ("beam_center_y", 1.0),
        "distance": ("distance", 1.0),
        "position_x": ("position_x", 1e-3),
        "position_y": ("position_y", 1e-3),
        "beam_center_position_x": ("beam_center_position_x", 1e-3),
        "beam_center_position_y": ("beam_center_position_y", 1e-3),
        "swing_horizontal": ("flightpath_swing", 1.0),
        "swing_vertical": ("flightpath_swing_vertical", 1.0),
    }

    for geometry_key, override in geometry_to_override.items():
        if geometry_key not in geometry:
            continue

        override_name, scale = override
        values[override_name] = resolve_value(geometry[geometry_key]) * scale

    if "pixel_size" in geometry:
        pixel_size = resolve_value(geometry["pixel_size"])
        values["x_pixel_size"] = pixel_size
        values["y_pixel_size"] = pixel_size

    # Everything above is keyed by short name; the caller wants NeXus paths.
    overrides = {}

    for name, value in values.items():
        overrides[OVERRIDE_PATHS[name]] = value

    return overrides


@contextmanager
def swapped_registers(leg):
    """Point the shared run-state fields at one leg, then put them back.

    dm_run_job() reads det_name/qmap_file/analysis_type/workflow_name from `expt`
    and has no override argument, and create_runtime_metadata_dict() needs
    det_name to resolve to a real device_position.yaml key before its
    additional_metadata override is applied. Both calls happen after the parallel
    window, one leg at a time, so the values are never contended -- and they are
    restored afterwards so a parallel run leaves no trace.

    Snapshot and restore both on `expt`. Until 2026-09-06 this saved from the
    EPICS registers, set on expt, and restored to the registers -- so the restore
    never touched what the readers actually read, and after a parallel run every
    later measurement in the session was stamped with the last leg's detector,
    qmap and analysis type.
    """
    names = ["det_name", "qmap_file", "analysis_type", "workflow_name"]

    # snapshot_run(), not getattr(): in a parallel run these fields have no global
    # value between legs, and a plain getattr would raise AttributeError on the
    # very first leg. The sentinel lets restore_run() put "absent" back.
    saved = expt.snapshot_run(names)

    try:
        expt.det_name = leg["device"]
        expt.qmap_file = leg["qmap_file"]
        expt.analysis_type = leg["analysis_type"]

        if leg.get("workflow_name"):
            expt.workflow_name = leg["workflow_name"]

        yield
    finally:
        expt.restore_run(saved)


def move_leg_motors(leg):
    """Move whatever the leg's `motors` block names. Absent axes are left alone.

    huber.delta and huber.nu are refused here even though master_plan already
    rejects them at validation time -- multi_acq_series() can be driven directly,
    bypassing that.
    """
    for dotted, position in (leg.get("motors") or {}).items():
        if dotted in FORBIDDEN_MOTORS:
            raise ValueError(
                f"Leg '{leg['label']}': '{dotted}' cannot be moved by a parallel acquisition. "
                f"Both huber axes are positioned once before acquisition by "
                f"master_plan.setup_huber_for_multi()."
            )

        motor = get_ophyd_object(dotted)
        print(f"Moving {dotted} to {position}")
        motor.move(float(position), wait=True)


def write_leg_metadata(leg):
    """Write one leg's NeXus metadata, then clear the path so it is written only once.

    Self-contained on purpose -- it does its own swapped_registers(), so
    cleanup_multi() can call it without the caller having set anything up first.
    The swap and the per-leg overrides are both required: without them the file is
    stamped with whatever detector expt.det_name happens to name, which in a
    parallel run is another leg most of the time. That is also why cleanup cannot
    write these files through cleanup_acquisition's metadata_fname argument, which
    knows nothing about legs.

    Does nothing once the metadata is on disk, so an abort later in the same
    repeat cannot overwrite a complete file.
    """
    metadata_fname = leg.get("metadata_fname")

    if not metadata_fname or leg.get("det") is None:
        return

    print(f"{timestamp()}, Writing metadata, {leg['file_name']}")

    with swapped_registers(leg):
        create_nexus_format_metadata(
            metadata_fname,
            det=leg["det"],
            additional_metadata=metadata_overrides(leg),
        )

    leg["metadata_fname"] = None


def cleanup_multi(legs):
    """Close the shutter, stop every leg, write metadata for whatever was in flight.

    Safe to call from a partially-started state.
    """
    try:
        blockbeam()
    except Exception as exc:
        print(f"Could not block beam during cleanup: {exc}")

    for leg in legs:
        cleanup_acquisition(leg.get("det"), leg.get("mode_info"))

    # Only after every leg is stopped, so no HDF file is still open while it is
    # read. One leg's metadata failing must not cost another leg its own.
    for leg in legs:
        try:
            write_leg_metadata(leg)
        except Exception as exc:
            print(f"Could not write metadata for leg {leg.get('label')}: {exc}")


# -----------------------------------------------------------------------------
# Leg preparation
# -----------------------------------------------------------------------------


def prepare_legs(leg_specs):
    """Resolve devices and mode metadata for each leg.

    Everything that can fail without touching hardware fails here, before any motion.
    """
    if not leg_specs:
        raise ValueError("A parallel acquisition needs at least one detector leg.")

    legs = []

    for spec in leg_specs:
        leg = dict(spec)

        device = leg["device"]
        mode = leg["mode"]

        if device not in ACQ_MODES or mode not in ACQ_MODES[device]:
            raise ValueError(f"Invalid detector/mode combination: {device} / {mode}")

        # Before anything is resolved or moved. master_plan.validate_leg() makes
        # the same call, but multi_acq_series() can be driven directly.
        assert_parallel_safe(device, mode, where=f"Leg '{leg.get('label', device)}'")

        mode_info = ACQ_MODES[device][mode]

        for device_name in mode_info["required_devices"]:
            get_connected_device(device_name)

        leg["mode_info"] = mode_info
        leg["behaviour"] = multi_leg_behaviour(device, mode)
        leg["setup"] = mode_info["setup"]
        leg["det"] = get_connected_device(mode_info.get("hardware_device", device))
        leg["stop_timeout"] = float(leg.get("stop_timeout", DEFAULT_STOP_TIMEOUT))

        legs.append(leg)

    return legs


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------


def multi_acq_series(leg_specs, num_repeats=1, wait_time=0.0, cam_timeout=None):
    """Run several detectors in one beam window, repeated num_repeats times in lockstep.

    Args:
        leg_specs: One dict per detector. Required keys: device, mode, label,
            acq_time, num_frames, qmap_file, analysis_type. Optional: geometry
            (only needed when the detector is not at the mount
            device_position.yaml describes), motors, select_device,
            stop_timeout, workflow_name. The start and HDF-drain allowances are
            fixed: START_TIMEOUT and HDF_TIMEOUT.
        num_repeats: Repeats for the set. Every detector runs every repeat; the
            faster ones idle until the slowest finishes before the next starts.
        wait_time: Delay before each repeat.
        cam_timeout: Override the computed per-repeat cam allowance, in seconds.

    Shutter contract -- the reason this exists rather than calling the serial
    acquire functions once per detector:

        * The shutter opens before any detector is armed and closes only after
          EVERY cam reports done. showbeam()/blockbeam() appear exactly once per
          repeat, outside the per-leg loop -- no leg-level code touches the shutter.
        * A leg that stages before it reports -- the Rigaku -- must confirm it has
          STARTED before the wait for everyone to finish begins. Without that a
          staging Rigaku (detector_state Idle) reads as finished, and the shutter
          closes on an acquisition that has not begun.
        * DM submission waits for every HDF plugin, not just the leg that finished
          first.

    File naming, with the measurement number incremented ONCE for the whole set so
    every detector carries the same run number, and one folder per detector:

        A0012_Bismuth_a0010_f001000_rigaku3M_r00001
        A0012_Bismuth_a0010_f001000_eiger4M_r00001
        A0012_Bismuth_a0010_f001000_lambda2M_r00001
    """
    # Enforced here rather than at import: this module loads at startup, and schema
    # drift must not be able to break a session that only ever runs
    # single-detector acquisitions.
    _assert_override_paths()
    assert_multi_legs_sane()

    legs = prepare_legs(leg_specs)

    try:
        post_align()
        shutteroff()

        for leg in legs:
            if leg.get("select_device"):
                select_device(leg["device"])

        # The huber was already positioned by master_plan.setup_huber_for_multi(),
        # before this function was called. Nothing from here on may touch
        # huber.delta or huber.nu.
        for leg in legs:
            move_leg_motors(leg)

        workflow_proc_api, dmuser = dm_setup()

        # Once for the set, so every detector shares a run number and
        # measurement_num is only bumped once.
        folder_prefix = gen_folder_prefix()

        for leg in legs:
            leg["file_header"] = f"{folder_prefix}_f{int(leg['num_frames']):06d}_{leg['label']}"

        # Published to the `persistent:` block of state/run_state.yml
        # (expt.file_name), so a GUI or a shell script can see which set is
        # running. The legs' names agree right up to the detector label
        # (..._f003000_eiger4M_r00001 / ..._f003000_rigaku3M_r00001), so the shared
        # stem names the set without naming any one detector. It stops before the
        # label, so the _rNNNNN repeat suffix -- which comes after it -- is not part
        # of the stem. num_frames is per leg and nothing requires the legs to agree,
        # so fall back to the run prefix when they differ rather than publishing a
        # truncated frame count that matches no file.
        frame_counts = set()

        for leg in legs:
            frame_counts.add(int(leg["num_frames"]))

        if len(frame_counts) == 1:
            # One distinct value, so popping it just reads the only member.
            only_frame_count = frame_counts.pop()
            shared_name = f"{folder_prefix}_f{only_frame_count:06d}"
        else:
            shared_name = folder_prefix

        expt.file_name = shared_name

        if cam_timeout is not None:
            effective_cam_timeout = float(cam_timeout)
        else:
            effective_cam_timeout = cam_timeout_for(legs)

        # Legs whose "armed" and "busy" states are distinguishable, and so have a
        # confirmed-start phase of their own. Usually just the Rigaku.
        staging = [leg for leg in legs if leg["behaviour"].get("needs_start_confirm")]

        # One repeat, start to finish: stop any live/TV mode -> set every leg up
        # -> open the shutter -> arm every leg -> wait for the staging legs to
        # confirm they started -> wait for every cam -> close the shutter -> wait
        # for every HDF plugin -> write metadata and submit the DM job, one leg at
        # a time.
        for rep in range(int(num_repeats)):
            ttime.sleep(wait_time)

            sample_mesh_move()

            # Before setup, never after: a running camera silently swallows the
            # acquire_time write that setup is about to make.
            stop_live_mode(legs)

            for leg in legs:
                leg["file_name"] = f"{leg['file_header']}_r{rep + 1:05d}"
                leg["metadata_fname"] = leg["setup"](
                    acq_time=leg["acq_time"],
                    num_frames=int(leg["num_frames"]),
                    file_header=leg["file_header"],
                    file_name=leg["file_name"],
                )

            names = ", ".join(leg["file_name"] for leg in legs)
            print(f"\n{timestamp()}, Starting parallel measurement: {names}")

            # ---- the shutter window ------------------------------------------
            # Opens once, here, and closes once, below. Nothing between these two
            # points touches the shutter.
            showbeam()

            for leg in legs:
                leg["behaviour"]["arm"](leg["det"])
                print(f"{timestamp()}, Armed {leg['label']}")

            # Only the legs that stage before they report -- the Rigaku. See the
            # needs_start_confirm note above MULTI_LEGS for why the ADCore legs
            # neither need this nor can be asked. Fault-aware, so a leg that fails
            # to arm raises within one poll instead of holding the beam open for
            # the whole START_TIMEOUT.
            if staging:
                labels = ", ".join(leg["label"] for leg in staging)
                wait_all(staging, "started", START_TIMEOUT, f"{labels} to start acquiring", want=True)
                print(f"{timestamp()}, {labels} is acquiring")

            wait_all(
                legs,
                "cam_busy",
                effective_cam_timeout,
                "all detectors to finish acquiring",
                settle=ARM_SETTLE,
            )

            blockbeam()
            # ---- end of the shutter window -----------------------------------

            print(f"{timestamp()}, All detectors done acquiring, waiting for HDF plugins")

            wait_all(legs, "hdf_busy", HDF_TIMEOUT, "all HDF plugins to finish writing")

            # Every file is closed, so metadata and DM can run one leg at a time.
            for leg in legs:
                write_leg_metadata(leg)

                with swapped_registers(leg):
                    print(f"{timestamp()}, Submitting to DM, {leg['file_name']}")
                    dm_run_job(workflow_proc_api, dmuser, leg["file_name"])

            print(f"{timestamp()}, Complete measurement, repeat {rep + 1} of {int(num_repeats)}")

    except KeyboardInterrupt as exc:
        cleanup_multi(legs)
        raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).") from exc
    except Exception:
        cleanup_multi(legs)
        raise
