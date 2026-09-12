"""
Consolidated acquisition script for 8-ID detectors.

Per-detector setup/acquire functions and their mode tables are just data and
callables until this module assembles them -- they live in sibling
*_modes.py files (see each file's docstring for what each mode does and when
to use it):
    eiger4m_modes.py  : eiger4M       -- "Internal Series", "Internal Enable",
                                          "External Series", "External Enable"
    lambda2m_modes.py : lambda2M      -- "Internal", "External"
    rigaku3m_modes.py : rigaku3M      -- "ZDT2bit", "ZDT4bit", "ZDT8bit"
                        rigaku3M_ftf   -- "ZDT2bit", "ZDT4bit", "ZDT8bit"
                        rigaku3M_epics -- "EPICS"

This module assembles those into ACQ_MODES and provides the two acquisition
loops that run against it:

    det_acq_series()    one detector at a time. The setup function arms the
                        hardware and the mode's own acquire function opens the
                        shutter, waits, and closes it again.

    multi_acq_series()  several detectors in ONE beam window. Reuses the setup
                        half of ACQ_MODES verbatim and replaces the acquire
                        half, which cannot be run twice at once because it
                        owns the shutter. Was trio_acq_rigaku3m_eiger4m_lambda2m.py
                        until 2026-09-11; see the section header further down
                        for why it folded in here.

Both are driven from master_plan.py, which picks between them per protocol:
a protocol with a scalar ``detector:`` takes the first, one with a
``detectors:`` list takes the second.

scan_8id.py and scan_8id_dev.py do `from id8_common.plans.acquire.ad_acq
import *` and call setup_eiger_internal(...) etc. directly, so this module
re-exports everything from the per-detector files the same way -- no __all__
needed: unlike named imports, ruff's --fix leaves `import *` alone (it can't
tell what's "unused" behind a star import), so nothing here is at risk of
being silently stripped on the next pre-commit run. The per-detector files
have no __all__ either, so all of their public names come through; acq_helpers
does have one, and it is what decides which of ITS names (gen_folder_prefix,
get_ophyd_object, sample_mesh_move, ...) reach here and, through startup.py,
the prompt.
"""

import os
import time as ttime
import warnings
from contextlib import contextmanager
from datetime import datetime

from id8_common.expt_config import expt
from id8_common.plans.acquire import acq_helpers
from id8_common.plans.acquire.acq_helpers import *
from id8_common.plans.acquire.acq_wait import DEFAULT_STOP_TIMEOUT
from id8_common.plans.acquire.acq_wait import STATE_IDLE
from id8_common.plans.acquire.acq_wait import acq_timeout
from id8_common.plans.acquire.acq_wait import cam_fault
from id8_common.plans.acquire.acq_wait import stop_acquiring
from id8_common.plans.acquire.acq_wait import wait_until
from id8_common.plans.acquire.eiger4m_modes import *
from id8_common.plans.acquire.lambda2m_modes import *
from id8_common.plans.acquire.rigaku3m_modes import *
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

# =============================================================================
# Mode table
# =============================================================================

ACQ_MODES = {
    "eiger4M": EIGER4M_MODES,
    "lambda2M": LAMBDA2M_MODES,
    "rigaku3M": RIGAKU3M_BIN_MODES,
    "rigaku3M_ftf": RIGAKU3M_FTF_MODES,
    "rigaku3M_epics": RIGAKU3M_EPICS_MODES,
}


# =============================================================================
# Cleanup helper
# =============================================================================

def cleanup_acquisition(det=None, mode_info=None, metadata_fname=None):
    """
    Close the shutter, stop softglue and the detector, write metadata.

    metadata_fname:
        NeXus metadata path of the measurement that was still in flight. The partial
        dataset left behind by the abort keeps its metadata this way. Pass None (the
        default) once the metadata has already been written, so an abort later in the
        same repeat does not overwrite a complete file.
    """
    # First, and in its own try, so no later failure can leave the beam on the sample.
    try:
        blockbeam()
    except Exception as e:
        print(f"Could not block beam during cleanup: {e}")

    if mode_info is not None and "softglue" in mode_info.get("required_devices", []):
        try:
            softglue = get_connected_device("softglue")
            softglue.stop_pulses.put("1!")
        except Exception:
            pass

    if det is not None:
        try:
            det.cam.acquire.put(0)
        except Exception:
            pass
        if hasattr(det, "hdf1"):
            try:
                det.hdf1.capture.put(0)
            except Exception:
                pass

    # Last, so the detector is stopped and the data file is closed before it is read.
    if metadata_fname is not None and det is not None:
        try:
            print(f"Writing metadata for interrupted measurement: {metadata_fname}")
            create_nexus_format_metadata(metadata_fname, det=det)
        except Exception as e:
            print(f"Could not write metadata {metadata_fname}: {e}")


# =============================================================================
# Main user-facing acquisition function
# =============================================================================

def det_acq_series(wait_time=0, hooks=None):
    """Run repeated detector acquisitions from the current run state.

    Parameters come from ``expt`` (id8_common/expt_config.py): the protocol
    values are loaded from measurement_info.yaml by
    master_plan.run_measurement(), the per-experiment ones from
    configs/experiment.yml, and measurement_num from 8ideSoft:Reg1.
    ``expt.file_name`` is an output, not an input -- this function builds it
    per repeat and writes it back. See docs/running-measurements.md for
    sample-motion behavior, hooks, and file-naming rules.

    To drive this directly, without master_plan, set the run state first::

        expt.set_measurement({"detector": "eiger4M", "mode": "Internal Series",
                              "acq_time": 0.01, "acq_period": 0.01,
                              "num_frames": 100, "num_repeats": 1,
                              "sample_move": "no", "qmap_file": "..."},
                             {"header": "A", "sample_name": "Test"})
        det_acq_series()

    Error handling, which matters to anything calling this in a loop: Ctrl+C is
    turned into a RuntimeError (so an upstream ``except KeyboardInterrupt``
    will NOT catch it), while any other exception is cleaned up after, printed,
    and swallowed -- this function then returns normally. A caller therefore
    cannot tell a failed measurement from a successful one by exception alone.
    Either way cleanup_acquisition() runs first and attempts to close the
    shutter before anything else.
    """
    det = None
    mode_info = None
    metadata_fname = None
    try:
        # The [:] overwrites the contents of the existing list instead of
        # rebinding the name. acq_helpers is star-imported here and, through
        # startup.py, into the interactive namespace, so several names are
        # already bound to this one list object; a plain assignment would
        # replace only acq_helpers' own attribute and leave the rest pointing
        # at the old, empty list.
        acq_helpers.active_hooks[:] = load_hooks(hooks)

        post_align()
        shutteroff()

        workflowProcApi, dmuser = dm_setup()

        detector = expt.det_name
        mode = expt.det_mode

        acq_time = expt.acq_time
        acq_period = expt.acq_period
        num_frames = expt.num_frames
        num_reps = expt.num_repeats
        num_segments = expt.num_segments
        trigger_period = expt.trigger_period

        mode_info = ACQ_MODES[detector][mode]

        for device_name in mode_info["required_devices"]:
            get_connected_device(device_name)

        # Most ACQ_MODES keys are themselves device names, but a few are aliases
        # for different ways of driving one detector -- rigaku3M_ftf and
        # rigaku3M_epics both drive the rigaku3M device -- and those table
        # entries carry a "hardware_device" naming the real one.
        det = get_connected_device(mode_info.get("hardware_device", detector))
        setup_func = mode_info["setup"]
        acquire_func = mode_info["acquire"]

        folder_prefix = gen_folder_prefix()
        file_header = f"{folder_prefix}_f{num_frames:06d}"

        for rep in range(num_reps):
            ttime.sleep(wait_time)

            sample_mesh_move()

            file_name = f"{file_header}_r{rep + 1:05d}"

            # Persistent session state: written to the `persistent:` block of
            # state/run_state.yml, so a GUI or a shell script can see which
            # measurement is running without a live session. It used to be
            # 8ideSoft:StrReg8 -- that Component still exists on pv_registers,
            # but nothing writes it any more, so don't caget it.
            # Bare name, no extension and no path -- the folder is rebuildable
            # from mount_point/cycle_name/experiment_name/use_subfolder, and
            # each detector appends its own suffix (.h5, .bin, _metadata.hdf).
            expt.file_name = file_name

            # Only the parameters a mode declares are passed, so a setup
            # function takes exactly the arguments it uses. Add a
            # "needs_<name>" flag to a mode table entry to extend this.
            setup_kwargs = {
                "acq_time": acq_time,
                "num_frames": num_frames,
                "file_header": file_header,
                "file_name": file_name,
            }

            # needs_acq_period is spelled with plain [...] because every mode
            # entry must define it; the two below use .get(..., False) because
            # they are optional and most entries leave them out.
            if mode_info["needs_acq_period"]:
                setup_kwargs["acq_period"] = acq_period

            if mode_info.get("needs_num_segments", False):
                setup_kwargs["num_segments"] = num_segments

            if mode_info.get("needs_trigger_period", False):
                setup_kwargs["trigger_period"] = trigger_period

            metadata_fname = setup_func(**setup_kwargs)

            time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"\n{time_now}, Starting measurement {file_name}")

            acquire_func()

            time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{time_now}, Complete measurement {file_name}")

            print(f"{time_now}, Writing metadata, {file_name}")

            # One writer, no env-var switch: utils/nexus_writer.py, which is
            # Miaoqi Chu's nexus_xpcs_aps for the schema and the HDF5 write,
            # over our utils/nexus_runtime.py for the EPICS-signal mapping he
            # does not model. The abort path above uses the same one, so an
            # aborted run and a good run produce the same shape of file.
            create_nexus_format_metadata(metadata_fname, det=det)

            # Complete, so an abort during DM submission must not rewrite it.
            metadata_fname = None

            print(f"{time_now}, Submitting to DM, {file_name}")
            dm_run_job(workflowProcApi, dmuser, file_name)

    except KeyboardInterrupt:
        cleanup_acquisition(det, mode_info, metadata_fname)
        raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).")
    except Exception as e:
        # DELIBERATE, do not "fix" this into a raise. A measurement_info.yaml run
        # is often an overnight queue of dozens of measurements; one that fails --
        # a detector hiccup, a full disk, a DM submission that bounces -- must not
        # take the rest of the night with it. So the error is reported and this
        # measurement is abandoned, and master_plan.run_measurement() goes on to
        # the next one.
        #
        # The cost, which is real: master_plan cannot tell a failed measurement
        # from a successful one, so a run can "complete" with gaps in it. Read the
        # printed errors before trusting a night's output. Decided by the owner
        # 2026-09-07 with that tradeoff understood.
        cleanup_acquisition(det, mode_info, metadata_fname)
        print(f"Error occurred during measurement: {e}")
    finally:
        pass

# =============================================================================
# Parallel multi-detector acquisition
# =============================================================================
# Runs several detectors in the same beam window instead of back to back, so a
# set of long measurements costs roughly the slowest one rather than the sum of
# all of them. The detectors are NOT frame-synced -- they merely overlap. Each
# writes its own file into its own folder, its own NeXus metadata, and submits
# its own DM job.
#
# This lived in trio_acq_rigaku3m_eiger4m_lambda2m.py until 2026-09-11, next to a
# trio_master_plan_* front end that duplicated master_plan.py's validation and run
# expansion. It folded in here once setup_rigaku_epics() stopped driving the fast
# shutter (softglue.enable_rigaku '1' + trigger_mode 'Start with Trigger'). That
# was what forced the split: with the Rigaku gating the beam through softglue,
# the shutter belonged to one nominated leg, every other leg had to be armed
# after that leg confirmed it was running, and none of it resembled the serial
# path. Now EVERY supported mode is internally timed and none of them touches
# softglue, so a parallel measurement is the serial one with the arm and the wait
# done over a list:
#
#     open the shutter  ->  arm every detector  ->  wait for the Rigaku to
#       confirm it started  ->  wait for every cam  ->  close the shutter
#       ->  wait for every HDF plugin  ->  metadata + DM, one leg at a time
#
# There are NO THREADS here and nothing runs out of order. Arming a detector is
# an EPICS put that tells the hardware to start and returns immediately, so the
# loop below arms one leg, then the next, then the next, and from that moment
# every detector is integrating while the Python side sits in a single polling
# loop (wait_all). Read it top to bottom like any other plan: the only thing
# happening concurrently is the hardware.
#
# Supported (device, mode) pairs are the keys of MULTI_LEGS. Anything else is
# rejected.

# Wall-clock allowance for a cam to finish, as a multiple of the expected
# acquisition time plus a fixed floor. Generous on purpose: the point of a
# timeout here is to stop one hung detector from stranding the others' data, not
# to police timing.
CAM_TIMEOUT_FACTOR = 3.0
CAM_TIMEOUT_FLOOR = 60.0

#: Seconds to give a detector to confirm it has started, once armed. The Rigaku
#: is the reason this is not instant: it stages internally for a few seconds
#: after cam.acquire.put(1) before it reports anything.
DEFAULT_START_TIMEOUT = 30.0

#: Seconds to give the HDF plugins to flush after the cams stop.
DEFAULT_HDF_TIMEOUT = 300.0

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
        f"ad_acq: NeXus override paths missing from xpcs_schema: {missing_override_paths()}. "
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

        mode_info = ACQ_MODES[device][mode]

        for device_name in mode_info["required_devices"]:
            get_connected_device(device_name)

        leg["mode_info"] = mode_info
        leg["behaviour"] = multi_leg_behaviour(device, mode)
        leg["setup"] = mode_info["setup"]
        leg["det"] = get_connected_device(mode_info.get("hardware_device", device))
        leg["start_timeout"] = float(leg.get("start_timeout", DEFAULT_START_TIMEOUT))
        leg["stop_timeout"] = float(leg.get("stop_timeout", DEFAULT_STOP_TIMEOUT))
        leg["hdf_timeout"] = float(leg.get("hdf_timeout", DEFAULT_HDF_TIMEOUT))

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
            start_timeout, stop_timeout, hdf_timeout, workflow_name.
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

        start_timeout = max((leg["start_timeout"] for leg in staging), default=0.0)
        hdf_timeout = max(leg["hdf_timeout"] for leg in legs)

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
            # the whole start_timeout.
            if staging:
                labels = ", ".join(leg["label"] for leg in staging)
                wait_all(staging, "started", start_timeout, f"{labels} to start acquiring", want=True)
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

            wait_all(legs, "hdf_busy", hdf_timeout, "all HDF plugins to finish writing")

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


# =============================================================================
# Example use cases -- see docs/running-measurements.md for
# run-state-driven det_acq_series() examples (Eiger/Lambda/Rigaku, with and
# without sample motion).
# =============================================================================
