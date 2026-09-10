"""
Parallel two-detector acquisition for 8-ID.

Runs Eiger and Rigaku in the same beam window instead of back to back, so a pair of long
measurements costs roughly the slower one rather than the sum of both. The two detectors are
NOT frame-synced -- they just overlap. Each writes its own file, its own NeXus metadata, and
submits its own DM job.

This module changes nothing in the serial acquisition path. It reuses the setup half of
ad_acq.ACQ_MODES verbatim and replaces only the acquire half, which in the serial code
interleaves showbeam()/blockbeam() with a blocking wait and so cannot be run twice at once.

Normally you drive this through dual_master_plan.run_dual_measurement_info(), which reads
user_plans/dual_measurement_info.yaml. dual_acq_series() below is the lower-level entry point.

Shutter contract (the reason this module exists rather than calling the serial acquire
functions twice):

    * The shutter opens before either detector is armed, and closes only after BOTH cams
      report done. showbeam()/blockbeam() appear exactly once per repeat, outside the per-leg
      loop -- no leg-level code touches the shutter.
    * One leg is the shutter owner (the Rigaku). It is armed first, and the other legs are
      armed only after it confirms it is acquiring. The Rigaku's internal staging means the
      shutter has then been open for several seconds already, so no detector can ever
      integrate against a closed shutter.
    * DM submission waits for both HDF plugins, not just the leg that finished first.

Supported (device, mode) pairs are the keys of DUAL_LEGS. Anything else is rejected.
"""

import time as ttime
import warnings
from contextlib import contextmanager
from datetime import datetime

from apsbits.core.instrument_init import oregistry
from id8_common.plans.acquire.ad_acq import ACQ_MODES
from id8_common.plans.acquire.ad_acq import cleanup_acquisition
from id8_common.plans.acquire.ad_acq import gen_folder_prefix
from id8_common.plans.acquire.ad_acq import get_connected_device
from id8_common.plans.acquire.ad_acq import get_ophyd_object
from id8_common.plans.acquire.ad_acq import sample_mesh_move
from id8_common.plans.set.select_device import select_device
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import post_align
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.utils.dm_util import dm_run_job
from id8_common.utils.dm_util import dm_setup
from id8_common.utils.nexus_utils import create_nexus_format_metadata
from id8_common.utils.xpcs_schema import xpcs_schema

pv_registers = oregistry["pv_registers"]

POLL_INTERVAL = 0.1

# Wall-clock allowance for a cam to finish, as a multiple of the expected acquisition time
# plus a fixed floor. Generous on purpose: the point of a timeout here is to stop one hung
# detector from stranding the other one's data, not to police timing.
CAM_TIMEOUT_FACTOR = 3.0
CAM_TIMEOUT_FLOOR = 60.0

DEFAULT_START_TIMEOUT = 30.0
DEFAULT_HDF_TIMEOUT = 300.0

# huber axes that nothing inside an acquisition may drive.
#
# Both are positioned exactly once per measurement, before acquisition starts, by
# dual_master_plan.setup_huber_for_dual(). From that point on nothing may touch them: not a
# leg's `motors:` block, not the sample mesh. Both routes are refused at validation time and
# again here at run time.
FORBIDDEN_MOTORS = ("huber.delta", "huber.nu")


# =============================================================================
# NeXus metadata override paths
# =============================================================================
# create_runtime_metadata_dict() derives detector geometry from device_position.yaml via
# pv_registers.det_name, which can only describe one detector at a time. It applies its
# additional_metadata argument LAST, so passing the per-leg values through that existing hook
# overrides every geometry field without touching nexus_utils.py.
#
# The cost of that is these path strings duplicating ones that live in nexus_utils.py. If a
# path is renamed there, an override would silently stop applying and the file would get the
# wrong detector's geometry. _assert_override_paths() below turns that into an import-time
# error instead. update_schema_at_runtime() walks xpcs_schema, so xpcs_schema is what a path
# has to exist in.

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
            "dual_acq.py writes NeXus fields that no longer exist in xpcs_schema: "
            f"{missing}. The schema changed under this module -- update OVERRIDE_PATHS "
            "to match id8_common/utils/nexus_utils.py before running a dual acquisition."
        )


# Warn at import, enforce at acquisition. This module is imported by startup.py, so a hard
# failure here would take down the whole session -- including single-detector runs, which do
# not use any of this. Schema drift must not be able to do that.
if missing_override_paths():
    warnings.warn(
        f"dual_acq: NeXus override paths missing from xpcs_schema: {missing_override_paths()}. "
        "Dual acquisition will refuse to run until OVERRIDE_PATHS is updated. "
        "Single-detector acquisition is unaffected.",
        RuntimeWarning,
        stacklevel=2,
    )


# =============================================================================
# Per-leg arm and poll behaviour
# =============================================================================
# The serial acquire functions in ad_acq.py interleave showbeam()/blockbeam() with a blocking
# wait, so they cannot be used here -- whichever finished first would close the shutter on the
# other. These callables are the same detector handling with the shutter and the waiting
# lifted out, which is what lets the shutter contract live in one place.
#
# started vs cam_busy are genuinely different for the Rigaku: detector_state stays at 0 for
# seconds after cam.acquire.put(1) while the detector stages, so "not busy" and "not started
# yet" are indistinguishable without a separate confirmed-start phase.

DUAL_LEGS = {
    ("eiger4M", "Internal Series"): {
        "arm": lambda d: (d.hdf1.capture.put(1), d.cam.acquire.put(1)),
        "started": lambda d: d.cam.acquire.get() == 1,
        "cam_busy": lambda d: d.cam.acquire.get() == 1,
        "hdf_busy": lambda d: d.hdf1.capture.get() == 1,
    },
    ("rigaku3M_epics", "EPICS"): {
        "arm": lambda d: (d.hdf1.capture.put(1), d.cam.acquire.put(1)),
        "started": lambda d: d.cam.detector_state.get() == 1,
        # != 0 rather than == 1 so state 2 (readout) still counts as busy.
        "cam_busy": lambda d: d.cam.detector_state.get() != 0,
        "hdf_busy": lambda d: d.hdf1.capture.get() == 1,
    },
}


def dual_leg_behaviour(device, mode):
    """Look up the arm/poll callables for a (device, mode) pair."""
    key = (device, mode)

    if key not in DUAL_LEGS:
        supported = ", ".join(f"{d}/{m}" for d, m in sorted(DUAL_LEGS))
        raise ValueError(f"'{device}' in mode '{mode}' is not supported for dual acquisition. Supported: {supported}")

    return DUAL_LEGS[key]


# =============================================================================
# Helpers
# =============================================================================


def timestamp():
    """Log timestamp, matching the format ad_acq.py prints."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def resolve_value(value):
    """Resolve a geometry entry that is either a literal number or a dotted ophyd path.

    'huber.nu' reads that motor's current position; 5.0 is used as-is. Lets a leg report a
    live motor position in its metadata without needing a device_position.yaml entry.
    """
    if isinstance(value, str):
        return float(get_ophyd_object(value).position)

    return float(value)


def wait_until(predicate, timeout, description):
    """Poll predicate() until it is true. Raise TimeoutError if it never becomes true.

    The serial code uses unbounded `while True` loops. That is worse here: a hung detector
    would strand the other detector's data too, since the shutter and the HDF waits are shared.
    """
    deadline = ttime.time() + timeout

    while not predicate():
        if ttime.time() > deadline:
            raise TimeoutError(f"Timed out after {timeout:.1f} s waiting for {description}.")
        ttime.sleep(POLL_INTERVAL)


def wait_all(legs, key, timeout, description):
    """Wait until `key` reads false for EVERY leg.

    AND, not OR: the shutter cannot close and DM cannot start while any leg is still going.
    """
    deadline = ttime.time() + timeout

    while True:
        busy = [leg for leg in legs if leg["behaviour"][key](leg["det"])]

        if not busy:
            return

        if ttime.time() > deadline:
            names = ", ".join(leg["label"] for leg in busy)
            raise TimeoutError(f"Timed out after {timeout:.1f} s waiting for {description}. Still busy: {names}.")

        ttime.sleep(POLL_INTERVAL)


def cam_timeout_for(legs):
    """Allowance for the slowest leg's cam to finish."""
    expected = max(leg["acq_time"] * leg["num_frames"] for leg in legs)

    return max(CAM_TIMEOUT_FLOOR, CAM_TIMEOUT_FACTOR * expected + CAM_TIMEOUT_FLOOR)


def metadata_overrides(leg):
    """Per-leg NeXus values, passed through create_nexus_format_metadata's additional_metadata.

    Applied after everything create_runtime_metadata_dict() derived from device_position.yaml,
    so whatever is returned here wins.

    Only what the leg actually states is overridden. A leg sitting at its normal
    device_position.yaml mount needs no `geometry` block at all -- the values that file already
    holds are correct, and this returns just the two things that are per-leg by definition
    (which detector wrote the file, and which qmap goes with it). Override a geometry field
    only when this detector is somewhere device_position.yaml does not describe, e.g. an Eiger
    remounted on the huber arm.

    Lengths are millimetres in the YAML and metres in the file, matching what nexus_utils.py
    does with the same fields.
    """
    geometry = leg.get("geometry") or {}

    values = {
        "detector_name": leg["label"],
        "qmap_file": leg["qmap_file"],
    }

    # geometry key -> (override name, scale applied to the YAML value)
    direct = {
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

    for key, (name, scale) in direct.items():
        if key in geometry:
            values[name] = resolve_value(geometry[key]) * scale

    if "pixel_size" in geometry:
        pixel_size = resolve_value(geometry["pixel_size"])
        values["x_pixel_size"] = pixel_size
        values["y_pixel_size"] = pixel_size

    return {OVERRIDE_PATHS[name]: value for name, value in values.items()}


@contextmanager
def swapped_registers(leg):
    """Point the shared registers at one leg, then put them back.

    dm_run_job() reads det_name/qmap_file/analysis_type/workflow_name from pv_registers and
    has no override argument, and create_runtime_metadata_dict() needs det_name to resolve to
    a real device_position.yaml key before its additional_metadata override is applied. Both
    calls happen after the parallel window, one leg at a time, so the registers are never
    contended -- and they are restored afterwards so a dual run leaves no trace.
    """
    names = ["det_name", "qmap_file", "analysis_type", "workflow_name"]
    saved = {name: getattr(pv_registers, name).get() for name in names}

    try:
        pv_registers.det_name.put(leg["device"])
        pv_registers.qmap_file.put(leg["qmap_file"])
        pv_registers.analysis_type.put(leg["analysis_type"])

        if leg.get("workflow_name"):
            pv_registers.workflow_name.put(leg["workflow_name"])

        yield
    finally:
        for name in names:
            getattr(pv_registers, name).put(saved[name])


def move_leg_motors(leg):
    """Move whatever the leg's `motors` block names. Absent axes are left alone.

    huber.delta and huber.nu are refused here even though dual_master_plan already rejects
    them at validation time -- dual_acq_series() can be driven directly, bypassing that.
    """
    for dotted, position in (leg.get("motors") or {}).items():
        if dotted in FORBIDDEN_MOTORS:
            raise ValueError(
                f"Leg '{leg['label']}': '{dotted}' cannot be moved by a dual acquisition. "
                f"Both huber axes are positioned once before acquisition by "
                f"dual_master_plan.setup_huber_for_dual()."
            )

        motor = get_ophyd_object(dotted)
        print(f"Moving {dotted} to {position}")
        motor.move(float(position), wait=True)


def cleanup_dual(legs):
    """Close the shutter and stop every leg. Safe to call from a partially-started state."""
    try:
        blockbeam()
    except Exception as exc:
        print(f"Could not block beam during cleanup: {exc}")

    for leg in legs:
        cleanup_acquisition(leg.get("det"), leg.get("mode_info"))


# =============================================================================
# Leg preparation
# =============================================================================


def prepare_legs(leg_specs):
    """Resolve devices and mode metadata for each leg, and pick the shutter owner.

    Everything that can fail without touching hardware fails here, before any motion.
    """
    if not leg_specs:
        raise ValueError("A dual acquisition needs at least one detector leg.")

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
        leg["behaviour"] = dual_leg_behaviour(device, mode)
        leg["setup"] = mode_info["setup"]
        leg["det"] = get_connected_device(mode_info.get("hardware_device", device))
        leg["start_timeout"] = float(leg.get("start_timeout", DEFAULT_START_TIMEOUT))
        leg["hdf_timeout"] = float(leg.get("hdf_timeout", DEFAULT_HDF_TIMEOUT))

        legs.append(leg)

    owners = [leg for leg in legs if leg.get("shutter_owner")]

    if len(legs) == 1 and not owners:
        owners = legs
        legs[0]["shutter_owner"] = True

    if len(owners) != 1:
        labels = ", ".join(leg["label"] for leg in owners) or "none"
        raise ValueError(f"Exactly one leg must set shutter_owner. Got {len(owners)}: {labels}.")

    return legs, owners[0]


# =============================================================================
# Main entry point
# =============================================================================


def dual_acq_series(leg_specs, num_repeats=1, wait_time=0.0, cam_timeout=None):
    """Run several detectors in parallel, repeated num_repeats times in lockstep.

    Args:
        leg_specs: One dict per detector. Required keys: device, mode, label, acq_time,
            num_frames, qmap_file, analysis_type. Optional: geometry (only needed when the
            detector is not at the mount device_position.yaml describes), motors,
            select_device, shutter_owner, start_timeout, hdf_timeout, workflow_name.
        num_repeats: Repeats for the pair. Both detectors run every repeat; the faster one
            idles until the slower one finishes before the next repeat starts.
        wait_time: Delay before each repeat.
        cam_timeout: Override the computed per-repeat cam allowance, in seconds.

    File naming, with the measurement number incremented ONCE for the whole pair so both
    detectors carry the same run number:

        A0012_Bismuth_a0010_f001000_eiger4M_r00001
        A0012_Bismuth_a0010_f003000_rigaku3M_r00001
    """
    # Enforced here rather than at import: this module loads at startup, and schema drift must
    # not be able to break a session that only ever runs single-detector acquisitions.
    _assert_override_paths()

    legs, owner = prepare_legs(leg_specs)
    others = [leg for leg in legs if leg is not owner]

    try:
        post_align()
        shutteroff()

        for leg in legs:
            if leg.get("select_device"):
                select_device(leg["device"])

        # The huber was already positioned by dual_master_plan.setup_huber_for_dual(), before
        # this function was called. Nothing from here on may touch huber.delta or huber.nu.
        for leg in legs:
            move_leg_motors(leg)

        workflow_proc_api, dmuser = dm_setup()

        # Once for the pair, so both detectors share a run number and measurement_num is
        # only bumped once.
        folder_prefix = gen_folder_prefix()

        for leg in legs:
            leg["file_header"] = f"{folder_prefix}_f{int(leg['num_frames']):06d}_{leg['label']}"

        effective_cam_timeout = float(cam_timeout) if cam_timeout is not None else cam_timeout_for(legs)

        for rep in range(int(num_repeats)):
            ttime.sleep(wait_time)

            sample_mesh_move()

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

            # ---- the shutter window -------------------------------------------------
            # Opens once, here, and closes once, below. Nothing between these two points
            # touches the shutter.
            showbeam()

            owner["behaviour"]["arm"](owner["det"])
            print(f"{timestamp()}, Armed shutter owner {owner['label']}, waiting for it to start")

            wait_until(
                lambda: owner["behaviour"]["started"](owner["det"]),
                timeout=owner["start_timeout"],
                description=f"{owner['label']} to start acquiring",
            )
            print(f"{timestamp()}, {owner['label']} is acquiring")

            # The owner is confirmed running, so the shutter has been open for a while.
            # Only now is it safe to start everyone else.
            for leg in others:
                leg["behaviour"]["arm"](leg["det"])
                print(f"{timestamp()}, Started {leg['label']}")

            wait_all(legs, "cam_busy", effective_cam_timeout, "all detectors to finish acquiring")

            blockbeam()
            # ---- end of the shutter window ------------------------------------------

            print(f"{timestamp()}, All detectors done acquiring, waiting for HDF plugins")

            hdf_timeout = max(leg["hdf_timeout"] for leg in legs)
            wait_all(legs, "hdf_busy", hdf_timeout, "all HDF plugins to finish writing")

            # Both files are closed, so metadata and DM can run one leg at a time.
            for leg in legs:
                print(f"{timestamp()}, Writing metadata, {leg['file_name']}")

                with swapped_registers(leg):
                    create_nexus_format_metadata(
                        leg["metadata_fname"],
                        det=leg["det"],
                        additional_metadata=metadata_overrides(leg),
                    )

                    print(f"{timestamp()}, Submitting to DM, {leg['file_name']}")
                    dm_run_job(workflow_proc_api, dmuser, leg["file_name"])

            print(f"{timestamp()}, Complete measurement, repeat {rep + 1} of {int(num_repeats)}")

    except KeyboardInterrupt as exc:
        cleanup_dual(legs)
        raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).") from exc
    except Exception:
        cleanup_dual(legs)
        raise
