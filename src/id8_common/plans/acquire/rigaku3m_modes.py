"""
Rigaku3M mode definitions: setup/acquire functions plus RIGAKU3M_BIN_MODES,
RIGAKU3M_FTF_MODES and RIGAKU3M_EPICS_MODES, the tables ad_acq.py assembles
into ACQ_MODES. Nothing here runs on its own: det_acq_series() in ad_acq.py
calls these, and trio_acq_series() runs setup_rigaku_epics() for the
("rigaku3M_epics", "EPICS") trio leg.

The three tables are separate ACQ_MODES *detector* keys rather than one table,
because BIN and FTF offer the same three mode names (ZDT2bit / ZDT4bit /
ZDT8bit) and ACQ_MODES is looked up as ACQ_MODES[det_name][det_mode]. The
det_name in the run state (expt.det_name, set per-protocol from
measurement_info.yaml) is what picks the output format:

    rigaku3M ("bin", RIGAKU3M_BIN_MODES)
        ZDT2bit / ZDT4bit / ZDT8bit -- zero-deadtime sparsified output at
        2/4/8 bit resolution, written directly to a .bin file by the
        detector's own fast-file mechanism (not areaDetector HDF5).
        output_control = "Sparsified".

    rigaku3M_ftf ("fast transfer", RIGAKU3M_FTF_MODES)
        ZDT2bit / ZDT4bit / ZDT8bit -- as the bin modes above but with
        output_control = "FastTransfer - HDF5", so the detector's fast-file
        mechanism emits HDF5 instead of the sparsified format, written as
        .h5 rather than .bin. Same trigger, timing and resolution.

    rigaku3M_epics (RIGAKU3M_EPICS_MODES)
        EPICS -- standard areaDetector-style acquisition through the HDF1
        plugin, for when neither fast-file format is wanted.

Both fast-file families share acquire_rigaku_zdt(); only EPICS needs its own
acquire, because it also has to arm and drain the HDF1 plugin.
"""

import os
import time as ttime

from id8_common.plans.acquire.acq_helpers import get_connected_device
from id8_common.plans.acquire.acq_helpers import get_rigaku_file_path
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import showbeam

# =============================================================================
# Detector setup functions
#
# All of these configure only -- nothing starts until the matching acquire_*
# function runs -- and each returns the NeXus metadata file path, which the
# caller (det_acq_series) fills in once the run is over. The six fast-file
# setups are one-line wrappers around _setup_rigaku_fast_file(); only the
# EPICS one has a body of its own, because it is the only mode that writes
# through the areaDetector HDF1 plugin.
# =============================================================================

def _setup_rigaku_fast_file(acq_time, num_frames, file_header, file_name,
                            output_resolution, output_control, file_extension):
    """Shared body of the six fast-file (bin/FTF) setup functions.

    The bin and FTF modes differ in ``output_control`` and in the extension
    they ask the IOC to write, and the three resolutions only in
    ``output_resolution``, so the named setup functions below are thin wrappers
    that pin those three strings. The PV write order is unchanged from the
    original per-resolution functions -- it has not been established which
    orderings the IOC tolerates, so don't reshuffle it without testing on the
    detector.

    ``file_extension`` is ``"bin"`` for the sparsified modes and ``"h5"`` for
    fast transfer, whose output really is HDF5 (verified by its \x89HDF magic
    on 2026-09-03, back when it was still being written under a .bin name).
    The IOC appends a per-module suffix to whatever it is given, so the files
    land as ``<file_name>.<ext>.000`` .. ``.005``.
    """
    rigaku3M = get_connected_device("rigaku3M")
    softglue = get_connected_device("softglue")

    softglue.enable_rigaku.put('1')

    rigaku3M.cam.trigger_mode.put('Start with Trigger')

    file_path, full_path = get_rigaku_file_path(file_header, file_name)

    rigaku3M.cam.acquire_time.put(acq_time)
    rigaku3M.cam.acquire_period.put(acq_time)

    rigaku3M.cam.fast_file_name.put(f"{file_name}.{file_extension}")
    rigaku3M.cam.fast_file_path.put(file_path)

    rigaku3M.cam.num_images.put(num_frames)
    rigaku3M.cam.output_control.put(output_control)
    rigaku3M.cam.image_mode.put("Zero-Deadtime")
    rigaku3M.cam.output_resolution.put(output_resolution)

    os.makedirs(full_path, mode=0o770, exist_ok=True)

    metadata_fname = f"{full_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_rigaku_zdt2bit(acq_time, num_frames, file_header, file_name):
    return _setup_rigaku_fast_file(
        acq_time, num_frames, file_header, file_name,
        output_resolution="2 Bit", output_control="Sparsified",
        file_extension="bin",
    )


def setup_rigaku_zdt4bit(acq_time, num_frames, file_header, file_name):
    return _setup_rigaku_fast_file(
        acq_time, num_frames, file_header, file_name,
        output_resolution="4 Bit", output_control="Sparsified",
        file_extension="bin",
    )


def setup_rigaku_zdt8bit(acq_time, num_frames, file_header, file_name):
    return _setup_rigaku_fast_file(
        acq_time, num_frames, file_header, file_name,
        output_resolution="8 Bit", output_control="Sparsified",
        file_extension="bin",
    )


def setup_rigaku_ftf2bit(acq_time, num_frames, file_header, file_name):
    return _setup_rigaku_fast_file(
        acq_time, num_frames, file_header, file_name,
        output_resolution="2 Bit", output_control="FastTransfer - HDF5",
        file_extension="h5",
    )


def setup_rigaku_ftf4bit(acq_time, num_frames, file_header, file_name):
    return _setup_rigaku_fast_file(
        acq_time, num_frames, file_header, file_name,
        output_resolution="4 Bit", output_control="FastTransfer - HDF5",
        file_extension="h5",
    )


def setup_rigaku_ftf8bit(acq_time, num_frames, file_header, file_name):
    return _setup_rigaku_fast_file(
        acq_time, num_frames, file_header, file_name,
        output_resolution="8 Bit", output_control="FastTransfer - HDF5",
        file_extension="h5",
    )


def setup_rigaku_epics(acq_time, num_frames, file_header, file_name):
    """Configure the Rigaku for plain areaDetector acquisition. Return the metadata path.

    The one mode that bypasses the detector's own fast-file writer: full 16-bit
    frames go through the areaDetector pipeline and out via the HDF1 plugin, so
    unlike the ZDT modes this one has an hdf1 block to fill in (and an hdf1 to
    arm and drain in acquire_rigaku_epics).
    """
    rigaku3M = get_connected_device("rigaku3M")
    softglue = get_connected_device("softglue")

    softglue.enable_rigaku.put('1')

    rigaku3M.cam.trigger_mode.put('Start with Trigger')

    # Only the second value is wanted here. The first is the IOC-relative
    # path, which only the fast-file modes use (cam.fast_file_path); the HDF1
    # plugin is given the absolute one instead.
    _, full_path = get_rigaku_file_path(file_header, file_name)

    rigaku3M.cam.acquire_time.put(acq_time)
    rigaku3M.cam.acquire_period.put(acq_time)

    rigaku3M.cam.num_images.put(num_frames)
    rigaku3M.cam.output_control.put("areaDetector")
    rigaku3M.cam.image_mode.put("Standard")
    rigaku3M.cam.output_resolution.put("16 Bit")

    rigaku3M.hdf1.file_name.put(file_name)
    rigaku3M.hdf1.file_path.put(full_path)
    rigaku3M.hdf1.num_capture.put(num_frames)

    os.makedirs(full_path, mode=0o770, exist_ok=True)

    metadata_fname = f"{full_path}/{file_name}_metadata.hdf"

    return metadata_fname


# =============================================================================
# Detector acquire functions
# =============================================================================

# ADCore ADStatus enum on 8idRigaku3m:cam1:DetectorState_RBV. Acquire/Readout/
# Saving/Waiting are all normal mid-acquisition states; the last three are
# terminal and will never become Idle on their own, so a wait loop that only
# tests "!= Idle" spins forever on them.
STATE_IDLE = 0
STATE_ACQUIRE = 1
STATE_FAILED = {6: "Error", 9: "Disconnected", 10: "Aborted"}

# Idle is the ONLY completion signal available for the fast-file modes -- see
# the note in acquire_rigaku_zdt() about why the frame counter cannot be used.

# Generous: covers the file-writing tail, which for a large sparsified run is
# far longer than the acquisition itself. Only there so a wedged detector
# cannot hold the shutter open indefinitely.
ACQUIRE_TIMEOUT = 3600.0
ARM_TIMEOUT = 60.0


def _wait_for_state(cam, want, timeout, what):
    """Block until detector_state reaches `want`, or raise.

    Raises on a terminal failure state or on timeout rather than looping
    forever -- the loop this replaces had no exit for either, so a detector
    that errored mid-acquisition hung the plan with the shutter still open.
    """
    deadline = ttime.time() + timeout
    while True:
        state = cam.detector_state.get()

        if state == want:
            return

        if state in STATE_FAILED:
            raise RuntimeError(
                f"rigaku3M went to {STATE_FAILED[state]} while {what}. "
                f"Detector message: {cam.status_message.get()!r}"
            )

        if ttime.time() > deadline:
            raise RuntimeError(
                f"rigaku3M timed out after {timeout:g}s while {what} "
                f"(detector_state={state})."
            )

        ttime.sleep(0.1)


def acquire_rigaku_zdt():
    """Run one zero-deadtime acquisition. Serves all six bin and FTF modes.

    Nothing here depends on the resolution or the output format -- those were
    already pushed to the detector by the setup function -- so one acquire
    covers both fast-file families.
    """
    rigaku3M = get_connected_device("rigaku3M")

    # blockbeam() in a finally: it used to sit after two unbounded wait loops,
    # so any detector fault left the beam on the sample until someone noticed.
    try:
        showbeam()
        ttime.sleep(0.1)

        # Start. The detector writes straight to disk itself in these modes,
        # so there is no HDF1 plugin to arm or drain.
        rigaku3M.cam.acquire.put(1)

        # Two stages: first confirm it actually started, then wait out the
        # acquisition and the file-writing tail that follows it.
        _wait_for_state(rigaku3M.cam, STATE_ACQUIRE, ARM_TIMEOUT, "arming")
        _wait_for_state(rigaku3M.cam, STATE_IDLE, ACQUIRE_TIMEOUT, "acquiring")
    finally:
        blockbeam()

    # Deliberately no frame-count check here. cam.num_images_counter is the
    # areaDetector array counter, and Zero-Deadtime output never goes through
    # the areaDetector pipeline -- the detector's own fast-file writer takes it
    # straight to disk -- so the counter barely moves regardless of how the run
    # went (measured 2026-09-03: it read 3 after a healthy 100000-frame ZDT
    # acquisition). Returning to Idle is the completion signal for these modes;
    # don't reintroduce a counter comparison, it aborts good measurements.


def acquire_rigaku_epics():
    """Run one areaDetector acquisition, writing through the HDF1 plugin.

    Same waits as acquire_rigaku_zdt(), plus the drain at the end: this is the
    only Rigaku mode with an HDF plugin to finish writing.
    """
    rigaku3M = get_connected_device("rigaku3M")

    try:
        showbeam()
        ttime.sleep(0.1)

        # Arm the file writer before the detector, so the plugin is already
        # capturing when the first frame arrives.
        rigaku3M.hdf1.capture.put(1)
        rigaku3M.cam.acquire.put(1)

        _wait_for_state(rigaku3M.cam, STATE_ACQUIRE, ARM_TIMEOUT, "arming")
        _wait_for_state(rigaku3M.cam, STATE_IDLE, ACQUIRE_TIMEOUT, "acquiring")
    finally:
        blockbeam()

    # Drain: capture stays at 1 until the plugin has written every frame.
    # Bounded like the state waits above and for the same reason -- a stuck
    # HDF1 plugin should surface as an error, not as a plan that never returns.
    deadline = ttime.time() + ACQUIRE_TIMEOUT
    while rigaku3M.hdf1.capture.get() == 1:
        if ttime.time() > deadline:
            raise RuntimeError(
                f"rigaku3M hdf1 still capturing after {ACQUIRE_TIMEOUT:g}s "
                f"({rigaku3M.hdf1.num_captured.get()} of "
                f"{rigaku3M.hdf1.num_capture.get()} frames written)."
            )
        ttime.sleep(0.1)


# =============================================================================
# Mode tables
#
# ad_acq.py mounts these as ACQ_MODES["rigaku3M"], ["rigaku3M_ftf"] and
# ["rigaku3M_epics"], and det_acq_series() looks up
# ACQ_MODES[expt.det_name][expt.det_mode] to find out how to run one
# measurement. What the keys mean:
#
#   setup             function that configures the hardware
#   acquire           function that then runs it
#   needs_acq_period  pass acq_period to setup as well? False throughout here:
#                     every Rigaku setup is given only acq_time, num_frames,
#                     file_header and file_name
#   required_devices  all must connect before the run starts
#   hardware_device   the oregistry device name, needed whenever the
#                     ACQ_MODES key is not itself a device ("rigaku3M_ftf",
#                     "rigaku3M_epics"); it defaults to the ACQ_MODES key
#   min_acq_time      floor validators.py enforces on acq_time, in seconds.
#                     The ZDT figures scale with bit depth -- twice the bits,
#                     twice the readout time per frame
#
# Adding a mode means adding a row here plus its setup_/acquire_ pair above.
# =============================================================================

RIGAKU3M_BIN_MODES = {
    "ZDT2bit": {
        "setup": setup_rigaku_zdt2bit,
        "acquire": acquire_rigaku_zdt,
        "needs_acq_period": False,
        "required_devices": ["rigaku3M"],
        "min_acq_time": 20e-6,
    },
    "ZDT4bit": {
        "setup": setup_rigaku_zdt4bit,
        "acquire": acquire_rigaku_zdt,
        "needs_acq_period": False,
        "required_devices": ["rigaku3M"],
        "min_acq_time": 40e-6,
    },
    "ZDT8bit": {
        "setup": setup_rigaku_zdt8bit,
        "acquire": acquire_rigaku_zdt,
        "needs_acq_period": False,
        "required_devices": ["rigaku3M"],
        "min_acq_time": 80e-6,
    },
}

# Same mode names and timing limits as RIGAKU3M_BIN_MODES -- the tables are
# kept apart because they are two ACQ_MODES detector keys, and "rigaku3M_ftf"
# is not an oregistry device name, so these need hardware_device.
RIGAKU3M_FTF_MODES = {
    "ZDT2bit": {
        "setup": setup_rigaku_ftf2bit,
        "acquire": acquire_rigaku_zdt,
        "needs_acq_period": False,
        "required_devices": ["rigaku3M"],
        "hardware_device": "rigaku3M",
        "min_acq_time": 20e-6,
    },
    "ZDT4bit": {
        "setup": setup_rigaku_ftf4bit,
        "acquire": acquire_rigaku_zdt,
        "needs_acq_period": False,
        "required_devices": ["rigaku3M"],
        "hardware_device": "rigaku3M",
        "min_acq_time": 40e-6,
    },
    "ZDT8bit": {
        "setup": setup_rigaku_ftf8bit,
        "acquire": acquire_rigaku_zdt,
        "needs_acq_period": False,
        "required_devices": ["rigaku3M"],
        "hardware_device": "rigaku3M",
        "min_acq_time": 80e-6,
    },
}

RIGAKU3M_EPICS_MODES = {
    "EPICS": {
        "setup": setup_rigaku_epics,
        "acquire": acquire_rigaku_epics,
        "needs_acq_period": False,
        "required_devices": ["rigaku3M"],
        "hardware_device": "rigaku3M",
        "min_acq_time": 0.01,
    },
}
