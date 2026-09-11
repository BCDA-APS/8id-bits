"""
Eiger4M mode definitions: setup/acquire functions plus EIGER4M_MODES, the
table ad_acq.py assembles into ACQ_MODES. Nothing here runs on its own. The
callers are det_acq_series() in ad_acq.py; dual_acq_series(), for the
("eiger4M", "Internal Series") dual leg; and scan_8id.py / scan_8id_dev.py,
which call setup_eiger_internal() directly.

Modes (matching the TriggerMode enum on 8idEiger4M:cam1):

    Internal Series (trigger_mode="Internal Series")
        No external hardware. One cam.acquire arms the whole series: the
        detector internally times and captures num_images frames back-to-back
        at acquire_time/acquire_period. Use for a plain, unsynchronized
        fixed-count exposure series.

    Internal Enable (trigger_mode="Internal Enable")
        No external hardware. manual_trigger is set to "Enable", so the
        detector waits for a software trigger (cam.special_trigger_button --
        the same signal the pre-armed step-scans in scan_8id.py/ophyd_scan.py
        use) before each exposure. num_images is ignored; the acquisition
        completes once num_triggers software triggers have been issued.
        acquire_time/acquire_period (on the cam) set both the exposure and
        the pacing this code uses between triggers. Use for software-paced
        single exposures with no external pulse generator.

    External Series (trigger_mode="External Series")
        Uses the same softglue pulse train as External Enable, but here the
        Eiger triggers on each pulse's rising edge and times the frames itself
        from acquire_time/acquire_period -- the pulse width does not set the
        exposure the way it does in External Enable. Use for externally
        synchronized acquisition timed by the Eiger rather than by softglue.

        The pulse width is still the SHUTTER open time, though, so it is set to
        one segment (num_frames * acq_period): beam is on while a segment is
        being collected and off in the dead time between segments. The shutter
        cannot instead be held open by shutteroff()/showbeam() the way Internal
        Series does it -- shutter logic "Normal" gates off the entire softglue
        output, so the detector would receive no triggers at all (measured
        2026-09-02: 0 pulses on "Normal" vs 3 on "Override").

        This is the only mode that takes num_segments (expt.num_segments) and
        trigger_period (expt.trigger_period), both run-state values set
        per-protocol from measurement_info.yaml. softglue sends num_segments pulses spaced
        trigger_period apart, and each pulse fires a segment of num_frames
        images, so num_frames is the count PER SEGMENT and the total written to
        HDF is num_frames * num_segments:

            cam.num_images     = num_frames     (frames in one segment)
            cam.num_triggers   = num_segments   (how many segments)
            hdf1.num_capture   = num_frames * num_segments
            softglue.acq_period = trigger_period (gap BETWEEN segments)

        Three timings, easy to confuse:

            acq_time       exposure of one frame
            acq_period     frame-to-frame spacing inside a segment, so a
                           segment lasts num_frames * acq_period
            trigger_period segment-to-segment spacing. The pulse width holds
                           the shutter open, so each segment is two shutter
                           movements and validate_timing() requires both to be
                           >= 0.1 s:
                               open   = num_frames * acq_period
                               closed = trigger_period - num_frames*acq_period
                           The closed-time rule also keeps trigger_period
                           longer than a segment -- otherwise the next pulse
                           lands mid-segment, the detector drops it, and then
                           hangs on triggers that were already spent.

        num_segments defaults to 1 when a protocol omits it; trigger_period is
        required. No other mode reads either, and num_frames stays the plain
        total everywhere else.

    External Enable (trigger_mode="External Enable")
        Uses softglue to generate one pulse per frame; softglue's acq_time is
        what actually sets the exposure length (the same value is also
        written to the Eiger's own acquire_time to avoid confusion, but the
        Eiger ignores its own copy in this mode). Use for externally
        synchronized acquisition timed by softglue.
"""

import time as ttime

from id8_common.plans.acquire.acq_helpers import get_common_file_path
from id8_common.plans.acquire.acq_helpers import get_connected_device
from id8_common.plans.acquire.acq_wait import STATE_FAILED
from id8_common.plans.acquire.acq_wait import STATE_IDLE
from id8_common.plans.acquire.acq_wait import DetectorWaitError
from id8_common.plans.acquire.acq_wait import cam_fault
from id8_common.plans.acquire.acq_wait import hdf_frames_written
from id8_common.plans.acquire.acq_wait import hdf_progress
from id8_common.plans.acquire.acq_wait import wait_until
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.plans.set.shutter_att import shutteron

#: Seconds to let the HDF plugin finish writing after the cam stops. 60 s is
#: the old 600 x 0.1 s poll budget, now enforced rather than abandoned.
HDF_DRAIN_TIMEOUT = 60.0

#: ADStatus, for messages. Only 6/9/10 are terminal (acq_wait.STATE_FAILED);
#: the rest reach Idle by themselves and must be waited out, not interrupted.
STATE_NAMES = {0: "Idle", 1: "Acquire", 2: "Readout", 3: "Correct", 4: "Saving",
               5: "Aborting", 6: "Error", 7: "Waiting", 8: "Initializing",
               9: "Disconnected", 10: "Aborted"}

#: How long to let a transient state finish before treating it as stuck.
EIGER_SETTLE_TIMEOUT = 20.0

#: Longest to spend on each step of returning the camera to Idle.
EIGER_DISARM_TIMEOUT = 5.0
EIGER_RECOVER_TIMEOUT = 15.0

#: Exposure for the throwaway frame in step 2 below. Short enough to cost
#: nothing, long enough that the DCU treats it as a real series.
EIGER_RECOVER_EXPOSURE = 0.01


def recover_eiger_idle(eiger4M=None):
    """Return the Eiger to Idle before arming it. Raises if it cannot.

    Ctrl+C out of an acquisition leaves the camera in Aborted (ADStatus 10).
    Aborted is TERMINAL -- it never becomes Idle on its own -- and arming over it
    is what hangs the next measurement: ``cam.acquire`` goes to 1 and never comes
    back, and the exposure waits in this module are deliberately unbounded, so the
    plan sits there for ever instead of raising.

    Why not just treat Aborted as Idle and carry on? Because that changes what we
    BELIEVE about the detector without changing the detector. The DCU still holds
    the abandoned series, so the next arm hangs exactly as it did. The state has to
    be cleared, not reinterpreted.

    Two steps, cheapest first:

    1. **Disarm** -- ``cam.acquire.put(0)``. No beam, no frames, no time. On its own
       this is usually enough, so step 2 is rarely reached.
    2. **One throwaway frame** -- a single short exposure with the HDF plugin off,
       so nothing is written, and the beam blocked first, so the sample sees
       nothing. Completing a trivial series is what walks the DCU back to Idle when
       a bare disarm will not.

**Nothing is restored.** Every parameter step 2 touches is rewritten by the
    calling setup_eiger_*() on the very next lines, so putting the old value back
    here is not merely redundant, it is harmful: the restore and the setup are two
    rapid writes to the same DCU-backed parameter, and the earlier one can be the
    one that sticks. That showed up as an exposure time visibly flashing back to
    the previous protocol's value and then staying there, ignoring
    measurement_info.yaml (reported 2026-09-10, from an earlier version of this
    function that saved and restored six parameters).

    This is safe only because every setup_eiger_*() is self-contained about the
    camera parameters. setup_eiger_internal_enable() did not write num_images and
    was made to, rather than restoring it here -- one write beats two.

    Returns True if it had to do anything, False if the camera was already Idle.
    """
    if eiger4M is None:
        eiger4M = get_connected_device("eiger4M")
    cam = eiger4M.cam

    def state():
        try:
            return int(cam.detector_state.get())
        except Exception:
            return None

    now = state()
    if now == STATE_IDLE:
        return False

    # Only Error, Disconnected and Aborted are terminal. Everything else --
    # Acquire, Readout, Correct, Saving, Aborting, Waiting, Initializing --
    # reaches Idle on its own, and intervening is actively harmful: this runs at
    # the top of every setup, so in a repeat series it can catch the PREVIOUS
    # repeat still reading out or saving, and cam.acquire.put(0) there truncates
    # that repeat's file. The first version tested only "not Idle" and so took a
    # throwaway frame after a perfectly healthy acquisition (reported 2026-09-10).
    if now not in STATE_FAILED:
        try:
            wait_until(lambda: state() == STATE_IDLE,
                       timeout=EIGER_SETTLE_TIMEOUT,
                       what=f"eiger4M finishing {STATE_NAMES.get(now, now)}")
            return False          # settled by itself; nothing to clear
        except DetectorWaitError:
            print(f"eiger4M stuck in {STATE_NAMES.get(state(), state())} for "
                  f"{EIGER_SETTLE_TIMEOUT:g}s -- clearing it before arming.")
    else:
        print(f"eiger4M is in {STATE_NAMES.get(now, now)}, not Idle -- "
              f"clearing it before arming.")

    # --- 1. disarm -------------------------------------------------------
    try:
        cam.acquire.put(0)
    except Exception as exc:
        print(f"  disarm write failed ({exc}); trying a throwaway frame")
    else:
        try:
            wait_until(lambda: state() == STATE_IDLE,
                       timeout=EIGER_DISARM_TIMEOUT, what="eiger4M disarm")
            print("  Idle after disarm.")
            return True
        except DetectorWaitError:
            print("  still not Idle after disarm; taking one throwaway frame")

    # --- 2. throwaway frame, no beam, no file ----------------------------
    blockbeam()
    try:
        eiger4M.hdf1.capture.put(0)
        cam.manual_trigger.put("Disable")   # else it waits for a trigger that never comes
        cam.trigger_mode.put("Internal Series")
        cam.num_images.put(1)
        cam.num_triggers.put(1)
        cam.acquire_time.put(EIGER_RECOVER_EXPOSURE)
        cam.acquire_period.put(EIGER_RECOVER_EXPOSURE)
        cam.acquire.put(1)
        wait_until(lambda: state() == STATE_IDLE,
                   timeout=EIGER_RECOVER_TIMEOUT, what="eiger4M recovery frame")
    except DetectorWaitError as exc:
        raise RuntimeError(
            f"eiger4M will not leave {STATE_NAMES.get(state(), state())}: {exc}. "
            f"Arming anyway is the hang this "
            f"check exists to prevent. Restart the Eiger IOC, or put cam.acquire to 0 "
            f"by hand and confirm DetectorState_RBV reads Idle."
        ) from exc
    finally:
        try:
            cam.acquire.put(0)
        except Exception:
            pass

    print("  Idle after a throwaway frame.")
    return True

# =============================================================================
# Detector setup functions
#
# All four have the same shape: look up the devices, build the file path, then
# push this mode's settings to the detector (and to softglue, for the external
# modes). They configure only -- nothing starts until the matching acquire_*
# function runs. Each returns the NeXus metadata file path, which the caller
# (det_acq_series) fills in once the run is over.
# =============================================================================

def setup_eiger_internal(acq_time, num_frames, file_header, file_name):
    """Configure the Eiger for "Internal Series". Return the metadata file path."""
    eiger4M = get_connected_device("eiger4M")
    # Ctrl+C out of the previous run leaves the camera Aborted; arming over
    # that hangs. See recover_eiger_idle().
    recover_eiger_idle(eiger4M)
    file_path = get_common_file_path(file_header, file_name)

    # This mode has no separate acq_period (needs_acq_period is False in the
    # table below, so det_acq_series never passes one), so the period is the
    # exposure: frames run back-to-back with no gap.
    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_time)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)

    eiger4M.cam.trigger_mode.put("Internal Series")
    eiger4M.cam.num_images.put(num_frames)
    eiger4M.cam.num_triggers.put(1)
    # Explicit, not inherited: setup_eiger_internal_enable() turns manual
    # triggering ON and nothing used to turn it back off, so any mode that ran
    # after it sat waiting for a software trigger that never came -- detector
    # Done, HDF Capturing, 0 frames, forever. Each mode now states what it
    # needs. (Found 2026-09-04.)
    eiger4M.cam.manual_trigger.put("Disable")

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_eiger_internal_enable(acq_time, acq_period, num_frames, file_header, file_name):
    """Configure the Eiger for "Internal Enable". Return the metadata file path."""
    eiger4M = get_connected_device("eiger4M")
    # Ctrl+C out of the previous run leaves the camera Aborted; arming over
    # that hangs. See recover_eiger_idle().
    recover_eiger_idle(eiger4M)
    file_path = get_common_file_path(file_header, file_name)

    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_period)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)

    eiger4M.cam.trigger_mode.put("Internal Enable")
    # The only mode that wants this on. Every other eiger setup function
    # explicitly puts it back to "Disable" -- see the note there.
    eiger4M.cam.manual_trigger.put("Enable")
    eiger4M.cam.num_triggers.put(num_frames)
    # One image per software trigger, num_triggers of them. The detector ignores
    # num_images in this mode, but it is written anyway so this setup states the
    # full camera configuration on its own: recover_eiger_idle() may have left it
    # at 1, and it restores nothing (see the note there about racing writes).
    eiger4M.cam.num_images.put(1)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_eiger_external_series(
    acq_time, acq_period, num_frames, num_segments, trigger_period, file_header, file_name
):
    """Configure the Eiger and softglue for "External Series". Return the metadata path.

    One softglue pulse per SEGMENT of num_frames images -- see the module
    docstring for how acq_period and trigger_period differ.
    """
    eiger4M = get_connected_device("eiger4M")
    # Ctrl+C out of the previous run leaves the camera Aborted; arming over
    # that hangs. See recover_eiger_idle().
    recover_eiger_idle(eiger4M)
    softglue = get_connected_device("softglue")
    # Not referenced below. get_connected_device() does connect it, so the only
    # effect of this line is that a dead mz2 softglue fails here at setup
    # rather than part way through the acquisition.
    softglue_8id_mz2 = get_connected_device("softglue_8id_mz2")

    file_path = get_common_file_path(file_header, file_name)

    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_period)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    # num_frames is per segment here (see module docstring), so the file holds
    # one segment's worth of frames for each of the num_segments pulses.
    eiger4M.hdf1.num_capture.put(num_frames * num_segments)

    eiger4M.cam.trigger_mode.put("External Series")
    eiger4M.cam.num_images.put(num_frames)
    eiger4M.cam.num_triggers.put(num_segments)
    # Explicit, not inherited: setup_eiger_internal_enable() turns manual
    # triggering ON and nothing used to turn it back off, so any mode that ran
    # after it sat waiting for a software trigger that never came -- detector
    # Done, HDF Capturing, 0 frames, forever. Each mode now states what it
    # needs. (Found 2026-09-04.)
    eiger4M.cam.manual_trigger.put("Disable")

    # The softglue pulse WIDTH is the shutter open time, so it spans one whole
    # segment. The Eiger triggers on the pulse's rising edge and times the
    # frames itself, so the width does not affect exposure here -- but the
    # shutter is AND-gated by this same pulse whenever shutter logic is on
    # "Override" (which acquire_eiger_external's shutteron() sets, and which is
    # required: with logic on "Normal" the whole softglue output is blocked and
    # the detector receives no triggers at all -- measured on 2026-09-02,
    # 0 pulses vs 3). A negligible 1e-5 s width would therefore open the
    # shutter for 10 us per segment and the entire dataset would be dark.
    #
    # Width < trigger_period is guaranteed by validate_timing(), so consecutive
    # pulses stay separated and each segment still gets its own rising edge.
    softglue.acq_time.put(acq_period * num_frames)

    # One pulse per SEGMENT, not per frame, and the spacing is its own number:
    # trigger_period, not acq_period. acq_period paces frames *inside* a
    # segment (so a segment lasts num_frames * acq_period); trigger_period is
    # the gap between segments and has to be longer than that, or the next
    # pulse arrives while the detector is still busy and is silently dropped.
    # master_plan.validate_timing() enforces it. (Observed on 2026-09-02 with
    # the spacing left at acq_period: softglue emitted all 5 pulses in 8 s, the
    # detector consumed the first, spent 20 s on segment 1, and hung waiting
    # for 4 triggers already spent -- 10 of 50 frames in the file.) External
    # Enable is one pulse per frame, which is why it passes acq_period through.
    softglue.acq_period.put(trigger_period)
    softglue.num_triggers.put(num_segments)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_eiger_external(acq_time, acq_period, num_frames, file_header, file_name):
    """Configure the Eiger and softglue for "External Enable". Return the metadata path.

    One softglue pulse per FRAME, and softglue's acq_time -- not the Eiger's --
    is what sets the exposure length in this mode.
    """
    eiger4M = get_connected_device("eiger4M")
    # Ctrl+C out of the previous run leaves the camera Aborted; arming over
    # that hangs. See recover_eiger_idle().
    recover_eiger_idle(eiger4M)
    softglue = get_connected_device("softglue")
    # Not referenced below; see the same line in setup_eiger_external_series().
    softglue_8id_mz2 = get_connected_device("softglue_8id_mz2")

    file_path = get_common_file_path(file_header, file_name)

    eiger4M.cam.acquire_time.put(acq_time)
    eiger4M.cam.acquire_period.put(acq_period)

    eiger4M.hdf1.file_name.put(file_name)
    eiger4M.hdf1.file_path.put(file_path)
    eiger4M.hdf1.num_capture.put(num_frames)

    eiger4M.cam.num_triggers.put(num_frames)
    # One image per external pulse, num_triggers of them. Written explicitly for
    # the same reason as manual_trigger below: this setup states the full camera
    # configuration rather than inheriting it. recover_eiger_idle() may have left
    # num_images at 1, and it restores nothing (see the note there).
    eiger4M.cam.num_images.put(1)
    eiger4M.cam.trigger_mode.put("External Enable")
    # Explicit, not inherited: setup_eiger_internal_enable() turns manual
    # triggering ON and nothing used to turn it back off, so any mode that ran
    # after it sat waiting for a software trigger that never came -- detector
    # Done, HDF Capturing, 0 frames, forever. Each mode now states what it
    # needs. (Found 2026-09-04.)
    eiger4M.cam.manual_trigger.put("Disable")

    softglue.acq_time.put(acq_time)
    softglue.acq_period.put(acq_period)
    softglue.num_triggers.put(num_frames)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


# =============================================================================
# Detector acquire functions
#
# Each takes no arguments: the matching setup_* function has already pushed
# every parameter to the hardware. They run the same stages in the same order
# -- open beam, arm the file writer, arm/start the detector, wait for the
# frames, close beam, wait for the file to finish -- and the order is the
# hardware protocol, so do not reshuffle the puts.
# =============================================================================

def acquire_eiger_internal():
    """Run one "Internal Series": the detector times the whole series itself."""
    eiger4M = get_connected_device("eiger4M")

    # Beam on for the whole series. shutteroff() takes the fast shutter out of
    # softglue's control so showbeam() drives it directly -- this mode has no
    # pulse train, so nothing else would open it.
    shutteroff()
    showbeam()
    ttime.sleep(0.1)

    # Arm the file writer before the detector, so the plugin is already
    # capturing when the first frame arrives.
    eiger4M.hdf1.capture.put(1)
    eiger4M.cam.acquire.put(1)
    ttime.sleep(0.5)

    # Wait out the exposure series. Unbounded on purpose for now: a detector
    # fault here hangs the plan instead of raising. (acq_wait.wait_until() is
    # the bounded version, used below in acquire_eiger_external.)
    while eiger4M.cam.acquire.get() == 1:
        ttime.sleep(0.1)

    blockbeam()

    # Drain: the HDF plugin keeps writing for a moment after the last frame.
    while eiger4M.hdf1.capture.get() == 1:
        ttime.sleep(0.1)


def acquire_eiger_internal_enable():
    """Fire num_triggers software triggers, one exposure per trigger (see
    module docstring: Internal Enable)."""
    eiger4M = get_connected_device("eiger4M")

    # Beam on for the whole series -- see acquire_eiger_internal().
    shutteroff()
    showbeam()
    ttime.sleep(0.1)

    # Arm both. The detector now sits waiting for software triggers rather
    # than exposing, because setup put manual_trigger on "Enable".
    eiger4M.hdf1.capture.put(1)
    eiger4M.cam.acquire.put(1)
    ttime.sleep(0.5)

    # Read the pacing back off the cam rather than taking arguments, so this
    # stays a zero-argument acquire like the others.
    num_triggers = int(eiger4M.cam.num_triggers.get())
    acq_period = eiger4M.cam.acquire_period.get()

    # One software trigger per frame, paced by acq_period. wait=False, so the
    # put returns at once and the sleep below is what sets the spacing.
    for _ in range(num_triggers):
        eiger4M.cam.special_trigger_button.put(1, wait=False)
        ttime.sleep(acq_period)

    # Unbounded on purpose for now -- see acquire_eiger_internal().
    while eiger4M.cam.acquire.get() == 1:
        ttime.sleep(0.1)

    blockbeam()

    # Drain the HDF plugin.
    while eiger4M.hdf1.capture.get() == 1:
        ttime.sleep(0.1)


def acquire_eiger_external():
    """Run one softglue-triggered acquisition. Serves both external modes.

    "External Series" and "External Enable" differ only in what their setup_*
    function programmed into the cam and softglue, so the run sequence here --
    clear softglue, arm the detector, start the pulse train, wait -- is shared.
    """
    eiger4M = get_connected_device("eiger4M")
    softglue = get_connected_device("softglue")

    # Reset softglue's pulse counters, so this run starts from zero rather than
    # continuing whatever the previous run left behind.
    softglue.pv_clear1.put("1!")
    softglue.pv_clear2.put("1!")

    # Beam under softglue control: each pulse gates the shutter open, instead
    # of it staying open for the whole run the way the internal modes do. The
    # module docstring explains why this cannot be shutteroff()/showbeam().
    shutteron()
    showbeam()
    ttime.sleep(0.1)

    # Arm the file writer, then arm the detector. Neither does anything yet --
    # the detector is waiting for its first external pulse.
    eiger4M.hdf1.capture.put(1)
    eiger4M.cam.acquire.put(1)
    ttime.sleep(1.0)

    # Trigger: softglue now emits the pulse train the setup programmed.
    softglue.start_pulses.put("1!")

    # Wait for the frames. Unbounded on purpose for now; only the HDF drain
    # below is bounded. acquire_busy (not acquire) is the signal that the
    # detector is still working through its triggers.
    while True:
        #### QZ on 2026/01/06 ####
        # without the 0.5 s wait time, the repeating acqs go out of sync.
        # Don't know why and maybe the 0.5 s can be made shorter
        #### QZ on 2026/01/06 ####
        ttime.sleep(0.5)
        det_status = eiger4M.cam.acquire_busy.get()
        if det_status == 1:
            ttime.sleep(0.1)
        if det_status == 0:
            break
    blockbeam()

    # Raises if the plugin does not finish. Until 2026-09-06 this loop polled 600
    # times and then simply carried on, force-closing a SHORT file -- and
    # det_acq_series() went straight on to write NeXus metadata and submit a DM
    # job for it, with nothing printed. A truncated dataset that looks complete
    # is worse than a failed one.
    try:
        # hdf_frames_written() returns (frames written, frames expected), so
        # [0] >= [1] means the plugin has written everything it was asked for.
        wait_until(
            lambda: hdf_frames_written(eiger4M.hdf1)[0] >= hdf_frames_written(eiger4M.hdf1)[1],
            timeout=HDF_DRAIN_TIMEOUT,
            what=hdf_progress(eiger4M.hdf1, "eiger4M"),
            faults=(cam_fault(eiger4M.cam, "eiger4M"),),
        )
    finally:
        eiger4M.hdf1.capture.put(0)


# =============================================================================
# Mode table
#
# ad_acq.py mounts this as ACQ_MODES["eiger4M"], and det_acq_series() looks up
# ACQ_MODES[expt.det_name][expt.det_mode] to find out how to run one
# measurement. What the keys mean:
#
#   setup                 function that configures the hardware
#   acquire               function that then runs it
#   needs_acq_period      pass acq_period to setup as well? Every setup is
#                         always given acq_time, num_frames, file_header and
#                         file_name; the needs_* flags add the rest, so a
#                         setup function takes exactly what it uses
#   needs_num_segments    pass num_segments to setup as well?
#   needs_trigger_period  pass trigger_period to setup as well?
#   required_devices      all must connect before the run starts. Listing
#                         "softglue" also tells cleanup_acquisition() to stop
#                         the pulse train if the plan is aborted
#
# Adding a mode means adding a row here plus its setup_/acquire_ pair above.
# needs_acq_period must be present; the other needs_* flags default to False.
# =============================================================================

EIGER4M_MODES = {
    "Internal Series": {
        "setup": setup_eiger_internal,
        "acquire": acquire_eiger_internal,
        "needs_acq_period": False,
        "required_devices": ["eiger4M"],
    },
    "Internal Enable": {
        "setup": setup_eiger_internal_enable,
        "acquire": acquire_eiger_internal_enable,
        "needs_acq_period": True,
        "required_devices": ["eiger4M"],
    },
    "External Series": {
        "setup": setup_eiger_external_series,
        "acquire": acquire_eiger_external,
        "needs_acq_period": True,
        "needs_num_segments": True,
        "needs_trigger_period": True,
        "required_devices": ["eiger4M", "softglue", "softglue_8id_mz2"],
    },
    "External Enable": {
        "setup": setup_eiger_external,
        "acquire": acquire_eiger_external,
        "needs_acq_period": True,
        "required_devices": ["eiger4M", "softglue", "softglue_8id_mz2"],
    },
}
