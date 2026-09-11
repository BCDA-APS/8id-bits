"""
TV mode: a live view on the Rigaku, the Eiger and the Lambda at once.

What people have always called "TV mode" at 8-ID is a detector left free-running
from its own GUI so you can watch the image update. This is that, scripted, and
on all three detectors together instead of one at a time:

    tv_mode()                                # 5000 frames of 1 s, all three
    tv_mode(acq_time=0.1, num_frames=200)    # something shorter
    tv_mode(detectors=["lambda2M"])          # just one

WRITES NOTHING. No HDF file, no NeXus metadata, no DM job, and -- because it
never calls gen_folder_prefix() -- no measurement number is consumed. The HDF
plugins are explicitly DISARMED rather than merely left alone, so a plugin an
aborted measurement left capturing cannot quietly append 5000 live-view frames
to that measurement's file.

Blocks the session. tv_mode() returns only when all three detectors have
finished their frames, which at the defaults is 5000 s -- roughly 83 minutes.
That is deliberate: the number is a ceiling, not a target. The expected way to
end a TV session is Ctrl+C, and Ctrl+C presses Stop on every detector and closes
the shutter before it returns. So does any other failure: the cleanup is in a
`finally`, so there is no path out of this function that leaves a detector
running or the beam on the sample.

BEAM IS ON for the whole run. tv_mode() calls shutteroff() then showbeam()
before arming, exactly as the serial acquire_*() functions do, and blockbeam()
on the way out. Set the attenuation you want with att() BEFORE starting it --
the session is blocked once it is running, so you cannot change it partway.
At the default 5000 frames that is up to 83 minutes of continuous beam on one
spot; Ctrl+C when you have seen enough.

The three detectors are NOT synchronised. Each is armed with its own EPICS put,
one after the other, and each times its own frames. They overlap; they do not
share a clock.
"""

import time as ttime
from datetime import datetime

from id8_common.plans.acquire.acq_wait import acq_time_mismatch
from id8_common.plans.acquire.acq_wait import acq_timeout
from id8_common.plans.acquire.acq_wait import cam_fault
from id8_common.plans.acquire.acq_wait import stop_acquiring
from id8_common.plans.acquire.acq_wait import wait_until
from id8_common.plans.acquire.eiger4m_modes import recover_eiger_idle
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.registry import get_connected_device

#: Defaults, as specified for the beamline's TV mode.
TV_ACQ_TIME = 1.0
TV_NUM_FRAMES = 5000

#: Seconds to wait for a detector to confirm it has started after being armed.
TV_START_TIMEOUT = 30.0

#: Seconds between polls of the "are they still going?" loop, and between the
#: progress lines it prints. The poll is slower than an acquisition's because
#: nothing downstream is waiting on the answer.
TV_POLL_INTERVAL = 0.5
TV_PROGRESS_INTERVAL = 30.0

#: Allowance on top of acq_time * num_frames before the wait gives up. Generous:
#: its job is to stop a wedged detector from blocking the session for ever, not
#: to police timing.
TV_TIMEOUT_FACTOR = 2.0
TV_TIMEOUT_FLOOR = 120.0

#: How long to let an acquire_time readback catch up with the value just written.
TV_ACQ_TIME_SETTLE = 10.0


def timestamp():
    """Log timestamp, matching the format the acquisition modules print."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# =============================================================================
# Per-detector configuration
# =============================================================================
# Free-running internal timing on each detector, and nothing else: cam
# parameters only, no file plugin, no softglue, no shutter. Each takes the
# device and returns nothing.
#
# These are NOT the ACQ_MODES setup functions. Those exist to prepare a
# measurement -- they build a file path, fill in the HDF plugin and hand back a
# NeXus metadata path, all of which would need a folder prefix and so would burn
# a measurement number for a live view that saves nothing.


def configure_rigaku_tv(rigaku3M, acq_time, num_frames):
    """Rigaku, internally timed, full frames through the areaDetector pipeline."""
    cam = rigaku3M.cam

    # 'Fixed Time' is the Rigaku's internally-timed trigger mode: it clocks its
    # own frames and needs nothing from softglue, which is what makes it the
    # analogue of the Eiger's 'Internal Series' and the Lambda's 'Internal'.
    #
    # Confirmed 2026-09-11 to be a legal value of 8idRigaku3m:cam1:TriggerMode
    # (the enum is Fixed Time / Continuous with Trigger / Start with Trigger /
    # Gated Trigger / Trigger Sync / Fixed, Start with Trigger / Burst Mode) --
    # but NOT yet run on the detector. If it turns out not to free-run, the
    # proven fallback is the acquisition path's combination: 'Start with
    # Trigger' plus softglue.enable_rigaku = '1'.
    cam.trigger_mode.put("Fixed Time")
    cam.image_mode.put("Standard")
    cam.output_control.put("areaDetector")
    cam.output_resolution.put("16 Bit")

    cam.acquire_time.put(acq_time)
    cam.acquire_period.put(acq_time)
    cam.num_images.put(num_frames)


def configure_eiger_tv(eiger4M, acq_time, num_frames):
    """Eiger, 'Internal Series' -- one arm, num_frames back-to-back exposures."""
    cam = eiger4M.cam

    # Ctrl+C out of a previous run leaves the camera Aborted, which is terminal,
    # and arming over it hangs. Same call setup_eiger_internal() makes, and
    # doubly relevant here because Ctrl+C is the EXPECTED way to end TV mode.
    recover_eiger_idle(eiger4M)

    cam.trigger_mode.put("Internal Series")
    # Explicit: 'Internal Enable' leaves manual triggering on, and a TV mode
    # that inherited it would sit waiting for a software trigger for ever.
    cam.manual_trigger.put("Disable")
    cam.num_triggers.put(1)

    cam.acquire_time.put(acq_time)
    cam.acquire_period.put(acq_time)
    cam.num_images.put(num_frames)


def configure_lambda_tv(lambda2M, acq_time, num_frames):
    """Lambda, 'Internal' -- free-running internal timing."""
    cam = lambda2M.cam

    cam.trigger_mode.put("Internal")
    # 'Single' would give one frame however large num_images is.
    cam.image_mode.put("Multiple")

    cam.acquire_time.put(acq_time)
    cam.acquire_period.put(acq_time)
    cam.num_images.put(num_frames)


#: Registry name -> configure function. Insertion order is arm order.
TV_DETECTORS = {
    "rigaku3M": configure_rigaku_tv,
    "eiger4M": configure_eiger_tv,
    "lambda2M": configure_lambda_tv,
}


# =============================================================================
# Helpers
# =============================================================================


def tv_progress(cams):
    """One line of 'how far has each detector got', for the periodic report."""
    parts = []

    for cam in cams:
        try:
            done = int(cam["det"].cam.num_images_counter.get())
            parts.append(f"{cam['name']} {done}/{cam['num_frames']}")
        except Exception:
            parts.append(f"{cam['name']} ?")

    return ", ".join(parts)


def stop_tv_mode(cams):
    """Close the shutter and press Stop on every detector. Safe from any state.

    This is what Ctrl+C runs. Every step is in its own try, so one detector that
    will not answer cannot stop the others from being stopped -- and the beam is
    blocked first, before anything that could fail.
    """
    try:
        blockbeam()
    except Exception as exc:
        print(f"Could not block beam while stopping TV mode: {exc}")

    for cam in cams:
        try:
            cam["det"].cam.acquire.put(0)
            print(f"{timestamp()}, Stopped {cam['name']}")
        except Exception as exc:
            print(f"Could not stop {cam['name']}: {exc}")

        # Belt and braces: TV mode never arms the file plugin, but an earlier
        # aborted measurement may have left it capturing.
        try:
            if hasattr(cam["det"], "hdf1"):
                cam["det"].hdf1.capture.put(0)
        except Exception:
            pass


def resolve_tv_detectors(names, acq_time, num_frames):
    """Look up each requested detector and check it is one TV mode knows how to drive."""
    cams = []

    for name in names:
        if name not in TV_DETECTORS:
            raise ValueError(f"'{name}' is not a TV-mode detector. Known: {sorted(TV_DETECTORS)}")

        cams.append(
            {
                "name": name,
                "det": get_connected_device(name),
                "configure": TV_DETECTORS[name],
                "acq_time": float(acq_time),
                "num_frames": int(num_frames),
            }
        )

    return cams


# =============================================================================
# Main entry point
# =============================================================================


def tv_mode(acq_time=TV_ACQ_TIME, num_frames=TV_NUM_FRAMES, detectors=None):
    """Run a live view on all three detectors until the frames run out or you Ctrl+C.

    Args:
        acq_time: Exposure per frame, seconds. Also the period -- frames run
            back-to-back.
        num_frames: Frames per detector. A ceiling, not a target: Ctrl+C is the
            normal way to end a TV session.
        detectors: Registry names to run, defaulting to all of TV_DETECTORS
            (rigaku3M, eiger4M, lambda2M). Pass a shorter list to view one or two.

    Raises:
        RuntimeError: on Ctrl+C, matching det_acq_series() and trio_acq_series().
            Note this means an upstream ``except KeyboardInterrupt`` will not
            catch it.
        DetectorWaitError: if a detector will not stop, will not start, will not
            take the exposure it was given, or reaches a terminal state.

    Whatever the outcome, stop_tv_mode() has run before this returns or raises:
    beam blocked, Stop pressed on every detector.
    """
    names = list(detectors) if detectors else list(TV_DETECTORS)
    cams = resolve_tv_detectors(names, acq_time, num_frames)

    print(f"\n{timestamp()}, TV mode: {', '.join(names)} at {float(acq_time):g} s x {int(num_frames)} frames")
    print(f"{timestamp()}, Writing no files. Ctrl+C stops every detector and closes the shutter.")

    try:
        # 1. Nothing may already be running. This is both the "clear a previous
        #    TV session" step and the guard for the bug that a write to
        #    acquire_time on a RUNNING camera does not stick -- see
        #    acq_wait.stop_acquiring(). It has to happen before step 3.
        for cam in cams:
            if stop_acquiring(cam["det"].cam, cam["name"]):
                print(f"{timestamp()}, {cam['name']} was already acquiring -- stopped it")

        # 2. TV mode saves nothing, so make sure nothing is armed to save.
        for cam in cams:
            if hasattr(cam["det"], "hdf1"):
                cam["det"].hdf1.capture.put(0)

        # 3. Free-running internal timing on each.
        for cam in cams:
            cam["configure"](cam["det"], cam["acq_time"], cam["num_frames"])

        # 4. Prove the exposure took, before any beam is spent on it.
        confirm_tv_acq_time(cams)

        # ---- beam on ---------------------------------------------------------
        # shutteroff() takes the fast shutter out of softglue's control so
        # showbeam() drives it directly: none of these modes emits a pulse
        # train, so nothing else would open it.
        shutteroff()
        showbeam()
        ttime.sleep(0.1)

        # 5. Arm. Each put returns immediately, so by the end of this loop all
        #    three are integrating -- offset from each other by the puts, which
        #    is the only sense in which they are "unsynchronised".
        for cam in cams:
            cam["det"].cam.acquire.put(1)
            print(f"{timestamp()}, Armed {cam['name']}")

        for cam in cams:
            wait_until(
                lambda cam=cam: cam["det"].cam.acquire.get() == 1,
                timeout=TV_START_TIMEOUT,
                what=f"{cam['name']} to start acquiring",
                faults=(cam_fault(cam["det"].cam, cam["name"]),),
            )

        print(f"{timestamp()}, All detectors live. Ctrl+C to stop.")

        # 6. Hold the session until every detector has finished its frames.
        wait_for_tv_frames(cams)

        print(f"{timestamp()}, TV mode complete: {tv_progress(cams)}")

    except KeyboardInterrupt as exc:
        # Converted, not propagated, the same way det_acq_series() and
        # trio_acq_series() do it -- so a caller cannot swallow it with a bare
        # `except KeyboardInterrupt` and carry on as though nothing happened.
        # The `finally` below has already stopped everything by this point.
        raise RuntimeError("\n TV mode stopped by user (Ctrl+C).") from exc
    finally:
        # Runs on every path out: normal finish, Ctrl+C, detector fault, timeout.
        stop_tv_mode(cams)


def confirm_tv_acq_time(cams, timeout=TV_ACQ_TIME_SETTLE):
    """Check every camera is holding the exposure it was just given.

    Same check trio_acq_series() makes, for the same reason: a camera that was
    still acquiring when acquire_time was written keeps its old exposure, and
    without this the live view would silently run at the previous session's
    value. Polls rather than reading once, since a put and its readback are not
    simultaneous.
    """

    def mismatches():
        found = []

        for cam in cams:
            reason = acq_time_mismatch(cam["det"].cam, cam["acq_time"], cam["name"])

            if reason:
                found.append(reason)

        return found

    wait_until(
        lambda: not mismatches(),
        timeout=timeout,
        what=lambda: "exposure time did not take on " + "; ".join(mismatches()),
    )


def wait_for_tv_frames(cams):
    """Block until every detector has finished, printing progress as it goes.

    Bounded, like every other wait in acq_wait: a wedged detector must surface
    as an error rather than as a session that never comes back.
    """
    slowest = max(cam["acq_time"] * cam["num_frames"] for cam in cams)
    timeout = acq_timeout(slowest, 1, factor=TV_TIMEOUT_FACTOR, floor=TV_TIMEOUT_FLOOR)

    last_report = [ttime.time()]

    def every_detector_finished():
        # The progress line lives inside the predicate because wait_until()
        # already polls once per TV_POLL_INTERVAL and TV mode blocks the
        # session -- this is the only sign of life the user gets for what may
        # be over an hour.
        now = ttime.time()

        if now - last_report[0] >= TV_PROGRESS_INTERVAL:
            last_report[0] = now
            print(f"{timestamp()}, {tv_progress(cams)}")

        for cam in cams:
            if cam["det"].cam.acquire.get() == 1:
                return False

        return True

    return wait_until(
        every_detector_finished,
        timeout=timeout,
        what=lambda: f"all detectors to finish ({tv_progress(cams)})",
        faults=tuple(cam_fault(cam["det"].cam, cam["name"]) for cam in cams),
        poll=TV_POLL_INTERVAL,
    )


__all__ = [
    "TV_ACQ_TIME",
    "TV_DETECTORS",
    "TV_NUM_FRAMES",
    "stop_tv_mode",
    "tv_mode",
]
