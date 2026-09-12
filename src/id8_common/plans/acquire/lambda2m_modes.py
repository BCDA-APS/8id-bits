"""
Lambda2M mode definitions: setup/acquire functions plus LAMBDA2M_MODES, the
table ad_acq.py assembles into ACQ_MODES. Nothing here runs on its own. The
callers are det_acq_series() in ad_acq.py, and trio_acq_series() for the
("lambda2M", "Internal") trio leg.

Modes:

    Internal (trigger_mode="Internal")
        No external hardware. Free-running internal timing. The only mode a
        parallel run can use: it needs nothing from softglue, so it can share
        one beam window with a Rigaku-owned shutter.

    External (trigger_mode="External_ImagePer")
        Uses softglue to generate one pulse per frame, same pattern as
        Eiger's External Enable. NOT usable as a trio leg -- softglue gates the
        shutter per pulse, which cannot coexist with a single shutter window
        held open by another detector.
"""

import time as ttime

from id8_common.plans.acquire import acq_helpers
from id8_common.plans.acquire.acq_helpers import get_common_file_path
from id8_common.plans.acquire.acq_helpers import get_connected_device
from id8_common.plans.acquire.acq_wait import cam_fault
from id8_common.plans.acquire.acq_wait import hdf_frames_written
from id8_common.plans.acquire.acq_wait import hdf_progress
from id8_common.plans.acquire.acq_wait import stop_acquiring
from id8_common.plans.acquire.acq_wait import wait_until
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.plans.set.shutter_att import shutteron

#: Seconds to let the HDF plugin finish writing after the cam stops.
HDF_DRAIN_TIMEOUT = 60.0

#: Seconds to give the Lambda to drop out of live/TV mode before it is armed.
LIVE_STOP_TIMEOUT = 5.0


def stop_lambda_live(lambda2M=None):
    """Disarm the Lambda so the next arm is not swallowed. Raises if it will not stop.

    Thin Lambda-flavoured wrapper over acq_wait.stop_acquiring(), which carries
    the full explanation of why a camera left in live/TV mode has to be stopped
    before its exposure time is written.

    Both acquire_* functions below already do this immediately before arming, so
    for the serial path this is a harmless second write. It exists as its own
    function because the TRIO path never calls them: trio_acq_series() takes only
    the setup half of the mode table and arms the detector itself.

    Called from setup_lambda_internal(), which runs BEFORE the shutter opens --
    the only safe place for a wait that is allowed to fail.
    """
    if lambda2M is None:
        lambda2M = get_connected_device("lambda2M")

    return stop_acquiring(lambda2M.cam, "lambda2M", timeout=LIVE_STOP_TIMEOUT)


# =============================================================================
# Detector setup functions
#
# Both have the same shape: look up the devices, build the file path, then push
# this mode's settings to the detector (and to softglue, for External). They
# configure only -- nothing starts until the matching acquire_* function runs.
# Each returns the NeXus metadata file path, which the caller (det_acq_series)
# fills in once the run is over.
# =============================================================================

def setup_lambda_internal(acq_time, num_frames, file_header, file_name):
    """Configure the Lambda for "Internal". Return the metadata file path."""
    lambda2M = get_connected_device("lambda2M")
    # Live/TV mode has to be off before the detector is armed. The trio path
    # never reaches acquire_lambda_internal(), which is where that used to be
    # the only time it happened -- see stop_lambda_live().
    stop_lambda_live(lambda2M)
    file_path = get_common_file_path(file_header, file_name)

    lambda2M.hdf1.enable.put(1)
    lambda2M.cam.trigger_mode.put("Internal")

    # No separate acq_period in this mode (needs_acq_period is False in the
    # table below), so the period is the exposure: frames run back-to-back.
    lambda2M.cam.acquire_time.put(acq_time)
    lambda2M.cam.acquire_period.put(acq_time)

    lambda2M.hdf1.file_name.put(file_name)
    lambda2M.hdf1.file_path.put(file_path)

    lambda2M.cam.num_images.put(num_frames)
    lambda2M.hdf1.num_capture.put(num_frames)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


def setup_lambda_external(acq_time, acq_period, num_frames, file_header, file_name):
    """Configure the Lambda and softglue for "External". Return the metadata path.

    One softglue pulse per frame, the same pattern as the Eiger's
    "External Enable".
    """
    lambda2M = get_connected_device("lambda2M")
    softglue = get_connected_device("softglue")

    file_path = get_common_file_path(file_header, file_name)

    lambda2M.hdf1.enable.put(1)

    lambda2M.cam.acquire_time.put(acq_time)
    lambda2M.cam.acquire_period.put(acq_period)

    lambda2M.hdf1.file_name.put(file_name)
    lambda2M.hdf1.file_path.put(file_path)
    lambda2M.hdf1.num_capture.put(num_frames)

    lambda2M.cam.num_images.put(num_frames)
    lambda2M.cam.trigger_mode.put("External_ImagePer")

    softglue.acq_time.put(acq_time)
    softglue.acq_period.put(acq_period)
    softglue.num_triggers.put(num_frames)

    metadata_fname = f"{file_path}/{file_name}_metadata.hdf"

    return metadata_fname


# =============================================================================
# Detector acquire functions
#
# Both take no arguments: the matching setup_* function has already pushed
# every parameter to the hardware. They run the same stages -- open beam, arm
# the file writer, start the detector, wait for the frames, close beam, wait
# for the file to finish -- and the order is the hardware protocol, so do not
# reshuffle the puts.
# =============================================================================

def acquire_lambda_internal():
    """Run one "Internal" series: the detector times the whole series itself."""
    lambda2M = get_connected_device("lambda2M")
    # Stop any TV/live mode that might be running before we arm.
    lambda2M.cam.acquire.put(0)

    # Beam on for the whole series. shutteroff() takes the fast shutter out of
    # softglue's control so showbeam() drives it directly -- this mode has no
    # pulse train, so nothing else would open it.
    shutteroff()
    showbeam()

    # Arm the file writer before the detector, so the plugin is already
    # capturing when the first frame arrives.
    lambda2M.hdf1.capture.put(1)
    lambda2M.cam.acquire.put(1)

    #  '''keithley voltage sequence'''
    # while lambda2M.cam.acquire.get() == 1:
    #     volt_cycle_single(voltage_file = np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_single.txt'))
        # volt_cycle_series(voltage_file_series=np.loadtxt('/home/beams10/8IDIUSER/bluesky/src/id8_common/plans/set/voltage_program_series.txt')

    # Wait out the exposure series. Unbounded on purpose for now: a detector
    # fault here hangs the plan instead of raising. (acq_wait.wait_until() is
    # the bounded version, used below for the HDF drain in External mode.)
    ttime.sleep(0.5)
    while lambda2M.cam.acquire.get() == 1:
        ttime.sleep(0.05)

    blockbeam()

    # Drain: the HDF plugin keeps writing for a moment after the last frame.
    while lambda2M.hdf1.capture.get() == 1:
        ttime.sleep(0.05)


def acquire_lambda_external():
    """Run one softglue-triggered acquisition, one pulse per frame."""
    lambda2M = get_connected_device("lambda2M")
    softglue = get_connected_device("softglue")
    softglue_8id_mz2 = get_connected_device("softglue_8id_mz2")
    # dpKeysight = get_connected_device("dpKeysight")

    # Stop any TV/live mode that might be running before external acquisition.
    lambda2M.cam.acquire.put(0)

    # Reset softglue's pulse counters, so this run starts from zero rather than
    # continuing whatever the previous run left behind.
    softglue.pv_clear1.put("1!")
    softglue.pv_clear2.put("1!")

    # Beam under softglue control: each pulse gates the shutter open, instead
    # of it staying open for the whole run the way Internal mode does it.
    shutteron()
    showbeam()
    ttime.sleep(0.1)

    # Arm the file writer, then arm the detector. Neither does anything yet --
    # the detector is waiting for its first external pulse.
    lambda2M.hdf1.capture.put(1)
    lambda2M.cam.acquire.put(1)
    ttime.sleep(1.0)

    # Trigger: softglue now emits one pulse per frame.
    softglue.start_pulses.put("1!")
    softglue_8id_mz2.load.put("1!")

    # Wait for the frames, running any hooks the protocol loaded (the `hooks:`
    # block of measurement_info.yaml) on every pass. With no hooks loaded
    # run_hooks() returns immediately, so this becomes a tight poll of the
    # acquire PV with no sleep in it. Unbounded on purpose for now.
    while lambda2M.cam.acquire.get() == 1:
        acq_helpers.run_hooks(acq_helpers.active_hooks)

    # A second, slower wait on the same PV. The loop above has already seen
    # acquire go to 0, so unless it comes back this only adds one 0.5 s pause
    # before the shutter closes.
    while True:
        ttime.sleep(0.5)
        det_status = lambda2M.cam.acquire.get()
        if det_status == 1:
            ttime.sleep(0.1)
        if det_status == 0:
            break
    blockbeam()

    # Nothing below reads either of these; the drain check uses
    # hdf_frames_written() instead.
    frame_num_set = lambda2M.hdf1.num_capture.get()
    count = 0

    # dpKeysight.output.put(0)


    # Raises if the plugin does not finish -- see the same fix in eiger4m_modes.py.
    try:
        # hdf_frames_written() returns (frames written, frames expected), so
        # [0] >= [1] means the plugin has written everything it was asked for.
        wait_until(
            lambda: hdf_frames_written(lambda2M.hdf1)[0] >= hdf_frames_written(lambda2M.hdf1)[1],
            timeout=HDF_DRAIN_TIMEOUT,
            what=hdf_progress(lambda2M.hdf1, "lambda2M"),
            faults=(cam_fault(lambda2M.cam, "lambda2M"),),
        )
    finally:
        lambda2M.hdf1.capture.put(0)


# =============================================================================
# Mode table
#
# ad_acq.py mounts this as ACQ_MODES["lambda2M"], and det_acq_series() looks up
# ACQ_MODES[expt.det_name][expt.det_mode] to find out how to run one
# measurement. What the keys mean:
#
#   setup             function that configures the hardware
#   acquire           function that then runs it
#   needs_acq_period  pass acq_period to setup as well? Every setup is always
#                     given acq_time, num_frames, file_header and file_name;
#                     this flag adds acq_period, so a setup function takes
#                     exactly what it uses
#   required_devices  all must connect before the run starts. Listing
#                     "softglue" also tells cleanup_acquisition() to stop the
#                     pulse train if the plan is aborted
#
# Adding a mode means adding a row here plus its setup_/acquire_ pair above.
# =============================================================================

# Two flags on every entry decide whether a mode may share a beam window with
# another detector (see parallel_objection() in ad_acq.py, which reads them):
#
#   drives_shutter  the mode routes the fast shutter through softglue, so the
#                   beam is gated by the detector's own trigger path instead of
#                   by the plan's showbeam()/blockbeam(). Such a mode cannot run
#                   alongside another detector -- the single window around the
#                   set would fight it.
#   self_paced      one arm runs the whole acquisition: the detector clocks every
#                   frame itself off acquire_time/acquire_period, with no
#                   external pulse train and no per-frame software trigger. Only
#                   a self-paced mode can be armed and left while the plan waits
#                   on all the other legs.
#
# Both must be stated. A mode that omits them is refused from parallel runs
# rather than assumed safe -- see parallel_objection().

LAMBDA2M_MODES = {
    "Internal": {
        "setup": setup_lambda_internal,
        "acquire": acquire_lambda_internal,
        "needs_acq_period": False,
        "required_devices": ["lambda2M"],
        # shutteroff() + showbeam() for the whole series, detector clocks itself.
        "drives_shutter": False,
        "self_paced": True,
    },
    "External": {
        "setup": setup_lambda_external,
        "acquire": acquire_lambda_external,
        "needs_acq_period": True,
        "required_devices": ["lambda2M", "softglue", "softglue_8id_mz2"],
        # One softglue pulse per frame, gating the shutter -- the reason
        # lambda2M/External is not a parallel-capable pair.
        "drives_shutter": True,
        "self_paced": False,
    },
}
