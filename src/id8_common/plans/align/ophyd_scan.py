"""
Ophyd-only scanning with live CSV output for the 8-ID beamlines.

This is ``scan_8id.py`` (Sam Marks' code) with Bluesky taken out. It is kept
deliberately close to the original -- same function names, same argument order
and defaults, same three detector branches in the same order, same nested
``step_*``/``inner_*`` helpers, same ``try: ... finally:`` bodies, same variable
names, same docstrings -- so anyone who knows ``scan_8id.py`` can read, edit and
maintain this file. ``scan_8id.py`` itself is untouched.

The translation rule is mechanical:

===============================================  ==============================
scan_8id.py                                      here
===============================================  ==============================
``yield from bps.mv(sig, x)``                    ``sig.put(x)``
``yield from bps.mv(sig1, a, sig2, b)``          one ``.put()`` per line
``yield from bps.mv(m1, s1, m2, s2)``            ``_return_motors([...])``
``yield from bps.move_per_step({m: p}, cache)``  ``m.move(p, wait=True)``
``yield from bps.sleep(t)``                      ``time.sleep(t)``
``create / read / read / ... / save``            ``scan.add_point(pos)``
``bpp.run_wrapper(..., md=md)``                  ``scan_csv.open_scan(...)``
``bpp.stage_wrapper(...)``                       gone -- nothing stages
===============================================  ==============================

A few multi-signal ``bps.mv()`` calls that Sam wrote once but that appear
identically in all six scans are called through one named helper here --
``_arm_eiger(det, num_pts, count_time)`` in place of his seven-signal
``bps.mv``. One call in, one call out, on the same line of the same branch: the
line-by-line correspondence is preserved and six copies become one. Nothing
else is factored out. There is no scan driver, no detector class, and the
position tables, ``step``/``inner`` bodies and ``try/finally`` skeletons stay
written out in each function, because that is what tells one scan from another.

Anything else that differs from ``scan_8id.py`` is marked with a
``# CHANGED (n):`` comment saying why. The full index:

 (1) ``gen_folder_prefix()`` is called exactly ONCE per scan. His
     ``save_images()`` calls it a second time, overwriting the caller's prefix:
     the ``.h5`` then gets a different name from its ``.csv`` and every scan
     burns two measurement numbers.
 (2) ``hdf1.enable`` / ``auto_save`` / ``file_write_mode`` are set for EVERY
     area detector, not just the eiger. ``bpp.stage_wrapper`` used to do this;
     nothing stages now. Without it a lambda scan writes no file at all and
     then spins forever in the frame wait -- and dmesh/mesh default to lambda2M.
 (3) ``setup_eiger_internal()`` is not called: it re-points hdf1 at a different
     folder, separating the ``.h5`` from its ``.csv``.
 (4) With ``save_img=0`` nothing captures, so ``num_captured`` never advances
     and his lambda frame wait spins forever. Wait out the exposure instead.
 (5) ``motor.stop()`` before the return move. See ``_return_motors``.
 (6) The eiger HDF drain is abort-aware -- it waits for the frames actually
     triggered, not for ``num_pts``, so an aborted scan does not sit out the
     whole timeout. Timeout is his ``num_pts * count_time + 10``.
 (7) ``blockbeam()`` is replaced by ``_blockbeam_verified()``, which confirms
     the blade against ``state_rbv``.
 (8) SIGINT is deferred for the duration of cleanup, so a second Ctrl+C cannot
     strand a moving axis with the shutter open. See ``_ScanGuard``.
 (9) ``bps.mv(m1, s1, m2, s2)`` moved both axes concurrently; two sequential
     ``move(wait=True)`` calls would not. See ``_return_motors``.
(10) ``time.sleep(count_time)`` in the eiger step of dmesh/mesh/d2scan/ascan/
     a2scan. Only his dscan had it; the other five read the ROI counters
     immediately after firing the trigger, i.e. before the frame was exposed,
     so every point recorded the PREVIOUS point's value. Fixing this makes
     those scans slower and moves their numbers -- see the note below.
(11) ``d2scan``'s lambda branch referenced an undefined ``frame_num``: a
     guaranteed ``NameError`` on the first point. It now counts frames the way
     dmesh does. ``ascan``/``a2scan`` had dropped the frame wait entirely.
(12) The tetramm arm/disarm is guarded by ``save_img``, matching his dmesh/mesh
     rather than his dscan (which armed HDF capture even with ``save_img=0``).
(13) The ``det.cam.stage_sigs`` save/patch/restore blocks are dropped. They are
     only ever applied by ``stage()``, and nothing stages without
     ``bpp.stage_wrapper``, so both halves were dead code.
(14) ``det=eiger4M`` becomes ``det=None`` + ``get_connected_device()``. A
     default argument is evaluated once, at import, so the old form captured
     ``None`` when ``safe_make_devices()`` had skipped that detector as
     offline, and failed later as an ``AttributeError`` on ``None``.
(15) The tetramm branch closes the shutter when it finishes, which his did not.
     Safe-direction, and it makes all three branches agree. It does NOT open
     the shutter -- that stays external, exactly as he had it.
(16) The public names are also exported with an ``_ophyd`` suffix. See the
     alias block at the bottom of this file for why that is not cosmetic.
(17) The return move waits for the axis to finish DECELERATING (DMOV back to 1)
     instead of sleeping a flat 0.2 s. Measured on hardware: with 0.2 s the
     record ignored the new target AND ophyd reported the move successful, so
     an aborted scan left the motor parked mid-scan while printing that it had
     gone home. See ``_settle``.

NOT PORTED, and why (all still present and working in ``scan_8id.py``):

  * ``dscan`` at :148 -- a dead twin of the one at :756 that shadows it.
  * ``dscan_test`` :1773 -- a ``bp.rel_scan`` wrapper with a debug print.
  * ``dscan_auto`` :1807 -- ``rate_limit`` is bound nowhere, so it raises
    ``NameError`` right after ``showbeam()``; its eiger branch never blocks the
    beam. Keep the idea (per-point attenuation) as a future ``att_map=`` kwarg.
  * the lup family (``x_lup``, ``y_lup``, ``huber_x_lup``, ``huber_y_lup``,
    ``rheo_x_lup``, ``rheo_y_lup``, ``rheo_set_x_lup``) :1500-1670 -- every one
    is a ``yield from bp.rel_scan(...)`` wrapper, i.e. a whole Bluesky plan
    driven by the RunEngine, with no line-by-line Ophyd equivalent. What they
    do is a relative scan of one axis, which is exactly ``dscan``:
    ``x_lup(-3, 3, 60)`` is ``dscan(sample.x, -3, 3, 60, <dwell>)``. Porting
    them meant inventing a dwell time ``bp.rel_scan`` never had, and made every
    alignment scan advance the shared measurement counter. Use ``dscan`` (or
    ``ascan``) directly. They remain available in a Bluesky session, where
    ``scan_8id.py`` still defines them.
  * the local ``att()`` :43 -- shadows ``shutter_att.att`` and raises
    ``TypeError`` when called with no argument.
  * ``from matplotlib.pylab import det`` :10 -- that binds ``numpy.linalg.det``
    as the bare name ``det`` in the session.

Each scan writes ONE .csv, named after that scan's .h5::

    .../data/bluesky/A0201_Test_a0010.h5     images
    .../data/bluesky/A0201_Test_a0010.csv    motor positions and counters

The NNNN in that name comes from ``gen_folder_prefix()``, which reads and then
advances ``expt.measurement_num`` -- the EPICS register ``8ideSoft:Reg1``, not a
per-session counter (see ``PV_FIELDS`` in ``expt_config.py`` for why it lives
there). ``det_acq_series()`` advances the SAME register, so scans and
acquisitions draw their numbers from one sequence and each moves the other's
numbering along. They do NOT share a folder, though: acquisitions write under
``data/`` and these scans under ``data/bluesky/``, so the highest number
already used may be under either one.

The .csv is closed after every point, so a viewer can poll it while the scan
runs. Two-motor scans (d2scan, a2scan, dmesh, mesh) write a column per motor,
so a mesh can be drawn as a map; see ``configs/scan_csv_template.yml``.

usage (call it directly -- do NOT wrap it in ``RE()``)::

    dscan(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
    dmesh(huber.nu, -0.5, 0.5, 11, huber.z, -0.5, 0.5, 11, 1.0, det=lambda2M)

Any positioner works -- ``huber.delta``, ``sample.x``, ``rheometer.y`` -- as
long as it has ``.name`` / ``.position`` / ``.move()`` / ``.stop()``.

This module resolves ``softglue`` and ``softglue_8id_acq`` from the
``oregistry`` at import time, so it must be imported after devices have been
loaded. The detectors are not -- they are looked up per call, see CHANGED (14).
"""

import os
import signal
import time

import numpy as np

from id8_common.expt_config import expt
from id8_common.plans.acquire.ad_acq import gen_folder_prefix
from id8_common.plans.align import scan_csv
from id8_common.plans.set.shutter_att import PIND_status
from id8_common.plans.set.shutter_att import att
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import pre_align
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.plans.set.shutter_att import shutteron
from id8_common.registry import get_connected_device
from id8_common.registry import oregistry

# softglue is present in every session, so it stays resolved at import time.
# The detectors are not: any of them can be skipped as offline by
# safe_make_devices(), so they are resolved per call by get_connected_device()
# -- see CHANGED (14). They are still available as bare names at the prompt,
# put there by the device loader, so `dscan(..., det=lambda2M)` works as before.
softglue = oregistry.get("softglue")
softglue_8id_acq = oregistry.get("softglue_8id_acq")

#: Which detector each scan reaches for when the caller names none. These are
#: Sam's defaults: his dscan/d2scan/ascan/a2scan default to eiger4M, his
#: dmesh/mesh to lambda2M.
DEFAULT_DETECTOR = "eiger4M"
DEFAULT_MESH_DETECTOR = "lambda2M"

#: Ceiling on a single return move, in seconds. Not optional: PositionerBase
#: builds its MoveStatus with self._timeout, which is None unless somebody set
#: it, so .wait() on an axis whose DMOV never clears would block forever -- and
#: with SIGINT deferred that is an unkillable session.
RETURN_TIMEOUT = 120.0

#: Longest to wait for a motor record to finish its deceleration ramp after a
#: stop(), before issuing the return move anyway. See :func:`_settle`.
SETTLE_TIMEOUT = 10.0


def data_folder():
    """``<mount_point><cycle>/<experiment>/data/bluesky``, created if needed.

    The leading tree comes from ``expt.mount_point`` and is NOT fixed to one
    station -- 8-ID-I and 8-ID-E have separate ones (``/gdata/dm/8ID/8IDI/`` and
    ``/gdata/dm/8ID/8IDE/``), and it has to match the station the current
    experiment belongs to.

    ``det_acq_series()`` writes under ``data/`` instead of ``data/bluesky``,
    but both name their files from the same measurement counter -- see the
    module docstring.

    Raises:
        RuntimeError: if the folder can't be created, naming the three
            ``configs/experiment.yml`` settings it was built from.
    """
    mount = expt.mount_point
    cycle = expt.cycle_name
    experiment = expt.experiment_name
    folder = f"{mount}{cycle}/{experiment}/data/bluesky"

    # Create it here rather than leave it to the first scan_csv write. These
    # three settings are set independently and drift apart -- a stale
    # mount_point names a tree that either doesn't exist or isn't ours, and DM
    # owns the cycle directories, so makedirs then fails with a bare
    # PermissionError from inside csv writing, several frames deep and naming
    # none of the settings responsible. Fail here instead, before the scan moves
    # anything, and say which one to look at.
    try:
        os.makedirs(folder, mode=0o770, exist_ok=True)
    except OSError as e:
        raise RuntimeError(
            f"Cannot use the scan data folder {folder!r}: {e.strerror}.\n"
            f"It is built from three settings in configs/experiment.yml -- "
            f"mount_point={mount!r}, cycle_name={cycle!r}, "
            f"experiment_name={experiment!r}.\n"
            f"Check that experiment {experiment!r} really lives under "
            f"{mount}{cycle}/ -- 8-ID-I and 8-ID-E have separate trees "
            f"(/gdata/dm/8ID/8IDI/ and /gdata/dm/8ID/8IDE/) and mount_point is "
            f"the usual thing left pointing at the previous experiment's one."
        ) from e

    return folder


def save_images(det, save_img, num_pts, num_frames=1, file_path=None, folder_prefix=None):
    """
    Toggle saving (1) or not saving (0) images.

    With save_img=0 this returns having touched nothing -- it is setup only, and
    does NOT disarm an HDF plugin that an earlier scan left armed.

    det: detector instance (eiger4M, lambda2M, or tetramm)
    save_img: 1 save, 0 don't save
    num_pts: how many frames the HDF plugin should expect (hdf1.num_capture)
    num_frames: number of frames to capture per point (default 1)
    file_path: override the folder (default: data_folder())
    folder_prefix: base file name for the .h5, from gen_folder_prefix().
        Required in practice when save_img=1: CHANGED (1) removed the fallback
        that used to fill it in, so None now raises TypeError on the print below.
    """
    if save_img not in (0, 1):
        raise ValueError("save_img must be 1 or 0 (to save or not to save)")

    if save_img == 1:
        # CHANGED (1): scan_8id.py calls gen_folder_prefix() again further down,
        # overwriting whatever the caller passed in. That is fatal here: the
        # caller has already named the .csv after its prefix, so the .h5 would
        # get a different name and the pair would no longer match. It also
        # burned two measurement numbers per scan. The `if folder_prefix is
        # None` guard he had here is deleted with it, so the only calls left in
        # this file are the six in the scans -- one each, and nowhere else.
        if file_path is None:
            file_path = data_folder()

        is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
        is_tetramm = "tetramm" in det.name.lower()

        def has(obj, attr):
            return getattr(obj, attr, None) is not None

        file_name = folder_prefix

        print("Scan folder created: " + folder_prefix)
        print("File path: ", file_path)

        if is_tetramm:
            det.hdf1.file_path.put(file_path)
            det.hdf1.file_name.put(file_name)
            det.hdf1.num_capture.put(num_pts)
            det.hdf1.file_write_mode.put(2)  # Stream mode
            return

        # for eiger4m, lambda2m
        if has(det, "cam"):
            if has(det.cam, "fw_enable"):
                det.cam.fw_enable.put(1)
            if has(det.cam, "save_files"):
                det.cam.save_files.put(1)

        if has(det, "hdf1"):
            if has(det.hdf1, "num_capture"):
                det.hdf1.num_capture.put(num_pts)
            if has(det.hdf1, "file_name"):
                det.hdf1.file_name.put(file_name)
            if has(det.hdf1, "file_path"):
                det.hdf1.file_path.put(file_path)

            # CHANGED (2): bpp.stage_wrapper used to enable the HDF plugin and
            # put it in Stream mode (apstools AD_EpicsHdf5FileName). Nothing
            # stages any more, so without these three the scan runs happily and
            # writes no file at all.
            #
            # These are OUTSIDE the `if is_eiger` below on purpose. In
            # scan_8id.py they only ever ran for the eiger, which means a
            # lambda2M scan with save_img=1 armed no HDF plugin, wrote nothing,
            # and then spun forever in the CHANGED (4) frame wait waiting for
            # num_captured on a plugin that was never capturing -- and dmesh and
            # mesh both DEFAULT to lambda2M.
            det.hdf1.enable.put(1)
            det.hdf1.auto_save.put(1)
            det.hdf1.file_write_mode.put(2)  # Stream

        if has(det, "cam"):
            if is_eiger:
                if has(det.cam, "trigger_mode"):
                    det.cam.trigger_mode.put("Internal Enable")
                if has(det.cam, "num_triggers"):
                    det.cam.num_triggers.put(1)
                if has(det.cam, "num_images"):
                    det.cam.num_images.put(1)

                # CHANGED (3): scan_8id.py calls setup_eiger_internal(...) here.
                # It re-points hdf1.file_path/file_name at get_common_file_path(),
                # a different folder, which separates the .h5 from its .csv and
                # breaks the viewer's "find the live scan" lookup. Everything
                # else it sets (acquire_time/period, trigger_mode, num_images,
                # num_triggers) the scan overwrites moments later anyway.
            else:
                if has(det.cam, "num_images"):
                    det.cam.num_images.put(num_frames)


# ===========================================================================
# ARMING -- each of these is ONE of Sam's bps.mv() calls. They are called from
# the same line of the same branch he wrote them on: one call in, one call out.
# ===========================================================================


def _arm_tetramm(det, save_img):
    """HDF capture on. His dmesh/mesh line; CHANGED (12) applies it everywhere.

    scan_8id.py's dscan armed capture even with save_img=0, which parks the
    plugin in Stream mode waiting for frames nobody asked to keep. His dmesh and
    mesh already guard this with save_img; that is the version used here.
    """
    if save_img == 1:
        det.hdf1.enable.put(1)
        det.hdf1.capture.put(1)
        print(f"TetrAMM HDF capture armed: {det.hdf1.file_name.get()}")


def _arm_lambda(det, num_pts, count_time):
    """The nine-signal bps.mv at scan_8id.py:202 / 424 / 634 / 810 / 1015 / 1200 / 1386."""
    det.cam.operating_mode.put(3)  # 24-bit dual threshold mode
    det.cam.trigger_mode.put("External_ImagePer")
    det.cam.acquire_time.put(count_time)
    det.cam.acquire_period.put(count_time)
    det.cam.num_images.put(num_pts)
    det.hdf1.num_capture.put(num_pts)
    softglue.num_triggers.put(1)  # one pulse per "1!"; else stale value from run_measurement
    softglue.acq_time.put(count_time)
    softglue.acq_period.put(count_time)


def _arm_eiger(det, num_pts, count_time):
    """The seven-signal bps.mv at scan_8id.py:265 / 483 / 693 / 873 / 1070 / 1252 / 1439.

    CHANGED (13): his cam.stage_sigs save/patch/restore block sat right after
    this and is dropped -- stage_sigs are only applied by stage(), and nothing
    stages without bpp.stage_wrapper, so both halves were dead code. The .put()
    calls here are what actually arms the detector.
    """
    det.cam.acquire_time.put(count_time)
    det.cam.acquire_period.put(count_time)
    det.cam.trigger_mode.put("Internal Series")
    det.cam.manual_trigger.put("Enable")
    det.cam.num_images.put(1)
    det.cam.num_triggers.put(num_pts)
    det.hdf1.num_capture.put(num_pts)
    # No hdf1.enable/hdf1.capture here, unlike _arm_tetramm() and the
    # inner_lambda() branches, and that is correct: arming the plugin is
    # save_images()' job and it only does it for save_img=1. A scan that is not
    # keeping images has no reason to arm the writer. Do not "even this up" with
    # the lambda path.


# ===========================================================================
# TEARDOWN -- one uninterrupted run from one of Sam's `finally:` blocks.
#
# One rule, and it decides the order everywhere below: everything after the
# verified shutter close is dose-free and may take as long as it needs;
# everything before it must be a non-blocking .put(). Each of these calls
# _enter_cleanup() on the way in, so Sam's finally blocks carry no interrupt
# bookkeeping of their own.
# ===========================================================================


def _stop_lambda(det):
    """softglue off, then acquire off, then capture off -- Sam's order, which matters.

    stop_pulses first: softglue drives the shutter during external-trigger
    acquisition, so a blockbeam() issued before the pulses stop is undone by the
    next pulse.
    """
    _enter_cleanup()
    _safe("stop softglue pulses", softglue.stop_pulses.put, "1!")
    _safe("stop lambda acquire", det.cam.acquire.put, 0)
    _safe("stop lambda capture", det.hdf1.capture.put, 0)


def _drain_eiger(det, num_pts, count_time, triggered=None):
    """Wait for the HDF plugin to flush, then stop acquiring.

    CHANGED (6): ``triggered`` is how many frames were actually fired. On a
    clean run that is num_pts and this behaves exactly as his loop did. On an
    abort it is however far the scan got, so cleanup does not sit out the full
    timeout waiting for frames that were never triggered -- with the shutter
    still open, because on this branch the close happens after the drain.

    Timeout is his dscan's ``num_pts * count_time + 10`` (his other five used
    ``* 3 + 10``; the tighter one is used everywhere, and it only ever bounds a
    wait that has already gone wrong).
    """
    _enter_cleanup()
    target = num_pts if triggered is None else triggered

    def wait_for_frames():
        t0 = time.time()
        timeout = num_pts * count_time + 10
        while det.hdf1.num_captured.get() < target:
            time.sleep(0.05)
            if time.time() - t0 > timeout:
                print("WARNING: HDF write timeout -- not all frames saved.")
                break

    _safe("drain eiger HDF", wait_for_frames)
    _safe("stop eiger acquire", det.cam.acquire.put, 0)


# ===========================================================================
# RESTORE -- Sam's outer `finally:` block, one call per detector.
# ===========================================================================


def _disarm_tetramm(det, save_img):
    """His `det.hdf1.capture.put(0)` plus the message. CHANGED (12)."""
    _enter_cleanup()
    if save_img == 1:
        _safe("stop tetramm capture", det.hdf1.capture.put, 0)
        print("TetrAMM HDF capture stopped.")
    # CHANGED (15): his tetramm branch left the shutter wherever it was. It
    # never opens it either, so this only ever closes something someone else
    # opened -- safe-direction, and it makes all three branches agree.
    _blockbeam_verified()


def _disarm_lambda(det):
    """His outer lambda finally: mode, trigger, preset, beam off, logic off, count."""
    _enter_cleanup()
    _safe("restore lambda operating_mode", det.cam.operating_mode.put, 3)
    _safe("restore lambda trigger_mode", det.cam.trigger_mode.put, 0)
    _safe("restore softglue preset", softglue_8id_acq.preset.put, 50)
    _blockbeam_verified()
    _safe("shutteroff", shutteroff)
    _safe("report frames", print, "# images captured: ", det.hdf1.num_captured.get())


def _disarm_eiger(det):
    """His outer eiger finally: back to Internal Enable, manual trigger off, beam off."""
    _enter_cleanup()
    _safe("restore eiger trigger_mode", det.cam.trigger_mode.put, "Internal Enable")
    _safe("restore eiger manual_trigger", det.cam.manual_trigger.put, "Disable")
    _blockbeam_verified()
    _safe("report frames", print, "# images captured: ", det.hdf1.num_captured.get())


# ===========================================================================
# CHANGED (8): SIGINT is deferred for the duration of cleanup.
#
# This is the one abstraction here with no counterpart in scan_8id.py, and it
# is not elegance -- ophyd is the proof that a bare try/finally cannot do the
# job. EpicsMotor.move() is
#
#     status = super().move(position, **kwargs)
#     self.user_setpoint.put(position, wait=False)     # <-- BEFORE the try
#     try:
#         if wait:
#             status_wait(status)
#     except KeyboardInterrupt:
#         self.stop()
#         raise                                        # <-- re-raises at you
#
# so the slowest and LAST step of cleanup -- sending the motor home -- is a call
# ophyd deliberately aborts and re-raises out of on a second ^C, skipping
# everything after it. On the eiger branch that leaves the shutter open. And the
# setpoint is commanded one line BEFORE that try, so a ^C landing in that window
# leaves the axis accelerating with no stop() issued at all. No arrangement of
# try/finally fixes this, because the problem is that the cleanup code is itself
# interruptible.
# ===========================================================================

# How the guard behaves is kept in plain module-level variables rather than in
# an object: exactly one scan runs at a time at the prompt, and every variable
# here is read and written only by the handful of small functions below.

#: How many scans have installed the guard. Normally 0 or 1. The counter only
#: matters if a scan is ever started from inside another one, so that the inner
#: one does not hand SIGINT back while the outer one is still running.
_guard_depth = 0

#: True once we have actually taken SIGINT over (see _install_guard).
_guard_installed = False

#: What SIGINT did before we took it over, so it can be put back exactly.
_guard_previous_handler = None

#: False while the scan is running -- a ^C raises KeyboardInterrupt, exactly as
#: it always has. True once hardware cleanup has started -- a ^C is then counted
#: and reported instead of obeyed.
_guard_in_cleanup = False

#: How many times ^C has been pressed since cleanup started, and when the first
#: of those presses arrived.
_guard_ctrl_c_count = 0
_guard_first_ctrl_c_time = 0.0

#: The escape hatch needs BOTH of these: three presses AND five seconds.
#: A count on its own is wrong -- a held-down ^C repeats about thirty times a
#: second, so three presses arrive within a tenth of a second and would abandon
#: a cleanup that was about to succeed. A timer on its own is wrong -- one
#: stray press should arm nothing.
_GUARD_PRESSES_TO_GIVE_UP = 3
_GUARD_SECONDS_TO_GIVE_UP = 5.0

#: Written with os.write, never print(). ScanCsv.add_point() print()s on EVERY
#: scan point, and a signal handler that re-enters the buffered-writer lock held
#: by the line it interrupted deadlocks -- in the one code path whose job is to
#: shut the shutter. A plain bytes constant, so there is no f-string to build,
#: no logging lock to take and no EPICS call to make from inside a handler.
_GUARD_MESSAGE = (
    b"\n^C  cleaning up: beam off, motors home. Hold on.\n"
    b"    (^C 3 more times over 5 s abandons the return move -- "
    b"the beam is already off by then.)\n"
)


def _handle_ctrl_c(signum, frame):
    """What Ctrl+C does while a scan owns it."""
    global _guard_ctrl_c_count, _guard_first_ctrl_c_time

    if not _guard_in_cleanup:
        # Still scanning. Do exactly what Python normally does, so that every
        # `except KeyboardInterrupt:` in this file keeps working unchanged.
        signal.default_int_handler(signum, frame)
        return

    # Cleaning up. Count the press, say so on screen, and carry on.
    _guard_ctrl_c_count += 1
    if _guard_ctrl_c_count == 1:
        _guard_first_ctrl_c_time = time.time()
    os.write(2, _GUARD_MESSAGE)


def _install_guard():
    """Take Ctrl+C over for one scan. ALWAYS pair with _restore_guard().

    Called before the scan touches any hardware, and put back in a `finally:`
    so that it is handed back even when the scan raises. Installing it later --
    inside the `except KeyboardInterrupt:` block, say -- would leave a window in
    which a second ^C is still delivered to the old handler and escapes.
    """
    global _guard_depth, _guard_installed, _guard_previous_handler
    global _guard_in_cleanup, _guard_ctrl_c_count, _guard_first_ctrl_c_time

    _guard_depth += 1
    if _guard_depth > 1:
        return  # an outer scan already owns SIGINT; leave its handler alone

    _guard_in_cleanup = False
    _guard_ctrl_c_count = 0
    _guard_first_ctrl_c_time = 0.0
    _guard_previous_handler = signal.getsignal(signal.SIGINT)
    try:
        signal.signal(signal.SIGINT, _handle_ctrl_c)
        _guard_installed = True
    except ValueError:
        # Not the main thread -- a queueserver worker, say. Python only delivers
        # ^C to the main thread at all, so such a scan simply runs to
        # completion; _safe() below still protects each cleanup step.
        _guard_installed = False


def _restore_guard():
    """Give Ctrl+C back to whoever had it. Safe even if we never took it."""
    global _guard_depth, _guard_installed, _guard_in_cleanup

    _guard_depth -= 1
    if _guard_depth > 0:
        return  # an outer scan still owns it
    _guard_depth = 0
    _guard_in_cleanup = False

    if not _guard_installed:
        return
    _guard_installed = False

    previous = _guard_previous_handler
    if previous is None:
        # getsignal() returns None when the handler it replaced came from C, and
        # None is not a legal argument to signal.signal(). Fall back to Python's
        # own handler rather than die on the last line of cleanup.
        previous = signal.default_int_handler
    try:
        signal.signal(signal.SIGINT, previous)
    except ValueError:
        pass


def _enter_cleanup():
    """Hardware cleanup starts here: stop raising on ^C, start counting instead.

    Called at the top of every teardown helper, so that the scan functions
    themselves carry no interrupt bookkeeping.
    """
    global _guard_in_cleanup
    _guard_in_cleanup = True


def _user_gave_up():
    """True once ^C has been pressed enough times, over enough seconds, to mean it.

    Only _return_motors() asks. Every other cleanup step is a quick .put() that
    is not worth abandoning; the return move is the one slow step, and by the
    time it runs the beam is already off.
    """
    if _guard_ctrl_c_count < _GUARD_PRESSES_TO_GIVE_UP:
        return False
    return time.time() - _guard_first_ctrl_c_time >= _GUARD_SECONDS_TO_GIVE_UP


def _safe(what, action, *args, **kwargs):
    """Run one cleanup step. Never let it stop the steps that follow.

    `what` is a short description used only in the warning if it fails::

        _safe("stop lambda acquire", det.cam.acquire.put, 0)

    calls ``det.cam.acquire.put(0)``. A step that raises is reported and skipped,
    so one dead PV cannot leave the shutter open or the motor stranded.
    """
    try:
        action(*args, **kwargs)
    except KeyboardInterrupt:
        # Only reachable when _install_guard() could not install (a worker
        # thread). The scan's own `except KeyboardInterrupt` already recorded
        # the abort, so just note it and keep cleaning up.
        os.write(2, b"\n^C  still cleaning up.\n")
    except Exception as exc:
        # Deliberately not BaseException: SystemExit must still exit.
        print(f"WARNING: cleanup step {what} failed: {exc!r}", flush=True)


def _blockbeam_verified(timeout=2.0):
    """blockbeam(), then make sure the blade really closed.  CHANGED (7).

    `operation` is only the request; `state_rbv` is where the blade actually is,
    and FastShutter documents in its own source that the two disagree while
    softglue holds the override -- which is every external-trigger acquisition.
    A pyepics put to a dropped IOC is also a silent no-op. A cleanup that prints
    "beam blocked" over an open shutter is the single failure this whole design
    exists to prevent, so this one step is verified and escalated. Nothing else
    is.

    Bounded at 2 s: this runs BEFORE the escape hatch is reachable, so it must
    not be able to trap the user. Poll first and escalate ONCE -- hammering
    shutteroff()/blockbeam() every 100 ms would drop and re-request the softglue
    lock twenty times on the safety path.
    """
    _enter_cleanup()
    closed = False

    def close_and_confirm():
        nonlocal closed
        shutter = get_connected_device("shutter_8ide")
        blockbeam()
        t0 = time.time()
        escalated = False
        while time.time() - t0 < timeout:
            if str(shutter.state_rbv.get(as_string=True)).lower().startswith("clos"):
                closed = True
                return
            if not escalated and time.time() - t0 > timeout / 2:
                shutteroff()  # drop softglue's override, then ask once more
                blockbeam()
                escalated = True
            time.sleep(0.05)
        print(
            "\033[91m*** SHUTTER DID NOT CONFIRM CLOSED -- CLOSE IT BY HAND ***\033[0m",
            flush=True,
        )

    _safe("block beam", close_and_confirm)
    return closed


def _settle(motor, timeout=SETTLE_TIMEOUT):
    """Wait for an EPICS motor record to finish decelerating (DMOV back to 1).

    CHANGED (17). This is not optional, and a fixed sleep does not do it.
    Measured on huber.nu (8ideSoft:CR8-E1:m4, VELO 0.4 deg/s) on 2026-09-06,
    mid-move, with an abort:

        stop(); sleep(0.2); move(home)  ->  DMOV was still 0 when the move was
            issued. status.wait() returned success=True. The motor never moved:
            it sat at the scan point, and was still there a second later.
        stop(); sleep(1.0); move(home)  ->  DMOV was 1 when the move was issued.
            The motor went home.

    Two things go wrong at once when the ramp has not finished. The record will
    not honour a new target; AND ophyd's MoveStatus completes on the DMOV 0->1
    edge produced by the DECELERATION rather than by the new move, so .wait()
    reports success for a move that never happened. A cleanup that prints
    "back at start" over a motor parked mid-scan is precisely the failure this
    module exists to prevent, so wait for the real thing instead of guessing a
    sleep that is right for one axis and wrong for the next.

    A SoftPositioner has no DMOV and nothing to settle, so it returns True.
    """
    done = getattr(motor, "motor_done_move", None)
    if done is None:
        return True
    t0 = time.time()
    while time.time() - t0 < timeout:
        try:
            if int(done.get()) == 1:
                return True
        except Exception:
            return True  # can't read DMOV -- don't block cleanup on it
        time.sleep(0.05)
    return False


def _return_motors(pairs, timeout=RETURN_TIMEOUT):
    """Send every motor home at once -- what bps.mv(m1, s1, m2, s2) did.  CHANGED (9).

    Two sequential move(wait=True) calls would double the return time and change
    the path through space on d2scan/dmesh/mesh, so this is the translation-table
    entry for a multi-signal bps.mv, not a consolidation. Started non-blocking so
    the axes travel together and one dead axis cannot hold up the others; waited
    separately so one failure is reported and the rest still finish.

    After a stop() the status finishes UnknownStatusFailure, so .wait() RAISING
    is the normal outcome for an axis the user escaped out of, not a bug.
    """
    _enter_cleanup()
    started = []
    for motor, start in pairs:
        if _user_gave_up():
            print(
                f"^C^C^C  leaving {motor.name} at {motor.position} (start was {start}). "
                f"The beam is off. Move it back at the prompt.",
                flush=True,
            )
            continue
        # CHANGED (5): stop() before the return move. NOT because of "another
        # set() is still in progress" -- that text is Signal.set and is
        # unreachable from EpicsMotor.move, whose PositionerBase.move() finishes
        # any stale status on its first two lines. The three real reasons:
        # (a) EpicsMotor commands user_setpoint.put() one line BEFORE the try
        #     that catches ^C, so only our own stop() halts an axis interrupted
        #     there; (b) a soft-positioner-backed axis does not self-stop, and
        #     there that message IS what you hit; (c) the EPICS motor record
        #     needs its deceleration ramp to finish before it will honour a new
        #     target -- which is what the _settle() calls below wait for. (They
        #     replaced a flat 0.2 s sleep; see CHANGED (17).)
        # EpicsMotor.stop carries @raise_if_disconnected, hence the wrapper.
        _safe(f"stop {motor.name}", motor.stop)
    # CHANGED (17), replacing a flat time.sleep(0.2): wait for each axis to
    # actually finish decelerating. See _settle() for the measurement that
    # showed 0.2 s silently loses the return move on huber.nu.
    for motor, _start in pairs:
        if _user_gave_up():
            continue

        def settle_and_warn(m=motor):
            if not _settle(m):
                print(
                    f"WARNING: {m.name} still reports DMOV=0 after {SETTLE_TIMEOUT}s; "
                    f"issuing the return move anyway, but it may not take -- check the axis.",
                    flush=True,
                )

        _safe(f"settle {motor.name}", settle_and_warn)
    for motor, start in pairs:
        if _user_gave_up():
            continue
        try:
            started.append((motor, start, motor.move(start, wait=False, timeout=timeout)))
        except Exception as exc:
            print(f"WARNING: {motor.name} would not accept the return move: {exc!r}", flush=True)
    for motor, start, status in started:
        _safe(f"return {motor.name}", status.wait, timeout=timeout + 5)
        # The status alone is not proof -- see _settle(). Let the ramp finish so
        # the position printed below is the one the axis actually stopped at.
        _safe(f"settle {motor.name} after return", _settle, motor)

        # Both numbers, no tolerance test: there is no epsilon that is right for
        # huber.delta and sample.x alike, and a wrong hard-coded one trains
        # people to ignore cleanup warnings.
        def report(m=motor, s=start):
            print(f"{m.name}: back at {m.position} (start {s})", flush=True)

        _safe(f"report {motor.name}", report)


def _move_motors(pairs):
    """Move every motor to its point at once -- what bps.move_per_step() did.

    ``bps.move_per_step({motor1: p1, motor2: p2}, cache)`` expands to a single
    ``mv()`` over both axes, so they travel together. Two sequential
    ``move(wait=True)`` calls would instead drive one axis and then the other,
    which on a d2scan/a2scan is a visibly different path through space -- an L
    instead of a diagonal. Not a consolidation: this is the translation-table
    entry for move_per_step with more than one motor.

    Unlike _return_motors this is NOT wrapped in _safe(): a motor that will not
    move during the scan is a real failure and must abort it, not be logged and
    skipped. It is only forgiving during cleanup.
    """
    statuses = [motor.move(position, wait=False) for motor, position in pairs]
    for status in statuses:
        status.wait()


# ===========================================================================
# The scans, in the order they appear in scan_8id.py.
# ===========================================================================


def dmesh(
    motor1,
    rel_begin1,
    rel_end1,
    num1,
    motor2,
    rel_begin2,
    rel_end2,
    num2,
    count_time,
    det=None,
    att_ratio=7,
    save_img=1,
    comment="",
):
    """
    Relative 2D raster (mesh) scan with per-step triggering for lambda2M, eiger4M, and tetramm.

    Begin/end values are relative to each motor's current position. motor1 is the
    outer (slow) axis, motor2 is the inner (fast) axis. For every position of
    motor1 the plan steps motor2 across its full range, counting the detector at
    each of the num1 * num2 grid points. Both motors return to their starting
    positions when the scan finishes.

    usage:
        dmesh(huber.x, -0.5, 0.5, 11, huber.z, -0.5, 0.5, 11, 1)

    args:
        motor1: outer/slow ophyd positioner (e.g. huber.x)
        rel_begin1, rel_end1: start/end for motor1 relative to its current position
        num1: number of points along motor1
        motor2: inner/fast ophyd positioner (e.g. huber.z)
        rel_begin2, rel_end2: start/end for motor2 relative to its current position
        num2: number of points along motor2
        count_time: detector acquisition time per point (s)
        det: detector (lambda2M, eiger4M, or tetramm1); default lambda2M
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
        comment: free-text note written near the top of the .csv

    returns:
        the :class:`scan_csv.ScanCsv` that was written.
    """
    if det is None:
        det = get_connected_device(DEFAULT_MESH_DETECTOR)  # CHANGED (14)

    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name.lower()
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

    num_pts = num1 * num2
    # One measurement number per scan whether or not images are saved, so every
    # .csv has a unique name that matches its .h5 when there is one. His
    # `gen_folder_prefix() if save_img == 1 else ""` left every save_img=0 scan
    # trying to write to the same nameless file.
    folder_prefix = gen_folder_prefix()
    file_path = data_folder()
    h5_file = f"{file_path}/{folder_prefix}.h5" if save_img == 1 else ""
    csv_file = f"{file_path}/{folder_prefix}.csv"
    if save_img != 1:
        print(f"save_img=0: no .h5 will be written; scan file is {folder_prefix}.csv")

    save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    # This is what run_wrapper(md=...) used to do. Resolved before anything
    # moves, so a typo in the template raises here and not halfway through a
    # scan with the beam on.
    command = (
        f"dmesh({motor1.name}, {rel_begin1}, {rel_end1}, {num1}, "
        f"{motor2.name}, {rel_begin2}, {rel_end2}, {num2}, {count_time}, det={det.name})"
    )
    start1 = motor1.position
    start2 = motor2.position
    positions1 = np.linspace(start1 + rel_begin1, start1 + rel_end1, num1)
    positions2 = np.linspace(start2 + rel_begin2, start2 + rel_end2, num2)

    scan = scan_csv.open_scan(
        csv_file,
        [motor1, motor2],
        det=det,
        scan_type="dmesh",
        command=command,
        num_points=num_pts,
        count_time=count_time,
        h5_file=h5_file,
        comment=comment,
        shape=(num1, num2),
        # The COMMANDED first and last position of each axis, so a viewer can
        # rebuild the exact grid instead of inferring it from the readbacks it
        # happens to have. Inferring is not good enough: a real motor's readback
        # differs at every visit to the same cell, and a scan stopped early
        # spans less than it was told to -- reconstructing from the measured
        # span then puts the measured rows in the WRONG cells, which looks
        # plausible and is wrong. Written before anything moves.
        extra={
            "motor1_start": repr(float(positions1[0])),
            "motor1_stop": repr(float(positions1[-1])),
            "motor2_start": repr(float(positions2[0])),
            "motor2_stop": repr(float(positions2[-1])),
        },
    )
    scan.write_header()
    expt.file_name = folder_prefix

    # CHANGED (8). Installed before anything moves, and handed back in the
    # finally below so it is given up even if the scan raises.
    _install_guard()
    try:
        if is_tetramm:
            _arm_tetramm(det, save_img)

            def inner_tetramm():
                try:
                    for p1 in positions1:
                        for p2 in positions2:
                            _move_motors([(motor1, p1), (motor2, p2)])
                            # The TetrAMM free-runs; trigger() only pokes Acquire
                            # and stamps a datum, so this is what
                            # trigger_and_read did.
                            det.acquire.put(1, wait=False)
                            time.sleep(count_time)
                            scan.add_point(p1, p2)
                finally:
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_tetramm()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_tetramm(det, save_img)
            return scan

        if is_lambda:
            _arm_lambda(det, num_pts, count_time)

            def step_lambda(p1, p2, frame_num):
                _move_motors([(motor1, p1), (motor2, p2)])
                softglue.start_pulses.put("1!")
                if save_img == 1:
                    # wait for this frame to be captured before reading stats
                    while det.hdf1.num_captured.get() < frame_num:
                        time.sleep(0.005)
                else:
                    # CHANGED (4): with save_img=0 nothing is capturing, so
                    # num_captured never advances and his wait above spins for
                    # ever. Wait out the exposure instead.
                    time.sleep(count_time)
                scan.add_point(p1, p2)  # was: create / read motors+stats1..3 / save

            def inner_lambda():
                det.cam.acquire.put(1)
                det.hdf1.capture.put(1)
                shutteron()
                showbeam()
                frame_num = 0
                try:
                    for p1 in positions1:
                        for p2 in positions2:
                            frame_num += 1
                            step_lambda(p1, p2, frame_num)
                finally:
                    _stop_lambda(det)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_lambda()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_lambda(det)
            return scan

        # eiger4M
        if is_eiger:
            _arm_eiger(det, num_pts, count_time)

            fired = 0

            def step(detectors, p1, p2):
                """Move motors, fire software trigger, wait for frame, read."""
                nonlocal fired
                _move_motors([(motor1, p1), (motor2, p2)])
                det.cam.special_trigger_button.put(1, wait=False)
                fired += 1
                # CHANGED (10): his dmesh read the ROI counters immediately after
                # firing the trigger, i.e. before the frame had been exposed, so
                # every point recorded the PREVIOUS point's value (and point 0 a
                # pre-scan one). Only his dscan waited. This is his dscan's line.
                time.sleep(count_time)
                scan.add_point(p1, p2)  # was: create / read motors+stats1..3 / save

            def inner():
                """Start pre-armed acquisition (accepts num_pts software triggers)."""
                showbeam()
                det.cam.acquire.put(1, wait=False)
                try:
                    for p1 in positions1:
                        for p2 in positions2:
                            step([det], p1, p2)
                finally:
                    _drain_eiger(det, num_pts, count_time, triggered=fired)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_eiger(det)
            return scan

        scan.close("error")
        raise ValueError(f"Unrecognized detector {det.name!r} (expected eiger, lambda or tetramm)")
    finally:
        _restore_guard()


def mesh(
    motor1,
    abs_begin1,
    abs_end1,
    num1,
    motor2,
    abs_begin2,
    abs_end2,
    num2,
    count_time,
    det=None,
    att_ratio=7,
    save_img=1,
    comment="",
):
    """
    Absolute 2D raster (mesh) scan with per-step triggering for lambda2M, eiger4M, and tetramm.

    motor1 is the outer (slow) axis, motor2 is the inner (fast) axis. For every position of
    motor1 the plan steps motor2 across its full range, counting the detector at
    each of the num1 * num2 grid points. Both motors return to their starting
    positions when the scan finishes.

    usage:
        mesh(huber.x, -0.5, 0.5, 11, huber.z, -0.5, 0.5, 11, 1)

    args:
        motor1: outer/slow ophyd positioner (e.g. huber.x)
        abs_begin1, abs_end1: absolute start/end for motor1
        num1: number of points along motor1
        motor2: inner/fast ophyd positioner (e.g. huber.z)
        abs_begin2, abs_end2: absolute start/end for motor2
        num2: number of points along motor2
        count_time: detector acquisition time per point (s)
        det: detector (lambda2M, eiger4M, or tetramm1); default lambda2M
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
        comment: free-text note written near the top of the .csv

    returns:
        the :class:`scan_csv.ScanCsv` that was written.
    """
    if det is None:
        det = get_connected_device(DEFAULT_MESH_DETECTOR)  # CHANGED (14)

    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name.lower()
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

    num_pts = num1 * num2
    # One measurement number per scan whether or not images are saved, so every
    # .csv has a unique name that matches its .h5 when there is one. His
    # `gen_folder_prefix() if save_img == 1 else ""` left every save_img=0 scan
    # trying to write to the same nameless file.
    folder_prefix = gen_folder_prefix()
    file_path = data_folder()
    h5_file = f"{file_path}/{folder_prefix}.h5" if save_img == 1 else ""
    csv_file = f"{file_path}/{folder_prefix}.csv"
    if save_img != 1:
        print(f"save_img=0: no .h5 will be written; scan file is {folder_prefix}.csv")

    save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    # This is what run_wrapper(md=...) used to do. Resolved before anything
    # moves, so a typo in the template raises here and not halfway through a
    # scan with the beam on.
    command = (
        f"mesh({motor1.name}, {abs_begin1}, {abs_end1}, {num1}, "
        f"{motor2.name}, {abs_begin2}, {abs_end2}, {num2}, {count_time}, det={det.name})"
    )
    start1 = motor1.position
    start2 = motor2.position
    positions1 = np.linspace(abs_begin1, abs_end1, num1)
    positions2 = np.linspace(abs_begin2, abs_end2, num2)

    scan = scan_csv.open_scan(
        csv_file,
        [motor1, motor2],
        det=det,
        scan_type="mesh",
        command=command,
        num_points=num_pts,
        count_time=count_time,
        h5_file=h5_file,
        comment=comment,
        shape=(num1, num2),
        # The COMMANDED first and last position of each axis, so a viewer can
        # rebuild the exact grid instead of inferring it from the readbacks it
        # happens to have. Inferring is not good enough: a real motor's readback
        # differs at every visit to the same cell, and a scan stopped early
        # spans less than it was told to -- reconstructing from the measured
        # span then puts the measured rows in the WRONG cells, which looks
        # plausible and is wrong. Written before anything moves.
        extra={
            "motor1_start": repr(float(positions1[0])),
            "motor1_stop": repr(float(positions1[-1])),
            "motor2_start": repr(float(positions2[0])),
            "motor2_stop": repr(float(positions2[-1])),
        },
    )
    scan.write_header()
    expt.file_name = folder_prefix

    # CHANGED (8). Installed before anything moves, and handed back in the
    # finally below so it is given up even if the scan raises.
    _install_guard()
    try:
        if is_tetramm:
            _arm_tetramm(det, save_img)

            def inner_tetramm():
                try:
                    for p1 in positions1:
                        for p2 in positions2:
                            _move_motors([(motor1, p1), (motor2, p2)])
                            # The TetrAMM free-runs; trigger() only pokes Acquire
                            # and stamps a datum, so this is what
                            # trigger_and_read did.
                            det.acquire.put(1, wait=False)
                            time.sleep(count_time)
                            scan.add_point(p1, p2)
                finally:
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_tetramm()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_tetramm(det, save_img)
            return scan

        if is_lambda:
            _arm_lambda(det, num_pts, count_time)

            def step_lambda(p1, p2, frame_num):
                _move_motors([(motor1, p1), (motor2, p2)])
                softglue.start_pulses.put("1!")
                if save_img == 1:
                    # wait for this frame to be captured before reading stats
                    while det.hdf1.num_captured.get() < frame_num:
                        time.sleep(0.005)
                else:
                    # CHANGED (4): with save_img=0 nothing is capturing, so
                    # num_captured never advances and his wait above spins for
                    # ever. Wait out the exposure instead.
                    time.sleep(count_time)
                scan.add_point(p1, p2)  # was: create / read motors+stats1..3 / save

            def inner_lambda():
                det.cam.acquire.put(1)
                det.hdf1.capture.put(1)
                shutteron()
                showbeam()
                frame_num = 0
                try:
                    for p1 in positions1:
                        for p2 in positions2:
                            frame_num += 1
                            step_lambda(p1, p2, frame_num)
                finally:
                    _stop_lambda(det)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_lambda()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_lambda(det)
            return scan

        # eiger4M
        if is_eiger:
            _arm_eiger(det, num_pts, count_time)

            fired = 0

            def step(detectors, p1, p2):
                """Move motors, fire software trigger, wait for frame, read."""
                nonlocal fired
                _move_motors([(motor1, p1), (motor2, p2)])
                det.cam.special_trigger_button.put(1, wait=False)
                fired += 1
                # CHANGED (10): his dmesh read the ROI counters immediately after
                # firing the trigger, i.e. before the frame had been exposed, so
                # every point recorded the PREVIOUS point's value (and point 0 a
                # pre-scan one). Only his dscan waited. This is his dscan's line.
                time.sleep(count_time)
                scan.add_point(p1, p2)  # was: create / read motors+stats1..3 / save

            def inner():
                """Start pre-armed acquisition (accepts num_pts software triggers)."""
                showbeam()
                det.cam.acquire.put(1, wait=False)
                try:
                    for p1 in positions1:
                        for p2 in positions2:
                            step([det], p1, p2)
                finally:
                    _drain_eiger(det, num_pts, count_time, triggered=fired)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_eiger(det)
            return scan

        scan.close("error")
        raise ValueError(f"Unrecognized detector {det.name!r} (expected eiger, lambda or tetramm)")
    finally:
        _restore_guard()


def dscan(motor, rel_begin, rel_end, num_pts, count_time, det=None, att_ratio=1e6, save_img=1, comment=""):
    """
    Pre-armed software-trigger scan for fast acquisitions.

    args:
        motor: ophyd positioner -- huber.delta, sample.x, rheometer.y, ...
        rel_begin, rel_end: relative start/end (motor units)
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (eiger4M, lambda2M, or tetramm1); default eiger4M
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
        comment: free-text note written near the top of the .csv, so you can
            grep the data folder for a keyword and land on the right scan.
            Newlines are flattened to spaces; empty means no comment line.

    returns:
        the :class:`scan_csv.ScanCsv` that was written -- ``.path`` is the file,
        ``.data`` is ``{column_label: [values]}``.

    On Ctrl+C the motor stops and returns to where it started, acquisition
    halts, the beam is blocked, and the .csv closes with ``#exit_status,aborted``.
    Repeated Ctrl+C during that cleanup is absorbed, not obeyed.
    """
    if det is None:
        det = get_connected_device(DEFAULT_DETECTOR)  # CHANGED (14)

    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name.lower()
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

    # One measurement number per scan whether or not images are saved, so every
    # .csv has a unique name that matches its .h5 when there is one.
    folder_prefix = gen_folder_prefix()
    file_path = data_folder()
    h5_file = f"{file_path}/{folder_prefix}.h5" if save_img == 1 else ""
    csv_file = f"{file_path}/{folder_prefix}.csv"
    if save_img != 1:
        print(f"save_img=0: no .h5 will be written; scan file is {folder_prefix}.csv")

    save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    # This is what run_wrapper(md=...) used to do. Resolved before anything
    # moves, so a typo in the template raises here and not halfway through a
    # scan with the beam on.
    command = f"dscan({motor.name}, {rel_begin}, {rel_end}, {num_pts}, {count_time}, det={det.name})"
    scan = scan_csv.open_scan(
        csv_file,
        motor,
        det=det,
        scan_type="dscan",
        command=command,
        num_points=num_pts,
        count_time=count_time,
        h5_file=h5_file,
        comment=comment,
    )
    scan.write_header()
    # One bare name for both files -- the .h5 and .csv differ only by extension,
    # and the folder is rebuildable from the other settings.
    expt.file_name = folder_prefix

    # CHANGED (8). Installed before anything moves, and handed back in the
    # finally below so it is given up even if the scan raises.
    _install_guard()
    try:
        if is_tetramm:
            _arm_tetramm(det, save_img)

            start_pos = motor.position
            positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)

            def inner_tetramm():
                try:
                    for pos in positions:
                        motor.move(pos, wait=True)
                        # The TetrAMM free-runs; trigger() only pokes Acquire and
                        # stamps a datum, so this is what trigger_and_read did.
                        det.acquire.put(1, wait=False)
                        time.sleep(count_time)
                        scan.add_point(pos)
                finally:
                    _return_motors([(motor, start_pos)])

            exit_status = "success"
            try:
                inner_tetramm()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_tetramm(det, save_img)
            return scan

        if is_lambda:
            _arm_lambda(det, num_pts, count_time)

            start_pos = motor.position
            positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)

            def step_lambda(detectors, pos, frame_num):
                motor.move(pos, wait=True)
                softglue.start_pulses.put("1!")
                if save_img == 1:
                    while det.hdf1.num_captured.get() < frame_num:
                        time.sleep(0.005)
                else:
                    # CHANGED (4): with save_img=0 nothing is capturing, so
                    # num_captured never advances and his wait above spins for
                    # ever. Wait out the exposure instead.
                    time.sleep(count_time)
                scan.add_point(pos)  # was: create / read motor+stats1..3 / save

            def inner_lambda():
                det.cam.acquire.put(1)
                det.hdf1.capture.put(1)
                shutteron()
                showbeam()
                try:
                    for ii, pos in enumerate(positions):
                        step_lambda([det], pos, ii + 1)
                finally:
                    _stop_lambda(det)
                    _blockbeam_verified()
                    _return_motors([(motor, start_pos)])

            exit_status = "success"
            try:
                inner_lambda()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_lambda(det)
            return scan

        # eiger4M
        if is_eiger:
            _arm_eiger(det, num_pts, count_time)

            start_pos = motor.position
            positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)

            fired = 0

            def step(detectors, pos):
                """Move motor, fire software trigger, wait for frame, read motor."""
                nonlocal fired
                motor.move(pos, wait=True)
                det.cam.special_trigger_button.put(1, wait=False)
                fired += 1
                time.sleep(count_time)
                scan.add_point(pos)  # was: create / read motor+stats1..3 / save

            def inner():
                """Start pre-armed acquisition (accepts num_pts software triggers)"""
                det.cam.acquire.put(1, wait=False)
                showbeam()
                try:
                    for pos in positions:
                        step([det], pos)
                finally:
                    _drain_eiger(det, num_pts, count_time, triggered=fired)
                    # His dscan blocked the beam in the OUTER finally, i.e. after
                    # the return move; his other five block here, before it.
                    # Theirs is used -- the return move is the long step and
                    # there is no reason to spend it with the shutter open.
                    _blockbeam_verified()
                    _return_motors([(motor, start_pos)])

            exit_status = "success"
            try:
                inner()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                # Return detector to normal state: Internal Enable, manual trigger off
                _disarm_eiger(det)
            return scan

        scan.close("error")
        raise ValueError(f"Unrecognized detector {det.name!r} (expected eiger, lambda or tetramm)")
    finally:
        _restore_guard()


def d2scan(
    motor1,
    rel_begin1,
    rel_end1,
    motor2,
    rel_begin2,
    rel_end2,
    num_pts,
    count_time,
    det=None,
    att_ratio=7,
    save_img=1,
    comment="",
):
    """
    Two-motor relative scan with per-step triggering for eiger4M, lambda2M, and tetramm.

    Both motors sweep together along a single line: point i puts motor1 at
    positions1[i] AND motor2 at positions2[i]. This is a 1D trajectory, not a
    raster -- use dmesh/mesh for a grid, and plot the result in the 1D viewer
    against whichever motor you like.

    usage:
        d2scan(m1, s1, f1, m2, s2, f2, npts, ct)

    args:
        motor1, motor2: ophyd positioners (e.g. huber.eta, huber.chi)
        rel_begin1, rel_end1: relative start/end for motor1
        rel_begin2, rel_end2: relative start/end for motor2
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (eiger4M, lambda2M, or tetramm1); default eiger4M
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
        comment: free-text note written near the top of the .csv

    returns:
        the :class:`scan_csv.ScanCsv` that was written.
    """
    if det is None:
        det = get_connected_device(DEFAULT_DETECTOR)  # CHANGED (14)

    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name.lower()
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

    folder_prefix = gen_folder_prefix()
    file_path = data_folder()
    h5_file = f"{file_path}/{folder_prefix}.h5" if save_img == 1 else ""
    csv_file = f"{file_path}/{folder_prefix}.csv"
    if save_img != 1:
        print(f"save_img=0: no .h5 will be written; scan file is {folder_prefix}.csv")

    save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    command = (
        f"d2scan({motor1.name}, {rel_begin1}, {rel_end1}, "
        f"{motor2.name}, {rel_begin2}, {rel_end2}, {num_pts}, {count_time}, det={det.name})"
    )
    # No shape=: both motors sweep one line, so this must NOT be gridded.
    scan = scan_csv.open_scan(
        csv_file,
        [motor1, motor2],
        det=det,
        scan_type="d2scan",
        command=command,
        num_points=num_pts,
        count_time=count_time,
        h5_file=h5_file,
        comment=comment,
    )
    scan.write_header()
    expt.file_name = folder_prefix

    start1 = motor1.position
    start2 = motor2.position
    positions1 = np.linspace(start1 + rel_begin1, start1 + rel_end1, num_pts)
    positions2 = np.linspace(start2 + rel_begin2, start2 + rel_end2, num_pts)

    # CHANGED (8). Installed before anything moves, and handed back in the
    # finally below so it is given up even if the scan raises.
    _install_guard()
    try:
        if is_tetramm:
            _arm_tetramm(det, save_img)

            def inner_tetramm():
                try:
                    for p1, p2 in zip(positions1, positions2, strict=True):
                        _move_motors([(motor1, p1), (motor2, p2)])
                        det.acquire.put(1, wait=False)
                        time.sleep(count_time)
                        scan.add_point(p1, p2)
                finally:
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_tetramm()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_tetramm(det, save_img)
            return scan

        if is_lambda:
            _arm_lambda(det, num_pts, count_time)

            def step_lambda(p1, p2, frame_num):
                _move_motors([(motor1, p1), (motor2, p2)])
                softglue.start_pulses.put("1!")
                if save_img == 1:
                    # CHANGED (11): his d2scan waited on an undefined `frame_num`
                    # -- a guaranteed NameError on the very first point, so this
                    # branch has never run. Counted the way his dmesh counts.
                    while det.hdf1.num_captured.get() < frame_num:
                        time.sleep(0.005)
                else:
                    time.sleep(count_time)  # CHANGED (4)
                scan.add_point(p1, p2)

            def inner_lambda():
                det.cam.acquire.put(1)
                det.hdf1.capture.put(1)
                shutteron()
                showbeam()
                frame_num = 0
                try:
                    for p1, p2 in zip(positions1, positions2, strict=True):
                        frame_num += 1
                        step_lambda(p1, p2, frame_num)
                finally:
                    _stop_lambda(det)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_lambda()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_lambda(det)
            return scan

        # eiger4M
        if is_eiger:
            _arm_eiger(det, num_pts, count_time)

            fired = 0

            def step(detectors, p1, p2):
                nonlocal fired
                _move_motors([(motor1, p1), (motor2, p2)])
                det.cam.special_trigger_button.put(1, wait=False)
                fired += 1
                time.sleep(count_time)  # CHANGED (10)
                scan.add_point(p1, p2)

            def inner():
                showbeam()
                det.cam.acquire.put(1, wait=False)
                try:
                    for p1, p2 in zip(positions1, positions2, strict=True):
                        step([det], p1, p2)
                finally:
                    _drain_eiger(det, num_pts, count_time, triggered=fired)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_eiger(det)
            return scan

        scan.close("error")
        raise ValueError(f"Unrecognized detector {det.name!r} (expected eiger, lambda or tetramm)")
    finally:
        _restore_guard()


def ascan(motor, abs_begin, abs_end, num_pts, count_time, det=None, att_ratio=7, save_img=1, comment=""):
    """
    Absolute single-motor scan with per-step triggering for eiger4M, lambda2M, and tetramm.

    args:
        motor: ophyd positioner -- huber.delta, sample.x, rheometer.y, ...
        abs_begin, abs_end: absolute start/end (motor units)
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (eiger4M, lambda2M, or tetramm1); default eiger4M
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
        comment: free-text note written near the top of the .csv, so you can
            grep the data folder for a keyword and land on the right scan.
            Newlines are flattened to spaces; empty means no comment line.

    returns:
        the :class:`scan_csv.ScanCsv` that was written -- ``.path`` is the file,
        ``.data`` is ``{column_label: [values]}``.

    On Ctrl+C the motor stops and returns to where it started, acquisition
    halts, the beam is blocked, and the .csv closes with ``#exit_status,aborted``.
    Repeated Ctrl+C during that cleanup is absorbed, not obeyed.
    """
    if det is None:
        det = get_connected_device(DEFAULT_DETECTOR)  # CHANGED (14)

    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name.lower()
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

    # One measurement number per scan whether or not images are saved, so every
    # .csv has a unique name that matches its .h5 when there is one.
    folder_prefix = gen_folder_prefix()
    file_path = data_folder()
    h5_file = f"{file_path}/{folder_prefix}.h5" if save_img == 1 else ""
    csv_file = f"{file_path}/{folder_prefix}.csv"
    if save_img != 1:
        print(f"save_img=0: no .h5 will be written; scan file is {folder_prefix}.csv")

    save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    # This is what run_wrapper(md=...) used to do. Resolved before anything
    # moves, so a typo in the template raises here and not halfway through a
    # scan with the beam on.
    command = f"ascan({motor.name}, {abs_begin}, {abs_end}, {num_pts}, {count_time}, det={det.name})"
    scan = scan_csv.open_scan(
        csv_file,
        motor,
        det=det,
        scan_type="ascan",
        command=command,
        num_points=num_pts,
        count_time=count_time,
        h5_file=h5_file,
        comment=comment,
    )
    scan.write_header()
    # One bare name for both files -- the .h5 and .csv differ only by extension,
    # and the folder is rebuildable from the other settings.
    expt.file_name = folder_prefix

    # CHANGED (8). Installed before anything moves, and handed back in the
    # finally below so it is given up even if the scan raises.
    _install_guard()
    try:
        if is_tetramm:
            _arm_tetramm(det, save_img)

            start_pos = motor.position
            positions = np.linspace(abs_begin, abs_end, num_pts)

            def inner_tetramm():
                try:
                    for pos in positions:
                        motor.move(pos, wait=True)
                        # The TetrAMM free-runs; trigger() only pokes Acquire and
                        # stamps a datum, so this is what trigger_and_read did.
                        det.acquire.put(1, wait=False)
                        time.sleep(count_time)
                        scan.add_point(pos)
                finally:
                    _return_motors([(motor, start_pos)])

            exit_status = "success"
            try:
                inner_tetramm()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_tetramm(det, save_img)
            return scan

        if is_lambda:
            _arm_lambda(det, num_pts, count_time)

            start_pos = motor.position
            positions = np.linspace(abs_begin, abs_end, num_pts)

            def step_lambda(detectors, pos, frame_num):
                motor.move(pos, wait=True)
                softglue.start_pulses.put("1!")
                if save_img == 1:
                    while det.hdf1.num_captured.get() < frame_num:
                        time.sleep(0.005)
                else:
                    # CHANGED (4): with save_img=0 nothing is capturing, so
                    # num_captured never advances and his wait above spins for
                    # ever. Wait out the exposure instead.
                    time.sleep(count_time)
                scan.add_point(pos)  # was: create / read motor+stats1..3 / save

            def inner_lambda():
                det.cam.acquire.put(1)
                det.hdf1.capture.put(1)
                shutteron()
                showbeam()
                try:
                    for ii, pos in enumerate(positions):
                        step_lambda([det], pos, ii + 1)
                finally:
                    _stop_lambda(det)
                    _blockbeam_verified()
                    _return_motors([(motor, start_pos)])

            exit_status = "success"
            try:
                inner_lambda()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_lambda(det)
            return scan

        # eiger4M
        if is_eiger:
            _arm_eiger(det, num_pts, count_time)

            start_pos = motor.position
            positions = np.linspace(abs_begin, abs_end, num_pts)

            fired = 0

            def step(detectors, pos):
                """Move motor, fire software trigger, wait for frame, read motor."""
                nonlocal fired
                motor.move(pos, wait=True)
                det.cam.special_trigger_button.put(1, wait=False)
                fired += 1
                time.sleep(count_time)
                scan.add_point(pos)  # was: create / read motor+stats1..3 / save

            def inner():
                """Start pre-armed acquisition (accepts num_pts software triggers)"""
                det.cam.acquire.put(1, wait=False)
                showbeam()
                try:
                    for pos in positions:
                        step([det], pos)
                finally:
                    _drain_eiger(det, num_pts, count_time, triggered=fired)
                    # His dscan blocked the beam in the OUTER finally, i.e. after
                    # the return move; his other five block here, before it.
                    # Theirs is used -- the return move is the long step and
                    # there is no reason to spend it with the shutter open.
                    _blockbeam_verified()
                    _return_motors([(motor, start_pos)])

            exit_status = "success"
            try:
                inner()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                # Return detector to normal state: Internal Enable, manual trigger off
                _disarm_eiger(det)
            return scan

        scan.close("error")
        raise ValueError(f"Unrecognized detector {det.name!r} (expected eiger, lambda or tetramm)")
    finally:
        _restore_guard()


def a2scan(
    motor1,
    abs_begin1,
    abs_end1,
    motor2,
    abs_begin2,
    abs_end2,
    num_pts,
    count_time,
    det=None,
    att_ratio=7,
    save_img=1,
    comment="",
):
    """
    Two-motor absolute scan with per-step triggering for eiger4M, lambda2M, and tetramm.

    Both motors sweep together along a single line: point i puts motor1 at
    positions1[i] AND motor2 at positions2[i]. This is a 1D trajectory, not a
    raster -- use dmesh/mesh for a grid, and plot the result in the 1D viewer
    against whichever motor you like.

    usage:
        a2scan(m1, s1, f1, m2, s2, f2, npts, ct)

    args:
        motor1, motor2: ophyd positioners (e.g. huber.eta, huber.chi)
        abs_begin1, abs_end1: absolute start/end for motor1
        abs_begin2, abs_end2: absolute start/end for motor2
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (eiger4M, lambda2M, or tetramm1); default eiger4M
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
        comment: free-text note written near the top of the .csv

    returns:
        the :class:`scan_csv.ScanCsv` that was written.
    """
    if det is None:
        det = get_connected_device(DEFAULT_DETECTOR)  # CHANGED (14)

    pre_align()
    att(att_ratio)
    PIND_status(0)

    is_tetramm = "tetramm" in det.name.lower()
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    is_lambda = ("lambda" in det.name.lower()) or ("lambda" in det.prefix.lower())

    folder_prefix = gen_folder_prefix()
    file_path = data_folder()
    h5_file = f"{file_path}/{folder_prefix}.h5" if save_img == 1 else ""
    csv_file = f"{file_path}/{folder_prefix}.csv"
    if save_img != 1:
        print(f"save_img=0: no .h5 will be written; scan file is {folder_prefix}.csv")

    save_images(det, save_img, num_pts, folder_prefix=folder_prefix)

    command = (
        f"a2scan({motor1.name}, {abs_begin1}, {abs_end1}, "
        f"{motor2.name}, {abs_begin2}, {abs_end2}, {num_pts}, {count_time}, det={det.name})"
    )
    # No shape=: both motors sweep one line, so this must NOT be gridded.
    scan = scan_csv.open_scan(
        csv_file,
        [motor1, motor2],
        det=det,
        scan_type="a2scan",
        command=command,
        num_points=num_pts,
        count_time=count_time,
        h5_file=h5_file,
        comment=comment,
    )
    scan.write_header()
    expt.file_name = folder_prefix

    start1 = motor1.position
    start2 = motor2.position
    positions1 = np.linspace(abs_begin1, abs_end1, num_pts)
    positions2 = np.linspace(abs_begin2, abs_end2, num_pts)

    # CHANGED (8). Installed before anything moves, and handed back in the
    # finally below so it is given up even if the scan raises.
    _install_guard()
    try:
        if is_tetramm:
            _arm_tetramm(det, save_img)

            def inner_tetramm():
                try:
                    for p1, p2 in zip(positions1, positions2, strict=True):
                        _move_motors([(motor1, p1), (motor2, p2)])
                        det.acquire.put(1, wait=False)
                        time.sleep(count_time)
                        scan.add_point(p1, p2)
                finally:
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_tetramm()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_tetramm(det, save_img)
            return scan

        if is_lambda:
            _arm_lambda(det, num_pts, count_time)

            def step_lambda(p1, p2, frame_num):
                _move_motors([(motor1, p1), (motor2, p2)])
                softglue.start_pulses.put("1!")
                if save_img == 1:
                    # CHANGED (11): his d2scan waited on an undefined `frame_num`
                    # -- a guaranteed NameError on the very first point, so this
                    # branch has never run. Counted the way his dmesh counts.
                    while det.hdf1.num_captured.get() < frame_num:
                        time.sleep(0.005)
                else:
                    time.sleep(count_time)  # CHANGED (4)
                scan.add_point(p1, p2)

            def inner_lambda():
                det.cam.acquire.put(1)
                det.hdf1.capture.put(1)
                shutteron()
                showbeam()
                frame_num = 0
                try:
                    for p1, p2 in zip(positions1, positions2, strict=True):
                        frame_num += 1
                        step_lambda(p1, p2, frame_num)
                finally:
                    _stop_lambda(det)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner_lambda()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_lambda(det)
            return scan

        # eiger4M
        if is_eiger:
            _arm_eiger(det, num_pts, count_time)

            fired = 0

            def step(detectors, p1, p2):
                nonlocal fired
                _move_motors([(motor1, p1), (motor2, p2)])
                det.cam.special_trigger_button.put(1, wait=False)
                fired += 1
                time.sleep(count_time)  # CHANGED (10)
                scan.add_point(p1, p2)

            def inner():
                showbeam()
                det.cam.acquire.put(1, wait=False)
                try:
                    for p1, p2 in zip(positions1, positions2, strict=True):
                        step([det], p1, p2)
                finally:
                    _drain_eiger(det, num_pts, count_time, triggered=fired)
                    _blockbeam_verified()
                    _return_motors([(motor1, start1), (motor2, start2)])

            exit_status = "success"
            try:
                inner()
            except KeyboardInterrupt:
                exit_status = "aborted"
                print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
            finally:
                scan.close(exit_status)
                _disarm_eiger(det)
            return scan

        scan.close("error")
        raise ValueError(f"Unrecognized detector {det.name!r} (expected eiger, lambda or tetramm)")
    finally:
        _restore_guard()


def auto_att(
    det,
    pilot_exptime: float = 0.05,
    rate_limit: float = 1e5,
    filter_factor: float = 5.0,
    retry_max: int = 10,
    grace_factor: float = 0.25,
):
    """Find the optimal attenuation using short pilot exposures.

    Args:
        det:            eiger4M or lambda2M
        pilot_exptime:  duration of each test frame (s)
        rate_limit:     max acceptable count rate (max pixel cts/s)
        filter_factor:  transmission multiplier per step, must be > 1
        retry_max:      max iterations before giving up
        grace_factor:   lower rate bound = rate_limit * grace_factor

    Example:

        auto_att(eiger4M, pilot_exptime=0.05, rate_limit=4e5)

    Copied from scan_8id.py:1671 essentially unchanged -- it contains no
    `yield`, so it was already Ophyd-only. The one change is that filter_8ide is
    resolved here instead of at import time; see CHANGED (14).
    """
    filter_beam = get_connected_device("filter_8ide")

    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    low_rate = rate_limit * grace_factor

    orig_acq_time = det.cam.acquire_time.get()
    orig_acq_period = det.cam.acquire_period.get()

    det.cam.acquire_time.put(pilot_exptime)
    det.cam.acquire_period.put(pilot_exptime)

    if is_eiger:
        det.cam.trigger_mode.put("Internal Series")
        det.cam.num_images.put(1)
        det.cam.num_triggers.put(1)
        det.cam.manual_trigger.put("Disable")
    else:
        # lambda2M
        det.cam.trigger_mode.put("Internal")
        det.cam.num_images.put(1)

    det.stats1.enable.put(1)
    det.stats1.compute_statistics.put(1)

    # Start from maximum attenuation (minimum transmission) for safety
    filter_beam.transmission.move(1e-10)
    time.sleep(0.5)

    # The beam is opened here, so the same interrupt guard applies: a ^C during
    # a pilot exposure must not leave the shutter open.
    # CHANGED (8). Installed before anything moves, and handed back in the
    # finally below so it is given up even if the scan raises.
    _install_guard()
    try:
        showbeam()
        try:
            for attempt in range(retry_max):
                det.cam.acquire.put(1)
                t0 = time.time()
                timeout = pilot_exptime * 5 + 2
                while det.cam.acquire.get() == 1:
                    time.sleep(0.02)
                    if time.time() - t0 > timeout:
                        print("  WARNING: pilot frame timed out")
                        break

                max_cts = det.stats1.max_value.get()
                rate = max_cts / pilot_exptime
                current_trans = filter_beam.transmission.readback.get()

                print(f"Attempt {attempt + 1}: trans={current_trans:.4f}  max_cts={max_cts:.0f}  rate={rate:.0f} cts/s")

                if rate > rate_limit:
                    new_trans = current_trans / filter_factor
                    print(f"Rate too high. Reducing transmission to {new_trans:.6f}")
                    filter_beam.transmission.move(new_trans)

                elif rate < low_rate:
                    if rate > 0:
                        new_trans = current_trans * (0.75 * rate_limit / rate)
                    else:
                        # rate=0: jump to a coarse fraction of max rather than
                        # tiny multiplier steps
                        new_trans = min(current_trans * 1000, 0.01)
                    new_trans = min(new_trans, 1.0)
                    if new_trans >= current_trans * 0.999:
                        if current_trans >= 0.999:
                            print("WARNING: already at max transmission, beam too weak.")
                            break
                    print(f"Rate too low. Raising transmission to {new_trans:.4f}")
                    filter_beam.transmission.move(new_trans)
                else:
                    print(f"    Rate in [{low_rate:.0f}, {rate_limit:.0f}] cts/s -- converged.")
                    break
            else:
                print(f"WARNING: auto_att did not converge in {retry_max} attempts")
        except KeyboardInterrupt:
            print("\n^C  auto_att aborted.", flush=True)
        finally:
            _blockbeam_verified()  # CHANGED (7)
            _enter_cleanup()
            _safe("restore acquire_time", det.cam.acquire_time.put, orig_acq_time)
            _safe("restore acquire_period", det.cam.acquire_period.put, orig_acq_period)
    finally:
        _restore_guard()

    trans = filter_beam.transmission.readback.get()
    atten = filter_beam.attenuation.readback.get()
    print(f"  Final: transmission={trans:.4f}  attenuation={atten}")


# ===========================================================================
# CHANGED (16): the suffixed aliases exist ONLY so that a Bluesky session, which
# does `from .plans.align.scan_8id import *` (startup.py:218), cannot end up
# with one dscan shadowing the other.
#
# The failure that prevents is not a clean error. Python evaluates arguments
# before it calls anything, so if the Ophyd dscan won that name,
#
#     RE(dscan(huber.eta, -0.5, 0.5, 40, 0.1, lambda2M, att_ratio=10))
#
# would run the ENTIRE Ophyd scan -- move the motor, open the shutter, expose,
# write the .h5, burn a measurement number -- and only THEN raise, because RE()
# got a ScanCsv instead of a generator. The user reads "error", runs it again,
# and double-exposes the sample. That exact call form is on record in
# src/user_plans/Archive/zhou202607.notes:14.
#
# __all__ restricts `import *` only, so the plain names below stay available to
# an explicit `from ... import dscan` -- which is what startup_ophyd.py does,
# and it never imports scan_8id at all.
# ===========================================================================
dmesh_ophyd, mesh_ophyd, dscan_ophyd = dmesh, mesh, dscan
d2scan_ophyd, ascan_ophyd, a2scan_ophyd = d2scan, ascan, a2scan
auto_att_ophyd = auto_att

__all__ = [
    "dmesh_ophyd",
    "mesh_ophyd",
    "dscan_ophyd",
    "d2scan_ophyd",
    "ascan_ophyd",
    "a2scan_ophyd",
    "auto_att_ophyd",
]
