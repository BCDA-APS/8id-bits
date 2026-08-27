"""
Ophyd-only scanning with live CSV output for the 8-ID beamlines.

This is ``save_images()`` and ``dscan()`` from ``scan_8id.py`` (Sam's code) with
Bluesky taken out. It is kept deliberately close to the original -- same
function shape, same three detector branches, same nested ``step_*``/``inner_*``
helpers, same ``try: ... finally:`` bodies, same variable names -- so anyone who
knows ``scan_8id.py`` can read, edit and maintain this file.

The translation rule was mechanical:

===============================================  ==============================
scan_8id.py                                      here
===============================================  ==============================
``yield from bps.mv(sig, x)``                    ``sig.put(x)``
``yield from bps.mv(sig1, a, sig2, b)``          one ``.put()`` per line
``yield from bps.move_per_step({m: p}, cache)``  ``m.move(p, wait=True)``
``yield from bps.sleep(t)``                      ``time.sleep(t)``
``create / read / read / ... / save``            ``scan.add_point(pos)``
``bpp.run_wrapper(..., md=md)``                  ``scan_csv.open_scan(...)``
``bpp.stage_wrapper(...)``                       gone -- nothing stages
===============================================  ==============================

Anything else that differs from ``scan_8id.py`` is a bug that was fatal here,
and every one of them is marked with a ``# CHANGED:`` comment saying why. There
are five. Nothing else was "improved", tidied or restructured.

Each scan writes ONE .csv, named after that scan's .h5::

    .../data/bluesky/A0201_Test_a0010.h5     images
    .../data/bluesky/A0201_Test_a0010.csv    motor positions and counters

The .csv is closed after every point, so a viewer can poll it while the scan
runs. Both names go to ``pv_registers.scan_h5_file`` / ``scan_csv_file`` so a
GUI can find the live scan.

What goes into the .csv is not hard-coded here -- the header lines and the data
columns come from ``configs/scan_csv_template.yml``; see :mod:`scan_csv`.

usage (call it directly -- do NOT wrap it in ``RE()``)::

    dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
    dscan_ophyd(sample.x, -0.1, 0.1, 21, 0.5, det=lambda2M, comment="spot 5")

Any positioner works -- ``huber.delta``, ``sample.x``, ``rheometer.y`` -- as
long as it has ``.name`` / ``.position`` / ``.move()``.

Like ``scan_8id.py``, this module resolves devices from the ``oregistry`` at
import time, so it must be imported after ``make_devices()`` in ``startup.py``.
"""

import time

import numpy as np
from apsbits.core.instrument_init import oregistry
from id8_common.plans.acquire.ad_acq import gen_folder_prefix
from id8_common.plans.align import scan_csv
from id8_common.plans.set.shutter_att import PIND_status
from id8_common.plans.set.shutter_att import att
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import pre_align
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.plans.set.shutter_att import shutteron

pv_registers = oregistry["pv_registers"]
softglue = oregistry["softglue"]
softglue_8id_acq = oregistry["softglue_8id_acq"]
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]


def data_folder():
    """``/gdata/dm/8ID/8IDE/<cycle>/<experiment>/data/bluesky``.

    scan_8id.py builds this same string with the run cycle written out in full
    ("2026-2"); here it comes from ``pv_registers.cycle_name`` so the .csv and
    the .h5 cannot end up in different folders after a cycle rolls over. Today
    the two produce exactly the same path.
    """
    mount = pv_registers.mount_point.get().strip()
    cycle = pv_registers.cycle_name.get().strip()
    experiment = pv_registers.experiment_name.get().strip()
    return f"{mount}{cycle}/{experiment}/data/bluesky"


def save_images(det, save_img, num_pts, num_frames=1, file_path=None, folder_prefix=None):
    """
    Toggle saving (1) or not saving (0) images.

    det: detector instance (eiger4M, lambda2M, or tetramm)
    save_img: 1 save, 0 don't save
    num_frames: number of frames to capture per point (default 1)
    file_path: override base path (default uses a safe local path)
    """
    if save_img not in (0, 1):
        raise ValueError("save_img must be 1 or 0 (to save or not to save)")

    if save_img == 1:
        if folder_prefix is None:
            folder_prefix = gen_folder_prefix()

        if file_path is None:
            file_path = data_folder()

        is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
        is_tetramm = "tetramm" in det.name.lower()

        def has(obj, attr):
            return getattr(obj, attr, None) is not None

        # CHANGED (1): scan_8id.py calls gen_folder_prefix() again on this line,
        # overwriting whatever the caller passed in. That is fatal here: the
        # caller has already named the .csv after its prefix, so the .h5 would
        # get a different name and the pair would no longer match. It also
        # burned two measurement numbers per scan.
        file_name = folder_prefix
        file_header = folder_prefix

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

        if has(det, "cam"):
            if is_eiger:
                if has(det.cam, "trigger_mode"):
                    det.cam.trigger_mode.put("Internal Enable")
                if has(det.cam, "num_triggers"):
                    det.cam.num_triggers.put(1)
                if has(det.cam, "num_images"):
                    det.cam.num_images.put(1)

                # CHANGED (2): bpp.stage_wrapper used to enable the HDF plugin
                # and put it in Stream mode (apstools AD_EpicsHdf5FileName).
                # Nothing stages any more, so without these three the scan runs
                # happily and writes no file at all.
                det.hdf1.enable.put(1)
                det.hdf1.auto_save.put(1)
                det.hdf1.file_write_mode.put(2)  # Stream

                # CHANGED (3): scan_8id.py calls setup_eiger_internal(...) here.
                # It re-points hdf1.file_path/file_name at get_common_file_path(),
                # a different folder, which separates the .h5 from its .csv and
                # breaks the viewer's "find the live scan" lookup. Everything
                # else it sets (acquire_time/period, trigger_mode, num_images,
                # num_triggers) dscan_ophyd overwrites moments later anyway.
                #     setup_eiger_internal(
                #         acq_time=float(det.cam.acquire_time.get()) if has(det.cam, "acquire_time") else 1.0,
                #         num_frames=num_frames,
                #         file_header=file_header,
                #         file_name=file_name,
                #     )
            else:
                if has(det.cam, "num_images"):
                    det.cam.num_images.put(num_frames)


def dscan_ophyd(motor, rel_begin, rel_end, num_pts, count_time, det=eiger4M, att_ratio=1e6, save_img=1, comment=""):
    """
    Pre-armed software-trigger scan for fast acquisitions.

    args:
        motor: ophyd positioner -- huber.delta, sample.x, rheometer.y, ...
        rel_begin, rel_end: relative start/end (motor units)
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: detector (eiger4M, lambda2M, or tetramm1)
        att_ratio: attenuation ratio
        save_img: 1 save, 0 don't save
        comment: free-text note written near the top of the .csv, so you can
            grep the data folder for a keyword and land on the right scan.
            Newlines are flattened to spaces; empty means no comment line.

    returns:
        the :class:`scan_csv.ScanCsv` that was written -- ``.path`` is the file,
        ``.data`` is ``{column_label: [values]}``.

    On Ctrl+C the motor stops and returns to where it started, acquisition
    halts, the beam is blocked, and the .csv closes with ``#END,aborted``.
    """
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
    command = f"dscan_ophyd({motor.name}, {rel_begin}, {rel_end}, {num_pts}, {count_time}, det={det.name})"
    scan = scan_csv.open_scan(
        csv_file,
        motor,
        det=det,
        scan_type="dscan_ophyd",
        command=command,
        num_points=num_pts,
        count_time=count_time,
        h5_file=h5_file,
        comment=comment,
    )
    scan.write_header()
    pv_registers.scan_h5_file.put(h5_file)
    pv_registers.scan_csv_file.put(csv_file)

    if is_tetramm:
        det.hdf1.enable.put(1)
        det.hdf1.capture.put(1)
        print(f"TetrAMM HDF capture armed: {det.hdf1.file_name.get()}")

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
                motor.stop()  # CHANGED (5) -- see the note in inner_lambda
                time.sleep(0.2)
                motor.move(start_pos, wait=True)

        exit_status = "success"
        try:
            inner_tetramm()
        except KeyboardInterrupt:
            exit_status = "aborted"
            print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
        finally:
            scan.close(exit_status)
            if save_img == 1:
                det.hdf1.capture.put(0)
        return scan

    if is_lambda:
        det.cam.operating_mode.put(3)  # 24-bit dual threshold mode
        det.cam.trigger_mode.put("External_ImagePer")
        det.cam.acquire_time.put(count_time)
        det.cam.acquire_period.put(count_time)
        det.cam.num_images.put(num_pts)
        det.hdf1.num_capture.put(num_pts)
        softglue.num_triggers.put(1)  # one pulse per "1!"; else stale value from run_measurement
        softglue.acq_time.put(count_time)
        softglue.acq_period.put(count_time)

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
                # num_captured never advances and scan_8id.py's wait above spins
                # for ever. Wait out the exposure instead.
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
                softglue.stop_pulses.put("1!")
                det.cam.acquire.put(0)
                det.hdf1.capture.put(0)
                blockbeam()
                # CHANGED (5): stop() before the return move. A ^C leaves the
                # interrupted move's MoveStatus unfinished, and ophyd then
                # refuses the next move() ("another set() is still in
                # progress") -- so scan_8id.py's line below silently does
                # nothing and the motor is stranded mid-scan.
                motor.stop()
                time.sleep(0.2)
                motor.move(start_pos, wait=True)

        exit_status = "success"
        try:
            inner_lambda()
        except KeyboardInterrupt:
            exit_status = "aborted"
            print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
        finally:
            scan.close(exit_status)
            det.cam.operating_mode.put(3)
            det.cam.trigger_mode.put(0)
            softglue_8id_acq.preset.put(50)
            blockbeam()
            shutteroff()
            print("# images captured: ", det.hdf1.num_captured.get())
        return scan

    # eiger4M
    if is_eiger:
        det.cam.acquire_time.put(count_time)
        det.cam.acquire_period.put(count_time)
        det.cam.trigger_mode.put("Internal Series")
        det.cam.manual_trigger.put("Enable")
        det.cam.num_images.put(1)
        det.cam.num_triggers.put(num_pts)
        det.hdf1.num_capture.put(num_pts)

        # scan_8id.py saves and patches det.cam.stage_sigs here, then restores
        # it in the outer finally. stage_sigs are only applied by stage(), and
        # nothing stages without bpp.stage_wrapper, so both halves would be
        # dead code. The .put() calls just above are what arms the detector.

        start_pos = motor.position
        positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)

        def step(detectors, pos):
            """Move motor, fire software trigger, wait for frame, read motor."""
            motor.move(pos, wait=True)
            det.cam.special_trigger_button.put(1, wait=False)
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
                t0 = time.time()
                timeout = num_pts * count_time + 10
                while det.hdf1.num_captured.get() < num_pts:
                    time.sleep(0.05)
                    if time.time() - t0 > timeout:
                        print("WARNING: HDF write timeout -- not all frames saved.")
                        break
                det.cam.acquire.put(0)
                motor.stop()  # CHANGED (5) -- see the note in inner_lambda
                time.sleep(0.2)
                motor.move(start_pos, wait=True)  # return to start position

        exit_status = "success"
        try:
            inner()
        except KeyboardInterrupt:
            exit_status = "aborted"
            print(f"\n^C  Scan aborted at point {scan.num_points}/{num_pts}.", flush=True)
        finally:
            scan.close(exit_status)
            # Return detector to normal state: Internal Enable, manual trigger off
            det.cam.trigger_mode.put("Internal Enable")
            det.cam.manual_trigger.put("Disable")
            blockbeam()
            print("# images captured: ", det.hdf1.num_captured.get())
        return scan

    scan.close("error")
    raise ValueError(f"Unrecognized detector {det.name!r} (expected eiger, lambda or tetramm)")
