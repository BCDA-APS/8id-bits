"""
Ophyd-only scanning with live SPEC-format output for the 8-ID beamlines.

**No Bluesky.** This module contains no RunEngine, no plans, no generators, no
``plan_stubs``, and emits no documents. Everything is plain Python driving Ophyd
signals with ``.put()`` / ``.get()`` / ``.move()``.

The scan writes one row to a SPEC-format ``.spec`` file per point, closing the
file after every row, so an external plotting tool can poll the file while the
scan runs. Because the file is append-only and never held open, a reader can
never collide with the writer -- which is the whole reason for doing it this way
rather than reading the detector ``.h5`` mid-scan.

**What goes into the file is not hard-coded here.** The header lines, the
``#O``/``#P`` positioner snapshot and the data columns all come from
``configs/spec_template.yml``; edit that file to add a column or drop a header
entry, no Python change needed. Resolution lives in :mod:`ophyd_spec_config`, and
the writing itself in :mod:`ophyd_spec_writer` (which imports no hardware at all).

usage (call it directly -- do NOT wrap it in ``RE()``)::

    from id8_common.plans.align.ophyd_scan import dscan_ophyd
    dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)

Any positioner works -- ``huber.delta``, ``sample.x``, ``rheometer.y`` -- as long
as it has ``.name`` / ``.position`` / ``.move()``.

The companion live viewer is ``specr_py`` (``~/Documents/specr_py``).

Import order note: like ``scan_8id.py``, this module resolves devices from the
``oregistry`` at import time, so it must be imported after ``make_devices()``
has run in ``startup.py``.
"""

import os
import time

import numpy as np
from apsbits.core.instrument_init import oregistry

# Everything imported below is a plain function (not a Bluesky plan) -- verified:
# neither shutter_att.py nor ad_acq.py imports bluesky, so this chain stays clean.
from id8_common.plans.acquire.ad_acq import gen_folder_prefix
from id8_common.plans.align import ophyd_spec_config
from id8_common.plans.align.ophyd_scan_utils import ScanResult
from id8_common.plans.align.ophyd_scan_utils import disable_ctrl_c
from id8_common.plans.align.ophyd_scan_utils import restore_ctrl_c
from id8_common.plans.align.ophyd_scan_utils import safe_call
from id8_common.plans.align.ophyd_scan_utils import table_header
from id8_common.plans.align.ophyd_scan_utils import table_row
from id8_common.plans.align.ophyd_spec_writer import SpecFile
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


# =============================================================================
# Where the files go
# =============================================================================
#
# Every path is built from pv_registers, so the files follow the run cycle
# without anything being hard-coded.


def image_file_path():
    """Directory the detector HDF5 files go to.

    Matches ``save_images()`` in scan_8id.py, which builds
    ``/gdata/dm/8ID/8IDE/<cycle>/<experiment>/data/bluesky`` -- reconstructed
    here from the registers rather than hard-coding the cycle.
    """
    mount = pv_registers.mount_point.get().strip()
    cycle = pv_registers.cycle_name.get().strip()
    experiment = pv_registers.experiment_name.get().strip()
    return f"{mount}{cycle}/{experiment}/data/bluesky"


def default_spec_path():
    """``<mount><cycle>/<experiment>/data/bluesky/<spec_file_name>.spec``.

    The SPEC file sits in the same directory as the detector ``.h5`` files, and
    its name comes from ``pv_registers.spec_file_name``. Writing that PV is how
    you switch SPEC files on demand -- the next scan appends to the named file,
    or creates it (``#F`` header block and all) if it does not exist yet.

    ``.spec`` is appended when the register does not already end in it, because
    the viewer's file dialog filters on that extension.

    Note this directory is on GPFS, which is not mounted on every analysis
    workstation. Run the live viewer somewhere that can see ``/gdata``, or pass
    ``spec_path=`` to put a copy somewhere shared.
    """
    name = pv_registers.spec_file_name.get().strip()
    if not name:
        # A blank register would otherwise give a file called ".spec".
        name = pv_registers.experiment_name.get().strip()
        print(f"WARNING: pv_registers.spec_file_name is empty; using {name!r}.")
    if not name.lower().endswith(".spec"):
        name = name + ".spec"
    return os.path.join(image_file_path(), name)


def scan_number_from_prefix(folder_prefix, fallback=0):
    """``"A0092_G10_a0007"`` -> ``92``.

    Reuses the number ``gen_folder_prefix()`` already burned, so SPEC ``#S 92``
    lines up with ``A0092_*.h5``. Do not read ``measurement_num`` again here --
    ``gen_folder_prefix()`` increments it as a side effect.
    """
    try:
        digits = "".join(c for c in folder_prefix.split("_")[0] if c.isdigit())
        return int(digits)
    except (ValueError, IndexError):
        return fallback


# =============================================================================
# Detector plumbing
# =============================================================================


def detector_kind(det):
    """Return 'eiger', 'lambda', or 'tetramm'."""
    name = (det.name or "").lower()
    prefix = (getattr(det, "prefix", "") or "").lower()
    if "eiger" in name or "eiger" in prefix:
        return "eiger"
    if "lambda" in name or "lambda" in prefix:
        return "lambda"
    if "tetramm" in name:
        return "tetramm"
    raise ValueError(f"Unrecognized detector {det.name!r}")


def arm_hdf(det, file_path, file_name, num_capture, kind):
    """
    Prepare the HDF5 plugin so images are written exactly as ``dscan()`` does.

    The two detectors need different handling because ``dscan()`` in scan_8id.py
    treats them differently:

    * **lambda2M** -- ``dscan`` stages only the motor, so the HDF plugin is
      never staged. Everything it gets comes from ``save_images()``:
      ``num_capture``, ``file_name``, ``file_path``. Capture is switched on
      later, after ``cam.acquire``. Nothing else is touched, so the IOC's own
      ``FileWriteMode``/``AutoSave`` settings are left exactly as the beamline
      configured them. Replicated faithfully here.
    * **eiger4M** -- ``dscan`` stages the detector, and Bluesky staging applies
      ``file_write_mode="Stream"`` and ``capture=1`` (verified in
      ``apstools.devices.area_detector_support.AD_EpicsHdf5FileName``). Without
      Bluesky nothing applies those, so they are set explicitly -- omit them and
      the scan runs but writes no file.
    """
    det.hdf1.num_capture.put(num_capture)
    det.hdf1.file_name.put(file_name)
    det.hdf1.file_path.put(file_path)

    if kind == "eiger":
        # Replace what stage_wrapper would have done.
        det.hdf1.enable.put(1)
        det.hdf1.file_write_mode.put(2)  # Stream
        det.hdf1.auto_save.put(1)
        det.hdf1.capture.put(1)  # must be last


def report_frames_captured(det):
    """Print how many frames the HDF plugin actually wrote."""
    print(f"# images captured: {det.hdf1.num_captured.get()}")


def disarm_hdf(det):
    """Stop capture so the IOC closes the file."""
    try:
        det.hdf1.capture.put(0)
    except Exception as exc:
        print(f"WARNING: could not stop HDF capture: {exc}")


def wait_for_frames(det, n_expected, timeout):
    """Block until the HDF plugin reports n_expected frames, or timeout."""
    t0 = time.time()
    while det.hdf1.num_captured.get() < n_expected:
        time.sleep(0.005)
        if time.time() - t0 > timeout:
            print(
                f"WARNING: timeout waiting for frame {n_expected} "
                f"(captured {det.hdf1.num_captured.get()})."
            )
            return False
    return True


# =============================================================================
# The scan
# =============================================================================


def dscan_ophyd(
    motor,
    rel_begin,
    rel_end,
    num_pts,
    count_time,
    det=None,
    comment="",
    att_ratio=1e6,
    save_img=1,
    spec_path=None,
    spec_template=None,
    set_attenuation=True,
    beam_control=True,
    verbose=True,
):
    """
    Relative single-motor scan, Ophyd only, with live SPEC output.

    Mirrors the hardware choreography of ``dscan`` in scan_8id.py (pre-armed
    software trigger for the Eiger, softglue-pulsed external trigger for the
    Lambda) but contains no Bluesky.

    Call directly::

        dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
        dscan_ophyd(sample.x, -0.1, 0.1, 21, 0.5, det=lambda2M)

    args:
        motor: any ophyd positioner -- huber.delta, sample.x, rheometer.y, ...
        rel_begin, rel_end: start/end relative to the current position
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: eiger4M (default) or lambda2M
        comment: free-text note describing this scan -- letters, digits, spaces
            and punctuation are all fine. Written as the **first** ``#MD`` line
            of the scan block so you can Ctrl+F the ``.spec`` file for a keyword
            and land on the scan you want::

                dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M,
                            comment="3x3 grid spot 5, after realigning KB")

            Any newlines are flattened to spaces so the file structure cannot
            break. Leave it empty and no ``#MD comment`` line is written.
        att_ratio: attenuation ratio, applied only if ``set_attenuation``
        save_img: 1 to write detector images, 0 to scan without saving
        spec_path: override the SPEC file location
        spec_template: override the SPEC template (default: configs/spec_template.yml,
            or ``$ID8_SPEC_TEMPLATE``). Controls every header line and column.
        set_attenuation: if False, leave the attenuator exactly as it is. Use
            when the current attenuation is already correct and you do not want
            the scan changing the beam condition.
        beam_control: if False, skip ``pre_align``/``PIND_status`` and the
            shutter open/close. Use for a dry run that touches only the scanned
            motor and the detector.
        verbose: print a table row per point as the scan runs -- motor readback,
            elapsed time, and every counter -- instead of only seeing the values
            at the end. Set False for a quiet scan.

    returns:
        ``ScanResult``, which unpacks as ``(positions, {label: [values]})`` but
        prints as a one-line summary rather than dumping every value.

    On **Ctrl+C** the scan stops the motor, halts acquisition, closes the image
    file, blocks the beam, drives the motor back to where it started, and closes
    the SPEC block with ``exit_status = aborted`` so the partial scan still
    parses. That teardown is protected from further interrupts, so a second
    Ctrl+C will not strand the motor. The ``KeyboardInterrupt`` is re-raised
    afterwards, so an enclosing loop still stops.
    """
    if det is None:
        det = eiger4M
    kind = detector_kind(det)
    if kind == "tetramm":
        raise NotImplementedError(
            "tetramm is not supported here. Its trigger path and HDF stream mode "
            "differ from the area detectors and were not verified; use dscan() in "
            "scan_8id.py for tetramm scans."
        )

    if beam_control:
        pre_align()
        PIND_status(0)
    if set_attenuation:
        att(att_ratio)

    folder_prefix = gen_folder_prefix() if save_img == 1 else ""
    scan_num = scan_number_from_prefix(folder_prefix)

    # --- SPEC file -------------------------------------------------------
    # Resolved before anything moves: a broken template raises here rather than
    # partway through a scan.
    rendered = ophyd_spec_config.render(
        motor,
        det=det,
        scan_type="dscan_ophyd",
        num_points=num_pts,
        count_time=count_time,
        image_file=f"{folder_prefix}.h5" if folder_prefix else "",
        comment=comment,
        template=spec_template,
    )
    spec = SpecFile(spec_path or default_spec_path())
    spec.write_file_header(rendered.positioner_names)

    command = (
        f"dscan_ophyd({motor.name}, {rel_begin}, {rel_end}, {num_pts}, "
        f"{count_time}, det={det.name})"
    )
    spec.start_scan(
        scan_num,
        command,
        rendered.labels,
        metadata=rendered.metadata,
        motor_positions=rendered.positioner_positions,
        comments=rendered.comments,
    )
    print(f"SPEC file: {spec.path}   (#S {scan_num})")

    # --- detector / file setup -------------------------------------------
    if save_img == 1:
        file_path = image_file_path()
        print(f"Scan folder created: {folder_prefix}")
        print(f"File path: {file_path}")
        arm_hdf(det, file_path, folder_prefix, num_pts, kind)

    start_pos = motor.position
    positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
    t_start = time.time()
    counter_labels = rendered.counter_labels
    measured, columns = [], {lab: [] for lab in counter_labels}
    status = "success"
    aborted = False

    if verbose:
        print(table_header(motor.name, counter_labels), flush=True)

    def record(setpoint):
        """Read every template column for one point and append a SPEC row.

        ``setpoint`` is where the motor was told to go; the motor column is
        where it actually is -- logging only the commanded value would draw a
        perfect ramp even if the motor never moved.
        """
        elapsed = time.time() - t_start
        actual = motor.position
        values = rendered.read(actual, setpoint, elapsed)
        spec.add_point(values)
        measured.append(actual)
        counters = [values[i] for i in rendered.counter_indices]
        for lab, value in zip(counter_labels, counters, strict=False):
            columns[lab].append(value)
        if verbose:
            print(table_row(len(measured), motor.name, actual, elapsed,
                             counter_labels, counters), flush=True)

    try:
        if kind == "eiger":
            det.cam.acquire_time.put(count_time)
            det.cam.acquire_period.put(count_time)
            det.cam.trigger_mode.put("Internal Series")
            det.cam.manual_trigger.put("Enable")
            det.cam.num_images.put(1)
            det.cam.num_triggers.put(num_pts)
            det.hdf1.num_capture.put(num_pts)

            det.cam.acquire.put(1, wait=False)  # pre-arm for num_pts triggers
            if beam_control:
                showbeam()
            for pos in positions:
                motor.move(pos, wait=True)
                det.cam.special_trigger_button.put(1, wait=False)
                time.sleep(count_time)
                record(pos)

        else:  # lambda
            det.cam.operating_mode.put(3)  # 24-bit dual threshold
            det.cam.trigger_mode.put("External_ImagePer")
            det.cam.acquire_time.put(count_time)
            det.cam.acquire_period.put(count_time)
            det.cam.num_images.put(num_pts)
            det.hdf1.num_capture.put(num_pts)
            softglue.num_triggers.put(1)  # one pulse per "1!"
            softglue.acq_time.put(count_time)
            softglue.acq_period.put(count_time)

            det.cam.acquire.put(1)
            if save_img == 1:
                det.hdf1.capture.put(1)  # after acquire, exactly as dscan() does
            if beam_control:
                shutteron()
                showbeam()
            for i, pos in enumerate(positions):
                motor.move(pos, wait=True)
                softglue.start_pulses.put("1!")
                if save_img == 1:
                    # scan_8id.py spins here with no timeout; bounded here so a
                    # missed trigger cannot hang the session indefinitely.
                    wait_for_frames(det, i + 1, count_time * 5 + 10)
                else:
                    time.sleep(count_time)
                record(pos)

    except KeyboardInterrupt:
        status, aborted = "aborted", True
        print(f"\n^C  Scan aborted at point {len(measured)}/{num_pts} -- "
              "stopping motor and detector...", flush=True)
        raise
    except Exception:
        status = "error"
        raise
    finally:
        # One shared teardown for both detector branches, so they cannot drift.
        # Ctrl+C is switched off for the duration so cleanup always finishes;
        # the inner try/finally guarantees it gets switched back on.
        saved_ctrl_c = disable_ctrl_c()
        try:
            if aborted:
                # Stop motion first: this also clears the interrupted MoveStatus
                # that would otherwise make the return move complain that
                # another set() is still in progress.
                safe_call(f"stop {motor.name}", motor.stop)
                time.sleep(0.2)  # let the motor record settle before re-commanding

            # Only wait for frames on a normal finish. On an abort the remaining
            # frames are never coming, and waiting would stall teardown for up to
            # num_pts * count_time + 10 s with the detector still live.
            if save_img == 1 and not aborted:
                safe_call("wait for frames", wait_for_frames,
                          det, num_pts, num_pts * count_time + 10)

            if kind == "lambda":
                safe_call("stop softglue pulses", softglue.stop_pulses.put, "1!")
            safe_call("stop detector acquisition", det.cam.acquire.put, 0)
            if save_img == 1:
                safe_call("close the HDF file", disarm_hdf, det)
            if beam_control:
                safe_call("block the beam", blockbeam)

            if kind == "eiger":
                safe_call("restore trigger mode",
                          det.cam.trigger_mode.put, "Internal Enable")
                safe_call("restore manual trigger",
                          det.cam.manual_trigger.put, "Disable")
            else:
                safe_call("restore operating mode", det.cam.operating_mode.put, 3)
                safe_call("restore trigger mode", det.cam.trigger_mode.put, 0)
                safe_call("restore softglue preset", softglue_8id_acq.preset.put, 50)
                if beam_control:
                    safe_call("switch the shutter off", shutteroff)

            if aborted:
                print(f"    returning {motor.name} to {start_pos:.5g} ...", flush=True)
            safe_call(f"return {motor.name} to {start_pos:.5g}",
                      motor.move, start_pos, wait=True)

            safe_call("close the SPEC scan block", spec.end_scan, status)

            if save_img == 1:
                safe_call("read the captured-frame count", report_frames_captured, det)
            print(f"SPEC scan {scan_num} finished ({status}), {len(measured)} points.")
        finally:
            restore_ctrl_c(saved_ctrl_c)

    return ScanResult(np.array(measured), columns,
                      scan_num=scan_num, motor_name=motor.name)
