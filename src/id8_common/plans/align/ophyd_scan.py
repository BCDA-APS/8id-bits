"""
Ophyd-only scanning with live SPEC-format output for the 8-ID beamlines.

**No Bluesky.** This module contains no RunEngine, no plans, no generators, no
``plan_stubs``, and emits no documents. Everything is plain Python driving Ophyd
signals with ``.put()`` / ``.get()`` / ``.move()``.

The scan writes one row to a SPEC-format ``.dat`` file per point, closing the file
after every row, so an external plotting script can poll the file while the scan
runs. Because the file is append-only and never held open, a reader can never
collide with the writer -- which is the whole reason for doing it this way rather
than reading the detector ``.h5`` mid-scan.

usage (call it directly -- do NOT wrap it in ``RE()``)::

    from id8_common.plans.align.ophyd_scan import dscan_ophyd
    dscan_ophyd(huber.x, -0.5, 0.5, 41, 1.0, det=lambda2M)

The companion reader is a separate process; see ``spec_writer_plan.md``.

Import order note: like ``scan_8id.py``, this module resolves devices from the
``oregistry`` at import time, so it must be imported after ``make_devices()`` has
run in ``startup.py``.
"""

import math
import os
import socket
import getpass
import time
import datetime

import numpy as np

from apsbits.core.instrument_init import oregistry

# All of these are plain functions (not Bluesky plans) -- verified: neither
# shutter_att.py nor ad_acq.py imports bluesky, so this import chain stays clean.
from id8_common.plans.set.shutter_att import showbeam
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import shutteron
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.plans.set.shutter_att import pre_align
from id8_common.plans.set.shutter_att import PIND_status
from id8_common.plans.set.shutter_att import att
from id8_common.plans.acquire.ad_acq import gen_folder_prefix

pv_registers = oregistry["pv_registers"]
softglue = oregistry["softglue"]
softglue_8id_acq = oregistry["softglue_8id_acq"]
eiger4M = oregistry["eiger4M"]
lambda2M = oregistry["lambda2M"]

SPEC_TIME_FORMAT = "%a %b %d %H:%M:%S %Y"


# =============================================================================
# SPEC file writer
# =============================================================================


def _fmt(value):
    """Format one data cell so that ``float()`` will always parse it.

    The spec2nexus parser does ``float(val)`` on every cell, so a non-finite
    value or a string would break the whole scan block. Substitute 0 for
    nan/inf rather than emitting a token some parsers reject.
    """
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "0"
    if not math.isfinite(v):
        return "0"
    return f"{v:.12g}"


class SpecFile:
    """
    Minimal SPEC data file writer.

    Standard library only -- no apstools, no bluesky.

    The format rules below are not stylistic; they come from the ``spec2nexus``
    parser and getting any of them wrong produces a file that looks correct and
    silently mis-parses:

    * ``#L`` labels are separated by **two** spaces. The parser splits on a
      ``\\s{2,}`` regex; a single space only works via a fallback that requires
      ``#N`` to match the column count exactly.
    * ``#N`` must equal ``len(labels)``.
    * ``#S`` needs the trailing space -- scan detection is
      ``line.startswith("#S ")``.
    * A blank line must precede each ``#S``.
    * Every data cell must parse as a float.

    Every write opens the file, appends one or more complete newline-terminated
    lines, and closes it. A single append-mode write of a short line is atomic
    for a regular file, so a polling reader never observes a torn row.
    """

    def __init__(self, path):
        self.path = str(path)
        self.labels = []
        self._points = 0

    # -- internals ---------------------------------------------------------

    def _append(self, lines):
        """Append complete lines. Never raises into the caller."""
        try:
            os.makedirs(os.path.dirname(self.path) or ".", exist_ok=True)
            with open(self.path, "a") as f:
                f.write("".join(line + "\n" for line in lines))
        except Exception as exc:  # never break a scan over a log file
            print(f"WARNING: SPEC write failed ({exc}); continuing scan.")

    def _write_file_header(self):
        """Write ``#F``/``#E``/``#D``/``#C`` once, only for a brand-new file."""
        if os.path.exists(self.path) and os.path.getsize(self.path) > 0:
            return
        now = time.time()
        stamp = datetime.datetime.fromtimestamp(now).strftime(SPEC_TIME_FORMAT)
        user = getpass.getuser() or "8idUser"
        host = socket.gethostname() or "localhost"
        self._append(
            [
                f"#F {self.path}",
                f"#E {int(now)}",
                f"#D {stamp}",
                f"#C Ophyd (no bluesky)  user = {user}  host = {host}",
            ]
        )

    # -- public API --------------------------------------------------------

    def start_scan(self, scan_num, command, labels, comments=()):
        """Open a new ``#S`` block. ``labels`` defines the columns."""
        self._write_file_header()
        self.labels = list(labels)
        self._points = 0
        stamp = datetime.datetime.now().strftime(SPEC_TIME_FORMAT)
        lines = [
            "",  # blank line before #S is required
            f"#S {int(scan_num)} {command}",
            f"#D {stamp}",
        ]
        lines += [f"#C {c}" for c in comments]
        lines += [
            "#N " + str(len(self.labels)),
            "#L " + "  ".join(self.labels),  # TWO spaces -- see class docstring
        ]
        self._append(lines)

    def add_point(self, values):
        """Append one data row. Length must match the labels given to start_scan."""
        vals = list(values)
        if len(vals) != len(self.labels):
            print(
                f"WARNING: SPEC row has {len(vals)} values but {len(self.labels)} "
                "labels; row skipped."
            )
            return
        self._append([" ".join(_fmt(v) for v in vals)])
        self._points += 1

    def end_scan(self, status="success"):
        """Close out the scan block."""
        stamp = datetime.datetime.now().strftime(SPEC_TIME_FORMAT)
        self._append(
            [
                f"#C {stamp}.  num_points = {self._points}",
                f"#C {stamp}.  exit_status = {status}",
            ]
        )


# =============================================================================
# Paths, scan numbering, detector plumbing
# =============================================================================


def default_spec_path():
    """``<mount><cycle>/<experiment>/data/bluesky/<experiment>.dat``.

    Mirrors the layout used by ``get_common_file_path`` in ad_acq.py, so the
    SPEC file lands beside the detector ``.h5`` files.
    """
    mount = pv_registers.mount_point.get().strip()
    cycle = pv_registers.cycle_name.get().strip()
    exp = pv_registers.experiment_name.get().strip()
    return f"{mount}{cycle}/{exp}/data/bluesky/{exp}.dat"


def scan_number_from_prefix(folder_prefix, fallback=0):
    """``"A0092_G10_a0007"`` -> ``92``.

    Reuses the number ``gen_folder_prefix()`` already burned, so SPEC ``#S 92``
    lines up with ``A0092_*.h5``. Do not read ``measurement_num`` again here --
    ``gen_folder_prefix()`` increments it as a side effect.
    """
    try:
        return int("".join(c for c in folder_prefix.split("_")[0] if c.isdigit()))
    except (ValueError, IndexError):
        return fallback


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


def signal_map(det):
    """[(column_label, Signal)] to record at each point."""
    if detector_kind(det) == "tetramm":
        pairs = [("sum_all", det.sum_all.mean_value)]
        pairs += [(f"current{n}", getattr(det, f"current{n}").mean_value) for n in (1, 2, 3, 4)]
        return [(f"{det.name}_{lab}_mean_value", sig) for lab, sig in pairs]
    return [
        (f"{det.name}_stats{n}_total", getattr(det, f"stats{n}").total)
        for n in (1, 2, 3)
    ]


def arm_hdf(det, file_path, file_name, num_capture, start_capture=True):
    """
    Prepare the HDF5 plugin for capture.

    This replaces what Bluesky's ``stage_wrapper`` would have done. Verified in
    ``apstools.devices.area_detector_support.AD_EpicsHdf5FileName``: staging sets
    ``file_write_mode="Stream"`` and ``capture=1`` (plus ``enable``/``auto_save``
    from the ophyd FileStore mixins). Without Bluesky nothing applies those, so
    they are set explicitly here -- omit them and the scan runs but writes no file.

    ``start_capture`` exists because the two detectors want opposite ordering, and
    scan_8id.py is the reference for both. The Eiger is staged (capture on) before
    ``cam.acquire``; the Lambda sets ``cam.acquire`` first and only then turns
    capture on. Pass False to leave capture off and switch it on at the call site.
    """
    det.hdf1.enable.put(1)
    det.hdf1.file_path.put(file_path)
    det.hdf1.file_name.put(file_name)
    det.hdf1.file_write_mode.put(2)  # Stream
    det.hdf1.auto_save.put(1)
    det.hdf1.auto_increment.put(1)
    det.hdf1.num_capture.put(num_capture)
    det.cam.array_callbacks.put(1)  # plugins receive frames
    if start_capture:
        det.hdf1.capture.put(1)  # must be last


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
    att_ratio=1e6,
    save_img=1,
    spec_path=None,
):
    """
    Relative single-motor scan, Ophyd only, with live SPEC output.

    Mirrors the hardware choreography of ``dscan`` in scan_8id.py (pre-armed
    software trigger for the Eiger, softglue-pulsed external trigger for the
    Lambda) but contains no Bluesky.

    Call directly::

        dscan_ophyd(huber.x, -0.5, 0.5, 41, 1.0, det=lambda2M)

    args:
        motor: ophyd positioner (e.g. huber.x)
        rel_begin, rel_end: start/end relative to the current position
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: eiger4M (default) or lambda2M
        att_ratio: attenuation ratio
        save_img: 1 to write detector images, 0 to scan without saving
        spec_path: override the SPEC file location

    returns:
        (positions, {label: [values]}) for the points actually measured.
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

    pre_align()
    att(att_ratio)
    PIND_status(0)

    folder_prefix = gen_folder_prefix() if save_img == 1 else ""
    scan_num = scan_number_from_prefix(folder_prefix)

    # --- SPEC file -------------------------------------------------------
    spec = SpecFile(spec_path or default_spec_path())
    sigs = signal_map(det)
    labels = [motor.name, "Epoch", "Epoch_float"] + [lab for lab, _ in sigs]
    command = (
        f"dscan_ophyd({motor.name}, {rel_begin}, {rel_end}, {num_pts}, {count_time}, "
        f"det={det.name})"
    )
    comments = [f"detector = {det.name}", f"attenuation_ratio = {att_ratio}"]
    if folder_prefix:
        comments.append(f"image_file = {folder_prefix}.h5")
    spec.start_scan(scan_num, command, labels, comments)
    print(f"SPEC file: {spec.path}   (#S {scan_num})")

    # --- detector / file setup -------------------------------------------
    if save_img == 1:
        file_path = os.path.dirname(spec.path)
        print(f"Scan folder created: {folder_prefix}")
        print(f"File path: {file_path}")
        # Eiger: capture on before cam.acquire. Lambda: after -- see arm_hdf().
        arm_hdf(det, file_path, folder_prefix, num_pts, start_capture=(kind == "eiger"))

    start_pos = motor.position
    positions = np.linspace(start_pos + rel_begin, start_pos + rel_end, num_pts)
    t_start = time.time()
    measured, columns = [], {lab: [] for lab, _ in sigs}
    status = "success"

    def record(pos):
        """Read the signals for one point and append a SPEC row."""
        now = time.time()
        elapsed = now - t_start
        values = [sig.get() for _, sig in sigs]
        spec.add_point([pos, round(elapsed), elapsed] + values)
        measured.append(pos)
        for (lab, _), v in zip(sigs, values):
            columns[lab].append(v)

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
            showbeam()
            try:
                for pos in positions:
                    motor.move(pos, wait=True)
                    det.cam.special_trigger_button.put(1, wait=False)
                    time.sleep(count_time)
                    record(pos)
            finally:
                if save_img == 1:
                    wait_for_frames(det, num_pts, num_pts * count_time + 10)
                det.cam.acquire.put(0)
                blockbeam()
                motor.move(start_pos, wait=True)

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
                det.hdf1.capture.put(1)
            shutteron()
            showbeam()
            try:
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
            finally:
                softglue.stop_pulses.put("1!")
                det.cam.acquire.put(0)
                blockbeam()
                motor.move(start_pos, wait=True)

    except KeyboardInterrupt:
        status = "aborted"
        raise
    except Exception:
        status = "error"
        raise
    finally:
        if save_img == 1:
            disarm_hdf(det)
        if kind == "eiger":
            det.cam.trigger_mode.put("Internal Enable")
            det.cam.manual_trigger.put("Disable")
        else:
            det.cam.operating_mode.put(3)
            det.cam.trigger_mode.put(0)
            softglue_8id_acq.preset.put(50)
            shutteroff()
        blockbeam()
        spec.end_scan(status)
        if save_img == 1:
            print(f"# images captured: {det.hdf1.num_captured.get()}")
        print(f"SPEC scan {scan_num} finished ({status}), {len(measured)} points.")

    return np.array(measured), columns
