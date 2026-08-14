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

The output follows the layout produced by ``apstools`` ``SpecWriterCallback2``
(``#F/#E/#D/#C``, ``#O``/``#o``, ``#S/#D/#C/#MD/#P/#N/#L``, data rows, closing
``#C`` lines) so that existing SPEC readers -- ``spec2nexus``, PyMca, silx, and
the beamline's MATLAB ``specr`` GUI -- parse it unchanged.

usage (call it directly -- do NOT wrap it in ``RE()``)::

    from id8_common.plans.align.ophyd_scan import dscan_ophyd
    dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)

The companion live viewer is ``specr_py`` (see ``~/bluesky/specr_py/``).

Import order note: like ``scan_8id.py``, this module resolves devices from the
``oregistry`` at import time, so it must be imported after ``make_devices()``
has run in ``startup.py``.
"""

import datetime
import getpass
import math
import os
import socket
import time
import uuid

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

#: Diffractometer axes recorded in the SPEC ``#O``/``#P`` control lines.
HUBER_AXES = ("nu", "delta", "mu", "eta", "chi", "phi", "x", "y", "z")


# =============================================================================
# SPEC file writer
# =============================================================================


def _fmt(value):
    """Format one data cell so that ``float()`` will always parse it.

    The spec2nexus parser does ``float(val)`` on every cell, so a non-finite
    value or a string would break the whole scan block. Substitute 0 for
    nan/inf rather than emitting a token some parsers reject.
    """
    if isinstance(value, bool):
        return "1" if value else "0"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "0"
    if not math.isfinite(v):
        return "0"
    return f"{v:.12g}"


def _now():
    return datetime.datetime.now().strftime(SPEC_TIME_FORMAT)


class SpecFile:
    """
    Minimal SPEC data file writer.

    Standard library only -- no apstools, no bluesky.

    The format rules below are not stylistic; they come from the ``spec2nexus``
    parser and the MATLAB ``specr`` reader. Getting any of them wrong produces a
    file that looks correct and silently mis-parses:

    * ``#L`` labels are separated by **two** spaces. Both parsers split on a
      ``\\s{2,}`` pattern, which is what allows a label to contain a single
      space (e.g. ``SMS granite`` in the reference file).
    * ``#N`` must equal ``len(labels)``.
    * ``#S`` needs the trailing space -- scan detection is
      ``line.startswith("#S ")`` in both readers.
    * A blank line must precede each ``#S``.
    * Every data cell must parse as a float.

    Every write opens the file, appends complete newline-terminated lines, and
    closes it. A single append-mode write of a short line is atomic for a
    regular file, so a polling reader never observes a torn row.
    """

    def __init__(self, path):
        self.path = str(path)
        self.labels = []
        self._points = 0

    # -- internals ---------------------------------------------------------

    def _append(self, lines):
        """Append complete lines. Never raises into the caller."""
        try:
            directory = os.path.dirname(self.path)
            if directory:
                os.makedirs(directory, exist_ok=True)
            with open(self.path, "a") as f:
                f.write("".join(line + "\n" for line in lines))
        except Exception as exc:  # never break a scan over a log file
            print(f"WARNING: SPEC write failed ({exc}); continuing scan.")

    def exists(self):
        return os.path.exists(self.path) and os.path.getsize(self.path) > 0

    # -- public API --------------------------------------------------------

    def write_file_header(self, motor_names=(), force=False):
        """Write the ``#F``/``#E``/``#D``/``#C``/``#O``/``#o`` block.

        Skipped when the file already has content, so one file accumulates
        scans under a single header (unless ``force``).
        """
        if self.exists() and not force:
            return
        now = time.time()
        user = getpass.getuser() or "8idUser"
        host = socket.gethostname() or "localhost"
        lines = [
            f"#F {self.path}",
            f"#E {int(now)}",
            f"#D {datetime.datetime.fromtimestamp(now).strftime(SPEC_TIME_FORMAT)}",
            f"#C Bluesky  user = {user}  host = {host}",
        ]
        lines += _numbered_lines("#O", motor_names)
        lines += _numbered_lines("#o", motor_names)
        self._append(lines)

    def start_scan(
        self,
        scan_num,
        command,
        labels,
        metadata=None,
        motor_positions=(),
        comments=(),
    ):
        """Open a new ``#S`` block. ``labels`` defines the columns."""
        self.labels = list(labels)
        self._points = 0
        stamp = _now()
        lines = [
            "",  # blank line before #S is required
            f"#S {int(scan_num)}  {command}",
            f"#D {stamp}",
        ]
        lines += [f"#C {stamp}.  {c}" for c in comments]
        for key, value in (metadata or {}).items():
            lines.append(f"#MD {key} = {value}")
        lines += _numbered_lines("#P", [_fmt(v) for v in motor_positions]) or ["#P0 "]
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
        """Write the closing ``#C`` lines for the scan block."""
        stamp = _now()
        self._append(
            [
                f"#C {stamp}.  num_events_primary = {self._points}",
                f"#C {stamp}.  exit_status = {status}",
            ]
        )


def _numbered_lines(tag, items, per_line=8):
    """SPEC wraps #O/#o/#P at 8 entries per line: ``#O0 a  b``, ``#O1 c  d``."""
    items = list(items)
    if not items:
        return []
    out, row = [], 0
    while items:
        out.append(f"{tag}{row} " + "  ".join(str(v) for v in items[:per_line]))
        items = items[per_line:]
        row += 1
    return out


# =============================================================================
# Paths, scan numbering, detector plumbing
# =============================================================================


#: Where the SPEC files are written. This is deliberately NOT the GPFS data
#: directory: /gdata is mounted on the acquisition host (amber) but not on the
#: analysis workstations (e.g. kouga), whereas ~/bluesky is shared by both. The
#: detector .h5 files still go to GPFS -- only the small text file lives here,
#: so the live viewer can be run from anywhere.
SPEC_DIR = os.path.expanduser("~/bluesky")


def default_spec_path():
    """``~/bluesky/{experiment_name}.spec``."""
    exp = pv_registers.experiment_name.get().strip()
    return os.path.join(SPEC_DIR, f"{exp}.spec")


def image_file_path():
    """Directory the detector HDF5 files go to.

    Matches ``save_images()`` in scan_8id.py, which builds
    ``/gdata/dm/8ID/8IDE/<cycle>/<experiment>/data/bluesky`` -- reconstructed
    here from the registers rather than hard-coding the cycle.
    """
    mount = pv_registers.mount_point.get().strip()
    cycle = pv_registers.cycle_name.get().strip()
    exp = pv_registers.experiment_name.get().strip()
    return f"{mount}{cycle}/{exp}/data/bluesky"


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


def huber_motor_table():
    """Return ``(names, positions)`` for the diffractometer ``#O``/``#P`` lines.

    Read-only: touches only ``.position`` (the motor ``.RBV``). Degrades to
    empty lists if the device is absent or a read fails, so a missing IOC can
    never stop a scan.
    """
    names, positions = [], []
    try:
        huber = oregistry["huber"]
    except Exception:
        return names, positions
    for axis in HUBER_AXES:
        try:
            motor = getattr(huber, axis)
            positions.append(motor.position)
            names.append(motor.name)
        except Exception:
            continue
    return names, positions


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
    """[(column_label, Signal)] to record at each point.

    Ordered per SPEC convention: the most important detector goes last, so a
    reader defaulting the Y axis to the final column plots stats1.
    """
    if detector_kind(det) == "tetramm":
        pairs = [(f"current{n}", getattr(det, f"current{n}").mean_value) for n in (2, 3, 4)]
        pairs += [("sum_all", det.sum_all.mean_value), ("current1", det.current1.mean_value)]
        return [(f"{det.name}_{lab}_mean_value", sig) for lab, sig in pairs]
    # stats1 last: "first detector moved to last column per SPEC convention"
    return [
        (f"{det.name}_stats{n}_total", getattr(det, f"stats{n}").total)
        for n in (2, 3, 1)
    ]


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
    set_attenuation=True,
    beam_control=True,
):
    """
    Relative single-motor scan, Ophyd only, with live SPEC output.

    Mirrors the hardware choreography of ``dscan`` in scan_8id.py (pre-armed
    software trigger for the Eiger, softglue-pulsed external trigger for the
    Lambda) but contains no Bluesky.

    Call directly::

        dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)

    args:
        motor: ophyd positioner (e.g. huber.delta)
        rel_begin, rel_end: start/end relative to the current position
        num_pts: number of points
        count_time: detector acquisition time per point (s)
        det: eiger4M (default) or lambda2M
        att_ratio: attenuation ratio, applied only if ``set_attenuation``
        save_img: 1 to write detector images, 0 to scan without saving
        spec_path: override the SPEC file location
        set_attenuation: if False, leave the attenuator exactly as it is. Use
            when the current attenuation is already correct and you do not want
            the scan changing the beam condition.
        beam_control: if False, skip ``pre_align``/``PIND_status`` and the
            shutter open/close. Use for a dry run that touches only the scanned
            motor and the detector.

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

    if beam_control:
        pre_align()
        PIND_status(0)
    if set_attenuation:
        att(att_ratio)

    folder_prefix = gen_folder_prefix() if save_img == 1 else ""
    scan_num = scan_number_from_prefix(folder_prefix)
    uid = str(uuid.uuid4())

    # --- SPEC file -------------------------------------------------------
    spec = SpecFile(spec_path or default_spec_path())
    motor_names, motor_positions = huber_motor_table()
    spec.write_file_header(motor_names)

    sigs = signal_map(det)
    # Readback first (that is what SPEC readers plot), commanded value second.
    # Logging only the commanded value would draw a perfect ramp even if the
    # motor never moved, so both are recorded -- same convention as apstools.
    labels = [motor.name, f"{motor.name}_setpoint", "Epoch", "Epoch_float"]
    labels += [lab for lab, _ in sigs]
    command = (
        f"dscan_ophyd({motor.name}, {rel_begin}, {rel_end}, {num_pts}, "
        f"{count_time}, det={det.name})"
    )
    metadata = {
        "uid": uid,
        "beamline_id": "8-ID-E",
        "login_id": f"{getpass.getuser()}@{socket.gethostname()}",
        "pid": os.getpid(),
        "scan_type": "dscan_ophyd",
        "detectors": [det.name],
        "motors": [motor.name],
        "num_points": num_pts,
        "count_time": count_time,
        "acquisition": "ophyd (no bluesky)",
    }
    if folder_prefix:
        metadata["image_file"] = f"{folder_prefix}.h5"
    spec.start_scan(
        scan_num,
        command,
        labels,
        metadata=metadata,
        motor_positions=motor_positions,
        comments=[f"plan_type = function", f"uid = {uid}"],
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
    measured, columns = [], {lab: [] for lab, _ in sigs}
    status = "success"

    def record(setpoint):
        """Read the signals for one point and append a SPEC row.

        ``setpoint`` is where the motor was told to go; the first column is
        where it actually is.
        """
        elapsed = time.time() - t_start
        actual = motor.position
        values = [sig.get() for _, sig in sigs]
        spec.add_point([actual, setpoint, round(elapsed), elapsed] + values)
        measured.append(actual)
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
            if beam_control:
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
                if beam_control:
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
                det.hdf1.capture.put(1)  # after acquire, exactly as dscan() does
            if beam_control:
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
                if beam_control:
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
            if beam_control:
                shutteroff()
        if beam_control:
            blockbeam()
        spec.end_scan(status)
        if save_img == 1:
            print(f"# images captured: {det.hdf1.num_captured.get()}")
        print(f"SPEC scan {scan_num} finished ({status}), {len(measured)} points.")

    return np.array(measured), columns
