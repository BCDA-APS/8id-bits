"""
Write a scan to a plain, flat CSV file -- one file per scan.

This module knows about *file format only*. It has no scan logic and moves no
hardware, so the motor-control code in :mod:`ophyd_scan` stays readable.

What the file looks like::

    #start_time,2026-08-26 14:32:07
    #h5_file,/gdata/dm/8ID/8IDE/2026-3/comm202609/data/bluesky/S01459_HuberDelta_Lambda2M-a1345678-1s.h5
    ... more one-per-line entries ...
    #DATA
    ,huber_delta,elapsed_time,lambda2M_stats1_total
    0,31.0,0.0,1234
    1,31.025,1.05,1301
    ...
    #END
    #exit_status,success
    #points_written,41

* Everything **above** the ``#DATA`` marker is a ``#label,value`` line: things
  that do not change during the scan (start time, the ``.h5`` file name, motor
  positions at the start, ROI settings, ...). Which lines appear is set by
  ``configs/scan_csv_template.yml``, not by this module.
* The line **after** the marker holds the column names, led by an empty one.
* Everything after that is one row per scan point, led by the point index.
* A bare ``#END`` closes the table; how the scan ended follows it.

**This is BLUETELLA's (9-ID's) own extended CSV format, deliberately.** It was
adopted on 2026-09-07 so the 9-ID viewers -- ``scanviewer.py`` (1D) and
``meshviewer.py`` (2D) -- can open an 8-ID scan with no 8-ID-specific code in
them, which is what makes those viewers' 8-ID support upstreamable instead of a
fork we have to maintain. Two things changed to get there, and only two: the
header labels gained a leading ``#``, and ``#END,<status>,<n>`` became a bare
``#END`` with the outcome on the lines after it. No header line was dropped,
moved, or sent to a sidecar file. See ``ScanCsv.close`` for why the terminator
mattered and the ``#`` did not.

The first column is the **point index**, 0-based, and its name is deliberately
blank -- that is what ``pandas.to_csv()`` writes for an index, so ``read_csv``
with ``index_col=0`` picks it up as the index rather than as data. See
:data:`INDEX_LABEL`.

Reading one back is three lines of Python::

    import pandas
    rows = pandas.read_csv(path, skiprows=lambda i: i <= marker_line,
                           comment="#", index_col=0)

or, without pandas, see ``read_scan_csv()`` at the bottom of this file.

**What goes in the file is not hard-coded here.** The header lines and the data
columns both come from ``configs/scan_csv_template.yml``. To add a column or
record another PV in the header, edit that YAML file -- no Python change needed.

The file is opened, appended to, and closed once per scan point, so an external
viewer can poll it while the scan runs without ever colliding with the writer.
That is the whole reason a scan writes a text file instead of having the viewer
try to read the detector ``.h5`` mid-acquisition.

Note on loading: this uses ``yaml.safe_load`` rather than
``apsbits.utils.config_loaders.load_config``. That helper merges whatever it
reads into a process-global ``_iconfig`` and overwrites ``ICONFIG_PATH`` /
``INSTRUMENT_PATH``, so pointing it at this template would corrupt the session's
real instrument config.
"""

import csv
import numbers
import os
import time
from pathlib import Path

import yaml

from id8_common.registry import oregistry

#: Shipped default. The repo is an editable install, so users edit it in place.
DEFAULT_TEMPLATE = Path(__file__).resolve().parents[2] / "configs" / "scan_csv_template.yml"

#: How many scanned motors a template may name. Four is well past anything at
#: 8-ID (a mesh moves two) and keeps ``{motor5}`` an honest error rather than a
#: silently blank column.
MAX_MOTORS = 4

#: Column sources the scan computes itself instead of reading from a PV.
#: ``motor`` and ``motor_setpoint`` are aliases for ``motor1``, so a template
#: written when only one-motor scans existed still means exactly what it meant.
SPECIAL_SOURCES = (
    ("motor", "motor_setpoint", "elapsed", "epoch")
    + tuple(f"motor{i}" for i in range(1, MAX_MOTORS + 1))
    + tuple(f"motor{i}_setpoint" for i in range(1, MAX_MOTORS + 1))
)
# With MAX_MOTORS = 4 that is: "motor", "motor_setpoint", "elapsed", "epoch",
# "motor1".."motor4", and "motor1_setpoint".."motor4_setpoint".

#: Name of the leading point-index column in the data block. Blank on purpose:
#: it is the pandas index convention, so ``read_csv(..., index_col=0)`` treats
#: it as the index. Not a template entry -- every scan file gets it, and the
#: template's ``columns:`` stay the real measured data.
INDEX_LABEL = ""

#: How wide one column is *on the screen*. The .csv file is unaffected.
SCREEN_WIDTH = 13

#: Decimal places *on the screen*. The .csv keeps full precision.
SCREEN_DECIMALS = 3


# =============================================================================
# Loading the template
# =============================================================================


def template_path(explicit=None):
    """Which template to use: argument, then $ID8_SCAN_CSV_TEMPLATE, then default."""
    if explicit:
        return Path(os.path.expanduser(str(explicit)))
    from_env = os.environ.get("ID8_SCAN_CSV_TEMPLATE")
    if from_env:
        return Path(os.path.expanduser(from_env))
    return DEFAULT_TEMPLATE


def load_template(path=None):
    """Read the YAML template. Returns (config_dict, resolved_path)."""
    resolved = template_path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"CSV template not found: {resolved}")
    with open(resolved, "r") as handle:
        config = yaml.safe_load(handle) or {}
    if not isinstance(config, dict):
        raise ValueError(f"CSV template {resolved} must be a mapping at the top level")
    return config, resolved


# =============================================================================
# Turning template text into live signals and values
# =============================================================================


def _split_motor_index(source):
    """``"motor2_setpoint"`` -> ``("motor_setpoint", 1)``; ``"motor"`` -> ``("motor", 0)``.

    Indexed names, NOT bracket syntax: ``str.format`` already owns brackets, so
    ``"{motor[1]}".format(motor="huber_x")`` silently yields ``"u"`` instead of
    erroring -- a whole column of single characters.
    """
    base, _, tail = str(source).partition("_")
    if not base.startswith("motor"):
        return source, 0
    digits = base[len("motor") :]
    if digits.isdigit():
        index = int(digits) - 1  # "motor2" is the second motor, i.e. index 1
    else:
        index = 0  # a bare "motor" means the first one
    if tail == "setpoint":
        special = "motor_setpoint"
    else:
        special = "motor"
    return special, index


def resolve_signal(dotted, motors, det):
    """ "det.stats1.total" or "tetramm1.current1.mean_value" -> the ophyd object.

    ``det.`` and ``motor.`` mean the objects passed to *this* scan, which is what
    lets one template serve eiger4M, lambda2M and tetramm1.

    ``motors`` is the tuple of scanned positioners. ``motor.`` is the first one;
    ``motor1.``/``motor2.``/... name them individually, so a header line can
    record e.g. ``motor2.user_offset``.
    """
    head, _, rest = str(dotted).partition(".")
    if head == "det":
        if det is None:
            raise KeyError(f"{dotted!r} needs a detector but this scan has none")
        obj = det
    elif head == "motor" or (head.startswith("motor") and head[len("motor") :].isdigit()):
        # "motor" is the first scanned motor, "motor2" the second, and so on.
        digits = head[len("motor") :]
        if digits.isdigit():
            index = int(digits) - 1
        else:
            index = 0
        if index >= len(motors):
            raise KeyError(f"{dotted!r}: this scan has only {len(motors)} motor(s)")
        obj = motors[index]
    else:
        try:
            obj = oregistry[head]
        except Exception as exc:
            raise KeyError(f"no device {head!r} in the oregistry ({dotted!r})") from exc
    for part in rest.split("."):
        if not part:
            continue
        try:
            obj = getattr(obj, part)
        except AttributeError as exc:
            raise KeyError(f"{dotted!r}: {obj.name!r} has no attribute {part!r}") from exc
    return obj


def substitute(text, context):
    """Apply the {motor} / {det} / {command} ... replacements, with a clear error."""
    try:
        return str(text).format(**context)
    except KeyError as exc:
        raise KeyError(f"unknown substitution {exc} in {text!r}; available: {', '.join(sorted(context))}") from exc


def read_value(obj):
    """Positioners have .position; plain signals only have .get()."""
    if hasattr(obj, "position"):
        return obj.position
    return obj.get()


def cell(value):
    """Format one value for a CSV cell, keeping integers as integers."""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, numbers.Integral):
        return str(int(value))
    if isinstance(value, numbers.Real):
        return repr(float(value))  # shortest text that reads back exactly
    return str(value)


def one_line(text):
    """Flatten newlines so a multi-line value cannot break the file structure."""
    return " ".join(str(text).split())


def show(value):
    """Format one value for the *screen*.

    Three decimals is all anyone reads off a scrolling terminal, and the full
    ``repr(float)`` -- ``29.500603937599992`` -- is what made the columns run
    into each other. The file still gets every digit; see :func:`cell`.
    """
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, numbers.Integral):
        return str(int(value))
    if isinstance(value, numbers.Real):
        return f"{float(value):.{SCREEN_DECIMALS}f}"
    return one_line(value)


def short_label(label):
    """Fit a column name into :data:`SCREEN_WIDTH`, keeping the *end*.

    ``lambda2M_stats1_total`` -> ``.stats1_total``: the tail is the part that
    tells two columns apart. The leading dot means "trimmed". The .csv keeps
    the full name.
    """
    if len(label) <= SCREEN_WIDTH:
        return label
    return "." + label[-(SCREEN_WIDTH - 1) :]


# =============================================================================
# One column of the data table
# =============================================================================


class Column:
    """A data column: either something the scan computes, or a signal to read."""

    def __init__(self, label, special=None, signal=None, motor_index=0):
        """Store the column name and where its value comes from."""
        self.label = label
        self.special = special
        self.signal = signal
        #: which scanned motor this column follows; 0 for "motor"/"motor1"
        self.motor_index = motor_index

    def value(self, positions, setpoints, elapsed):
        """Return this column's value for the point being written.

        ``positions`` and ``setpoints`` are one entry per scanned motor, in the
        order they were given to :func:`open_scan`.
        """
        if self.special == "motor":
            return positions[self.motor_index]
        if self.special == "motor_setpoint":
            return setpoints[self.motor_index]
        if self.special == "elapsed":
            return elapsed
        if self.special == "epoch":
            return time.time()
        return self.signal.get()


# =============================================================================
# The file
# =============================================================================


class ScanCsv:
    """One scan, one CSV file. Built by :func:`open_scan`, never directly."""

    def __init__(self, path, header_lines, marker, columns, motors=(), verbose=True):
        """Store the layout. Nothing is written until write_header() is called."""
        self.path = str(path)
        self.header_lines = header_lines  # list of [label, v1, v2, ...]
        self.marker = marker
        self.columns = columns  # list of Column
        self.labels = [c.label for c in columns]
        self.verbose = verbose
        self.t0 = time.time()
        self.num_points = 0
        #: the scanned positioners, in the order given to open_scan(). Held here
        #: rather than attached by open_scan() afterwards, so a ScanCsv built
        #: directly is usable -- it used to raise AttributeError on _motor at the
        #: first add_point().
        self._motors = tuple(motors)
        #: every value read so far, {label: [...]} -- handy at the prompt
        self.data = {label: [] for label in self.labels}

    @property
    def _motor(self):
        """The first scanned motor. Kept so one-motor code reads as it did."""
        return self._motors[0] if self._motors else None

    # -- writing ----------------------------------------------------------

    def _append(self, rows):
        """Append rows and close immediately, so a reader can never collide."""
        folder = os.path.dirname(self.path)
        if folder:
            os.makedirs(folder, exist_ok=True)
        with open(self.path, "a", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            for row in rows:
                writer.writerow(row)

    def write_header(self):
        """Header lines, the marker, and the column-name line.

        Refuses to write into a file that already exists. _append() opens in
        mode "a" -- which is right for the per-point appends that let a viewer
        poll a running scan, but means that without this check a second scan
        landing on the same name silently CONCATENATES itself onto the first:
        two headers, two ``#END`` lines, one file. Nothing errors, and every
        reader is then wrong -- ``read_scan_csv`` and BLUETELLA both stop at the
        first ``#END``, so they report the OLD scan's data and the OLD scan's
        exit status for a file the new scan is still writing.

        A duplicate name means the measurement counter has repeated, and the
        ``.h5`` is about to collide too, so this is worth stopping for. It runs
        before the scan has moved anything or opened the shutter.

        That counter (``expt.measurement_num``, the EPICS register
        ``8ideSoft:Reg1``) is shared with ``det_acq_series``, which names files
        under ``data/`` rather than the ``data/bluesky/`` these scans write to.
        So when working out which number is already taken, look in both.
        """
        if os.path.exists(self.path):
            raise FileExistsError(
                f"{self.path} already exists -- refusing to append a second scan to it.\n"
                f"The measurement number has repeated, so the .h5 for this scan would "
                f"collide too. Check expt.measurement_num (caget 8ideSoft:Reg1) before "
                f"scanning again; move or delete the old file if it is not wanted."
            )
        # Each header label carries a leading "#", which is what makes this
        # file the SAME extended format BLUETELLA (9-ID) already writes and
        # already reads: "#label,value" lines, "#DATA", the table, "#END".
        # Nothing is lost or moved by the "#" -- every line stays in this file,
        # in this order -- and a reader that skips to "#DATA" never looks at
        # them either way. It is what lets the 9-ID viewers open an 8-ID scan
        # with no 8-ID-specific code in them.
        rows = [["#" + str(line[0])] + list(line[1:]) for line in self.header_lines]
        rows.append([self.marker])
        rows.append([INDEX_LABEL] + self.labels)
        self._append(rows)
        if self.verbose:
            print(f"Scan file: {self.path}")
            print("   # " + " ".join(f"{short_label(x):>{SCREEN_WIDTH}s}" for x in self.labels), flush=True)

    def add_point(self, *setpoints):
        """Read every column once and append a row. Returns the values.

        One setpoint per scanned motor, in the order given to :func:`open_scan`:
        ``scan.add_point(pos)`` for a dscan or ascan, ``scan.add_point(p1, p2)``
        for a d2scan, a2scan, dmesh or mesh.
        """
        elapsed = time.time() - self.t0
        if len(setpoints) != len(self._motors):
            raise ValueError(
                f"add_point() got {len(setpoints)} setpoint(s) but this scan has "
                f"{len(self._motors)} motor(s): {[m.name for m in self._motors]}"
            )
        # Where each motor actually ended up, which is never exactly the
        # setpoint it was sent to. A None motor is only possible for a ScanCsv
        # built by hand -- open_scan() always has real positioners -- and then
        # the setpoint stands in for the missing readback.
        positions = []
        for motor, setpoint in zip(self._motors, setpoints, strict=True):
            if motor is not None:
                positions.append(read_value(motor))
            else:
                positions.append(setpoint)
        values = [c.value(positions, setpoints, elapsed) for c in self.columns]
        # num_points is still the count *before* this row, so the index is
        # 0-based like a pandas index. (The screen counter below is 1-based --
        # it says how far along the scan is, which reads oddly starting at 0.)
        self._append([[self.num_points] + [cell(v) for v in values]])
        self.num_points += 1
        for index, column in enumerate(self.columns):
            self.data[column.label].append(values[index])
        if self.verbose:
            # Separated by a real space, so a value wider than the column
            # pushes the row out of line instead of merging with its neighbour.
            print(
                f"{self.num_points:>4d} " + " ".join(f"{show(v):>{SCREEN_WIDTH}s}" for v in values),
                flush=True,
            )
        return values

    def close(self, status="success"):
        """Write the closing lines. Safe to call twice.

        The terminator is a BARE ``#END``, with the outcome on separate lines
        after it. That is not cosmetic: ``#END`` is the last line a reader
        parses as part of the table, and the extra fields of the older
        ``#END,success,41`` form were read as one more data point -- a phantom
        row whose motor column held the string "success", which flipped that
        column to text and broke centre-of-mass and the derivative view in any
        reader that matched the marker by whole line (BLUETELLA's does).

        Everything after ``#END`` is ignored by those readers, so the status and
        the real point count live there safely. Named ``points_written`` rather
        than ``num_points`` on purpose: the header already carries
        ``num_points``, the count the scan was ASKED for, and an aborted scan
        must not silently overwrite it with the smaller number it achieved.
        """
        if getattr(self, "_closed", False):
            return
        self._closed = True
        self._append(
            [
                ["#END"],
                ["#exit_status", status],
                ["#points_written", self.num_points],
            ]
        )
        if self.verbose:
            print(f"Scan file closed ({status}, {self.num_points} points): {self.path}")


# =============================================================================
# Building a ScanCsv from the template
# =============================================================================


def open_scan(
    path,
    motor,
    det=None,
    scan_type="",
    command="",
    num_points=0,
    count_time=0.0,
    h5_file="",
    template=None,
    extra=None,
    verbose=True,
    shape=None,
):
    """Resolve the template against the live hardware and return a :class:`ScanCsv`.

    Everything is resolved *before* the scan starts, so a typo in the template
    raises here rather than halfway through a scan with the beam on.

    args:
        path: where to write the .csv
        motor: the scanned positioner, or a sequence of them for a two-motor
            scan (d2scan, a2scan, dmesh, mesh)
        det: detector for ``det.*`` sources, or None
        scan_type, command, num_points, count_time, h5_file:
            values the template can substitute as ``{scan_type}`` etc.
        template: path override; otherwise $ID8_SCAN_CSV_TEMPLATE, else default
        extra: more substitutions available to the template as ``{name}``
        verbose: print a line per point as the scan runs
        shape: ``(num1, num2)`` for a raster (dmesh/mesh); None for everything
            else. A viewer should test for the *presence* of the resulting
            ``shape`` header line to decide raster-vs-trajectory -- a d2scan
            sweeps both motors along one line and must NOT be gridded.
    """
    config, source = load_template(template)
    if source != DEFAULT_TEMPLATE:
        # Silent when it is the shipped default; loud when it is not, so nobody
        # wonders why their columns changed.
        print(f"NOTE: using CSV template {source}")

    # One motor or several. A bare positioner becomes a one-element tuple and
    # everything below treats it uniformly, so a dscan's .csv stays byte-identical
    # to the ones written before two-motor scans existed.
    motors = tuple(motor) if isinstance(motor, (list, tuple)) else (motor,)
    if not motors:
        raise ValueError("open_scan needs at least one motor")
    names = [m.name for m in motors]
    if len(set(names)) != len(names):
        # Two columns would render to the same label and be rejected downstream
        # with a confusing message; and a scan that moves one axis against itself
        # is a typo, not a scan.
        raise ValueError(f"the same motor is listed twice in this scan: {names}")

    context = {
        "motor": motors[0].name,
        "motors": " ".join(names),
        "det": det.name if det is not None else "",
        "scan_type": scan_type,
        "command": command,
        "num_points": num_points,
        "count_time": count_time,
        "h5_file": h5_file,
        # Bare file name as well as the full path. Both, because they answer
        # different questions: the path says where the images went, the name is
        # the stem this .csv shares with them and is what you type when looking
        # a measurement up. Empty for save_img=0, when there is no .h5 at all.
        "h5_name": os.path.basename(h5_file) if h5_file else "",
        "start_time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "epoch": f"{time.time():.3f}",
    }
    # {motor1}, {motor2}, ... one per scanned motor; blank for the ones this scan
    # does not have, so a header line or column for motor2 drops out through
    # skip_if_empty / optional rather than raising. {motor} stays an alias for
    # {motor1}.
    for index in range(1, MAX_MOTORS + 1):
        context[f"motor{index}"] = motors[index - 1].name if index <= len(motors) else ""
    # scan_8id.py's dmesh put "shape": [num1, num2] into its run_wrapper md dict.
    # This is the same number, in the one place a viewer can still read it.
    context["num1"] = shape[0] if shape else ""
    context["num2"] = shape[1] if shape else ""
    context["shape"] = f"{shape[0]}x{shape[1]}" if shape else ""
    # Default the commanded-grid keys to blank so a template line naming them
    # renders (and is then dropped by skip_if_empty) on a scan that is not a
    # raster, instead of raising "unknown substitution". A raster overwrites
    # them from extra= just below.
    for _key in ("motor1_start", "motor1_stop", "motor2_start", "motor2_stop"):
        context.setdefault(_key, "")
    context.update(extra or {})

    header_lines, skipped_lines = _build_header_lines(config, motors, det, context)
    columns, skipped_columns = _build_columns(config, motors, det, context)

    if not columns:
        raise ValueError(f"CSV template {source} defines no usable columns")

    # The template decides what goes in the file -- except that a scan must be
    # plottable against every axis it actually moved. A .csv with no column for
    # motor2 cannot be drawn as a map, and nothing downstream can reconstruct
    # where that axis was. Catch it here, not in analysis a week later.
    recorded = set()
    for column in columns:
        if column.special in ("motor", "motor_setpoint"):
            recorded.add(column.motor_index)
    missing = []
    for index, scanned_motor in enumerate(motors):
        if index not in recorded:
            missing.append(scanned_motor.name)
    if missing:
        raise ValueError(
            f"CSV template {source} records no column for scanned motor(s) {missing}. "
            'Add e.g.   - {label: "{motor2}", source: motor2, optional: true}   '
            "to its 'columns:' list."
        )
    if skipped_lines:
        print(f"NOTE: header lines skipped (not available): {', '.join(skipped_lines)}")
    if skipped_columns:
        print(f"NOTE: columns skipped (not available): {', '.join(skipped_columns)}")

    scan = ScanCsv(
        path,
        header_lines,
        str(config.get("marker", "#DATA")),
        columns,
        motors=motors,
        verbose=verbose,
    )
    return scan


def _entry_values(entry, motors, det, context):
    """One header entry -> the list of cell values after its label."""
    if "values" in entry:
        # Several PVs on one line, e.g. an ROI as min_x,min_y,size_x,size_y.
        values = []
        for dotted in entry["values"]:
            signal = resolve_signal(dotted, motors, det)
            values.append(read_value(signal))
        return values
    if "source" in entry:
        # Name the individual signal, not a parent device: a bare device's
        # .get() returns a namedtuple of every component, not a value.
        return [read_value(resolve_signal(entry["source"], motors, det))]
    return [substitute(entry.get("value", ""), context)]


def _build_header_lines(config, motors, det, context):
    """The ``label,value...`` lines above the marker."""
    lines, skipped = [], []
    for entry in config.get("lines", []) or []:
        if not isinstance(entry, dict) or "label" not in entry:
            print(f"WARNING: skipping malformed CSV header entry {entry!r}.")
            continue
        raw_label = str(entry["label"])
        try:
            # substitute() inside the try: a "{motor2}" line is exactly as absent
            # from a one-motor scan as a det.* line is from a tetramm, so
            # optional: true has to cover a failed label render too. It used to
            # sit outside, where nothing could rescue it.
            label = substitute(raw_label, context)
            values = _entry_values(entry, motors, det, context)
        except Exception as exc:
            if entry.get("optional"):
                skipped.append(raw_label)
                continue
            raise RuntimeError(f"CSV header line {raw_label!r}: {exc}") from exc
        cells = [one_line(cell(v)) for v in values]
        if entry.get("skip_if_empty") and not any(c.strip() for c in cells):
            continue
        lines.append([label] + cells)
    return lines, skipped


def _build_columns(config, motors, det, context):
    """The data columns, in template order."""
    columns, skipped = [], []
    for entry in config.get("columns", []) or []:
        if not isinstance(entry, dict) or "source" not in entry:
            print(f"WARNING: skipping malformed CSV column entry {entry!r}.")
            continue
        source = str(entry["source"])
        raw_label = str(entry.get("label", source))
        try:
            # See the note in _build_header_lines: substitute() belongs inside
            # the try so `optional: true` covers a label that cannot render.
            label = substitute(raw_label, context)
            if source in SPECIAL_SOURCES:
                special, index = _split_motor_index(source)
                if special in ("motor", "motor_setpoint") and index >= len(motors):
                    raise KeyError(f"source {source!r} needs {index + 1} motors, this scan has {len(motors)}")
                columns.append(Column(label, special=special, motor_index=index))
                continue
            signal = resolve_signal(source, motors, det)
            probe = signal.get()  # fail now, not at point 1 with the beam on
        except Exception as exc:
            if entry.get("optional"):
                skipped.append(raw_label)
                continue
            raise RuntimeError(f"CSV column {raw_label!r}: {exc}") from exc
        if not isinstance(probe, numbers.Real):
            print(
                f"WARNING: CSV column {label!r} reads {type(probe).__name__}, "
                "not a single number; it will be written as text."
            )
        columns.append(Column(label, signal=signal))
    seen, repeated = set(), []
    for column in columns:
        if column.label in seen:
            repeated.append(column.label)
        seen.add(column.label)
    if repeated:
        raise ValueError(
            f"CSV template has duplicate column labels: {sorted(set(repeated))}. "
            "Every column needs its own name -- most often this means two "
            "entries render to the same text once {det} or {motor} is filled "
            "in. Rename one of them in the template's 'columns:' list."
        )
    return columns, skipped


# =============================================================================
# Reading one back (no dependencies beyond the standard library)
# =============================================================================


def read_scan_csv(path):
    """Read a scan file. Returns ``(header, labels, columns, status)``.

    ``header`` is ``{label: [values]}`` for the lines above the marker,
    ``columns`` is ``{label: [floats]}`` for the data. Missing/partial trailing
    rows are ignored, so this is safe to call on a scan that is still running.

    The leading point-index column (see :data:`INDEX_LABEL`) is dropped: it is
    just the row number, and leaving it out keeps ``labels``/``columns`` the
    same shape they had before that column existed. Files written before it was
    added have no such column and read exactly as they always did.
    """
    with open(path, "r", newline="") as handle:
        rows = list(csv.reader(handle))
    marker = "#DATA"

    header = {}
    labels = []
    columns = {}
    status = "running"  # stays "running" when the file has no #END line yet
    seen_marker = False
    seen_end = False
    has_index = False
    for row in rows:
        if not row:
            continue
        key = row[0].strip()
        if key == "#END":
            # Bare marker: the table ends here, but the outcome lines follow, so
            # keep reading rather than break.
            seen_end = True
            status = "unknown"
            continue
        if seen_end:
            # Trailing "#exit_status,..." / "#points_written,..." and anything
            # else a future writer adds after the table.
            if key == "#exit_status":
                status = row[1] if len(row) > 1 else "unknown"
            else:
                header[key.lstrip("#")] = row[1:]
            continue
        if not seen_marker:
            # The marker is the exact "#DATA" line. Testing for a leading "#"
            # would stop at the FIRST header line now that they all carry one.
            if key == marker:
                seen_marker = True
            else:
                # "#label,value" (this writer) or bare "label,value" (files
                # written before 2026-09-07). Strip the "#" so both give the
                # same keys.
                header[key.lstrip("#")] = row[1:]
            continue
        if not labels:
            has_index = row[0] == INDEX_LABEL
            labels = row[1:] if has_index else row
            columns = {label: [] for label in labels}
            continue
        if has_index:
            row = row[1:]
        if len(row) != len(labels):
            continue  # a half-written trailing row; it will be complete next poll
        for index, label in enumerate(labels):
            text = row[index]
            try:
                columns[label].append(float(text))
            except ValueError:
                columns[label].append(text)
    return header, labels, columns, status
