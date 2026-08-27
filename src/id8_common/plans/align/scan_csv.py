"""
Write a scan to a plain, flat CSV file -- one file per scan.

This module knows about *file format only*. It has no scan logic and moves no
hardware, so the motor-control code in :mod:`ophyd_scan` stays readable.

What the file looks like::

    start_time,2026-08-26 14:32:07
    h5_file,/gdata/dm/8ID/8IDE/2026-2/pope202607/data/bluesky/A0201_Test_a0010.h5
    comment,3x3 grid spot 5 after realigning KB
    ... more one-per-line entries ...
    #DATA
    huber_delta,elapsed_time,lambda2M_stats1_total
    31.0,0.0,1234
    31.025,1.05,1301
    ...
    #END,success,41

* Everything **above** the ``#DATA`` marker is a ``label,value`` line: things
  that do not change during the scan (start time, the ``.h5`` file name, the
  comment, motor positions at the start, ROI settings, ...).
* The line **after** the marker holds the column names.
* Everything after that is one row per scan point.
* The last line starts with ``#`` and says how the scan ended.

Reading one back is three lines of Python::

    import pandas
    rows = pandas.read_csv(path, skiprows=lambda i: i <= marker_line, comment="#")

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
from apsbits.core.instrument_init import oregistry

#: Shipped default. The repo is an editable install, so users edit it in place.
DEFAULT_TEMPLATE = Path(__file__).resolve().parents[2] / "configs" / "scan_csv_template.yml"

#: Column sources the scan computes itself instead of reading from a PV.
SPECIAL_SOURCES = ("motor", "motor_setpoint", "elapsed", "epoch")

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


def resolve_signal(dotted, motor, det):
    """ "det.stats1.total" or "tetramm1.current1.mean_value" -> the ophyd object.

    ``det.`` and ``motor.`` mean the objects passed to *this* scan, which is what
    lets one template serve eiger4M, lambda2M and tetramm1.
    """
    head, _, rest = str(dotted).partition(".")
    if head == "det":
        if det is None:
            raise KeyError(f"{dotted!r} needs a detector but this scan has none")
        obj = det
    elif head == "motor":
        obj = motor
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
    """Apply the {motor} / {det} / {comment} ... replacements, with a clear error."""
    try:
        return str(text).format(**context)
    except KeyError as exc:
        raise KeyError(f"unknown substitution {exc} in {text!r}; " f"available: {', '.join(sorted(context))}") from exc


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
    """Flatten newlines so a pasted comment cannot break the file structure."""
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

    def __init__(self, label, special=None, signal=None):
        """Store the column name and where its value comes from."""
        self.label = label
        self.special = special
        self.signal = signal

    def value(self, position, setpoint, elapsed):
        """Return this column's value for the point being written."""
        if self.special == "motor":
            return position
        if self.special == "motor_setpoint":
            return setpoint
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

    def __init__(self, path, header_lines, marker, columns, verbose=True):
        """Store the layout. Nothing is written until write_header() is called."""
        self.path = str(path)
        self.header_lines = header_lines  # list of [label, v1, v2, ...]
        self.marker = marker
        self.columns = columns  # list of Column
        self.labels = [c.label for c in columns]
        self.verbose = verbose
        self.t0 = time.time()
        self.num_points = 0
        #: every value read so far, {label: [...]} -- handy at the prompt
        self.data = {label: [] for label in self.labels}

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
        """Header lines, the marker, and the column-name line."""
        rows = list(self.header_lines)
        rows.append([self.marker])
        rows.append(self.labels)
        self._append(rows)
        if self.verbose:
            print(f"Scan file: {self.path}")
            print("   # " + " ".join(f"{short_label(x):>{SCREEN_WIDTH}s}" for x in self.labels), flush=True)

    def add_point(self, setpoint):
        """Read every column once and append a row. Returns the values."""
        elapsed = time.time() - self.t0
        position = read_value(self._motor) if self._motor is not None else setpoint
        values = [c.value(position, setpoint, elapsed) for c in self.columns]
        self._append([[cell(v) for v in values]])
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
        """Write the closing line. Safe to call twice."""
        if getattr(self, "_closed", False):
            return
        self._closed = True
        self._append([["#END", status, self.num_points]])
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
    comment="",
    template=None,
    extra=None,
    verbose=True,
):
    """Resolve the template against the live hardware and return a :class:`ScanCsv`.

    Everything is resolved *before* the scan starts, so a typo in the template
    raises here rather than halfway through a scan with the beam on.

    args:
        path: where to write the .csv
        motor: the scanned positioner
        det: detector for ``det.*`` sources, or None
        scan_type, command, num_points, count_time, h5_file, comment:
            values the template can substitute as ``{scan_type}`` etc.
        template: path override; otherwise $ID8_SCAN_CSV_TEMPLATE, else default
        extra: more substitutions (a future d2scan passes ``motor2`` here)
        verbose: print a line per point as the scan runs
    """
    config, source = load_template(template)
    if source != DEFAULT_TEMPLATE:
        # Silent when it is the shipped default; loud when it is not, so nobody
        # wonders why their columns changed.
        print(f"NOTE: using CSV template {source}")

    context = {
        "motor": motor.name,
        "det": det.name if det is not None else "",
        "scan_type": scan_type,
        "command": command,
        "num_points": num_points,
        "count_time": count_time,
        "h5_file": h5_file,
        "comment": one_line("" if comment is None else comment),
        "start_time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "epoch": f"{time.time():.3f}",
    }
    context.update(extra or {})

    header_lines, skipped_lines = _build_header_lines(config, motor, det, context)
    columns, skipped_columns = _build_columns(config, motor, det, context)

    if not columns:
        raise ValueError(f"CSV template {source} defines no usable columns")
    if skipped_lines:
        print(f"NOTE: header lines skipped (not available): {', '.join(skipped_lines)}")
    if skipped_columns:
        print(f"NOTE: columns skipped (not available): {', '.join(skipped_columns)}")

    scan = ScanCsv(path, header_lines, str(config.get("marker", "#DATA")), columns, verbose=verbose)
    scan._motor = motor
    return scan


def _entry_values(entry, motor, det, context):
    """One header entry -> the list of cell values after its label."""
    if "values" in entry:
        # Several PVs on one line, e.g. an ROI as min_x,min_y,size_x,size_y.
        return [read_value(resolve_signal(p, motor, det)) for p in entry["values"]]
    if "source" in entry:
        # Name the individual signal, not a parent device: a bare device's
        # .get() returns a namedtuple of every component, not a value.
        return [read_value(resolve_signal(entry["source"], motor, det))]
    return [substitute(entry.get("value", ""), context)]


def _build_header_lines(config, motor, det, context):
    """The ``label,value...`` lines above the marker."""
    lines, skipped = [], []
    for entry in config.get("lines", []) or []:
        if not isinstance(entry, dict) or "label" not in entry:
            print(f"WARNING: skipping malformed CSV header entry {entry!r}.")
            continue
        label = substitute(entry["label"], context)
        try:
            values = _entry_values(entry, motor, det, context)
        except Exception as exc:
            if entry.get("optional"):
                skipped.append(label)
                continue
            raise RuntimeError(f"CSV header line {label!r}: {exc}") from exc
        cells = [one_line(cell(v)) for v in values]
        if entry.get("skip_if_empty") and not any(c.strip() for c in cells):
            continue
        lines.append([label] + cells)
    return lines, skipped


def _build_columns(config, motor, det, context):
    """The data columns, in template order."""
    columns, skipped = [], []
    for entry in config.get("columns", []) or []:
        if not isinstance(entry, dict) or "source" not in entry:
            print(f"WARNING: skipping malformed CSV column entry {entry!r}.")
            continue
        source = str(entry["source"])
        label = substitute(entry.get("label", source), context)
        if source in SPECIAL_SOURCES:
            columns.append(Column(label, special=source))
            continue
        try:
            signal = resolve_signal(source, motor, det)
            probe = signal.get()  # fail now, not at point 1 with the beam on
        except Exception as exc:
            if entry.get("optional"):
                skipped.append(label)
                continue
            raise RuntimeError(f"CSV column {label!r}: {exc}") from exc
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
    """
    with open(path, "r", newline="") as handle:
        rows = list(csv.reader(handle))

    header, labels, columns, status = {}, [], {}, "running"
    seen_marker = False
    for row in rows:
        if not row:
            continue
        if row[0].startswith("#END"):
            status = row[1] if len(row) > 1 else "unknown"
            break
        if not seen_marker:
            if row[0].startswith("#"):
                seen_marker = True
            else:
                header[row[0]] = row[1:]
            continue
        if not labels:
            labels = row
            columns = {label: [] for label in labels}
            continue
        if len(row) != len(labels):
            continue  # a half-written trailing row; it will be complete next poll
        for index, label in enumerate(labels):
            text = row[index]
            try:
                columns[label].append(float(text))
            except ValueError:
                columns[label].append(text)
    return header, labels, columns, status
