"""
SPEC-format file writer.

**Standard library only** -- no ophyd, no apsbits, no bluesky. That is deliberate:
this module is the piece every Ophyd scan shares, so keeping it free of hardware
imports means it can be exercised offline and reused by any future scan
(``d2scan``, ``mesh``, ...) without dragging in the beamline stack.

The format rules enforced here are not stylistic. They come from the
``spec2nexus`` parser and from ``specfile.py`` in the beamline's ``specr_py``
viewer. Getting any of them wrong produces a file that *looks* correct and
silently mis-parses:

* ``#L`` labels are separated by **two** spaces. Both readers split on a
  ``\\s{2,}`` pattern, which is what allows a label to contain a single space.
  A label that itself contains two spaces therefore becomes two columns --
  :func:`validate_labels` rejects that up front.
* ``#N`` must equal ``len(labels)``.
* ``#S`` needs the trailing space -- scan detection is ``line.startswith("#S ")``.
* A blank line must precede each ``#S``.
* Every data cell must parse as a float.

Every write opens the file, appends complete newline-terminated lines, and closes
it. A single append-mode write of a short line is atomic for a regular file, so a
polling reader never observes a torn row.
"""

import datetime
import getpass
import math
import os
import re
import socket
import time

SPEC_TIME_FORMAT = "%a %b %d %H:%M:%S %Y"

#: Both readers split ``#L``/``#O`` payloads on runs of two or more spaces.
_MULTISPACE = re.compile(r"\s{2,}")

#: Columns already reported as non-scalar, so the warning prints once per label
#: rather than once per scan point.
_NONSCALAR_WARNED = set()


def _fmt(value, label=None):
    """Format one data cell so that ``float()`` will always parse it.

    Readers do ``float(cell)`` on every token, so a non-finite value or a string
    would break the whole scan block. Substitute 0 rather than emitting a token
    some parsers reject.

    A non-scalar (an array-valued signal, say) would format to a token containing
    spaces and shift every later column, so it is reported by ``label`` -- once --
    and written as 0. Without the warning the corruption is invisible.
    """
    if isinstance(value, bool):
        return "1" if value else "0"
    try:
        v = float(value)
    except (TypeError, ValueError):
        if label and label not in _NONSCALAR_WARNED:
            _NONSCALAR_WARNED.add(label)
            print(
                f"WARNING: SPEC column {label!r} is not a scalar "
                f"({type(value).__name__}); writing 0. Choose a scalar signal in "
                "the SPEC template."
            )
        return "0"
    if not math.isfinite(v):
        return "0"
    return f"{v:.12g}"


def _now():
    return datetime.datetime.now().strftime(SPEC_TIME_FORMAT)


def _one_line(text):
    """Flatten text so it cannot break the one-record-per-line structure.

    Free-text header values (a user's scan comment, most obviously) may contain
    anything. A newline in the middle of a ``#MD`` line would end that record
    early and leave the remainder looking like a data row -- which a reader would
    either reject or, worse, silently append to the scan. Collapsing all
    whitespace runs to single spaces removes that class of corruption; tabs and
    carriage returns go the same way.
    """
    return " ".join(str(text).split())


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


def _numbered_payload(line, tag):
    """``"#O2 a  b"`` -> ``(2, "a  b")``; ``(None, "")`` when it is not ours.

    Mirrors ``SpecDataFile._numbered_payload`` in the viewer so that what we
    write back is read the same way.
    """
    body = line[len(tag):]
    digits = ""
    while body and body[0].isdigit():
        digits += body[0]
        body = body[1:]
    if not digits:
        return None, ""
    return int(digits), body.strip()


def _split_labels(text):
    """Split an ``#L``/``#O`` payload on two-or-more spaces."""
    return [p for p in _MULTISPACE.split(text.strip()) if p]


def validate_labels(labels):
    """Raise ``ValueError`` if these column labels cannot survive a round trip.

    Called before a scan touches hardware, so a bad template fails immediately
    instead of producing a file that mis-parses hours later.
    """
    problems = []
    for label in labels:
        if not str(label).strip():
            problems.append("empty label")
        elif _MULTISPACE.search(str(label)):
            problems.append(
                f"{label!r} contains two or more consecutive spaces, which "
                "readers would split into separate columns"
            )
    seen, dupes = set(), []
    for label in labels:
        if label in seen and label not in dupes:
            dupes.append(label)
        seen.add(label)
    for label in dupes:
        problems.append(
            f"{label!r} appears more than once; readers index columns by label "
            "and would always return the first"
        )
    if problems:
        raise ValueError(
            "invalid SPEC column labels -- fix the SPEC template:\n  - "
            + "\n  - ".join(problems)
        )


class SpecFile:
    """Minimal SPEC data file writer. One instance per scan is fine."""

    def __init__(self, path):
        self.path = str(path)
        self.labels = []
        self._points = 0

    # -- internals ---------------------------------------------------------

    def _append(self, lines):
        """Append complete lines. Never raises into the caller.

        A failure here must not kill a scan that is otherwise fine -- the images
        are the primary data; this file is a convenience.
        """
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

    def current_motor_names(self):
        """Motor names in effect from the file's last ``#O`` block.

        Returns ``None`` when the file is absent or unreadable. Used to detect a
        stale ``#O`` -- see :meth:`write_file_header`.
        """
        if not self.exists():
            return None
        names = []
        try:
            with open(self.path, "r", errors="replace") as f:
                for line in f:
                    line = line.rstrip("\n").rstrip("\r")
                    if line.startswith("#F"):
                        names = []
                    elif line.startswith("#O"):
                        idx, payload = _numbered_payload(line, "#O")
                        if idx is None:
                            continue
                        if idx == 0:
                            names = []
                        names.extend(_split_labels(payload))
        except OSError:
            return None
        return names

    # -- public API --------------------------------------------------------

    def write_file_header(self, motor_names=(), force=False):
        """Write the ``#F``/``#E``/``#D``/``#C``/``#O``/``#o`` block if needed.

        Normally skipped when the file already has content, so one file
        accumulates scans under a single header.

        The exception matters. ``#O`` (motor *names*) is written once per file,
        but ``#P`` (motor *positions*) once per scan, and a reader binds the
        names in force at each ``#S``. So if the template's positioner list
        changes mid-experiment, later scans would pair new positions with old
        names -- wrong labels, and no error anywhere. When that mismatch is
        detected we emit a fresh ``#F`` block; readers reset their motor-name
        table on ``#F``, which rebinds the names correctly from that point on.
        """
        motor_names = [str(n) for n in motor_names]
        reopening = False

        if self.exists() and not force:
            existing = self.current_motor_names()
            if existing == motor_names:
                return
            if existing is not None and not motor_names:
                # Nothing to declare; leave whatever is already in force.
                return
            reopening = True
            print(
                "NOTE: SPEC positioner list changed; starting a new #F block so "
                "#O names stay matched to #P values."
            )

        now = time.time()
        user = getpass.getuser() or "8idUser"
        host = socket.gethostname() or "localhost"
        lines = [""] if reopening else []
        lines += [
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
        """Open a new ``#S`` block. ``labels`` defines the columns.

        Raises ``ValueError`` for labels that cannot round-trip. This runs before
        the caller moves anything, so a bad template stops the scan at the start.
        """
        labels = [str(lab) for lab in labels]
        validate_labels(labels)

        self.labels = labels
        self._points = 0
        stamp = _now()
        lines = [
            "",  # blank line before #S is required
            f"#S {int(scan_num)}  {_one_line(command)}",
            f"#D {stamp}",
        ]
        lines += [f"#C {stamp}.  {_one_line(c)}" for c in comments]
        # Metadata order is the caller's; the scan comment is put first so it is
        # the first thing you see (and the first thing Ctrl+F lands on).
        for key, value in (metadata or {}).items():
            lines.append(f"#MD {_one_line(key)} = {_one_line(value)}")
        lines += _numbered_lines("#P", [_fmt(v) for v in motor_positions]) or ["#P0 "]
        lines += [
            "#N " + str(len(self.labels)),
            "#L " + "  ".join(self.labels),  # TWO spaces -- see module docstring
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
        self._append([" ".join(_fmt(v, lab) for lab, v in zip(self.labels, vals, strict=False))])
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
