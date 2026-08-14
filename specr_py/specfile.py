"""
Incremental SPEC data-file reader.

No GUI dependencies -- importable and testable headless.

Designed for *live* reading: the file is opened, the bytes appended since the
last call are consumed, and the file is closed again. Only complete lines (up
to the last newline) are parsed, so a reader polling a file that is being
appended to never sees a torn row.

Handles the SPEC dialect that apstools / Bluesky and ``ophyd_scan.py`` emit:

    #F <path>                       file header (may repeat per session)
    #E <epoch>
    #D <date>
    #C Bluesky  user = ...  host = ...
    #O0 mot1  mot2  ...             motor names, 8 per line, TWO-space separated
    #o0 mot1  mot2  ...             mnemonics
                                    (blank line)
    #S <n>  <command>               start of a scan
    #D <date>
    #C <date>.  key = value         comments
    #MD key = value                 metadata
    #P0 <pos> <pos> ...             motor positions matching #O
    #N <ncols>
    #L lab1  lab2  ...              column labels, TWO-space separated
    <numeric rows>
    #C <date>.  exit_status = ...

Label splitting follows the same rule the MATLAB ``specr`` reader and
``spec2nexus`` use: split on runs of two or more spaces, because a label may
itself contain a single space (e.g. ``SMS granite``).
"""

import os
import re

import numpy as np

_MULTISPACE = re.compile(r"\s{2,}")
_BOOLEANS = {"true": 1.0, "false": 0.0, "on": 1.0, "off": 0.0}


def _to_float(token):
    """Data cells are usually numeric, but Bluesky writes bare True/False too."""
    try:
        return float(token)
    except ValueError:
        return _BOOLEANS.get(token.strip().lower(), float("nan"))


def split_labels(text, ncols=None):
    """Split a ``#L``/``#O`` payload on two-or-more spaces.

    Falls back to single-space splitting only when that yields exactly the
    column count declared by ``#N`` -- the same guard spec2nexus applies.
    """
    parts = [p for p in _MULTISPACE.split(text.strip()) if p]
    if ncols and len(parts) != ncols:
        single = text.split()
        if len(single) == ncols:
            return single
    return parts


class SpecScan:
    """One ``#S`` block."""

    def __init__(self, index, number, command, motor_names):
        self.index = index
        self.number = number
        self.command = command
        self.motor_names = list(motor_names)
        self.motor_positions = []
        self.date = ""
        self.labels = []
        self.ncols = None
        self.rows = []
        self.comments = []
        self.meta = {}
        self.exit_status = None

    # -- accessors ---------------------------------------------------------

    @property
    def npoints(self):
        return len(self.rows)

    @property
    def header(self):
        return f"#S {self.number}  {self.command}"

    def column(self, label):
        """Values for one column label as a float array."""
        if label not in self.labels:
            return np.empty(0)
        i = self.labels.index(label)
        return np.array([r[i] for r in self.rows if len(r) > i], dtype=float)

    def as_dict(self):
        return {lab: self.column(lab) for lab in self.labels}

    def motor_table(self):
        """[(name, position)] pairs, tolerant of a length mismatch."""
        n = min(len(self.motor_names), len(self.motor_positions))
        return list(zip(self.motor_names[:n], self.motor_positions[:n]))

    def __repr__(self):
        return f"<SpecScan #{self.number} {self.npoints}pts {self.command[:40]!r}>"


class SpecDataFile:
    """Incrementally-parsed SPEC file.

    ``refresh()`` consumes only the newly appended bytes and returns True when
    something changed, so it is cheap to call on a timer.
    """

    def __init__(self, path):
        self.path = str(path)
        self.scans = []
        self.file_headers = []
        self._offset = 0
        self._motor_names = []
        self._motor_mnemonics = []
        self._current = None

    # -- internals ---------------------------------------------------------

    def _reset(self):
        self.scans = []
        self.file_headers = []
        self._offset = 0
        self._motor_names = []
        self._motor_mnemonics = []
        self._current = None

    def _numbered_payload(self, line, tag):
        """``#O2 a  b`` -> (2, 'a  b'); returns (None, '') if it isn't ours."""
        body = line[len(tag):]
        digits = ""
        while body and body[0].isdigit():
            digits += body[0]
            body = body[1:]
        if not digits:
            return None, ""
        return int(digits), body.strip()

    def _feed(self, line):
        line = line.rstrip("\r")

        if line.startswith("#S "):
            payload = line[3:].strip()
            bits = payload.split(None, 1)
            try:
                number = int(bits[0])
            except (ValueError, IndexError):
                return
            command = bits[1].strip() if len(bits) > 1 else ""
            self._current = SpecScan(
                len(self.scans), number, command, self._motor_names
            )
            self.scans.append(self._current)
            return

        if line.startswith("#F"):
            self.file_headers.append(line[2:].strip())
            self._motor_names = []
            self._motor_mnemonics = []
            return

        if line.startswith("#O"):
            idx, payload = self._numbered_payload(line, "#O")
            if idx is not None:
                if idx == 0:
                    self._motor_names = []
                self._motor_names.extend(split_labels(payload))
            return

        if line.startswith("#o"):
            idx, payload = self._numbered_payload(line, "#o")
            if idx is not None:
                if idx == 0:
                    self._motor_mnemonics = []
                self._motor_mnemonics.extend(split_labels(payload))
            return

        scan = self._current
        if scan is None:
            return

        if line.startswith("#D "):
            if not scan.date:
                scan.date = line[3:].strip()
            return

        if line.startswith("#P"):
            idx, payload = self._numbered_payload(line, "#P")
            if idx is not None:
                if idx == 0:
                    scan.motor_positions = []
                scan.motor_positions.extend(_to_float(t) for t in payload.split())
            return

        if line.startswith("#N "):
            try:
                scan.ncols = int(line[3:].split()[0])
            except (ValueError, IndexError):
                pass
            return

        if line.startswith("#L "):
            scan.labels = split_labels(line[3:], scan.ncols)
            return

        if line.startswith("#MD "):
            key, _, value = line[4:].partition("=")
            scan.meta[key.strip()] = value.strip()
            return

        if line.startswith("#C"):
            text = line[2:].strip()
            scan.comments.append(text)
            if "exit_status" in text:
                scan.exit_status = text.split("=", 1)[-1].strip()
            return

        if line.startswith("#"):
            return

        if not line.strip():
            return

        if scan.labels:
            row = [_to_float(t) for t in line.split()]
            if len(row) == len(scan.labels):
                scan.rows.append(row)
            elif len(row) > 1:
                # tolerate a ragged row rather than dropping the scan
                row = (row + [float("nan")] * len(scan.labels))[: len(scan.labels)]
                scan.rows.append(row)

    # -- public API --------------------------------------------------------

    def refresh(self):
        """Consume newly appended complete lines. True if anything changed.

        Two details matter for reading a file that another host is appending to
        over NFS:

        * The size is taken with ``seek(0, 2)`` on a freshly opened handle, not
          with ``os.stat``. NFS close-to-open consistency revalidates on open,
          whereas a cached ``stat`` can report a stale size.
        * ``_offset`` advances **only** past the last complete line. A partial
          trailing line is simply left unread and picked up next time. An
          earlier version also kept the partial text in a ``_carry`` buffer
          *and* rewound the offset over it, so the fragment was counted twice
          and rows got concatenated ("dedef" from "de" + "def"). There is no
          carry now -- the file itself is the only buffer.
        """
        try:
            with open(self.path, "rb") as f:
                f.seek(0, os.SEEK_END)
                size = f.tell()
                if size < self._offset:  # truncated or replaced -> start over
                    self._reset()
                if size == self._offset:
                    return False
                f.seek(self._offset)
                raw = f.read()
        except OSError:
            return False
        if not raw:
            return False

        cut = raw.rfind(b"\n")
        if cut == -1:
            return False  # nothing complete yet; leave the offset alone
        self._offset += cut + 1

        for line in raw[:cut].decode("utf-8", errors="replace").split("\n"):
            self._feed(line)
        return True

    def scan_headers(self):
        return [s.header for s in self.scans]

    def last(self):
        return self.scans[-1] if self.scans else None

    def by_number(self, number):
        for s in self.scans:
            if s.number == number:
                return s
        return None

    def __len__(self):
        return len(self.scans)

    def __getitem__(self, i):
        return self.scans[i]


# =============================================================================
# Scan statistics -- ports params.m from the MATLAB specr
# =============================================================================


def scan_stats(x, y):
    """Return ``(peak_x, peak_y, com, fwhm, fwhm_center)``.

    Mirrors ``params.m``: peak by maximum, centre of mass by intensity
    weighting, FWHM by interpolating the half-maximum crossings either side of
    the peak, with the same fallbacks when a crossing is not bracketed.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if x.size == 0:
        return (np.nan,) * 5

    ipeak = int(np.argmax(y))
    peak_x, peak_y = x[ipeak], y[ipeak]

    total = float(np.sum(y))
    com = float(np.sum(x * y) / total) if total != 0 else 0.0

    half = peak_y / 2.0
    left = right = None
    try:
        below = np.where(y[: ipeak + 1] < half)[0]
        if below.size:
            seg = slice(below[-1], ipeak + 1)
            xs, ys = x[seg], y[seg]
            order = np.argsort(ys)
            left = float(np.interp(half, ys[order], xs[order]))
    except Exception:
        left = None
    try:
        below = np.where(y[ipeak:] < half)[0]
        if below.size:
            seg = slice(ipeak, ipeak + below[0] + 1)
            xs, ys = x[seg], y[seg]
            order = np.argsort(ys)
            right = float(np.interp(half, ys[order], xs[order]))
    except Exception:
        right = None

    if left is not None and right is not None:
        center, fwhm = (left + right) / 2.0, right - left
    elif left is not None:
        center, fwhm = left, float(x[-1]) - left
    elif right is not None:
        center, fwhm = right, right - float(x[0])
    else:
        center = (float(x[-1]) + float(x[0])) / 2.0
        fwhm = float(x[-1]) - float(x[0])

    return peak_x, peak_y, com, fwhm, center
