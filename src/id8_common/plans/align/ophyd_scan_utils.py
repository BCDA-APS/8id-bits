"""
Small helpers the Ophyd scans use: safe cleanup, the result object, the live table.

Plain Python only -- this module imports no ophyd, no apsbits, no EPICS. Nothing
here knows what a motor is. It lives apart from ``ophyd_scan.py`` so that the
scan file contains the hardware and the scan itself, and nothing else.

Deliberately written with basic Python: plain functions, one simple class, no
decorators and no context managers.
"""

import signal
import threading

# =============================================================================
# Turning Ctrl+C off while a scan cleans up
# =============================================================================
#
# Cleanup must not be abandoned halfway. Leaving the detector acquiring, or the
# motor parked in the middle of a scan, is worse than making the user wait a
# couple of seconds. So Ctrl+C is switched off while the cleanup runs and
# switched back on afterwards.
#
# Use them in a pair, like this:
#
#     saved = disable_ctrl_c()
#     try:
#         ...cleanup steps...
#     finally:
#         restore_ctrl_c(saved)
#
# The try/finally matters: it guarantees Ctrl+C is switched back on even if a
# cleanup step raises.


def _ignore_ctrl_c(signal_number, frame):
    """Called instead of raising KeyboardInterrupt while cleanup is running."""
    print("    (cleanup in progress -- interrupt ignored)", flush=True)


def disable_ctrl_c():
    """Switch Ctrl+C off. Pass the return value to :func:`restore_ctrl_c`.

    Returns ``(changed, previous_handler)``. ``changed`` is False when we could
    not switch it off -- Python only allows this on the main thread -- and in
    that case the scan simply behaves as it always did.
    """
    if threading.current_thread() is not threading.main_thread():
        return (False, None)
    try:
        previous_handler = signal.signal(signal.SIGINT, _ignore_ctrl_c)
    except (ValueError, OSError):
        return (False, None)
    return (True, previous_handler)


def restore_ctrl_c(saved):
    """Switch Ctrl+C back on. ``saved`` is what :func:`disable_ctrl_c` returned."""
    changed, previous_handler = saved
    if not changed:
        return
    try:
        signal.signal(signal.SIGINT, previous_handler)
    except (ValueError, OSError):
        pass


def safe_call(description, function, *args, **kwargs):
    """Run one cleanup step. If it fails, print a warning and carry on.

    A dead IOC partway through cleanup must not stop us parking the motor or
    closing the SPEC block, so every teardown step goes through here::

        safe_call("block the beam", blockbeam)
        safe_call("stop the detector", det.cam.acquire.put, 0)
        safe_call("park the motor", motor.move, 30.0, wait=True)

    ``description`` comes first so the call site reads as a sentence, and it is
    what gets printed if the step fails.
    """
    try:
        function(*args, **kwargs)
    except (Exception, KeyboardInterrupt) as exc:
        print(f"WARNING: cleanup step '{description}' failed: {exc}")


# =============================================================================
# What a scan hands back
# =============================================================================


class ScanResult:
    """The value a scan returns.

    Holds the measured motor positions and the counter columns. It prints as a
    single summary line, so an interactive ``dscan_ophyd(...)`` does not dump
    every number as the cell's ``Out[N]`` echo -- the per-point table has
    already shown them.

    Both of these work::

        result = dscan_ophyd(...)
        result.positions            # array of measured positions
        result.columns              # {column label: [values]}

        positions, columns = dscan_ophyd(...)      # unpacking still works
    """

    def __init__(self, positions, columns, scan_num=None, motor_name="motor"):
        self.positions = positions
        self.columns = columns
        self.scan_num = scan_num
        self.motor_name = motor_name

    def __iter__(self):
        """Allows ``positions, columns = dscan_ophyd(...)``."""
        return iter((self.positions, self.columns))

    def __repr__(self):
        count = len(self.positions)
        if count:
            span = f"{self.positions[0]:.5g} -> {self.positions[-1]:.5g}"
        else:
            span = "no points"
        peaks = []
        for label, values in self.columns.items():
            if values:
                peaks.append(f"{label.split('_', 1)[-1]} max={max(values):g}")
            else:
                peaks.append(f"{label} empty")
        return (f"<ScanResult #S {self.scan_num}: {count} pts, "
                f"{self.motor_name} {span} | {', '.join(peaks)}>")


# =============================================================================
# The table printed while a scan runs
# =============================================================================


def table_header(motor_name, labels):
    """Column header for the per-point live table."""
    columns = [("#", 4), (motor_name, max(len(motor_name), 13)), ("time[s]", 8)]
    for label in labels:
        columns.append((label, max(len(label), 12)))
    heading = "  ".join(name.rjust(width) for name, width in columns)
    underline = "  ".join("-" * width for _, width in columns)
    return heading + "\n" + underline


def table_row(index, motor_name, position, elapsed, labels, values):
    """One line of the per-point live table."""
    cells = [
        str(index).rjust(4),
        f"{position:.6f}".rjust(max(len(motor_name), 13)),
        f"{elapsed:.1f}".rjust(8),
    ]
    for label, value in zip(labels, values, strict=False):
        try:
            text = f"{float(value):g}"
        except (TypeError, ValueError):
            text = str(value)
        cells.append(text.rjust(max(len(label), 12)))
    return "  ".join(cells)
