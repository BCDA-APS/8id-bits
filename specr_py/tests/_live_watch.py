"""Drive the real SpecrWindow headlessly against a growing SPEC file.

Proves the GUI's own monitor path updates live: each tick calls the same
_on_monitor_tick() the Qt timer calls, and we record how many points are
actually on the plotted line.

    QT_QPA_PLATFORM=offscreen python _live_watch.py <specfile> <seconds>
"""
import os
import sys
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from PyQt5 import QtWidgets  # noqa: E402
import specr  # noqa: E402

path = sys.argv[1]
duration = float(sys.argv[2]) if len(sys.argv) > 2 else 180.0

deadline = time.time() + duration
while not os.path.exists(path) and time.time() < deadline:
    time.sleep(0.5)
if not os.path.exists(path):
    print("FAIL: spec file never appeared")
    sys.exit(1)

app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
win = specr.SpecrWindow(path=path)
win.erase_action.setChecked(True)   # newest scan only, as during a live scan
print(f"GUI opened {path}: {len(win.specfile)} scan(s), "
      f"showing scan {win._primary_scan().number}", flush=True)

# Only track a scan that appears AFTER we start watching -- the file may
# already end with a completed scan, which tells us nothing about liveness.
baseline = win._primary_scan().number
print(f"baseline scan at startup: {baseline}; waiting for a newer one...", flush=True)

samples, last_n, tracked = [], None, None
t0 = time.time()
while time.time() < deadline:
    win._on_monitor_tick()
    app.processEvents()
    scan = win._primary_scan()
    if tracked is None:
        if scan.number > baseline:
            tracked = scan.number
            print(f"  new scan {tracked} appeared at t={time.time() - t0:.1f}s", flush=True)
        else:
            time.sleep(0.5)
            continue
    lines = win.axes.get_lines()
    n = len(lines[0].get_xdata()) if lines else 0
    if n != last_n:
        elapsed = time.time() - t0
        try:
            with open(path, "rb") as fh:
                fh.seek(0, 2)
                size = fh.tell()
        except OSError:
            size = -1
        print(f"  [{time.strftime('%H:%M:%S')}] t={elapsed:6.1f}s  scan {scan.number}: "
              f"{n} points plotted  (file {size} B)", flush=True)
        samples.append((round(elapsed, 1), scan.number, n))
        last_n = n
    if scan.number == tracked and scan.exit_status and n > 0:
        print(f"  scan {scan.number} closed: exit_status = {scan.exit_status}", flush=True)
        break
    time.sleep(1.0)

png = os.path.join(HERE, "_live_watch_result.png")
win.figure.savefig(png, dpi=90)
print(f"\ndistinct GUI updates observed: {len(samples)}")
print(f"point counts seen: {[n for _, _, n in samples]}")
print(f"final render: {png}")
grew = len(samples) >= 3 and samples[-1][2] > samples[0][2]
print("LIVE UPDATE CONFIRMED" if grew else "LIVE UPDATE NOT OBSERVED")
sys.exit(0 if grew else 1)
