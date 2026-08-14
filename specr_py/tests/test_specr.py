"""Offline tests: SPEC parser, ophyd_scan writer round-trip, headless GUI live update.

No EPICS, no hardware, no display. Run with:

    QT_QPA_PLATFORM=offscreen python tests/test_specr.py
"""
import os
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)

from specfile import SpecDataFile, scan_stats, split_labels  # noqa: E402

FAILS = []


def check(label, cond, extra=""):
    print(("  PASS  " if cond else "  FAIL  ") + label + (f"   {extra}" if extra else ""))
    if not cond:
        FAILS.append(label)


# =============================================================================
print("\n=== 1. parse the reference example ===")
ref = SpecDataFile(os.path.join(HERE, "reference_example.spe"))
ref.refresh()
check("4 scans found", len(ref) == 4, [s.number for s in ref.scans])
check("scan numbers", [s.number for s in ref.scans] == [3757, 3758, 3762, 3763])
check("two-space #S -> command intact",
      ref[1].command == "my_plan()", ref[1].command)
check("motor names from #O", ref[0].motor_names == ["attenE"], ref[0].motor_names)
check("motor position from #P", ref[0].motor_positions == [-0.01], ref[0].motor_positions)
check("empty scan has no rows", ref[0].npoints == 0)
check("exit_status abort", ref[0].exit_status == "abort", ref[0].exit_status)
check("exit_status fail", ref[2].exit_status == "fail", ref[2].exit_status)
check("exit_status success", ref[3].exit_status == "success", ref[3].exit_status)

s = ref[3]
check("#N 12 == len(#L)", s.ncols == len(s.labels) == 12, (s.ncols, len(s.labels)))
check("labels with embedded single space preserved",
      "SMS granite" in s.labels and "OSS underhang" in s.labels, s.labels[3:6])
check("10 data rows", s.npoints == 10, s.npoints)
check("True -> 1.0", list(s.column("tc32E_ch_connected")) == [1.0] * 10)
check("numeric column intact", abs(s.column("AIR")[0] - 22.08344841003418) < 1e-9,
      s.column("AIR")[0])
check("negative sentinel intact", s.column("SMS granite")[0] == -9999.0)
check("long #MD iconfig captured", "iconfig" in ref[0].meta)
check("#MD uid captured", ref[0].meta.get("uid", "").startswith("7797fce2"))

# =============================================================================
print("\n=== 2. split_labels rules ===")
check("splits on 2+ spaces", split_labels("a  b  c") == ["a", "b", "c"])
check("keeps single space inside label",
      split_labels("Epoch  SMS granite  x") == ["Epoch", "SMS granite", "x"])
check("single-space fallback when #N matches",
      split_labels("a b c", ncols=3) == ["a", "b", "c"])
check("no fallback when #N disagrees",
      split_labels("a b  c", ncols=2) == ["a b", "c"])

# =============================================================================
print("\n=== 3. ophyd_scan SpecFile writer round-trips ===")
# Import ophyd_scan's SpecFile without pulling in EPICS/oregistry.
src = open(os.path.join(PKG, "..", "src", "id8_common", "plans", "align",
                        "ophyd_scan.py")).read()
start = src.index("def _fmt(")
end = src.index("# =============================================================================\n"
                "# Paths, scan numbering")
mod = types.ModuleType("ophyd_spec_writer")
mod.__dict__.update(
    {"math": __import__("math"), "os": os, "time": __import__("time"),
     "datetime": __import__("datetime"), "getpass": __import__("getpass"),
     "socket": __import__("socket"),
     "SPEC_TIME_FORMAT": "%a %b %d %H:%M:%S %Y"}
)
exec(compile(src[start:end], "ophyd_scan_excerpt", "exec"), mod.__dict__)

tmp = tempfile.mkdtemp()
path = os.path.join(tmp, "pope202607.spec")
SpecFile = mod.SpecFile
motors = ["huber_nu", "huber_delta", "huber_mu", "huber_eta", "huber_chi",
          "huber_phi", "huber_x", "huber_y", "huber_z"]
positions = [0.0, 30.00018, 0.1, 1.2, 2.3, 3.4, 4.5, 5.6, 6.7]
labels = ["huber_delta", "Epoch", "Epoch_float",
          "lambda2M_stats2_total", "lambda2M_stats3_total", "lambda2M_stats1_total"]

sf = SpecFile(path)
sf.write_file_header(motors)
sf.start_scan(173, "dscan_ophyd(huber_delta, -0.5, 0.5, 5, 1.0, det=lambda2M)",
              labels, metadata={"uid": "abc-123", "image_file": "A0173_x_a0001.h5"},
              motor_positions=positions,
              comments=["plan_type = function", "uid = abc-123"])
for i in range(5):
    sf.add_point([29.5 + 0.25 * i, i, i * 1.01, 100 + i, 200 + i, 300 + i])
sf.end_scan("success")

# a second scan appended to the same file must NOT repeat the file header
sf2 = SpecFile(path)
sf2.write_file_header(motors)
sf2.start_scan(174, "dscan_ophyd(huber_delta, -0.1, 0.1, 3, 1.0, det=lambda2M)",
               labels, motor_positions=positions)
for i in range(3):
    sf2.add_point([29.9 + 0.1 * i, i, i * 1.01, 10 + i, 20 + i, 30 + i])
sf2.end_scan("success")

w = SpecDataFile(path)
w.refresh()
check("exactly one #F block", len(w.file_headers) == 1, w.file_headers)
check("two scans", len(w) == 2, [x.number for x in w.scans])
check("scan numbers 173, 174", [x.number for x in w.scans] == [173, 174])
check("labels round-trip", w[0].labels == labels, w[0].labels)
check("#N matches", w[0].ncols == len(labels))
check("9 motor names via #O0/#O1 wrap", w[0].motor_names == motors, w[0].motor_names)
check("9 motor positions via #P0/#P1 wrap",
      len(w[0].motor_positions) == 9 and abs(w[0].motor_positions[1] - 30.00018) < 1e-9,
      w[0].motor_positions)
check("5 rows then 3 rows", (w[0].npoints, w[1].npoints) == (5, 3))
check("data values intact", list(w[0].column("lambda2M_stats1_total")) ==
      [300, 301, 302, 303, 304])
check("exit_status parsed", w[0].exit_status == "success")
check("image_file metadata present", w[0].meta.get("image_file") == "A0173_x_a0001.h5")
check("command intact", w[0].command.startswith("dscan_ophyd(huber_delta"), w[0].command)

# =============================================================================
print("\n=== 4. incremental (live) parsing ===")
live = os.path.join(tmp, "live.spec")
lf = SpecFile(live)
lf.write_file_header(motors)
lf.start_scan(200, "dscan_ophyd(live)", labels, motor_positions=positions)
reader = SpecDataFile(live)
reader.refresh()
check("scan visible before any point", len(reader) == 1 and reader[0].npoints == 0)

seen = []
for i in range(6):
    lf.add_point([i, i, i * 1.0, i * 2, i * 3, i * 4])
    changed = reader.refresh()
    seen.append((changed, reader[0].npoints))
check("each append is seen exactly once",
      seen == [(True, 1), (True, 2), (True, 3), (True, 4), (True, 5), (True, 6)], seen)
check("no-op refresh returns False", reader.refresh() is False)

# partial trailing line must not be parsed until its newline arrives
with open(live, "a") as f:
    f.write("999 999 999 999 999")   # deliberately unterminated
reader.refresh()
check("torn row ignored until newline", reader[0].npoints == 6, reader[0].npoints)
with open(live, "a") as f:
    f.write(" 999\n")
reader.refresh()
check("row appears once completed", reader[0].npoints == 7, reader[0].npoints)

lf.end_scan("success")
reader.refresh()
check("exit_status seen live", reader[0].exit_status == "success")

# Regression: a poll landing mid-line must not double-count the fragment.
# The first version rewound the offset over the partial line AND kept it in a
# _carry buffer, so "de" + "def" parsed as "dedef" -- rows silently corrupted
# whenever the reader polled between two writes. Likely over NFS at 1 Hz.
frag = os.path.join(tmp, "frag.spec")
ff = SpecFile(frag)
ff.write_file_header(motors)
ff.start_scan(300, "fragment probe", labels, motor_positions=positions)
fr = SpecDataFile(frag)
fr.refresh()
raw = open(frag, "a")
def row(i):
    return f"{10 + i} {20 + i} {30 + i} {40 + i} {50 + i} {60 + i}"
# The corrupting case is a single read that SPANS a newline and then ends
# mid-line -- i.e. the reader is behind and picks up a complete row plus the
# start of the next. (A read with no newline at all takes a different branch
# that happened to be self-consistent, which is why this went unnoticed.)
for i in range(0, 5, 2):
    raw.write(row(i) + "\n"); raw.flush()
    if i + 1 < 5:
        raw.write(row(i + 1)[:6]); raw.flush()   # ...plus half of the next row
    fr.refresh()                                  # <-- poll spans \n, ends mid-line
    if i + 1 < 5:
        raw.write(row(i + 1)[6:] + "\n"); raw.flush()
        fr.refresh()
raw.close()
fr.refresh()
got = [list(r) for r in fr[0].rows]
want = [[10.0 + i, 20.0 + i, 30.0 + i, 40.0 + i, 50.0 + i, 60.0 + i] for i in range(5)]
check("mid-line polls produce exactly 5 rows", len(got) == 5, len(got))
check("no fragment duplication in any row", got == want,
      got[0] if got else "no rows")

# =============================================================================
print("\n=== 5. scan_stats (params.m port) ===")
import numpy as np  # noqa: E402
x = np.linspace(-5, 5, 101)
y = np.exp(-(x ** 2) / (2 * 1.0 ** 2))          # sigma=1 -> FWHM = 2.3548
px, py, com, fwhm, center = scan_stats(x, y)
check("peak at 0", abs(px) < 1e-9, px)
check("peak value 1", abs(py - 1.0) < 1e-12, py)
check("COM at 0", abs(com) < 1e-9, com)
check("FWHM ~ 2.355", abs(fwhm - 2.3548) < 0.02, fwhm)
check("FWHM centre at 0", abs(center) < 1e-9, center)
check("empty input is safe", np.isnan(scan_stats([], [])[0]))

# =============================================================================
print("\n=== 6. headless GUI: load, plot, live update ===")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PyQt5 import QtWidgets  # noqa: E402
import specr  # noqa: E402

app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
win = specr.SpecrWindow(path=path)
check("GUI loaded file", win.specfile is not None and len(win.specfile) == 2)
check("newest scan selected by default", win.selection == [1], win.selection)
check("X defaults to first column", win.x_combo.currentText() == "huber_delta",
      win.x_combo.currentText())
check("Y defaults to last column (stats1 per SPEC convention)",
      win.y_combo.currentText() == "lambda2M_stats1_total", win.y_combo.currentText())
check("a line was plotted", len(win.axes.get_lines()) == 1)
check("xlabel carries Peak/COM/FWHM", "FWHM" in win.axes.get_xlabel())

# now drive a live append and tick the monitor by hand
win.open_file(live)
n_before = len(win.axes.get_lines()[0].get_xdata())
lf2 = SpecFile(live)
lf2.start_scan(201, "dscan_ophyd(second live scan)", labels, motor_positions=positions)
for i in range(4):
    lf2.add_point([i, i, i * 1.0, i * 5, i * 6, i * 7])
win._on_monitor_tick()
check("monitor picked up the new scan", win.selection == [1], win.selection)
check("new scan plotted with 4 points",
      len(win.axes.get_lines()[0].get_xdata()) == 4,
      len(win.axes.get_lines()[0].get_xdata()))
lf2.add_point([9, 9, 9.0, 45, 54, 63])
win._on_monitor_tick()
check("incremental point appears on the plot",
      len(win.axes.get_lines()[0].get_xdata()) == 5,
      len(win.axes.get_lines()[0].get_xdata()))

png = os.path.join(tmp, "live_plot.png")
win.figure.savefig(png, dpi=80)
check("canvas renders to PNG", os.path.getsize(png) > 5000, os.path.getsize(png))
print(f"        (rendered {png})")

# style switching must not raise
for i, style in enumerate(specr.PLOT_STYLES):
    win.style_combo.setCurrentIndex(i)
check("all plot styles applied", win.style_combo.currentText() == "logxy")
win.style_combo.setCurrentIndex(0)

win.show_scan_text.__self__  # attribute exists
check("motor table available to GUI", len(win._primary_scan().motor_table()) == 9,
      len(win._primary_scan().motor_table()))

print("\n" + ("ALL CHECKS PASSED" if not FAILS else f"FAILURES: {FAILS}"))
sys.exit(1 if FAILS else 0)
