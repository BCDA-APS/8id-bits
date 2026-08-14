# Ophyd scans + live SPEC viewer — user guide

Two pieces that work together at 8-ID-E:

| | what it is | where |
|---|---|---|
| `dscan_ophyd()` | a relative motor scan written in **plain Ophyd** — no Bluesky RunEngine, plans, or documents. Writes a SPEC-format text file, one row per point. | [src/id8_common/plans/align/ophyd_scan.py](src/id8_common/plans/align/ophyd_scan.py) |
| `specr_py` | a Python/Qt rewrite of the MATLAB `specr` GUI that plots those SPEC files and updates **while the scan runs**. | [specr_py/](specr_py/) |

The scan appends one complete row per point and closes the file each time. The
viewer only ever reads complete lines. That is why it is safe to watch a scan as
it happens — and why the scan writes a small text file instead of having the
viewer read the detector `.h5`, which the IOC still holds open.

---

## Quick start

**On amber** (where the IOCs and `/gdata` live):

```bash
./start_bluesky.sh                     # or: conda activate 8id_bits; ipython -i -c "from id8_common.startup import *"
```

```python
from id8_common.plans.align.ophyd_scan import dscan_ophyd
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M, att_ratio=5)
```

**On kouga** (or any workstation with a display):

```bash
conda activate 8id_bits
python ~/bluesky/specr_py/specr.py ~/bluesky/<experiment>.spec --watch
```

`--watch` opens the file and starts the monitor immediately. The plot fills in
point by point.

> **Call it directly. Do not wrap it in `RE()`.** `dscan_ophyd` is an ordinary
> function, not a Bluesky plan. `RE(dscan_ophyd(...))` will not work.

---

## The scan

```python
dscan_ophyd(motor, rel_begin, rel_end, num_pts, count_time,
            det=None, att_ratio=1e6, save_img=1, spec_path=None,
            set_attenuation=True, beam_control=True)
```

| argument | meaning |
|---|---|
| `motor` | any Ophyd positioner, e.g. `huber.delta`, `huber.eta`, `huber.x` |
| `rel_begin`, `rel_end` | start/end **relative to the current position** |
| `num_pts` | number of points (inclusive of both ends) |
| `count_time` | detector exposure per point, seconds |
| `det` | `eiger4M` (default) or `lambda2M` |
| `att_ratio` | attenuation ratio, same meaning as in `dscan()` |
| `save_img` | `1` writes detector images and burns a measurement number; `0` scans without saving |
| `spec_path` | override the SPEC file location |
| `set_attenuation` | `False` leaves the attenuator exactly as-is |
| `beam_control` | `False` skips `pre_align()`, `PIND_status()` and the shutter |

Returns `(positions, columns)` — `positions` is the array of **measured**
readbacks, `columns` is `{column_label: [values]}`.

```python
pos, cols = dscan_ophyd(huber.eta, -0.2, 0.2, 21, 0.5, det=lambda2M, att_ratio=7)
import numpy as np
print("peak at", pos[np.argmax(cols["lambda2M_stats1_total"])])
```

### Common variants

```python
# quiet dry run: only the scanned motor and the detector move.
# no attenuator, no PIND, no shutter.
dscan_ophyd(huber.delta, -0.1, 0.1, 11, 0.1, det=lambda2M,
            set_attenuation=False, beam_control=False)

# no images: nothing written to /gdata, measurement_num NOT incremented
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M, save_img=0)

# send the SPEC output somewhere else
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M,
            spec_path="/home/beams/8IDIUSER/bluesky/scratch.spec")
```

`tetramm` raises `NotImplementedError` — its trigger path and HDF stream mode
differ and were never verified here. Use `dscan()` in `scan_8id.py` for tetramm.

### What it does to the hardware

With defaults (`set_attenuation=True, beam_control=True`) it follows `dscan()`
in [scan_8id.py](src/id8_common/plans/align/scan_8id.py) exactly: `pre_align()`,
`att(att_ratio)`, `PIND_status(0)`, arm the detector, open the shutter, step and
trigger, then close the shutter, restore the trigger mode, and return the motor
to its starting position. Both flags exist so you can opt out of everything
except the motor and the detector.

---

## Where files go

| | path |
|---|---|
| SPEC file | `~/bluesky/{experiment_name}.spec` |
| Detector images | `{mount_point}{cycle_name}/{experiment_name}/data/bluesky/{prefix}.h5` |

e.g. `~/bluesky/pope202607.spec` and
`/gdata/dm/8ID/8IDE/2026-2/pope202607/data/bluesky/A0177_HEA-15GPa-3x3Grid_a0005.h5`.

**The SPEC file is deliberately not next to the images.** `/gdata` is mounted on
amber but **not** on kouga, whereas `~/bluesky` is shared by both. Putting the
small text file in the home directory is what lets you run the viewer on a
workstation while the scan runs on amber. To change it, edit `SPEC_DIR` at the
top of `ophyd_scan.py`.

Both paths come from `pv_registers` (`mount_point`, `cycle_name`,
`experiment_name`), so they follow the run cycle automatically.

### Scan numbering

The SPEC scan number is taken from the image prefix, so they always line up:

```
#S 177  dscan_ophyd(huber_delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
        ^^^                    A0177_HEA-15GPa-3x3Grid_a0005.h5
```

`gen_folder_prefix()` increments `pv_registers.measurement_num` once per scan
(only when `save_img=1`). Scans append to one SPEC file per experiment; the
`#F` header block is written only when the file is created.

---

## Columns

```
#L huber_delta  huber_delta_setpoint  Epoch  Epoch_float  lambda2M_stats2_total  lambda2M_stats3_total  lambda2M_stats1_total
```

| column | meaning |
|---|---|
| `<motor>` | **measured readback** at the time of the point |
| `<motor>_setpoint` | where the motor was told to go |
| `Epoch` / `Epoch_float` | seconds since scan start, rounded / full precision |
| `<det>_stats1..3_total` | ROI totals from the detector's stats plugins |

Two deliberate choices:

- **Readback and setpoint are both recorded.** Logging only the commanded value
  would draw a perfect ramp even if the motor never moved. Comparing the two
  columns is the fastest way to confirm a scan really happened.
- **`stats1` is written last**, following SPEC's "primary detector goes in the
  last column" convention. Readers that default the Y axis to the final column —
  including this GUI and the MATLAB one — then plot the right thing with no
  clicking.

---

## The viewer

```bash
python ~/bluesky/specr_py/specr.py                        # then File > Open
python ~/bluesky/specr_py/specr.py <file>.spec            # open a file
python ~/bluesky/specr_py/specr.py <file>.spec --watch    # open and start monitoring
```

Needs a display — run it on a workstation, or `ssh -X`.

| feature | how |
|---|---|
| **Scan Monitor** | `Ctrl+M` or View ▸ Scan Monitor. Polls the file and redraws as points land. |
| **Monitor period** | Tools ▸ Settings, default 2 s. If a refresh takes longer than the interval the interval backs off, so polls never queue up. |
| **Erase mode** | View ▸. On: monitor shows only the newest scan. Off: accumulates scans sharing the same columns. |
| **Select Scan** | Multi-select list. Scans must share a column layout to overlay. |
| **Step scans** | `Ctrl+←` / `Ctrl+→` |
| **X / Y columns** | Dropdowns. Default X = first column, Y = last. |
| **Peak / COM / FWHM** | Under the X axis, recomputed on every redraw. |
| **Plot styles** | linear / logx / logy / logxy |
| **Show Current Scan** | Tools ▸ — the raw text of the scan block. |
| **Show Motor Positions** | Tools ▸ — the `#O`/`#P` table (all 9 huber axes at scan start). |
| **Save Displayed Data** | File ▸ `Ctrl+S` — two-column text of the current X/Y. |
| Grid, legend, axis inversion | View ▸ and Tools ▸ |

### Typical session

1. Start the viewer on kouga with `--watch` pointed at this experiment's `.spec`.
2. Run scans on amber. Each new scan appears automatically and plots live.
3. Turn erase mode **off** to overlay successive scans and compare.

---

## Reading SPEC files from your own scripts

`specfile.py` has no GUI dependencies, so it is useful on its own:

```python
import sys; sys.path.insert(0, "/home/beams/8IDIUSER/bluesky/specr_py")
from specfile import SpecDataFile, scan_stats

sd = SpecDataFile("/home/beams/8IDIUSER/bluesky/pope202607.spec")
sd.refresh()                       # parse; call again later to pick up new data

print(len(sd), "scans:", [s.number for s in sd.scans])

s = sd.by_number(177)              # or sd[-1], or sd.last()
x = s.column("huber_delta")
y = s.column("lambda2M_stats1_total")

peak_x, peak_y, com, fwhm, center = scan_stats(x, y)
print(f"peak {peak_y:g} @ {peak_x:g}, FWHM {fwhm:g}")

print(s.command, s.date, s.exit_status)
print(s.meta["image_file"])        # links the scan to its .h5
print(s.motor_table())             # [(name, position), ...] from #O/#P
```

`refresh()` is incremental — it consumes only bytes appended since the last
call and returns `True` if anything changed, so it is cheap to poll.

---

## Troubleshooting

**The curve is flat and `<motor>` barely changes while `<motor>_setpoint` sweeps.**
The motor is not physically moving. Ophyd will *not* raise — a disabled drive
accepts the setpoint and reports `DMOV=1` immediately. Check:

```python
from epics import caget
caget("8ideSoft:CR8-E1:m5.CNEN", as_string=True)   # 'Disable' -> drive is off
caget("8ideSoft:CR8-E1:m5.MSTA")                    # want the POSITION bit set
```

This has happened in practice. Comparing the readback and setpoint columns is
the quickest check, and the `.h5` `Huber_*` NDAttributes are an independent
second opinion — they will also be constant if the motor did not move.

**The GUI shows nothing new during a scan.** Confirm it is pointed at
`~/bluesky/<experiment>.spec`, not a copy, and that the monitor is on
(`Ctrl+M`). The title bar shows the file and the status bar shows
`[running]` or the exit status.

**"No such file or directory" for the images on kouga.** Expected — `/gdata` is
not mounted there. Only the SPEC file is readable from kouga; open the `.h5` on
amber.

**Counts are all zero.** Check the ring current and the shutter before
suspecting the software:

```python
caget("S:SRcurrentAI")
```

**A scan aborts partway.** The scan block is still closed off with
`#C ... exit_status = aborted` (or `error`), so the file stays parseable and the
partial data is kept.

---

## Tests

```bash
cd ~/bluesky/specr_py
QT_QPA_PLATFORM=offscreen python tests/test_specr.py
```

Runs headless, needs no hardware. Covers parsing a reference file with the
awkward parts of real Bluesky SPEC output (repeated `#F` blocks, labels
containing single spaces, bare `True` in a numeric column, `#N 0` scans),
round-tripping what `ophyd_scan.py` writes, incremental live reads including a
torn trailing row, the Peak/COM/FWHM port, and a headless GUI live-update run.

Hardware scaffolding, for re-running the beamline checks:

| script | what it does |
|---|---|
| `tests/_stage_a_check.py` | brings up the session and prints the resolved config, motor limits and detector state. **Moves nothing.** |
| `tests/_stage_b_scan.py` | the test scan, with a limits assertion that refuses to run out of range |
| `tests/_live_watch.py` | drives the real GUI headlessly against a growing file and reports whether it updated live |

---

## Notes and limits

- **Not Bluesky.** No RunEngine, no documents, nothing in the Tiled catalog. If
  you need the run in databroker/Tiled, use `dscan()` in `scan_8id.py`.
- `ophyd_scan.py` resolves devices from `oregistry` at import time, so import it
  after `make_devices()` has run — i.e. after `from id8_common.startup import *`.
- The SPEC writer never raises into a scan. A full disk or bad path logs a
  warning and the scan continues.
- One SPEC file per experiment, appended across sessions. Scan numbers come from
  `measurement_num`, so they are unique as long as that register is not reset.
- The Lambda frame wait is bounded (`count_time * 5 + 10` s). `dscan()` in
  `scan_8id.py` spins forever if a trigger is missed; this one warns and moves on.
