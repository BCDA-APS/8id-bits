# Ophyd scans and the CSV file template

How to run an alignment scan without Bluesky, and how to control what the scan
writes into its data file — by editing a YAML template, not Python.

Two audiences:

- **Just running scans?** Read [Quick start](#quick-start) and
  [`dscan_ophyd()`](#dscan_ophyd). That is all you need.
- **Want a different column, or another PV recorded?** Read
  [Defining the CSV file structure](#defining-the-csv-file-structure).

## What is here

| file | what it is | edit it when |
|---|---|---|
| [ophyd_scan.py](ophyd_scan.py) | the scan itself — moves the motor, arms the detector, reads counters | you are adding a new scan (`d2scan`, `mesh`, …) |
| [scan_csv.py](scan_csv.py) | the CSV writer and a reader. Knows the file *format*, moves no hardware | almost never |
| [../../configs/scan_csv_template.yml](../../configs/scan_csv_template.yml) | what goes in the file | you want a different column or another PV in the header |

The split is deliberate: `ophyd_scan.py` is about hardware, `scan_csv.py` is
about a file. A future `d2scan` gets its file format for free, and someone
reading the scan code is not reading string formatting.

`scan_csv.py` imports `oregistry` (to look up devices named in the template) but
nothing from Bluesky. `ophyd_scan.py` imports no Bluesky either — that is the
point of the module.

## Quick start

```python
# start the beamline session on Pearl with:  ~/bin/start_bluesky.sh
# (it activates 8id_bits, then `from id8_common.startup import *`)
#
# NB *not* start_bluesky_8ide.sh -- that one is stale. It still does
# `from id8_e.startup import *`, and id8_e was merged into id8_common.
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
```

That writes, in the experiment's `data/bluesky` folder:

```
A0201_Test_a0010.h5     the detector images
A0201_Test_a0010.csv    motor positions and counters, one row per point
```

**One CSV per scan**, named after that scan's `.h5`. Nothing is appended to a
previous scan's file, so there is no scan numbering to keep track of and no way
for two scans to interfere.

The file is closed after every point, so you can open it in a viewer, in Excel,
or in the specr_py viewer *while the scan is still running*. That is the whole
reason a scan writes a text file at all — reading a detector `.h5` mid-write is
not safe, reading an append-only text file is.

Call it **directly**. Do not wrap it in `RE()`; there is no RunEngine involved.

## `dscan_ophyd()`

```python
dscan_ophyd(motor, rel_begin, rel_end, num_pts, count_time,
            det=eiger4M, att_ratio=1e6, save_img=1, comment="")
```

| argument | meaning |
|---|---|
| `motor` | any ophyd positioner: `huber.delta`, `sample.x`, `rheometer.y`, … |
| `rel_begin`, `rel_end` | scan range **relative to where the motor is now** |
| `num_pts` | number of points (inclusive of both ends) |
| `count_time` | seconds per point |
| `det` | `eiger4M`, `lambda2M`, or `tetramm1` |
| `att_ratio` | attenuation ratio passed to `att()` |
| `save_img` | `1` writes the detector `.h5`; `0` scans without it (the CSV is still written) |
| `comment` | free text describing this scan — see below |

Returns the `ScanCsv` object: `.path` is the file written, `.data` is
`{column_label: [values]}` if you want the numbers in the session.

The motor is returned to its starting position when the scan finishes — or when
you interrupt it.

**Any positioner works.** The scan only uses `.name`, `.position` and `.move()`.
A motor does *not* have to be listed in the template to be scanned; the template's
motor list only says whose *start position* gets recorded in the header.

### Annotating a scan

```python
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M,
            comment="3x3 grid spot 5, after realigning KB")
```

lands on the third line of the file, where it is easy to find:

```
start_time,2026-08-26 14:32:07
h5_file,/gdata/dm/8ID/8IDE/2026-2/pope202607/data/bluesky/A0201_Test_a0010.h5
comment,"3x3 grid spot 5, after realigning KB"
```

Commas and quotes inside a comment are safe — the writer quotes the field the
way any CSV reader expects. Newlines are flattened to spaces so a pasted
paragraph cannot break the file's structure. Leave `comment` out and the line is
omitted entirely.

Because it is plain text, `grep -l "spot 5" *.csv` in the data folder finds the
scan.

### What it does to the hardware

Before the loop: `pre_align()`, `att(att_ratio)`, `PIND_status(0)`, then the
detector is armed. During the loop: move, trigger, wait `count_time`, read.
After: acquisition stops, the image file closes, the beam is blocked, trigger
modes are restored, `shutteroff()`, and the motor drives back to start.

### Interrupting with Ctrl+C

```
 12            30.15            30.15  0.5012  ...
^C  Scan aborted at point 12/41.
Scan file closed (aborted, 12 points): /gdata/.../A0201_Test_a0010.csv
# images captured:  12
```

The teardown halts acquisition, closes the image file, blocks the beam, stops
the motor, drives it back to where it started, and closes the CSV with
`#END,aborted,12` so the partial scan is still a valid file. The
`KeyboardInterrupt` is caught, so you get the message above rather than a
traceback, and `dscan_ophyd` still returns its `ScanCsv`.

That teardown is Sam's `finally:` block from `scan_8id.py`, with two lines added
in front of the return move:

```python
motor.stop()                        # added
time.sleep(0.2)                     # added
motor.move(start_pos, wait=True)    # Sam's line
```

**Why.** A ^C interrupts a move and leaves its `MoveStatus` unfinished; ophyd
then refuses the next `move()` on that motor because "another set() is still in
progress". That is what the abort test on 2026-08-26 ran into: everything else
worked — the motor stopped (`DMOV=1`, `SPMG=Go`, no limit violation), the
detector went idle, the beam was blocked, and the file closed as
`#END,aborted,7` — but the motor stayed at 30.176 instead of driving back to
31.001. `stop()` clears the stale status so the return move is accepted.

> **⚠ Diagnosed, not yet retested on hardware.** The two lines above are the
> fix for that failure, but no one has given the beamline another ^C since.
> **Check the motor position after a ^C before starting the next scan.**

Two things this file deliberately does *not* do, because `scan_8id.py` does not
either:

- **A second Ctrl+C during cleanup is not caught.** It interrupts the teardown,
  possibly before the return move.
- **The frame wait still runs after an abort** (eiger branch), so an abort there
  can take up to `num_pts * count_time + 10` s to return.

---

# Defining the CSV file structure

## Anatomy of the file

```
start_time,2026-08-26 14:32:07                        <- lines:
h5_file,/gdata/.../A0201_Test_a0010.h5                <- lines:
comment,"3x3 grid spot 5"                             <- lines:  (skipped if empty)
command,"dscan_ophyd(huber_delta, -0.5, 0.5, 41, 1.0, det=lambda2M)"
scan_type,dscan_ophyd
motor,huber_delta
detector,lambda2M
num_points,41
count_time,1.0
beamline,8-ID-E
epoch,1787793756.041
attenuation,10                                        <- lines:  (source: a PV)
huber_nu,-0.00049                                     <- lines:  (source: a PV)
huber_delta,30.0004
...
roi1,730,983,100,10                                   <- lines:  (values: 4 PVs)
#DATA                                                 <- marker:
huber_delta,huber_delta_setpoint,elapsed_time,...     <- columns:  NAMES
30.0,30.0,0.0,1234                                    <- columns:  one row per point
30.025,30.025,1.05,1301
#END,success,41                                       <- automatic
```

Two parts, two template keys:

- **`lines:`** — things that *do not* change during the scan. Written once, above
  the marker, as `label,value`.
- **`columns:`** — things that *do* change. One column each, one row per point,
  below the marker.

`#END,<status>,<points>` is written by the scan, not the template. Status is
`success`, `aborted`, or `error`. Its absence means the scan is still running —
which is how a live viewer knows to keep polling.

## The three ways to give a line a value

```yaml
lines:
  - {label: beamline,    value: "8-ID-E"}                    # beamline,8-ID-E
  - {label: attenuation, source: filter_8ide.attenuation.readback}   # attenuation,10
  - {label: roi1, values: [det.roi1.min_xyz.min_x, det.roi1.min_xyz.min_y,
                           det.roi1.size.x, det.roi1.size.y]}        # roi1,730,983,100,10
```

| key | takes | produces |
|---|---|---|
| `value:` | text, with `{...}` substitutions | `label,text` |
| `source:` | **one** dotted path to a PV | `label,value` |
| `values:` | a **list** of dotted paths | `label,v1,v2,v3` |

All are read once, at scan start.

## Dotted paths

```
det.stats1.total       the detector passed to THIS scan (eiger4M / lambda2M / tetramm1)
motor.velocity         the motor being scanned
huber.delta            any device in the oregistry
sample.x
filter_8ide.attenuation.readback
tetramm1.current1.mean_value
```

`det.` and `motor.` are why one template covers every detector and every motor.

**Always name the individual signal, never a parent device.** `lambda2M.roi1`
reads back a namedtuple of every field; `lambda2M.roi1.size.x` is a number. A
non-numeric value is written as-is and a one-line warning names the column, so
the mistake shows up rather than silently producing junk.

## `{...}` substitutions

Usable in any `label:` and in any `value:`.

| | |
|---|---|
| `{start_time}` | `2026-08-26 14:32:07` |
| `{epoch}` | the same instant as a unix timestamp |
| `{h5_file}` | **full path** of this scan's detector `.h5` (see below); empty when `save_img=0` |
| `{comment}` | the `comment="..."` argument |
| `{command}` | `dscan_ophyd(huber_delta, -0.5, 0.5, 41, 1.0, det=lambda2M)` |
| `{scan_type}` | `dscan_ophyd` |
| `{motor}` | `huber_delta` |
| `{det}` | `lambda2M` |
| `{num_points}`, `{count_time}` | as passed |

### Where `{h5_file}` comes from

It is **not** an argument to `dscan_ophyd()`. The scan builds it:

1. `gen_folder_prefix()` (in `plans/acquire/ad_acq.py`) combines
   `pv_registers.header` + `measurement_num` + `sample_name` + the current
   attenuation into e.g. `A0201_Test_a0010`, incrementing `measurement_num`.
2. That prefix names both files, in the folder built from
   `pv_registers.mount_point` + `cycle_name` + `experiment_name` + `/data/bluesky`.

So the CSV and the `.h5` always share a name, and each scan gets exactly one
measurement number.

## Per-entry flags

```yaml
  - {label: comment,     value: "{comment}", skip_if_empty: true}
  - {label: rheometer_x, source: rheometer.x, optional: true}
```

- **`optional: true`** — if the device or signal does not exist, skip the entry
  and print one note. Use it for anything that is not on every station or not on
  every detector. Everything in the shipped template's `columns:` is optional,
  which is how one template serves eiger, lambda and tetramm.
- **`skip_if_empty: true`** — leave the line out when the value comes back blank.

An entry **not** marked `optional` whose source is missing is an error raised
**before the scan starts** — not halfway through with the beam on.

## Columns

```yaml
columns:
  - {label: "{motor}",            source: motor}
  - {label: "{motor}_setpoint",   source: motor_setpoint}
  - {label: elapsed_time,         source: elapsed}
  - {label: tetramm1_current1,    source: tetramm1.current1.mean_value, optional: true}
  - {label: "{det}_stats1_total", source: det.stats1.total, optional: true}
```

Four sources are computed by the scan rather than read from a PV:

| source | value |
|---|---|
| `motor` | where the motor actually is (readback) |
| `motor_setpoint` | where the motor was told to go |
| `elapsed` | seconds since the scan started |
| `epoch` | unix timestamp of the point |

Anything else is a dotted path, read fresh at every point.

**Order matters for plotting.** Viewers conventionally default X to the first
column and Y to the last, so keep the scanned motor first and your primary
counter last — which is why the shipped template ends with `stats1`.

**Labels must be unique.** A duplicate makes the file ambiguous, so it is
rejected before the scan starts. The usual cause is `{det}` or `{motor}`
rendering two entries to the same text: with `det=tetramm1`, a
`{det}_current1` column would collide with a `tetramm1_current1` column. The
error names the offending labels.

## Worked edits

**Record another PV in the header** — say the monochromator energy:

```yaml
lines:
  ...
  - {label: energy, source: mono.energy, optional: true}
```

(`mono.energy` is an `EpicsMotor`, and the writer reads `.position` from
anything that has one, so naming the positioner itself is right here.)

**Add a counter column** — put it before the last one so `stats1` stays the
default Y axis:

```yaml
columns:
  ...
  - {label: "{det}_stats2_total",  source: det.stats2.total, optional: true}
  - {label: tetramm2_sum_all,      source: tetramm2.sum_all.mean_value, optional: true}
  - {label: "{det}_stats1_total",  source: det.stats1.total, optional: true}
```

**Record a whole ROI on one line:**

```yaml
  - {label: roi2, values: [det.roi2.min_xyz.min_x, det.roi2.min_xyz.min_y,
                           det.roi2.size.x, det.roi2.size.y], optional: true}
```

**Stop recording the rheometer** — delete its three lines. The template is
re-read at the start of every scan, so the next scan picks the change up. No
restart, no `reload`.

## Trying a template without touching the default

```bash
cp ~/bluesky/src/id8_common/configs/scan_csv_template.yml ~/my_template.yml
# edit ~/my_template.yml
export ID8_SCAN_CSV_TEMPLATE=~/my_template.yml     # before starting bluesky
```

Resolution order is: the `template=` argument to `scan_csv.open_scan()`, then
`$ID8_SCAN_CSV_TEMPLATE`, then the packaged default. A scan using anything other
than the packaged default prints `NOTE: using CSV template ...`, so there is no
mystery about which one is in force.

---

## Where files go

```
<mount_point><cycle_name>/<experiment_name>/data/bluesky/
    A0201_Test_a0010.h5      images
    A0201_Test_a0010.csv     this scan
    A0202_Test_a0010.h5      next scan
    A0202_Test_a0010.csv
```

built from `pv_registers.mount_point`, `.cycle_name`, `.experiment_name` — the
same folder `save_images()` in `scan_8id.py` writes to.

### Finding the running scan from anywhere

Both names are published to EPICS as the scan starts:

```bash
caget -S 8ideSoft:StrReg21     # the .h5  file of the current/last scan
caget -S 8ideSoft:StrReg22     # the .csv file of the current/last scan
```

```python
pv_registers.scan_h5_file.get()
pv_registers.scan_csv_file.get()
```

These are 256-**character waveform** PVs, so a full path fits — which is also
why `caget` needs `-S`, or it prints the path as a list of character codes.
`scan_h5_file` is empty when `save_img=0`.

The specr_py viewer uses this for its `--live` flag and File ▸ Follow Live Scan:
read the register, open what it names, start monitoring. It is a convenience
only — the viewer's actual live update comes from polling the data folder, so it
works on machines with no channel access to `8ideSoft:`.

## Reading the file back

With pandas — skip everything down to and including the marker:

```python
import pandas

def read_scan(path, marker="#DATA"):
    with open(path) as f:
        skip = next(i for i, line in enumerate(f) if line.startswith(marker)) + 1
    return pandas.read_csv(path, skiprows=skip, comment="#")

df = read_scan("/gdata/.../A0201_Test_a0010.csv")
df.plot(x=df.columns[0], y=df.columns[-1])
```

Without pandas, `scan_csv.py` ships a reader that imports only the stdlib:

```python
from id8_common.plans.align.scan_csv import read_scan_csv

header, labels, columns, status = read_scan_csv(path)
header["h5_file"]      # ['/gdata/.../A0201_Test_a0010.h5']
header["roi1"]         # ['730', '983', '100', '10']
columns["huber_delta"] # [30.0, 30.025, ...]
status                 # 'success' | 'aborted' | 'error' | 'running'
```

`status == "running"` means no `#END` line yet — call it again in a second and
you have a live view. A half-written trailing row is ignored until it is
complete, so polling a scan in progress is safe.

## Adding another scan (d2scan, mesh, …)

Copy `dscan_ophyd`, change the loop, and reuse the same four lines:

```python
scan = scan_csv.open_scan(
    csv_file, motor, det=det, scan_type="d2scan_ophyd", command=command,
    num_points=num_pts, count_time=count_time, h5_file=h5_file, comment=comment,
)
scan.write_header()
# per point:
scan.add_point(setpoint)
# at the end (in a finally:):
scan.close(status)
```

For a two-motor scan, pass the second motor's position through
`extra={"motor2": motor2.name}`, and the template can then use `{motor2}` and a
`motor2.user_readback` column.

## Troubleshooting

| symptom | cause |
|---|---|
| `no device 'foo' in the oregistry` before the scan starts | typo in a template `source:`, or the device is not loaded. Add `optional: true` if it is genuinely not always present. |
| `CSV template has duplicate column labels` | two entries render to the same text once `{det}`/`{motor}` is filled in. Rename one. |
| `unknown substitution '{...}'` | a `{name}` that is not in the list above. The error prints the available ones. |
| `WARNING: column X is not a number` | the dotted path points at a device rather than a signal. Name the individual signal. |
| CSV written but no `.h5` | `save_img=0` (the scan says so on its first line), or the HDF plugin is not enabled. With `save_img=1` the scan prints `Scan folder created:` and the path it armed. |
| no `#END` line | the session was killed outright. Everything up to the last complete row is still readable. |
| `another set() is still in progress` on the next scan | should not happen — abort calls `motor.stop()` first. Report it. |

## Notes and limits

- **No RunEngine.** Nothing here writes to Tiled or databroker, and none of the
  data is in the catalog. The CSV *is* the record.
- **Not queueserver-safe.** `dscan_ophyd` is a plain function, not a plan, so
  the QS cannot run it. Use it interactively.
- The template is re-read at the start of every scan; the module-level device
  handles (`eiger4M`, `lambda2M`, `pv_registers`, …) are resolved once at import,
  so `ophyd_scan` must be imported after `make_devices()` — as `startup.py` does.
