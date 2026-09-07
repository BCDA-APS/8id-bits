# Ophyd scans and the CSV file template

How to run an alignment scan without Bluesky, and how to control what the scan
writes into its data file — by editing a YAML template, not Python.

Two audiences:

- **Just running scans?** Read [Quick start](#quick-start), [The six
  scans](#the-six-scans) and [Interrupting with Ctrl+C](#interrupting-with-ctrlc).
- **Want a different column, or another PV recorded?** Read
  [Defining the CSV file structure](#defining-the-csv-file-structure).

## What is here

| file | what it is | edit it when |
|---|---|---|
| [ophyd_scan.py](ophyd_scan.py) | the scan itself — moves the motor, arms the detector, reads counters | you are adding a new scan (it already has `dscan`, `ascan`, `d2scan`, `a2scan`, `dmesh`, `mesh` and the lups) |
| [scan_csv.py](scan_csv.py) | the CSV writer and a reader. Knows the file *format*, moves no hardware | almost never |
| [../../configs/scan_csv_template.yml](../../configs/scan_csv_template.yml) | what goes in the file | you want a different column or another PV in the header |

The split is deliberate: `ophyd_scan.py` is about hardware, `scan_csv.py` is
about a file. A new scan gets its file format for free, and someone
reading the scan code is not reading string formatting.

`scan_csv.py` imports `oregistry` (to look up devices named in the template) but
nothing from Bluesky. `ophyd_scan.py` imports no Bluesky either — that is the
point of the module.

## Quick start

Two sessions, both on Pearl, and **the scans have different names in each**:

```bash
~/bin/start_ophyd.sh      # Ophyd only.   dscan(...)        <- plain names
~/bin/start_bluesky.sh    # full Bluesky. dscan_ophyd(...)  <- suffixed names
```

```python
# in the Ophyd-only session
dscan(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)

# the same scan in the Bluesky session
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
```

See [Two sessions, two sets of names](#two-sessions-two-sets-of-names) for why.
Everything else on this page is identical in both — except the lups, which only
the Ophyd-only session has. The plain names are used from here on.

> NB **not** `start_bluesky_8ide.sh` — that one is stale. It still does
> `from id8_e.startup import *`, and `id8_e` was merged into `id8_common`.

A scan writes, in the experiment's `data/bluesky` folder:

```
A0113_Test_a1010041.h5     the detector images
A0113_Test_a1010041.csv    motor positions and counters, one row per point
```

**One CSV per scan**, named after that scan's `.h5`. Nothing is appended to a
previous scan's file, so there is no scan numbering to keep track of and no way
for two scans to interfere — a scan that finds its name already taken raises
`FileExistsError` before it moves anything, rather than concatenating itself
onto the older file.

The file is closed after every point, so you can open it in the BLUETELLA scan
viewer, or in Excel, *while the scan is still running*. That is the whole reason
a scan writes a text file at all — reading a detector `.h5` mid-write is not
safe, reading an append-only text file is.

Call it **directly**. Do not wrap it in `RE()`; there is no RunEngine involved.

## The six scans

```python
dscan (motor,  rel_begin,  rel_end,  num_pts, count_time, det=None, att_ratio=1e6, save_img=1, comment="")
ascan (motor,  abs_begin,  abs_end,  num_pts, count_time, det=None, att_ratio=7,   save_img=1, comment="")

d2scan(motor1, rel_begin1, rel_end1,
       motor2, rel_begin2, rel_end2, num_pts, count_time, det=None, att_ratio=7,   save_img=1, comment="")
a2scan(motor1, abs_begin1, abs_end1,
       motor2, abs_begin2, abs_end2, num_pts, count_time, det=None, att_ratio=7,   save_img=1, comment="")

dmesh (motor1, rel_begin1, rel_end1, num1,
       motor2, rel_begin2, rel_end2, num2,   count_time, det=None, att_ratio=7,   save_img=1, comment="")
mesh  (motor1, abs_begin1, abs_end1, num1,
       motor2, abs_begin2, abs_end2, num2,   count_time, det=None, att_ratio=7,   save_img=1, comment="")
```

| scan | axes | range | shape | default `det` |
|---|---|---|---|---|
| `dscan` | one motor | **relative** to where it is now | line | `eiger4M` |
| `ascan` | one motor | **absolute** | line | `eiger4M` |
| `d2scan` | two motors | **relative** | one **trajectory** — both axes sweep together, point *i* is `(positions1[i], positions2[i])` | `eiger4M` |
| `a2scan` | two motors | **absolute** | one **trajectory** | `eiger4M` |
| `dmesh` | two motors | **relative** | **raster**, `num1 * num2` points | `lambda2M` |
| `mesh` | two motors | **absolute** | **raster**, `num1 * num2` points | `lambda2M` |

On a raster **`motor1` is the outer (slow) axis and `motor2` the inner (fast)
one**: for each position of `motor1` the scan steps `motor2` across its whole
range. `num_pts` is not an argument to a raster — it is `num1 * num2`.

`d2scan`/`a2scan` are *not* rasters. They drive both motors along a single line
and deliberately write no grid metadata, so a viewer knows not to grid them.
Plot one in the 1-D viewer against whichever motor you like.

Arguments common to all six:

| argument | meaning |
|---|---|
| `motorN` | any ophyd positioner: `huber.delta`, `huber.nu`, `sample.x`, `rheometer.y`, … |
| `num_pts` / `num1`, `num2` | number of points (inclusive of both ends) |
| `count_time` | seconds per point |
| `det` | `eiger4M`, `lambda2M`, or `tetramm1`. `None` looks the default up **at call time**, so a detector that was offline at startup is not baked in |
| `att_ratio` | attenuation ratio passed to `att()`. Note `dscan`'s default is `1e6`, everything else `7` — Sam's defaults, kept |
| `save_img` | `1` writes the detector `.h5`; `0` scans without it (the CSV is still written, and still gets its own measurement number) |
| `comment` | free text describing this scan — see below |

Each returns the `ScanCsv` object: `.path` is the file written, `.data` is
`{column_label: [values]}` if you want the numbers in the session.

Every motor is returned to its starting position when the scan finishes — or
when you interrupt it. On `ascan`/`a2scan`/`mesh` that is where the axis stood
when you typed the command, not the first point of the absolute range.

**Any positioner works.** The scan only uses `.name`, `.position`, `.move()` and
`.stop()`. A motor does *not* have to be listed in the template to be scanned;
the template's motor list only says whose *start position* gets recorded in the
header. It does have to have a **column** — see
[Two-motor and raster scans](#two-motor-and-raster-scans).

### The lups

Six one-line wrappers around `dscan` for the axes people align most often. They
resolve the stage from the oregistry at call time and always pass `save_img=0`,
so they write a `.csv` and no `.h5`:

```python
x_lup      (rel_begin=-3,   rel_end=3,   num_pts=60, att_ratio=7, det=None, count_time=0.1, comment="")
y_lup      (rel_begin=-3,   rel_end=3,   num_pts=60, att_ratio=7, det=None, count_time=0.1, comment="")
huber_x_lup(rel_begin=-0.3, rel_end=0.3, num_pts=60, att_ratio=1, det=None, count_time=0.1, comment="")
huber_y_lup(rel_begin=-0.3, rel_end=0.3, num_pts=60, att_ratio=1, det=None, count_time=0.1, comment="")
rheo_x_lup (rel_begin=-3,   rel_end=3,   num_pts=60, att_ratio=7, det=None, count_time=0.1, comment="")
rheo_y_lup (rel_begin=-3,   rel_end=3,   num_pts=60, att_ratio=7, det=None, count_time=0.1, comment="")
```

`x_lup`/`y_lup` scan `sample.x`/`sample.y`, `huber_*_lup` the huber axes, and
`rheo_*_lup` the rheometer. All six default to **`tetramm1`**, not an area
detector — which is why the template's `tetramm1_*` columns matter (see the note
in [scan_csv_template.yml](../../configs/scan_csv_template.yml)).

**The lups above exist only in the Ophyd-only session.** `scan_8id.py` defines
`x_lup`, `y_lup`, `huber_x_lup`, `huber_y_lup`, `rheo_x_lup` and `rheo_y_lup`
too, as Bluesky generators, and `startup.py` star-imports them — so in the
Bluesky session those plain names are Sam's generators, to be run as
`RE(x_lup())`. There is no `x_lup_ophyd`: `startup.py` imports none of the lups
from `ophyd_scan.py`. To use these, start the Ophyd-only session, or call
`dscan_ophyd(sample.x, -3, 3, 60, 0.1, det=tetramm1, save_img=0)` directly.

### `auto_att()`

```python
auto_att(det, pilot_exptime=0.05, rate_limit=1e5, filter_factor=5.0,
         retry_max=10, grace_factor=0.25)
```

Not a scan and it writes no file. It takes short pilot exposures on `eiger4M` or
`lambda2M`, starting from maximum attenuation, and steps the transmission until
the peak count rate lands in `[rate_limit * grace_factor, rate_limit]`. The beam
is opened, so it runs under the same interrupt guard as a scan: a Ctrl+C blocks
the beam and restores `acquire_time`/`acquire_period` before returning.

In the Bluesky session it is `auto_att_ophyd`. Plain `auto_att` exists there
too, but it is `scan_8id.py`'s own copy — same signature, no interrupt guard.

## Two sessions, two sets of names

| session | started by | scan names |
|---|---|---|
| Ophyd only | `~/bin/start_ophyd.sh` | `dscan`, `ascan`, `d2scan`, `a2scan`, `dmesh`, `mesh`, `auto_att`, the lups |
| Bluesky | `~/bin/start_bluesky.sh` | `dscan_ophyd`, `ascan_ophyd`, `d2scan_ophyd`, `a2scan_ophyd`, `dmesh_ophyd`, `mesh_ophyd`, `auto_att_ophyd` — **no lups**, see [The lups](#the-lups) |

Both names are the **same function object**; the aliases are assigned at the
bottom of [ophyd_scan.py](ophyd_scan.py).

The reason for the suffix is that `startup.py` also does
`from .plans.align.scan_8id import *`, and Sam's module defines `dscan`,
`ascan`, `d2scan`, `a2scan`, `dmesh` and `mesh` too — as Bluesky **generators**,
which is what `RE()` expects. If the Ophyd versions won those names in that
session, then

```python
RE(dscan(huber.eta, -0.5, 0.5, 40, 0.1, lambda2M, att_ratio=10))
```

would not fail cleanly. Python evaluates arguments before it calls anything, so
the **entire Ophyd scan would run** — motor moved, shutter opened, sample
exposed, `.h5` written, a measurement number burned — and only then would `RE()`
raise, because it got a `ScanCsv` instead of a generator. The user reads
"error", runs it again, and double-exposes the sample.

`ophyd_scan.py`'s `__all__` therefore exports only the suffixed names, so a
`import *` cannot shadow Sam's. The plain names stay available to an explicit
`from ... import dscan` — which is what `startup_ophyd.py` does, and that session
never imports `scan_8id` at all.

The `.csv` records the **bare** name either way: `scan_type` is `dscan`, and
`command` reads `dscan(...)`, never `dscan_ophyd(...)`.

## Annotating a scan

```python
dscan(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M,
      comment="3x3 grid spot 5, after realigning KB")
```

lands near the top of the file, where it is easy to find:

```
#start_time,2026-09-07 15:21:49
#h5_file,/gdata/dm/8ID/8IDE/2026-3/comm202609/data/bluesky/A0113_Test_a1010041.h5
#h5_name,A0113_Test_a1010041.h5
#comment,"3x3 grid spot 5, after realigning KB"
```

Commas and quotes inside a comment are safe — the writer quotes the field the
way any CSV reader expects. Newlines are flattened to spaces so a pasted
paragraph cannot break the file's structure. Leave `comment` out and the line is
omitted entirely.

Because it is plain text, `grep -l "spot 5" *.csv` in the data folder finds the
scan.

## What it does to the hardware

Before the loop: `pre_align()`, `att(att_ratio)`, `PIND_status(0)`, then the
detector is armed. During the loop: move, trigger, wait `count_time`, read.
After: acquisition stops, the image file closes, the beam is blocked **and the
close is confirmed against the shutter readback**, trigger modes are restored,
`shutteroff()`, and every motor drives back to start.

The order is fixed, and one rule sets it: everything before the verified shutter
close must be a non-blocking `.put()`; everything after it is dose-free and may
take as long as it needs. That is why the motors go home *after* the beam is
off, never during.

## Interrupting with Ctrl+C

```
  12        30.150         0.501       1234567
^C
huber_delta: back at 31.001 (start 31.001)
^C  Scan aborted at point 12/41.
Scan file closed (aborted, 12 points): /gdata/.../A0113_Test_a1010041.csv
# images captured:  12
```

The hardware cleanup runs in the scan's `finally:`, so its lines come out
*before* the "Scan aborted" message, not after it.

What happens, in order:

1. Acquisition halts and the image file closes. On the eiger branch the HDF
   drain waits only for the frames actually **triggered**, not for `num_pts`, so
   an abort does not sit out the full `num_pts * count_time + 10` s timeout.
2. `blockbeam()`, then the blade is polled against `state_rbv` until it reads
   closed. If it does not confirm within 2 s the scan escalates once
   (`shutteroff()` to drop softglue's override, then `blockbeam()` again), and if
   it still does not confirm it prints
   `*** SHUTTER DID NOT CONFIRM CLOSED -- CLOSE IT BY HAND ***` in red. The
   request alone is not trusted: `operation` and `state_rbv` disagree while
   softglue holds the override, which is every external-trigger acquisition, and
   a pyepics put to a dropped IOC is a silent no-op.
3. Every motor is stopped, then **settled** — the scan waits for `DMOV` to go
   back to 1, i.e. for the deceleration ramp to finish, before commanding the
   return.
4. All motors are sent home **at once** (`move(wait=False)`, then wait on each
   status), so a two-motor scan returns along the same kind of path it scanned
   and one dead axis cannot hold up the others. Each axis prints where it ended
   up and where it started.
5. The CSV is closed with a bare `#END`, then `#exit_status,aborted` and
   `#points_written,12`, so the partial scan is still a valid file.

The `KeyboardInterrupt` is caught, so you get the message above rather than a
traceback, and the scan still returns its `ScanCsv`.

### Why the settle, and why a sleep is not enough

The return move used to be `stop()`, `sleep(0.2)`, `move(start)`. Measured on
`huber.nu` (`8ideSoft:CR8-E1:m4`, VELO 0.4 deg/s) on 2026-09-06, mid-move:

| cleanup | `DMOV` when the move was issued | result |
|---|---|---|
| `stop(); sleep(0.2); move(home)` | 0 | motor **never moved** — and `status.wait()` returned `success=True` |
| `stop(); sleep(1.0); move(home)` | 1 | motor went home |

Two things go wrong at once when the ramp has not finished: the record will not
honour a new target, *and* ophyd's `MoveStatus` completes on the `DMOV` 0→1 edge
produced by the **deceleration** rather than by the new move, so `.wait()`
reports success for a move that never happened. A cleanup that prints "back at
start" over a motor parked mid-scan is exactly the failure this module exists to
prevent, so it waits for the real thing (up to `SETTLE_TIMEOUT`, 10 s) instead of
guessing a sleep that is right for one axis and wrong for the next.

### Repeated Ctrl+C

**A second Ctrl+C during cleanup is absorbed, not obeyed.** SIGINT is deferred
for the whole of the teardown: the handler counts the press, writes one line to
stderr with `os.write` (never `print`, which would deadlock on the buffered
writer it interrupted), and returns. During the scan proper the handler *is*
`signal.default_int_handler`, so the first ^C raises exactly as it always has.

This is not defensive decoration. `EpicsMotor.move()` puts the setpoint one line
**before** the `try` that catches `KeyboardInterrupt`, and re-raises out of that
`try` on the way past — so the slowest, last step of cleanup is a call ophyd
deliberately aborts out of, skipping everything after it. On the eiger branch
that used to leave the shutter open. No arrangement of `try`/`finally` fixes it,
because the cleanup code is itself interruptible.

**The escape hatch needs three presses *and* five seconds.** Count alone is
wrong — a held key repeats at ~30 Hz and reaches three presses in ~100 ms, so a
key bounce would abandon a cleanup that was about to succeed. Time alone is
wrong — one stray press should arm nothing. When it does fire it abandons only
the **return move**, and it prints where each motor was left:

```
^C^C^C  leaving huber_nu at 0.4998 (start was 0.1999). The beam is off. Move it back at the prompt.
```

Nothing else is skippable: the beam is already off by the time the hatch can arm,
and every other cleanup step is a fast `.put()`.

Two limits worth knowing:

- The guard only works on the main thread. `signal.signal()` raises `ValueError`
  anywhere else, and the scan then runs with no interrupt guard — which is
  harmless, because CPython only delivers ^C to the main thread anyway.
- Each individual cleanup step is wrapped so a failure is printed and the
  remaining steps still run. A warning line naming a step is not a reason to
  stop reading; the lines after it say what else happened.

## Two-motor and raster scans

A two-motor scan writes **one column per motor**, in the order the scan was
given them, both before `elapsed_time`:

```
,huber_nu,huber_delta,elapsed_time,tetramm1_current1,...,lambda2M_stats1_total
```

That is a template entry (`{motor2}`, `source: motor2`, `optional: true`), and
it is checked before the scan starts: **a scan refuses to run if the template
has no column for an axis it is about to move**, because the resulting file
could not be plotted against that axis and nothing downstream could reconstruct
where it was.

A **raster** (`dmesh`/`mesh`) additionally writes the grid, which a trajectory
(`d2scan`/`a2scan`) does not:

```
#scan_type,dmesh
#motor,huber_nu
#motor2,huber_delta
#num_points,9
#shape,3x3
#num1,3
#num2,3
#motor1_start,0.1999
#motor1_stop,0.7999
#motor2_start,29.702
#motor2_stop,30.302
```

- **`shape` is the flag.** A viewer should test for the *presence* of that line
  to decide grid-versus-scatter, never string-match `scan_type`.
- **`motorN_start`/`motorN_stop` are the COMMANDED first and last position** of
  each axis, written before anything moves. A viewer needs them because the
  readbacks are not enough: a real stage repeats only to its retry deadband, so
  a 5×5 raster has 25 distinct readbacks and not 5 per axis, and an aborted mesh
  spans less than it was told to. Reconstructing the grid from the measured span
  puts the measured rows in the *wrong cells*, which looks plausible and is
  wrong.
- `num_points` stays the count the scan was **asked** for (`num1 * num2`). What
  it achieved is `points_written`, after the `#END`.

### Viewing one

```bash
start_scanviewer.sh              # 1D  (dscan, ascan, d2scan, a2scan, the lups)
start_scanviewer.sh --mesh       # 2D  (dmesh, mesh)
start_scanviewer.sh --dir <folder>
```

`~/bin/start_scanviewer.sh` launches BLUETELLA (`~/Documents/BLUETELLA_9ID`) on
this experiment's `data/bluesky` folder, and adds `--scan-prefix A` unless you
pass your own — 8-ID names scans `A####`, BLUETELLA names its own `S#####`.

That command is a thin wrapper around
[`scripts/start_scanviewer.sh`](../../../../scripts/start_scanviewer.sh) in this
repo, which is the copy to edit; the one in `~/bin` only `exec`s it.

`--mesh` selects `meshviewer.py`: X vs Y coloured by Z. A file carrying a
declared `shape` is drawn as a filled image, binned by nearest **commanded**
position; anything without one stays a coloured scatter. Cells that were never
measured stay blank, not zero. There is a "Draw raster as grid" checkbox if you
want the dots back. Without `--mesh` you get `scanviewer.py`, the 1-D viewer,
which is also the right one for `d2scan`/`a2scan`.

> **Under heavy attenuation the area-detector columns read zero.** At
> `att_ratio=1e6` the `lambda2M_stats*_total` columns are all 0, so a viewer
> defaulting Z (or Y) to the last column shows a blank map. Pick
> `tetramm1_sum_all` instead. This is the first thing that confuses people.

---

# Defining the CSV file structure

## Anatomy of the file

```
#start_time,2026-09-07 15:21:49                       <- lines:
#h5_file,/gdata/.../A0113_Test_a1010041.h5            <- lines:
#h5_name,A0113_Test_a1010041.h5                       <- lines:
#comment,"3x3 raster with images"                     <- lines:  (skipped if empty)
#command,"dmesh(huber_nu, -0.3, 0.3, 3, huber_delta, -0.3, 0.3, 3, 0.1, det=lambda2M)"
#scan_type,dmesh
#motor,huber_nu
#motor2,huber_delta                                   <- lines:  (two-motor scans only)
#detector,lambda2M
#num_points,9
#shape,3x3                                            <- lines:  (rasters only)
#num1,3
#num2,3
#motor1_start,0.1999
#motor1_stop,0.7999
#motor2_start,29.702
#motor2_stop,30.302
#count_time,0.1
#beamline,8-ID-E
#epoch,1789137709.041
#attenuation,1010041                                  <- lines:  (source: a PV)
#huber_nu,0.1999                                      <- lines:  (source: a PV)
#huber_delta,29.702
...
#roi1,730,983,100,10                                  <- lines:  (values: 4 PVs)
#DATA                                                 <- marker:
,huber_nu,huber_delta,elapsed_time,...,lambda2M_stats1_total   <- columns: NAMES
0,0.1999,29.702,2.52,...,0                            <- columns: one row per point
1,0.1999,30.002,3.61,...,0
...
#END                                                  <- automatic
#exit_status,success                                  <- automatic
#points_written,9                                     <- automatic
```

Two parts, two template keys:

- **`lines:`** — things that *do not* change during the scan. Written once, above
  the marker, as `#label,value`.
- **`columns:`** — things that *do* change. One column each, one row per point,
  below the marker.

Every header label carries a leading **`#`**. Nothing is dropped and nothing
moves to a sidecar — all ~40 template lines are in this one file, in template
order.

**This is BLUETELLA's (9-ID's) own extended CSV format**, adopted on 2026-09-07
so the 9-ID viewers can open an 8-ID scan with no 8-ID-specific code in them.
Two things changed to get there, and only two: the header labels gained the
`#`, and `#END,<status>,<n>` became a bare `#END` with the outcome after it.

The leading column — 0-based point index, **blank name** — is not a template
entry either. Every scan file gets it, and the blank name is the pandas index
convention, so `read_csv(..., index_col=0)` treats it as the index rather than
as a data column. That is also why it does not become the default X axis: the
scanned motor is still the first *named* column.

### The terminator

The three closing lines are written by the scan, not the template.

| line | meaning |
|---|---|
| `#END` | end of the table. **Bare** — no fields |
| `#exit_status,<status>` | `success`, `aborted`, or `error` |
| `#points_written,<n>` | rows actually written |

The bare `#END` is not cosmetic. `#END` is the last line a reader parses as part
of the table, so the extra fields of the older `#END,success,41` form were read
as one more data point: a phantom row whose motor column held the string
`success`, which flipped that column to text and broke centre-of-mass and the
derivative view in any reader that matches the marker by whole line (BLUETELLA's
does). A real 21-point file read back as 22 rows.

`points_written`, not `num_points`, for a related reason: the header already
carries `num_points`, the count the scan was **asked** for, and an aborted scan
must not silently overwrite it with the smaller number it achieved.

The absence of `#END` means the scan is still running — which is how a live
viewer knows to keep polling.

## What you see on the screen

The same columns are echoed to the terminal as the scan runs, but **only as a
summary** — it is meant to be read at a glance while the beam is on:

```
Scan file: /gdata/dm/8ID/8IDE/2026-3/comm202609/data/bluesky/A0113_Test_a1010041.csv
   #   huber_delta  elapsed_time .stats4_total .stats3_total .stats2_total .stats1_total
   1        29.501         1.050         0.000         0.000         0.000       1234567
   2        29.750         2.101         0.000         0.000         0.000       1234567
Scan file closed (success, 5 points): /gdata/.../A0113_Test_a1010041.csv
```

- Numbers are shown to **3 decimals**; integers are left alone. The `.csv` still
  gets every digit — `29.500603937599992`, not `29.501`.
- A column name longer than 13 characters is trimmed **from the left**, with a
  leading `.` to say so: `lambda2M_stats4_total` shows as `.stats4_total`. The
  end is the part that tells two columns apart, and the file keeps the full name.
- The first column is the point number, so you can see how far along the scan is.

Both settings live at the top of `scan_csv.py` as `SCREEN_DECIMALS` and
`SCREEN_WIDTH`. Neither touches the file. Pass `verbose=False` to `open_scan()`
to print nothing at all.

## The three ways to give a line a value

```yaml
lines:
  - {label: beamline,    value: "8-ID-E"}                    # #beamline,8-ID-E
  - {label: attenuation, source: filter_8ide.attenuation.readback}   # #attenuation,10
  - {label: roi1, values: [det.roi1.min_xyz.min_x, det.roi1.min_xyz.min_y,
                           det.roi1.size.x, det.roi1.size.y]}        # #roi1,730,983,100,10
```

| key | takes | produces |
|---|---|---|
| `value:` | text, with `{...}` substitutions | `#label,text` |
| `source:` | **one** dotted path to a PV | `#label,value` |
| `values:` | a **list** of dotted paths | `#label,v1,v2,v3` |

All are read once, at scan start. The `#` is added by the writer — do **not**
put one in the template's `label:`, or the file gets `##beamline`.

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
| `{start_time}` | `2026-09-07 15:21:49` |
| `{epoch}` | the same instant as a unix timestamp |
| `{h5_file}` | **full path** of this scan's detector `.h5` (see below); empty when `save_img=0` |
| `{h5_name}` | just the **name**, `A0113_Test_a1010041.h5`; empty when `save_img=0` |
| `{comment}` | the `comment="..."` argument |
| `{command}` | `dscan(huber_delta, -0.5, 0.5, 41, 1.0, det=lambda2M)` |
| `{scan_type}` | `dscan` — the bare name, never the `_ophyd` alias |
| `{motor}` | `huber_delta` |
| `{motor1}` … `{motor4}`, `{motors}` | one per scanned motor, blank past the ones this scan has; `{motor}` is an alias for `{motor1}` |
| `{num1}`, `{num2}`, `{shape}` | the raster grid, from `shape=(num1, num2)`; blank on a scan that is not a raster |
| `{motor1_start}`, `{motor1_stop}`, `{motor2_start}`, `{motor2_stop}` | the commanded ends of each raster axis; blank on a scan that is not a raster |
| `{det}` | `lambda2M` |
| `{num_points}`, `{count_time}` | as passed |

`{h5_file}` and `{h5_name}` are both there because they answer different
questions: the path says where the images went, the name is the stem the `.csv`
shares with them and is what you type when looking a measurement up.

### Where `{h5_file}` comes from

It is **not** an argument to the scan. The scan builds it:

1. `gen_folder_prefix()` (defined in `plans/acquire/acq_helpers.py`, and
   re-exported by `ad_acq.py`) combines `header` + `expt.measurement_num` +
   `sample_name` + the current attenuation into
   `{header}{measurement_num:04d}_{sample_name}_a{attenuation:04d}`, e.g.
   `A0113_Test_a1010041`, then advances `measurement_num`.
2. That prefix names both files, in the folder built from `expt.mount_point` +
   `cycle_name` + `experiment_name` + `/data/bluesky`.

So the CSV and the `.h5` always share a name, and each scan gets exactly one
measurement number.

**Where the name comes from on a standalone scan.** `header` and `sample_name`
are read from run state (`expt.header` / `expt.sample_name`) when a measurement
is in flight — during a multi-sample run that is the only thing that knows which
sample is actually in the beam. With no measurement running, which is the normal
case for an alignment scan, they fall back to the
`sample_{expt.sample_index}` block of `sample_info.yaml` — the same source
`run_measurement()` itself reads, and `sample_index` is what `select_sample(<n>)`
sets. Before that fallback existed an align scan on a fresh session died with
`'header' has not been set yet`, *after* `pre_align()` and `att()` had already
changed beamline state, and the workaround was to set the two by hand, which
invited invented names that then live in file names for ever.

If neither source works the scan says so before moving anything, and names both
ways out: pick a sample with `select_sample(<n>)`, or set `expt.header` /
`expt.sample_name` yourself.

**That counter is shared with acquisitions.** `measurement_num` is the EPICS
register `8ideSoft:Reg1`, and `det_acq_series()` calls the same
`gen_folder_prefix()` — so a scan advances the acquisition numbering and an
acquisition advances the scans'. The two do not share a *folder*, though:
acquisitions write under `data/`, these scans under `data/bluesky/`, so the
highest number already used may be under either one. See
[docs/configuration.md](../../../../docs/configuration.md#the-measurement-counter-stays-in-epics)
for why the counter is kept in EPICS rather than in `state/run_state.yml`.

## Per-entry flags

```yaml
  - {label: comment,     value: "{comment}", skip_if_empty: true}
  - {label: rheometer_x, source: rheometer.x, optional: true}
```

- **`optional: true`** — if the device or signal does not exist, skip the entry
  and print one note. Use it for anything that is not on every station or not on
  every detector. Every *detector* column in the shipped template is optional,
  which is how one template serves eiger, lambda and tetramm; the first motor
  column and `elapsed_time` are deliberately not.
- **`skip_if_empty: true`** — leave the line out when the value comes back blank.

An entry **not** marked `optional` whose source is missing is an error raised
**before the scan starts** — not halfway through with the beam on.

## Columns

```yaml
columns:
  - {label: "{motor1}",           source: motor1}
  - {label: "{motor2}",           source: motor2, optional: true}
  - {label: elapsed_time,         source: elapsed}
  - {label: tetramm1_current1,    source: tetramm1.current1.mean_value, optional: true}
  - {label: "{det}_stats1_total", source: det.stats1.total, optional: true}
```

These sources are computed by the scan rather than read from a PV:

| source | value |
|---|---|
| `motor1` … `motor4` | where that scanned motor actually is (readback) |
| `motor1_setpoint` … `motor4_setpoint` | where it was told to go |
| `elapsed` | seconds since the scan started |
| `epoch` | unix timestamp of the point |

A bare `motor` / `motor_setpoint` means the first one, so a beamline-local
template written before two-motor scans existed still means what it meant.
`MAX_MOTORS` is 4, which keeps `{motor5}` an honest error rather than a silently
blank column.

Anything else is a dotted path, read fresh at every point.

**Order matters for plotting.** `scanviewer.py` defaults X to the first column;
`meshviewer.py` defaults X to the first, Y to the second and Z to the **last**.
So: the scanned motors first, in the order the scan moved them, and your main
counter last — which is why the shipped template ends with `stats1`, and why
`tetramm1_sum_all` sits after `tetramm3_sum_all` (on a scan whose detector *is* a
tetramm, every `det.stats*` column is skipped and the last column standing has to
be the picoammeter being scanned).

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
    A0113_Test_a1010041.h5      images
    A0113_Test_a1010041.csv     this scan
    A0114_Test_a1010041.csv     next scan; no .h5 if it ran with save_img=0
```

built from `expt.mount_point`, `expt.cycle_name` and `expt.experiment_name` —
three settings in `configs/experiment.yml` — the same folder `save_images()` in
`scan_8id.py` writes to. The folder is created, if needed, **before** the scan
moves anything; if that fails the error names those three settings, because a
`mount_point` left pointing at the previous experiment's station tree is the
usual cause and the bare `PermissionError` from inside `os.makedirs` says
nothing useful.

### Finding the running scan from anywhere

The name shared by both files is recorded as the scan starts, right after the
header is written:

```python
expt.file_name          # 'A0113_Test_a1010041'
```

`file_name` is persistent session state, so it is on disk as well as in the
session — in the `persistent:` block of `state/run_state.yml` in the checkout
the session is running from. A GUI or a shell script running **as the account
that owns the session** can read it there with no channel access at all:

```bash
grep file_name ~/bluesky/state/run_state.yml
```

The file is written mode `0600`, so another account cannot read it — unlike the
registers below, which anything with channel access could `caget`.

It is the **bare name — no extension, no path**. Append `.h5` or `.csv`, and
rebuild the folder from `mount_point` + `cycle_name` + `experiment_name` +
`/data/bluesky` (the same three `configs/experiment.yml` settings the scan
itself uses).

The name is written whether or not `save_img=1` — the `.csv` is always
produced, and only the `.h5` is conditional.

`det_acq_series()` in `plans/acquire/ad_acq.py` sets this **same** field with
its measurement name, so one look finds whichever of the two is running. That
also means it holds the most recent of either — a scan started after an
acquisition overwrites it.

> **Changed 2026-09-06.** This used to be published to EPICS. First as two
> registers holding full paths, `scan_h5_file` (`StrReg21`) and `scan_csv_file`
> (`StrReg22`) — both commented out in `registers_device.py` on 2026-09-02 —
> and then as the bare name in `file_name` (`StrReg8`). Nothing writes
> `StrReg8` any more: the Component is still declared on
> `EpicsPvStorageRegisters`, so a `caget` on it returns whatever was left there
> before the move, which is worse than nothing. Anything that polled those
> registers — including the specr_py viewer below — needs to read
> `state/run_state.yml` instead. `measurement_num` (`8ideSoft:Reg1`) is the one
> register that did *not* move; see
> [docs/configuration.md](../../../../docs/configuration.md#the-measurement-counter-stays-in-epics).

The specr_py viewer (`~/Documents/specr_py`, `~/bin/start_specr_py.sh`) does
this at startup, and again on File ▸ Follow Live Scan: read the published name,
open what it names, start monitoring. That lookup is the part the note above
affects — it still reads `8ideSoft:StrReg8`. It is a convenience only: the
viewer's actual live update comes from polling the data folder, so it works on
machines with no channel access to `8ideSoft:`.

BLUETELLA (`start_scanviewer.sh`, see
[Viewing one](#viewing-one)) reads no registers and no run state at all — it
lists and polls the folder — so the note above does not affect it.

## Reading the file back

With pandas — skip everything down to and including the marker:

```python
import pandas

def read_scan(path, marker="#DATA"):
    with open(path) as f:
        skip = next(i for i, line in enumerate(f) if line.startswith(marker)) + 1
    return pandas.read_csv(path, skiprows=skip, comment="#", index_col=0)

df = read_scan("/gdata/.../A0113_Test_a1010041.csv")
df.plot(x=df.columns[0], y=df.columns[-1])
```

`index_col=0` consumes the blank-named point-index column, so `df.columns[0]`
is the scanned motor and the plot above is unchanged. Drop `index_col=0` on a
file written before that column existed.

`comment="#"` is what drops the `#END` / `#exit_status` / `#points_written`
lines. Do not remove it, and do not read a file whose terminator is the old
`#END,success,41` form without it — those extra fields parse as one more data
point.

Without pandas, `scan_csv.py` ships a reader that imports only the stdlib:

```python
from id8_common.plans.align.scan_csv import read_scan_csv

header, labels, columns, status = read_scan_csv(path)
header["h5_file"]        # ['/gdata/.../A0113_Test_a1010041.h5']
header["shape"]          # ['3x3']   -- absent on a non-raster
header["roi1"]           # ['730', '983', '100', '10']
header["points_written"] # ['9']     -- absent while the scan is running
columns["huber_nu"]      # [0.1999, 0.1999, ...]
status                   # 'success' | 'aborted' | 'error' | 'running' | 'unknown'
```

Header keys are stored with the `#` stripped, so the same code reads files
written before 2026-09-07, when the labels carried none. It strips the index
column too, so `labels` and `columns` are the same whether or not the file has
one.

`status == "running"` means no `#END` line yet — call it again in a second and
you have a live view. `"unknown"` means there is an `#END` but no
`#exit_status` after it, which is a file truncated between the two. A
half-written trailing row is ignored until it is complete, so polling a scan in
progress is safe.

## Adding another scan

Copy the closest existing scan, change the loop, and reuse the same four lines:

```python
scan = scan_csv.open_scan(
    csv_file, motor, det=det, scan_type="d3scan", command=command,
    num_points=num_pts, count_time=count_time, h5_file=h5_file, comment=comment,
)
scan.write_header()
# per point:
scan.add_point(setpoint)
# at the end (in a finally:):
scan.close(status)
```

Multi-motor scans need no special handling: pass a **list** of positioners as
`motor`, and one setpoint per motor to `add_point(p1, p2)`. `{motor1}`,
`{motor2}` and the `motor1`/`motor2` column sources then resolve on their own —
that is what `d2scan`, `a2scan`, `dmesh` and `mesh` already do. Up to
`scan_csv.MAX_MOTORS` (4) of them. A raster additionally passes
`shape=(num1, num2)`, which is how a viewer tells a grid from a trajectory —
`d2scan` sweeps both axes along one line and deliberately passes no `shape` — and
the commanded ends of each axis through `extra=`:

```python
    shape=(num1, num2),
    extra={
        "motor1_start": repr(float(positions1[0])),
        "motor1_stop":  repr(float(positions1[-1])),
        "motor2_start": repr(float(positions2[0])),
        "motor2_stop":  repr(float(positions2[-1])),
    },
```

Anything in `extra=` becomes a `{name}` the template can substitute, so a new
scan can publish its own metadata without a change to `scan_csv.py`.

A template that has no column for a motor the scan actually moved is rejected
before the scan starts, so a new scan cannot quietly write an unplottable file.

Wrap the loop in `_scan_guard()` and end the cleanup with `_blockbeam_verified()`
and `_return_motors([...])` rather than open-coding a `finally:` — that is what
makes a new scan behave like the six under Ctrl+C.

## Troubleshooting

| symptom | cause |
|---|---|
| `no device 'foo' in the oregistry` before the scan starts | typo in a template `source:`, or the device is not loaded. Add `optional: true` if it is genuinely not always present. |
| `CSV template has duplicate column labels` | two entries render to the same text once `{det}`/`{motor}` is filled in. Rename one. |
| `unknown substitution '{...}'` | a `{name}` that is not in the list above. The error prints the available ones. |
| `WARNING: column X is not a number` | the dotted path points at a device rather than a signal. Name the individual signal. |
| CSV written but no `.h5` | `save_img=0` — the scan says so before it starts (`save_img=0: no .h5 will be written; scan file is ...`) — or the HDF plugin is not enabled. With `save_img=1` the scan prints `Scan folder created:` and the path it armed. |
| no `#END` line | the session was killed outright. Everything up to the last complete row is still readable. |
| `FileExistsError: ... refusing to append a second scan to it` | the measurement counter has repeated, so the `.h5` would collide too. Check `expt.measurement_num` (`caget 8ideSoft:Reg1`) against **both** `data/` and `data/bluesky/`; move or delete the old file if it is not wanted. |
| `Cannot work out the scan name` | no measurement running and `sample_{expt.sample_index}` in `sample_info.yaml` could not be read, or has no `header`/`sample_name`. Run `select_sample(<n>)`, or set `expt.header` / `expt.sample_name`. |
| `Cannot use the scan data folder ...` | one of `mount_point` / `cycle_name` / `experiment_name` in `configs/experiment.yml` is stale. The error names all three; a `mount_point` left on the other station's tree (`/gdata/dm/8ID/8IDI/` vs `8IDE/`) is the usual one. |
| `*** SHUTTER DID NOT CONFIRM CLOSED ***` | the blade did not read closed within 2 s, even after one `shutteroff()` + `blockbeam()` escalation. **Close it by hand**, then look at the shutter IOC. |
| `WARNING: <motor> still reports DMOV=0 after 10s` | the axis never finished decelerating, so the return move may not have taken. Check where it actually is before the next scan. |
| a scan left a motor away from its start | the ^C escape hatch fired (3 presses over 5 s). It prints `^C^C^C leaving <motor> at ...`; the beam is off, move the axis back at the prompt. |
| map is blank in the 2-D viewer | the `{det}_stats*` columns are all zero under heavy attenuation. Set Z to `tetramm1_sum_all`. |

## Notes and limits

- **No RunEngine.** Nothing here writes to Tiled or databroker, and none of the
  data is in the catalog. The CSV *is* the record.
- **Not queueserver-safe.** These are plain functions, not plans, so the QS
  cannot run them. Use them interactively. The interrupt guard also needs the
  main thread; on a worker it installs nothing and says nothing, which is
  harmless because ^C is never delivered there anyway.
- **Every scan burns a measurement number**, `save_img=0` included, so the `.csv`
  always has a unique name. That counter is shared with `det_acq_series()`.
- The template is re-read at the start of every scan, and the detectors and
  template devices are looked up per scan, so a device that was offline at
  startup fails at the call that needs it rather than poisoning the import.
  `ophyd_scan` still binds `softglue` and `softglue_8id_acq` at module scope,
  so it must be imported after `make_devices()` — as `startup.py` does.
- **`scan_8id.py` is untouched** and still works. The Bluesky session imports
  both; Sam's generators keep the plain names there. Everything `ophyd_scan.py`
  does differently is marked `# CHANGED (n):` in its source, with the reason.

## Verified on hardware

All six scans were run on pearl on 2026-09-07 against `lambda2M` and two real
motors — `huber.nu` (limits ±1.1) and `huber.delta` (kept inside 29–31) — in
`/gdata/dm/8ID/8IDE/2026-3/comm202609/data/bluesky/`:

| file | scan | points | shape |
|---|---|---|---|
| `A0107_Test_a1010041.csv` | `dscan` (huber.nu) | 5 | 1-D |
| `A0108_Test_a1010041.csv` | `ascan` (huber.delta) | 5 | 1-D |
| `A0109_Test_a1010041.csv` | `d2scan` | 5 | trajectory, no shape → scatter |
| `A0110_Test_a1010041.csv` | `a2scan` | 5 | trajectory, no shape → scatter |
| `A0111_Test_a1010041.csv` | `dmesh` | 25 | 5×5 raster → grid |
| `A0112_Test_a1010041.csv` | `mesh` | 16 | 4×4 raster → grid |
| `A0113_Test_a1010041.csv` | `dmesh`, `save_img=1` | 9 | 3×3 raster → grid, plus `.h5` |

Ctrl+C was exercised on all three detector branches.
