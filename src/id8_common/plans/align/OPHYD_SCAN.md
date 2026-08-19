# Ophyd scans and the SPEC file template

How to run `dscan_ophyd()`, and how to control exactly what ends up in the
`.spec` file it writes — without editing any Python.

## What is here

| file | what it is | imports hardware? |
|---|---|---|
| `ophyd_scan.py` | the scans themselves (`dscan_ophyd`), detector plumbing, interrupt-safe teardown | yes |
| `ophyd_spec_config.py` | reads the YAML template and resolves it against the live `oregistry` | yes |
| `ophyd_spec_writer.py` | formats and appends SPEC lines; `SpecFile` | **no** — stdlib only |
| `/home/beams/8IDIUSER/bluesky/src/id8_common/configs/spec_template.yml` | **the file you edit** to change the SPEC layout | — |

`ophyd_spec_writer.py` is deliberately free of ophyd/apsbits/bluesky imports so it
can be exercised offline and reused by any future scan.

**No Bluesky anywhere.** No RunEngine, no plans, no generators, no documents —
just `.put()` / `.get()` / `.move()`. Nothing lands in databroker or Tiled. If you
need the run in the catalog, use `dscan()` in `scan_8id.py` instead.

---

## Quick start

```python
# in the beamline session (start_bluesky.sh)
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M, att_ratio=10)
```

Call it **directly** — do *not* wrap it in `RE()`. It is an ordinary function.

Watch it live from a workstation with the `specr_py` viewer:

```bash
conda activate specr_py
python ~/Documents/specr_py/specr.py /gdata/dm/8ID/8IDE/2026-2/pope202607/data/bluesky/<spec_file_name>.spec --watch
```

This document is the reference for the **scan and the file format**. For the
viewer — install, live mode, overlays, reading `.spec` files from your own
scripts — see `~/Documents/specr_py/README.md` and its `OPHYD_SCAN_GUIDE.md`
(<https://github.com/qzhang234/specr_py>).

---

## `dscan_ophyd()`

```python
dscan_ophyd(motor, rel_begin, rel_end, num_pts, count_time,
            det=None, comment="", att_ratio=1e6, save_img=1, spec_path=None,
            spec_template=None, set_attenuation=True, beam_control=True,
            verbose=True)
```

| argument | meaning |
|---|---|
| `motor` | **any** Ophyd positioner — `huber.delta`, `sample.x`, `rheometer.y`, ... |
| `rel_begin`, `rel_end` | start/end **relative to the current position** |
| `num_pts` | number of points, inclusive of both ends |
| `count_time` | detector exposure per point, seconds |
| `det` | `eiger4M` (default) or `lambda2M` |
| `comment` | free-text note, written as the first `#MD` line — see [below](#annotating-a-scan) |
| `att_ratio` | attenuation ratio, same meaning as in `dscan()` |
| `save_img` | `1` writes images and burns a measurement number; `0` scans without saving |
| `spec_path` | override the SPEC file location |
| `spec_template` | override the template for this one scan |
| `set_attenuation` | `False` leaves the attenuator exactly as-is |
| `beam_control` | `False` skips `pre_align()`, `PIND_status()` and the shutter |
| `verbose` | `False` suppresses the per-point table |

Nothing in the scan is specific to the diffractometer — any object with `.name`,
`.position` and `.move()` works:

```python
dscan_ophyd(sample.x,     -0.1,  0.1, 21, 0.5, det=lambda2M)
dscan_ophyd(rheometer.y,  -0.05, 0.05, 11, 2.0, det=eiger4M)
dscan_ophyd(huber.eta,    -0.2,  0.2, 21, 0.5, det=lambda2M, save_img=0)
```

Returns a `ScanResult`, which unpacks as `(positions, columns)` but prints as one
line instead of dumping every value:

```python
pos, cols = dscan_ophyd(huber.eta, -0.2, 0.2, 21, 0.5, det=lambda2M)
import numpy as np
print("peak at", pos[np.argmax(cols["lambda2M_stats1_total"])])
```

`positions` is the array of **measured** readbacks, not the commanded ones.

### Annotating a scan

Pass `comment=` to record why you ran a scan. It is written as the **first**
`#MD` line of the scan block, so you can `Ctrl+F` the `.spec` file for a keyword
and land straight on the scan you want:

```python
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M,
            comment="3x3 grid spot 5, after realigning KB")
```

```
#S 196  dscan_ophyd(huber_delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
#D Tue Aug 18 17:31:23 2026
#C ... plan_type = function
#MD comment = 3x3 grid spot 5, after realigning KB      <- first, easy to find
#MD beamline_id = 8-ID-E
...
```

Letters, digits, spaces and punctuation are all fine, including `=`, `#`, quotes
and non-ASCII — a reader splits a `#MD` line on the *first* `=` only, so the rest
of your text survives untouched. Newlines and tabs are flattened to single
spaces, because a line break mid-record would end the `#MD` line early and leave
the remainder looking like a data row. Leave `comment` empty (the default) and no
`#MD comment` line is written at all.

It is an ordinary template entry — the first one in `metadata:`:

```yaml
- {key: comment, value: "{comment}", skip_if_empty: true}
```

so you can move it, rename the key, or delete it like anything else. If you
delete it but still pass `comment=`, the note is written first anyway rather than
being silently dropped.

> Not to be confused with the template's `comments:` list, which writes
> timestamped `#C` lines and is the same for every scan. `comment=` is per-scan.

### What it does to the hardware

With defaults it follows `dscan()` in `scan_8id.py`: `pre_align()`,
`att(att_ratio)`, `PIND_status(0)`, arm the detector, open the shutter, step and
trigger, then close the shutter, restore the trigger mode and return the motor to
its starting position. The Eiger is pre-armed for `num_pts` software triggers; the
Lambda is externally triggered by a softglue pulse per point.

### Interrupting with Ctrl+C

Safe at any point. The scan stops the motor, halts acquisition, closes the image
file, blocks the beam, restores the trigger configuration, drives the motor **back
to where it started**, and closes the SPEC block with `exit_status = aborted` so
the partial scan still parses and still plots.

```
^C  Scan aborted at point 12/41 -- stopping motor and detector...
    returning huber_delta to 30.0004 ...
SPEC scan 187 finished (aborted), 12 points.
```

- **A second Ctrl+C during cleanup is ignored** (`(cleanup in progress -- interrupt
  ignored)`), so the motor cannot be stranded mid-scan.
- **Aborting is fast.** A normal finish waits for the detector to flush every
  frame; an abort skips that wait, since the remaining frames are never coming.
  Measured: a 41-point 1 s scan tears down in ~2 s rather than the ~51 s the frame
  wait would take.
- The `KeyboardInterrupt` is re-raised after cleanup, so a loop over samples stops
  instead of silently continuing.

---

## Defining the SPEC file structure

Everything the file contains is declared in:

```
/home/beams/8IDIUSER/bluesky/src/id8_common/configs/spec_template.yml
```

Edit it and the next scan picks it up. No Python change, no session restart.

### Anatomy: template key → SPEC output

```
#F /gdata/.../data/bluesky/pope202607.spec            <- automatic
#O0 huber_nu  huber_delta  ... huber_y                <- positioners:  NAMES, first 8
#O1 huber_z  sample_x  ... rheometer_z                <- positioners:  NAMES, next 7
#o0 / #o1                                             <- same, mnemonics (a copy)

#S 186  dscan_ophyd(huber_delta, -0.5, 0.5, 41, 1.0, det=lambda2M)
#D Mon Aug 17 23:05:12 2026                           <- automatic
#C ... plan_type = function                           <- comments:
#MD comment = 3x3 grid spot 5                         <- metadata:  first entry, from comment=
#MD beamline_id = 8-ID-E                              <- metadata:  (fixed text)
#MD roi1 = min_x=730 min_y=983 size_x=100 size_y=10   <- metadata:  (sources: several PVs)
#P0 -0.00049  30.0004  ... 21.1945                    <- positioners:  VALUES, first 8
#P1 -4.7409  173.654  ... 0                           <- positioners:  VALUES, next 7
#N 6                                                  <- automatic (len of columns)
#L huber_delta  huber_delta_setpoint  ...             <- columns:
29.500655 29.5 -0.400127 0 0 0                        <- columns:, one row per point
#C ... exit_status = success                          <- automatic
```

The trailing digit on `#O`/`#o`/`#P` is **only line wrapping** — SPEC puts 8
entries per line, so 15 positioners give `#O0`+`#O1` and `#P0`+`#P1`. Both come
from the single `positioners:` list; see [below](#the-op-snapshot).

### Where values come from — `source:`

| `source:` | value written |
|---|---|
| `motor` | scanned motor **readback** at the point |
| `motor_setpoint` | where the motor was told to go |
| `epoch` | whole seconds since the scan started |
| `epoch_float` | the same, with fractions |
| `det.<path>` | attribute path on the detector **passed to this scan** — so one template serves eiger4M and lambda2M (`det.stats1.total`, `det.roi1.size.x`) |
| `<device>.<path>` | any device in the `oregistry` (`tetramm1.current1.mean_value`) |

Every source must name a signal that returns a **single value** — see
[the gotcha below](#gotcha-bare-device-paths-are-not-scalars).

### Text substitutions

Usable in any `label:` or `value:`:

| token | value |
|---|---|
| `{motor}` | name of the scanned motor, e.g. `huber_delta` |
| `{det}` | name of the detector, e.g. `lambda2M` |
| `{scan_type}` | the scan function, e.g. `dscan_ophyd` |
| `{num_points}`, `{count_time}` | exactly as passed to `dscan_ophyd()` |
| `{image_file}` | name of the detector HDF5 file — **not** an argument, see below |
| `{comment}` | the scan's `comment=` text (used by the first `metadata:` entry) |

`{image_file}` is **not** something you pass in. `dscan_ophyd()` builds it by
calling `gen_folder_prefix()` (`plans/acquire/ad_acq.py`), which assembles

```
<header><measurement_num:04d>_<sample_name>_a<attenuation:04d>.h5
```

from `pv_registers.header`, `pv_registers.measurement_num`,
`pv_registers.sample_name`, and the live `filter_8ide` attenuation — giving e.g.
`A0186_HEA-15GPa-3x3Grid_a0010.h5`. This is also where the SPEC scan number comes
from, which is why `#S 186` always matches `A0186_*.h5`.

`measurement_num` increments **once per scan, and only when `save_img=1`**. With
`save_img=0` there is no image file, `{image_file}` is empty, and the
`skip_if_empty: true` flag on that entry drops the `#MD` line entirely.

### Per-entry flags

| flag | effect |
|---|---|
| `optional: true` | device missing → skip with a warning instead of failing the scan |
| `skip_if_empty: true` | (metadata only) omit the line entirely when the value is blank |

### Several PVs on one header line: `sources:`

A header entry can read more than one PV and join them into a single `#MD` line.
Write `sources:` (plural) as a mapping and the names are kept:

```yaml
- {key: roi1, optional: true, sources: {min_x: det.roi1.min_xyz.min_x, min_y: det.roi1.min_xyz.min_y, size_x: det.roi1.size.x, size_y: det.roi1.size.y}}
```
```
#MD roi1 = min_x=730 min_y=983 size_x=100 size_y=10
```

Write it as a plain list instead to get bare values in order:

```yaml
- {key: roi1, sources: [det.roi1.min_xyz.min_x, det.roi1.size.x]}   ->  "730 100"
```

Nothing in the code knows what an ROI is — the template names the PVs, so the
same mechanism works for any group of related readings. Readers split a `#MD`
line only on the **first** `=`, so the `name=value` pairs survive intact.

`sources:` is **metadata only**. A data cell must be a single number, so a
`sources:` on a column is a hard error — a joined value would contain spaces and
shift every column after it.

### Worked edits

**Add a counter column.** Insert before the last entry — the last column is the
viewer's default Y axis, so the primary counter stays there:

```yaml
columns:
  # ...
  - {label: tetramm2_current1, source: tetramm2.current1.mean_value, optional: true}
  - {label: "{det}_stats1_total", source: det.stats1.total}   # keep last
```

**Add a time column.** Elapsed time is not recorded by default. Two sources are
available, both counting from the start of the scan (despite the
SPEC-conventional `Epoch` name) — they are present but commented out in the
template:

```yaml
columns:
  - {label: "{motor}", source: motor}
  - {label: Epoch, source: epoch}              # whole seconds
  - {label: Epoch_float, source: epoch_float}  # fractional seconds
  - {label: "{det}_stats1_total", source: det.stats1.total}   # keep last
```

**Drop a column.** Delete its entry. Nothing depends on any particular column
except the ordering rule below.

**Add a header line from a device** — one PV with `source:`, or several joined
onto one line with `sources:` (see [below](#several-pvs-on-one-header-line-sources)):

```yaml
metadata:
  - {key: sample_temperature, source: lakeshore1.readback_ch1, optional: true}
  - {key: slit4, optional: true, sources: {h: sl4.h.size.readback, v: sl4.v.size.readback}}
```

**Add a fixed or substituted header line:**

```yaml
metadata:
  - {key: operator, value: "night shift"}
  - {key: exposure, value: "{count_time} s x {num_points}"}
```

**Change the context snapshot** — edit `positioners:`. Safe mid-experiment; see
below.

### Rules the template must respect

All three are checked **before the scan moves anything**, so a mistake fails
immediately rather than producing a file that silently mis-parses hours later.

1. **Scanned motor first, primary counter last.** `specr_py` and the MATLAB
   `specr` both default X to the first column and Y to the last. The writer warns
   if the last column is the motor or its setpoint.
2. **No label may contain two consecutive spaces.** `#L` is split on runs of two or
   more spaces, so such a label would silently become two columns. Hard error.
3. **No duplicate labels.** Readers index columns by name and would only ever
   return the first. Hard error.

A missing device in a **non-optional** entry raises before the scan starts. A
missing device in an `optional:` entry is skipped with a warning.

### Gotcha: bare device paths are not scalars

**Always name the individual PV, never a parent device.** `Device.get()` returns a
*namedtuple* of every component, which is not a value anyone wants in a column or
a header:

```yaml
- {label: eta, source: huber.eta}                  # WRONG -> warns, writes 0
- {label: eta, source: huber.eta.user_readback}    # right -> the .RBV float
```

The same applies to the ROIs. `det.roi1` is a plugin, not a value; its geometry
lives in four separate PVs, which is how the template records them:

| template `source:` | equivalent in the session |
|---|---|
| `det.roi1.min_xyz.min_x` | `lambda2M.roi1.min_xyz.min_x.get()` |
| `det.roi1.min_xyz.min_y` | `lambda2M.roi1.min_xyz.min_y.get()` |
| `det.roi1.size.x` | `lambda2M.roi1.size.x.get()` |
| `det.roi1.size.y` | `lambda2M.roi1.size.y.get()` |

The template names each PV explicitly — there is no ROI-aware code — while
`sources:` keeps the four of them on one `#MD` line.

For a column, the writer prints `WARNING: SPEC column 'eta' is not a scalar (...)`
once per label if you get this wrong, so it is visible rather than silent.

### The `#O`/`#P` snapshot

**The `positioners:` list is the sole source of every `#O`, `#o` and `#P` line.**
It records where everything *not* being scanned was sitting, so the geometry can
be reconstructed later.

| line | contains | written |
|---|---|---|
| `#O0`, `#O1`, … | motor **names** | once per **file**, in the `#F` header |
| `#o0`, `#o1`, … | mnemonics — a copy of the names | once per **file** |
| `#P0`, `#P1`, … | their **positions** | once per **scan**, in the `#S` block |

The trailing digit is nothing but line wrapping: SPEC allows 8 entries per line,
so the 15 default positioners produce `#O0` (first 8) and `#O1` (remaining 7),
and the matching `#P0`/`#P1`. Ordering follows the `positioners:` list exactly, so
the Nth name across the `#O` lines pairs with the Nth value across the `#P` lines:

```
#O0 huber_nu  huber_delta  huber_mu  huber_eta  huber_chi  huber_phi  huber_x  huber_y
#O1 huber_z  sample_x  sample_y  sample_z  rheometer_x  rheometer_y  rheometer_z
#P0 -0.00049  30.0004  -0.00039  -0.00027  90.00002  0.00033  -0.12761  21.1945
#P1 -4.7409  173.654  24.9957  3.38223  0  0  0
```

Add a 17th positioner and you get a `#P2`; trim the list to 8 and `#O1`/`#P1`
disappear. Anything unavailable is skipped with a warning. `specr_py` reads these
back under **Show Motor Positions**.

#### This list does *not* restrict what you can scan

`positioners:` only controls the **context snapshot**. You can scan **any**
positioner in the session, listed or not:

```python
dscan_ophyd(mono.bragg, -0.01, 0.01, 21, 1.0, det=lambda2M)   # not in the list
```

The scanned motor is **appended automatically** to `#O`/`#P` when it is absent, so
it is always recorded. Nothing else is needed to support a new motor — no code
change, no template change.

One side effect is worth knowing. Appending changes the length of the positioner
list, so **alternating between a listed and an unlisted scanned motor starts a new
`#F` block on every switch** (the writer will not let new `#P` values pair with
stale `#O` names). Three such scans give three `#F` blocks. Everything still
parses and plots — but if you scan a motor regularly, add it to `positioners:` so
the list stays fixed and the file stays tidy.

Note that **`huber` and `psic` drive the same physical motors** (both on
`8ideSoft:CR8-E1:`; `huber.delta` and `psic.delta` are the same PV `m5`). Listing
both records every angle twice. The default lists `huber`, which is the superset —
`psic` has no `x`/`y`/`z`.

Changing this list mid-experiment is safe. `#O` is written once per *file* but
`#P` once per *scan*, and readers bind names at each `#S`, so a changed list would
otherwise pair new positions with stale names. The writer detects the change and
starts a fresh `#F` block instead:

```
NOTE: SPEC positioner list changed; starting a new #F block so #O names stay
matched to #P values.
```

Readers reset their motor-name table on `#F`, so both halves of the file stay
correct.

### Trying a template without touching the default

```python
dscan_ophyd(..., spec_template="~/my_template.yml")   # one scan
```
```bash
export ID8_SPEC_TEMPLATE=~/my_template.yml            # whole session
```

Resolution order: the `spec_template=` argument, then `$ID8_SPEC_TEMPLATE`, then
the default at
`/home/beams/8IDIUSER/bluesky/src/id8_common/configs/spec_template.yml`.

(In code that default is derived from the module location rather than hard-coded,
so a different checkout still finds its own template. The scan prints the
resolved absolute path if you ever need to confirm which one was used.)

### Effect on the viewer

Editing the template does not rewrite earlier scans, so one `.spec` file can hold
scans with different column layouts. Each plots normally, but `specr_py` will not
**overlay** across the boundary — it reports *"Multi-selections have to be scans of
the same type"*. That is expected, not a bug.

---

## Where files go

| | path |
|---|---|
| SPEC file | `<mount_point><cycle>/<experiment>/data/bluesky/<spec_file_name>.spec` |
| Detector images | `<mount_point><cycle>/<experiment>/data/bluesky/*.h5` |

The SPEC file sits in the **same directory as the images**, e.g.
`/gdata/dm/8ID/8IDE/2026-2/pope202607/data/bluesky/pope202607.spec`. Every part of
the path comes from `pv_registers`, so it follows the run cycle automatically.

### Switching SPEC files on demand

The file name comes from **`pv_registers.spec_file_name`**, so you change it at
any time without touching code:

```python
pv_registers.spec_file_name.put("alignment_run3")
```

The next scan writes to `alignment_run3.spec` in the same directory. `.spec` is
appended when the register does not already end in it (the viewer's file dialog
filters on that extension). If the register is empty the scan falls back to
`experiment_name` and says so.

**New files are detected automatically.** If the target file does not exist, the
scan writes the full `#F`/`#E`/`#D`/`#C`/`#O`/`#o` header block before the first
`#S`; if it already exists, the scan just appends and no second header is written.
Switching back to an earlier file appends to it as normal. The directory is
created if it is not there yet.

> **The viewer must be able to see `/gdata`.** This directory is on GPFS, which is
> not mounted on every analysis workstation — an earlier version of this module
> deliberately put the SPEC file in `~/bluesky` for exactly that reason. Run
> `specr_py` on a host that mounts `/gdata`, or pass `spec_path=` to write a copy
> somewhere shared.

The SPEC scan number comes from the image prefix, so `#S 186` lines up with
`A0186_*.h5`. `gen_folder_prefix()` increments `pv_registers.measurement_num` once
per scan, and only when `save_img=1`.

---

## Adding another scan (d2scan, mesh, ...)

Scans live in `ophyd_scan.py`; the SPEC machinery is already factored out, so a new
scan reuses it rather than copying it:

```python
rendered = ophyd_spec_config.render(
    motor, det=det, scan_type="d2scan_ophyd",
    num_points=num_pts, count_time=count_time,
    image_file=f"{folder_prefix}.h5" if folder_prefix else "",
    extra={"motor2": motor2.name},        # adds {motor2} as a substitution
)
spec = SpecFile(spec_path or default_spec_path())
spec.write_file_header(rendered.positioner_names)
spec.start_scan(scan_num, command, rendered.labels,
                metadata=rendered.metadata,
                motor_positions=rendered.positioner_positions,
                comments=rendered.comments)
# per point:
spec.add_point(rendered.read(motor.position, setpoint, elapsed))
```

Also reusable from `ophyd_scan.py`: `protect_cleanup()` and `safe_step()` for
interrupt-safe teardown, and `detector_kind()`, `arm_hdf()`, `disarm_hdf()`,
`wait_for_frames()` for the detector.

`extra=` supplies **substitutions** (label text), not new column sources. To record
a second motor's readback today, use a dotted path — `huber.eta.user_readback`. To
make `motor2` a first-class source, add it to `_SPECIAL_SOURCES` and to
`_Column.value()` in `ophyd_spec_config.py`.

---

## Troubleshooting

| symptom | cause |
|---|---|
| `invalid SPEC column labels` at scan start | duplicate label, or a label with two consecutive spaces |
| `no device 'x' in the oregistry` at scan start | typo in a `source:`, or the device is not loaded — add `optional: true` if it is genuinely absent sometimes |
| `WARNING: SPEC column '...' is not a scalar` | a bare device path used as a column; append the signal (`.user_readback`, `.total`, `.mean_value`) |
| A column is silently missing | it was `optional:` and its device did not resolve — look for the warning at scan start |
| Viewer plots the wrong Y axis | the primary counter is no longer the last entry in `columns:` |
| Viewer refuses to overlay two scans | their column layouts differ, i.e. the template changed between them |
| Motor positions look wrong in the viewer | should not happen — the writer re-emits `#F` when `positioners:` changes; report it if you see it |
| Scan runs but no `.h5` appears | `save_img=0`, or the HDF plugin was not armed — check the printed file path |

The plot showing a perfect ramp with no motion is why both readback and setpoint
are recorded by default: compare the two columns to confirm the motor really moved.
