# Running measurements

[← index](../README.md)

Two ways in. Almost always use the first.

| | Entry point | You edit | Use when |
|---|---|---|---|
| **High level** | `run_measurement_info()` | `measurement_info.yaml` | normal operation |
| **Low level** | `det_acq_series()` | Python at the prompt | scripting something the YAML cannot express |

## The normal workflow

```python
dry_run_measurement_info()    # expand + validate + preview. Moves nothing.
run_measurement_info()        # go
```

`dry_run` expands the YAML and checks every measurement it produces: the
detector and mode names, the timing rules, the frame and repeat counts, the
analysis type, and any detector-position override. Use it after every edit —
it costs seconds.

**⚠ The serial dry run does not touch hardware.** Device connectivity and the
sample mesh are checked by `run_measurement_info()` itself, per measurement,
immediately before that measurement runs — so a missing device or an
unreachable mesh motor surfaces when its turn comes, part way down a long list,
rather than up front. (The trio path gates this differently — see
[Trio-detector acquisition](#trio-detector-acquisition).)

## Anatomy of `measurement_info.yaml`

Two sections. `runs:` says *what to execute*; `protocols:` says *how*.

```yaml
loop_order: sample_major

runs:
 - name: mode_test
   samples: [1]                       # indices into sample_info.yaml
   protocols:
     - eiger_internal_series          # executed in this order
     - rigaku_epics
   repeats: 1                         # repeats of the whole run block

protocols:

  eiger_internal_series:
    detector: eiger4M                 # must be a key of ACQ_MODES
    mode: Internal Series             # must be a mode of that detector
    att_level: 2
    acq_time: 0.01
    num_frames: 100
    num_repeats: 3                    # repeats of THIS protocol
    wait_time: 0
    sample_move: no
    position_reset: no
    analysis_type: Multitau
    qmap_file: eiger4m_qmap_default.hdf
```

**`repeats` vs `num_repeats`** — `repeats` in the run block re-runs the whole
list of protocols; `num_repeats` inside a protocol repeats that one measurement,
each getting its own `_r0000N` file. They multiply.

**⚠ A timing key the mode does not use is rejected, not ignored.**
`num_segments` or `trigger_period` on a mode whose table does not declare it
raises at validation. This is deliberate: a silently ignored timing parameter
is how you get data that is not what you asked for. Those two are the only
keys checked this way — a misspelt field name is still ignored.

See [Detector modes](reference/detector-modes.md) for which fields each mode
requires, and the timing rules.

## What happens per measurement

```
run_measurement_info()
   └─ expand runs × samples × protocols  →  a list of measurements
       └─ for each:  validate  →  load into expt  →  att()  →  select_device()
           └─ det_acq_series()
               ├─ gen_folder_prefix()      A0061_Test_a0002   (bumps the counter)
               ├─ for each repeat:
               │    ├─ sample_mesh_move()  (only if sample_move: yes)
               │    ├─ setup_<mode>()      arm the detector
               │    ├─ acquire_<mode>()    open shutter, collect, close shutter
               │    ├─ NeXus metadata      → <name>_metadata.hdf
               │    └─ DM job submitted
```

File naming: `A0061_Test_a0002_f000100_r00001` —
`<header><counter>_<sample>_a<attenuation>_f<frames>_r<repeat>`.

**⚠ The `a0002` comes from the attenuator *readback*, not from `att_level`.** If
you ask for 20 and it reads back 22, the folder says `a0022`.

The `0061` is `expt.measurement_num`, and its store is the EPICS register
`8ideSoft:Reg1` — not `state/run_state.yml`, which keeps only a mirror of the
last value this checkout saw. `gen_folder_prefix()` reads the register, builds
the name from it, and writes it back one higher, so the counter advances once
per measurement (the repeats of one measurement share a prefix and differ only
in `_rNNNNN`). A trio measurement bumps it once too, and both legs share the
number.

```bash
caget 8ideSoft:Reg1        # the number the next measurement will get
```

**⚠ One counter, two naming streams.** The alignment scans call
`gen_folder_prefix()` as well, but write into `data/bluesky/`, not `data/`, so
the highest `NNNN` on disk may be under `data/bluesky/`. How far a shift of
alignment scans drags the counter depends on which module you used.
`ophyd_scan.py` (`dscan_ophyd` and the rest) bumps it exactly once per scan,
saved images or not. `scan_8id.py` bumps it **twice** for a `save_img=1` scan —
the caller and `save_images()` each call `gen_folder_prefix()`, so the number in
the run metadata's `image_file` is one behind the number in the `.h5` name it
actually wrote — and not at all for `save_img=0`.

[Configuration → the measurement counter](configuration.md#the-measurement-counter-stays-in-epics)
covers what happens when the register goes backwards, and what a session does
if `pv_registers` is not connected.

## Low level: `det_acq_series()`

Set the run state yourself, then call it:

```python
expt.measurement_num = 12          # EPICS write -- 8ideSoft:Reg1

expt.set_measurement(
    measurement={
        "detector": "eiger4M",
        "mode": "Internal Series",
        "acq_time": 0.1,
        "acq_period": 0.1,
        "num_frames": 1000,
        "num_repeats": 5,
        "sample_move": "no",
        "qmap_file": "eiger4m_qmap_default.hdf",
    },
    sample={"header": "A", "sample_name": "G10"},
)

det_acq_series(wait_time=0)
```

Individual fields work too — `expt.num_frames = 2000`.

`expt.measurement_num` is the odd one out: assigning to it is a real write to
`8ideSoft:Reg1`, shared with every other session and with the alignment scans,
and the value you set is not checked against what is already on disk. Lower it
only if you mean to overwrite.

With sample motion, the sample index and the per-sample mesh position are
persistent state (the `persistent:` block of `state/run_state.yml`) rather than
run state:

```python
expt.sample_index = 1
expt.sample_position(1)               # where sample 1's mesh is now; -1 = start
expt.set_sample_position(1, -1)       # restart the mesh for sample 1
```

then add `"sample_move": "yes"` plus `inner_motor` / `outer_motor` /
`inner_center` / `outer_center` / `inner_range` / `outer_range` / `inner_pts` /
`outer_pts` to the `sample=` dict.

## Trio-detector acquisition

Runs the Rigaku, the Eiger and the Lambda **in the same beam window** instead of
back to back. They are *not* frame-synced — the acquisitions merely overlap, so a
set of long measurements costs roughly the slowest one rather than the sum.

"Trio" is the original name, not a limit. Any number of legs works, and the
supported `(device, mode)` pairs are the keys of `TRIO_LEGS` in
`plans/acquire/trio_acq_rigaku3m_eiger4m_lambda2m.py`:

| device | mode |
|---|---|
| `rigaku3M_epics` | `EPICS` |
| `eiger4M` | `Internal Series` |
| `lambda2M` | `Internal` |

`lambda2M` `External` is deliberately absent: it gates the shutter through
softglue, one pulse per frame, which cannot coexist with a Rigaku-owned shutter
window. Only the Lambda's `Internal` mode can share one.

```python
dry_run_trio_measurement_info(check_hardware=True)
run_trio_measurement_info()
```

**`check_hardware` defaults to `False`** — pass `True` on a live session, as
above. It gates the whole per-leg check: mode name, the `acq_time` floor,
`qmap_file`, `analysis_type`, the `geometry:` block, a leg's `motors:` block,
and every device the mode drives. Left off you still get run expansion, the
required protocol fields, the duplicate-label and shutter-owner rules, and the
sample-mesh check — which resolves the mesh motors either way, so a
`sample_move: yes` protocol still needs a live session.

Both front ends now run the same field checks, from
`plans/acquire/validators.py`. That closed a real gap on this side: the shared
`require_mode_devices()` checks a mode's `required_devices` as well as its
`hardware_device`, so a trio protocol with an `eiger4M` `External Series` leg is
now rejected at dry-run time when `softglue` is missing. It used to pass
validation and fail only once the run was underway.

The config is `trio_measurement_info.yaml`, in the same folder as the others. A
protocol carries a `detectors:` **list** instead of a scalar `detector:`:

```yaml
  trio_att2:
    att_level: 2
    num_repeats: 1
    sample_move: no
    position_reset: yes
    detectors:
      - device: rigaku3M_epics
        mode: EPICS
        acq_time: 1
        num_frames: 3000
        qmap_file: rigaku3m_qmap_default.hdf
        analysis_type: Multitau
        select_device: yes
        shutter_owner: yes          # exactly one leg must own the shutter
        start_timeout: 30
        hdf_timeout: 600
      - device: eiger4M
        mode: Internal Series
        acq_time: 1
        num_frames: 3000
        qmap_file: eiger4m_qmap_default.hdf
        analysis_type: Multitau
        hdf_timeout: 600
      - device: lambda2M
        mode: Internal
        acq_time: 1
        num_frames: 3000
        qmap_file: lambda2m_qmap_default.hdf
        analysis_type: Multitau
        hdf_timeout: 600
```

Each step produces **one output per leg, all sharing one run number**, each with
its own folder, its own `_metadata.hdf` and its own DM job. The names agree up to
the detector label:

```
A0058_HEA_a0022_f003000_rigaku3M_r00001
A0058_HEA_a0022_f003000_eiger4M_r00001
A0058_HEA_a0022_f003000_lambda2M_r00001
```

`use_subfolder` decides where those land. Under `no` — what
`configs/experiment.yml` currently sets — each is a folder directly under
`data/`. Under `yes` each sits one level deeper, inside a per-measurement
`data/A0058_HEA_a0022_f003000_<label>/` folder.

**A one-leg `detectors:` list is legal**, and that leg becomes the shutter owner
automatically. That is the right way to smoke-test each detector before running
them together.

### Motion in a trio run

Moves once per measurement, in this order:

| | |
|---|---|
| `filter_8ide.attenuation` | to `att_level`, in `run_trio_measurement()` |
| `huber.delta → 10`, `huber.nu → 0` | `setup_huber_for_trio()`, still before any acquisition |
| `detector.x`, `detector.y` | via `select_device()` on each leg with `select_device: yes` — **but only while that detector's `allow_motion` is true in `device_position.yaml`, and `eiger4M`, `rigaku3M` and `lambda2M` are all currently `false`, so this row moves nothing today** |
| whatever a leg's `motors:` block names | `move_leg_motors()`, after the `select_device` pass and before the shutter window |

The first two are in `run_trio_measurement()`, the last two inside
`trio_acq_series()`.

**⚠ `setup_huber_for_trio()`'s motion is currently commented out** for testing.
Re-enable the commented lines in `trio_master_plan_rigaku3m_eiger4m_lambda2m.py` before a
real trio run, or the leg geometry in the metadata will not describe the true
beam path. They reference `oregistry`, which that module does not currently
import, so add the import at the same time.

Never moves during acquisition, and refuses if asked: once the huber is set,
`huber.delta` and `huber.nu` are off limits — rejected in a leg's `motors:`
block at validation *and* at run time, and rejected as `sample_info.yaml`'s
`inner_motor`/`outer_motor`. To change where a trio run acquires, edit
`TRIO_HUBER_DELTA` / `TRIO_HUBER_NU`, not the YAML.

`setup_huber_for_trio()` is independent of `master_plan.py`'s
`placeholder_rigaku3M()`. Changing the serial placeholder does not move the trio
geometry, or vice versa.

### Per-leg metadata

In a trio run each leg needs its own detector name, qmap and geometry. That is
handled by a scoped swap around each leg's metadata write, so the two
`_metadata.hdf` files are stamped correctly and independently.

With no `geometry:` block, a leg's metadata comes from `device_position.yaml` as
a single-detector run would — but in a trio run at most one detector can be at
its calibrated preset, so add a `geometry:` block to the other leg. Every field
is optional and takes a literal number or a dotted ophyd path read live:

```yaml
        geometry:
          db_x: 716
          db_y: 815
          distance: 11.5            # metres
          pixel_size: 75.0e-6       # metres
          position_x: 256.0         # MILLIMETRES, written as metres
          position_y: 77.0
          swing_horizontal: huber.nu
```

## TV mode: a live view on all three detectors

`tv_mode()` free-runs the Rigaku, the Eiger and the Lambda together so you can
watch the images. It is the scripted version of leaving a detector running from
its own GUI, which is what "TV mode" has always meant here.

```python
tv_mode()                              # 5000 frames of 1 s, all three
tv_mode(acq_time=0.1, num_frames=200)  # something shorter
tv_mode(detectors=["lambda2M"])        # just one
```

| | |
|---|---|
| Writes | **nothing** — no HDF file, no NeXus metadata, no DM job, and no measurement number consumed |
| Beam | `shutteroff()` + `showbeam()` before arming, `blockbeam()` on the way out. **Set `att()` before you start it** — the session is blocked once it is running |
| Blocks | until every detector has finished its frames: 5000 s ≈ 83 min at the defaults |
| Ctrl+C | presses Stop on all three detectors and closes the shutter, then raises `RuntimeError` |

The 5000 is a ceiling, not a target — Ctrl+C is the normal way to end a TV
session. Cleanup is in a `finally`, so no path out of `tv_mode()` leaves a
detector running or the beam on the sample. A progress line prints every 30 s.

The HDF plugins are explicitly **disarmed**, not merely left alone: a plugin
that an aborted measurement left capturing would otherwise quietly append
live-view frames to that measurement's file.

Trigger modes are the free-running one on each detector — Eiger `Internal
Series`, Lambda `Internal`, Rigaku `Fixed Time`. The Rigaku has no mode called
"internal series"; `Fixed Time` is its internally-timed one, needing nothing
from softglue. **That choice is confirmed as a legal enum value but has not yet
been run on the detector** — if it turns out not to free-run, the proven
fallback is the acquisition path's `Start with Trigger` plus
`softglue.enable_rigaku = '1'`.

## ⚠ Live/TV mode blocks an exposure change

A camera left free-running ignores writes to `acquire_time`. The setup function
returns as though it worked, and the measurement then runs at the live view's
exposure instead of the protocol's — frames written, metadata recording the
*requested* value, and the exposure the data was actually taken at simply lost.

`trio_acq_series()` guards against this in two steps per repeat, both **before**
`showbeam()` so a failure costs no beam:

1. `stop_live_mode()` presses Stop on every leg and waits for it, before any
   leg's setup runs.
2. `confirm_acq_time()` then reads each camera's `acquire_time` **readback** and
   raises if it is not what the protocol asked for (5% tolerance, for the
   detectors' own quantisation). This checks the symptom, so it catches any
   other reason a write did not stick, not just live mode.

`tv_mode()` makes the same two checks. A leg can raise its stop allowance with
`stop_timeout:` in the protocol; the default is 10 s.

## ⚠ The qmap must already be on disk

A plan names a qmap by bare file name, and it is resolved against the
experiment's `data/` directory — the same place DM looks for it. Since
2026-09-08 both the dry run and the real run check the file is there and refuse
to start if it is not, naming what *is* available:

```
qmap 'rigaku3m_qmap_Sq360_Dq18_Sphi16_Dphi1_lin.hdf' not found at
  /gdata/dm/8ID/8IDE/2026-3/comm202609/data/…
  Available in that directory: ['eiger4m_qmap_default.hdf', 'rigaku3m_qmap_default.hdf']
```

Before that check existed, only the name's non-emptiness was validated, so a
typo or a qmap that was never copied into the experiment stayed invisible until
the analysis stage — after the data was taken.

## Related

* [Detector modes](reference/detector-modes.md)
* [Configuration](configuration.md)
* [Troubleshooting](reference/troubleshooting.md)
