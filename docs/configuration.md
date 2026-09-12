# Configuration: where every setting lives

[← index](../README.md)

A handful of YAML files describe an experiment. One object, `expt`, is what the
code actually reads. This page explains the path from one to the other.

## The whole flow

```
 ┌─ EDITED BY HAND ─────────────────────────────────────────────────────────┐
 │                                                                          │
 │  configs/experiment.yml            ← changes a few times per cycle       │
 │      cycle_name, mount_point, experiment_name,                           │
 │      analysis_machine, workflow_name, use_subfolder                      │
 │                    │                                                     │
 │                    ├── selects ──────────┐                               │
 │                    │                     ▼                               │
 │                    │   user_plans/<cycle_name>/<experiment_name>/        │
 │                    │       sample_info.yaml       ← what/where the       │
 │                    │       measurement_info.yaml     samples are, and    │
 │                    │       trio_measurement_info.yaml  what to measure   │
 │                    │                                                     │
 │  plans/set/device_position.yaml    ← per-detector geometry:              │
 │      motors, db_x/db_y, distance, pixel_size, allow_motion               │
 └──────────────────────────────────────────────────────────────────────────┘
                      │
                      ▼
 ┌─ expt  (src/id8_common/expt_config.py) ──────────────────────────────────┐
 │                                                                          │
 │   static      from experiment.yml, read once at import                   │
 │   run state   set per measurement, from measurement_info + sample_info   │
 │   persistent  sample_index, file_name, mesh positions                    │
 │   EPICS       measurement_num  ->  8ideSoft:Reg1                         │
 │                                                                          │
 └──────────────────────────────────────────────────────────────────────────┘
          │                                    │
          │ read by everything                 │ mirrored to
          ▼                                    ▼
   master_plan → ad_acq → *_modes         state/run_state.yml
   nexus_utils (NeXus metadata)
   dm_util     (DM job submission)
```

## The four buckets in `expt`

`expt` is one attribute namespace over four different backends. You never need
to know which is which to *read* a value — `expt.cycle_name` and `expt.acq_time`
look identical — but it matters for knowing where to change something.

| bucket | source | lifetime | example |
|---|---|---|---|
| **static** | `configs/experiment.yml` | a cycle | `expt.cycle_name`, `expt.mount_point` |
| **run state** | `measurement_info.yaml` + `sample_info.yaml`, per measurement | one measurement | `expt.acq_time`, `expt.det_name` |
| **persistent** | `state/run_state.yml`, read back at startup | between sessions | `expt.sample_index`, `expt.file_name` |
| **EPICS** | `pv_registers` → `8ideSoft:Reg1` | outlives the checkout | `expt.measurement_num` |

```python
expt.cycle_name          # '2026-3'                → edit experiment.yml
expt.acq_time            # 0.01                    → edit measurement_info.yaml
expt.measurement_num     # 61                      → 8ideSoft:Reg1, automatic
expt.user_plan_dir       # .../user_plans/2026-3/comm202609
```

Two behaviours worth knowing:

* **An unset run field raises, it does not return `None`.** Before any
  measurement is configured, `expt.acq_time` raises `AttributeError` naming the
  cause. A typo like `expt.cycl_name` raises too. This is deliberate — a silent
  `None` used to surface much later as an unrelated crash.
* **YAML `no` is normalised.** A bare `sample_move: no` parses as the *boolean*
  `False`, which would never compare equal to `"no"`. `expt` coerces booleans
  back to `"yes"`/`"no"`, so hand-edited YAML behaves the way it reads.

## Switching to a different experiment

Edit two lines in `configs/experiment.yml`:

```yaml
cycle_name: "2026-3"
experiment_name: "comm202609"
```

Everything follows. The plan files are found at
`user_plans/<cycle_name>/<experiment_name>/`, resolved **at call time** — so
nothing is pinned to a path and no Python needs editing. A session that is
already open needs `expt.reload()` (or a restart) to see the edit —
`experiment.yml` is read at import. Both `run_measurement_info()` and
`run_trio_measurement_info()` print the directory they read from, so a wrong
setting is visible immediately:

```
Reading plans from /home/beams10/8IDIUSER/bluesky/src/user_plans/2026-3/comm202609
```

If the folder does not exist you get the path *and* the two settings that built
it, rather than a bare `FileNotFoundError`.

## What lives where — quick lookup

| Setting | File |
|---|---|
| cycle, experiment name, data root | `configs/experiment.yml` |
| where analysis runs, DM workflow | `configs/experiment.yml` |
| sample names, headers, mesh centres | `user_plans/<cycle>/<expt>/sample_info.yaml` |
| protocols: detector, mode, timings, frames | `user_plans/<cycle>/<expt>/measurement_info.yaml` |
| multi-detector (parallel) protocols | either file — conventionally `user_plans/<cycle>/<expt>/trio_measurement_info.yaml` |
| detector distance, beam centre, **pixel size** | `plans/set/device_position.yaml` |
| which devices exist at all | `configs/devices.yml`, `ad_devices.yml` |
| sample index, mesh positions | `state/run_state.yml` (managed, do not hand-edit) |
| the measurement counter | `8ideSoft:Reg1` (managed — hand-set only to recover a lost counter) |

**⚠ Pixel size is per detector, not per experiment.** eiger4M 75 µm, rigaku3M
76 µm, lambda2M 55 µm — all in `device_position.yaml`. It used to be a single
global, which meant every Rigaku file recorded the Eiger's value.

## `state/run_state.yml`

Three sections under a `_note:` header, and they are **not** equivalent:

```yaml
persistent:            # ← READ BACK at startup. Do not hand-edit while running.
  sample_index: 1
  file_name: A0060_Test_a0002_f000100_r00001
  measurement_num: 61  #   ...except this one, which is only a MIRROR of Reg1
  sample_positions: {1: -1, 2: 0, 3: 3}
static: {...}          # ← output only, a mirror of experiment.yml
run:    {...}          # ← output only, the measurement last configured
```

`static` and `run` exist so a failed run leaves behind what it was trying to do,
and so a GUI can see the current measurement without a live session. Nothing
reads them back. `persistent` is the real store for the sample index, the name
of the measurement in flight and the per-sample mesh positions — those used to
be `8ideSoft:` Reg6/StrReg8/Reg16-41.

**⚠** Deleting `state/run_state.yml` forgets every mesh position, and it is
gitignored, so a checkout does not restore it.

## The measurement counter stays in EPICS

`measurement_num` is the `NNNN` in every folder and file name
(`A0061_Test_a0002_f000100_r00001`). It only counts up, and **nothing downstream
checks whether a name is already taken** — so a counter that restarts at 0
silently overwrites existing measurements.

That is why it is the one value that did not move into `run_state.yml`. Its
store is `8ideSoft:Reg1`:

```bash
caget 8ideSoft:Reg1        # what the next measurement will be numbered
```

| | `run_state.yml` | `8ideSoft:Reg1` |
|---|---|---|
| survives `git clean -x`, a fresh clone, an accidental `rm` | ✗ | ✓ |
| shared between two sessions / two machines | ✗ | ✓ |
| visible to `caget` and non-Bluesky tools | ✗ | ✓ |

Every read of it is a real Channel Access get (`use_monitor=False`), not ophyd's
monitor cache, so a number another session or another process has just taken is
visible here immediately.

`run_state.yml` still records `measurement_num`, but only as a **mirror** — the
last value this checkout saw. It is read for exactly one purpose: if the
register comes back *lower* than the mirror (a soft-IOC restart that lost Reg1,
or someone clearing it by hand), the next read pushes the register back up and
prints a warning, rather than handing out numbers that already name files on
disk. The correction only ever goes upwards; a register that is ahead of the
mirror is always taken as correct.

If `pv_registers` is not connected, the session falls back to the mirror and
says so in red at first use; with no mirrored value either it raises rather than
guess a number. Check for existing files before acquiring in that state — the
counter is not shared with anything else while it lasts.

**⚠ Two naming streams share one counter.** `gen_folder_prefix()`
(`plans/acquire/acq_helpers.py`) increments `Reg1` once per call, and it is
called from both sides of the session: `det_acq_series()` names its output
`<mount_point><cycle>/<experiment>/data/A####…`, while the align scans
(`dscan_ophyd`, `dmesh_ophyd`, `ascan_ophyd`, … in `plans/align/ophyd_scan.py`,
and the Bluesky `plans/align/scan_8id.py`) name theirs `…/data/bluesky/A####…`.
So if you are working out what the counter should be from what is already on
disk, the highest `NNNN` may be under `data/bluesky/`, not `data/` — look in
both before setting `Reg1` by hand.

## Related

* [Running measurements](running-measurements.md) — protocol syntax and examples
* [Devices](devices.md) — `devices.yml` and friends
