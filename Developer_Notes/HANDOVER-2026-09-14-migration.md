# Migration handover — 2026-09-14

Three capabilities from `~/ophyd` (`dev`), rebuilt on `~/bluesky`'s own
architecture. **Nothing here has been run against hardware.** No session was
started, no detector armed, no file written, no DM job submitted.

## Where things are

| path | branch | what it is |
|---|---|---|
| `~/bluesky` | `main` @ `d5fed29` | **Untouched.** What `~/bin/start_bluesky.sh` starts, exactly as before |
| `~/bluesky_trio` | `migrate/ophyd-capabilities-20260914` @ `6626d41` | The migrated code. Started by `~/bin/start_bluesky_trio.sh` |
| `~/ophyd` | `dev` @ `49d3ec0` | **Untouched**, including its 231 uncommitted files |

`~/bluesky_trio` is a git worktree of the same repo, so it has its own branch but
shares `.git`. Branch pushed to `origin`.

```bash
~/bin/start_bluesky.sh        # unchanged: ~/bluesky, main. For Sam and Suresh
~/bin/start_bluesky_trio.sh   # new: ~/bluesky_trio, the migrated code
```

The only difference between the two launchers is `PYTHONPATH`. Do not run both
at once against the same detector.

## What was migrated, and what was not

| # | Capability | State |
|---|---|---|
| 1 | `sample_info` / `measurement_info` as CSV | Ported, verified offline |
| 2 | Metadata after Ctrl+C | **Not ported** — dropped 2026-09-14 |
| 3 | Eiger External Series, Rigaku Fast Transfer | **Not ported** — dropped 2026-09-14 |
| 4 | Startup does not disturb a running acquisition | Ported, verified offline |
| 5 | Trio detector acquisition | Ported, verified up to the hardware boundary |

No acquisition mode was added or changed. `dual_leg_behaviour()` still rejects
`eiger4M/External Series`, which is the check that 3 stayed out.

## The three changes

### 4. Startup no longer fires an Eiger exposure

`ad_setup()` called apstools' `AD_prime_plugin2()`, which fires a **real
exposure**. It was gated on `AD_plugin_primed()`, which compares
`cam.data_type` with `hdf1.data_type` — on the eiger4M those are `Int8` and
`UInt32`, permanently unequal (re-verified against live PVs 2026-09-14). So the
gate never passed and **every session start fired an exposure**. That is what
crashed measurement G0209 on 2026-09-12.

Priming is not needed: eiger4M, lambda2M and rigaku3M all run hdf1
`LazyOpen=Yes` in `Stream` mode (verified against live PVs), which per apstools'
own `AD_plugin_primed` docstring removes the need. `ALLOW_AREA_DETECTOR_WARMUP`
is now a dead key, set `False`.

`utils/startup_guard.py` is the backstop: it refuses `EpicsSignal` `put`/`set`
attempted while the startup module is on the stack, prints what it blocked, and
**never raises**. Writes typed at the prompt afterwards are unaffected.

### 1. Plans can be written as CSV

`master_plan.read_yaml()` dispatches on the `.csv` suffix; the `plan_csv` import
is inside that branch, so a fault in it cannot reach a YAML session. Switch a run
by passing a `.csv` path as `measurement_info_file` / `sample_info_file` — no
code edit, reversible by passing the `.yaml` back.

The converter emits whichever of this tree's two schemas fits:

* **one** detector → the serial shape `master_plan.py` reads (`detector:` /
  `mode:` singular at protocol level)
* **two or more** → the `detectors:` list `dual_master_plan.py` reads
* `parallel,yes` in `#MEASUREMENT` forces the list for a single leg

### 5. Trio acquisition

`dual_acq.py` and `dual_master_plan.py` were already written against
`measurement["detectors"]` as a list of *any* length, so this adds only:

* the `("lambda2M", "Internal")` entry in `DUAL_LEGS` — an **existing serial
  mode** whose setup half, `setup_lambda_internal()`, is reused verbatim
* `run_trio_measurement_info()` / `dry_run_trio_measurement_info()`, thin
  wrappers over the dual functions

**The shutter contract is unchanged.** Exactly one leg is the shutter owner (the
Rigaku, which still gates the beam via softglue in `Start with Trigger`), it is
armed first, and the other two follow once it confirms it is acquiring. `~/ophyd`
removed the owner by switching the Rigaku to `Fixed Time`; that is an
acquisition-behaviour change and was deliberately not brought over. **So this
tree's `trio_measurement_info.yaml` needs `shutter_owner: yes` on the Rigaku leg,
where `~/ophyd`'s rejects it.** That is the one format difference between them.

**A trio run moves no huber axis.** `setup_huber_for_trio()` reports the position
and moves nothing, matching the live decision in `~/ophyd` where the equivalent
motion is commented out: three detectors have different presets in
`device_position.yaml` and at most one can be satisfied. Position all three
yourself first. The dual path keeps `setup_huber_for_dual()` and its existing
motion to delta 10 / nu 0, via a defaulted argument — **dual runs behave exactly
as before.** The attenuator is the only thing a trio run drives by itself.

## What has been verified, and how

Everything below ran on kouga with no EPICS and no beam.

* `python scripts/check_plan_csv.py` → **3/3 matched**, exit 0. Each pair
  compares the converter against a YAML written **by hand from the format spec**,
  not dumped from the converter. Covers the sample table and both protocol shapes.
* With `oregistry` stubbed (plan modules only *look up* devices at import, they
  do not instantiate), the **real** validators ran on the trio plan from both
  YAML and CSV: 3 legs `rigaku3M`/`eiger4M`/`lambda2M`, exactly one shutter
  owner, 3 leg specs, 10 s each.
* `dual_leg_behaviour()` resolves all three trio pairs and still rejects
  `eiger4M/External Series`.
* `arm_startup_guard()` / `report_startup_writes()` patch and restore
  `EpicsSignal.put`/`.set` cleanly, and refuse a second arm.

## Suggested test order

Beam was down on 2026-09-14 (ring current 0.007 mA), so the first three steps
work as dark frames; only the data is meaningless, not the plumbing.

1. **Startup.** `~/bin/start_bluesky_trio.sh`. Watch for
   `[startup_guard] disarmed: no EPICS writes attempted during startup`. If it
   names a blocked write, that write is a real startup-time hardware touch —
   report it rather than allow-listing it blindly.
   Confirm the Eiger took no frame: read `8idEiger4m:cam1:ArrayCounter_RBV`
   before and after; it must not change.
2. **Dry run.** `dry_run_trio_measurement_info(check_hardware=True)`. This is the
   first thing that checks the three detectors are actually connected.
3. **One leg at a time.** A `detectors:` list of length one is legal and becomes
   its own shutter owner — smoke-test the Rigaku, then the Eiger, then the
   Lambda, before running all three.
4. **The trio.** `run_trio_measurement_info()`. 10 frames × 1 s.
5. **Analysis.** Confirm three folders under
   `/gdata/dm/8ID/8IDE/2026-3/comm202609/data/` sharing one run number, each with
   its own `_metadata.hdf`; then `dmjob.sh <UUID>` per leg, and results under
   `.../comm202609/analysis/`.

### Two things to know before step 4

* **`lambda2m_qmap_default.hdf` was copied into `comm202609/data/` on
  2026-09-14** (from `pope202609`). It was missing, and DM stalls at
  `02-WAIT-QMAP` without it. It is the generic default, not a qmap measured for
  this geometry — fine for a smoke test, replace it for real data.
* **Beam-centre metadata for the non-Rigaku legs is probably wrong.** With no
  `geometry:` block a leg's metadata comes from `device_position.yaml`, which was
  calibrated with each detector at its own preset. In a trio at most one is.
  Raw frames are unaffected; only `_metadata.hdf`. Add a `geometry:` block per
  leg to correct it.

## Undoing all of it

```bash
cd ~/bluesky && git diff checkpoint/bluesky-main-20260914   # expect: empty
git worktree remove ~/bluesky_trio                          # drop the migrated tree
rm ~/bin/start_bluesky_trio.sh
```

`~/bluesky` was never changed, so there is nothing to revert there. Restore
points, both pushed:

* `checkpoint/bluesky-main-20260914` → `d5fed29`, `~/bluesky` main
* `backup/ophyd-worktree-20260914` → `8dba7b0`, an exact snapshot of `~/ophyd`'s
  231 uncommitted files, taken with a temporary index so that tree was never
  touched
