# Demoing the rolled-back `main` branch

**Audience:** whoever runs the Friday demo (and the colleagues taking over
support on 2026-09-16, if they want to see what `main` does before touching it).

`main` was rolled back on 2026-09-10 to the 2026-08-04 tree, so that support
staff find the code they already know. Everything written since then lives on
`dev`.

**Since the 2026-09-11 split, each branch has a checkout of its own**, so
running `main` no longer needs a throwaway clone or a special launcher. It is
simply the other session.

---

## What exists

| | path | branch | who uses it |
|---|---|---|---|
| Ophyd checkout | `~/ophyd` | `dev` | the beamline, every real measurement |
| Bluesky checkout | `~/bluesky` | `main` @ `d5fed29` | the August tree |
| Ophyd launcher | `~/bin/start_ophyd.sh` | → `~/ophyd` | normal operations |
| Bluesky launcher | `~/bin/start_bluesky.sh` | → `~/bluesky` | the August tree |

Both are independent `git clone`s — a `.git` each, so deleting or moving one
cannot break the other. Each launcher pins `PYTHONPATH` to its own `src`,
because `cd` alone is not enough: see the sys.path warning below.

**`~/bluesky_main` and `~/bin/start_bluesky_main.sh` are obsolete.** They were
the throwaway clone this page used to describe, from before the split. Nothing
points at them and deleting them costs nothing.

The August tree gives you: `scan_8id.py` for all scans (Sam's code, driven
through `RE(...)`), the in-tree `nexus_utils.py` NeXus writer, `pv_registers`
instead of `expt_config`, and `dual_master_plan.py`. It does **not** contain
`ophyd_scan.py`, `expt_config.py` or `startup_ophyd.py`, and it knows nothing of
the `trio_*` modules or `tv_mode()`, which are dev-only.

---

## Before you start

1. **No live session.** Both trees drive the same motors, the same shutter and
   the same detectors. Two clients on that hardware is the real risk here, not
   anything about the code.

   ```bash
   ssh pearl "pgrep -af ipython"
   ```

   If that prints anything, stop and find out whose it is.

2. **Run it on pearl or amber.** Not kouga — no PVs there.

3. **Know that `measurement_num` is shared.** It is the one value still in
   EPICS (`8ideSoft:Reg1`), and both trees read and increment the same
   register. A demo acquisition consumes a real `A####` number out of the live
   sequence. That is survivable — nothing breaks — but the gap in the numbering
   is permanent, so say so rather than letting someone wonder later.

---

## Running it

```bash
~/bin/start_bluesky.sh
```

### ⚠ Which tree a session actually imports

`cd` does not decide this. `id8_bits` is installed **editable** in
`8ide_bits_test`, and `__editable__.id8_bits-0.0.1.pth` hardcodes
`/home/beams10/8IDIUSER/bluesky/src`, so `import id8_common` goes there from any
working directory. Both launchers therefore pin `PYTHONPATH` to their own `src`,
which lands earlier in `sys.path` than site-packages and wins.

`id8_common` is also a **namespace package** (no `__init__.py`), which means two
`src` roots on `sys.path` *merge* rather than one shadowing the other: a module
missing from the pinned tree silently falls back to the other checkout. As of
2026-09-11 the only files on `main` and not on `dev` are

    plans/acquire/dual_acq.py      utils/default_metadata.py
    plans/acquire/dual_master_plan.py   utils/nexus_utils.py
    user/ad_acq_zhou202607.py

and no module either startup imports pulls any of them in, so nothing on the
startup path can pick up August code. The visible consequence is small but worth
knowing: from an **Ophyd** session, `import id8_common.plans.acquire.dual_acq`
still succeeds — serving `main`'s August module.

The retired `start_bluesky_main.sh` solved this properly, by removing the other
entry from `sys.path` before `id8_common` is first imported and then printing
proof of which tree won. Porting that check into both launchers is the right
follow-up; until then the pin is what you have.

To leave, `exit` the ipython session as usual. Nothing to clean up.

---

## What to demo, in order of safety

**Safe — moves nothing:**

```python
dry_run_measurement_info()
dry_run_dual_measurement_info()
```

These validate the whole plan structure and move no hardware. This is the right
thing to show for "here is how a measurement is described."

**Needs staging first** (see the trap below):

```python
run_measurement_info()
run_dual_measurement_info()
```

**Do not run in a demo:**

```python
select_device("rigaku3M")   # or eiger4M, lambda2M
```

Reason under "Traps", below. It is a real detector move to a stale preset.

---

## Traps

These are all cases where the *code* is fine and the *state it reads* is not.
The August tree was correct against the world as it stood on 2026-08-04; parts
of that world have since moved.

### 1. The plan files are not where the August code looks

`master_plan.py` hardcodes:

```python
USER_PLAN_DIR = Path("/home/beams10/8IDIUSER/bluesky/src/user_plans")
SAMPLE_INFO_FILE       = USER_PLAN_DIR / "sample_info.yaml"
MEASUREMENT_INFO_FILE  = USER_PLAN_DIR / "measurement_info.yaml"
```

Note that this points into whichever tree the session imported — since the
split that is **`~/bluesky`, the `main` checkout** — not into the dev tree at
`~/ophyd`. And those flat files no longer exist: the current layout is
`user_plans/<cycle>/<experiment>/`, i.e.
`user_plans/2026-3/pope202609/sample_info.yaml`.

So `run_measurement_info()` and `dry_run_measurement_info()` will raise
`FileNotFoundError` out of the box. This is a *good* failure — loud and
immediate — but you have to plan for it.

To demo them, stage copies at the flat path:

```bash
cd ~/ophyd/src/user_plans/2026-3/pope202609   # the dev tree holds the real plans
cp sample_info.yaml measurement_info.yaml dual_measurement_info.yaml ../../
```

This is safe: the `dev` tree in `~/ophyd` resolves plan files through
`expt_config.USER_PLANS_ROOT / <cycle> / <experiment>`, so it never reads the
flat copies, and they are being written into the *other* checkout anyway. Remove them afterwards anyway, so nobody finds two versions later:

```bash
rm ~/bluesky/src/user_plans/{sample,measurement,dual_measurement}_info.yaml
```

Whether the August parser accepts the *current* YAML schema is untested —
`dry_run_*` is exactly the cheap way to find out, and it is the reason to run
the dry run before the real one.

### 2. `select_device()` would drive the detector to an August position

`device_position.yaml` in the clone carries the presets as they were before the
recalibration:

| | August (`~/bluesky_main`) | current (`~/bluesky`) |
|---|---|---|
| `rigaku3M` `detector.x` | **−250.0** | −275.0 |
| `rigaku3M` `db_x` | **1378** | 1052 |

`select_device()` moves horizontal and vertical translation and writes the
beam-centre registers. Calling it from the demo session physically moves the
detector 25 mm to a stale position and overwrites the live beam centre.

**Don't call `select_device()` from the demo session.** If someone does, the
recovery is to re-run it from a normal `~/bin/start_bluesky.sh` session, which
will restore both the position and the registers from the current YAML.

### 3. `StrReg5` (analysis machine) is stale

The August code reads the analysis machine from EPICS. Live value:

```
StrReg1  2026-3               cycle_name        current
StrReg2  WA-XPCS              geometry          stale-looking
StrReg3  /gdata/dm/8ID/8IDE/  mount_point       current
StrReg4  pope202609           experiment_name   current
StrReg5  polaris              analysis_machine  STALE
StrReg6  xpcs8-boost-corr     workflow_name     current
StrReg7  no                   use_subfolder     current
```

`dev` reads `analysis_machine` from `configs/experiment.yml`, which now says
`local` (resolved to adamite or amazonite at submit time). The register still
says `polaris`. A demo acquisition that submits a DM job would submit it to polaris — which, per
`project_gladier_boost_corr_bug`, silently produces nothing while DM reports
SUCCEEDED.

If the demo goes as far as submitting a job, either accept that the job is
meaningless, or `caput 8ideSoft:StrReg5 adamite` first **and put it back
afterwards** — the live `dev` session ignores this register, but the next
person to read it should not be misled.

---

## "If it starts up, will it run?"

Short answer: **no, not reliably — and the three traps above are why.**

A clean startup is a genuine and useful signal. It proves:

- every module imports (no syntax or missing-dependency breakage from the rollback)
- the IOCs are up and the devices connect
- `ad_setup()` succeeded and the area detectors are wired
- the pre-flight check passed, so you really are on the August tree

That is most of what could go wrong *with the code*, and it is why the launcher
does its checking up front rather than partway through a plan.

But startup exercises almost nothing that reads **state**. Every trap above is
invisible at startup and only appears when a plan actually runs:

- the plan-file paths are `Path` objects at import time — no I/O happens until
  `run_measurement_info()` is called, so a missing file cannot be detected at
  startup
- `device_position.yaml` is read inside `select_device()`, not at import
- the EPICS string registers are read per-plan

There is also a category startup cannot speak to at all: whether the August
code paths still *agree* with things outside the repo that have moved since —
the `nexus_xpcs_aps` install, the DM workflow definitions, the current YAML
schema. The August NeXus writer is in-tree (`nexus_utils.py`) so it is
self-consistent, which helps; the DM side is not.

**So the useful rule for the demo:** a clean startup means the session is real
and connected, and you can show devices, read positions, and talk through the
code with confidence. Before showing anything that *runs*, use
`dry_run_measurement_info()` — it is the step that actually exercises the state,
and it moves nothing.

---

## If you need to abandon the demo

```bash
rm -f ~/bluesky/src/user_plans/{sample,measurement,dual_measurement}_info.yaml
```

`~/ophyd` — the dev tree the beamline runs — is untouched throughout: the only
thing this procedure writes is those staged copies, and they land in the `main`
checkout.
