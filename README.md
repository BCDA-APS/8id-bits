# 8-ID XPCS instrument

Beamline control for the APS 8-ID XPCS instrument. The working code is
`src/id8_common`; this page is the way in.

Start here and follow the links. Every page is plain Markdown, so it reads on
GitHub, in an editor, or with `less`.

## Starting

On a beamline host that can see the EPICS network (`pearl` or `amber`):

```bash
~/bin/start_ophyd.sh      # Ophyd-only session -- devices and plans, no RunEngine
~/bin/start_bluesky.sh    # the same, plus Bluesky (RunEngine, databroker)
```

Either one activates the `8id_bits` conda environment and sources the APS Data
Management setup. `start_ophyd.sh` is the normal way in; see
[Starting a session](docs/starting-a-session.md) for what it loads and how to
read the output.

## I want to…

| … | Go to |
|---|---|
| start a session at the beamline | [Starting a session](docs/starting-a-session.md) |
| understand what got loaded, or add/remove a device | [Devices](docs/devices.md) |
| know where a setting lives, or change the experiment | [Configuration](docs/configuration.md) |
| read or set `expt` at the prompt, or work out why a field raised | [Using `expt`](docs/using-expt.md) |
| add a new device's reading to the NeXus metadata file | [Adding metadata fields](docs/adding-metadata.md) |
| check, debug or re-run a DM analysis job | [Data Management](docs/reference/data-management.md) |
| know why we use our NeXus writer and not Miaoqi's | [How the NeXus file is written](docs/reference/nexus-writer.md) |
| check the 2026-09-08 six-mode verification run | [Verification 2026-09-08](docs/reference/verification-2026-09-08.md) |
| find where session state is kept | [Configuration → `state/run_state.yml`](docs/configuration.md#staterun_stateyml) |
| know why the counter lives in EPICS | [Configuration → the measurement counter](docs/configuration.md#the-measurement-counter-stays-in-epics) |
| define a protocol and run a measurement | [Running measurements](docs/running-measurements.md) |
| plot a scan, live or afterwards — and know when to ask Peco | [Viewing scans](docs/viewing-scans.md) |
| set up or debug the scan viewer itself — checkout, conda env, what we changed | [The scan viewer](docs/reference/scan-viewer.md) |
| look up a detector mode or its timing rules | [Detector modes](docs/reference/detector-modes.md) |
| work out why something failed | [Troubleshooting](docs/reference/troubleshooting.md) |

## The 30-second version

```
   configs/experiment.yml          which cycle, which experiment, where data goes
   user_plans/<cycle>/<expt>/
       sample_info.yaml            what the samples are and where they sit
       measurement_info.yaml       what to measure and how
                 │
                 ▼
            expt  (id8_common/expt_config.py)     ← one object, everything reads it
                 │
                 ▼
        run_measurement_info()      expands YAML into measurements, validates, runs
                 │
                 ▼
        det_acq_series()            drives the detector, writes NeXus metadata,
                                    submits the DM analysis job
```

Every folder and file is named `A0061_<sample>_a0007_…`. That `0061` is
`expt.measurement_num`, and its store is the EPICS register `8ideSoft:Reg1`, not
a file in this repository — so a second session, or a fresh clone, carries on
from the same number instead of restarting at 0 over existing data.
Generated session state (sample index, mesh positions) lands in the gitignored
`state/run_state.yml`. Both are covered in [Configuration](docs/configuration.md).

Two commands cover most days:

```python
dry_run_measurement_info()   # validate and preview -- moves nothing
run_measurement_info()       # go
```

## Known issues and to-do

Ordered by what will bite you soonest. Each line says where the problem is, so
it can be picked up without re-deriving it.

### Will bite you during a run

* **⚠ `eiger4M` Internal Enable can hang indefinitely.** The wait at
  `plans/acquire/eiger4m_modes.py:337` (`acquire_eiger_internal_enable`) has no
  deadline and no fault check, so a trigger the detector never services stalls
  the session with the shutter open. Seen 2026-09-08: 12 minutes before a manual
  `^C`, then the same protocol ran in 21 s. Intermittent, not reproducible on
  demand. `plans/acquire/acq_wait.py` already has the bounded, fault-aware
  helper — this call site was never migrated.
* **Five more unbounded waits** in the same family: cam and HDF drains at
  `eiger4m_modes.py:299,305,343`, `lambda2m_modes.py:135,141`, plus
  `lambda2m_modes.py:180`, which busy-spins with no `sleep()` at all.
  `rigaku3m_modes.py:207` has its own private `_wait_for_state` rather than
  using `acq_wait`.
* **Two commented-out metadata lines are booby-trapped.** Uncommenting
  `/entry/instrument/incident_beam/fofb_s09_horizontal` or `…_vertical` in
  `utils/nexus_runtime.py` raises `KeyError` — there is no schema node — and
  `det_acq_series()` swallows it, so the measurement finishes and writes **no
  metadata file at all**. Add the schema node first. The other 28 commented
  lines are safe.

### This beamtime

* **The Alicat pressure field is not in the NeXus file yet.** `pcd1`/`pcd2` load
  and read, but nothing writes them to metadata. Two schema nodes plus two
  runtime lines — see
  [Adding metadata fields](docs/adding-metadata.md#case-2--a-field-he-does-not-model-eg-a-pressure-controller).
* **No Eiger measurement has been verified with beam.** The 2026-09-08 six-mode
  run was taken at 1.09 mA, so all four Eiger `g2` arrays are `NaN` — zero
  counts, nothing to correlate. The pipeline is proven; the numbers are not.
  Repeat at least Internal Series with beam.
* **A 100 000-frame Rigaku ZDT dataset is ~118 GB** and takes 2–7 minutes to
  correlate, loading each 35 GB module through CPU RAM. Fine on adamite (1 TB);
  it would fail on a smaller box. Budget for it before queueing many.

### Waiting on someone else

* **[AZjk/nexus_xpcs_aps#1](https://github.com/AZjk/nexus_xpcs_aps/pull/1)** —
  when it merges, delete the five local patches listed in
  [How the NeXus file is written](docs/reference/nexus-writer.md#local-patches-we-carry),
  and the scratch clone at `~/pr_work/nexus_xpcs_aps`.
* **`keysight_freq` and `keysight_amp` write `units="any"`** because their schema
  nodes *declare* `NX_ANY`. Raised in the PR as a question; not fixed there,
  because it may be deliberate.

### Cleanup, whenever

* **`conda activate 8id_bits` runs `8ide_bits_test`'s interpreter** — both
  `8id_bits/bin/pip` and `bin/ipython` have shebangs pointing there. Harmless
  today (the package is installed in both), but testing an import in "8id_bits"
  tests the wrong environment.
* **A stale second clone** of Miaoqi's package at
  `~/Documents/Miaoqi/nexus_xpcs_aps` (`main` @ 73d0be4, Aug 2025, with
  uncommitted edits). Nothing should point at it; the live one is
  `nexus_xpcs_aps_95ab368`.
* **`utils/Archive/*.txt`** — the retired NeXus writer, kept for reference.
  Deletable once the new path has a cycle behind it.
* **Three modules still bracket-index the registry at module scope** —
  `plans/set/hooks.py`, `volt_seq.py`, `ur5_pipetting.py`. None is imported at
  startup, so they are latent rather than dangerous, but each raises `KeyError`
  on import today.
* **Rarely-used devices with real bugs**, deferred by decision:
  `devices/lakeshore.py:set_temp_lakeshore()` raises `NameError`;
  `utils/misc.py:temp2str()` raises `IndexError` on an integer temperature.
* **`id8_common_dev/`** is well behind `id8_common` — treat as stale, not as a
  sibling to keep in sync.

## Conventions in these pages

* `expt` is the single configuration object. If a page says "`expt.cycle_name`",
  you can type exactly that at the session prompt.
* Paths are relative to the top of this repository — except the shorthands
  used throughout: `configs/…` and `plans/…` are inside `src/id8_common/`, and
  `user_plans/…` is inside `src/`. The pages themselves live in `docs/`.
* Anything marked **⚠** is a known trap, not a suggestion.
