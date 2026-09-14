# Evaluation — 2026-09-13: stack consolidation (8-ID + 9-ID), and the YAML→CSV question

Written so the thread survives losing the chat it came from. Three separate
questions got answered in one session; they are only loosely related and can be
acted on independently.

| § | Question | Answer |
|---|---|---|
| 1 | Is 8-ID "on Bluesky"? | No, and it never was — only the align scans were |
| 2 | Should 8-ID and 9-ID merge code bases? | Yes, and it is ~2/3 done already. Target Jan 2027 |
| 3 | Can hklpy2's math be used without ophyd? | Yes for *use*, no for *install*. See the measured answer |
| 4 | Convert `sample_info` / `measurement_info` to CSV before Wed 2026-09-16? | **No** — but for schedule/handover reasons, not format ones. See the §4.1 correction |

Repos involved, easy to confuse:

| repo | path | branch | role |
|---|---|---|---|
| `BCDA-APS/8id-bits` | `~/ophyd` | `dev` | live development; ophyd-only session |
| `BCDA-APS/8id-bits` | `~/bluesky` | `main` | rolled back to the 2026-08-04 tree for the 09/16 support handover |
| `AdvancedPhotonSource/BLUETELLA_9ID` | `~/Documents/BLUETELLA_9ID` | — | 9-ID's pyepics stack (Peco Myint) |

All measurements below were taken on **kouga** (no PV access) on 2026-09-13
against `dev` at `49d3ec0` and `main` at its rollback point `6afa83b`.

---

## 1. What 8-ID actually runs on

**The XPCS acquisition path has never used Bluesky — on either branch.** On
`main` at `6afa83b` (the "real Bluesky" tree), `ad_acq.py`, `master_plan.py`,
`dual_acq.py`, `dual_master_plan.py` and `tetramm_acq.py` contain **zero**
`bp.` / `bps.` / `yield from` / `RE(` constructs. All 5,304 lines are direct
ophyd puts and gets. The only Bluesky in the live science path is
`plans/align/scan_8id.py`.

### The "layers cause instability" claim does not survive the history

| | `main` (Bluesky era) | `dev` (Ophyd era) |
|---|---|---|
| except clauses in `id8_common` | 32 / 15,353 lines = **2.1 per kLOC** | 108 / 21,255 = **5.1 per kLOC** |
| align scan module | `scan_8id.py`: 2,044 lines, **0** except | `ophyd_scan.py`: 2,352 lines, **30** except |
| net of comments/blanks | 1,708 lines | 1,679 lines |
| writes its own CSV/SPEC output | no (databroker callbacks did it) | yes — 54 file-writing references |

Defensive code **more than doubled** after Bluesky was removed. The two scan
modules are the same size net of comments, but the ophyd one also absorbed the
data-capture responsibility the framework used to provide.

The stability fixes in the log are overwhelmingly **not** framework failures:
"Clear the Eiger out of Aborted before arming", "Eiger recovery: wait out
transient states", "Survive a dead 8idAlicat IOC", "Guard the QNW cells",
"scans hanging right after external mode". These are EPICS/IOC/detector
state-machine problems. They hit BLUETELLA identically — same IOCs.

### Where the criticism *is* correct

`AD_prime_plugin2` fired a real exposure on every session start, from inside
apstools, gated by a check (`AD_plugin_primed`) that can never pass on these
detectors. A layer nobody here wrote, doing an unrequested hardware write,
invisible to a grep of our own source. It crashed measurement G0209 on
2026-09-12. Fixed 2026-09-13 — see `STARTUP_HARDWARE_SAFETY.md`.

`registry.py` existing at all (325 hand-written lines, replacing
`apsbits.core.instrument_init`) is the same story: apsbits' loader let one dead
IOC starve every device declared after it.

**Conclusion.** The Bluesky-vs-not argument is retrospective; it was settled in
practice months ago. The substantive open question is whether to drop **ophyd**,
which is a bigger change and deserves to be argued on its own terms rather than
under a Bluesky banner.

---

## 2. Merging with 9-ID / BLUETELLA

### Already converged — two of three components

| Component | Owner | State |
|---|---|---|
| NeXus writer | Miaoqi Chu | **Done.** External package `nexus_xpcs_aps`; our in-tree writer deleted in `8b83c67`, retired to `utils/Archive` |
| Scan format + viewers | Peco Myint | **Mostly done.** 8-ID writes BLUETELLA's extended CSV and its `S` prefix since 2026-09-12; stock 9-ID viewers read our scans unmodified |
| Acquisition / device layer | Q. Zhang | **Not started.** This is the work |

BLUETELLA is already shaped for a merge: `tellalib.py` (1,660-line engine —
`MyMotor`, `MyDet`, `MyCounter`, `MyFilter`, `mv`/`umv`/`wm`,
`ascan`/`dscan`/`lscan`/`tscan`) contains **2** hardcoded 9-ID references.
`tellaconfig.py` contains **422**. The library/config split — the hard
architectural decision — has been made on both sides.

### Measured port cost, 8-ID → a BLUETELLA-style pyepics stack

| Layer | Size | ophyd coupling | Port cost |
|---|---|---|---|
| `plans/acquire/` | 5,304 lines | **zero.** 0 `stage()`, 0 `stage_sigs`, 0 `.set()`, 0 `.kind`, 0 `.read()`, 0 `.connected`. Just 120 `.get()` + 157 `.put()` on attribute chains | **Near-free.** `det.cam.acquire.put(1)` → `det.cam.put('Acquire', 1)`. Mechanical |
| Motors | 120 `EpicsMotor`, 88 `.move()`, 90 `.position` | `EpicsMotor` | **Near-free.** `MyMotor` wraps `epics.Motor` — pyepics' own motor-record class. Like-for-like |
| `devices/` | 2,990 lines | 499 `Component(` + 127 `Cpt(` + 12 `FormattedComponent`, 47 classes | **Bulk, low risk.** Becomes `tellaconfig`-style dicts. Data translation |
| Area detectors | — | 21 `ADComponent` + ophyd's plugin classes + apstools `ad_creator` | **The real engineering risk.** `MyDet` already drives eiger16m/pilatus/mmpad, so the pattern exists |
| `registry.py` | 325 lines | 17 `connected` / `wait_for_connection` | **Rebuild on pyepics.** This is `safe_make_devices` — why the beamline boots with dead IOCs. No BLUETELLA equivalent found |
| hklpy2 / `psic` | 4 files, 14 `add_reflection` | Bluesky `Movable`/`Readable` | **See §3.** BLUETELLA has zero diffractometer/UB support |

Headline: **the acquisition path is already written in a pyepics-compatible
dialect.** Most of the risk is gone before the port starts.

Remaining portability debt on our side: **134 hardcoded `8id` references in
`plans/`** (vs 49 in `devices/`). Plans are where beamline-specificity leaked.

### Caveats worth arguing about in January

* **Merging three single-author code bases yields three single-author
  subsystems.** Peco wrote BLUETELLA alone, Miaoqi the NeXus writer alone, and
  all 71 commits on `dev` since the rollback are by one author. Co-location does
  not create review. What creates review is a rule that every non-trivial change
  has a reviewer who is not its author — cheap, already working informally via
  PRs into BLUETELLA, and does not require a merge to start.
* **Release discipline is a prerequisite, not a follow-up.** `main` had to be
  rolled back two months to hand over support for *one* beamline, on an
  NFS-live checkout where an edit anywhere is instantly live on `pearl`. Shared
  files across two beamlines without tagged releases and per-beamline pinning
  converts a bus-factor problem into a blast-radius problem.
* **Dropping ophyd means BCDA can no longer help with anything.** Counterweight:
  this week's bug was in the layer they maintain, hidden where our grep could
  not reach.
* **Scope it to XPCS, not to APS.** Two beamlines, same detector class, same
  qmaps, same boost_corr/DM pipeline, same NeXus schema, same three people.
  That is a coherent domain library, not a framework competing with Bluesky.

### Zero-risk work that can happen before January

1. Parameterise the 134 hardcoded `8id` references in `plans/`. Worth doing
   regardless; also the best proxy for real portability.
2. Confirm `epics.Motor` error semantics match what our scans assume —
   `MyMotor.move()` defaults `wait=False`, and `epics.Motor.move` returns a
   status code rather than raising. Different contract from `EpicsMotor`.
3. Settle the hklpy2-without-ophyd question (§3 below — largely done).
4. Tagged releases + per-beamline pinning. Needed either way.

---

## 3. hklpy2: can we use only the math?

**Measured on 2026-09-13** in `8ide_bits_test`, hklpy2 with `hkl_soleil` 5.1.3.

### Yes for use — the solver runs with plain numbers, no ophyd device

```python
from hklpy2.backends.hkl_soleil import HklSolver
s = HklSolver(geometry="E6C", engine="hkl", mode="constant_phi_vertical")
s.real_axis_names    # ['mu','omega','chi','phi','gamma','delta']
s.pseudo_axis_names  # ['h','k','l']
```

Module-by-module, what is and is not entangled:

| Module | Imports |
|---|---|
| `backends/base.py`, `backends/no_op.py`, `backends/th_tth_q.py` | pure python |
| `backends/hkl_soleil.py`, `backends/hkl_soleil_utils.py` | `gi` (PyGObject → libhkl) only |
| `blocks/lattice.py`, `reflection.py`, `sample.py`, `constraints.py`, `configure.py` | pure python |
| `ops.py` (the core operations layer) | pure python |
| `diffract.py` | **bluesky + ophyd** |
| `incident.py` | **ophyd** |
| `misc.py`, `user.py`, `blocks/zone.py` | **bluesky** (and ophyd in `misc`) |

### No for install — `__init__.py` drags ophyd and bluesky in regardless

`hklpy2/__init__.py:60-61` does `from .diffract import DiffractometerBase` and
`creator`. Any `import hklpy2.<anything>` executes `__init__.py` first, so:

```
ophyd  pulled in by the import?  True
bluesky pulled in by the import? True
```

They must remain **installed**. They are never **used** at runtime by the math
path.

### What this does and does not mean for the environment

It does **not** mean keeping `8ide_bits_test` (2.4 GB, 281 packages). The
diffractometer math needs roughly: `hklpy2`, `pygobject` + `libhkl`, `ophyd`,
`bluesky`, `numpy`, `pyepics`. Not apsbits, databroker, tiled, apstools,
guarneri, or the queueserver stack. Two escape hatches if even that is
unacceptable:

* Import `hklpy2/backends/hkl_soleil.py` by file path via `importlib`, bypassing
  the package `__init__`. Fragile across hklpy2 releases.
* Call libhkl directly through `gi.repository.Hkl` — which is all
  `hkl_soleil.py` (~470 lines) does. No hklpy2 at all. This is the honest option
  if BLUETELLA gains diffractometer support.

**Recommendation:** keep hklpy2 and accept ophyd+bluesky as inert install-time
dependencies. They cost disk, not stability. Revisit only if a clean-env build
proves difficult.

---

## 4. YAML → CSV for `sample_info` / `measurement_info` before 2026-09-16

**Recommendation: do the trio-plan migration on Tue 2026-09-15. Do NOT do the
CSV conversion the same day.**

Context: beam drops Mon 2026-09-14 08:00; operation starts Wed 2026-09-16;
Tuesday 09-15 is the only working day, remote, with Sam and possibly Suresh.
Migrating the trio-detector acquisition plan from `~/ophyd` (`dev`) to
`~/bluesky` (`main`) is a hard requirement. The CSV conversion is discretionary.

### 4.1 The field-coverage gap — corrected 2026-09-13 after review

**An earlier draft of this note claimed the CSV could not express nine per-leg
fields. That was wrong; it is one.** Checked against the code:

| Field | Verdict |
|---|---|
| `start_timeout` | **Droppable.** `ad_acq.py:338` — `DEFAULT_START_TIMEOUT = 30.0`. The trio's `start_timeout: 30` is literally the default, so the line is redundant |
| `hdf_timeout` | **Droppable.** `ad_acq.py:341` — `DEFAULT_HDF_TIMEOUT = 300.0`. The trio's 600 doubles it; hardcode 600 as the default and the field disappears |
| `stop_timeout` | **Droppable.** Same pattern, never set in any live protocol |
| `num_segments` / `trigger_period` | **Not needed for the trio** — External Series only, and the trio is Internal Series / Internal / EPICS. But they are *measurement* parameters, not DM ones: `eiger4m_modes.py:335,339` put `num_frames * num_segments` into `hdf1.num_capture` and `num_segments` into `cam.num_triggers` |
| `geometry`, `motors`, `label`, `workflow_name` | **Confirmed unused.** No live leg in `measurement_info.yaml` or either `trio_measurement_info.yaml` sets them; every match is in a reference block |
| `acq_time`, `num_frames` | **Per-leg required** (`master_plan.py:83`) but shared in the CSV. Happens to be fine for the trio — all three legs are 1 s / 3000 frames. A latent trap, not a blocker |
| **`select_device`** | **The one that genuinely survives.** See below |

#### `select_device` is NOT neutered by `allow_motion: false`

`allow_motion: false` suppresses **motion**, not **registers**. In
`select_device()`, the registers loop runs *before* the `_motion_allowed(cfg)`
check:

```python
for reg_path, value in cfg.get("registers", {}).items():
    _resolve(reg_path).put(value)          # always runs
...
if _motion_allowed(cfg):                    # only motion is gated
```

And `device_position.yaml` gives each detector a *different* value of the same
register:

| entry | `allow_motion` | `registers` |
|---|---|---|
| `rigaku3M` | false | `softglue.enable_rigaku: '1'` |
| `rigaku3M_epics` | false | `softglue.enable_rigaku: '0'` |
| `eiger4M` | false | `softglue.enable_rigaku: '0'` |
| `lambda2M` | false | `softglue.enable_rigaku: '0'` |

So `select_device: yes` on the trio's `rigaku3M_epics` leg writes
`softglue.enable_rigaku = 0`. Drop it and the register keeps whatever the
previous run left — and a preceding `rigaku3M` run leaves it at `1`, i.e.
softglue enabled during a run that must not have it. This is exactly why
`master_plan.py:116` notes that `rigaku3M_epics` stopped being an alias of
`rigaku3M` on 2026-09-11: "it needs its own `softglue.enable_rigaku` value".

`trio_measurement_info.yaml` already documents this correctly in its own
reference block: *"With allow_motion: false on all three entries it currently
moves nothing either way; it still writes the entry's registers: block and sets
expt.det_name."*

**Net:** the CSV needs one new per-leg column, not nine. The technical objection
to the format is therefore much weaker than first stated. The reasons below are
what the recommendation now rests on.

### 4.2 It defeats the stated purpose of the `main` rollback

Commit `0376e6f`: *"Restores main to the tree of 6afa83b (2026-08-04). Requested
by the beamline scientist so that colleagues taking over support on 2026-09-16
work from the code base they know."*

The handover date **is** the Wednesday the operation starts. Handing colleagues
a rolled-back tree *and* a config file format nobody has seen contradicts the
reason the rollback happened.

### 4.3 It destroys operational knowledge at the worst moment

`measurement_info.yaml` and `measurement_info_2.yaml` carry roughly 60 lines of
comments that exist nowhere else:

* Internal Enable hung at `acq_period` 0.02 s — the Eiger serviced 12 of 100
  triggers (observed 2026-09-04).
* External Enable: one softglue pulse holds the shutter, so `acq_time` and
  `acq_period` both have a 0.1 s floor.
* External Series: `trigger_period` ≥ segment + 0.1 s, or the next pulse lands
  mid-segment, is dropped, and the acquisition hangs.
* ZDT floors: 2e-5 s for 2 bit, 4e-5 for 4 bit, 8e-5 for 8 bit.
* Folder names come from the attenuator **readback**, not `att_level` — 20 read
  back as 22, so expect `a0022`.
* `att10`'s stuck-filter history (E0160 ran at the wrong attenuation).

CSV has no comment mechanism; `#` rows are already structural in the drafted
format. This knowledge would be lost in the same week its author leaves for
months.

### 4.4 The "easier for non-programmers" argument inverts here

CSV-in-Excel is *more* dangerous for a non-Python-savvy editor, not less. Excel
silently reformats quoted comma lists (`"rigaku3M_epics, eiger4M"`), coerces
numbers to dates, strips leading zeros, and appends trailing commas. YAML fails
loudly at parse time; a mangled CSV parses fine and runs the wrong acquisition.

Also, the drafted `measurement_info.csv` is ambiguous where the YAML is explicit:

```csv
sample,protocol
1,test_protocol
2,
3,
1,test_protocol_2
,test_protocol_2
```

Blank `protocol` on rows 2-3 and blank `sample` on the last row have no stated
meaning — inherit from above? skip? The YAML `runs:` block states
`name` / `samples` / `protocols` / `repeats` outright.

### 4.5 The steelman, and the right version of this idea

The ergonomic win is real: a spreadsheet batch table is genuinely easier for a
beamline scientist than nested YAML, and Suresh should not have to count
indentation. The safe form is a **CSV importer, not a CSV runtime**: keep YAML
as the source of truth (comments and all), and add a converter that reads the
spreadsheet and *emits* the YAML, which `dry_run_measurement_info()` then
validates exactly as today. Runtime path unchanged, one new file, throwaway if
it does not work out. Build it in January with tests, not on 09-15.

### 4.6 Recommended plan for Tue 2026-09-15

1. Migrate the trio acquisition plan `dev` → `main` in `~/bluesky`. Nothing else.
2. `dry_run_trio_measurement_info()` and `dry_run_measurement_info(check_hardware=True)`
   until both are clean.
3. Walk Sam through the existing YAML as-is — that is the handover artifact.
4. Leave the CSV idea in this note for January.

### 4.7 Two discrepancies spotted in the attached files — check before Wednesday

Neither is caused by anything above; both look like drift.

* **`measurement_info.yaml` `att1` has `num_repeats: 441`** while `att40`/`att20`/
  `att10`/`att5`/`att2` all have `1`. At 3,000 frames × 1 s that is ~367 hours.
  `441 = 21²`, so it looks like a leftover from a 21×21 mesh plan.
* **The header comment is stale relative to `sample_info.yaml`.** The comment
  describes "a 1 x 1 mesh … total_pts is 1 … huber.x -0.2179 / huber.y 21.3250".
  `sample_info.yaml` `sample_1` now reads `inner_center: -0.142`,
  `outer_center: 21.18`, `inner_pts: 3`, `outer_pts: 3` — 9 points at a
  different position. With `sample_move: yes` and `position_reset: yes` the run
  will not do what the comment says.

---

## Provenance

Everything above is measured from the trees named at the top, on 2026-09-13, on
kouga. No claim here was verified against live PVs — kouga has no PV access.
The hklpy2 results in §3 are the exception: those were executed in
`8ide_bits_test` and are reproducible with the snippet shown.
