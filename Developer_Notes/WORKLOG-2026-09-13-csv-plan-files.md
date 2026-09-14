# Work log — 2026-09-13: reading measurement plans from CSV

Written so the thread survives losing the chat it came from. Companion to
`EVAL-2026-09-13-stack-consolidation-and-config-format.md` §4, which argued the
case; this is what was actually built.

**Status: built and tested off-beamline on kouga. Not used by anything yet.**
Every session today still reads YAML and is completely unaffected — see
*Safety* below, which is the part to check before believing that.

## Why

Users who set up runs are beamline scientists, not programmers. YAML
indentation is a real barrier and a *silent* one: a key indented two spaces too
far is still valid YAML, it just lands on the wrong parent. Demonstrated:

```
detectors:                   detectors:
  - device: eiger4M            - device: eiger4M
    acq_time: 1                  acq_time: 1
    num_frames: 3000         num_frames: 3000      <- two spaces left
```

Both parse. The left leg has `['acq_time','device','num_frames']`; the right leg
has `['acq_time','device']` and `num_frames` silently became a protocol key.

## What was added

| File | Lines | Role |
|---|---|---|
| `src/id8_common/plans/acquire/plan_csv.py` | new | CSV → the dict `yaml.safe_load` would have returned |
| `src/id8_common/plans/acquire/validators.py` | +18 | `read_yaml()` dispatches on suffix |
| `scripts/check_plan_csv.py` | new | equivalence checker, runnable anywhere |
| `Developer_Notes/csv_examples/` | new | 5 example CSVs + 2 hand-written expected YAMLs |

### The contract

```
read_plan_csv(x.csv)  ==  yaml.safe_load(x.yaml)
```

The converter does **shape and type only**. It returns the same nested dict, so
`expand_measurements()`, `REQUIRED_LEG_FIELDS`, `validate_acq_time` and
`validate_qmap_exists` all run unchanged and do every semantic check. Nothing
downstream knows which format was on disk. **No file is ever written** — there
is no generated YAML to go stale.

### Where it hooks

`read_yaml` in `validators.py` has exactly four callers, all in `master_plan.py`
(1060/1061 in `run_measurement_info`, 1104/1105 in `dry_run_measurement_info`).
Hooking it there means **the dry run validates precisely the structure the run
will use** — no second parser in the critical path.

### How to switch a run to CSV

Point `sample_info_file` / `measurement_info_file` in `configs/experiment.yml`
at the `.csv`. That is the whole change: no code edit, per-experiment, and
reversible by editing the line back. Both formats work side by side
indefinitely — Sam can run YAML while Suresh runs CSV.

## Safety — why this cannot disturb a running beamline

1. **The YAML path is byte-identical.** Verified against all three live
   pope202609 files: `read_yaml(f) == yaml.safe_load(open(f))` → True for
   `measurement_info.yaml`, `sample_info.yaml`, `trio_measurement_info.yaml`.
2. **The `plan_csv` import is INSIDE the `.csv` branch, deliberately.** Verified
   by replacing `plan_csv.py` with a syntax error: YAML still read fine, only
   the CSV path raised. A bug in the new module cannot reach a YAML session.
3. **Nothing else changed.** No device, plan, config or `user_plans` file was
   touched. No EPICS, no registry, no hardware in any new code path.
4. **Nothing opts in by itself.** Until `experiment.yml` names a `.csv`, none of
   this executes.
5. Already-running sessions are unaffected regardless — Python caches imported
   modules, so an edit to `validators.py` reaches only sessions started after it.

## Format

A row whose first cell starts with `#` is one of three things, told apart by
name and shape:

* a name in `KNOWN_SECTIONS` (`DEFAULTS SAMPLES DETECTOR MEASUREMENT BATCH
  REPEAT LOOP_ORDER`) opens a section; `#REPEAT,1` also carries its value;
* not a known section **but has a value** → a data row whose key is written with
  a `#`, e.g. `#device,"rigaku3M_epics, eiger4M"`. The `#` is stripped. The
  drafted files do this, and a parser that treated every `#` row as a section
  would silently drop the detectors;
* anything else → **a free-text comment, ignored.** Put notes in the spreadsheet
  this way; they cost nothing.

Wholly blank rows are dropped, so Excel's trailing empties are harmless.

`#DETECTOR` values are comma-separated, one per leg; a single value broadcasts
to all legs, any other length mismatch is an error rather than a guess.
`#MEASUREMENT` keys are split by `LEG_FIELDS`: a per-leg field (`acq_time`,
`num_frames`, …) is copied into every leg because `REQUIRED_LEG_FIELDS` demands
it per leg; the rest stay protocol-level. `#BATCH` in a protocol file yields one
protocol per row — `test_protocol` + `att_level` 1,2,3 → `test_protocol_att1/2/3`.
Always emits the parallel `detectors:` shape, even for one detector.

## How it was tested

**Equivalence (the real test).** `scripts/check_plan_csv.py` deep-compares the
converter's dict against a YAML **a human wrote from the format spec**. A golden
file dumped from the converter would only prove the converter equals itself.

```
$ python scripts/check_plan_csv.py
=== sample_info.csv vs sample_info.expected.yaml ===        MATCH
=== measurement_info.csv vs measurement_info.expected.yaml === MATCH
2/2 pair(s) matched
```

Runs anywhere — no EPICS, no `/gdata`, no beam. Exits non-zero on mismatch, so
it can gate a commit. Point it at your own pair:
`python scripts/check_plan_csv.py mine.csv mine.yaml`.

It compares through `read_yaml`, not `read_plan_csv`, so it also proves the
suffix dispatch routes the way a real run will.

**Dry-run parity.** `dry_run_measurement_info()` on the CSV pair reaches exactly
the same point the YAML does on kouga — `validate_qmap_exists` failing because
`/gdata` is not mounted here. Everything before that validated identically, and
`Total measurements planned: 15` matches the expected 5 rows × 3 att levels. On
pearl or amber this goes all the way through.

**Still to do on the beamline:** `dry_run_measurement_info(check_hardware=True)`
against a CSV, then one short real acquisition.

## Findings the dry run produced

Both are about the *drafted CSVs*, not the converter, and both were caught by
checks that already existed — which is the architecture working as intended.

1. **`sample_name` is required and was missing.** `validators.require_fields`
   raises `missing sample field: 'sample_name'`. The column is now in the
   example. Any real `sample_info.csv` needs it.
2. **The drafted `measurement_info.csv` is rejected as written.** Its last row
   has a blank `sample`, which forward-fills to 1 — giving sample 1 /
   `test_protocol_2` twice. `check_duplicate_assignments` refuses: *"Use
   runs[].repeats instead of duplicating run blocks."* Correct behaviour, and it
   means **the blank-cell rule needs a decision** (see below).
   `measurement_info_valid.csv` is the same file with that row naming sample 2.

## Open decisions before this is used in anger

* **Blank cells forward-fill from the row above.** That is spreadsheet
  intuition, and every fill is printed to the terminal so a stray blank is
  visible rather than silent. But it is a guess about intent — confirm it is
  what users expect, or require every cell.
* **Run-block names** are `<protocol>_sample<N>` and may repeat. Harmless:
  `expand_run_block` uses only `samples`, `protocols`, `repeats`, `loop_order`;
  `name` is cosmetic.
* **Protocol naming** strips a trailing `_level`, so `att_level: 20` →
  `test_protocol_att20`, matching the existing `att20` convention. Any other
  batch column keeps its full key.
* **Comments.** Durable reference material (the 0.02 s Internal Enable hang, the
  0.1 s shutter floor, the `trigger_period` derivation, the ZDT floors) should
  move to `docs/running-measurements.md`, which the YAML headers already name as
  its home — it is currently duplicated across three YAML files and drifting.
  Short per-row notes can live as `#` comment rows in the CSV.

## Recommended sequencing

Unchanged from the EVAL note: on Tue 2026-09-15, **migrate the trio acquisition
plan to `~/bluesky` first and finish it.** This CSV work is additive and
reversible and can be tried afterwards, or not at all. Do not make CSV the only
format before the 09-16 handover.

Note this was built in `~/ophyd` (`dev`). `~/bluesky` (`main`) does not have it;
porting is three files and one 18-line edit, but do it after the trio migration.
