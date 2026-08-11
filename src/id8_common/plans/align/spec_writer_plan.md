# Hand-written SPEC file writer for live scan visualization

Status: proposed, not implemented. Written 2026-08-10.

Supersedes [live_scan_feed_plan.md](src/id8_common/plans/align/live_scan_feed_plan.md) (the UDP-feed
proposal). Cost analysis of the file-per-point alternative is in
[live_viz_overhead_report.md](src/id8_common/plans/align/live_viz_overhead_report.md).

## Context

Alignment scans in [scan_8id.py](src/id8_common/plans/align/scan_8id.py) need a live curve
(ROI/stats intensity vs motor position). Reading the AreaDetector `.h5` from a second process while
the IOC writes it is the hazard that started this; writing a **SPEC-format text file, one row per
scan point**, avoids it entirely — the file is append-only, closed between writes, and a reader can
never collide with the writer.

SPEC is a standard format at APS and other US light sources, so a polling reader is not a bespoke
protocol: `spec2nexus`, PyMca, and silx already parse it.

**Decision:** write the SPEC file directly from the plan with plain Python + Ophyd. No `apstools`
dependency, no RunEngine subscriber. This was chosen deliberately for full control after reviewing
the alternative.

**The alternative, for the record.** `apstools.callbacks.SpecWriterCallback2` already does
incremental per-point writing — verified in the installed source: `event()` calls
`write_scan_data_row(doc)` on every event, and `_write_lines_(mode="a+")` opens/writes/**closes** per
row. Wiring it is ~4 lines in `startup.py` plus copying the existing
[demo_spec_callback.py](src/legacy/id8_e/callbacks/demo_spec_callback.py). Replaying real `dscan`
documents through it produces a valid file with all the needed columns. If the hand-written writer
becomes a maintenance burden, that is the fallback.

Two related facts worth recording:

- `SPEC_DATA_FILES: ENABLE: true` in [iconfig.yml](src/id8_common/configs/iconfig.yml#L50) is **dead
  config today** — nothing in [startup.py](src/id8_common/startup.py) reads it, and the legacy block
  at [legacy/id8_e/startup.py:76](src/legacy/id8_e/startup.py#L76) is commented out. No SPEC file is
  currently being written.
- The `fw_specfile_example` documentation page uses `SpecWriterCallback` **v1**, which batches at
  `stop()`. It does not demonstrate the incremental path. Only v2 writes per point.

---

## Format contract

These come from reading the `spec2nexus` parser, not from convention. Getting any of them wrong
produces a file that looks fine and silently mis-parses.

| Rule | Why |
|---|---|
| `#L` labels separated by **two or more spaces** | `split_column_labels()` splits on a `\s{2,}` regex. Single-space only works as a fallback when `#N` exactly matches the column count. |
| `#N <n>` must equal the number of `#L` labels | It is the fallback's only guard; a wrong value corrupts the fallback path too. |
| `#S ` requires the **trailing space** | Scan detection is `line.startswith("#S ")` (spec.py:159, 386). |
| `#S <int> <command>` | The integer is `scanNum`; the rest is the command text. |
| Blank line before each `#S` | Sections are delimited by `#E`/`#F`/`#S`. |
| Every data cell must `float()`-parse | Parser does `buf[self.L[col]] = float(val)`. No strings, no empty cells, no `None`. |
| `#F` or `#S` must appear | That is what marks the file as SPEC (spec.py:150). |

Number formatting: use `f"{v:.12g}"`. Guard non-finite values — replace `nan`/`inf` with `0` and note
it in a `#C` line rather than emitting a token that some parsers reject.

---

## Changes

### 1. New module: `src/id8_common/utils/spec_writer.py`

A plain class. Standard library only — no `apstools`, no `bluesky` imports.

```python
class SpecFile:
    def __init__(self, path): ...          # remembers path; does not create the file
    def start_scan(self, scan_num, command, labels, comments=None): ...
    def add_point(self, values): ...       # one row; values must align with labels
    def end_scan(self, status="success"): ...
```

Behaviour:

- **File header** (`#F`, `#E`, `#D`, `#C`) written once, only when the file does not yet exist.
  Appending to an existing file skips straight to the scan header, so one file accumulates scans.
- **`start_scan`** writes a blank line, then `#S`, `#D`, `#C uid/comments`, `#N`, `#L`.
- **`add_point`** formats one row and appends it.
- **Every write is `open(path, "a")` → one complete line ending in `\n` → `close()`.** A single
  append-mode write of a short line is atomic for regular files, so a polling reader never sees a
  torn row. This matches what `SpecWriterCallback2` does and is the property the whole design rests
  on.
- No exception may escape into the RunEngine. Wrap each public method so a full disk or a bad path
  degrades to a logged warning, never a failed scan.

### 2. Call it from `dscan` first — prove it before rolling out

The `bps.create` / `bps.read` × 4 / `bps.save` block is copy-pasted at **16 sites** across
`dscan`/`dmesh`/`mesh`/`d2scan`/`ascan`/`a2scan`/`dscan_auto` × the lambda and eiger branches
(verified: `grep -c 'bps.create("primary")'`). Do **not** edit all 16 up front.

**Phase 1 — `dscan` eiger + lambda branches only.** Add `spec.start_scan(...)` next to the existing
`md` construction, `spec.add_point(...)` inside the per-step alongside `bps.save()`, and
`spec.end_scan()` in the existing `finally:`. Confirm live plotting works against a real scan.

**Phase 2 — collapse to one helper and roll out.** Only once Phase 1 is proven:

```python
def _read_point(motors, det, spec=None):
    yield from bps.create("primary")
    for m in motors:
        yield from bps.read(m)
    readings = [(yield from bps.read(getattr(det, f"stats{n}"))) for n in (1, 2, 3)]
    yield from bps.save()
    if spec is not None:
        spec.add_point([m.position for m in motors] + [...])   # reuse readings, no extra CA calls
```

Reuse the values `bps.read()` already returned rather than issuing fresh `.get()` calls — that avoids
4 extra CA round-trips per point. Confirm the reading-dict shape when wiring it up.

### 3. File location and scan numbering

- **Path:** `<...>/data/bluesky/<experiment_name>.dat`, alongside the `.h5`. Reuse the path
  construction already in `save_images` at
  [scan_8id.py:76](src/id8_common/plans/align/scan_8id.py#L76) rather than rebuilding it, and
  `pv_registers.experiment_name` for the stem.
- **Scan number:** derive from the `folder_prefix` the plan already has. `gen_folder_prefix()`
  ([ad_acq.py:100](src/id8_common/plans/acquire/ad_acq.py#L100)) returns e.g. `A0092_G10_a0007` and
  increments `pv_registers.measurement_num` as a side effect — so parse the number out of the prefix
  rather than reading the register again. This makes `#S 92` line up with `A0092_*.h5`.
- **Link to the images:** emit `#C image_file = A0092_G10_a0007.h5` in the scan header. This is what
  lets a reader get from a SPEC scan to the frames.

### 4. Reader / plotter (separate process, separate file)

Poll the `.dat`, track file size, parse only appended bytes. Two rules that matter:

- **Ignore a trailing incomplete line** — only parse up to the last `\n`.
- Reset the plot on a new `#S`; read axis labels from `#L`.

Polling every 1 s is fine — the file is small and the read is microseconds. 5 s also works; a
41-point scan at ~2 s/point runs ~80 s.

---

## Known issue that will bite during implementation

[scan_8id.py:145](src/id8_common/plans/align/scan_8id.py#L145) and
[scan_8id.py:753](src/id8_common/plans/align/scan_8id.py#L753) are **byte-identical `dscan`
definitions**. Line 753 shadows line 145, so edits to the first have no effect. Delete lines
145–328 before starting, or Phase 1 will appear to do nothing.

---

## Verification

1. **Format compliance — the main risk of hand-writing.** Parse the output with the installed
   parser and confirm the columns round-trip:
   ```python
   from spec2nexus.spec import SpecDataFile
   s = SpecDataFile(path).getScan(N); s.interpret()
   assert s.L == expected_labels and len(s.data[s.L[0]]) == num_points
   ```
   This is the single most important check. Run it before trusting anything else.
2. **Differential test against the reference implementation.** Run the same document stream through
   `SpecWriterCallback2` and diff against the hand-written output. Divergence is not automatically a
   bug, but each difference should be understood rather than discovered later by a downstream tool.
3. **Atomicity.** Run a scan while a reader polls at 10 Hz. It must never parse a partial row and
   never raise. Assert the reader's point count matches `num_points` exactly at the end.
4. **Append across scans.** Run three scans into one file; confirm `getScanNumbers()` returns three
   and that scan numbers match the `A####` prefixes of the corresponding `.h5` files.
5. **Writer cannot break a scan.** Point the writer at an unwritable path and run a `dscan`. It must
   complete normally with a logged warning.
6. **No file contention.** `lsof` on the reader during a scan must show no handle on any `.h5`, and
   no persistent handle on the `.dat`.
7. **Phase 2 is behaviour-preserving.** After the refactor, confirm Tiled event counts and column
   names are unchanged from a pre-refactor run of each plan type.

There is no unit-test suite in this repo — the `test_*.py` files under `plans/` and `user/` are
beamline scripts. Verification is by running the plans.

---

## Not doing

- Wiring `SpecWriterCallback2` as an RE subscriber. Considered and declined; kept above as the
  documented fallback.
- Removing the existing `bps.create/read/save`. It stays; the SPEC write is additional.
- The `#MD iconfig` bloat, detector-vs-unassigned classification, and empty `#O`/`#P` lines. These
  are artifacts of the apstools writer and mostly do not apply to a hand-written one; revisit only
  if a downstream tool complains.
