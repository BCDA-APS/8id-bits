# The scan viewer (BLUETELLA): checkout, environment, and what we changed

[← index](../../README.md) · [Viewing scans](../viewing-scans.md)

The maintainer's page for the viewer 8-ID plots its scans with: where it is
checked out, how its environment is built, what we changed in it, and why those
changes are upstreamable rather than a fork.

Day-to-day use — which viewer for which scan, the flags, the attenuation trap,
who to contact — is [Viewing scans](../viewing-scans.md). This page does not
repeat it.

**BLUETELLA is Peco Myint's (9-ID) code.** We run it and contribute back. The
string "8-ID" appears nowhere in either viewer, and that is deliberate: see
[The format decision](#the-format-decision).

## Where it lives

| | |
|---|---|
| checkout | `~/Documents/BLUETELLA_9ID` — i.e. `/home/beams10/8IDIUSER/Documents/BLUETELLA_9ID` |
| branch | `8id_test`; local `main` tracks Peco's `main` |
| remote | `github-qzhang234:AdvancedPhotonSource/BLUETELLA_9ID.git` |
| upstream | <https://github.com/AdvancedPhotonSource/BLUETELLA_9ID> — **private** |
| launcher, tracked | `~/bluesky/scripts/start_scanviewer.sh` — **edit this one** |
| launcher, on `$PATH` | `~/bin/start_scanviewer.sh`, a wrapper that does nothing but `exec` the tracked one (`BLUESKY_DIR` overrides the checkout) |
| conda environment | `bluetella_viewer` |
| launcher overrides | `VIEWER_DIR` (checkout), `VIEWER_ENV` (conda env), `EXPT_YML` (which `experiment.yml` the data folder is read from), `BLUESKY_DIR` (wrapper only) |

`github-qzhang234` is not a host — it is an alias in `~/.ssh/config` pointing at
`github.com` with `IdentityFile ~/.ssh/id_ed25519_qzhang234`. Clone with that
alias, or the private repo will refuse the default key.

`/home/beams/8IDIUSER` is a symlink to `/home/beams10/8IDIUSER`, so those two
paths are one directory and `$HOME` resolves there on every 8-ID workstation —
a path quoted either way is the same file. `/home/beams10` is NFS, so this
**one** checkout is what both kouga and pearl run. There is no deploy step: an
edit on kouga is live on pearl immediately. That matters because of the 8-ID
network split (`getent hosts kouga pearl`):

| host | network | GitHub | beamline PVs | `/gdata` | use it for |
|---|---|---|---|---|---|
| kouga | 164.54.116.55 | yes | no | **no** | clone, fetch, commit, push, run the tests, rebuild the conda env |
| pearl | 10.54.116.66 | no | yes | yes | **run the viewer** |

So do the git work on kouga and run the GUI on pearl (`ssh -Y`, since the viewer
needs an X display). Same files either way.

Before editing, `git pull`: CI fixes land on `8id_test` upstream while the pull
request is open, so the local branch is routinely behind.

Run `scripts/install-git-hooks.sh` once per clone. It sets
`core.hooksPath=.githooks`, whose `pre-push` refuses any push that updates
`refs/heads/main`. That hook plus a CI warning is the whole guard — GitHub
branch protection is not available for this private repo under the
organization's plan, as `CONTRIBUTING.md` says. Push feature branches only.

## What is in the repository

| file | what it is | used at 8-ID |
|---|---|---|
| `scanviewer.py` | 1-D scan browser and plotter: any column against any other, multi-scan overlay, Peak / COM / FWHM, dy/dx, "Follow latest" | **yes** |
| `meshviewer.py` | 2-D raster viewer: X vs Y coloured by Z, for `dmesh` / `mesh` | **yes**, via `--mesh` |
| `tests/test_scanviewer_csv.py` | 19 tests — CSV loading, prefix parsing, statistics | yes |
| `tests/test_meshviewer_csv.py` | 38 tests — extended-format loading, gridding, orientation | yes |
| `.github/workflows/ci.yml` | test matrix, plus a warning when `main` is updated without a PR | yes (it gates our PR) |
| `.githooks/pre-push`, `scripts/install-git-hooks.*` | local guard against pushing to `main` | yes, install it |
| `requirements-ci.txt` | pinned `pandas` / `scipy` / `matplotlib` for CI; NumPy is pinned per matrix leg | CI only |
| `tellalib.py` | 9-ID's instrument-control library (EPICS, scanning, file I/O) — their equivalent of our `id8_common` | no |
| `tellaconfig.py` | 9-ID's PV/device dictionary | no |
| `timeIOC.py` | small caproto soft IOC serving a `9idclock` timestamp PV | no |
| `user_config_template.py` | template for a 9-ID user's per-experiment scan hooks | no |

Both viewers are standalone — they import no EPICS and no `tellalib`, so they
run anywhere the data files are visible. Nothing at 8-ID imports the 9-ID
instrument code, and it is not meaningful here.

## The `bluetella_viewer` conda environment

The viewer has its own environment and **deliberately does not share one with
Bluesky**. It is a read-only plotting tool; it must never become the reason to
upgrade `pandas` or `matplotlib` inside an environment that runs the instrument.
`8id_bits`, `8ide_bits_test` and the rest are left alone — nothing here installs
into them or imports from them. Do not point `VIEWER_ENV` at one.

Five dependencies: `pandas`, `numpy`, `scipy`, `matplotlib`, `tkinter`.

Rebuild it from a **164.\*** host such as kouga — conda needs the internet and
the 10.\* beamline machines do not have it. Conda environments live under the
shared home, so pearl picks the result up with no further step:

```bash
ssh 8idiuser@kouga
conda create -y -n bluetella_viewer python=3.11 pandas numpy scipy matplotlib-base tk
```

`matplotlib-base` rather than `matplotlib` keeps Qt out of an environment that
only ever uses the Tk backend.

Fallback if the environment is missing mid-experiment: the **system `python3`**
on both hosts has all five, so the viewer runs without conda at all. Pass the
prefix by hand, because only the launcher supplies it:

```bash
python3 ~/Documents/BLUETELLA_9ID/scanviewer.py --dir <folder> --scan-prefix A
python3 ~/Documents/BLUETELLA_9ID/meshviewer.py --dir <folder> --scan-prefix A
```

## Running the tests

```bash
cd ~/Documents/BLUETELLA_9ID && python3 -m unittest discover -s tests -q
```

57 tests, in well under a second. They are pure file-format tests: no EPICS,
no hardware, no display, no `/gdata` — so they run on kouga, and they are the
first thing to run after touching either viewer.

Peco's CI runs the same command on two legs, Python 3.10 / NumPy 1.26.4 and
Python 3.12 / NumPy 2.4.6, after `python -m compileall`. Locally the same suite
has been run on NumPy 1.24.4, 1.26.4 and 2.4.6. Two NumPy majors is not
pedantry — see change 1 below.

## How the viewer decides which scan is "current"

Entirely from **file names on disk**; it never reads EPICS. With "Follow latest"
ticked, `_poll_follow()` runs once a second (`POLL_INTERVAL_MS = 1000`),
re-lists the directory, parses the leading number out of each name, sorts by
`(scan_number, mtime)` and follows the highest. A name that does not match is
ignored entirely. Before switching, it proves the new file actually loads and
caches, so a half-written file cannot blank an already-good plot; only the
followed file is re-read on each tick.

The pattern comes from the prefix the viewer was given. BLUETELLA's own default
is glob `S*.csv` / `S*.json` with number regex `^S(\d+)_`; `--scan-prefix A`
makes that `A*` and `^A(\d+)_`. 8-ID names its scans
`A0113_Test_a1010041.csv` — from `gen_folder_prefix()` in
[`plans/acquire/acq_helpers.py`](../../src/id8_common/plans/acquire/acq_helpers.py),
`{header}{measurement_num:04d}_{sample_name}_a{attenuation:04d}` — so
`start_scanviewer.sh` adds `--scan-prefix A` unless you pass your own. **That is
the only local adjustment 8-ID needs.** Without it the browser globs `S*` and an
8-ID folder lists as empty.

One prefix is active at a time; the viewers never glob two. `scanviewer.py` also
takes `--meta-prefix` (default `M`) for its Meta source mode, which 8-ID does
not use; `meshviewer.py` has no Meta mode and so only `--scan-prefix`.

Two consequences worth knowing:

* The viewer only ever watches the one directory it was given. After a cycle or
  experiment rollover, restart it — the launcher rebuilds the folder from
  `configs/experiment.yml` each time.
* It assumes measurement numbers increase. They do: `gen_folder_prefix()`
  increments a shared EPICS counter. But if that counter is reset, or an older
  file with a higher number is copied in, "Follow latest" will sit on the wrong
  file.

## The format decision

**8-ID writes BLUETELLA's own extended CSV format.** Adopted 2026-09-07 in
[`plans/align/scan_csv.py`](../../src/id8_common/plans/align/scan_csv.py): a
`#label,value` preamble, `#DATA`, the column-name line, the rows, a **bare**
`#END`, then `#exit_status` and `#points_written`.

This is the load-bearing decision on this page. Because the file is already in
his format, pristine upstream `scanviewer.py` parses an 8-ID scan with **no
8-ID-specific code at all** — which is what makes our viewer changes
upstreamable instead of a fork we maintain for ever. Only two things changed on
our side to get there: header labels gained a leading `#`, and `#END,<status>,<n>`
became a bare `#END` with the outcome on the lines after it. No header line was
dropped, moved, or sent to a sidecar; all ~40 template lines are still in the
one file, in order.

What a viewer relies on:

| line | meaning |
|---|---|
| bare `#END` | the terminator, and nothing else on the line. The old `#END,success,21` form sat *inside* the region pandas reads, so it parsed as a phantom data point whose motor value was the string `success` — flipping that column to text and taking COM and the derivative view with it. A 21-point file read as 22 rows |
| `#exit_status`, `#points_written` | after `#END`, outside the parsed region. `points_written` is what was written; the header's `num_points` is what was asked for, and an aborted scan must not overwrite it |
| `#shape`, `#num1`, `#num2` | written **only** by a raster. A `d2scan`/`a2scan` writes none, which is how a viewer knows not to grid it |
| `#motor`, `#motor2` | the scan's own motor names — the grid check compares them against the chosen X and Y columns |
| `#motor1_start/stop`, `#motor2_start/stop` | the **commanded** grid, not readbacks |

Which lines appear is controlled by
[`configs/scan_csv_template.yml`](../../src/id8_common/configs/scan_csv_template.yml)
(40 header entries, 13 column entries today), not by Python. A scan also refuses
to write into a file that already exists (`FileExistsError`) rather than append
a second scan to it.

### Approaches that were tried and deleted

Nothing in either viewer implements these any more. Do not re-invent them.

| deleted approach | why it went |
|---|---|
| teach the viewer 8-ID's own dialect — match `#END` as the first comma-separated field to cope with `#END,success,21` | the trailing fields were the bug, not the reader. Fixed on the **writer** side instead |
| hard-code `A` alongside `S`/`M` in the filename regex | replaced by the generic `--scan-prefix`, which no site has to patch |
| a 9-ID/8-ID "format" layer, with a `--format` flag and a sniffer to guess the dialect | deleted once 8-ID adopted the 9-ID format. There is one dialect, and neither viewer has any notion of a second |

## What we contributed

Four separable changes, none containing beamline-specific code. They can be
merged independently, in any order.

| # | file | change | kind |
|---|---|---|---|
| 1 | `scanviewer.py` | `np.trapz` → a version-agnostic `_trapezoid`, resolved once at import | **9-ID bug fix** |
| 2 | `scanviewer.py`, `meshviewer.py` | `--scan-prefix` / `--meta-prefix`, defaulting to `S` / `M` | new generic option |
| 3 | `meshviewer.py` | read the extended `#DATA` / `#END` format | **9-ID bug fix** |
| 4 | `meshviewer.py` | draw a raster that declares its shape as a filled image | new generic feature |

1. **`np.trapz` → `_trapezoid`.** NumPy 2.0 renamed `np.trapz` to
   `np.trapezoid` and kept the old name as a deprecated alias; **2.4 removed
   `np.trapz` outright**, so the COM line in `compute_scan_stats()` raised
   `AttributeError` and took the whole Peak/COM/FWHM panel down with it.
   `np.trapezoid` does not exist before 2.0, so neither spelling alone covers
   the interpreters in use. Nothing to do with 8-ID — it surfaced only because
   the viewer got its own, newer environment.

2. **`--scan-prefix` / `--meta-prefix`.** The prefixes were hard-coded in the
   glob and in the number regex. They are arguments now, defaulting to `S` and
   `M`, so 9-ID behaviour is byte-for-byte unchanged and any site with another
   convention configures rather than patches. This is the *whole* of what 8-ID
   needs on the viewer side.

3. **`meshviewer.py` reads the extended format.** It previously handled plain
   CSVs only — it could not open 9-ID's *own* extended `#DATA`/`#END` files,
   which `scanviewer.py` has always read. A 9-ID fix that happens to also let it
   open ours.

4. **Declared rasters drawn as a grid.** A file carrying a shape is rendered
   with `imshow`, with a "Draw raster as grid" checkbox to turn it off;
   anything without a shape keeps the original scatter, unchanged. Cells are
   binned by nearest **commanded** position, never by reshaping rows — a real
   stage repeats only to its retry deadband, so a 5×5 raster has 25 distinct
   readbacks and a plain reshape corrupts silently on an aborted scan or a
   dropped point. Unmeasured cells stay NaN and render blank, never zero. X/Y/Z
   also now default to the first, second and last column rather than all three
   to the first.

Diff against `origin/main`: `scanviewer.py` +47/−7, `meshviewer.py` +322/−23,
plus 554 lines of new tests.

### Peco's review, and what changed for it

All four corrections are in commit `21c26b5`.

| he asked for | what it means in the code |
|---|---|
| literal filename prefixes, not a character class | `_number_re()` builds `^` + `re.escape(prefix)` + `(\d+)_`. A two-letter prefix `AB` now means files beginning `AB`, not "A or B" |
| grid only when X and Y are the two motors the scan declared | `grid_orientation()` accepts `X=motor1, Y=motor2` or the swap, and returns `None` otherwise. Previously a declared shape was enough, so plotting `elapsed_time` against `motor2` produced a plausible, meaningless filled map |
| bulk cell assignment, not a per-point `argmin` loop | `_nearest_index()` uses `np.searchsorted`, flipping and re-mapping a descending axis (`dscan(m, +1, -1, …)` is legitimate) |
| regression tests | the suite went to 57 tests |

### Where it stands

**Pull request #1 is OPEN, not merged:**
<https://github.com/AdvancedPhotonSource/BLUETELLA_9ID/pull/1>

Peco has approved the four fixes in a comment; the PR needs @ennogra's approval
or beamline testing before it merges. His CI passes on both legs.

The two 8-ID site documents that were on that branch — `8ID_SCANVIEWER.md` and
`USING_SCANVIEWER.md` — were deleted at his request and ignored in that clone,
since his repository should describe the package, not one beamline's
deployment. **This page and [Viewing scans](../viewing-scans.md) are what became
of them.** Do not put 8-ID documentation back in his repository.

## Hardware verification, 2026-09-08

8-ID-E, on pearl: all six scan types, against **two real motors** — `huber.nu`
and `huber.delta` — with the Lambda2M as detector. Files are in
`/gdata/dm/8ID/8IDE/2026-3/comm202609/data/bluesky/`:

| file | scan | what it exercises |
|---|---|---|
| `A0114` | `dscan` | 1-D |
| `A0115` | `ascan` | 1-D |
| `A0116` | `d2scan` | two-motor trajectory, no `#shape` → stays a scatter |
| `A0117` | `a2scan` | two-motor trajectory, no `#shape` → stays a scatter |
| `A0118` | `dmesh` | 5×5 raster → grid |
| `A0119` | `mesh` | 4×4 raster → grid |

Open the first four with `start_scanviewer.sh`, the last two with
`start_scanviewer.sh --mesh`. An earlier, equivalent run on 2026-09-07 left
`A0107`–`A0113` in the same folder; `A0113` also ran with `save_img=1`, so it
has a `.h5` beside it.

Those runs were taken under heavy attenuation, so every `lambda2M_stats*` column
reads zero — pick `tetramm1_sum_all`. This is the single most common "the viewer
is broken" report; [Viewing scans](../viewing-scans.md) explains it.

Ctrl+C was exercised on all three detector branches in the 2026-09-07 run:
motors return to where the scan started, the shutter is closed and confirmed against its
readback, acquisition stops, and the file closes with `#END` /
`#exit_status,aborted` / a short `#points_written`. Such a file is still a valid
scan file — the 1-D viewer plots the points that were taken, and the mesh viewer
places them in the right cells and leaves the rest blank
(`test_partial_raster_places_cells_where_the_full_map_has_them` and
`test_a_partial_large_raster_leaves_the_rest_blank`).

## Relationship to `specr_py`

`~/Documents/specr_py` (<https://github.com/qzhang234/specr_py>, launched by
`~/bin/start_specr_py.sh`, conda environment `specr_py`) is 8-ID's own PyQt5
viewer for the same scan CSVs. It is not a competitor so much as a different
mechanism, and they fail in different places:

| | `scanviewer.py` (BLUETELLA) | `specr.py` (specr_py) |
|---|---|---|
| finds the live scan by | file names in one directory | EPICS: `caget` on `8ideSoft:StrReg8` (bare scan name) plus `StrReg3`/`StrReg1`/`StrReg4` for mount point, cycle, experiment |
| needs the beamline network | no — runs on kouga, or a laptop, on copied files | yes: `caget` and a route to the soft IOC (it says so and falls back to File ▸ Open) |
| follows a cycle/experiment change | no — restart it | yes, it re-reads the registers |
| can be fooled by a reset counter | yes, it follows the highest number | no, it is told the exact name |
| header preamble | skipped; it plots the data table | parsed (`csvfile.py`, `CsvScan`: `meta`, `positions`, `comments`), and the originating command goes in the plot title |
| extras | multi-scan overlay, dy/dx, Peak/COM/FWHM, 2-D raster viewer | SPEC-file support for pre-2026-08 data |
| toolkit | tkinter | PyQt5 |

Neither launcher sources `dm.setup.sh`: it overwrites `LD_LIBRARY_PATH`, which
neither GUI toolkit wants, and a viewer needs no Data Management.

**⚠ `specr_py` has not been updated for the current CSV format** — verified
2026-09-08 against a file in the current format. Its reader treats the first
line beginning with `#` as the `#DATA` marker, so now that every header label
carries one it takes the *second* header line as the column names and the rest
of the preamble as data; it also still expects `#END,<status>,<n>`. Reproduce it
against any current scan `.csv` (on pearl, or a copy of one anywhere):

```bash
cd ~/Documents/specr_py && python3 -c "
import csvfile
d = csvfile.CsvDataFile('<folder>/<a current scan>.csv'); d.refresh()
s = d[-1]; print(s.labels[:3], len(s.rows), repr(s.command), s.exit_status)"
```

A current file yields `labels` that are the two fields of the `#h5_file` line, a
row count that counts preamble lines as points, an empty `command`, no motor
positions, and `exit_status` of `unknown`. Fixing it means teaching
`csvfile.py` the `#`-prefixed preamble and the bare `#END` — a small change, in
*our* repository, on the reader side only. Until then, BLUETELLA is the viewer
that reads current scans.

## See also

* [Viewing scans](../viewing-scans.md) — how to use it, and when to contact Peco
* [Ophyd scans and the CSV file template](../../src/id8_common/plans/align/OPHYD_SCAN.md)
  — the scans themselves and what the `.csv` contains
* [`plans/align/scan_csv.py`](../../src/id8_common/plans/align/scan_csv.py) — the writer
* [`configs/scan_csv_template.yml`](../../src/id8_common/configs/scan_csv_template.yml) — which lines and columns a scan writes
* `scripts/start_scanviewer.sh` — the tracked launcher
