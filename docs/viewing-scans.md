# Viewing scans

[← index](../README.md)

Every scan in [`ophyd_scan.py`](../src/id8_common/plans/align/ophyd_scan.py)
writes one `.csv` next to its `.h5`, and BLUETELLA plots it — live, while the
scan is still running.

```bash
start_scanviewer.sh              # 1-D: dscan, ascan, d2scan, a2scan
start_scanviewer.sh --mesh       # 2-D: dmesh, mesh, drawn as a grid
```

That is the whole thing for normal use. With no arguments it opens **this
experiment's** scan folder, worked out from `configs/experiment.yml` the same
way the scans work it out, so it follows a cycle or experiment rollover with no
edit.

## Which viewer

Follow the scan, not the number of motors.

| you ran | viewer | why |
|---|---|---|
| `dscan`, `ascan` | 1-D | one axis |
| `d2scan`, `a2scan` | **1-D** | two motors, but along **one line** — there is no 2-D field to colour. Pick either motor as X |
| `dmesh`, `mesh` | 2-D (`--mesh`) | a real raster: `motor1` slow, `motor2` fast |

The 2-D viewer decides from the **file**, not the scan's name: it draws a grid
only when the header declares a shape (`#shape,5x5`). A `d2scan` writes no
shape, so it stays a scatter even if you open it there.

## Options

| flag | what it does |
|---|---|
| `--mesh` | the 2-D viewer instead of the 1-D one |
| `--dir <folder>` | some other folder; suppresses the automatic lookup |
| `--scan-prefix <letter>` | the launcher already passes `A`; override only if your files start with something else |
| `--help` | the viewer's own options |

Tick **Follow latest** to track the running scan. Both viewers are read-only:
no EPICS, no ophyd, no Data Management. They only ever read files.

## The first thing that confuses people

**Under high attenuation the `lambda2M_stats*` columns are all zero.** The plot
is then a flat line at 0.00 and Peak / COM / FWHM read `N/A`, which looks like a
broken viewer and is not.

Pick **`tetramm1_sum_all`** as Y (or Z, in the mesh viewer). The picoammeter
reads real current whatever the attenuator is doing. An alignment `dscan`
defaults to `att_ratio=1e6`, so this is the normal case, not the exception.

## Other things worth knowing

* **Blank cells in a 2-D grid mean "never scanned", not zero.** A mesh stopped
  with Ctrl+C shows the part it measured, in the right cells, and the status
  line says how many of the grid it got to.
* **The 2-D viewer defaults X / Y / Z to the first, second and last column** —
  the scanned motors and your main counter, because that is the order the scans
  write. Any pick of your own wins over that.
* **The 1-D viewer does not pick a Y for you.** X defaults to the first column,
  but Y is a multi-select list and starts empty, so a freshly opened file plots
  nothing until you choose one. That is not a broken file — click a Y (see the
  attenuation note above for which one).
* **The file list shows one prefix at a time.** Ours is `A`; the launcher passes
  it. Started without it the viewer looks for `S*` and an 8-ID folder lists as
  empty.
* **A live scan is safe to open.** The `.csv` is closed after every point, and
  the viewer tolerates catching it mid-write.

## Where the launcher lives

`~/bin/start_scanviewer.sh` is a **thin wrapper**. The real script is version
controlled in this repo:

```
scripts/start_scanviewer.sh
```

**Edit that one.** A change to the copy in `~/bin` is invisible to everyone
else and is lost on a re-clone. Set `BLUESKY_DIR` if your checkout is not at
`~/bluesky`.

The viewers themselves are a separate checkout, `~/Documents/BLUETELLA_9ID`, on
branch `8id_test`. They are 9-ID's code — see below before changing them.

## When to contact Peco Myint

**BLUETELLA is Peco Myint's (9-ID).** We use it; we do not own it. The 8-ID
changes on `8id_test` are deliberately written to be beamline-neutral so they
can go back upstream rather than fork — the string "8-ID" appears nowhere in
either viewer.

**Contact Peco when:**

| situation | why him |
|---|---|
| the viewer crashes, mis-plots, or a statistic looks wrong | it is his code, and the fix belongs upstream where 9-ID gets it too |
| you want a new viewer feature (a fit, an export, a second Y axis) | so it lands in `main` and survives our next pull |
| you are about to change `scanviewer.py` or `meshviewer.py` | agree the approach first; an 8-ID-only patch is a fork we then maintain for ever |
| a pull request from us is open and needs a decision | he is the reviewer and approver |
| you want to know what 9-ID's own file format does | he defined it |

**Do NOT contact him for** — these are ours, and he cannot fix them:

| situation | ours, in this repo |
|---|---|
| a column is missing, wrong, or you want another PV recorded | [`configs/scan_csv_template.yml`](../src/id8_common/configs/scan_csv_template.yml) |
| the `.csv` is malformed, or the header is wrong | [`plans/align/scan_csv.py`](../src/id8_common/plans/align/scan_csv.py) |
| the scan itself misbehaves — motion, shutter, Ctrl+C, detector | [`plans/align/ophyd_scan.py`](../src/id8_common/plans/align/ophyd_scan.py), and [Ophyd scans](../src/id8_common/plans/align/OPHYD_SCAN.md) |
| the launcher opens the wrong folder | `scripts/start_scanviewer.sh` and `configs/experiment.yml` |
| the file list is empty | almost always the prefix or the folder — see above |

**Before you write to him**, it is worth being able to say which side the
problem is on. The quickest test: open the `.csv` in a text editor. If the file
looks wrong, it is ours. If the file looks right and the viewer shows something
else, it is his.

Second useful test — the viewers run with no EPICS and no beamline, so you can
reproduce off the instrument and send him a file:

```bash
start_scanviewer.sh --dir /path/to/a/folder/with/one/csv
```

And run his tests, which need nothing but Python:

```bash
cd ~/Documents/BLUETELLA_9ID && python3 -m unittest discover -s tests -q
```

## See also

* [Ophyd scans and the CSV file template](../src/id8_common/plans/align/OPHYD_SCAN.md)
  — the scans, their arguments, and what the `.csv` contains
* `~/Documents/BLUETELLA_9ID/USING_SCANVIEWER.md` — the operator's guide to the
  window itself: panels, overlays, statistics
* `~/Documents/BLUETELLA_9ID/8ID_SCANVIEWER.md` — where BLUETELLA is installed,
  what we changed, and what should go upstream
