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
experiment's** scan folder, assembled from three settings in
`configs/experiment.yml` exactly the way the scans assemble it:

```
<mount_point><cycle_name>/<experiment_name>/data/bluesky
```

so it follows a cycle or experiment rollover with no edit.

**If the launcher says that folder does not exist, check `mount_point` first.**
A stale `mount_point` still naming the other station's tree — `8IDI` where this
experiment is `8IDE`, or the reverse — is the usual cause. The launcher prints
all three values so you can see which is wrong, then opens anyway on an empty
browser; fix `experiment.yml`, or pass `--dir` for now.

Run it where `/gdata` is mounted — pearl — and over `ssh -Y`, since it opens a
window. The launcher says so itself if either is missing.

```bash
ssh -Y 8idiuser@pearl
start_scanviewer.sh
```

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
| `--dir <folder>` | some other folder; suppresses the automatic lookup. A viewer started by hand with no `--dir` opens on the current directory instead — **Browse...** moves it |
| `--scan-prefix <letter>` | the launcher already passes `A`; override only if your files start with something else |
| `--meta-prefix <letter>` | 1-D viewer only; 9-ID's metadata files, unused at 8-ID |
| `--help` | the viewer's own options |

`--mesh` is the only flag the launcher consumes; everything else is handed
straight through. Tick **Follow latest** to track the running scan. Both viewers
are read-only: no EPICS, no ophyd, no Data Management. They only ever read
files.

## The first thing that confuses people

**Under high attenuation the `lambda2M_stats*` columns are all zero.** The plot
is then a flat line at 0.00 and Peak / COM / FWHM read `N/A`, which looks like a
broken viewer and is not.

Pick **`tetramm1_sum_all`** as Y (or Z, in the mesh viewer). The picoammeter
reads real current whatever the attenuator is doing. An alignment `dscan`
defaults to `att_ratio=1e6`, so this is the normal case, not the exception.

## The window

```
 Directory: [ /gdata/dm/8ID/8IDE/2026-3/comm202609/data/bluesky ]
            [Browse...] [Refresh]     (o) Scan  ( ) Meta   [x] Follow latest
 --------------------------------+------------------------------------------
  Files (multi-select: overlay)  |
    File    Number Modified Rows |
    A0114_...  114  14:21      5 |               the plot
    A0115_...  115  14:28      5 |
    A0116_...  116  14:34      5 |  [ matplotlib toolbar: pan / zoom / PNG ]
                                 +------------------------------------------
  X column:  [ combo         v ] |  Peak / COM / FWHM (raw linear data)
                                 |   File  X  Y  Peak X  Peak Y  COM  FWHM
  Y column(s) (multi-select)     +------------------------------------------
    elapsed_time                 |  Status:
    lambda2M_stats1_total        |   A0116...: missing column, curve skipped
    tetramm1_sum_all             |   (a running log -- scroll back through it)
                                 |
  [ ] Log X   [ ] Log Y          |
  [x] Show Peak/COM/FWHM overlay |
  [ ] Show derivative            |
```

Left is what you choose, right is what you get. The **Status** box at the
bottom right is a running log, not a single line, so an earlier message is never
overwritten. When something does not plot, the reason is almost always sitting
in it.

**The file list shows one prefix at a time.** It lists `A*.csv` and `A*.json` —
ours is `A`, and the launcher passes it. Started without it the viewer looks for
`S*` and an 8-ID folder lists as empty. The **Number** column is parsed from the
`A<digits>_` part of the name; a file that does not match shows `?` there and
sorts to the top, and **Follow latest** will never choose it.

**Rows is blank until the file has been read.** It fills in the moment you
select the file, and keeps up while a scan is running.

## Plotting a scan

1. **Pick a file** in the list on the left.
2. **Pick the X column** from the dropdown. It defaults to the first column of
   the file, which for a scan CSV is the moved motor.
3. **Pick one or more Y columns** from the list below it.

**The 1-D viewer does not pick a Y for you.** The Y list is a multi-select and
starts empty, so a freshly opened file plots nothing until you choose one. That
is not a broken file — click a Y (see the attenuation note above for which one).

The column choices come from whichever files are loaded, so they follow the
detector: a `lambda2M` scan offers `lambda2M_stats1_total` … `stats4`, a
`tetramm` scan offers its current channels. A pick you have already made is kept
when you move to another file that still has that column, so stepping through a
run of scans does not mean re-choosing every time.

**Several Y columns stack, they do not share an axis.** Each Y gets its own
panel, one above the other, sharing X. That is deliberate — a counter reading
1e-9 A and one reading 1e6 counts on the same axis would flatten one of them.

### Overlaying several scans

Ctrl-click or shift-click several **files** and they are drawn on the same axes,
one colour and one legend entry each, with the same X/Y choice. This is the
usual way to compare an alignment before and after a change.

A file that does not have the chosen X or Y column is skipped rather than
guessed at, with a line in **Status** naming it. Nothing is inferred from column
names: if two scans call the same counter two different things, they will not
overlay.

The mesh viewer is **single-select** instead — there is no meaningful way to
overlay two colour maps, so it plots one file at a time.

## Following a scan as it runs

Tick **Follow latest**. Once a second the viewer re-lists the directory,
switches to the **highest-numbered** file as soon as one appears, and re-reads
only that one file as it grows. Leave it open, start a `dscan_ophyd()` from your
Bluesky session (or plain `dscan()` from the Ophyd-only one), and watch the plot
fill in point by point. Nothing needs restarting between scans.

**A live scan is safe to open.** The `.csv` is closed after every point, and the
viewer tolerates catching it mid-write: a new file that is still header-only is
skipped until it reads cleanly, so a half-written scan cannot blank a good plot.

Two things to know:

* It picks the current scan **from file names**, not from EPICS — the highest
  measurement number wins, with modification time only as a tiebreak. Copy an
  old file with a higher number into the folder and it will follow that instead.
* Following works in **Scan** mode only, and only inside the one directory you
  opened. After a cycle or experiment rollover, restart it —
  `start_scanviewer.sh` picks up the new folder on its own.

## Reading the numbers

The table under the plot gives one row per drawn curve: three files by two Y
columns is six rows.

| column | meaning | how it is computed |
|---|---|---|
| **Peak X / Peak Y** | position and height of the maximum | plain `argmax` — the largest **measured point**, not a fit. On a noisy scan that is the noisiest point |
| **COM** | centre of mass | trapezoidal `∫x·y dx / ∫y dx` over the whole scan, so a background offset drags it toward the middle of the range |
| **FWHM center** | midpoint of the half-maximum crossings | linear interpolation between the points that bracket half the peak |
| **FWHM** | full width at half maximum | right crossing − left crossing |

**These are always computed on the raw linear data.** Log X and Log Y change
how the curve is drawn and nothing else; the numbers do not move.

**⚠ A half-maximum crossing that does not exist is replaced by the end of the
scan.** If the curve never comes back down below half its peak on one side — a
monotonic scan, or a peak the scan did not fully bracket — that side's crossing
becomes the last point of the scan, `FWHM center` collapses onto the crossing
that *was* found, and the reported FWHM is really "peak to the edge of the
scan": a lower bound, not a width. If neither side crosses, FWHM is the whole
scan range and the centre is its midpoint. Widen the scan rather than trust the
number.

A statistic that cannot be computed meaningfully is reported as `N/A` rather
than as a nonsense value, and the row is all-or-nothing: fewer than two points,
or any single value coming out NaN or infinite, and the whole row reads `N/A`.
An all-zero counter does exactly this — `COM` is `0/0` — which is why the
attenuation trap above shows five `N/A`s. During a live scan the row appears
once there are two points.

**Show Peak/COM/FWHM overlay** draws the same three things on the plot itself:
a star at the peak, the centre of mass, and the half-maximum span.

## The other controls

* **Log X / Log Y** — logarithmic axes. Non-positive values cannot be shown on a
  log axis, so if any are present the viewer says which column in **Status** and
  falls back to linear **for that plot only**, rather than silently dropping the
  points. This is the second way a heavily attenuated scan announces itself.
* **Show derivative** — adds a `dy/dx` panel under each curve. It is computed in
  **acquisition order**, not sorted by X, because for a beamline scan the order
  the points were taken in is part of the data. It is refused, with a note in
  Status, when X has duplicate or non-finite values — dividing by a zero step
  would manufacture a number rather than measure one.
* **Scan / Meta** — which family of files to browse: **Scan** lists the scan
  prefix (`A` here), **Meta** the meta prefix (`M`). One at a time, never both.
  Meta is 9-ID's metadata snapshots; at 8-ID you will only ever use Scan.
  Switching to Meta also jumps to a sibling `meta/` directory if one exists —
  ours does not, so it simply lists nothing. Following is disabled in Meta mode.
* **Refresh** — re-list the directory now instead of waiting for the next poll.
* **Browse...** — change directory.
* **matplotlib toolbar** — pan, zoom, and **save the current plot as a PNG**.
* **Status** — read it first. Skipped columns, log-scale fallbacks, refused
  derivatives, transient read errors and partial rasters are all reported here.

## Mesh scans

`meshviewer.py` is a separate window for two-motor rasters. Same layout, but
X, Y and Z are three dropdowns, the file list is single-select, and there is no
log, derivative or statistics table — just one extra checkbox, **Draw raster as
grid**. **Follow latest** and the **Status** log behave exactly as in the 1-D
viewer. Start it through the same launcher:

```bash
start_scanviewer.sh --mesh
start_scanviewer.sh --mesh --dir /gdata/dm/8ID/8IDE/2026-3/comm202609/data/bluesky
```

Run `meshviewer.py` by hand only off the beamline, in an environment you have
already set up.

### Picking X, Y and Z

They default to the **first**, **second** and **last** column of the file — for
a scan CSV that is the two moved axes and the main counter, which is usually
what you want. Nothing reads a column's name or contents to guess a meaning; the
defaults are purely positional. Your own pick always wins, and is kept when you
move to another file that has the same column.

If Z is flat and single-coloured, check the attenuation before you check the
beam.

### A raster is drawn as an image

A `dmesh` or `mesh` scan is drawn as a filled image, one coloured cell per grid
point, instead of scattered dots. The plot title carries the grid, for example
`A0118_Test_a1010041   (5x5 grid)`.

What decides this is the **file**, not the scan's name:

| you ran | header declares a shape | drawn as |
|---|---|---|
| `dmesh` / `mesh` | yes | filled image |
| `d2scan` / `a2scan` | no | scatter |

`d2scan` and `a2scan` write no shape on purpose. They sweep both motors along a
single line, so there is no 2-D field to colour — a grid would be almost all
empty and would imply an area nobody measured.

Two further conditions, both of which fall back to the scatter and say so in
Status:

* **X and Y must be the two motors the scan declared.** Put `elapsed_time` on an
  axis and the declared grid no longer describes what is being plotted, so it is
  not used. The status line names the two motors it expected.
* **Untick Draw raster as grid** to force the scatter back at any time.

Cells are matched to the nearest **commanded** position, never by reshaping the
rows in order. A stage only repeats a position to within its retry deadband, so
a 5x5 raster really does contain 25 slightly different readbacks, and reshaping
would put a point in the wrong cell the first time one was retried.

**A blank cell means never scanned, not zero.** A mesh stopped with Ctrl+C shows
the part it measured, in the right cells, and Status says how many of the grid's
cells are filled.

## If something looks wrong

**The file list is empty.** Almost always the prefix or the folder. The viewer
lists only names starting with the active prefix — `A0046_Test_a1010041.csv`
under `--scan-prefix A`, which is what `start_scanviewer.sh` passes. Started
without it the prefix is `S` and an 8-ID folder looks empty. Also check the
directory, and that you are in Scan mode rather than Meta.

**Nothing is drawn.** You need both an X column and at least one Y column. The Y
list is a multi-select and starts empty.

**"Follow latest" is greyed out.** You are in Meta mode. Switch to Scan.

**The plot is a flat line at zero.** Check the attenuation first — see the
section above. If `tetramm1_sum_all` is zero too, then the counter really is
zero: check the shutter and that there is beam. The viewer plots what the file
contains.

**A statistic reads `N/A`.** Fewer than two points so far, or the value came out
NaN or infinite (an all-zero counter does this). This is deliberate: `N/A` is
honest, a fabricated number is not.

**An FWHM looks far too wide.** The scan probably did not come back down below
half maximum — see the warning under "Reading the numbers".

**A curve is missing from an overlay.** That file does not have the chosen X or
Y column. Status names it.

**Check Status.** Every one of the above that the viewer can detect is reported
there.

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
* [Scan viewer reference](reference/scan-viewer.md) — where BLUETELLA is
  installed, the conda environment, what we changed and why, and what is open
  upstream
