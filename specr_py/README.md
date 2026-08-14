# specr_py — live SPEC file viewer

A Python/Qt rewrite of the MATLAB `specr` GUI (Zhang Jiang, APS), for watching
8-ID scans as they are written.

Pairs with `dscan_ophyd()` in
[ophyd_scan.py](../src/id8_common/plans/align/ophyd_scan.py), which appends one
complete SPEC row per scan point and closes the file each time. Because the file
is append-only and never held open by the writer, the viewer can poll it while a
scan runs without any risk of collision — this is the whole reason the scan
writes a text file rather than having the viewer read the detector `.h5`.

## Where things live

| | path | visible from |
|---|---|---|
| SPEC files | `~/bluesky/<experiment>.spec` | amber **and** kouga (NFS home) |
| Detector images | `/gdata/dm/8ID/8IDE/<cycle>/<exp>/data/bluesky/*.h5` | amber only (GPFS) |

`/gdata` is **not** mounted on kouga, which is why the SPEC file goes to
`~/bluesky` instead of next to the images. Scans run on amber; the viewer can run
anywhere that sees the home directory.

## Running it

```bash
conda activate 8id_bits

python ~/bluesky/specr_py/specr.py                              # then File > Open
python ~/bluesky/specr_py/specr.py ~/bluesky/pope202607.spec    # open a file
python ~/bluesky/specr_py/specr.py ~/bluesky/pope202607.spec --watch   # open + monitor
```

Needs a display: run it on a workstation, or `ssh -X`.

## What it does

- **Scan Monitor** (`Ctrl+M`, or View ▸ Scan Monitor) polls the file and redraws
  as points land. Period is configurable in Tools ▸ Settings (default 2 s), with
  the MATLAB auto-period behaviour — if a refresh takes longer than the interval,
  the interval backs off so polls never queue up.
- **Erase mode** (View ▸) — on, the monitor shows only the newest scan; off, it
  accumulates scans that share the same columns.
- **Select Scan** — multi-select list; scans must share a column layout to be
  overlaid. `Ctrl+←` / `Ctrl+→` step through scans.
- **X / Y** — any column. Defaults match the MATLAB reader: X is the first
  column, Y is the last. `ophyd_scan.py` deliberately writes the scanned motor
  first and `stats1` last (SPEC's "primary detector goes last" convention), so
  the default plot is the one you want.
- **Peak / COM / FWHM** — shown under the X axis, ported from `params.m`
  including its fallbacks when a half-maximum crossing is not bracketed.
- Plot styles (linear / logx / logy / logxy), grid, legend, axis inversion,
  Show Current Scan, Show Motor Positions (from `#O`/`#P`), Save Displayed Data.

## Layout

- `specfile.py` — incremental SPEC parser. No GUI dependencies, so it is usable
  on its own for analysis and is what the tests exercise.
- `specr.py` — the Qt front end.
- `tests/` — offline tests plus `reference_example.spe`, a fixture covering the
  awkward parts of real Bluesky-written SPEC files: repeated `#F` blocks, labels
  containing single spaces (`SMS granite`), bare `True` in a numeric column,
  `#N 0` scans with no columns, and `abort`/`fail`/`success` exit states.

```bash
cd ~/bluesky/specr_py
QT_QPA_PLATFORM=offscreen python tests/test_specr.py
```

Covers parsing the reference file, round-tripping what `ophyd_scan.py` writes,
incremental live reads (including that a torn trailing row is ignored until its
newline arrives), the statistics port, and a headless GUI live-update run.
