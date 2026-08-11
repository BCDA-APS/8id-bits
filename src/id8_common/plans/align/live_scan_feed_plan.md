# Local live-scan feed for 8-ID scan plans

Status: proposed, not implemented. Written 2026-08-10.

## Context

`dscan` and its siblings in [scan_8id.py](src/id8_common/plans/align/scan_8id.py) hand-assemble a
Bluesky event per scan point (`bps.create` / `bps.read` × 4 / `bps.save`). The goal is a live
visualization tool for alignment scans, built without adding a second copy of the scan record and
without risking concurrent access to the AreaDetector `.h5` file.

Three options were considered and rejected before this one:

- **A. One `.h5` per scan point.** Large time and space overhead, and it does not help — the pixels
  are not what the alignment curve plots.
- **B. A bespoke scan text file.** Duplicates a record that already exists in two places (the Tiled
  catalog and the SPEC `.dat` callback, both already enabled).
- **C. Write a databroker-shaped export file and feed it to gemviz.** Does not compose: gemviz has no
  file reader. Its entire data layer (`gemviz/tapi.py`) is `from_uri(uri, "dask")` — it is a Tiled
  HTTP client.

### Findings that shape this plan

- **Reading the `.h5` live is the wrong channel and is unnecessary.** A reader opened `'r'` cannot
  corrupt an HDF5 file — corruption needs a second *writer*. The real hazards are (a) POSIX advisory
  lock contention, which can make the **IOC's** write fail, and (b) a stale metadata view without
  SWMR. Both vanish if the live path never opens the file.
- **The `.h5` does not contain the plotted quantity.** [scan_structure.txt](src/id8_common/plans/align/scan_structure.txt)
  lists `Huber_*` positions but no `Stats*_Total` in NDAttributes. Reproducing `lambda2M_stats1_total`
  from the file means re-reducing every frame (41 × 1813 × 1558 for a typical scan). The per-step
  event is therefore *not* redundant with the file for the exact quantity being plotted.
- **The existing Tiled path is not viable for live viz.** `iconfig.yml` sets `TILED_PROFILE_NAME: 8ide`,
  so `apsbits.init_RE` wraps the catalog in `TiledWriter(cat, batch_size=1)` and subscribes it to the
  RunEngine — one synchronous HTTP POST to `albite.xray.aps.anl.gov:8000` per point, inside the scan
  loop. Fast alignment feedback should not sit behind that.
- **gemviz cannot read the current server.** gemviz pins tiled `0.1.0a116`; the server runs `0.2.3`.
  Three independent breakages: it expects `run[stream]["data"]` (the server serves a Composite with
  parts `['internal']`), it expects `stream_metadata()["descriptors"]` (absent from the new schema),
  and `plottable_signals()` requires `start["hints"]`, which the hand-built `md` dicts omit.
- **The plan already owns the correlation.** Independent CA monitors cannot reliably pair a motor
  position with the detector total for the same point — the motor has moved on by the time the stats
  callback fires, and two equal consecutive totals are indistinguishable from no new point. The
  per-step block can, and does.

**Outcome:** the per-step block additionally emits each `(position, signals)` pair over a
fire-and-forget local datagram. The viz tool is a separate process that listens. No file access, no
Tiled, no HTTP, no broker, no new persisted record.

## Design decision: UDP datagram to localhost

`socket.sendto()` on a connectionless UDP socket cannot block and does not raise when nothing is
listening. It therefore **cannot slow down or break a scan** — strictly safer than the synchronous
HTTP write already in the loop. It is cross-process and cross-machine (unlike an in-process
matplotlib buffer), and needs only stdlib `socket` + `json`.

Rejected alternatives:

- **EPICS waveform PVs** would fit the `pv_registers` convention
  ([registers_device.py](src/id8_common/devices/registers_device.py)) and CLAUDE.md's note that it is
  "shared with non-Bluesky tools" — but that device is all scalar `Reg`/`StrReg`. It would need new
  records in the `8ideSoft` IOC, outside this repo.
- **`bluesky.callbacks.zmq`** (`Publisher` + `RemoteDispatcher`) would reintroduce a broker process to
  run and babysit, and a Bluesky dependency in the viz tool.
- **Pure Ophyd CA monitors in the viz process** — no plan changes at all, but hits the correlation
  problem above. Would need `ArrayCounter_RBV` / `NDArrayUniqueId` as a sequence key and would still
  be racy at the edges.

## Recommended: keep `bps.create/read/save`

This is the one point where this plan argues against the original instinct. It is ~5 floats per point,
and it is the only queryable record of what an alignment scan did. Since the stats totals are absent
from the `.h5`, it is not redundant with the file. The new feed is *additional*, not a replacement.

If that turns out to be unconvincing in practice, deleting it later is one line per site — the shared
helper below makes that a single edit rather than sixteen.

## Changes

### 1. New module: `src/id8_common/utils/live_feed.py`

A `LiveFeed` class, ~40 lines, no dependencies beyond stdlib:

- `__init__(host="127.0.0.1", port=DEFAULT_PORT)` — creates one `SOCK_DGRAM` socket, reused for the run.
- `start(**md)` / `point(**kw)` / `stop(**kw)` — each builds
  `{"type": ..., "seq": n, "ts": time.time(), ...}`, `json.dumps`, `sendto`. Auto-incrementing `seq`
  so the receiver can detect drops. UDP does not guarantee delivery; for a localhost loopback datagram
  under ~64 KB this is reliable in practice, but the receiver should tolerate gaps.
- Every method wrapped in a bare `except Exception: pass`. **The feed must never be able to raise into
  the RunEngine.**
- Module-level `DEFAULT_PORT` so the viz tool imports the same constant rather than hardcoding it.

### 2. One shared per-step helper in `scan_8id.py`

The `bps.create` → `bps.read(motor)` → `bps.read(det.stats1/2/3)` → `bps.save()` block appears
**16 times** (verified by `grep -c 'bps.create("primary")'`) across `dscan` / `dmesh` / `mesh` /
`d2scan` / `ascan` / `a2scan` / `dscan_auto` × the lambda and eiger branches. Representative sites:
lines 223–228, 288–293, 439–445, 501–506, 1088–1093, 1954–1959.

Collapse to one module-level helper:

```python
def _read_point(motors, det, feed=None, extra=None):
    """Emit one primary event and mirror it to the live feed."""
    yield from bps.create("primary")
    for m in motors:
        yield from bps.read(m)
    for n in (1, 2, 3):
        yield from bps.read(getattr(det, f"stats{n}"))
    yield from bps.save()
    if feed is not None:
        feed.point(**{m.name: m.position for m in motors},
                   **{f"{det.name}_stats{n}_total": getattr(det, f"stats{n}").total.get()
                      for n in (1, 2, 3)},
                   **(extra or {}))
```

Adding the feed once instead of sixteen times is the whole reason to do this refactor now.

`bps.read()` returns the reading, so `feed.point()` should reuse those values rather than issuing a
second `.get()` per signal — that avoids 4 extra CA round-trips per point. Confirm the reading-dict
shape when wiring it up.

### 3. Instantiate the feed per run

In each plan, alongside the existing `folder_prefix` / `md` construction (e.g. lines 774–777):

```python
feed = LiveFeed()
feed.start(plan_name="dscan", motors=[motor.name], detector=det.name,
           num_points=num_pts, image_file=folder_prefix,
           signals=[f"{det.name}_stats{n}_total" for n in (1, 2, 3)])
```

and `feed.stop(reason=...)` in the existing `finally:` block. The `start` payload gives the viz tool
everything it needs to label axes, pre-size arrays, and reset between scans.

### 4. Viz tool (separate process, separate file)

Outside this repo's scope, but the receiving end is ~25 lines:
`socket.bind(("127.0.0.1", DEFAULT_PORT))` → on `start` clear the axes and read the labels → on `point`
append and redraw → on `stop` finalize. `matplotlib` + `FuncAnimation`, or PyQtGraph for a snappier
redraw.

To run it on another machine, change the destination host — UDP handles that unchanged.

## Fixes worth folding in

Small, independent, same file or its neighbours:

- **Duplicate `dscan`.** [scan_8id.py:145](src/id8_common/plans/align/scan_8id.py#L145) and
  [scan_8id.py:753](src/id8_common/plans/align/scan_8id.py#L753) are byte-identical definitions. Line
  753 shadows line 145, so edits to the first have no effect. Delete lines 145–328.
- **`check_h5_shape` leaks a file handle.** [check_file_dim.py](src/id8_common/utils/check_file_dim.py)
  does `h5py.File(path)` and never closes it. h5py 3.15 defaults to mode `'r'` so it will not corrupt
  anything, but the handle — and its read lock — stays open for the life of the Bluesky session, and it
  is imported into the session namespace by [startup.py:115](src/id8_common/startup.py#L115). This is a
  live instance of exactly the hazard that motivated this whole design. Use a `with` block.
- **Optional, cheap insurance: add `hints` to `md`.**
  `"hints": {"dimensions": [([motor.name], "primary")]}` — one entry per `md` dict. Standard Bluesky
  plans emit this; the hand-built dicts do not, which is why no generic client can pick the axes.
  Costs nothing and keeps the Tiled/gemviz door open without depending on it.
- **Stale AD path templates.** `write_path_template: /gdata/dm/8IDI/2025-2/` in
  [ad_devices.yml](src/id8_common/configs/ad_devices.yml) does not match the real path built in
  `save_images` (`/gdata/dm/8ID/8IDE/2026-2/<exp>/data/bluesky`). Harmless today because no
  Resource/Datum documents are emitted, but it will produce wrong paths the moment they are.
- **`scan_8id_dev.py`.** Contains 10 more copies of the same block and is imported nowhere. Decide
  whether to update it alongside or delete it; do not leave it silently diverging. Same reasoning as
  CLAUDE.md's note about `id8_common_dev`.

## Verification

1. **Feed in isolation, no beamline.** Run a receiver on `127.0.0.1:<DEFAULT_PORT>`; drive `LiveFeed`
   from a plain Python loop. Confirm ordering, `seq` continuity, and that `start`/`point`/`stop`
   round-trip.
2. **Feed cannot break a scan.** With *nothing* listening, run a `dscan`. It must complete normally.
   Then kill the receiver mid-scan and confirm the same.
3. **Simulated scan.** `RE(dscan(motor, -0.1, 0.1, 11, 0.1, det=..., save_img=0))` against a sim
   detector, receiver running. Curve should build point by point.
4. **Correlation check — the point of the whole exercise.** On a real `dscan` over a feature with
   structure, compare the live curve against `run["primary"]["internal"]["events"].read()` pulled from
   Tiled afterwards. The arrays must match element for element. If they do not, the pairing is wrong.
5. **Refactor is behaviour-preserving.** For one scan of each type, confirm the Tiled event count and
   column names are unchanged from a pre-refactor run of the same plan.
6. **No file contention.** Confirm the viz tool holds no handle on any `.h5`:
   `lsof -p <pid> | grep '\.h5'` must be empty during a scan.

There is no unit-test suite in this repo — the `test_*.py` files under `plans/` and `user/` are
beamline scripts. Verification is by running the plans.

## Not doing

- Reading the `.h5` from a second process, with or without SWMR. Nothing in the live-curve requirement
  needs it, and SWMR over GPFS would have to be validated separately.
- Emitting Resource/Datum, or porting gemviz to tiled 0.2.3. Both are real gaps; neither is needed for
  live curves. Worth revisiting if image display is ever wanted.
- Removing the existing `bps.create/read/save`. See the reasoning above; the helper makes it a
  one-line change later if wanted.