# Live scan visualization: measured cost of one-file-per-point (Plan A)

Measurement report. Written 2026-08-10. All numbers measured on `amber` against real
8-ID-E data in `/gdata/dm/8ID/8IDE/2026-2/`, not estimated.

Companion document: [live_scan_feed_plan.md](src/id8_common/plans/align/live_scan_feed_plan.md)
(the UDP-feed proposal this report cross-evaluates against).

---

## 1. The question

The alignment plans in [scan_8id.py](src/id8_common/plans/align/scan_8id.py) emit a Bluesky event
per scan point (`bps.create` / `bps.read` × 4 / `bps.save`, e.g. lines 288–293). A live
visualization tool is wanted. Reading the AreaDetector `.h5` from a second process while the IOC
writes it is the obvious hazard, so three options were on the table:

- **Plan A** — write each point to its own `.h5`. Each file is closed before the next opens, so a
  reader can safely open file *N* while the IOC writes *N+1*. A separate script polls for new files,
  computes ROI intensity from pixels, and plots it against motor position (which the file already
  carries as an NDAttribute).
- **Plan B** — a bespoke scan text file. Rejected: duplicates records that already exist.
- **UDP feed** — the plan emits each `(position, signals)` pair as a fire-and-forget datagram.

The concern about Plan A was "lots of overhead both in time and file size." This report measures
that overhead.

**Headline result: the overhead is far smaller than feared. Time and space are noise. The only
large number is file count.**

---

## 2. What the files actually are

14 most recent files in `/gdata/dm/8ID/8IDE/2026-2/pope202607/data/bluesky`:

```
file                               frames          HxW    dtype          chunk    rawMB     imgMB    fileMB   ovhdMB
A0156_G10_a0001.h5                      3    1813x1558   uint32 (1, 1813, 1558)     33.9      0.13      0.45    0.311
A0157_G10_a0001.h5                     81    1813x1558   uint32 (1, 1813, 1558)    915.2     84.80     85.12    0.318
A0158_G10_a0001.h5                      8    1813x1558   uint32 (1, 1813, 1558)     90.4      2.27      2.58    0.311
A0159_G10_a0001.h5                     81    1813x1558   uint32 (1, 1813, 1558)    915.2    116.87    117.19    0.318
A0160_G10_a0001.h5                     81    1813x1558   uint32 (1, 1813, 1558)    915.2     66.06     66.38    0.318
A0162_G10_a0001.h5                     81    1813x1558   uint32 (1, 1813, 1558)    915.2     65.81     66.13    0.318
A0163_G10_a0001.h5                     11    1813x1558   uint32 (1, 1813, 1558)    124.3      2.61      2.92    0.311
A0164_G10_a0001.h5                     81    1813x1558   uint32 (1, 1813, 1558)    915.2    121.06    121.37    0.318
A0165_G10_a17274.h5                    41    1813x1558   uint32 (1, 1813, 1558)    463.2      2.14      2.33    0.191
A0166_G10_a17274.h5                    41    1813x1558   uint32 (1, 1813, 1558)    463.2      2.15      2.34    0.191
A0167_G10_a17274.h5                    41    1813x1558   uint32 (1, 1813, 1558)    463.2      2.14      2.34    0.191
A0169_HEA-15GPa_a0001.h5               61    1813x1558   uint32 (1, 1813, 1558)    689.2     74.01     74.26    0.254
A0170_HEA-15GPa_a0001.h5               20    1813x1558   uint32 (1, 1813, 1558)    226.0      5.75      6.00    0.248
A0171_HEA-15GPa_a0001.h5               61    1813x1558   uint32 (1, 1813, 1558)    689.2     71.01     71.26    0.254
```

Key structural facts:

- `uint32`, 1813 × 1558 → **11.3 MB raw per frame**.
- **Chunk = `(1, 1813, 1558)` — exactly one frame.**
- Compression is HDF5 filter **32004 (LZ4)**, applied per chunk. `h5py` reports
  `compression: None` because it is a third-party filter; the 10–20× ratio between `rawMB` and
  `imgMB` is real.
- 83 NDAttribute datasets per file, each chunked at the frame count, 328 B stored for a 41-frame
  scan (22.0 KB total).

**The single most important consequence:** because the chunk is exactly one frame and LZ4 runs
per chunk, splitting a scan into one file per point **cannot change the compressed image bytes at
all**. Every byte of Plan A's cost is structural overhead plus per-file write calls. This is what
makes the cost bounded and measurable rather than a guess.

Full byte breakdown for `A0165_G10_a17274.h5` (41 frames):

| component | bytes |
|---|---:|
| `/entry/data/data` (LZ4) | 2,142,894 |
| 83 NDAttribute datasets | 22,550 |
| `performance/timestamp` | 1,640 |
| HDF5 structure (headers, B-trees, superblock) | 166,444 |
| **total file** | **2,333,528** |

---

## 3. Benchmark

Method: read the first 41 real frames from a source file, then write them two ways with an
identical AD-style layout (LZ4 image dataset chunked one frame per chunk, 83 chunked NDAttribute
datasets, performance table) — once as a single 41-frame file, once as 41 single-frame files.
Compare bytes and wall time. Then read each back and compute a 400 × 400 ROI sum per frame, both
with a handle held open and with open/close per frame.

Two regimes were tested because the overhead ratio depends entirely on how compressible the frames
are:

- **sparse** — `A0165_G10_a17274.h5`, heavily attenuated, ~52 KB/frame compressed. Representative
  of alignment scans.
- **dense** — `A0169_HEA-15GPa_a0001.h5`, ~1.2 MB/frame compressed. Representative of data
  collection.

Writes went to local `/tmp` on amber; file **sizes** are filesystem-independent, so those are exact.
Write **times** on GPFS will be somewhat higher, so the timing figures are a lower bound. Reads were
measured both on `/tmp` and against the real GPFS files.

Environment: amber, Python 3.9, h5py 3.14.0, hdf5plugin 5.0.0, 32 cores, `/gdata` = GPFS
(25.2 PB, 31% used).

### Results

```
sparse (alignment-like)  <- A0165_G10_a17274.h5
source: 41 frames (1813, 1558) uint32 | img 2.143 MB | file 2.334 MB | 83 attrs

  WRITE (local /tmp, lower bound for GPFS)
    1 file  x 41 frames :    2.182 MB     0.18 s   (4.4 ms/frame)
    41 files x  1 frame :    3.766 MB     0.40 s   (9.7 ms/frame)
    delta               :   +1.583 MB (+72.5%)   +0.22 s (+122.3%)
    per-file structural overhead: 37.7 KB

  READ + ROI sum, 41 frames
    open/close per frame   : 0.415 s  (10.1 ms/frame)
    one handle held open   : 0.394 s  ( 9.6 ms/frame)
    open/close penalty     : +0.5 ms/frame
    real file on GPFS      : 0.395 s  ( 9.6 ms/frame)

dense (data-like)  <- A0169_HEA-15GPa_a0001.h5
source: 61 frames (1813, 1558) uint32 | img 74.008 MB | file 74.262 MB | 83 attrs

  WRITE (local /tmp, lower bound for GPFS)
    1 file  x 41 frames :   49.218 MB     0.38 s   ( 9.4 ms/frame)
    41 files x  1 frame :   50.802 MB     0.60 s   (14.6 ms/frame)
    delta               :   +1.583 MB (+3.2%)   +0.22 s (+56.1%)
    per-file structural overhead: 37.7 KB

  READ + ROI sum, 41 frames
    open/close per frame   : 0.564 s  (13.8 ms/frame)
    one handle held open   : 0.541 s  (13.2 ms/frame)
    open/close penalty     : +0.6 ms/frame
    real file on GPFS      : 0.542 s  (13.2 ms/frame)
```

### Reading the results

- **Per-file structural overhead is 37.7 KB and is constant** across both regimes, as predicted
  from the per-chunk compression argument. The percentage figures differ (+72.5% vs +3.2%) only
  because the denominator differs.
- **The open/close penalty is 0.5–0.6 ms/frame.** Opening a fresh HDF5 file per point is
  effectively free. This was the cost most likely to sink Plan A, and it does not.
- **GPFS reads match warm local reads exactly** (9.6 and 13.2 ms/frame). The read is bound by LZ4
  decompression and the ROI sum, not by I/O. Note that a chunk is a whole frame, so requesting only
  an ROI still decompresses the full 11.3 MB frame — there is no way to make this cheaper by
  slicing.

---

## 4. At real 8-ID-E scale

Scan-size distribution from the Tiled catalog (`8ide` profile, `/raw`), last 600 runs:

| metric | value |
|---|---|
| plan mix | 573 `dscan`, 10 `count`, 9 `mesh`, 3 `rel_scan`, 3 `a2scan`, 2 `ascan` |
| points per run | min 1, **median 40**, mean 67, p90 100, **max 10,000** |
| total points | 40,050 |

Current file counts in cycle 2026-2: 6,001 `.h5` total; 1,066 across the eight
`*/data/bluesky` directories (55–422 per experiment).

Extrapolating Plan A across those 600 runs:

| | now | Plan A | delta |
|---|---:|---:|---|
| wall-clock write time | — | — | **+212 s (3.5 min) total**, ≈0.3% of scan time |
| disk | — | — | **+1.51 GB** |
| files | 600 | 40,050 | **67×** |

Context for each:

- **Time.** 3.5 minutes spread over 600 scans. A single `dscan` runs ~2 s/point (measured: 41 points
  in 81 s), dominated by motor motion. The extra ~5.3 ms/point is invisible.
- **Disk.** 1.51 GB against 17.5 PB free. Invisible.
- **Files.** This is the only figure that is genuinely large.

Worst case worth naming: one run in that sample had **10,000 points**. Under Plan A that is 10,000
files and 377 MB of pure structural overhead for a single scan.

---

## 5. Cross-evaluation

| | Plan A (file per point) | UDP feed |
|---|---|---|
| Live latency | ~13 ms decompress + file-ready detection | ~0.1 ms |
| Extra disk | +37.7 KB/point | 0 |
| Extra scan time | ~+5.3 ms/point | ~+0.05 ms/point |
| Concurrency risk | eliminated by construction | eliminated by construction |
| Depends on Tiled / doc stream | no | no (feed itself is plain sockets) |
| Gives pixels, not just a curve | **yes** | no (would need the PVA plugin) |
| ROI changeable after the fact | **yes** | no — baked in at scan time by `stats1/2/3` |
| Replayable if viz wasn't running | **yes** | no (ephemeral) |
| Plan changes needed | HDF plugin mode in `save_images` | ~4 lines + per-step refactor |
| Downstream impact | **67× file count** | none |

### Correction to the earlier argument against Plan A

The main technical objection raised in [live_scan_feed_plan.md](src/id8_common/plans/align/live_scan_feed_plan.md)
was that the stats totals are absent from the `.h5`
([scan_structure.txt](src/id8_common/plans/align/scan_structure.txt) has `Huber_*` positions but no
`Stats*_Total`), so the plotted curve is not reproducible from the file alone.

**That objection does not survive under Plan A.** If ROI intensity is computed from pixels, the
stored stats totals are not needed, and the motor positions are already NDAttributes. Under Plan A
the `.h5` set is a complete, self-sufficient record and the Bluesky event adds nothing that cannot
be recomputed.

The two rows where Plan A is strictly better — pixels available, and ROI changeable after the fact —
are real scientific advantages during alignment, not just conveniences.

---

## 6. What these numbers do not cover

Two risks that the benchmark cannot see. Both must be checked before committing to Plan A.

### 6.1 Detector re-arm — the one number that could change the conclusion

The benchmark timed `h5py` writing files, not the IOC. If one-file-per-point were implemented by
putting the **cam** in Single mode, the Eiger's disarm/re-arm cycle would likely dwarf everything
measured above and could double or triple alignment scan time.

The way to avoid it: leave the cam pre-armed exactly as `dscan` does today
(`trigger_mode="Internal Series"`, `manual_trigger="Enable"`, `num_triggers=num_pts`) and change
only the **HDF plugin** to `FileWriteMode=Single` with auto-increment. `NDFileHDF5` writes one file
per NDArray callback independently of cam arm state, so the existing pre-armed software-trigger
scheme in [scan_8id.py:753](src/id8_common/plans/align/scan_8id.py#L753) should carry over unchanged.

**Verify on hardware before committing.** This is the only measurement that could invalidate the
conclusion.

### 6.2 Knowing a file is complete

A file appearing in a directory listing does not mean the IOC has closed it. The polling script
needs a real completion signal. Options, in rough order of robustness:

1. Gate on `hdf1:NumCaptured_RBV` / `WriteFile_RBV` over Channel Access.
2. Use the ordering guarantee — file *N* is definitely closed once file *N+1* exists. Costs one
   point of lag.
3. `inotify` plus a settle delay. Simplest, least reliable.

This is the genuinely fragile part of Plan A. It is solvable, but it is real work and it puts EPICS
back in the reader's loop.

---

## 7. Conclusion

**The overhead concern about Plan A was unfounded, and it should not be the deciding factor.** Time
cost is ~0.3% of scan time; disk cost is 1.5 GB per 600 scans on a 17 PB-free filesystem; the
open/close penalty that seemed most threatening is 0.5 ms per frame.

The real trade is not overhead. It is that **every downstream consumer currently assumes one file per
scan** — DM workflow upload, the XPCS analysis pipeline,
[check_file_dim.py](src/id8_common/utils/check_file_dim.py), qmap generation. 67× the files does not
strain GPFS, but it touches all of that.

Two scoping measures make Plan A comfortable:

- **Apply it only on alignment scans** (gated on `save_img` in `save_images`), not data collection.
  Those are the scans actually watched live, and they are the sparse-frame case; this caps the
  file-count growth to where it buys something.
- **Guard the mesh case.** At 10,000 points the file count and the 377 MB overhead both stop being
  negligible. Consider a point-count threshold above which the plan falls back to a single file.

If those hold and §6.1 checks out on hardware, Plan A is a sound design and is preferable to the UDP
feed whenever live pixel-level ROI — with an ROI that can be changed mid-alignment — is what is
actually wanted. The UDP feed remains the better choice if only the stats-plugin curve is needed,
since it costs nothing and adds no files.

Unrelated to the choice, and worth fixing either way:
[check_file_dim.py](src/id8_common/utils/check_file_dim.py) does `h5py.File(path)` and never closes
it, leaking a read lock for the life of the Bluesky session. Under Plan A a reader script will be
opening files constantly, which makes handle hygiene matter more, not less.

---

## Appendix: reproducing the benchmark

Three scripts were used; all are read-only against `/gdata` and write only to `/tmp` on amber.

**A. File characterization** — shape, dtype, chunking, storage size, per-file overhead:

```python
import h5py, os, glob
d = "/gdata/dm/8ID/8IDE/2026-2/pope202607/data/bluesky"
for f in sorted(glob.glob(d + "/*.h5"), key=os.path.getmtime)[-14:]:
    with h5py.File(f, "r") as h:
        ds = h["/entry/data/data"]
        img = ds.id.get_storage_size() / 1e6
        tot = os.path.getsize(f) / 1e6
        print(os.path.basename(f), ds.shape, ds.dtype, ds.chunks,
              "img=%.2fMB file=%.2fMB ovhd=%.3fMB" % (img, tot, tot - img),
              "nattr=%d" % len(h["/entry/instrument/NDAttributes"]))
```

**B. Filter identification** — confirms LZ4 (32004) and per-chunk compression:

```python
plist = ds.id.get_create_plist()
for i in range(plist.get_nfilters()):
    print(plist.get_filter(i))
```

**C. Write/read benchmark** — the script that produced §3. Reads N real frames, writes them as one
file and as N files with an equivalent AD layout, then times reads with and without per-frame
open/close. Run with:

```bash
ssh 10.54.116.10 python3 - < h5bench.py
```

Key parameters to keep if re-running: chunk `(1, H, W)`, `hdf5plugin.LZ4()`, 83 NDAttribute datasets
chunked at the frame count, and a `performance/timestamp` dataset of shape `(n, 5)`. Dropping any of
these will understate the per-file overhead.
