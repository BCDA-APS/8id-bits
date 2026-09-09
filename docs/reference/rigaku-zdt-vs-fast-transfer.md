# Rigaku 3M: ZDT 2-bit vs Fast Transfer 2-bit

Measured 2026-09-08 on `pearl`. Compares the two ways of getting 2-bit data off
the Rigaku 3M, across four photon-counting lower thresholds.

| | mode string | detector name | files written |
|---|---|---|---|
| **ZDT** | `ZDT2bit` | `rigaku3M` | six `.bin.000` … `.bin.005` |
| **FTF** | `ZDT2bit` | `rigaku3M_ftf` | six `.h5.000` … `.h5.005` |

Both write six files, one per detector module, plus one `_metadata.hdf`.

## Bottom line

**Compression:** neither format wins outright — they cross over. ZDT is up to
**4× smaller** than FTF on sparse data (7 keV) and **12.6× larger** on dense
data (4 keV), where it exceeded even uncompressed storage by 10.6×. ZDT costs
64 bits per photon, so it beats dense storage only below **3.1 % occupancy**;
FTF costs a whole 786 KB frame whenever that frame has any count at all, which
bounds it between ~58 MB and 7.86 GB per 10 000 frames.

**Acquisition time:** FTF is faster and flatter everywhere — **2.7–3.0 s per
repeat regardless of data volume**, against ZDT's 2.9–12.2 s. At 10 000 frames
this is a comparison of per-measurement overhead, not throughput.

**Recommendation:** keep ZDT as the default for real XPCS measurements, which
live in the sparse regime where it is 4× better and analysable. Reach for FTF
only if a measurement is genuinely dense — and note that FTF output currently
cannot be analysed at all (see [below](#what-you-can-do-with-the-data-afterwards)),
so today the choice is really "ZDT, and keep the occupancy below a few percent".

## ⚠ Read this before using the numbers

**1. 10 000 frames per repeat, not the standard 100 000.** Chosen to keep the
low-threshold points from filling the disk — ZDT at 4 keV writes 84 GB per
repeat even at 10 000 frames, so 100 000 would have been ~840 GB per repeat and
8.4 TB for the ten repeats of that one point. Sizes scale with frame count and
can be multiplied by ten with confidence. **Times cannot**: at 20 µs per frame,
10 000 frames is 0.2 s of exposure inside a 2.7–12 s repeat, so what the timing
column measures is almost entirely per-measurement overhead (arm, write, close,
metadata), not streaming throughput. See [Acquisition efficiency](#acquisition-efficiency).

**2. There was no beam.** The counts are detector dark noise, so the lower
threshold acts as a count-rate knob rather than as a real photon-energy cut.
What that noise actually consists of — and the finding that ~188 pixels produce
most of it — is in [Rigaku 3M dark noise](rigaku-3m-dark-noise.md).
That is exactly what a compression comparison needs — a controlled sweep from
dense to sparse — but the absolute rates are not what a real sample would give,
and the threshold at which the two formats cross over will move with the actual
count rate.

## Conditions

```
detector      rigaku3M, 2-bit
acq_time      0.00002 s (20 us)
num_frames    10 000
num_repeats   10 per (threshold, mode)
att_level     2
sample_move   no
thresholds    4, 5, 6, 7 keV   (8idRigaku3m:cam1:LowerThreshold)
```

80 measurements, 931 GB written, plus a 20-measurement repeat of the 6 keV
point (1.2 GB) to confirm the timing result below.

## Results

Mean over the 10 repeats of each point. `sd` is the standard deviation of the
per-repeat size — reproducibility was excellent throughout (worst case 0.2 %).

| threshold | mode | MB / repeat | sd | s / repeat | sd | total |
|---|---|---|---|---|---|---|
| 4 keV | ZDT | 83 583.36 | 24.16 | 12.18 | 0.12 | 835.83 GB |
| 4 keV | FTF | 6 652.64 | 2.65 | 3.02 | 0.07 | 66.53 GB |
| 5 keV | ZDT | 1 807.50 | 3.96 | 2.85 | 0.06 | 18.08 GB |
| 5 keV | FTF | 899.10 | 1.54 | 2.66 | 0.04 | 8.99 GB |
| 6 keV | ZDT | 39.02 | 0.02 | 6.47 | 0.08 | 0.39 GB |
| 6 keV | FTF | 77.96 | 0.05 | 2.80 | 0.10 | 0.78 GB |
| 7 keV | ZDT | 14.48 | 0.01 | 6.53 | 0.08 | 0.14 GB |
| 7 keV | FTF | 58.27 | 0.01 | 2.76 | 0.09 | 0.58 GB |

### Which files are which threshold

All under `/gdata/dm/8ID/8IDE/2026-3/comm202609/data/`. Every measurement is ten
folders, `_r00001` … `_r00010`, each holding six data files plus one
`_metadata.hdf`. The measurement numbers are **not** in threshold order —
6 and 7 keV were measured first, then 4 and 5 keV were re-measured after the
first attempt was discarded.

| threshold | mode | folder | data files |
|---|---|---|---|
| 4 keV | ZDT | `A0167_Thr4keV_a0002_f010000_r000NN` | `.bin.000` … `.bin.005` |
| 4 keV | FTF | `A0168_Thr4keV_a0002_f010000_r000NN` | `.h5.000` … `.h5.005` |
| 5 keV | ZDT | `A0169_Thr5keV_a0002_f010000_r000NN` | `.bin.000` … `.bin.005` |
| 5 keV | FTF | `A0170_Thr5keV_a0002_f010000_r000NN` | `.h5.000` … `.h5.005` |
| 6 keV | ZDT | `A0163_Thr6keV_a0002_f010000_r000NN` | `.bin.000` … `.bin.005` |
| 6 keV | FTF | `A0164_Thr6keV_a0002_f010000_r000NN` | `.h5.000` … `.h5.005` |
| 7 keV | ZDT | `A0165_Thr7keV_a0002_f010000_r000NN` | `.bin.000` … `.bin.005` |
| 7 keV | FTF | `A0166_Thr7keV_a0002_f010000_r000NN` | `.h5.000` … `.h5.005` |
| 6 keV | ZDT — repeat | `A0171_Thr6keV_a0002_f010000_r000NN` | `.bin.000` … `.bin.005` |
| 6 keV | FTF — repeat | `A0172_Thr6keV_a0002_f010000_r000NN` | `.h5.000` … `.h5.005` |

For these ten, **the `ThrNkeV` in the folder name is the real threshold** — each
was recorded by a single process whose threshold readback was sampled every
0.5 s and never varied.

**`A0153`–`A0162` are the discarded first attempt and their names lie about the
threshold** — see [the appendix](#appendix-the-discarded-first-attempt). There
is no correct threshold to assign them; they are not in the table above and
should be deleted rather than re-labelled.

## Compression efficiency

The natural reference is the **dense** size — every pixel stored, no
compression. The FTF file's internal layout gives it directly:

```
entry/data/data    shape (131072, 10000)  uint8  chunks (131072, 1)  no compression
```

131 072 bytes per frame per module × 6 modules = **786 KB per frame**, so a
10 000-frame measurement is **7.86 GB** dense. Against that:

| threshold | ZDT | FTF |
|---|---|---|
| 4 keV | **10.6× larger than dense** | 0.85× dense |
| 5 keV | 4.4× smaller | 8.7× smaller |
| 6 keV | 202× smaller | 101× smaller |
| 7 keV | 543× smaller | 135× smaller |

**The two formats fail in opposite directions, and they cross over between 5
and 6 keV.**

*ZDT is an event list.* Its size is proportional to the number of photons, with
no upper bound. When the detector is sparse it is superb — 543× at 7 keV, four
times better than FTF. When the detector is flooded it is catastrophic: at
4 keV it wrote **ten times more than storing every pixel uncompressed**, because
each event costs more bits than the 2 bits a dense pixel would have taken.

### The break-even rule

Every ZDT file measured, at every threshold, has a size divisible by 8, so the
record is **8 bytes = 64 bits per event**. A dense 2-bit pixel costs 2 bits.
So ZDT breaks even against dense storage at an occupancy of

```
2 bits / 64 bits = 3.125 % of pixels hit per frame
```

Below that ZDT wins, above it ZDT is worse than storing everything. The 3M has
524 288 pixels per module × 6 = **3 145 728 pixels**, so:

| threshold | events / frame | occupancy | predicted vs dense | measured |
|---|---|---|---|---|
| 4 keV | 1 044 792 | 33.2 % | 10.6× larger | 10.6× |
| 5 keV | 22 594 | 0.72 % | 4.4× smaller | 4.4× |
| 6 keV | 488 | 0.0155 % | 202× smaller | 202× |
| 7 keV | 181 | 0.0058 % | 543× smaller | 543× |

This makes ZDT's behaviour predictable rather than empirical: **estimate the
occupancy and you have the file size.** Anything above ~3 % occupancy should not
be recorded in ZDT.

*FTF is chunked HDF5, one chunk per frame, chunks allocated only when written.*
That bounds it at both ends: it can never exceed the 7.86 GB dense size (at
4 keV it reached 0.85× of it, i.e. ~85 % of frames had at least one count), and
it has a floor of ~58 MB set by per-frame chunk overhead, which is why it is
only 135× at 7 keV where ZDT manages 543×.

The practical rule, independent of threshold:

> **Sparse data → ZDT. Dense data → FTF.** ZDT wins by up to 4× when the
> detector is quiet and loses by more than 12× when it is not. FTF's worst case
> is bounded; ZDT's is not.

For a real XPCS measurement — a speckle pattern, not a noise flood — the sparse
regime is the normal one, which is where ZDT is the better choice. **But there
is no protection against the dense case:** a mis-set threshold, a mis-set
attenuator, or a direct-beam hit will silently produce a file an order of
magnitude larger than dense storage, at a rate of ~7 GB/s. Nothing in the
acquisition path checks for this.

## Acquisition efficiency

**FTF is flat.** 2.66–3.02 s per repeat across a 7.4× range in file size — the
per-repeat cost is essentially independent of how much data it writes, over this
range.

**ZDT is not, and it is not monotonic either**: 12.18 s at 4 keV, 2.85 s at
5 keV, 6.47 s at 6 keV, 6.53 s at 7 keV. Only the 4 keV point is explained by
I/O (84 GB in 12.18 s ≈ 7 GB/s). The 6 and 7 keV points take more than twice as
long as 5 keV while writing 46× less data, which is backwards for anything I/O
bound.

That non-monotonicity is **real, not an artefact**. The 6 keV point was measured
twice, by two separate processes about an hour apart:

| | size | time |
|---|---|---|
| 6 keV ZDT, run 1 (`A0163`) | 39.02 MB | 6.48 ± 0.08 s |
| 6 keV ZDT, run 2 (`A0171`) | 39.05 MB | 6.49 ± 0.07 s |
| 6 keV FTF, run 1 (`A0164`) | 77.96 MB | 2.80 ± 0.10 s |
| 6 keV FTF, run 2 (`A0172`) | 77.94 MB | 2.78 ± 0.12 s |

Two mechanisms were checked and ruled out: the acquire loop polls at 0.1 s, too
fine to quantise a 3.6 s difference; and all six modules write a non-empty file
at every threshold, so it is not one module stalling. **The cause is not
known.** It costs ~3.6 s per repeat in the sparse regime — the regime real XPCS
measurements run in — so it is worth chasing: over a 1 000-repeat series it is
an hour.

Because 10 000 frames at 20 µs is only 0.2 s of exposure, 93–98 % of every
repeat here is fixed overhead. **This table therefore ranks per-measurement
overhead, not throughput**, and the ranking would compress at 100 000 frames
where the exposure is 2 s and the files are ten times larger.

## What you can do with the data afterwards

These are two separate problems, and it is worth not conflating them:

| | can `boost_corr` read the format? | does the DM job succeed? |
|---|---|---|
| ZDT `.bin.000` | **yes** | **no** — but for an unrelated, Polaris-side reason; the same file analyses fine locally |
| FTF `.h5.000` | **no** | n/a |

Tested directly on this data (`A0166` rep 1), rather than assumed:

```
TypeError: File type [.000] is not supported
  -> boost_corr.xpcs_aps_8idi.exceptions.DatasetError
```

The dispatch in `boost_corr/xpcs_aps_8idi/dataset/utils.py` has an explicit
`raw_fname.endswith(".bin.000")` branch and a `.tpx.000` branch, but for
`.h5.000` `os.path.splitext` yields `.000`, which is in none of the lists, so it
falls through to the `raise`. The file itself is perfectly good HDF5 — it
carries the `\x89HDF` magic number.

**This is more than a missing extension.** The `.bin.000` branch dispatches to
`Rigaku3MDataset`, which knows to find the other five module files and unpack
2-bit frames. Supporting FTF needs an equivalent reader for its layout (six
files, `entry/data/data` chunked one frame at a time), not just another
`elif`.

So FTF's compression advantage in the dense regime is currently unusable end to
end — the data is written but nothing will read it — and that is a code gap in
`boost_corr`, not a property of the format. ZDT's failure is a deployment
problem that does not affect the data. Both are tracked in
[Data Management](data-management.md) and in the README's Known issues.

## Reproducing this

Everything needed is in [`scripts/benchmarks/`](../../scripts/benchmarks). On a
beamline host, from `src/`, in the `8id_bits` environment:

```bash
python ../scripts/benchmarks/rigaku_zdt_vs_ftf.py 4 5 6 7   # ~8 min, ~930 GB
python ../scripts/benchmarks/analyse.py                     # rebuilds the table above
```

`analyse.py` with no arguments re-reads the published measurement numbers; pass
your own in `ZDT FTF` pairs to analyse a fresh sweep.

Two details in the driver are not incidental, and both exist because of how the
first attempt failed (see below):

1. **A threshold watchdog.** A background thread samples the threshold readback
   every 0.5 s and reports every distinct value seen during the block. A clean
   block prints a single value (`### RBV seen during block: {4.0: 331}`).
   Anything else means something is writing the PV underneath you, and that
   block's data is worthless — the file name will say one threshold and the
   contents will be another.
2. **Re-raise the interrupt.** `det_acq_series()` converts Ctrl+C into a
   `RuntimeError`, so a driver that catches broad exceptions will swallow the
   interrupt and march on to the next threshold.

## Appendix: the discarded first attempt

`A0153`–`A0162` are on disk, are labelled with thresholds, and are **wrong**.
Do not use them.

A first benchmark process did not die when it was killed; it kept running and
stepping its own threshold sweep while the replacement run started. Two
processes then drove the same detector and wrote the same
`LowerThreshold` PV. The result is data whose file name says one threshold and
whose contents are another: `A0157` is named `Thr4keV`, but repeat 7 is 83.6 GB
(the 5 keV size) and repeat 10 is 39 MB (the 6 keV size). 27 of those
directories contain only a `_metadata.hdf` and no data at all, because the other
process owned the detector at the time.

Two lessons, both cheap:

* **Directory mtimes are the audit trail.** The overlap was invisible in the
  log and obvious the moment the directories were sorted by time — `A0156_r00001`
  and `A0157_r00001` share a timestamp to the second.
* **A measurement number does not identify a run.** `expt.measurement_num` comes
  from a single EPICS register (`8ideSoft:Reg1`), so two concurrent processes
  interleave and take alternate numbers. `A0157`, `A0159`, `A0161` belonged to
  one process and `A0158`, `A0160` to the other.
