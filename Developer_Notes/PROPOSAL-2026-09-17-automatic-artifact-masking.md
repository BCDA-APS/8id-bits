# Proposal: automatic artifact masking in pySimpleMask

**For:** Miaoqi Chu · **From:** Q. Zhang · **Date:** 2026-09-17

Diamond-anvil-cell patterns at 8-ID-E carry narrow dark streaks crossing the
scattering — neither detected nor masked automatically today. Reference frame:

```
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/G0256_HEA-Dec9GPa_a0001_f003000_lambda2M_r00001/
```

Lambda2M, where the q map is trustworthy. Below: a method that finds them, and
what it recovers on that frame.

---

## 1. What they are

**Kossel lines from the diamond anvil.** Diffuse scattering from the sample acts
as a divergent source illuminating the downstream anvil, a large single crystal.
Directions meeting a Bragg condition get diffracted out of the transmitted cone,
leaving a narrow dark line where intensity is *missing*.

Two properties matter for the software:

- **They are dark, not bright** — subtractive. Every line found in this study was
  dark.
- **They are fixed in the lab frame, by the cell orientation**, not in detector
  pixels. They move drastically when the cell is rotated and barely at all when
  the sample is translated (operator experience). The confirmed case in §3 backs
  this from the data: the same pixels show a −23% line at one detector angle and
  nothing 4° away.

The bright spots are the complementary case: an anvil or coarse-grain reflection
landing *on* the detector.

The practical consequence: an artifact mask is valid for one cell orientation,
and stays valid across a translation mesh.

---

## 2. Proposed algorithm

Work in the **q-φ frame from pysimplemask itself**, never in a hand-rolled one.
`compute_transmission_qmap()` returns per-pixel `q`, `TTH` and `phi`; the swing
angles come from `flightpath_swing` (horizontal) and `flightpath_swing_vertical`
(vertical), and the beam centre is the direct-beam pixel at **zero swing on both**.
For the reference frame here that gives q ∈ 3.111–3.455 Å⁻¹ and φ ∈ 87.19–92.82°.

Two facts make the rest easy, and both come straight out of those maps:

```
dq/drow = -1.88e-4 A-1/px     dq/dcol = 0
dphi/dcol = -0.00344 deg/px   dphi/drow = 0
```

**q depends only on row, φ only on column.** A narrow q bin is therefore a band of
rows, and I-vs-φ within it is simply the intensity profile across columns.

### The method

```
0. restrict to the peak    only pixels inside the analysis ROI matter, and the
                           high intensity there makes a fractional dip far easier
                           to see than in the wings
1. narrow q bins           ~4 rows (~0.00075 A-1)
2. I versus phi per bin    powder scattering is flat in phi; a Kossel line is a
                           localised DIP at one phi
3. find the discontinuity  running median over phi as baseline, flag phi where
                           I falls > 3.5 robust sigma below it
4. LINK across q bins      a Kossel line is continuous, so its dip drifts smoothly
                           in phi. Track bin-to-bin, predicting the next phi by
                           extrapolating the slope so far
5. rasterise               interpolate each track and paint its measured width
```

**Step 4 is what makes it work.** Kossel lines are conic sections, so they are
curved and run at arbitrary angles. Tracking follows that curvature; anything that
assumes a fixed orientation does not.

### Three things that break it

**(a) Do not collapse a profile over the whole detector.** A column profile
averaged over all rows smears a drifting dip into nothing and returns a confident
null. The dip must be found *per narrow q bin*, then linked.

**(b) Do not use a straight-line matched filter.** A prototype scanned straight
kernels at 10° steps. Because the lines are curved a straight kernel only matches
short chords, so it fragmented long lines and missed shallow ones entirely: 7
fragments where there are 5 continuous streaks, and two never found at all.

**(c) Dead pixels and Kossel lines look identical to a threshold.** A dead pixel
reads *exactly zero*; a Kossel line is a fractional dip. Split on depth ratio —
this also guards against a stale blemish map, since unflagged dead pixels are
otherwise the detector's first find.

## 3. Result on real data

`G0256_HEA-Dec9GPa_a0001_f003000_lambda2M_r00001`, δ = 24.48°, 300-frame average,
restricted to the analysis ROI q = 3.2632–3.3126 Å⁻¹ (rows 740–1019).

178 raw dips across 69 q bins, linked into 11 tracks, merged into **5 distinct
streaks**:

| # | q bins | rows | cols | angle | depth |
|---|---|---|---|---|---|
| 1 | 6 | 914–970 | 220–247 | +115.7° | −27.8% |
| 2 | 5 | 806–854 | 278–304 | +118.4° | −27.2% |
| 3 | 7 | 782–850 | 544–566 | +72.1° | −27.8% |
| 4 | 66 | 742–1002 | 761–774 | +87.6° | −35.7% |
| 5 | 58 | 802–1014 | 1481–1555 | +109.2° | −40.8% |

**3949 px, 1.0% of the live pixels in the ROI.**

![Kossel lines on Lambda2M, G0256](figures/kossel_phi_final.png)

*Full detector; dashed lines mark the peak ROI. Left: original. Middle: residual
against I_ref(q) — the streaks are invisible in the raw frame and unmistakable
here. Right: masked pixels in magenta. Dark bands are module gaps and the five
dead chip-gap columns, which are detector structure, not artifacts.*

Five angles spanning 72° to 118° confirms these are not detector structure — no
detector artifact is slanted, and none is slanted five different ways.

The streaks overlap the sample peak at q = 3.288 Å⁻¹ used for the XPCS analysis of
this same dataset, contributing a few percent of its measured static baseline.
That is the concrete case for the feature.

## 4. Suggested interface

Three mask layers with different lifetimes, kept separate:

| layer | lifetime | source |
|---|---|---|
| blemish | permanent, per detector | site TIFF (maintained by the beamline) |
| **artifact** | per cell orientation | proposed auto-detection |
| user | per analysis | manual ROI |

A minimal CLI surface consistent with the existing `build` options:

```
--auto-artifact {off,streaks}     default off
--artifact-qrange QLO:QHI         restrict to the ROI that matters; the high
                                  intensity inside a peak is what makes a
                                  fractional dip detectable
--artifact-qbin N                 rows per narrow q bin (default ~4)
--artifact-nsig N                 dip threshold in robust sigma (default ~3.5)
--artifact-min-bins N             q bins a track must span to be kept (default ~5)
--artifact-dead-frac F            depth/ref beyond this = dead pixel, routed to
                                  the blemish layer (default 0.9)
--artifact-from FILE [FILE ...]   build from a high-statistics sum, apply to target
--output-artifact-mask FILE       write the layer separately
```

Two points of principle:

- **Report, never silently mask.** Print counts and total area removed, and write
  the layer as its own file. A mask that quietly deletes 5% of a detector is worse
  than the artifact when it is wrong.
- **Default off**, until it has a track record.

---

## 5. Suggested order

1. Add the per-pixel q/φ accessor if one is not already public —
   `compute_transmission_qmap()` already returns everything needed.
2. Implement narrow-q-bin dip detection with φ linking. **Validate against G0256**
   (§3): known answer, trustworthy geometry, five streaks from −27% to −41%.
3. Expose `--artifact-nsig` and a q-range restriction, and report detections
   rather than applying them silently.
4. Spot detection last — compact bright features are easier, and
   `--threshold-high` already covers part of it.

---

## Appendix — reproduction

Prototype scripts, on `amber`, all under
`/home/beams10/8IDIUSER/xpcs_contrast_check/`:

| file | does |
|---|---|
| `phi_track2.py` | narrow-q-bin dip detection + φ linking (the detector) |
| `phi_final.py` | merges tracks, rasterises the mask, writes the figure |
| `hunt_kossel.py` | earlier full-detector search (superseded) |
| `sum_list.py` | high-statistics sum over a scan |
| `psm_maps.npy` | per-pixel q / TTH / phi from pysimplemask |

Products, same directory:

| file | contents |
|---|---|
| `kossel_phi_final.png` | the figure above |
| `k2_img.npy` | 300-frame average of G0256 (1813×1558) |
| `phi_mask_final.npy` | the five detected streaks, bool (3949 px) |
| `phi_dips.npy` | raw per-bin dips before linking |
| `k2_good2.npy` | valid-pixel mask, blemish applied |

Data referenced:

```
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/G0256_HEA-Dec9GPa_a0001_f003000_lambda2M_r00001/
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/G025{3,4,5}_HEA-Dec9GPa_*_f003000_lambda2M_r00001/
/home/beams/8IDIUSER/Documents/areaDetectorBlemish/8idLambda2m/latest_blemish.tif
```

All read-only with respect to the beamline: file reads and compute only.
