# Proposal: automatic artifact masking in pySimpleMask

**For:** Miaoqi Chu · **From:** Q. Zhang · **Date:** 2026-09-17

Diamond-anvil-cell patterns at 8-ID-E carry narrow dark lines crossing the
scattering, plus localised bright spots — neither detected automatically today.
Example frame:

```
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/E0157_EHEA-Mesh_a0001_f000010_eiger4M_r00006/E0157_EHEA-Mesh_a0001_f000010_eiger4M_r00006.h5
```

Below: a way to find them, and what a prototype did on real data.

---

## 1. What they are

**Kossel lines from the diamond anvil.** Diffuse scattering from the sample acts
as a divergent source illuminating the downstream anvil, a large single crystal.
Directions meeting a Bragg condition get diffracted out of the transmitted cone,
leaving a narrow dark line where intensity is *missing*.

Two properties matter for the software:

- **They are dark, not bright** — subtractive. Every line found in this study was
  dark.
- **They are fixed by the cell orientation.** They move drastically when the cell
  is rotated and barely at all when the sample is translated (established from
  operator experience, not from a dedicated measurement).

The bright spots are the complementary case: an anvil or coarse-grain reflection
landing *on* the detector.

The practical consequence: an artifact mask is valid for one cell orientation,
and stays valid across a translation mesh.

---

## 2. Proposed algorithm

The assumption is that true I(q, φ) is smooth for a powder-like sample, so
anything abrupt is an artifact. Build the mask from a **high-statistics sum**, not
a single frame — the Eiger example averages 0.287 counts/pixel, and the lines are
visible only because the eye integrates along them. Artifacts are static within an
orientation, so summing a mesh is free and valid.

Apply the blemish layer before thresholding, eroded a few pixels so gap and border
neighbours go with it — pySimpleMask already loads one by default. The prototype
skipped this and immediately tripped: five of its six strongest "line" detections
were the detector's last column.

```
1. reference    I_ref(q) = median over φ within fine q bins
2. residual     R = I − I_ref(q)
3. integrate    R_s = smooth(R)              # lines are EXTENDED
4. threshold    σ = 1.4826 · MAD(R_s);  flag |R_s| > n·σ
5. classify     connected components → aspect ratio: high = line, low = spot
```

Use the **median**, not the mean, at step 1 — it does not chase the outliers being
hunted.

A note on using derivatives, which was the original suggestion: the instinct is
right, but at ~0.3 counts/pixel a literal finite difference is pure shot noise.
Steps 2–3 are the same idea implemented as a matched filter — difference from a
smooth reference at the feature's own scale — which is noise-optimal.

### Three things that break it

All three were hit while prototyping; they are the useful part of this proposal.

**(a) The reference must follow the true q contours.** A first attempt used
"median along detector rows". On Lambda2M that is nearly right, because the
virtual beam centre sits ~19,000 rows off-detector and rings are almost straight
horizontal lines. **On Eiger4M it fails completely** — the q gradient runs across
columns there, so a row-median averages straight across the signal and detects
nothing. Use the real q map.

So on Eiger, where the geometry is uncertain: the q map has to be *locally*
correct, not absolutely correct — a wrong beam centre distorts ring shape and
smears the reference. Where geometry is hopeless, see (c).

**(b) Dead pixels and physical artifacts look identical to a threshold.** They are
not the same thing and have different lifetimes:

```
depth / I_ref ≈ −100%        dead pixel      → permanent, belongs in blemish
depth / I_ref ≈ −10…−60%     Kossel line     → per orientation
```

This is the guard against a *stale* blemish map: dead pixels accumulate between
updates, and whatever is not yet flagged will be the detector's first find — which
is exactly what happened (§3).

**(c) Geometry-free fallback.** Where the q map cannot be trusted, a line is still
a line: detect narrow extended features directly in pixel space with a
Radon/Hough transform of the residual, or morphological opening with a line
structuring element at several orientations. Less sensitive, but needs no
geometry. Worth shipping alongside the q-φ path.

---

## 3. What the prototype found

Run over a sample of the 1398 Lambda2M datasets from pope202609, plus the Eiger
example.

**Kossel lines were not confirmed on Lambda2M.** After excluding detector edges
and dead pixels, no convincing physical line survived in the 16 datasets sampled.
Possible reasons, not distinguished: Lambda2M subtends only ~2.2° per frame at
2.2 m so intercepts few anvil cones; the sampling was sparse; sensitivity scales
with local intensity. The algorithm is therefore validated so far on *detector*
defects rather than on the anvil lines it targets, and the Eiger frame remains the
reference case.

**What it did find was dead pixels** — a band at rows 112–125, columns 712–828 on
Lambda2M reading 0.0045 counts against a detector mean of 1.5833, with 1479 pixels
reading exactly zero while flagged good in the current blemish TIFF. Probed at
fixed coordinates across eight mesh points, depth/reference was −97% to −100%
every time.

That is useful here only as evidence, in two directions: the detector works, and
**the dead-pixel split in §2(b) is not optional** — without it, defects swamp the
physical artifacts you are trying to find. The Lambda2M bad-pixel map itself is a
beamline responsibility and will be established separately, from a flat-scattering
standard such as NIST glassy carbon measured in transmission geometry. Nothing is
being asked of pySimpleMask on that front.

### A model-selection step worth adding

On the Eiger frame the recorded geometry is **falsified by the data**. Comparing
residual RMS under three reference models on the 1500-frame sum:

| reference model | residual RMS |
|---|---|
| I(row) | 0.42934 |
| I(2θ) from metadata geometry | 0.42664 |
| **I(col)** | **0.05179** |

The 2θ map built from the recorded centre, distance and delta is no better than
assuming rings run along rows — while an empirical I(col) reference fits 8× better.
The Eiger sits on a manual stage here, so its geometry is known to be unreliable;
the point is that **this test detects that condition cheaply**, before the
reference model silently ruins the detection. Recommend trying a small set of
reference models and picking by residual RMS rather than trusting the metadata.

### Current state of the prototype on Eiger

With the empirical I(col) reference at 4σ it finds **2 dark lines** (589 px at
row 507/col 816, 370 px at row 1495/col 895) and masks 0.040% of the live area.
That is a sane mask but an incomplete result — roughly four lines are visible by
eye in the frame. Raising the threshold to 7σ loses both.

**No validated artifact-removed dataset exists yet.** The limitation is that the
residual still contains real structure, so the MAD threshold has no clean
separation between artifact and signal. A better reference model — a 2-D smooth
surface fit rather than a 1-D profile — is the obvious next thing to try.

---

## 4. Suggested interface

Three mask layers with different lifetimes, kept separate:

| layer | lifetime | source |
|---|---|---|
| blemish | permanent, per detector | site TIFF (maintained by the beamline) |
| **artifact** | per cell orientation | proposed auto-detection |
| user | per analysis | manual ROI |

A minimal CLI surface consistent with the existing `build` options:

```
--auto-artifact {off,lines,spots,both}   default off
--artifact-nsig N                        default ~6
--artifact-min-aspect N                  line vs spot, default ~5
--artifact-dead-frac F                   depth/ref beyond this = dead pixel,
                                         routed to the blemish layer (default 0.9)
--artifact-from FILE [FILE ...]          build from a high-statistics sum,
                                         apply to the target
--artifact-exclude-border N              default ~8 px
--output-artifact-mask FILE              write the layer separately
```

Two points of principle:

- **Report, never silently mask.** Print counts and total area removed, and write
  the layer as its own file. A mask that quietly deletes 5% of a detector is worse
  than the artifact when it is wrong.
- **Default off**, until it has a track record.

---

## 5. Suggested order

1. Implement the q-φ detector with border exclusion and the dead-pixel split.
   Validate against the Eiger frame, where the lines are visible by eye.
2. Add the geometry-free fallback for detectors with poor q calibration.
3. Spot detection last — compact bright features are easier, and
   `--threshold-high` already covers part of it.

---

## Appendix — reproduction

Prototype scripts, on `amber`:

| file | does |
|---|---|
| `/home/beams10/8IDIUSER/xpcs_contrast_check/artifact_detect.py` | core detector: reference, residual, classify |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/fish_lambda.py` | q-binned reference from metadata; scans many datasets |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/line_persist.py` | tracks a detected line across mesh points |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/line_fixed_probe.py` | probes fixed pixels across mesh points (the §3 test) |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/sum_list.py` | high-statistics sum over a scan |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/eiger_clean.py` | reference-model comparison (the table above) |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/eiger_mask_out.py` | detection + writes the products below |

Products, for inspection — **diagnostics, not validated output**:

| file | contents |
|---|---|
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_before_after.png` | 3 panels: before / residual / after. **Start here** |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_sum.npy` | 1500-frame sum, float64 (2162×2068) |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_sum_bad.npy` | detector sentinel/gap mask, bool |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_ref.npy` | fitted reference I(col) |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_artifact_mask.npy` | detected artifacts, bool (1681 px) |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_artifact_mask.tif` | same, 8-bit TIFF |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_cleaned.npy` | sum with artifacts + dead pixels set to NaN |
| `/home/beams10/8IDIUSER/xpcs_contrast_check/e0157_tth.npy` | 2θ map from metadata (the one the test rejects) |

Data referenced:

```
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/E0157_EHEA-Mesh_a0001_f000010_eiger4M_r00006/
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/G0242_HEA-Mesh30_a0001_f000010_lambda2M_r00227/
/home/beams/8IDIUSER/Documents/areaDetectorBlemish/8idLambda2m/latest_blemish.tif
```

All read-only with respect to the beamline: file reads and compute only.
