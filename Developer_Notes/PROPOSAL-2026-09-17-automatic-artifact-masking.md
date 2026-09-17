# Proposal: automatic artifact masking in pySimpleMask

**For:** Miaoqi Chu · **From:** Q. Zhang · **Date:** 2026-09-17

Diamond-anvil-cell patterns at 8-ID-E carry narrow dark lines crossing the
scattering, plus localised bright spots — neither detected automatically today.
Example frame: `E0157_EHEA-Mesh_10001_f00010_eiger4M_r00006`. Below: a way to find
them, what a prototype did on real data, and one bug worth fixing regardless.

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
anything abrupt is an artifact.

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

### Four things that break it

All four were hit while prototyping; they are the useful part of this proposal.

**(a) The reference must follow the true q contours.** A first attempt used
"median along detector rows". On Lambda2M that is nearly right, because the
virtual beam centre sits ~19,000 rows off-detector and rings are almost straight
horizontal lines. **On Eiger4M it fails completely** — the q gradient runs across
columns there, so a row-median averages straight across the signal and detects
nothing. Use the real q map.

So on Eiger, where the geometry is uncertain: the q map has to be *locally*
correct, not absolutely correct — a wrong beam centre distorts ring shape and
smears the reference. Where geometry is hopeless, see (e).

**(b) Exclude detector borders and module-gap edges first.** Without this they
dominate: in a six-dataset trial, five of the six strongest "line" detections were
at column 1556–1557 of a 1558-wide detector. Erode the valid-pixel mask a few
pixels before thresholding.

**(c) Single frames lack the statistics.** The Eiger example averages **0.287
counts/pixel**, maximum 9 counts anywhere; the lines are visible only because the
eye integrates along them. Build the mask from a high-statistics sum — artifacts
are static within an orientation, so summing a whole mesh is free and valid.

**(d) Dead pixels and physical artifacts look identical to a threshold.** They are
not the same thing and have different lifetimes:

```
depth / I_ref ≈ −100%        dead pixel      → permanent, belongs in blemish
depth / I_ref ≈ −10…−60%     Kossel line     → per orientation
```

Without this split the detector's first find is unflagged dead pixels — which is
exactly what happened (§3).

**(e) Geometry-free fallback.** Where the q map cannot be trusted, a line is still
a line: detect narrow extended features directly in pixel space with a
Radon/Hough transform of the residual, or morphological opening with a line
structuring element at several orientations. Less sensitive, but needs no
geometry. Worth shipping alongside the q-φ path.

---

## 3. What the prototype found

Run over a sample of the 1398 Lambda2M datasets from pope202609, plus the Eiger
example.

### A gap in the Lambda2M bad-pixel map

The strongest, most repeatable detection was a band at **rows 112–125, columns
712–828** on Lambda2M:

```
raw mean in the band            0.0045 counts   (detector mean 1.5833)
pixels reading exactly zero     1564 of 1638  (95.5%)
blemish map says GOOD (b == 1)  1553
read zero AND flagged good      1479 pixels
```

Probed at fixed pixel coordinates across eight mesh points, depth/reference was
**−97% to −100% every time** — dead pixels, not a Kossel line. Those 1479 pixels
are absent from

```
/home/beams/8IDIUSER/Documents/areaDetectorBlemish/8idLambda2m/latest_blemish.tif
```

and therefore enter every qmap built on this detector today. There are also
all-zero row bands at **902–907** and **1552–1557**, beyond the two known module
gaps (516–646, 1163–1296).

**Worth fixing on its own**, independently of whether the masking feature is built.

### What was not established

**Kossel lines were not confirmed on Lambda2M.** After excluding dead pixels and
detector edges, no convincing physical line survived in the 16 datasets sampled.
Possible reasons, not distinguished: Lambda2M subtends only ~2.2° per frame at
2.2 m so intercepts few anvil cones; the sampling was sparse; sensitivity scales
with local intensity. So the algorithm is so far validated on *detector* defects
rather than on the anvil lines it targets.

---

## 4. Suggested interface

Three mask layers with different lifetimes, kept separate:

| layer | lifetime | source |
|---|---|---|
| blemish | permanent, per detector | site TIFF (needs the §3 correction) |
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

1. **Correct the blemish file** (§3) — independent of the feature, affects every
   analysis on Lambda2M.
2. Implement the q-φ detector with border exclusion and the dead-pixel split.
   Validate against the Eiger frame, where the lines are visible by eye.
3. Add the geometry-free fallback for detectors with poor q calibration.
4. Spot detection last — compact bright features are easier, and
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

Data referenced:

```
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/E0157_EHEA-Mesh_a0001_f000010_eiger4M_r00006/
/gdata/dm/8ID/8IDE/2026-3/pope202609/data/G0242_HEA-Mesh30_a0001_f000010_lambda2M_r00227/
/home/beams/8IDIUSER/Documents/areaDetectorBlemish/8idLambda2m/latest_blemish.tif
```

All read-only with respect to the beamline: file reads and compute only.
