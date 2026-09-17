# Proposal: automatic artifact masking in pySimpleMask

**For:** Miaoqi Chu (pySimpleMask), via Q. Zhang
**Date:** 2026-09-17
**Context:** diamond-anvil-cell WAXS/XPCS at 8-ID-E. Patterns carry narrow dark
lines crossing the scattering, plus localised bright spots. Neither is currently
detected automatically. Example: `E0157_EHEA-Mesh_10001_f00010_eiger4M_r00006`.

The proposal rests on one assumption, which the user supplied and which is sound:
**for a powder-like sample the true scattering I(q, φ) is smooth**, so anything
abrupt is an artifact. What follows is how to turn that into a detector, what
breaks if you do it naively, and what a prototype actually found on real data.

---

## 1. What the lines probably are — and they are not Kikuchi lines

**Kikuchi lines are an electron-diffraction phenomenon.** They require inelastically
scattered electrons inside the crystal to act as a divergent source which then
Bragg-diffracts. There is no X-ray Kikuchi effect.

The X-ray analogue is the **Kossel line**, and the geometry here matches it well:

- The sample scatters diffusely into a broad cone — a *divergent source*.
- That cone passes through the **downstream diamond anvil**, a large single crystal.
- Directions satisfying a Bragg condition for some diamond *hkl* are diffracted
  out of their path. The locus of such directions is a cone, which intersects a
  flat detector as a conic section — a narrow, gently curved **dark** line.

In the high-pressure community these are usually just called *diamond Bragg lines*
or *anvil shadows*. They are **subtractive** (intensity removed), which matches the
observation that every line found in this study was dark and none was bright. The
localised bright spots are the complementary effect: a diamond (or a coarse gasket
grain) Bragg reflection landing *on* the detector.

### Three tests that would settle it

| test | if diamond-anvil Bragg |
|---|---|
| Rotate the cell 1–2° about the vertical | lines **sweep rapidly** across the detector; the sample rings do not move |
| Translate the sample across a mesh | lines **stay put** — the anvil orientation is unchanged |
| Change pressure at fixed orientation | lines barely move; sample rings shift |

The first is decisive and nearly free: the ±5° rotations already performed for
sample centring would show it. **Do that before building any physics into the
masking** — if the lines move with cell angle, the mask is valid only for one
orientation, which changes how it must be stored (§4).

---

## 2. Proposed algorithm

### 2.1 The core, in five steps

```
1. build the reference   I_ref(q)  = median over φ within fine q bins
2. residual              R = I - I_ref(q)
3. integrate             R_s = smooth(R)          # lines are EXTENDED
4. robust threshold      σ = 1.4826 · MAD(R_s);  flag |R_s| > n·σ
5. classify              connected components -> aspect ratio
                           high aspect  -> line
                           low aspect   -> spot
```

**Use the median, not the mean,** at step 1. The median is insensitive to the
outliers being hunted, so the reference does not chase the artifact.

**On the user's derivative idea.** The instinct is right — "abrupt" means high
spatial frequency — but a literal finite difference is the wrong estimator here.
At ~0.3 counts/pixel (measured, §3) a derivative is pure shot noise. Step 2–3 is
the same idea implemented as a matched filter: difference from a *smooth reference*
at the feature's own scale, which is noise-optimal. Recommend this over derivatives.

### 2.2 Four failure modes, all hit while prototyping

**(a) The reference must follow the true q contours.** A first attempt used
"median along detector rows" as the reference. On Lambda2M this is nearly right —
the virtual beam centre is ~19,000 rows off-detector, so rings are almost straight
horizontal lines. **On Eiger4M it is completely wrong**: the q gradient runs across
columns there, so a row-median averages straight across the signal and the method
detects nothing. Use the actual q map.

This is also the honest answer to "will this work on Eiger, where the geometry is
uncertain?" — **the algorithm needs the q map to be locally correct, not
absolutely correct.** An error in the beam centre distorts ring shape and smears
the reference. A useful consequence: azimuthal uniformity of the residual is itself
an objective function for fitting the centre, so detection and centre-finding could
be run together. Where geometry is hopeless, fall back to §2.3.

**(b) Exclude detector borders and module-gap edges first.** Without this, edge
columns dominate the detection completely. In a 6-dataset trial, 5 of 6 "strongest
line detections" were at column 1556–1557 of a 1558-wide detector — the last
column, not a physical artifact. Erode the valid-pixel mask by a few pixels before
thresholding.

**(c) Single frames do not have the statistics.** The Eiger example averages
**0.287 counts/pixel** over 10 frames, with a maximum of 9 counts anywhere. The
lines are visible to the eye only because the eye integrates along them. Build the
mask from a **high-statistics sum**: the artifacts are static, so summing a whole
mesh (441 points × 10 frames here) is legitimate and costs nothing.

**(d) Separate dead pixels from physical artifacts — they look identical to a
threshold but are not.** A dead pixel reads **exactly zero**; a Bragg line removes
a *fraction* of the intensity. Test the depth ratio:

```
depth / I_ref  ≈ -100%   ->  dead pixel   (detector defect, permanent)
depth / I_ref  ~  -10..-60%  ->  physical  (anvil line, orientation-dependent)
```

This distinction is not cosmetic: the two have different lifetimes and belong in
different mask layers (§4). Without it, the first thing the detector finds is
unflagged dead pixels — which is exactly what happened (§3).

### 2.3 Geometry-free fallback

Where the q map cannot be trusted, a line is still a line. Detect narrow, extended
features directly in pixel space with a directional matched filter — a Radon/Hough
transform of the residual, or morphological opening with a line structuring element
at several orientations. Less sensitive than the q-φ method, but it needs no
geometry at all. Recommend shipping both, with the q-φ path as default.

---

## 3. What a prototype found on real data

Implemented as above and run over a sample of the 1398 Lambda2M datasets from
pope202609, plus the Eiger example.

### 3.1 A genuine gap in the Lambda2M bad-pixel map

The strongest and most consistent detection was a band at **rows 112–125,
columns 712–828** on Lambda2M:

```
raw mean intensity in the band   0.0045 counts   (detector mean 1.5833)
pixels reading exactly zero      1564 of 1638  (95.5%)
blemish map says GOOD (b == 1)   1553
READ ZERO *and* flagged good     1479 pixels     <-- unmasked dead region
```

Probed at fixed pixel coordinates across eight mesh points, the depth/reference
ratio was **−97% to −100% every time** — the signature of dead pixels, not of a
Bragg line. These 1479 pixels are **not** in
`areaDetectorBlemish/8idLambda2m/latest_blemish.tif` and are therefore included in
every qmap built today.

Also present: contiguous all-zero row bands at **902–907** and **1552–1557**,
beyond the two known module gaps (516–646, 1163–1296).

**This is worth acting on independently of the masking feature** — it is a
correction to the blemish file, and it silently affects every analysis on this
detector.

### 3.2 What was *not* established

**Diamond Bragg lines were not confirmed on Lambda2M.** The detector found dead
pixels and detector edges; after excluding both, no convincing physical line
survived in the datasets sampled. Possible reasons, not distinguished here:

- The Lambda2M sits at 2.2 m and subtends only ~2.2° per frame, so it sees a much
  smaller solid angle than the Eiger and may simply not intercept many anvil cones.
- The sampling was sparse (16 of 1398 datasets at 40 frames each).
- Detection sensitivity scales with local intensity; a fractional line is hardest
  to see exactly where the sample is faint.

So the Eiger observation stands on its own, and the algorithm is validated on
*detector* defects rather than on the anvil lines it was designed for. That gap
should be closed before the feature is called done — §5.

---

## 4. Suggested interface

Three mask layers with different lifetimes, kept separate:

| layer | lifetime | source |
|---|---|---|
| blemish | permanent, per detector | site TIFF (needs the §3.1 correction) |
| **artifact** | **per cell orientation** | **proposed auto-detection** |
| user | per analysis | manual ROI |

A minimal CLI surface consistent with the existing `build` options:

```
--auto-artifact {off,lines,spots,both}     default off
--artifact-nsig N                          default ~6
--artifact-min-aspect N                    line vs spot, default ~5
--artifact-dead-frac F                     depth/ref beyond this = dead pixel,
                                           routed to the blemish layer (default 0.9)
--artifact-from FILE [FILE ...]            build the mask from a high-statistics
                                           sum of these files, apply to the target
--artifact-exclude-border N                default ~8 px
--output-artifact-mask FILE                write the layer separately
```

Two points of principle worth preserving:

- **Report, never silently mask.** Print counts and total area removed, and write
  the layer as its own file. A mask that quietly deletes 5% of the detector is
  worse than the artifact when it is wrong.
- **Default off.** Anything that removes real data on the user's behalf should be
  opt-in until it has a track record.

---

## 5. Suggested order of work

1. **Fix the blemish file** (§3.1). Independent of the feature, affects everything,
   costs nothing.
2. **Run the cell-rotation test** (§1). Decides whether the artifact layer is
   per-orientation, which determines the interface.
3. Implement the q-φ detector with border exclusion and the dead-pixel split.
   Validate against the Eiger frame where the lines are visible by eye.
4. Add the geometry-free fallback for detectors with poor q calibration.
5. Only then consider spot detection — compact bright features are easier and
   partly covered by `--threshold-high` already.

---

## Appendix — reproduction

On `amber`, in `~/xpcs_contrast_check/`:

| script | does |
|---|---|
| `artifact_detect.py` | core detector: reference, residual, classify |
| `fish_lambda.py` | q-binned reference from metadata; scans many datasets |
| `line_persist.py` | tracks a detected line across mesh points |
| `line_fixed_probe.py` | probes fixed pixels across mesh points (§3.1 test) |
| `sum_list.py` | high-statistics sum over a scan |

All of it read-only with respect to the beamline: file reads and compute only.
No motor, attenuator, shutter or PV was touched and no ophyd session was started.
