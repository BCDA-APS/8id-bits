# Rigaku500k as a fixed wide-angle phase monitor at 8-ID-E

**Status:** design proposal, ready for mechanical implementation
**Date:** 2026-09-16
**Context:** high-pressure XPCS on a eutectic high-entropy alloy in a diamond-anvil
cell, following the pope202609 beamtime (2026-09-11 → 2026-09-15)
**Companion:** `EVAL-2026-09-13-DAC-WAXS-geometry-and-protocol.md` — the science
result, why the pattern could not be indexed, and the Lambda2M distance question.
That note and this one share an experiment; read §2.2 and §3 there before acting
on §3.5 here.

Every number in this note is either measured on pope202609 data or computed from
the beamline geometry; the distinction is marked throughout. Reproduction scripts
are listed in Appendix C.

---

## 1. What this detector is for

Three detectors, each with one job:

| detector | role |
|---|---|
| Rigaku3M | SA-XPCS |
| Lambda2M | precision XRD (delta scans) **and** WA-XPCS |
| **Rigaku500k** | **fixed wide-angle monitor: is a phase present, and roughly where** |

The Rigaku500k does not replace the Lambda2M delta scan. It is fixed, so it sees
every accessible reflection on every frame, with no scan and no motion. That is
its whole value: phase appearance and disappearance become visible *live*, during
an XPCS acquisition, instead of requiring a dedicated scan.

The BCC→FCC transition found during pope202609 — a line at d ≈ 2.006 Å present at
10.5 GPa and gone at 20.2 GPa — would have been visible on this detector as it
happened.

**It cannot be positioned accurately, and it will not be moved once mounted.**
Every pixel's q must therefore come from an in-place calibration against a known
standard. Section 4 is the core of this design.

---

## 2. Detector position

### 2.1 Specification for mechanical design

```
Sample-to-detector-centre distance : 350 mm, measured along the central ray
Central ray angle                  : 2θ = 23.0° from the direct beam
Detector face                      : NORMAL TO THE CENTRAL RAY (facing the sample)
                                     — not parallel to the beam axis
Long axis (1024 px)                : in the scattering plane
Plane                              : the same scattering plane as the Lambda2M delta arm

Centre position : 136.8 mm transverse from the beam axis
                  322.2 mm downstream of the sample
Active area     : 77.8 × 38.9 mm  (1024 × 512 px @ 76 µm)
Pixel           : 0.0124° in 2θ
```

### 2.2 Clearance

The detector must clear the Lambda2M flight path, which requires a keep-out cone
of at least 5° about the direct beam.

```
inner active edge   100.9 mm off axis, at 340 mm downstream
5° cone radius there 29.5 mm
CLEARANCE            71.4 mm available for housing, cabling and mount
```

The constraint is not close to binding. At this distance the detector centre
could drop to 2θ ≈ 11.5° before the inner edge reached the cone.

### 2.3 Tolerances — where precision is and is not needed

**Angular tolerance on the mount is not critical.** Because the detector already
views at 23° to the beam, a small additional tilt is degenerate with a change in
the fitted distance and centre angle: the calibration absorbs it. Simulation
(§4.3) shows an in-plane tilt of up to **5°** leaves a residual q error of
0.00008 Å⁻¹ — four hundred times smaller than the features being measured.

What the mount **does** need:

- **Rigidity.** The calibration is taken once, in place. Any creep, thermal walk
  or vibration after calibration is an uncorrectable error. This is where the
  engineering effort belongs.
- **A reproducible sample position.** See §4.4 — the dominant error in the whole
  scheme is not the detector at all, it is where the calibration standard sits
  relative to the sample.

---

## 3. Angular coverage and what it buys

### 3.1 Coverage

```
2θ  16.66 – 29.34°        q  2.233 – 3.904 Å⁻¹
```

> **Cell limit: 2θ ≤ 30°, q ≤ 3.990 Å⁻¹.** The anvil opening spans −30° to +30°,
> but exploiting the full 60° would mean sending the incident beam in along the
> −30° edge of the cone — grazing the anvil, long diamond path, asymmetric
> geometry. That is not a practical configuration. With the beam on the cell axis
> the usable maximum scattering angle is **30°**.
>
> This window reaches **within 0.66° of that limit** — essentially nothing the
> cell can deliver is left uncovered.

> ⚠ **This contradicts the companion note, which should be corrected.**
> `EVAL-2026-09-13-…` §2.2 states the cell allows 60°, lists "six accessible
> reflections" (111, 200, 220, 311, 222, 400), and concludes "the indexing was
> solvable the whole time and nobody looked." **That conclusion does not hold.**
> At 30° a 3.456 Å FCC cell gives only **two** reflections — 111 at 23.6° and
> 200 at 27.3°. The failure to index was an energy/aperture limitation, not a
> scan-range oversight, and the remedy is higher energy or a wider-aperture cell,
> not longer scans.

Nothing of interest lies below 16.7°.

### 3.2 Peaks covered

| feature | 2θ | margin to nearest edge |
|---|---|---|
| CeO2 200 *(calibration)* | 17.33° | 0.67° (54 px) |
| **Au 111 @ 0 GPa** | 19.94° | 3.28° (264 px) |
| **Au 111 @ 28.5 GPa** | 20.76° | 4.10° (330 px) |
| **HEA low-q peak** | 22.00° | 5.34° (430 px) |
| **vanished BCC line** (d = 2.006 Å) | 23.45° | 5.89° (474 px) |
| CeO2 220 *(calibration)* | 24.60° | 4.74° (381 px) |
| **HEA high-q peak** | 25.00° | 4.34° (349 px) |
| CeO2 311 *(calibration)* | 28.93° | 0.41° (33 px) |

Gold 111 positions are from a Birch–Murnaghan EOS (K₀ = 167 GPa, K₀′ = 5.79) on
a₀ = 4.0782 Å, and agree with the measured ladder in the companion note (19.962°
at ambient, 20.816° at 30.8 GPa). The gold marker never leaves the detector across
the full 0–30 GPa range, every science feature sits at least 3.3° from an edge,
and the two outermost calibration rings being near the edges is acceptable — §4.2
shows two rings suffice.

### 3.2a A real bonus: Au 111 and Au 200 in the same frame

This window also contains **Au 200 at 23.1°** alongside Au 111 at 19.9°. That is
worth more than it first appears.

For FCC the ratio q₂₀₀/q₁₁₁ = 2/√3 is exact and known a priori. With the distance
already fixed by CeO2, having both lines in one exposure gives **two independent
determinations of a_Au at the sample position, on every frame** — a redundant
pressure readout and a continuous consistency check.

It does **not** remove the need for CeO2. The gold angles are not known a priori
(a_Au is the unknown), and §4.0 shows the scale degeneracy is exact: gold alone
determines neither the distance nor the lattice parameter. The pair is a check on
the calibration, not a substitute for it.

The companion note (§3.2) identifies exactly this — Au 111 and Au 200 in a single
frame — as the main prize of relocating Lambda2M to 1.0 m, a route of unconfirmed
feasibility (§3.4). **The Rigaku500k at 350 mm delivers it without moving Lambda2M
at all**, and continuously rather than in dedicated scans. For a 441-point mesh
this turns one frame per point into a simultaneous phase map *and* pressure map —
with the accuracy caveat in §4.0: read absolute pressure from Lambda2M, and use
the Rigaku500k to map and track it.

> **Crowding warning.** Au 200 moves from 23.07° (ambient) to 24.03° (30 GPa),
> straight through the region occupied by the vanished BCC line (23.45°) and close
> to CeO2 220 (24.60°). At ambient, Au 200 and the BCC line are 0.38° apart —
> comparable to the 0.377° sample-peak FWHM, so they will overlap. They remain
> separable because gold is textured and the sample is not (spottiness 5.18 vs
> ≤0.14): resolve them azimuthally, not radially.

> Caveat: the two strong observed peaks (22.0° and 25.0°) do **not** index cleanly
> as FCC 111/200 — their q ratio is 1.134 against the required 1.1547, off by 1.8%.
> That remains unresolved (companion note §1.3). The window covers both regardless,
> and also covers where the paper's FCC cell (a = 3.456 Å at 30.6 GPa) puts 111
> and 200, at 23.57° and 27.29°.

### 3.3 Angular resolution: can it track the gold 111 shift?

**Yes, comfortably.** The gold peak width was measured directly, not assumed:
scan **S00247** (Lambda2M delta scan, gold 111 at ~10 GPa, peak at delta 20.312°
— matching the EOS prediction of 20.30°) gives

```
gold 111 FWHM = 0.167°       peak/background = 7.0
```

At 350 mm that peak spans **13.4 pixels**. Against it:

| pressure range | 5 GPa shift | pixels | vs FWHM | 1 GPa |
|---|---|---|---|---|
| low P (0–5 GPa) | 0.186° | 15.0 | 1.12 × FWHM | 3.0 px |
| high P (25–30 GPa) | 0.112° | 9.0 | 0.67 × FWHM | 1.8 px |

The gold peak moves by roughly its own width per 5 GPa step. Even 1 GPa is two to
three pixels. Combined with the systematic budget in §3.4, a 5 GPa step is
detected at about **13:1**.

> **Caveat.** The detector subtends only ~16° of azimuth. Gold is the textured
> phase — azimuthal spottiness measured 5.18 during pope202609, against ≤0.14 for
> the sample — so a short arc samples few grains and the centroid can be biased by
> which grains happen to diffract. Use the Rigaku500k to *track* the gold peak;
> use Lambda2M delta scans to *measure* it. The sample peaks, being smooth, are
> better served by this detector than gold is.

### 3.4 Why not move the detector further back

Backing off gives finer pixels, but pixel size is not the limiting term and the
trade is net negative:

| L | span | px across gold FWHM | azimuth | photon σ | displacement σ | total σ |
|---|---|---|---|---|---|---|
| 300 mm | 14.78° | 11.5 | 18.2° | 0.0021° | 0.0095° | 0.0098° |
| **350 mm** | **12.69°** | **13.4** | **15.7°** | **0.0022°** | **0.0082°** | **0.0085°** |
| 500 mm | 8.90° | 19.2 | 11.0° | 0.0027° | 0.0057° | 0.0063° |
| 700 mm | 6.36° | 26.8 | 7.9° | 0.0032° | 0.0041° | 0.0052° |

Two reasons the gain is illusory:

1. **13 pixels across a peak is already plenty** for a centroid fit. Subdividing
   further does not add information.
2. **Moving back collects fewer photons, not more.** The captured azimuthal arc
   shrinks as 1/L (15.7° → 7.9°), so total counts fall and centroid precision
   degrades as √L even while pixels get finer.

The only term that improves is sample-displacement immunity, and at 350 mm that
is already seven times below the signal. Against the 0.1116° worst-case 5 GPa
shift, 350 mm gives 13:1 and 700 mm gives 21:1 — bought by halving the angular
coverage and losing every calibration ring but one.

### 3.5 Withdrawn: there is no useful high-angle placement

An earlier draft of this note proposed an alternative placement at 2θ = 41° to
capture FCC 220, 311 and 222 for indexing. **That option does not exist.** It
rested on the 60° aperture figure; at the real 30° limit there is nothing above
30° to detect.

The consequence is worth stating plainly: **no additional detector and no longer
scan can solve the indexing problem.** The reflections required are not accessible
at 15.2 keV in this cell. Only higher energy or a wider-aperture cell can reach
them. Companion note §4 argues against raising the energy on other grounds; that
trade needs revisiting now that the scan-range escape route is closed.

---

## 4. Calibration with CeO2

This is the part the design lives or dies on. The detector cannot be positioned
accurately and cannot be moved, so **every pixel's q comes from an in-place fit to
a known standard**, taken once, at the science position.

### 4.0 Can the gold pair replace CeO2? No.

Since Au 111 and Au 200 both land on this detector (§3.2a), it is natural to ask
whether CeO2 is still needed. **It is.** The two are not interchangeable, and the
reason is a hard degeneracy rather than a question of precision.

CeO2 ring angles are **known absolutely** — the lattice parameter is certified and
fixed. Gold's are **not**: a_Au depends on the pressure, which is the unknown being
measured. Only the ratio q₂₀₀/q₁₁₁ = 2/√3 is known a priori, and that ratio is
scale-invariant in exactly the wrong way. Scaling the distance and the lattice
parameter together, (L, a) → (kL, ka), leaves the pattern on the detector
unchanged to first order in 2θ.

Simulated, fitting `(L, c, v0, a_Au)` to both gold conics:

```
correlation(L, a_Au) = +1.0000        fit runs away to L = 639 mm, a_Au = 7.35 Å
```

Perfectly degenerate. The tan/arcsin nonlinearity breaks it only at second order,
far below the noise at these angles. **Gold alone cannot determine either the
distance or the lattice parameter.** CeO2, with three rings of known angle,
determines L to 0.002%.

**What the gold pair is genuinely worth,** once L is fixed by CeO2:

1. **Relative pressure tracking — excellent.** Peak *shifts* are geometry-free;
   this is the detector's main job and needs no absolute scale.
2. **A coarse consistency check.** With L fixed, each gold line yields its own
   a_Au; they agree only if the geometry is right. But the test is weak — a 1 mm
   sample-position error produces a 111/200 split of just 0.00031 Å, about
   **0.06 pixel**, which is below the grain-sampling bias on a spotty gold ring
   (azimuthal spottiness 5.18). Treat it as an alarm for millimetre-scale
   blunders, not as a precision transfer standard.

**And a warning about using this detector for pressure.** Because δq/q = δz/L,
the short 350 mm arm is unforgiving: a **1 mm sample-position error biases a_Au by
0.0095 Å ≈ 1.2 GPa**. At Lambda2M's 2.2 m the same error costs ~0.2 GPa. Read
pressure from Lambda2M; use the Rigaku500k to watch it change.

### 4.1 What is fitted

Three parameters:

| parameter | meaning |
|---|---|
| `L` | effective sample-to-detector distance |
| `c` | 2θ of the detector centre |
| `v0` | out-of-plane offset of the detector centre from the scattering plane |

Energy is held fixed — a prior CeO2 measurement at 8-ID-E established it is known
well. Tilts are **not** fitted; §4.3 shows they are degenerate with `L` and `c`
and are absorbed harmlessly.

Because the detector faces the sample, it lies ~23° to the beam axis, so the
Debye–Scherrer rings land on it as strongly-curved **conic sections**, not
near-straight arcs. That curvature is what breaks the distance/centre degeneracy,
and it is the reason a short azimuthal arc is still enough to calibrate.

### 4.2 CeO2 rings and how many are needed

NIST SRM 674b, a = 5.411651 Å. At 15.209 keV (λ = 0.81521 Å):

| reflection | d (Å) | q (Å⁻¹) | 2θ | on this detector |
|---|---|---|---|---|
| 111 | 3.1244 | 2.0110 | 14.99° | no (below window) |
| **200** | 2.7058 | 2.3221 | **17.33°** | **yes** |
| **220** | 1.9133 | 3.2839 | **24.60°** | **yes** |
| **311** | 1.6317 | 3.8508 | **28.93°** | **yes** |
| 222 | 1.5622 | 4.0220 | 30.25° | no (above this window; the cell reaches it) |

Simulated fit precision, 0.2 px ring-finding scatter:

| rings used | L precision | resulting q error at q = 3.3 |
|---|---|---|
| 4 | 0.002% | 0.00005 Å⁻¹ |
| **3 (this design)** | **0.002%** | **0.00007 Å⁻¹** |
| 2 | 0.005% | 0.00018 Å⁻¹ |
| 1 | 0.204% | 0.0067 Å⁻¹ |

Three rings calibrate ~700× finer than the peak width being measured. **Two would
suffice**, which is why the window was placed for coverage rather than for ring
count. One ring is not enough — with a single ring the distance and centre become
degenerate and precision falls by two orders of magnitude.

Recompute this table if the energy changes. A 0.06% energy shift moves the rings
by ~0.015° (about one pixel) — negligible, but the ring table should still be
regenerated rather than assumed.

### 4.3 Error budget

Fitting `(L, c, v0)` only, and transferring the result to the sample:

| error source | residual q error | verdict |
|---|---|---|
| perfect | 0.00005 Å⁻¹ | — |
| in-plane detector tilt 0.5° | 0.00010 Å⁻¹ | absorbed |
| in-plane detector tilt 2.0° | 0.00005 Å⁻¹ | absorbed |
| in-plane detector tilt **5.0°** | 0.00008 Å⁻¹ | **absorbed** |
| out-of-plane tilt 1.0° | 0.0013 Å⁻¹ | negligible |
| **CeO2 0.5 mm off in z** | **0.0061 Å⁻¹** | marginal |
| **CeO2 1.0 mm off in z** | **0.0123 Å⁻¹** | ≈ the whole 5 GPa gold shift |
| **CeO2 2.0 mm off in z** | **0.0246 Å⁻¹** | half a peak width — unacceptable |

*(For scale: sample peak FWHM = 0.0495 Å⁻¹; 5 GPa gold shift = 0.0155 Å⁻¹ at q = 3.3.)*

### 4.4 The one error that matters, and why it is dangerous

**Everything hinges on where the CeO2 standard sits along the beam relative to the
DAC sample.**

A z-offset is not merely large — it is **invisible**. A 1 mm offset causes the fit
to return L low by 0.92 mm with a clean, convincing residual. Nothing in the
calibration output indicates a problem, and the q scale is then wrong by 0.4%
everywhere. A previous CeO2 calibration at this beamline came out ~1% low in q;
an undetected offset of this kind is the most likely explanation.

Two defences, both recommended:

1. **Co-locate the standard.** Mount CeO2 inside the DAC, or in a dummy cell of
   the same geometry, and record the sample-stage z for both. This should hold the
   offset inside 0.5 mm.
2. **Use gold as a transfer standard — but cross-referenced to Lambda2M, not on
   its own.** The gold is inside the cell at exactly the sample position, so it
   carries the right z. It cannot supply an absolute angle by itself (§4.0), so
   pin gold 111's 2θ with a Lambda2M delta scan — which the experiment does anyway
   — and apply a single scalar correction to the fitted `L` so the Rigaku500k
   agrees. This removes the z-offset, and it is the only check capable of catching
   it at the precision that matters.

Use CeO2 for the *shape* of the q map (centre, curvature, pixel→2θ) and gold for
the *absolute scale*. Neither alone is sufficient.

### 4.5 Procedure

1. Mount and align the detector. Do not move it again for the rest of the beamtime.
2. Insert the CeO2 standard at the sample position; record the stage coordinates.
3. Collect a CeO2 pattern with all three rings well exposed.
4. Find ring centroids azimuthally; fit `(L, c, v0)` against the known 2θ values.
   Inspect the residual — with three rings it is over-determined, so a
   structured residual means the model is wrong.
5. Install the DAC. Record the stage coordinates and compare to step 2.
6. Measure gold 111 on Lambda2M; scale `L` so the Rigaku500k agrees.
7. Write the calibration out **with the run**, alongside the geometry it is valid
   at. A q map that does not record the geometry it was generated under is the
   exact failure that cost time during pope202609.
8. Re-verify against gold after any intervention that could have disturbed the
   mount.

---

## 5. Measured XPCS contrast on Lambda2M

Recorded here because it sets expectations for what WA-XPCS on Lambda2M can
deliver, and because the analysis contains a trap worth documenting.

### 5.1 Result

**Dataset:** `G0256_HEA-Dec9GPa_a0001_f003000_lambda2M_r00001`
att 1, 3000 frames × 1 s, Lambda2M middle module, sample decompressed to ~9 GPa.
Sample peak (not gold) at **q = 3.288 Å⁻¹, 2θ = 24.65°**.

Analysis: one dynamic ROI spanning exactly one peak FWHM (q 3.2632–3.3126 Å⁻¹),
25 static ROIs, 1405 hot/dead pixels (0.175%) removed by an 8σ MAD filter.

```
g2(τ = 1 s)                 = 1.0191
static artifact floor       = 1.0091      <- NOT dynamics
BETA (real contrast)        = 0.0096 ± 0.0001
τc                          = 500 ± 29 s
γ                           = 1.62 ± 0.14
```

**β ≈ 0.010**, with slow compressed-exponential relaxation on a ~500 s timescale.

*Uncertainty.* The decay is incomplete within 3000 s, so β, τc and the baseline are
correlated. Fixing the floor at its independently measured value gives β = 0.0096;
floating it gives β = 0.0148 with a floor of 1.0042 — below the measured static
variance, which is unphysical. **Take β ≈ 0.010 with a systematic uncertainty of
roughly +0.005**, not the ±0.0001 statistical error.

### 5.2 The trap: fake contrast from intensity variation

`boost_corr` averages pixels within each **static** q bin *before* dividing
(`average_with_index_map`, in `compute_g2`). Therefore

> **g2(∞) = ⟨Ī²⟩ₚ / ⟨Ī⟩ₚ² = 1 + Var_p(Ī)/⟨Ī⟩²**

Any spatial variation of mean intensity inside a static bin lifts the g2 baseline
and reads as contrast. Predicted from the time-averaged image, then measured:

| ROI / binning | predicted baseline | measured |
|---|---|---|
| tight ROI + 0.002 Å⁻¹ static bins | 1.0091 | 1.0056 |
| tight ROI, one static bin | 1.0163 | 1.0130 |
| whole middle module, one bin | 1.0654 | 1.0620 |

Quoting `g2(0) − 1` from the whole-module ROI would give **0.074 — an 8.6×
overstatement**, almost entirely spatial structure rather than dynamics.

**Rules.** Keep the dynamic ROI to about one peak FWHM; use ≥10 static bins so
each is near-iso-intensity; remove hot pixels; and always quote the *decay
amplitude above an independently determined floor*, never `g2(0) − 1`.

### 5.3 Controls

- **Per-frame normalisation** (`-nf 1`) gave an identical curve — not beam drift
  (incident intensity moved 1.2% over the run).
- **Half-length rerun** (1500 frames) overlaid the full run *at the same τ*
  (1.0130 vs 1.0136 at τ = 256 s) rather than scaling with τ/T — so the decay is
  real decorrelation, not finite-run-length bias.
- Points at τ > 1000 s fall below the static floor, indicating residual downward
  bias; they were excluded from the fit.

### 5.4 Why β is only 1%

Coherence, not the sample. At 2θ = 24.65° the path-length spread across the sample
(~1.8 µm for 20 µm thickness) far exceeds the longitudinal coherence length
(~0.3 µm at ΔE/E ≈ 1.4 × 10⁻⁴), and the ~18 µm speckle is well under the 55 µm
pixel. Both are addressable — a thinner sample, a lower-q reflection, or finer
pixels — and should be weighed before committing to WA-XPCS on Lambda2M.

### 5.5 Output files

```
/gdata/dm/8ID/8IDE/2026-3/pope202609/analysis/Multitau/
    G0256_HEA-Dec9GPa_a0001_f003000_lambda2M_r00001_TightPeakROI_results.hdf
    qmap_G0256_TightPeakROI.hdf
```
Readable directly in pyXpcsViewer.

---

## Appendix A — measured quantities

| quantity | value | source |
|---|---|---|
| beam energy | 15.19991 keV (λ = 0.81569 Å) | G0256 metadata |
| Lambda2M distance / pixel | 2.2 m / 55 µm | G0256 metadata |
| Lambda2M beam centre | col 779, row 987 | G0256 metadata |
| Lambda2M module gaps | rows 516–646, 1163–1296 | G0256 average image |
| gold 111 FWHM | 0.167° | S00247 delta scan |
| sample peak, direct | q 3.288 Å⁻¹, 2θ 24.65°, FWHM 0.377° | G0256 average image |
| sample peak, via delta scan | FWHM 0.473° / 0.612° | S00234 |
| hot/dead pixels | 1405 of 803928 (0.175%) | G0256, 8σ MAD |
| azimuthal spottiness | gold 5.18, sample ≤0.14 | pope202609 |

**Methodological note.** The delta-scan FWHM (0.473°) is *broader* than the
detector-resolved width (0.377°) because the scan integrates a ROI while the arm
rotates, convolving the peak with the ROI acceptance. **For peak widths, resolve
the peak on the detector; use the delta scan for peak position.** Using the
rocking-curve width to size a q ROI would over-wide it by ~25%.

## Appendix B — what was *not* verified

- Physical clearance against the real flight-tube envelope and Rigaku3M. The 5°
  cone was the only constraint applied; the actual hardware model was not checked.
- Whether the Rigaku500k is a single continuous chip. A module seam was assumed
  absent; if one exists, keep it away from 2θ ≈ 20–25°.
- Gold EOS consistency: a BM3 EOS puts Au 111 at q = 2.75 at 20 GPa, while a
  pope202609 note recorded a gold feature at q = 2.856 there. Both lie inside the
  detector window so the geometry is unaffected, but the pressure scale deserves
  resolution.
- **Where the cell shadow actually cuts in.** 30° is the nominal usable maximum
  with the beam on the cell axis; the real edge, after seat and gasket, was not
  measured. One wide CeO2 exposure would locate it and is worth taking early —
  the outermost calibration ring (311 at 28.93°) sits only 1.1° inside it.

## Appendix C — reproduction

In this directory:

| script | produces |
|---|---|
| `rigaku500k_geometry.py` | coverage, CeO2 ring landings, error budget |
| `rigaku500k_optimum.py` | distance optimisation, 5 GPa resolution |
| `ceo2_calib_sim.py` | calibration fit precision vs ring count |
| `ceo2_systematics.py` | tilt and z-offset error budget (§4.3) |

On `amber`, in `~/xpcs_contrast_check/`: qmap builders, the three comparison
`boost_corr` runs, and the fake-contrast diagnostic of §5.2.

All pope202609 analysis was read-only — no motor, attenuator, shutter or PV was
touched, and no ophyd session was started (session startup fires an Eiger frame).
