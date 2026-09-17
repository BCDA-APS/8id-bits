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

- **Rigidity.** The calibration is taken once, in place, so any creep or thermal
  walk afterwards is uncorrectable and invisible. The planned mount — **bolted to
  a post** — is the right answer; this note assumes it and drops creep from the
  error budget accordingly (§4.0). Keep it that way: do not add adjusters that
  could drift or be nudged.
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

### 3.2a Marker choice — switch from gold to tantalum

An earlier draft of this note treated "Au 111 and Au 200 in the same frame" as a
bonus, and proposed separating Au 200 from the sample lines azimuthally on the
grounds that gold is textured and the sample is not. **Both claims are withdrawn.**

**Gold texture disappears once pressure is applied** (observed on pope202609
data). The ambient spottiness of 5.18 does not survive compression, so azimuthal
discrimination is not available under experimental conditions. Au 200 must
therefore be separated radially or not at all — and it cannot be:

| P (GPa) | Au 200 | sample FCC 111 | separation |
|---|---|---|---|
| 0 | 23.06° | 22.68° | 0.39° |
| 10 | 23.47° | 22.96° | 0.50° |
| 20 | 23.79° | 23.22° | 0.57° |
| 30 | 24.06° | 23.44° | 0.61° |

Against a sample-peak FWHM of 0.377°, the two overlap at half maximum across the
entire pressure range and **never separate**. This is not a crossing to be worked
around; it is a permanent blend. It corrupts the sample FCC 111 measurement, so it
is a science problem, not merely a calibration one.

#### Marker screen

Every in-window line of each candidate, tracked over 0–30 GPa against the sample
band (22.68–27.14°, plus a 0.4° guard):

| marker | clean line | interfering line | 5 GPa shift, worst case |
|---|---|---|---|
| Au gold | 111 (19.94–20.80°) | **200 fouls the sample** | 9.0 px |
| Pt platinum | 111 (20.73–21.35°) | **200 fouls the sample** | 7.0 px |
| **Ta tantalum** | **110 (20.08–20.94°)** | **none** | **9.9 px** |
| W tungsten | 110 (20.99–21.58°) | none | 7.0 px |
| Mo molybdenum | 110 (21.11–21.78°) | none | ~7 px |

**Recommendation: tantalum.** It is better than gold on every axis that matters
here:

- **No interfering line.** Ta is BCC, so its second reflection (200) lands at
  28.55–29.31° — inside the *upper* clean zone, not on the sample. Gold and
  platinum, both FCC, put their 200 squarely in the sample band.
- **Better pressure sensitivity than gold** at the demanding end: 9.9 px per
  5 GPa at high pressure against gold's 9.0.
- **A second clean line below 17.8 GPa.** Ta 200 gives an in-situ check on
  detector drift (§4.0b) — the largest remaining error term, and one nothing else
  can catch once the cell is installed. W and Mo have no second line in the window
  at all.
- Ta 110 (20.1–20.9°) is far from the sample BCC 110 (~23.2°); no confusion there.

W and Mo also avoid fouling but give one line only and ~30% less sensitivity.

**Costs, which are real:** switching requires a fresh cell loading; the Ta
pressure scale (Cynn & Yoo 1999, K₀ = 194 GPa, K₀′ = 3.5) is well established but
less canonical than the gold scale, so pressures will not be directly comparable
to Au-referenced literature without conversion; and Ta 200 leaves the window above
17.8 GPa, so the in-situ drift check of §4.0b is available only below that (Ta 110,
and therefore pressure measurement, is unaffected at all pressures). Chemical
compatibility with the helium medium is settled in §4.0c — the sealed gasket
volume cannot supply enough hydrogen to hydride Ta at any realistic gas purity.

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

### 4.0 CeO2 is measured at ambient, in a DAC-like holder

**CeO2 can only be used at ambient** — pressurised, it would compress and stop
being a standard. It is therefore loaded into a **DAC-like holder containing no
sample**, measured separately from the experiment.

Two consequences, both favourable:

- **The anvils are in the beam during calibration.** Attenuation, aperture and
  vignetting match the real experiment, so the cell shadow is characterised by the
  calibration rather than being a surprise afterwards.
- **The residual z-offset is small by construction**, and what remains is measured
  directly — see below.

**The marker still cannot replace CeO2** (§4.0a): the marker sits at the sample
position but its lattice parameter is unknown, while CeO2 has a certified lattice
parameter but is in a different holder. Neither alone suffices.

#### The z-offset is solved procedurally, not left as an error

Centre the object on the diffractometer rotation axis, then measure its apparent
transverse position at ω = −5°, 0°, +5°. For a displacement (dx, dz) from the axis,
x(ω) = x₀ + dx·cos ω + dz·sin ω, so

```
dz = [x(+5°) − x(−5°)] / (2 sin 5°)          lever = 5.74
```

The small angle works in your favour: **1 mm of dz shows up as a 174 µm transverse
split**, which is easy to see. Iterate to convergence.

| transverse centring | dz determined to | q error at 350 mm |
|---|---|---|
| 2 µm | 16 µm | 0.00015 Å⁻¹ |
| 5 µm | 41 µm | 0.00038 Å⁻¹ |
| 20 µm | 162 µm | 0.0015 Å⁻¹ |

Even sloppy centring leaves the offset negligible. **This removes what §4.3 called
the dominant systematic.** Applying the same procedure to the DAC and to the CeO2
holder, and differencing, gives the offset between them directly.

#### Error budget as built

The Rigaku500k is **bolted to a post**, so the creep term an earlier draft called
dominant does not apply. Diffractometer positioning is trusted: **0.03% precision
has been demonstrated** in practice with z alignment and huber-delta scans on
Lambda2M.

| source | q error at q = 3.3 | % of the 5 GPa signal |
|---|---|---|
| sample moves 50 µm between pressure points | 0.00114 Å⁻¹ | 7.3% |
| demonstrated alignment, 0.03% | 0.00099 Å⁻¹ | 6.4% |
| ring centroid, photon-limited | 0.00031 Å⁻¹ | 2.0% |
| CeO2 z-offset, 20 µm residual | 0.00019 Å⁻¹ | 1.2% |
| **quadrature total** | **0.00155 Å⁻¹** | **10.0%** |

**A 5 GPa step is measured at 10:1, and a 1 GPa step at 2:1.** The dominant term
is no longer anything about the detector or its calibration — it is the sample
moving when pressure changes, which is a property of the cell.

> **One caveat on transferring the 0.03%.** Whether that figure carries over to a
> 350 mm arm depends on what limits it. Energy and angular-readback errors are
> L-independent and transfer unchanged. An *absolute distance or sample-position*
> error scales as 1/L and would become **6.3× worse** — 0.03% at 2.2 m would
> become 0.19%, or 0.0062 Å⁻¹.
>
> The arithmetic says this is not a concern: 0.03% at 2.2 m corresponds to a
> 0.66 mm position error *if* it were distance-limited, and the ±5° procedure pins
> dz 4–30× better than that. So the 0.03% was almost certainly limited by energy,
> readback or peak fitting rather than by position, and should transfer intact.
> Worth confirming once against CeO2 on the 500k itself rather than assumed.

### 4.0b Is the second marker line required? No.

**Note on indices: tantalum is BCC, so there is no Ta 111** — body-centring
forbids reflections with h+k+l odd. The available lines are **Ta 110** (strong,
20.08–20.94° over 0–30 GPa) and **Ta 200** (28.55–29.31°).

An earlier draft said to "transfer early, before Ta 200 leaves the window."
**That instruction is withdrawn.** It assumed the two-line transfer was needed to
cancel the CeO2 position error; the ±5° centring procedure (§4.0) removes that
error directly, so **Ta 110 alone is sufficient** for pressure measurement and
peak tracking across the whole range.

What the second line is still worth:

| | |
|---|---|
| Ta 200 in the detector window | up to **17.8 GPa** |
| Ta 200 within the cell aperture | up to 36.2 GPa |

Below 17.8 GPa, Ta 110 and Ta 200 must give the same lattice parameter. If they
disagree, the geometry has moved — the only in-situ geometry check available once
the cell is in and CeO2 is gone. With the detector bolted to a post this is a
reassurance rather than a necessity, but it costs nothing and would catch a knock
or a bumped cable. Above 17.8 GPa you keep full pressure sensitivity on Ta 110 and
simply lose the check.

So: **not required.** It is a mild preference for Ta over W and Mo, whose second
lines never enter the window at all — but with creep off the table, the choice
between the three now rests on sensitivity (Ta 9.9 px vs W/Mo ~7.0) and on the
chemistry of §4.0c.

### 4.0c Helium pressure medium with tantalum — no objection

**Helium is not a problem, and is a positive.** It is a noble gas and forms no
compound with tantalum at DAC pressures and room temperature. He penetration is
documented for open-framework materials (silica, zeolites, ice), not for
close-packed or body-centred metals; and porosity is irrelevant to a diffraction
marker in any case, since the measurement is of the lattice, not bulk density.
Because He is the most hydrostatic medium available, it *reduces* the deviatoric
stress that would otherwise bias a relatively ductile marker like Ta. He and Ta
are a good pairing.

**The hydride concern is withdrawn.** Tantalum is a strong hydride former, so an
earlier draft flagged H₂ contamination as a risk. Quantifying it settles the
question: **the gasket hole is a sealed micro-volume, so the hydrogen inventory is
bounded at loading regardless of gas purity.**

A 100 µm × 35 µm hole is **0.275 nL** — about 1 × 10⁻¹¹ mol of gas. Against a
20 × 20 × 5 µm Ta flake (1.1 × 10¹⁴ atoms):

| H₂ in the He | H atoms available | H/Ta | |
|---|---|---|---|
| 1 ppm | 9.9 × 10⁹ | 9 × 10⁻⁵ | safe |
| 10 ppm | 9.9 × 10¹⁰ | 9 × 10⁻⁴ | safe |
| 100 ppm | 9.9 × 10¹¹ | 9 × 10⁻³ | safe |
| 1000 ppm | 9.9 × 10¹² | 9 × 10⁻² | marginal |

TaHx shifts the lattice measurably at H/Ta ≈ 0.1–1.0. Research-grade He is 5N–6N,
with H₂ a fraction of the ≤10 ppm total impurity, putting the ratio at ≤10⁻³ —
**two orders of magnitude below the threshold**. Only implausibly dirty gas, around
0.1% H₂, would approach it.

The physical reason: a furnace or an H₂ atmosphere is an effectively infinite
hydrogen reservoir and *can* hydride tantalum. A sealed sub-nanolitre bubble of
clean helium cannot — there simply are not enough atoms in it.

With a dedicated He gas-loading system, **use tantalum without reservation.**

> Two second-order notes. A smaller marker raises the ratio proportionally — a
> 10 × 10 × 3 µm flake at 100 ppm reaches 6 × 10⁻², marginal — so do not use a
> vanishingly small flake with questionable gas. And Ta carries a thin native
> Ta₂O₅ layer, which is amorphous and does not affect the diffraction lines.

**Tungsten remains an alternative**, not a fallback: no hydride chemistry to think
about at all, at the cost of sensitivity (7.0 vs 9.9 px per 5 GPa) and with no
second line in the window. Choose it only if W is already on hand.

### 4.0a Can the marker replace CeO2? No.

Since the marker sits at the sample position, it is natural to ask whether CeO2 is
needed at all. **It is.** The two are not interchangeable, and the reason is a hard
degeneracy rather than a question of precision.

CeO2 ring angles are **known absolutely** — the lattice parameter is certified and
fixed. The marker's are **not**: its lattice parameter depends on the pressure,
which is the unknown being measured. Only the ratio q₂₀₀/q₁₁₁ = 2/√3 is known a priori, and that ratio is
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

**What the marker is genuinely worth,** once L is fixed by CeO2:

1. **Relative pressure tracking — excellent.** Peak *shifts* are geometry-free;
   this is the detector's main job and needs no absolute scale.
2. **A two-line transfer that removes the CeO2 position error entirely** — but
   only where two clean marker lines are visible. Simulated: calibrating with CeO2
   displaced by up to 5 mm and then re-fitting `(L, centre)` against **both**
   marker lines recovers the true geometry to **0.00004 Å⁻¹**, i.e. completely.
   A **single** line cannot do this — a z-offset biases two parameters and one
   line gives one constraint; single-line correction of a 1 mm offset only
   improves the error from 0.0105 to 0.0085 Å⁻¹. This is the strongest argument
   for tantalum: Ta 110 and Ta 200 are both clean below ~15 GPa, whereas gold's
   second line is unusable at any pressure (§3.2a).

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

1. **Measure it with the ±5° rotation procedure (§4.0).** This is the primary
   defence and it is procedural, not mechanical: it pins dz to tens of microns
   even with sloppy centring, leaving ≤0.0015 Å⁻¹. Apply it to both the DAC and
   the CeO2 holder and difference the results. The DAC-like holder keeps the raw
   offset small to begin with, so the correction is small and the iteration
   converges quickly.
2. **Use the two marker lines as a cross-check, not as the primary fix.** Below
   17.8 GPa, Ta 110 and Ta 200 must yield the same lattice parameter; disagreement
   means the geometry moved. Simulation shows a two-line re-fit can absorb a CeO2
   offset of up to 5 mm to 0.00004 Å⁻¹, so it is a strong backstop if the centring
   is ever in doubt — and it is the only in-situ detector-creep check available
   once the cell is installed (§4.0b).

Use CeO2 for the *shape* of the q map (centre, curvature, pixel→2θ) and gold for
the *absolute scale*. Neither alone is sufficient.

### 4.5 Procedure

1. Mount and align the detector. Do not move it again for the rest of the beamtime.
2. Insert the CeO2 standard on its matched-z holder (§4.4); record stage
   coordinates. CeO2 is measured at ambient, outside the cell — the cell is not in
   the beam for this step, so the rings are unattenuated and unvignetted. Note
   where the cell shadow later cuts in.
3. Collect a CeO2 pattern with all three rings well exposed.
4. Find ring centroids azimuthally; fit `(L, c, v0)` against the known 2θ values.
   Inspect the residual — with three rings it is over-determined, so a
   structured residual means the model is wrong.
5. Install the DAC. Record the stage coordinates and compare to step 2.
6. **At low pressure, before compressing:** re-fit `(L, centre)` against both
   marker lines (Ta 110 and Ta 200). This cancels any residual CeO2 position
   error. With a single-line marker, pin it against a Lambda2M delta scan instead.
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
| azimuthal spottiness, **ambient** | gold 5.18, sample ≤0.14 | pope202609 |
| gold texture under pressure | **disappears** — no azimuthal discrimination | pope202609 |

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
  the outermost calibration ring (311 at 28.93°) sits only 1.1° inside it, and
  Ta 200 sits closer still.
- **Chemical compatibility of tantalum with the alloy itself.** The helium
  question is settled quantitatively in §4.0c, but Ta-versus-HEA contact chemistry
  was not examined. Low risk at room temperature, worth a moment's thought before
  loading.
- Literature confirmation of Ta as a pressure marker in He-loaded cells. The
  reasoning in §4.0c is self-contained, but search tools are blocked by policy
  here so no reference was consulted.
- Whether tantalum texture also vanishes under pressure. Assumed irrelevant, since
  the Ta recommendation does not rely on azimuthal discrimination — but the same
  surprise that invalidated the gold argument could apply.

## Appendix C — reproduction

In this directory:

| script | produces |
|---|---|
| `rigaku500k_geometry.py` | coverage, CeO2 ring landings, error budget |
| `rigaku500k_optimum.py` | distance optimisation, 5 GPa resolution |
| `ceo2_calib_sim.py` | calibration fit precision vs ring count |
| `ceo2_systematics.py` | tilt and z-offset error budget (§4.3) |
| `gold_only_calib.py` | the (L, a) degeneracy proof (§4.0a) |
| `ceo2_ambient_transfer.py` | one- vs two-line transfer from an off-position CeO2 |
| `marker_screen.py` | marker collision screen and sensitivity (§3.2a) |
| `ta_hydride_budget.py` | hydrogen inventory in the gasket hole (§4.0c) |

On `amber`, in `~/xpcs_contrast_check/`: qmap builders, the three comparison
`boost_corr` runs, and the fake-contrast diagnostic of §5.2.

All pope202609 analysis was read-only — no motor, attenuator, shutter or PV was
touched, and no ophyd session was started (session startup fires an Eiger frame).
