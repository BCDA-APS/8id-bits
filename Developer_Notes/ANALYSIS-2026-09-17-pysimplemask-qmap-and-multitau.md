# Analysis log — 2026-09-17: pysimplemask qmaps + boost_corr multitau on G0256

**Host:** `amber` (needs `/gdata`)
**Tools:** `pysimplemask` 0.4.0.post2, `boost_corr` 0.1.5.dev132+gc3c61567a
**Dataset:** `G0256_HEA-Dec9GPa_a0001_f003000_lambda2M_r00001` — Lambda2M,
att 1, 3000 frames × 1 s, HEA sample decompressed to ~9 GPa
**Outcome:** β = 0.0108 ± 0.0001, τ_c = 479 ± 24 s, γ = 1.52

> **On the CLI entry point.** The documented entry point is
> `launch_simplemask_dev build` (on PATH at `/home/beams/8IDIUSER/bin/`), which
> exposes two options the older `pysimplemask-build-qmap` lacks:
> `--metadata-fname` and `--use-groupindex-for-subpartition`. The qmap below was
> first built with the older command and then rebuilt with
> `launch_simplemask_dev build --metadata-fname <meta>`: **both produced an
> identical qmap** — 388,419 pixels, q 3.2642–3.3116 Å⁻¹, virtual centre 19199.0 —
> so the result is unaffected. Neither extra flag was needed here; auto-discovery
> already found the right metadata file.

All work was read-only with respect to the beamline: file reads plus compute
jobs. No motor, attenuator, shutter or PV was touched, and no ophyd session was
started.

---

## 1. Locating the tools

`launch_simplemask_dev` **is** on PATH; the underlying package lives in Miaoqi
Chu's XPCS environment:

```
/home/beams/8IDIUSER/bin/launch_simplemask_dev         # documented entry point
/home/beams/8IDIUSER/bin/launch_simplemask_dev_combine
/home/beams/MQICHU/.conda/envs/{p,d}2504_xpcs/bin/pysimplemask[-build-qmap]
/home/beams/8IDIUSER/bin/boost_corr_bin                # wrapper, same env
```

Raw data root is `/gdata/dm/8ID/8IDE/<cycle>/<proposal>/` — note **`8ID/8IDE`**,
not `8IDE` or `8IDI` directly.

---

## 2. qmap settings used

### 2.1 Primary — tight ROI on the sample peak

```bash
launch_simplemask_dev build "$RAW" \
  --metadata-fname "$META" \
  --no-find-center \
  --beamstop-diameter 0 \
  --param-constraint q:AND:3.2632:3.3126 \
  --threshold-high 20 \
  --mode q-phi \
  --dq-num 1 --sq-num 25 --dp-num 1 --sp-num 1 \
  --style linear \
  --output-qmap psm_tight.hdf --output-mask "" --report psm_tight.pdf
```

| flag | value | why |
|---|---|---|
| `--no-find-center` | set | This is a **wide-angle** frame at `huber_delta` = 24.48°. There is no direct beam on the detector, so `goto_max`/`find_center` has nothing to lock onto and would corrupt the geometry. Use the metadata centre. |
| `--beamstop-diameter` | `0` | No beamstop in the field at wide angle; the default 30 px would punch a hole in real data. |
| `--param-constraint` | `q:AND:3.2632:3.3126` | **One FWHM centred on the peak.** Peak measured at q = 3.2879 Å⁻¹ with FWHM 0.0495 Å⁻¹ from the scattering pattern itself (not the rocking curve — see §5). |
| `--threshold-high` | `20` | Hot-pixel cut in raw counts. Typical pixel here is ~0.5 counts/frame, so 20 is ~40× the norm and removes only genuine outliers. |
| `--dq-num` / `--sq-num` | `1` / `25` | One dynamic ROI spanning the whole peak; 25 static ROIs so each is near iso-intensity. This is the key setting — see §5. |
| `--dp-num` / `--sp-num` | `1` / `1` | No azimuthal subdivision; the ROI already spans only ~16° of arc. |
| `--style` | `linear` | Uniform q bins across a narrow range; logarithmic would be pointless over 0.05 Å⁻¹. |

Result: **388,419 pixels kept** (86.25% of the detector masked out),
`dqmap/sqmap consistency check: True`.

### 2.2 Control — deliberately sloppy

Identical except no `--param-constraint`, and `--sq-num 1`: the **whole detector
as a single static bin**. Built only to quantify the artifact in §4. Kept in
scratch, deliberately *not* written to the analysis folder.

### 2.3 Blemish handling

pysimplemask found the site blemish map on its own:

```
/home/beams/8IDIUSER/Documents/areaDetectorBlemish/8idLambda2m/latest_blemish.tif
```

This is better than the ad-hoc 8σ MAD filter used in the earlier hand-built
attempt, and is one reason the two runs differ slightly (§3).

---

## 3. Geometry check — pysimplemask handles `huber_delta` correctly

Worth verifying explicitly, because a **stale qmap that ignores delta** was a real
source of error during pope202609 (see `EVAL-2026-09-13-…` §2.3).

pysimplemask emits `beam_center_y = 19199.03`, a *virtual* centre far off the
1813-row detector. Decoding it:

```
2θ(row 889) = atan((19199.03 − 889) × 55 µm / 2.2 m) = 24.593°
independent calculation from metadata delta = 24.4798°:  24.62°
```

Agreement to 0.027°, about two pixels. **The delta is being folded into the
geometry properly** — the qmap is not stale.

---

## 4. Results

### 4.1 Contrast

`boost_corr` reports g2 per *static* bin, averaging pixels within each bin before
dividing, so g2(∞) = 1 + Var_p(Ī)/⟨Ī⟩². The static floor was therefore measured
independently from the time-averaged image over pysimplemask's own 25 static
bins, and held fixed in the fit:

```
static floor (25 bins)  = 1.00860
g2(τ = 1 s)             = 1.0195
beta                    = 0.0108 ± 0.0001
tau_c                   = 479 ± 24 s
gamma                   = 1.52
```

Fit range τ ≤ 1000 s; longer delays fall below the static floor, indicating
residual finite-length bias, and were excluded.

### 4.2 Cross-check against the hand-built qmap

The same dataset was analysed on 2026-09-16 with a qmap constructed by hand in
h5py. Two independently built qmaps, different masking, same physics:

| | hand-built | pysimplemask | agreement |
|---|---|---|---|
| static floor | 1.0091 | 1.00860 | — |
| g2(1 s) | 1.0191 | 1.0195 | 0.04% |
| **β** | **0.0096** | **0.0108** | **12%** |
| τ_c | 500 s | 479 s | 4% |
| γ | 1.62 | 1.52 | 6% |

The 12% spread on β exceeds the 1% statistical error, so it is systematic —
different pixel sets (site blemish map vs MAD filter) and slightly different ROI
edges. It sits well inside the ±0.005 systematic band quoted previously.

**Take β ≈ 0.010 for this dataset**, systematic-limited, not statistics-limited.

### 4.3 The fake-contrast artifact, re-measured

| qmap | g2(1 s) | naive "contrast" g2(0) − 1 | overstatement |
|---|---|---|---|
| tight, 25 static bins | 1.0195 | 0.0195 | 1.8× |
| **whole detector, 1 static bin** | **1.4649** | **0.4649** | **43×** |

The whole-detector control is far worse than the whole-*module* control tried
earlier (8.6×) because it pools three modules spanning q = 3.12–3.45 Å⁻¹ and all
azimuths into one bin. **Quoting g2(0) − 1 from an unconstrained ROI would have
reported 0.46 contrast where the truth is 0.010** — a factor of 43, essentially
all of it spatial intensity structure rather than dynamics.

---

## 5. Two settings that matter more than the rest

**Size the ROI from the scattering pattern, not the rocking curve.** The delta
scan S00234 gives FWHM 0.473°, but that is convolved with the stats ROI
acceptance. Resolving the peak directly on the detector gives **0.377°**
(0.0495 Å⁻¹). Using the rocking-curve width would have over-widened the ROI by
~25% and inflated the floor.

**Static bins are where the artifact lives.** `--dq-num` controls the output
binning; `--sq-num` controls what gets averaged before the division that forms
g2. Fine static bins suppress the fake baseline even when the dynamic ROI is
wide. 25 bins over one FWHM puts each bin at ~0.002 Å⁻¹, near iso-intensity.

---

## 6. Output files

Written to the analysis folder, both readable in pyXpcsViewer:

```
/gdata/dm/8ID/8IDE/2026-3/pope202609/analysis/Multitau/
    G0256_HEA-Dec9GPa_a0001_f003000_lambda2M_r00001_psm_tight_results.hdf   (1.3 MB)
    qmap_G0256_psm_tight.hdf                                                (227 kB)
```

The earlier hand-built pair is alongside, suffixed `_TightPeakROI`, for
comparison. Scratch (qmaps, the sloppy control, PDF reports, diagnostics) is in
`~/xpcs_contrast_check/` on `amber`.

---

## 7. Not done / open

- Only **G0256** was analysed. The rest of the G0252–G0256 attenuation series
  (att 10, 5, 2, 1) is untouched; running it would show whether β is
  flux-independent, which is the standard check that the contrast is real and
  not a detector artifact. Worth doing.
- The **gold peak** was not analysed for XPCS — this was the sample peak only.
- `launch_simplemask_dev_combine` was not needed and not exercised.
- `--use-groupindex-for-subpartition` was not needed here (one q region), but it
  is the right tool if several disjoint peaks are ever partitioned at once.
- Detector artifact masking: see
  `PROPOSAL-2026-09-17-automatic-artifact-masking.md`. Note it identifies **1479
  dead Lambda2M pixels missing from the blemish file**, which affects the qmaps
  built here (they are inside the masked-out region for this particular tight ROI,
  so this result is unaffected, but wider ROIs would include them).
- β ≈ 0.010 remains coherence-limited rather than sample-limited: at 2θ = 24.65°
  the path-length spread across the sample far exceeds the longitudinal coherence
  length. See `Rigaku500k_Experimental_Design.md` §5.4.
