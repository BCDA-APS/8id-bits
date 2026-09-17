# Evaluation — 2026-09-13: DAC/WAXS geometry and protocol for 2027-1

Written so the thread survives losing the chat it came from. Covers the
pope202609 diamond-anvil-cell beamtime (Pope/Vohra EHEA, BCC→FCC under
pressure), what the analysis could and could not establish, why, and what to
change for 2027-1 (LDRD + Christian Gutt).

| § | Question | Answer |
|---|---|---|
| 1 | What did we actually measure? | A clean pressure ladder and one line that vanishes between 11.5 and 22 GPa |
| 2 | Why could we not index the pattern? | Three avoidable reasons — **plus a hard aperture limit; §2.2 corrected 2026-09-16** |
| 3 | How far should Lambda2M sit? | ≤1.5 m to get both gold lines in one frame; **1.0 m** for everything at once. Route: reduce the detector block — but feasibility unconfirmed, see §3.4 |
| 4 | Should we raise the energy? | **No.** Three multiplicative losses. Fix the aperture/geometry instead |
| 5 | Protocol for next time | §5 — the ambient reference is the single highest-value item |
| 6 | Tooling | `analysis/mesh_acq/`, reusable as-is |

Science context: Pope et al., *Sci. Rep.* **14**, 16472 (2024) — 3-D printed
EHEA Ni40Co20Fe10Cr10Al18W2, nanolamellar BCC (a=2.871 Å) + FCC (a=3.591 Å) at
ambient, irreversible BCC→FCC completing at 9 GPa, single-phase FCC a=3.456 Å at
30.6 GPa. They used 0.4246 Å (29.2 keV) at HPCAT. We are at 15.2 keV.

---

## 1. What was measured

### 1.1 Pressure ladder — from Au 111 arm scans, not from any qmap

| scan | pcd2 | Au 111 2θ | a_Au (Å) | P (GPa) |
|---|---|---|---|---|
| B0140 | 0 | 19.962 | 4.0756 | ~0 |
| S00212 | 370 | 20.326 | 4.0035 | 11.5 |
| S00237 | 450 | 20.612 | 3.9485 | 22.0 |
| S00240 | 505 | 20.816 | 3.9102 | 30.8 |

Ambient point gives a₀ = 4.0822 Å against the literature 4.0782 — **0.10%**.
That is the validation of the whole 2θ reconstruction.

**`pcd1` is NOT the membrane.** It reads 0.8 with setpoint 0.0 on every scan.
**`pcd2`** is the live one. This wasted real time; log it in the note-to-self
for next run.

### 1.2 The phase result

A line at **d = 2.006 Å** is present at 11.5 GPa in two independent scans
(S00217, F0183 — different positions, different ranges, agreeing to 0.005 Å)
and absent at 22 GPa in two more (S00219 at 8.4σ, S00234 at 5.4σ). It cannot
have moved out of range (would need 12% compression) and it is not textured
away (azimuthal uniformity 0.08, i.e. fine-grained powder).

**Something transformed between 11.5 and 22 GPa.** That stands.

As BCC 110 that line gives a = 2.8373 Å, −1.2% vs the paper's ambient 2.871 —
exactly right for a compressed BCC. But since *nothing else* on the pattern
indexes on the paper's cell, treat the BCC label as suggestive, not established.

### 1.3 What remains unexplained

Two survivors at 22 GPa, both real smooth powder rings, neither indexed:

| line | d at 11.5 GPa | d at 22 GPa | as FCC 111 | verdict |
|---|---|---|---|---|
| low-Q, 21.95° | 2.1389 | 2.1419 | a = 3.72 Å (expanded) | unidentified |
| high-Q, 24.97° | 1.9116 | 1.8866 | a = 3.27 Å (→105 GPa) | unidentified |

Ruled out for both: gold (three independent tests), stainless-steel gasket
(every SS reflection would have to *expand*), and each other as a 111/200 pair
(d ratio 1.1373 vs the required 1.1547, ~5σ).

---

## 2. Why the pattern could not be indexed — four avoidable causes

**2.1 No ambient reference at the sample position.** This is the root cause.
B0139 is the only ambient sample scan: 0.5 s/pt and taken 68 µm from where the
sample sat. With a proper ambient pattern you index at ambient — where the
paper tells you the phases — and then follow each line *by continuity* as you
compress. No ratios, no extra reflections needed. Two lines would have sufficed.

**2.2 ~~Scans stopped at 2θ = 27.2° when the cell allows 60°.~~ — CORRECTED
2026-09-16, this cause is withdrawn.**

> The original text read: *"The DAC half-angle is ±30°, giving d_min = 0.816 Å at
> 15.2 keV — six accessible reflections for a ~3.5 Å cell (111, 200, 220, 311,
> 222, 400). Every scan of the beamtime stayed below 27.2°. The indexing was
> solvable the whole time and nobody looked."*
>
> **This is wrong.** The ±30° opening does not give 60° of usable 2θ. Reaching 60°
> would require the incident beam to enter along the −30° edge of the cone,
> grazing the anvil with a long diamond path — not a practical geometry. With the
> beam on the cell axis the usable maximum is **2θ = 30°**, q ≤ 3.990 Å⁻¹.
>
> At 30° a 3.456 Å FCC cell offers only **two** reflections: 111 at 23.6° and 200
> at 27.3°. Scans reaching 27.2° were therefore not leaving reflections unmeasured.
>
> **The indexing failure was an energy/aperture limitation, not a scan-range
> oversight.** No longer scan and no extra detector can fix it; only higher energy
> or a wider-aperture cell can. This reopens §4 (which argues against raising the
> energy) — that trade must be revisited now that the scan-range escape route is
> closed. See `Rigaku500k_Experimental_Design.md` §3.1 and §3.5.

The remaining causes below (2.1, 2.3, 2.4) stand as written.

**2.3 Stale qmap.** The lambda2M qmap is built from beam centre and distance and
carries no delta. It is valid only at the delta it was generated at (≈20.28).
Every dataset taken at another delta has a systematically wrong q axis. G0222
(delta 20.58) was analysed for hours before this surfaced.

**2.4 Coverage-edge bias.** A peak sitting near the end of a delta scan has a
badly biased centroid — measured: truncating S00234 at S00217's edge moved the
fitted centre by **0.35°**, twenty times the shift being looked for. Two separate
wrong conclusions came from this. Keep any peak of interest ≥0.5° inside the
scan range.

---

## 3. Lambda2M distance

### 3.1 The numbers (1813 × 55 µm = 99.7 mm tall)

| L (m) | 2θ span/frame | mdeg/px | azimuth @20° | Au 111+200 together? |
|---|---|---|---|---|
| **2.2 (current)** | 2.23° | 1.43 | **6°** | no |
| 1.5 | 3.27° | 2.10 | 9° | **yes, marginally** |
| **1.0** | **4.84°** | 3.15 | **13.5°** | yes, with margin |

Maximum distance that still captures both gold lines: **1.52 m at 20 GPa**,
1.48 m at 40 GPa — the Au 111/200 separation only grows 3.12° → 3.30° from
ambient to 40 GPa, so ~1.5 m covers any pressure we would use.

### 3.2 Why 1.0 m rather than 1.5 m

At 1.5 m one frame holds the two gold lines *or* the two sample lines. At 1.0 m
one frame (2θ 20.54–25.38) holds **all four at once**: Au 111, low-Q, Au 200,
high-Q. That is the qualitative change:

- **Pressure and phase in the same exposure, same spot, same instant.** This run
  we measured pressure on gold at (−0.1943, 21.2728) and phase on sample at
  (−0.2007, 21.3250), 55 µm and sometimes hours apart, and had to *assume* they
  corresponded. Never checked.
- **Diffraction becomes a readout, not an experiment.** A phase check costs
  25–35 min now, so it happened a handful of times in two days. At 1 m it is one
  exposure.
- **Spatial scans are where this pays.** 441 mesh points × one frame each = 441
  four-line patterns instead of 441 one-line slices → a **phase map and a
  pressure map**, since Au 111 + Au 200 in every frame gives a_Au from the exact
  2/√3 ratio at every point, with no calibration. Same acquisition cost as the
  mesh we already ran (~2 h).
- **"Absent" starts to mean absent.** At 2.2 m only ~1.7% of each Debye ring
  lands on the detector. That is why Au 200 never appeared on Lambda despite
  being in range, and why several "the line is gone" claims needed walking back.

Long XPCS runs do **not** need this — a few Lambda frames at spaced delta during
a 3000 s acquisition is sufficient (at 1 m, delta 22/27/32 covers 2θ 19.5–34.0;
at 2.2 m the same needs ~7 frames).

### 3.3 Hardware routes for 2027-1 (from the exchange with Suresh)

Open question to Suresh: **the 8-ID-E door poster shows Lambda without the
cylindrical flight path — what is the shortest sample-to-detector distance in
that configuration?** Working assumption 1 m.

| option | gains | costs |
|---|---|---|
| Remove cylindrical pipe only | 1 m, easy, was designed for it | none known; stopgap only |
| Push Lambda upstream, terminate the 8-ID-I entrance flight path with the 10 mm SiN window | 1 m without rail or 3rd flight path; WAXS retained; the 1 m can be left in air | **loses nu motion** |
| **Reduce the detector block — i.e. bring the detector in on the arm (preferred *if feasible*, see §3.4)** | see below | at the short distance, no XPCS on Lambda; reversible only if the arm can be extended without the large block |

**Preferred route: reduce the detector block.** Reducing the block and moving the
detector toward the sample are the same physical change, not two options. Three
reasons it beats the alternatives, and they compound:

1. **delta and nu are both retained.** The block still exists, just shorter, so
   the arm keeps its full angular freedom — unlike the "push Lambda upstream"
   route, which sacrifices nu.
2. **It allows a 10 mm SiN window instead of 15 mm, which is the real
   de-risking.** Required aperture scales linearly with distance from the
   sample. Rigaku's SA-XPCS needs q ≤ 0.0292 Å⁻¹ = 2θ 0.217°, so:

   | window at | diameter actually needed | margin with 10 mm |
   |---|---|---|
   | **1.0 m** | **3.8 mm** | **2.6×** |
   | 2.0 m | 7.6 mm | 1.3× |
   | 2.2 m | 8.3 mm | 1.2× |

   A 10 mm window at 1 m passes q ≤ 0.0770 Å⁻¹, well beyond what Rigaku uses.
   15 mm is only forced if the window sits at ~4 m. **We do not know whether a
   15 mm SiN window is even feasible** — they are fragile at that diameter — so
   designing the geometry to need only 10 mm removes a vacuum/safety and
   schedule risk rather than merely being convenient.
3. **It need not foreclose long-distance XPCS on Lambda, and might extend it.**
   *If* a longer rail can be mounted off the reduced block — and *if* the lower
   flight path and some of the 8-ID-I entrance flight path can be removed —
   Lambda could reach 2.2 m or **beyond** when an experiment wants XPCS rather
   than WAXS. The loss of Lambda XPCS at the short position is a consequence of
   *distance*, not of angular freedom, so it is reversible by extension rather
   than being designed out.

So one change buys the close WAXS geometry, keeps both arm angles, shrinks the
window to a size known to work, and leaves the long-distance XPCS option open.
The other two options each give up one of those.

### 3.4 Open engineering questions — none of this is confirmed

The preferred route rests on two assumptions that **have not been checked with
engineering**, and the whole §3.3 recommendation collapses to the fallback
options if either fails:

| assumption | status |
|---|---|
| The detector block can be reduced / the detector brought in on the arm at all | **unknown** |
| A longer rail can be mounted off the reduced block to get Lambda back out to ≥2.2 m | **unknown** |
| Shortest sample-to-detector distance in the poster (no cylindrical flight path) configuration | **unknown** — working assumption 1 m |
| Whether a 15 mm SiN window is feasible at all | **unknown** — which is precisely why the 10 mm route is worth the design effort |

If the block cannot be reduced, the fallback is "push Lambda upstream and
terminate with the 10 mm window", accepting the loss of nu. If neither is
available in time for 2027-1, removing the cylindrical pipe alone still gets
Lambda to ~1 m and delivers most of the §3.2 benefit for the spatial scans.

Going closer than 1 m by removing both flight paths is too disruptive for user
operation and is not the route being proposed here.

---

## 4. Energy: do not raise it

Tempting (29.2 keV is what the paper used) and wrong for this instrument. Three
multiplicative losses going 15.2 → 27 keV:

| penalty | factor |
|---|---|
| Si 350 µm QE on Rigaku/Lambda (∝E⁻³) | 5.6× |
| 1st → 3rd harmonic (1st absorbed by mono) | ≥2× |
| speckle 90 → 50 µm vs 76 µm Rigaku pixel | 2.3× |

Eiger is CdTe and would be fine; Rigaku and Lambda are not. 15.2 keV is also
about the ceiling of the 1st harmonic, and the mono needs enough incident angle
that Si 311 does not crack under heat load. The observed ~10% contrast on Rigaku
at 27 keV (11 m, 10 µm focus) is exactly what the three factors predict.

**The aperture, not the energy, is the lever** — and at ±30° we already have six
reflections. We simply never scanned out to them.

---

## 5. Protocol for next time

1. **Ambient pattern at the sample position before the first compression.** Full
   counting time, correct position, same geometry used for everything after.
   Highest-value single item in this list. Index here, then track by continuity.
2. **Confirm which pressure channel is live** (`pcd2`, not `pcd1`) and record it.
3. **Pressure from a fixed-ROI arm scan on Au 111, every time.** With a fixed ROI
   the sample-to-detector distance cancels exactly — peak appears at
   `delta = 2θ − φ`, φ constant. The qmap route cannot compete: holding 2θ to
   0.02° (≈0.5–1 GPa) needs L known to 2.3 mm at 2.2 m, 1.1 mm at 1 m. The delta
   motor gives ~0.001° ≈ 0.03 GPa. **Note this gets harder, not easier, as the
   detector moves closer** — the close detector is a survey instrument, the arm
   is the metrology one.
4. **One wide survey per pressure point, out to the cell limit** (2θ → ~58°).
   At 2.2 m: delta 21→58, step 0.5°, 20 s = 75 pts ≈ 25 min. Also measures the
   true aperture (where signal dies) and gives far better Δd/d, since precision
   improves as cot θ — a line at 48° is ~4× better than one at 24°.
5. **Re-find the sample after every pressure change.** It moves. The G0202 pick
   was stale by the time of G0222.
6. **Keep scan geometry identical between pressure points** you intend to
   compare, and keep peaks ≥0.5° inside the scan range (§2.4).
7. **Regenerate the Lambda qmap whenever delta moves**, or record delta in it —
   or, with the §5.3 split, stop relying on the Lambda qmap altogether.
8. **Spatial mesh with full patterns** once Lambda is ≤1.5 m. Gives phase and
   local pressure at every point, which also settles whether the gold spot and
   the sample spot are at the same pressure.

---

## 6. Tooling

In `/home/8-id-i/<cycle>/<expt>/analysis/mesh_acq/` (not in the repo):

| script | does |
|---|---|
| `reduce_delta_scan.py` | saved delta-scan frames → absolute 2θ powder profile. Each frame spans 2.6° and the scan slides it, so a scan is a heavily over-determined 2θ map |
| `reduce_spatial.py` | huber_x/y scan → position vs 2θ map |
| `spottiness.py` | azimuthal texture vs 2θ — separates textured gold (metric ~5) from fine-grained sample (~0.1) with no lattice assumption |
| `show_frames.py` | renders frames at chosen 2θ — ring vs spots by eye |
| `mesh_pick.py` / `gold_pick.py` | rank mesh points for gold-free or gold-bearing spots; `infer_mesh` now tolerates ≤2% off-grid readbacks |
| `audit_devices.py` | builds the registry **without** `ad_setup`, i.e. hardware-safe device access mid-acquisition |

Raw detector frames use the LZ4 filter — read them with the `8ide_bits_test`
conda python (`hdf5plugin`), not the system python3. `/gdata` is visible from
amber and pearl only.

Two gotchas baked into the scripts, worth not rediscovering:

- the q axis is `static_v_list_dim0[static_index_mapping]`, **0-based** — a
  detector with module gaps reports fewer bins than the v_list holds (lambda2M:
  318 of 360), and reading the v_list in order shifts the profile by a module.
- masked pixels carry ~2²⁴ sentinels; mask first, then clip.
