# DAC centre-aperture search — 2026-09-09, 8-ID-E, pope202609

Scans are on /gdata under `data/bluesky/`; this folder holds logs, maps and notes.

## Result

The transmitting region is a **diagonal band**, not a round hole:
roughly 1.8 mm long × 0.7 mm wide, running from about
(huber_x −2.4, huber_y 41.5) to (−1.1, 42.6). Zero transmission everywhere outside it.

| quantity | value |
|---|---|
| peak transmission (Lambda2M ROI1) | huber_x **−2.167**, huber_y **41.742** |
| top-hat fit of pind cut at y = 41.742 | centre x **−2.006**, edges −2.354 / −1.658, width **0.696 mm**, R² 0.9986 |
| geometric centre of the band | x ≈ −1.75, y ≈ 42.16 |
| contrast | pind 14.5× (1.198 → 17.36); Lambda ROI1 0 → 749 258 |

The starting position (x 1.250, y 42.950) was **3.4 mm** from the aperture in x —
outside the ±2 mm box searched first, which is why the first three scans were blank.

## Scans

| file | detector | range / step | result |
|---|---|---|---|
| A0025 | pind | 2D ±2 mm, 0.25 mm, 289 pts | flat 1.1958–1.1975 — aperture not in inner box |
| A0026 | pind | x cut ±5 mm, 0.25 mm | flat, span 0.0045 |
| A0027 | pind | y cut ±5 mm, 0.25 mm | flat, span 0.0007 |
| A0028 | pind | 2D ±5 mm, 0.5 mm (aborted 251/441) | **found it** — 18.63 at x −2.25, y 41.45 |
| A0029 | pind | fine 0.21 mm (aborted 26/169) | superseded by Lambda |
| A0030 | Lambda2M ROI1 | 2D 2.5×2.5 mm, 0.21 mm, delta=nu=0, att 1e4 | full map of the band |
| A0031 | pind | x cut ±1.5 mm at y 41.742, delta=10, att=1 | verification, top-hat R² 0.9986 |

## Worth a second look

A 1–2 mm circular aperture should give a round spot. A 1.8 × 0.7 mm band tilted
~45° in the (x, y) plane suggests either the cell is tilted so the bore is seen
obliquely, or the two stage axes are not orthogonal in the plane normal to the beam.

## Notes for whoever reads the logs

* `8ideSoft:fastshutter:State` and `State_RBV` carry **opposite enum labels**.
  `showbeam()` leaves State="Close"/State_RBV="Open". `_blockbeam_verified()`
  reads state_rbv, so it is correct; anything reading `State` gets the inverse.
* `pre_align()` inserts the pind (`8idiSoft:FLIGHT:bo1:8` → 0 → reads "IN");
  `post_align()` retracts it. Inserting it alone raises current3 by ~0.28
  (dark offset), which is easy to mistake for beam.

## Final position (set 2026-09-09 20:06)

Stage moved to the **geometric centre of the band**, chosen over the peak:

    huber.x = -1.7801
    huber.y = 42.1600
    huber.delta = 10.0   huber.nu = 0.0

Confirmed in place with the pind: **15.704 shutter-open vs 1.198 closed**, i.e.
14.51 above baseline against 16.17 at the peak — **90% of peak flux**, not the
77% the Lambda ROI1 map predicted. ROI1 is only 100x10 px and clips part of the
beam; the pind integrates all of it, so trust the pind for relative flux and the
Lambda for position.

Candidates considered:

| | huber.x | huber.y | flux |
|---|---|---|---|
| peak pixel (Lambda ROI1) | -2.1667 | 41.7416 | 100% |
| bbox midpoint of >50% pixels | -1.7500 | 42.1583 | |
| intensity-weighted centroid | -1.8025 | 42.1599 | |
| **chosen: mean of the two centre estimates** | **-1.78** | **42.16** | 90% (pind) |

Left with the beam blocked, attenuation 1, pind IN, Lambda out of the beam (delta=10).
