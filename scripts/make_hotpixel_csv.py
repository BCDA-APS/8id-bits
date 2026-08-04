#!/usr/bin/env python
"""Write the flagged Rigaku 3M hot-pixel list to CSV.

Detection and per-pixel g2 leverage both come from the *pre-fix* results file, so
the numbers describe the problem as it was; the last column re-reads the patched
qmap to confirm each pixel is now blanked.
"""

import csv
import glob
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from plot_unnormalized_g2 import find_hot_pixels, hot_pixel_table  # noqa: E402

BASE = "/home/8-id-e/2026-2/pope202607/analysis/Multitau/A0096_GlassyCarbon_a0002_f001000_r%05d_G2_results.hdf"
QMAP_GLOB = "/home/8-id-e/2026-2/pope202607/data/*hotpixfix*.hdf"
OUTPUTS = [
    "/home/8-id-e/2026-2/pope202607/data/rigaku3m_hotpixels_A0096_r00001.csv",
    "/home/beams10/8IDIUSER/A0096_rigaku3M_hot_pixels.csv",
]

HEADER = [
    "row", "col", "dynamic_roi", "static_roi", "q_inv_angstrom",
    "cts_per_frame_r1", "cts_per_frame_r2", "cts_per_frame_r3",
    "cts_per_frame_r4", "cts_per_frame_r5",
    "local_7x7_median", "ratio_to_local",
    "IP_tau768", "IF_tau768", "G2_tau768",
    "g2_static_bin", "g2_static_bin_without_pixel",
    "delta_g2_ROI_tau1", "delta_g2_ROI_tau768", "blanked_in_patched_qmap",
]


def main():
    f = h5py.File(BASE % 1, "r")
    saxs = f["xpcs/temporal_mean/scattering_2d"][0].astype(np.float64)
    mask = f["xpcs/qmap/mask"][()]
    hot, local, ratio = find_hot_pixels(saxs, mask)
    dset = f["xpcs/multitau/unnormalized_G2"]
    qlist = f["xpcs/qmap/dynamic_v_list_dim0"][()]

    rows_t0, _, _ = hot_pixel_table(f, *(dset[0, k].astype(np.float64) for k in range(3)), hot)
    rows_t33, _, _ = hot_pixel_table(f, *(dset[33, k].astype(np.float64) for k in range(3)), hot)
    lev0 = {(r["y"], r["x"]): r["lev"] for r in rows_t0}
    last = {(r["y"], r["x"]): r for r in rows_t33}

    pers = []
    for run in range(1, 6):
        with h5py.File(BASE % run, "r") as g:
            s = g["xpcs/temporal_mean/scattering_2d"][0]
            pers.append([float(s[y, x]) for y, x in hot])
    pers = np.array(pers)

    qmap = sorted(glob.glob(QMAP_GLOB))[0]
    with h5py.File(qmap, "r") as g:
        bl, mk = g["qmap/blemish"][()], g["qmap/mask"][()]
        dqn, sqn = g["qmap/dynamic_roi_map"][()], g["qmap/static_roi_map"][()]
    print(f"verified against {qmap}")

    recs = []
    for i, (y, x) in enumerate(hot):
        r = last[(y, x)]
        blanked = bl[y, x] == 0 and not mk[y, x] and dqn[y, x] == 0 and sqn[y, x] == 0
        recs.append([
            y, x, r["dq"], r["sq"], f"{qlist[r['dq'] - 1]:.5f}",
            *[f"{v:.4f}" for v in pers[:, i]],
            f"{local[y, x]:.5f}", f"{ratio[y, x]:.0f}",
            f"{r['IP']:.4f}", f"{r['IF']:.4f}", f"{r['G2']:.2f}",
            f"{r['g2s']:.4f}", f"{r['g2s_wo']:.4f}",
            f"{lev0[(y, x)]:+.4f}", f"{r['lev']:+.4f}",
            "yes" if blanked else "NO",
        ])
    recs.sort(key=lambda v: -float(v[18]))

    for out in OUTPUTS:
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(HEADER)
            w.writerows(recs)
        print(f"wrote {out}")
    print(f"rows: {len(recs)} | all blanked in patched qmap: {all(r[-1] == 'yes' for r in recs)}")


if __name__ == "__main__":
    main()
