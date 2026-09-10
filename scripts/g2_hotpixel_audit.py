#!/usr/bin/env python
"""Find the pixels responsible for inflated multitau g2 in a given dynamic ROI.

Reproduces the boost_corr reduction exactly:

    per static bin s:  G2_s = <G2>_s,  IP_s = <IP>_s,  IF_s = <IF>_s
                       g2_s = G2_s / (IP_s * IF_s)
    per dynamic bin d: g2_d = mean_{s in d} g2_s          (unweighted!)

so a single bright pixel enters the numerator of its static bin as I^2 but the
denominator only as I -> it inflates g2_s, and because dynamic bins average
static-bin *ratios* without weighting, one small static bin can drag the whole
ROI up.

Leverage of pixel j = g2_d(all pixels) - g2_d(pixel j dropped from its bin).
"""

import argparse
from pathlib import Path

import h5py
import numpy as np


def load(path, tau_index):
    with h5py.File(path, "r") as f:
        dq = f["xpcs/qmap/dynamic_roi_map"][()]
        sq = f["xpcs/qmap/static_roi_map"][()]
        G2, IP, IF = (f["xpcs/multitau/unnormalized_G2"][tau_index, k].astype(np.float64) for k in range(3))
        ref = f["xpcs/multitau/normalized_g2"][tau_index]
        qlist = f["xpcs/qmap/dynamic_v_list_dim0"][()]
        delays = f["xpcs/multitau/delay_list"][()]
        saxs = f["xpcs/temporal_mean/scattering_2d"][0].astype(np.float64)
    return dq, sq, G2, IP, IF, ref, qlist, delays, saxs


def bin_sums(sq_flat, nsq, *arrays):
    return np.stack([np.bincount(sq_flat, weights=a, minlength=nsq + 1) for a in arrays])


def g2_static(sums, counts):
    with np.errstate(invalid="ignore", divide="ignore"):
        avg = sums / counts
    ipif = avg[1] * avg[2]
    ipif = np.where(ipif == 0, 1e8, ipif)
    return avg[0] / ipif, avg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hdf", type=Path)
    ap.add_argument("--tau", type=int, default=-1, help="index along axis 0 (default: last)")
    ap.add_argument("--rois", type=int, nargs="+", default=[7, 8], help="dynamic ROI labels (1-based, as in dynamic_roi_map)")
    ap.add_argument("--top", type=int, default=15, help="pixels to list per ROI")
    args = ap.parse_args()

    dq2d, sq2d, G2, IP, IF, ref, qlist, delays, saxs = load(args.hdf, args.tau)
    ntau = h5py.File(args.hdf, "r")["xpcs/multitau/unnormalized_G2"].shape[0]
    ti = args.tau % ntau
    ny, nx = dq2d.shape
    dq, sq = dq2d.ravel(), sq2d.ravel()
    nsq = int(sq.max())
    counts = np.bincount(sq, minlength=nsq + 1).astype(np.float64)

    sums = bin_sums(sq, nsq, G2.ravel(), IP.ravel(), IF.ravel())
    g2s, avg = g2_static(sums, counts)

    print(f"file  : {args.hdf.name}")
    print(f"tau   : index {ti}  (delay = {delays[ti]} frames)")
    print(f"static bins {nsq}, dynamic bins {int(dq.max())}\n")

    for roi in args.rois:
        sel = dq == roi
        sbins = np.unique(sq[sel])
        sbins = sbins[sbins > 0]
        g2_roi = np.nanmean(g2s[sbins])
        print("=" * 96)
        print(f"dynamic ROI {roi}  (normalized_g2 column {roi - 1}, q = {qlist[roi - 1]:.5f} A^-1)")
        print(f"  pixels {sel.sum()},  static sub-bins {len(sbins)},  g2 = {g2_roi:.4f}  (file: {ref[roi - 1]:.4f})")

        # --- which static sub-bins carry the excess -------------------------
        order = np.argsort(g2s[sbins])[::-1]
        med = np.nanmedian(g2s[sbins])
        print(f"  median static-bin g2 = {med:.4f}")
        print(f"  worst static sub-bins (each one contributes g2_s/{len(sbins)} to the ROI mean):")
        print(f"    {'sbin':>6} {'npix':>6} {'g2_s':>10} {'<IP>':>10} {'maxIP':>10} {'maxIP/<IP>':>11} {'pull':>8}")
        for k in order[:8]:
            s = sbins[k]
            m = sq == s
            ipv = IP.ravel()[m]
            pull = (g2s[s] - med) / len(sbins)
            print(f"    {s:6d} {int(counts[s]):6d} {g2s[s]:10.3f} {avg[1][s]:10.4f} "
                  f"{ipv.max():10.4f} {ipv.max() / max(avg[1][s], 1e-12):11.1f} {pull:+8.4f}")

        # --- exact single-pixel leverage ------------------------------------
        # candidates: the brightest pixels of the worst bins
        cand = []
        for k in order[: max(8, len(sbins) // 4)]:
            s = sbins[k]
            idx = np.flatnonzero(sq == s)
            bright = idx[np.argsort(IP.ravel()[idx])[::-1][:6]]
            cand.extend(bright.tolist())
        cand = np.unique(cand)

        g2f, ipf, iff = G2.ravel(), IP.ravel(), IF.ravel()
        rows = []
        for j in cand:
            s = sq[j]
            n = counts[s]
            if n <= 1:
                continue
            a = (sums[0][s] - g2f[j]) / (n - 1)
            b = (sums[1][s] - ipf[j]) / (n - 1)
            c = (sums[2][s] - iff[j]) / (n - 1)
            new = a / (b * c) if b * c != 0 else np.nan
            delta = (g2s[s] - new) / len(sbins)   # effect on the ROI mean
            rows.append((delta, j, s, new))
        rows.sort(reverse=True)

        print(f"\n  top {args.top} single pixels by leverage on ROI g2:")
        print(f"    {'row':>5} {'col':>5} {'sbin':>6} {'IP':>10} {'IF':>10} {'G2':>12} "
              f"{'<I>_run':>9} {'I/nbr':>7} {'g2_s->':>9} {'dg2_ROI':>9}")
        for delta, j, s, new in rows[: args.top]:
            y, x = divmod(int(j), nx)
            y0, y1, x0, x1 = max(0, y - 3), min(ny, y + 4), max(0, x - 3), min(nx, x + 4)
            patch = saxs[y0:y1, x0:x1].copy()
            patch[y - y0, x - x0] = np.nan
            nbr = np.nanmedian(patch)
            ratio = saxs[y, x] / nbr if nbr > 0 else np.inf
            print(f"    {y:5d} {x:5d} {s:6d} {ipf[j]:10.4f} {iff[j]:10.4f} {g2f[j]:12.2f} "
                  f"{saxs[y, x]:9.4f} {ratio:7.1f} {new:9.3f} {delta:+9.4f}")

        # --- cumulative effect of dropping the worst pixels -----------------
        print("\n  cumulative: drop the N highest-leverage pixels and recompute the ROI g2")
        for n_drop in (0, 1, 2, 5, 10, 20, 50):
            drop = np.array([r[1] for r in rows[:n_drop]], dtype=np.int64)
            s2 = sums.copy()
            c2 = counts.copy()
            for j in drop:
                s2[0][sq[j]] -= g2f[j]
                s2[1][sq[j]] -= ipf[j]
                s2[2][sq[j]] -= iff[j]
                c2[sq[j]] -= 1
            gg, _ = g2_static(s2, c2)
            print(f"    drop {n_drop:3d} -> g2 = {np.nanmean(gg[sbins]):.4f}")
        print()


if __name__ == "__main__":
    main()
