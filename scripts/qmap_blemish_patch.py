#!/usr/bin/env python
"""Add hot pixels to a qmap's blemish mask and drop them from every ROI map.

In this qmap layout `blemish == 1` means *good*, and

    mask == (blemish == 1) & (dynamic_roi_map > 0)

holds exactly, so blanking a pixel means all four of:
blemish -> 0, mask -> False, dynamic_roi_map -> 0, static_roi_map -> 0.
Editing blemish alone would leave the pixel in the ROI maps and change nothing.

Hot pixels are detected from a results file's run-averaged image
(`xpcs/temporal_mean/scattering_2d`) as pixels far above their local
neighbourhood; pass --extra to add coordinates by hand.
"""

import argparse
import shutil
from pathlib import Path

import h5py
import numpy as np
from scipy.ndimage import median_filter

Q = "qmap/"


def detect(results_hdf, ratio_thresh, floor, box):
    with h5py.File(results_hdf, "r") as f:
        saxs = f["xpcs/temporal_mean/scattering_2d"][0].astype(np.float64)
        mask = f["xpcs/qmap/mask"][()]
    local = median_filter(saxs, size=box)
    hot = mask & (saxs > floor) & (saxs / np.maximum(local, 1e-9) > ratio_thresh)
    return np.argwhere(hot), saxs, local


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("qmap", type=Path)
    ap.add_argument("--from-results", type=Path, required=True, help="results hdf to detect hot pixels from")
    ap.add_argument("--ratio", type=float, default=20.0)
    ap.add_argument("--floor", type=float, default=0.05, help="min counts/frame to consider")
    ap.add_argument("--box", type=int, default=7)
    ap.add_argument("--extra", type=int, nargs="+", default=[], metavar="ROW COL",
                    help="extra pixels to blank, as flat row col row col ...")
    ap.add_argument("--backup", action="store_true", help="write <qmap>.bak before editing")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    hot, saxs, local = detect(args.from_results, args.ratio, args.floor, args.box)
    if args.extra:
        if len(args.extra) % 2:
            ap.error("--extra needs an even number of values (row col row col ...)")
        extra = np.array(args.extra, dtype=np.int64).reshape(-1, 2)
        hot = np.unique(np.vstack([hot, extra]), axis=0) if len(hot) else extra

    with h5py.File(args.qmap, "r") as f:
        bl = f[Q + "blemish"][()]
        mk = f[Q + "mask"][()]
        dq = f[Q + "dynamic_roi_map"][()]
        sq = f[Q + "static_roi_map"][()]

    assert np.array_equal(mk, (bl == 1) & (dq > 0)), "unexpected blemish/mask convention"
    y, x = hot[:, 0], hot[:, 1]
    new = int(np.count_nonzero(bl[y, x] == 1))

    print(f"qmap        : {args.qmap}")
    print(f"detected    : {len(hot)} hot pixels ({new} not already blemished)")
    print(f"blemish==0  : {int((bl == 0).sum())} -> {int((bl == 0).sum()) + new}")
    print(f"mask True   : {int(mk.sum())} -> {int(mk.sum()) - int(mk[y, x].sum())}")
    print(f"dq>0 pixels : {int((dq > 0).sum())} -> {int((dq > 0).sum()) - int((dq[y, x] > 0).sum())}")
    dqs, cnts = np.unique(dq[y, x][dq[y, x] > 0], return_counts=True)
    print("per dynamic ROI: " + ", ".join(f"{a}:{b}" for a, b in zip(dqs, cnts)))

    # a static bin must not be emptied
    sq_hit = sq[y, x]
    for s in np.unique(sq_hit[sq_hit > 0]):
        before = int((sq == s).sum())
        after = before - int((sq_hit == s).sum())
        if after < 10:
            print(f"  WARNING static bin {s}: {before} -> {after} pixels")

    if args.dry_run:
        print("\ndry run - nothing written")
        return

    if args.backup:
        bak = args.qmap.with_suffix(args.qmap.suffix + ".bak")
        if not bak.exists():
            shutil.copy2(args.qmap, bak)
            print(f"backup      : {bak}")

    bl[y, x] = 0
    mk[y, x] = False
    dq[y, x] = 0
    sq[y, x] = 0
    with h5py.File(args.qmap, "r+") as f:
        f[Q + "blemish"][...] = bl
        f[Q + "mask"][...] = mk
        f[Q + "dynamic_roi_map"][...] = dq
        f[Q + "static_roi_map"][...] = sq

    with h5py.File(args.qmap, "r") as f:
        bl2, mk2 = f[Q + "blemish"][()], f[Q + "mask"][()]
        dq2, sq2 = f[Q + "dynamic_roi_map"][()], f[Q + "static_roi_map"][()]
    ok = (not bl2[y, x].any() and not mk2[y, x].any()
          and not dq2[y, x].any() and not sq2[y, x].any()
          and np.array_equal(mk2, (bl2 == 1) & (dq2 > 0)))
    print(f"\nwritten, verified: {ok}")


if __name__ == "__main__":
    main()
