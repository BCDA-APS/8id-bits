"""Turn a rigaku_zdt_vs_ftf.py sweep into the table in the summary document.

    python analyse.py A0167 A0168 A0169 A0170 ...

Takes measurement numbers in (ZDT, FTF) pairs, or edit CLEAN below and run with
no arguments to reproduce the published table.

Size comes from the data files only -- the six .bin.NNN or .h5.NNN -- so the
~90 KB _metadata.hdf does not distort the sparse points, where it would
otherwise be several percent of the total.

Time comes from the gaps between directory mtimes, not from the log, because the
log only has one-second resolution and a repeat here is 2.7 s. Repeat 1 is
excluded: the interval before it includes one-off protocol setup (arming, file
templates), which is paid once rather than per repeat, so reps 2..N are the
steady-state acquire -> write -> metadata -> re-arm cycle.
"""

import glob
import os
import statistics
import sys

BASE = "/gdata/dm/8ID/8IDE/2026-3/comm202609/data"
FRAMES = 10000
# 131072 bytes per frame per module x 6 modules: what the same measurement would
# take stored densely at 2 bits per pixel, with no compression at all.
DENSE_BYTES = 131072 * 6 * FRAMES

# The published sweep. threshold keV -> {mode: measurement number}
CLEAN = {
    4: {"ZDT": "A0167", "FTF": "A0168"},
    5: {"ZDT": "A0169", "FTF": "A0170"},
    6: {"ZDT": "A0163", "FTF": "A0164"},
    7: {"ZDT": "A0165", "FTF": "A0166"},
}


def stats(measurement):
    """Per-repeat size and interval for one measurement number, or None."""
    dirs = sorted(glob.glob(f"{BASE}/{measurement}_*_r000*"))
    if not dirs:
        return None

    records = []
    for directory in dirs:
        data_files = [
            f
            for f in os.listdir(directory)
            if ".bin." in f or (".h5." in f and f.rsplit(".", 1)[-1].isdigit())
        ]
        records.append(
            (
                int(directory[-5:]),
                os.path.getmtime(directory),
                len(data_files),
                sum(os.path.getsize(os.path.join(directory, f)) for f in data_files),
            )
        )
    records.sort()

    sizes = [r[3] for r in records if r[2]]
    gaps = [records[i][1] - records[i - 1][1] for i in range(1, len(records))]
    return {
        "measurement": measurement,
        "repeats": len(records),
        "with_data": len(sizes),
        "mean": statistics.mean(sizes) if sizes else 0.0,
        "sd": statistics.stdev(sizes) if len(sizes) > 1 else 0.0,
        "total": sum(sizes),
        "gap": statistics.mean(gaps) if gaps else 0.0,
        "gap_sd": statistics.stdev(gaps) if len(gaps) > 1 else 0.0,
    }


def main():
    if len(sys.argv) > 1:
        names = sys.argv[1:]
        plan = {i: {"ZDT": names[i], "FTF": names[i + 1]} for i in range(0, len(names) - 1, 2)}
    else:
        plan = CLEAN

    print(
        f"{'thr':>4} {'mode':>4} {'meas':>6} {'reps':>5} {'mean/rep':>13} "
        f"{'sd':>8} {'vs dense':>10} {'s/rep':>7} {'sd':>5} {'total':>10}"
    )
    for threshold in sorted(plan):
        row = {}
        for mode in ("ZDT", "FTF"):
            s = stats(plan[threshold][mode])
            row[mode] = s
            if not s:
                print(f"{threshold:>4} {mode:>4}   -- not on disk --")
                continue
            # Near 1 the reciprocal form rounds away the whole result ("1/1"),
            # so only switch to it once the file is comfortably below dense.
            ratio = s["mean"] / DENSE_BYTES
            versus = f"{ratio:.2f}x" if ratio > 0.5 else f"1/{1 / ratio:.0f}"
            print(
                f"{threshold:>4} {mode:>4} {s['measurement']:>6} {s['with_data']:>5} "
                f"{s['mean'] / 1e6:>10.2f} MB {s['sd'] / 1e6:>7.3f} {versus:>10} "
                f"{s['gap']:>7.2f} {s['gap_sd']:>5.2f} {s['total'] / 1e9:>7.2f} GB"
            )
        if row.get("ZDT") and row.get("FTF"):
            print(
                f"{'':>4} {'-->':>4} FTF/ZDT size = {row['FTF']['mean'] / row['ZDT']['mean']:.2f}x   "
                f"FTF/ZDT time = {row['FTF']['gap'] / row['ZDT']['gap']:.2f}x"
            )


if __name__ == "__main__":
    main()
