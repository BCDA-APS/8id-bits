"""Extract all runs from the 8ide `raw` tiled catalog collected between SINCE and UNTIL

Writes, per run: a CSV of the primary stream table and a JSON of the run
metadata (start + stop docs). Also writes an index.csv summary.
"""

import datetime
import json
import os

import pandas as pd
from tiled.client import from_profile

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tiled_export_jun25-29")
SINCE = datetime.datetime(2026, 6, 25).timestamp()
UNTIL = datetime.datetime(2026, 6, 30).timestamp()  

def json_default(o):
    try:
        import numpy as np

        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
    except Exception:
        pass
    return str(o)

def read_primary(node):
    """Return the primary stream as a DataFrame, handling layout variants."""
    prim = node["primary"]
    for path in (["data"], ["internal"], ["internal", "events"]):
        cur = prim
        try:
            for p in path:
                cur = cur[p]
            return cur.read()
        except Exception:
            continue
    # last resort
    return prim.read()

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    cat = from_profile("8ide")["raw"]

    runs = []
    for uid, n in cat.items():
        st = n.metadata.get("start") or {}
        t = st.get("time")
        if t and SINCE <= t < UNTIL:
            runs.append((t, uid, n))
    runs.sort()
    print(f"Found {len(runs)} runs in window; writing to {OUTDIR}")

    index = []
    for t, uid, n in runs:
        st = n.metadata.get("start") or {}
        sp = n.metadata.get("stop") or {}
        sid = st.get("scan_id")
        stamp = datetime.datetime.fromtimestamp(t)
        base = f"scan_{sid:04d}_{stamp:%Y%m%d_%H%M%S}_{uid[:8]}"

        # metadata json
        with open(os.path.join(OUTDIR, base + ".meta.json"), "w") as f:
            json.dump({"uid": uid, "start": st, "stop": sp}, f, indent=2, default=json_default)

        # primary table csv
        nrows = 0
        try:
            df = read_primary(n)
            if isinstance(df, pd.DataFrame):
                df.to_csv(os.path.join(OUTDIR, base + ".csv"))
                nrows = len(df)
            else:
                df = df.to_dataframe() if hasattr(df, "to_dataframe") else pd.DataFrame(df)
                df.to_csv(os.path.join(OUTDIR, base + ".csv"))
                nrows = len(df)
        except Exception as e:
            print(f"  !! {base}: primary read failed: {e!r}")

        index.append(
            {
                "scan_id": sid,
                "uid": uid,
                "time": stamp.isoformat(),
                "plan_name": st.get("plan_name"),
                "motors": ",".join(st.get("motors", []) or []),
                "detectors": ",".join(st.get("detectors", []) or []),
                "num_points": st.get("num_points"),
                "exit_status": sp.get("exit_status"),
                "nrows": nrows,
                "csv": base + ".csv",
            }
        )
        print(f"  {base}  ({nrows} rows)")

    pd.DataFrame(index).to_csv(os.path.join(OUTDIR, "index.csv"), index=False)
    print(f"\nDone. {len(index)} runs -> {OUTDIR}/index.csv")


if __name__ == "__main__":
    main()
