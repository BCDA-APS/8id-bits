"""Rigaku 3M: ZDT 2-bit vs Fast Transfer 2-bit, swept over the lower threshold.

Results and interpretation: docs/reference/rigaku-zdt-vs-fast-transfer.md

Run it on a beamline host, from src/, inside the 8id_bits environment::

    python ../scripts/benchmarks/rigaku_zdt_vs_ftf.py 4 5 6 7

Each threshold runs mi_bench.yaml, which is two protocols (ZDT then Fast
Transfer) of ten repeats each. The threshold is written into the sample name via
si_thr<N>.yaml, so it ends up in every folder name.

Two details in here are not incidental -- both come from ways the first attempt
at this measurement produced wrong data, and both are cheap to keep:

1. A watchdog thread samples LowerThreshold_RBV every 0.5 s and prints every
   distinct value it saw during the block. A clean block prints exactly one
   value ({4.0: 331}). More than one means something else is writing the PV
   underneath you -- another session, or a previous run of this script that did
   not actually die -- and the data for that block is worthless, because the
   file name will say one threshold and the contents will be another.

2. The interrupt is re-raised. det_acq_series() converts Ctrl+C into a
   RuntimeError, so `except Exception` swallows it and the sweep marches on to
   the next threshold as if nothing happened. Stopping the sweep has to be
   explicit.

Analysis is disabled (analysis_machine = "none"): DM cannot read Fast Transfer
output at all, so submitting jobs would only produce failures to wade through.
"""

import json
import sys
import threading
import time
from pathlib import Path

# isort: off
# ORDER MATTERS -- do not let the formatter sort this block. startup_ophyd is
# what loads the devices, and the plan modules resolve device names at import
# time, so importing master_plan first would bind None for every device it
# looks up. Same rule as the startup scripts (see CLAUDE.md).
from id8_common.startup_ophyd import *  # noqa: F401,F403

from id8_common.registry import oregistry
from id8_common.expt_config import expt
from id8_common.plans.acquire import master_plan

# isort: on

HERE = Path(__file__).parent
RESTORE_TO = 6.0  # the threshold normal operation expects to find

thresholds = [float(a) for a in sys.argv[1:]] or [4.0, 5.0, 6.0, 7.0]

expt._static["analysis_machine"] = "none"
cam = oregistry["rigaku3M"].cam

seen = {}  # threshold -> {readback value: number of samples}
stop = threading.Event()
current = {"thr": None}


def watchdog():
    """Record every distinct LowerThreshold readback seen during each block."""
    while not stop.is_set():
        try:
            value = round(float(cam.lower_threshold.get()), 2)
            seen.setdefault(current["thr"], {}).setdefault(value, 0)
            seen[current["thr"]][value] += 1
        except Exception:
            pass  # a dropped read must not kill the sweep
        stop.wait(0.5)


threading.Thread(target=watchdog, daemon=True).start()

print(f"### threshold at launch = {cam.lower_threshold.get()} keV", flush=True)
wall = []
try:
    for thr in thresholds:
        cam.lower_threshold.put(thr)
        time.sleep(5)  # the detector needs a moment to apply it
        current["thr"] = thr
        got = cam.lower_threshold.get()
        print(f"\n### THRESHOLD {thr} keV (readback {got})", flush=True)

        started = time.time()
        try:
            master_plan.run_measurement_info(
                measurement_info_file=str(HERE / "mi_bench.yaml"),
                sample_info_file=str(HERE / f"si_thr{int(thr)}.yaml"),
            )
        except (KeyboardInterrupt, RuntimeError) as exc:
            if isinstance(exc, KeyboardInterrupt) or "stopped by user" in str(exc):
                print(f"### INTERRUPTED at {thr} keV -- stopping the sweep", flush=True)
                raise
            print(f"### THRESHOLD {thr} FAILED: {exc}", flush=True)

        wall.append(
            {
                "threshold": thr,
                "readback": got,
                "wall_s": round(time.time() - started, 1),
                "rbv_histogram": seen.get(thr, {}),
            }
        )
        print(f"### THRESHOLD {thr} keV took {wall[-1]['wall_s']} s", flush=True)
        print(f"### RBV seen during block: {seen.get(thr, {})}", flush=True)
finally:
    stop.set()
    cam.lower_threshold.put(RESTORE_TO)
    time.sleep(3)
    print(f"\n### restored lower_threshold = {cam.lower_threshold.get()} keV", flush=True)
    json.dump(wall, open("/tmp/bench_wall.json", "w"), indent=1)

print("### BENCH DONE", flush=True)
