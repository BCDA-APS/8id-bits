"""Stage B: the actual test scan.

    dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M, att_ratio=5)

Full dscan() behaviour (attenuator, pre_align, PIND, shutter) as requested.
Moves only huber.delta among the positioners.
"""
import time

print(">>> importing id8_common.startup ...", flush=True)
from id8_common.startup import *  # noqa: F401,F403

from apsbits.core.instrument_init import oregistry  # noqa: E402
from id8_common.plans.align.ophyd_scan import dscan_ophyd, default_spec_path  # noqa: E402

huber = oregistry["huber"]
lambda2M = oregistry["lambda2M"]
pv_registers = oregistry["pv_registers"]

start_pos = huber.delta.position
lo, hi = huber.delta.limits
span = (start_pos - 0.5, start_pos + 0.5)
print(f">>> huber.delta at {start_pos:.5f}, limits {lo}-{hi}, span {span[0]:.5f}..{span[1]:.5f}",
      flush=True)
assert lo <= span[0] and span[1] <= hi, "scan span outside motor limits -- refusing to run"

print(f">>> measurement_num before: {int(pv_registers.measurement_num.get())}", flush=True)
print(f">>> spec file: {default_spec_path()}", flush=True)
print(f">>> starting scan at {time.strftime('%H:%M:%S')}", flush=True)

t0 = time.time()
positions, columns = dscan_ophyd(
    huber.delta, -0.5, 0.5, 41, 1.0, det=lambda2M, att_ratio=5
)
elapsed = time.time() - t0

print(f"\n>>> scan finished at {time.strftime('%H:%M:%S')} in {elapsed:.1f} s", flush=True)
print(f">>> points returned      : {len(positions)}")
print(f">>> measurement_num after: {int(pv_registers.measurement_num.get())}")
print(f">>> huber.delta now      : {huber.delta.position:.5f} (started {start_pos:.5f})")
for label, values in columns.items():
    finite = [v for v in values if v is not None]
    print(f">>> {label}: n={len(finite)} min={min(finite)!r} max={max(finite)!r}")
print("STAGE B DONE")
