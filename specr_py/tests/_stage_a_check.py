"""Stage A: bring up the 8-ID-E session and verify the scan is wired, WITHOUT scanning.

Moves nothing. Reads only.
"""
print(">>> importing id8_common.startup ...", flush=True)
from id8_common.startup import *  # noqa: F401,F403

print(">>> importing dscan_ophyd ...", flush=True)
from id8_common.plans.align.ophyd_scan import (  # noqa: E402
    dscan_ophyd, default_spec_path, image_file_path, huber_motor_table,
    signal_map, detector_kind, SpecFile,
)

from apsbits.core.instrument_init import oregistry  # noqa: E402

huber = oregistry["huber"]
lambda2M = oregistry["lambda2M"]
pv_registers = oregistry["pv_registers"]

print("\n=== resolved configuration ===")
print("  spec path      :", default_spec_path())
print("  image path     :", image_file_path())
print("  spec exists?   :", __import__("os").path.exists(default_spec_path()))
print("  measurement_num:", int(pv_registers.measurement_num.get()))
print("  experiment     :", pv_registers.experiment_name.get())
print("  header/sample  :", pv_registers.header.get(), "/", pv_registers.sample_name.get())

print("\n=== motor to be scanned ===")
m = huber.delta
print("  name       :", m.name)
print("  prefix     :", m.prefix)
print("  position   :", m.position)
print("  limits     :", m.limits)
print("  egu        :", m.motor_egu.get())
target_lo, target_hi = m.position - 0.5, m.position + 0.5
print("  scan span  : %.5f .. %.5f" % (target_lo, target_hi))
lo, hi = m.limits
print("  within lim :", lo <= target_lo and target_hi <= hi)

print("\n=== detector ===")
print("  kind       :", detector_kind(lambda2M))
print("  columns    :", [lab for lab, _ in signal_map(lambda2M)])
print("  hdf path   :", lambda2M.hdf1.file_path.get())
print("  write mode :", lambda2M.hdf1.file_write_mode.get(as_string=True))
print("  capture    :", lambda2M.hdf1.capture.get(as_string=True))

names, positions = huber_motor_table()
print("\n=== #O / #P motor table (%d axes) ===" % len(names))
for n, p in zip(names, positions):
    print("  %-14s %.5f" % (n, p))

print("\nSTAGE A OK - nothing was moved.")
