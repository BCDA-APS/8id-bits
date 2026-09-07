# Troubleshooting

[← index](../../README.md)

Known failure modes, what they look like, and what to do.

## The session will not start

**A device is offline.** It should not stop startup — you get
`*** Devices not online: [...] ***` and the session continues. If startup
actually *aborts*, something else is wrong; read the traceback for the file, and
see [Devices](../devices.md).

**`experiment.yml not found`** — the file is required and holds the per-experiment
settings. See [Configuration](../configuration.md).

## A plan fails immediately

```
KeyError: 'eiger4M' is not in the device registry -- skipped at startup
```
The device did not connect. Check the red banner from startup, then `caget` the
real PV from pearl. This error is by design: it names the device instead of
surfacing later as `AttributeError: 'NoneType' ...`.

```
AttributeError: 'acq_time' has not been set yet
```
The run state has not been populated. Use `run_measurement_info()`, or call
`expt.set_measurement(...)` first for a manual `det_acq_series()`.

```
FileNotFoundError: .../user_plans/2026-3/comm202609 does not exist
```
`cycle_name` or `experiment_name` in `configs/experiment.yml` does not match a
real folder. The message names both settings.

```
RuntimeError: 'measurement_num' lives in pv_registers, which is not available
```
`pv_registers` did not connect *and* `state/run_state.yml` holds no mirrored
value to fall back on — a fresh checkout plus a dead soft IOC. Fix the IOC, or
set `expt.measurement_num` explicitly, past the highest number already on disk,
before acquiring.

## An acquisition hangs

**Eiger Internal Enable** — almost certainly `acq_period` too short; the detector
drops software triggers. See [Detector modes](detector-modes.md#eiger4m).

**Eiger External Series** — `trigger_period` shorter than one segment. The next
pulse lands mid-segment, is dropped, and the plan waits for a trigger already
spent.

**Anything else** — only some of the waits that used to abandon a run silently
are bounded. Bounded today:

* the Rigaku state waits (`_wait_for_state`) — raise on a timeout, or on the
  detector reaching Error/Disconnected/Aborted;
* the Rigaku HDF1 drain — raises on a timeout only, it has no fault check;
* the HDF drains in `acquire_eiger_external` and `acquire_lambda_external` —
  raise on a timeout or on a cam fault.

Everything else still spins, with no deadline and no fault check:

* the cam-busy loops in `acquire_eiger_internal`, `acquire_eiger_internal_enable`,
  `acquire_eiger_external` (both External modes) and both Lambda acquire
  functions — `acquire_lambda_external` has two back to back, and the first has
  no `sleep()` at all, it busy-spins on `run_hooks()`;
* the HDF drains in `acquire_eiger_internal`, `acquire_eiger_internal_enable`
  and `acquire_lambda_internal` — so an Internal-mode run stuck at the
  file-writing stage hangs there too.

**Recovering:** Ctrl+C. The cleanup path blocks the beam, stops softglue and the
detector, and writes metadata for the partial dataset.

## The measurement number looks wrong

`measurement_num` is the NNNN in every folder and file name (`A0061_Test_a0002_…`).
Its store is the EPICS register `8ideSoft:Reg1`, read through
`expt.measurement_num`; `state/run_state.yml` keeps only a mirror of the last
value this checkout saw. It counts up, and nothing downstream checks whether a
name is already taken.

**Red: `pv_registers is unavailable (…)`** — the register could not be reached,
so the session falls back to the `run_state.yml` mirror. That mirror is
per-checkout and gitignored, so the counter is **not** shared with another
session or with non-Bluesky tools this session: two of them will hand out the
same number. Check what is already on disk before acquiring, and fix the IOC.

**Red: `pv_registers.measurement_num reads N, but run_state.yml last saw M`** —
the register went backwards (soft-IOC restart, or Reg1 cleared by hand). The
first read of the session pushes it back up to M and says so. The correction
only ever goes upward; a register *ahead* of the mirror is always taken as
correct, so this warning does not appear for the normal case of another session
having advanced it.

**⚠ A stale counter overwrites data, and two naming streams share it.**
`det_acq_series()` writes under `data/`, while the align scans (`dscan`,
`ascan`, `d2scan`, `a2scan`, `mesh`, `dmesh`, the lups) write under
`data/bluesky/` — both trees name their files from `gen_folder_prefix()`, so the
highest number on disk may be under `data/bluesky/`, not `data/`.

How far one scan moves the counter depends on which session you are in, because
the same plain name means different code in each:

| you typed | `start_bluesky.sh` | `start_ophyd.sh` |
|---|---|---|
| `dscan`, `ascan`, `mesh`, `x_lup`, … | `scan_8id.py` (Bluesky generator) | `ophyd_scan.py` |
| `dscan_ophyd`, `ascan_ophyd`, … | `ophyd_scan.py` | same, as aliases |

`scan_8id.py`'s versions take **two** numbers when `save_img=1` and none when it
is 0. `ophyd_scan.py`'s take exactly one per call either way. So the same typed
command can move the counter by 0, 1 or 2 depending on the session and the
flag — check both trees before trusting it:

```bash
ls /gdata/dm/8ID/8IDE/<cycle>/<expt>/data/         | sort | tail
ls /gdata/dm/8ID/8IDE/<cycle>/<expt>/data/bluesky/ | sort | tail
caget 8ideSoft:Reg1
```

```python
expt.measurement_num = 1234    # writes Reg1 and the mirror
```

A scan at least stops itself: `scan_csv` refuses to append to a `.csv` that
already exists (`FileExistsError`), before the scan moves anything. An
acquisition has no such guard — it overwrites.

## Analysis produced nothing

**⚠ A green DM job does not mean the analysis ran.** DM can report
`status: done / stage: 12-DONE / gladierStatus: SUCCEEDED` with no
`*_results.hdf` produced.

The rest of this section describes code outside this repo — the Gladier compute
function on Polaris and the `boost_corr` CLI. It is what was observed at the
beamline, and nothing in this tree can confirm it.

Check for the actual output:

```bash
ls /gdata/dm/8ID/8IDE/<cycle>/<expt>/analysis/Multitau/*_results.hdf
```

If the folder holds only `boost_corr.log` (0 bytes), `boost_corr_err.log` and
`corr_metadata_output.json`, the correlation failed. Read the err log.

**Known upstream bug (as of 2026-09-06, reported):** the Gladier compute function
on Polaris launches `boost_corr` with **no arguments** —
`subprocess.run(cmd.split(" "), shell=True, ...)` runs `/bin/sh -c boost_corr`
and discards the rest. The error log contains only the argparse usage block. The
function records the failure in a field nothing inspects, so DM still reports
success. Fix is upstream, with DM support. The one-line cause is
`subprocess.run(<list>, shell=True)`, which on POSIX runs only `args[0]`.

**Workaround — run the correlation directly on amazonite:**

```bash
/home/beams/8IDIUSER/bin/boost_corr_bin \
    -r  /gdata/dm/8ID/8IDE/<cycle>/<expt>/data/<name>/<name>.h5 \
    -q  /gdata/dm/8ID/8IDE/<cycle>/<expt>/data/eiger4m_qmap_default.hdf \
    -o  /gdata/dm/8ID/8IDE/<cycle>/<expt>/analysis/<analysis type> \
    -t  Multitau -i 0
```

* `-o` is a **directory**; boost_corr names the file itself.
* `--meta-fname` is **not needed** — it finds `<name>_metadata.hdf` beside the raw file.
* `-i 0` selects GPU 0; `-i -1` forces CPU.
* `-u SUFFIX` distinguishes reanalysis variants.

**⚠ Check the qmap exists.** A protocol naming a qmap that is not in
`<expt>/data/` fails at analysis, not at acquisition.

## The metadata looks wrong

**Wrong pixel size** — pixel size is per detector in `device_position.yaml`
(`pixel_size:`), not global.

**Wrong detector name for a Rigaku run** — `detector_name` should be
`rigaku3M`, `rigaku3M_ftf` or `rigaku3M_epics`, matching the output format.

**Missing NeXus attributes** — if a metadata file has objects but none carry
attributes, that is the shallow-copy bug (fixed 2026-09-05). Files written before
that date, from the second measurement of a session onward, are affected: they
have no NeXus typing, and the damage propagates into `*_results.hdf`. Nothing
recovers them short of rewriting the metadata.

## Data volume is enormous

A Rigaku ZDT run at a low threshold records a noise flood — up to 1.5 TB for
100 000 frames. Raise the lower threshold. See
[Detector modes](detector-modes.md#-sparsified-size-depends-on-the-threshold).

## Checks worth running

```python
expt                                  # which experiment am I in
expt.user_plan_dir                    # which plan folder
expt.measurement_num                  # the shared counter (uncached read of Reg1)
oregistry.device_names                # what loaded
dry_run_measurement_info()            # validate without moving
```

```bash
caget 8ideSoft:Reg1                    # the measurement counter
caget 8idEiger4m:cam1:DetectorState_RBV
caget 8ideSoft:fastshutter:State_RBV   # the honest shutter readback
dmjob.sh <UUID>                        # DM job status
```

**⚠ `State` is what was last asked for, not where the blade is.** The two
disagree whenever softglue is driving the shutter (`Lock` on "Override"), which
is every external-trigger acquisition. `State_RBV` is the authoritative blade
position, and its enum is the honest one: "Open" / "Closed".
