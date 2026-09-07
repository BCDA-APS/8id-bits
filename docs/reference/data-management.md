# Data Management (DM) and analysis jobs

What gets submitted after a measurement, how to find a job again, and how to
tell a real success from a job that reported success and produced nothing.

## The 30-second version

```
  det_acq_series() finishes one measurement
          │
          ├─► writes  <name>.h5 / .bin.000       the data
          ├─► writes  <name>_metadata.hdf        the NeXus metadata
          │
          ▼
  dm_run_job()                        utils/dm_util.py
          │  builds argsDict from expt: experimentName, filePath,
          │  qmap, analysisMachine, gpuID, type, useSubdir
          │  ...plus four hardcoded: demand=True, saveG2=False,
          │  download=False, normalizeFrame=False
          ▼
  workflowProcApi.startProcessingJob(dmuser, workflow_name, argsDict)
          │
          │  prints:  Job e2ae89dd-01f8-45d8-a79c-8f73a4f7ff89
          │           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          │           the uuid -- printed only, NOT stored anywhere
          ▼
  DM runs the workflow (stages 01…NN) on analysisMachine
          │
          ▼
  <expt>/analysis/<type>/<name>_results.hdf        (note: resultS, plural)
```

## Finding the uuid again: `dm_jobs.log`

Every submitted job is appended to a plain text log at the experiment root,
beside `data/` and `analysis/`:

```
/gdata/dm/8ID/8IDE/<cycle>/<expt>/dm_jobs.log
```

```
# date       time      job_uuid  measurement  machine  workflow  data_file
2026-09-07 17:13:02  e2ae89dd-…-7ff89  A0101_Test_a0002_f000100_r00001  polaris  xpcs8-boost-corr  /gdata/…/A0101_Test.h5
```

```bash
grep A0101 dm_jobs.log                                # the job for one measurement
awk '!/^#/{print $3}' dm_jobs.log                     # every uuid
dmjob.sh $(tail -1 dm_jobs.log | awk '{print $3}')    # status of the most recent job
```

The uuid is the only handle DM gives you on a job, and before this it was
printed to the session and nowhere else — gone with the scrollback. It is a log
rather than a field in the NeXus file because the metadata file is written
*before* the job is submitted, so the uuid does not exist yet; and because one
greppable file beats opening sixty HDFs.

Writing the log never interrupts an acquisition: if the mount is full or
read-only the line is dropped with a warning and the run continues (the uuid is
still on screen).

If the log is missing a job — an older run, or a failed write — fall back to
`dm-list-processing-jobs` (below).

## Checking a job

The beamline wrappers live in `~/bin` and source the DM environment for you:

| command | what it tells you |
|---|---|
| `dmjob.sh <uuid>` | everything about one job, pretty-printed |
| `dmbatchjobstat.sh <uuid>` | just `countFiles` / `nCompleted` / `nFailed` |
| `dm_check_stage.sh <uuid> <stage-id>` | the actual shell command a stage ran. The argument is a stage **id string** such as `02-START` or `07-HPC`, not a number — read them off `dmjob.sh` output first |
| `dmjob_reprocess.sh <uuid>` | same as `dmjob.sh`, against the reprocess service |

Underneath, all of them are `dm-get-processing-job --id=<uuid>`.

To find a job when you no longer have the uuid:

```bash
source /home/dm_id/etc/dm.setup.sh
dm-list-processing-jobs --owner=$USER | tail -20
```

The CLI documents no sort order, so `tail` is a guess — match on the file name
in the output rather than assuming the last line is the newest.

**⚠ `dm_check_job.sh` and `dm_check_stage.sh` source `/home/dm/etc/dm.setup.sh`,
not `/home/dm_id/`.** `dmjob.sh` uses the `dm_id` one. If a wrapper returns
nothing or fails to authenticate, that is the first thing to check.

## What to read in the output

| field | meaning |
|---|---|
| `status` | `pending` → `running` → `finalizing` → `done` / `failed` (also `skipped`, `aborting`, `aborted`). Never the word SUCCEEDED — that comes from the Globus/Gladier layer inside a stage, which is why the two can disagree |
| `stage` | which workflow step it is on |
| `nCompletedChildProcesses` / `nFailedChildProcesses` | per **child process**, not per file. `dmbatchjobstat.sh` greps the shorter `nCompleted`/`nFailed`, which matches these by substring |
| `childProcesses` | the actual commands, with their stdout and stderr |
| `startTime` / `endTime` | wall clock; a job that "finished" in 2 s did nothing |

`childProcesses` is where the real answer usually is — but note that a failure
on the HPC side may leave nothing there at all, with the actual message in a log
on the compute node. A stage that "finished" in a couple of seconds is the
classic signature of a command that never really ran.

## ⚠ The trap: SUCCEEDED does not mean analysed

**DM can report a job as succeeded when the analysis produced no output at all.**
This has bitten this beamline: a bug in the Gladier compute function
(`xpcs_boost_corr.py`, `subprocess.run(cmd.split(" "), shell=True)` — with
`shell=True` and a list, every argument after the first is discarded) meant
`boost_corr` ran with no arguments. It exited non-zero with a two-line argparse
usage error, but that went to **stderr on the compute node**
(`boost_corr_err.log`), not into the DM child-process output, and the Gladier
layer still reported SUCCEEDED. No result file was written.

**So the only reliable check is the output file, not the job status:**

```bash
ls -l /gdata/dm/8ID/8IDE/<cycle>/<expt>/analysis/<type>/<name>_results.hdf
```

If the job says succeeded and that file is absent or implausibly small, read
`childProcesses` in `dmjob.sh` output and look at what was actually executed.

## Running the analysis by hand

To take DM out of the loop entirely and check whether the *analysis* works:

**⚠ `-o` is an output DIRECTORY, not a file.** boost_corr always names the file
itself, `<name>_results.hdf`. Passing a filename makes it `mkdir` a directory
with that name and put the real file inside — there is one of those in the
current experiment, `analysis/Multitau/A0023_…_r00001/`, left by exactly this
mistake.

```bash
# on an analysis machine (e.g. amazonite)
boost_corr_bin -r <data>.h5 -q <qmap>.hdf \
    -o /gdata/dm/8ID/8IDE/<cycle>/<expt>/analysis/Multitau \
    -i 0 -t Multitau
# writes <name>_results.hdf into that directory
```

If that produces a result and DM does not, the problem is in DM or the Gladier
function, not in the data or the qmap.

## Turning analysis off

Set `analysis_machine: none` in `configs/experiment.yml`. Data is written, no
job is submitted. Other accepted values: `polaris` (uses `gpuID: 0`), `local`
(picks a beamline analysis box via `get_machine_name()`), or any hostname typed
literally — both of the latter use `gpuID: -2`. **⚠ `local` is a fresh
`random.choice` on every submission** (currently adamite or califone, califone
twice so ~2/3 of the time), so consecutive measurements in one run can land on
different machines.

## Known rough edges in `dm_util.py`

* An unrecognised `expt.det_name` falls through the suffix `if/elif` chain and
  dies with `UnboundLocalError: cannot access local variable 'filepath'`. In
  practice you will not see it: `validate_detector_mode()` rejects an unknown
  detector with a message naming it long before submission. It only bites if you
  set `expt.det_name` by hand and call `det_acq_series()` directly.
* The matching `use_subfolder` branch is unreachable — `expt.use_subfolder` is
  coerced on load, so a bad value in `experiment.yml` raises
  `ValueError: expected yes or no` at `expt.reload()` instead.
* `rigaku3M_ftf` output is submitted as `<name>.h5.000`. The `.000` split is
  confirmed; that the IOC appends it to a `.h5` name the same way it did to
  `.bin` is inferred and not yet observed — check it on the first fast-transfer
  run. `boost_corr` also has no `.h5.000` branch, so fast-transfer data is not
  analysable today. See [detector modes](detector-modes.md).

## Related

* [Troubleshooting](troubleshooting.md) — including "analysis produced nothing"
* [Adding metadata fields](../adding-metadata.md) — what is in the file DM analyses
* [Configuration](../configuration.md) — `analysis_machine`, `workflow_name`, `use_subfolder`
