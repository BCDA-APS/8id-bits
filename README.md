# 8-ID XPCS instrument

Beamline control for the APS 8-ID XPCS instrument. The working code is
`src/id8_common`; this page is the way in.

Start here and follow the links. Every page is plain Markdown, so it reads on
GitHub, in an editor, or with `less`.

## Starting

On a beamline host that can see the EPICS network (`pearl` or `amber`):

```bash
~/bin/start_ophyd.sh      # Ophyd-only session -- devices and plans, no RunEngine
~/bin/start_bluesky.sh    # the same, plus Bluesky (RunEngine, databroker)
```

Either one activates the `8id_bits` conda environment and sources the APS Data
Management setup. `start_ophyd.sh` is the normal way in; see
[Starting a session](docs/starting-a-session.md) for what it loads and how to
read the output.

## I want to…

| … | Go to |
|---|---|
| start a session at the beamline | [Starting a session](docs/starting-a-session.md) |
| understand what got loaded, or add/remove a device | [Devices](docs/devices.md) |
| know where a setting lives, or change the experiment | [Configuration](docs/configuration.md) |
| read or set `expt` at the prompt, or work out why a field raised | [Using `expt`](docs/using-expt.md) |
| add a new device's reading to the NeXus metadata file | [Adding metadata fields](docs/adding-metadata.md) |
| check, debug or re-run a DM analysis job | [Data Management](docs/reference/data-management.md) |
| find where session state is kept | [Configuration → `state/run_state.yml`](docs/configuration.md#staterun_stateyml) |
| know why the counter lives in EPICS | [Configuration → the measurement counter](docs/configuration.md#the-measurement-counter-stays-in-epics) |
| define a protocol and run a measurement | [Running measurements](docs/running-measurements.md) |
| look up a detector mode or its timing rules | [Detector modes](docs/reference/detector-modes.md) |
| work out why something failed | [Troubleshooting](docs/reference/troubleshooting.md) |

## The 30-second version

```
   configs/experiment.yml          which cycle, which experiment, where data goes
   user_plans/<cycle>/<expt>/
       sample_info.yaml            what the samples are and where they sit
       measurement_info.yaml       what to measure and how
                 │
                 ▼
            expt  (id8_common/expt_config.py)     ← one object, everything reads it
                 │
                 ▼
        run_measurement_info()      expands YAML into measurements, validates, runs
                 │
                 ▼
        det_acq_series()            drives the detector, writes NeXus metadata,
                                    submits the DM analysis job
```

Every folder and file is named `A0061_<sample>_a0007_…`. That `0061` is
`expt.measurement_num`, and its store is the EPICS register `8ideSoft:Reg1`, not
a file in this repository — so a second session, or a fresh clone, carries on
from the same number instead of restarting at 0 over existing data.
Generated session state (sample index, mesh positions) lands in the gitignored
`state/run_state.yml`. Both are covered in [Configuration](docs/configuration.md).

Two commands cover most days:

```python
dry_run_measurement_info()   # validate and preview -- moves nothing
run_measurement_info()       # go
```

## Conventions in these pages

* `expt` is the single configuration object. If a page says "`expt.cycle_name`",
  you can type exactly that at the session prompt.
* Paths are relative to the top of this repository — except the shorthands
  used throughout: `configs/…` and `plans/…` are inside `src/id8_common/`, and
  `user_plans/…` is inside `src/`. The pages themselves live in `docs/`.
* Anything marked **⚠** is a known trap, not a suggestion.
