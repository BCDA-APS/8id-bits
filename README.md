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

Either one activates the `8ide_bits_test` conda environment (see the note
below) and sources the APS Data Management setup. `start_ophyd.sh` is the
normal way in; see [Starting a session](docs/starting-a-session.md) for what it
loads and how to read the output.

### Which conda environment

**`8ide_bits_test` is the environment that runs.** The name says "test"; it is
not. Three environments exist and only this one is real for our purposes:

| env | what it is |
|---|---|
| `8idi_bits_test` | the original, 2025-12-04. Used by `start_bluesky_8idi.sh`. |
| `8ide_bits_test` | `conda create --clone` of the above. **What every session actually runs.** |
| `8id_bits` | a filesystem copy of `8ide_bits_test`. Supplies a launcher script and nothing else. |

`8id_bits` was copied rather than cloned -- its `conda-meta/history` carries
only the parent's lineage, with no `--clone` line of its own -- so the copy kept
`8ide_bits_test`'s shebang:

```
$ conda activate 8id_bits
$ head -1 $(which ipython)
#!/home/beams/8IDIUSER/.conda/envs/8ide_bits_test/bin/python3.11
```

So until 2026-09-11 the launchers activated `8id_bits`, and then `ipython`
handed the session straight to `8ide_bits_test`'s interpreter and
site-packages. `python` and `ipython` in the same shell read different
`site-packages` directories.

Nothing was ever missing: the two are identical today -- 79 conda packages with
the same fingerprint, 277 pip dist-infos, zero difference either way, and both
`.pth` files pointing at the same `bluesky/src` and the same
`nexus_xpcs_aps_95ab368/src`. The hazard is the future one. A `pip install` into
`8id_bits` reports success and changes nothing about what runs, and the two
environments would then disagree with no symptom until an import fails.

`start_bluesky.sh`, `start_ophyd.sh` and `start_bluesky_main.sh` now activate
`8ide_bits_test` directly, which is self-consistent: `python` and `ipython`
resolve to the same interpreter and the same `site-packages`.
**Install into `8ide_bits_test`, and check `head -1 $(which ipython)` after any
`conda create --clone` -- a clone that keeps its parent's shebang is exactly how
this happened.**

## I want to…

| … | Go to |
|---|---|
| start a session at the beamline | [Starting a session](docs/starting-a-session.md) |
| see what loaded, or add and remove a device | [Devices](docs/devices.md) |
| change the experiment, or find where a setting or session state lives | [Configuration](docs/configuration.md) |
| read or set `expt` at the prompt | [Using `expt`](docs/using-expt.md) |
| define a protocol, run a measurement, or look up a detector mode | [Running measurements](docs/running-measurements.md) |
| add or remove a field in the NeXus metadata file | [Adding metadata fields](docs/adding-metadata.md) |
| run an alignment scan, and plot it live or afterwards | [Running and viewing scans](docs/running-and-viewing-scans.md) |
| work out why something failed, including an analysis job | [Troubleshooting](docs/reference/troubleshooting.md) |

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

## Where the code lives

Most of what runs at 8-ID is in this repository. Two pieces are not, and both
are installed from a checkout rather than from PyPI — so "upgrading" either one
means pulling in that directory, and it takes effect immediately.

| | what it does | on this machine | upstream |
|---|---|---|---|
| **this repo** | devices, plans, configuration | `~/bluesky` | [BCDA-APS/8id-bits](https://github.com/BCDA-APS/8id-bits) |
| **nexus_xpcs_aps** | writes the NeXus metadata file | `~/Documents/Miaoqi/nexus_xpcs_aps_95ab368`, `pip install -e` into `8ide_bits_test` | [AZjk/nexus_xpcs_aps](https://github.com/AZjk/nexus_xpcs_aps), branch `mc_refact` |
| **BLUETELLA** | the scan viewers | `~/Documents/BLUETELLA_9ID`, branch `8id_test` | 9-ID's, maintained by Peco Myint |

Both external packages are maintained by someone else. Before changing either,
read the "when to contact" section on the relevant page —
[NeXus writer](docs/reference/nexus-writer.md) or
[scan viewer](docs/reference/scan-viewer.md) — because a local-only patch is a
fork we then maintain for ever.

**⚠ Neither is version-pinned.** A `git pull` in one of those directories
changes what the beamline does, with no version bump and no warning. Treat it
as a change to this instrument, not a routine update, and re-run one
measurement per detector mode afterwards.

## Known issues

Short list of things that are known broken or unfinished. Details, with file
and line, are on the page named in each row.

| | impact | where |
|---|---|---|
| `eiger4M` Internal Enable can hang indefinitely | stalls the session with the shutter open; intermittent | [Troubleshooting → an acquisition hangs](docs/reference/troubleshooting.md#an-acquisition-hangs) |
| several other detector waits are unbounded | same failure mode, not yet seen | [Troubleshooting → an acquisition hangs](docs/reference/troubleshooting.md#an-acquisition-hangs) |
| Rigaku ZDT `.bin.000` fails in DM analysis | no result file; the same data analyses fine locally | [Data Management](docs/reference/data-management.md) |
| ~100 Rigaku 3M pixels are hot and unmasked | 88% of the dark count rate at 7 keV; raising the threshold does not remove them | [Rigaku 3M dark noise](docs/reference/rigaku-3m-dark-noise.md) |
| no Eiger measurement has been verified with beam | plumbing proven, numbers not | — |
| two commented-out metadata lines have no schema node | uncommenting either silently loses the whole metadata file | [Adding metadata fields](docs/adding-metadata.md) |

## Conventions in these pages

* `expt` is the single configuration object. If a page says "`expt.cycle_name`",
  you can type exactly that at the session prompt.
* Paths are relative to the top of this repository — except the shorthands
  used throughout: `configs/…` and `plans/…` are inside `src/id8_common/`, and
  `user_plans/…` is inside `src/`. The pages themselves live in `docs/`.
* Anything marked **⚠** is a known trap, not a suggestion.
