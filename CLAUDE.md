# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Bluesky/Ophyd instrument control for the APS 8-ID beamline (XPCS). Originally built on the [`apsbits`](https://BCDA-APS.github.io/BITS/) framework (Bluesky Instrument Template System), and steadily moving off it.

Packages under `src/`:

- `id8_common/` — **the live code.** Devices, plans, utilities, configs, and both session entry points (`startup.py`, `startup_ophyd.py`). Almost all work happens here.
- `user_plans/<cycle>/<experiment>/` — per-experiment YAML (`sample_info.yaml`, `measurement_info.yaml`, `trio_measurement_info.yaml`) and beamtime scripts. Which folder is live comes from `configs/experiment.yml`, not from a hardcoded path.
- `legacy/id8_i/`, `legacy/id8_e/` — the former per-instrument packages. Nothing live imports them. Don't edit them, and don't take them as a model.
- `id8_common_dev/` — parallel staging copy of `id8_common`, and by now well behind it. Treat it as stale rather than as a sibling to keep in sync; ask before touching it.

Two hand-written modules deliberately replace `apsbits` machinery — see them before assuming an `apsbits` answer applies:

- `id8_common/registry.py` — the device registry and YAML loader (`oregistry`, `safe_make_devices`, `get_connected_device`, `get_ophyd_object`), replacing `apsbits.core.instrument_init`.
- `id8_common/expt_config.py` — the `expt` session-state object, replacing `pv_registers`.

`apsbits` still provides the Guarneri device YAML shape, the RunEngine factory (`init_RE`), and queueserver scaffolding on the Bluesky path only.

## Architecture: how a session boots

**⚠ `src/id8_i/` and `src/id8_e/` no longer exist** — they were moved under `src/legacy/` and are not imported by anything live. The working code is `src/id8_common/`, which now owns both entry points:

* `id8_common/startup.py` — the Bluesky session (RunEngine, databroker, apsbits/guarneri). `~/bin/start_bluesky.sh`.
* `id8_common/startup_ophyd.py` — an Ophyd-only session: same devices and the same plans, no RunEngine, no apsbits. `~/bin/start_ophyd.sh`. Both populate the same `id8_common.registry.oregistry`, so plan modules work unchanged under either.

The rest of this section describes the Bluesky path; `startup_ophyd.py` does the same thing without steps 2–4. The order is:

1. Load `<package>/configs/iconfig.yml` — declares the databroker catalog, RunEngine defaults, baseline/BEC config, optional SPEC/Tiled callbacks, ophyd timeouts.
2. `init_instrument("guarneri")` → returns `(instrument, oregistry)`. The instrument is a Guarneri device manager. Devices built here are then re-registered into `id8_common.registry.oregistry`, which is what plan code looks up by name (`oregistry["eiger4M"]`) — it is the registry both startup paths share, and unlike the guarneri object it supports `in`, `len()` and iteration.
3. `init_bec_peaks`, `init_catalog`, `init_RE` → wire up the BestEffortCallback, the databroker catalog subscriber, and the `RE` itself.
4. Either import plans via `*` (queueserver mode, detected by `running_in_queueserver()`) or via prefixes (`bp`, `bps`) for interactive use.
5. `safe_make_devices(file="devices.yml", device_manager=instrument)` is called **once per YAML file** (typically `devices.yml`, `ad_devices.yml`, and `devices_aps_only.yml`). Each YAML maps a fully-qualified class path → list of instance dicts (`name`, `prefix`, and class-specific kwargs).
6. Area detectors need post-creation wiring: `ad_setup(oregistry["eiger4M"], iconfig)` from `id8_common.devices.area_detector`, plus `stream_rois(det)` from the per-instrument or common `utils/misc.py`.
7. Import the plans last, so device names are already in `oregistry` when plan modules read them at import time. Prefer `get_connected_device("name")` *inside* a plan function over a module-level `oregistry.get("name")` — a device skipped at startup binds `None` at import and then fails much later as an `AttributeError` on `None`, far from the cause.

Plans are organised as `id8_common/plans/{acquire,align,set}/`. Both startup files import the same ones; `startup_ophyd.py` imports the Ophyd-only subset (`ophyd_scan.dscan_ophyd` rather than `scan_8id`).

## Conventions to know before editing

- **Device classes live in `id8_common/devices/`**, not per-instrument. The per-instrument `devices/` dirs only contain a small `registers_device.py` shim. New hardware support belongs in `id8_common/devices/`, and the YAML entry that instantiates it goes in the relevant instrument's `configs/devices.yml`.
- **`expt`** (`id8_common/expt_config.py`) is the session-state store, and it replaced `pv_registers` on 2026-09-06. One attribute namespace over four backends: static settings from `configs/experiment.yml`, per-measurement run state loaded from `measurement_info.yaml`, persistent state (sample index, file name, mesh positions) in `state/run_state.yml`, and exactly one EPICS-backed field. Plans read `expt.acq_time`, `expt.cycle_name`, … and never touch registers directly. See `docs/configuration.md`.
- **`measurement_num` is the one value still in EPICS** (`8ideSoft:Reg1`, via `pv_registers.measurement_num` — see `PV_FIELDS` in `expt_config.py`). It is the `NNNN` in every file name, nothing checks whether a name is taken, and `state/run_state.yml` is gitignored and per-checkout — so a counter stored there could reset and silently overwrite data. `run_state.yml` keeps a mirror used only to push the register back up if it comes back lower. Don't move this one into Python state. Note the counter is shared by `data/A####` (acquisitions) *and* `data/bluesky/A####` (`dscan_ophyd`/`scan_csv`), so the highest number on disk may be under `data/bluesky/`.
- **The rest of `EpicsPvStorageRegisters` is dead.** The class still declares the old `StrRegN`/`RegN` Components and the device is still built from `devices.yml`, but nothing in `id8_common` reads or writes them any more.
- **`host_on_aps_subnet()`** from `apsbits.utils.aps_functions` gates loading of `devices_aps_only.yml`. Devices that talk to real APS hardware go there. Devices that work offline (sim, soft IOCs you bring up locally) go in `devices.yml`.
- **`safe_make_devices`** (`id8_common/registry.py`) is the default and only device loader on both startup paths — not opt-in. It builds each YAML entry in its own `try/except` and checks it with a plain `device.wait_for_connection()` through a bounded 4-worker pool; anything that fails to build or connect is skipped with a warning and named in the offline banner, so one dead IOC cannot starve the devices declared after it. It does **not** use `devices_heartbeat.yml` or `heartbeat_pv` — that mechanism is retired. `id8_common/utils/safe_devices.py` is a thin adapter for the Bluesky/guarneri path that also binds each device as a bare name in `__main__`.
- **Queueserver vs. interactive imports** are deliberately different. In QS the startup does `from bluesky.plans import *` (all plans must be importable by name through the QS permissions). In interactive mode it uses `bp` / `bps` prefixes. New plans must be importable cleanly in both paths.
- **Area-detector YAML (`ad_devices.yml`)** uses `apstools.devices.area_detector_factory.ad_creator` with per-plugin class overrides from `id8_common.devices.area_detector` (Eiger/Lambda variants of cam, codec, image, hdf1, overlay, process, pva, roi1-4, stats1-4, transform1). HDF5 `read_path_template`/`write_path_template` are real beamline paths under `/gdata/dm/8IDI/<cycle>/` and need updating per run cycle.

## Running it

**This code only runs on a host that can see the beamline PVs** — `pearl` or `amber`, not `kouga`. The checkout is NFS-shared, so an edit made anywhere is live on `pearl` immediately; check for a running session before changing anything mid-beamtime.

The two entry points, both using the `8ide_bits_test` conda environment and both sourcing `/home/dm_id/etc/dm.setup.sh` for APS Data Management. The environment is **`8ide_bits_test`** despite the name -- `8id_bits` is a filesystem copy whose `ipython` shebangs back into it, so installing there changes nothing about what runs. See README.md, "Which conda environment":

```bash
~/bin/start_bluesky.sh     # ipython -i -c "from id8_common.startup import *"
~/bin/start_ophyd.sh       # ipython -i -c "from id8_common.startup_ophyd import *"
```

Conda environment, if it ever has to be rebuilt:
```bash
conda create -y -n 8ide_bits_test python=3.11 pyepics apsu::aps-dm-api
conda activate 8ide_bits_test
pip install -e ."[all]"
```

**⚠ The queueserver is retired.** Its launchers now sit in `scripts/Archive/`
(`id8_i_qs_host.sh`, `id8_e_qs_host.sh`) and resolve `src/id8_i/configs` and
`src/id8_i/qserver`, which moved to `src/legacy/` — so they cannot work as
written either. Don't revive them on your own initiative; ask first.

Linting/formatting is enforced by pre-commit (`ruff` + `ruff-format`, line length 120, py311). CI runs only pre-commit (`.github/workflows/pre-commit.yml`). Install once with `pre-commit install`.

Tests: the README documents `pytest -vvv --lf ./src` and `pyproject.toml` configures pytest with `--import-mode=importlib -x`, but the `test_*.py` files in `plans/` and `user/` are beamline scripts, not unit tests — there is no real unit-test suite. Don't claim a change is verified just because `pytest` collected nothing.

## Things not to do silently

- Don't verify a change by claiming the tests pass. There is no unit-test suite; the real check is a session boot plus a dry run on `pearl` (`dry_run_measurement_info()` / `dry_run_trio_measurement_info()`), which validate everything and move nothing.
- Don't add your own try/except around device loading — `safe_make_devices` already isolates every entry, and a second layer just hides which device failed. A disconnected device is *expected* to load-and-skip: the beamline must start with hardware missing, and the error must surface when that specific device is used, not at startup for everything else. That is also why several entries in `devices.yml` are commented out rather than deleted.
- Don't reorder the startup blocks. Plan modules resolve devices at import time and rely on `safe_make_devices` having already run. The order of the star-imports also matters: it decides which module a name at the interactive prompt comes from.
- Don't push to a `*_dev` package thinking it's a feature branch — `id8_common_dev` is checked-in source code. It is now well behind `id8_common`; ask before touching it rather than mirroring changes into it by reflex.
- Don't hand-edit `state/run_state.yml` while a session is running, and don't reset the measurement counter. See the `expt` and `measurement_num` notes above.
