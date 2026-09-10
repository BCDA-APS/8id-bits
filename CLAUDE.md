# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Bluesky-based instrument control for the APS 8-ID beamline. The repo houses **two operational instruments** plus a **shared library**, all built on top of the [`apsbits`](https://BCDA-APS.github.io/BITS/) framework (Bluesky Instrument Template System).

Packages under `src/`:

- `id8_i/` — 8-ID-I (XPCS endstation)
- `id8_e/` — 8-ID-E
- `id8_common/` — devices, plans, and utilities shared by both
- `id8_common_dev/` — parallel staging copy of `id8_common`. Many files diverge between the two; when editing `id8_common`, check whether a sibling exists in `id8_common_dev` that should change too (or vice versa). There is no automatic sync.

The `apsbits` framework provides Guarneri-style declarative device YAML, the RunEngine factory (`init_RE`), and queueserver scaffolding. Treat it as upstream — most architectural decisions about *how* devices and plans are wired come from `apsbits`, not this repo.

## Architecture: how a session boots

Each instrument is invoked by importing its `startup.py` (e.g. `from id8_i.startup import *`). The startup module orchestrates everything in this order:

1. Load `<package>/configs/iconfig.yml` — declares the databroker catalog, RunEngine defaults, baseline/BEC config, optional SPEC/Tiled callbacks, ophyd timeouts.
2. `init_instrument("guarneri")` → returns `(instrument, oregistry)`. The instrument is a Guarneri device manager; `oregistry` is the global ophyd device registry that everything looks up by name (`oregistry["pv_registers"]`).
3. `init_bec_peaks`, `init_catalog`, `init_RE` → wire up the BestEffortCallback, the databroker catalog subscriber, and the `RE` itself.
4. Either import plans via `*` (queueserver mode, detected by `running_in_queueserver()`) or via prefixes (`bp`, `bps`) for interactive use.
5. `make_devices(file="devices.yml", device_manager=instrument)` is called **once per YAML file** (typically `devices.yml`, `ad_devices.yml`, and `devices_aps_only.yml`). Each YAML maps a fully-qualified class path → list of instance dicts (`name`, `prefix`, and class-specific kwargs).
6. Area detectors need post-creation wiring: `ad_setup(oregistry["eiger4M"], iconfig)` from `id8_common.devices.area_detector`, plus `stream_rois(det)` from the per-instrument or common `utils/misc.py`.
7. Import the per-instrument plans last, so device names are already in `oregistry` when plan modules read them at import time (several plans do `pv_registers = oregistry["pv_registers"]` at module level).

`id8_e` and `id8_i` startup files diverge in style — `id8_e` uses the structured `id8_common/plans/{acquire,align,set}/` layout, while `id8_i` still imports flat `id8_i/plans/*.py` modules. Don't try to unify them without checking with the user.

## Conventions to know before editing

- **Device classes live in `id8_common/devices/`**, not per-instrument. The per-instrument `devices/` dirs only contain a small `registers_device.py` shim. New hardware support belongs in `id8_common/devices/`, and the YAML entry that instantiates it goes in the relevant instrument's `configs/devices.yml`.
- **`pv_registers`** (`oregistry["pv_registers"]`, class `EpicsPvStorageRegisters`) is the project's session-state store, backed by `8ideSoft:` `StrRegN`/`RegN` EPICS PVs (spec filename, sample name, sample positions, acq parameters, etc.). Plans commonly read/write it directly. Don't replace it with Python state — it survives across processes and is shared with non-Bluesky tools.
- **`host_on_aps_subnet()`** from `apsbits.utils.aps_functions` gates loading of `devices_aps_only.yml`. Devices that talk to real APS hardware go there. Devices that work offline (sim, soft IOCs you bring up locally) go in `devices.yml`.
- **`safe_make_devices`** (`id8_common/utils/safe_devices.py`) is an opt-in alternative to `make_devices`. It reads `devices_heartbeat.yml`, runs `caget` on each entry's `heartbeat_pv` in parallel, and silently skips devices whose IOC is offline. The 8-ID-E startup has it commented out — currently both instruments use plain `make_devices`, so a missing IOC will raise at startup. If you need partial loading, switch to `safe_make_devices` rather than try/excepting around `make_devices`.
- **Queueserver vs. interactive imports** are deliberately different. In QS the startup does `from bluesky.plans import *` (all plans must be importable by name through the QS permissions). In interactive mode it uses `bp` / `bps` prefixes. New plans must be importable cleanly in both paths.
- **Area-detector YAML (`ad_devices.yml`)** uses `apstools.devices.area_detector_factory.ad_creator` with per-plugin class overrides from `id8_common.devices.area_detector` (Eiger/Lambda variants of cam, codec, image, hdf1, overlay, process, pva, roi1-4, stats1-4, transform1). HDF5 `read_path_template`/`write_path_template` are real beamline paths under `/gdata/dm/8IDI/<cycle>/` and need updating per run cycle.

## Running it

Conda environment (one-time):
```bash
export ENV_NAME=bs_8id_main
conda create -y -n $ENV_NAME python=3.11 pyepics apsu::aps-dm-api
conda activate $ENV_NAME
pip install -e ."[all]"
```

The 8-ID-I shell starter (`scripts/start_bluesky_8idi.sh`) hard-codes a different env name (`8idi_bits_test`) and sources `/home/dm_id/etc/dm.setup.sh` for APS Data Management. It also adds `164.54.116.40` to `EPICS_CA_ADDR_LIST` for the robocart IOC. Use it as-is on beamline workstations; for local dev, run `ipython -i -c "from id8_i.startup import *"` (or `id8_e.startup`) directly.

Queueserver (per instrument):
```bash
./scripts/id8_i_qs_host.sh {start|stop|restart|status|console|run}
./scripts/id8_e_qs_host.sh ...
```
The hostname check inside the script means it must run on the same host that will own the QS process. Config files (`qs-config.yml`, `user_group_permissions.yaml`) live in `src/id8_{i,e}/qserver/`. The startup module is wired via `startup.startup_module: id8_i.startup` in the QS config — don't rename `startup.py` without updating both.

Linting/formatting is enforced by pre-commit (`ruff` + `ruff-format`, line length 120, py311). CI runs only pre-commit (`.github/workflows/pre-commit.yml`). Install once with `pre-commit install`.

Tests: the README documents `pytest -vvv --lf ./src` and `pyproject.toml` configures pytest with `--import-mode=importlib -x`, but the `test_*.py` files in `plans/` and `user/` are beamline scripts, not unit tests — there is no real unit-test suite. Don't claim a change is verified just because `pytest` collected nothing.

## Things not to do silently

- Don't move a device class out of `id8_common/devices/` into a per-instrument package — both instruments share it.
- Don't introduce try/except around device loading to "be safe." If an IOC is expected to be offline, gate the YAML entry with `heartbeat_pv` and switch the startup to `safe_make_devices`. If it should always be present, let the failure surface.
- Don't reorder the startup blocks. The plan modules do `oregistry["..."]` at import time and rely on `make_devices` having already run.
- Don't push to a `*_dev` package thinking it's a feature branch — `id8_common_dev` is checked-in source code, parallel to `id8_common`. Either edit the right one, or edit both deliberately.
