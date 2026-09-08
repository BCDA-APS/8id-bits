# Starting a session

[← index](../README.md)

## Start it

```bash
~/bin/start_ophyd.sh
```

That is the whole thing. It sources the APS Data Management setup
(`/home/dm_id/etc/dm.setup.sh` — `utils/dm_util.py` imports the `dm` package to
submit analysis jobs), activates the `8id_bits` conda environment, adds the
robocart IOC to `EPICS_CA_ADDR_LIST`, and drops you into IPython with the
instrument loaded.

You must run it on a machine that can see the beamline PVs — **pearl** or
**amber**. Analysis machines and office workstations cannot; the devices will
simply fail to connect.

## What you should see

```
[startup_ophyd] Starting Ophyd-only session (.../configs/iconfig.yml)
[startup_ophyd] Connecting devices (devices.yml, ad_devices.yml, devices_aps_only.yml) ...
[startup_ophyd] 59 device(s) connected
[startup_ophyd] eiger4M area-detector plugins configured
[startup_ophyd] lambda2M area-detector plugins configured
[startup_ophyd] rigaku3M area-detector plugins configured (warmup skipped)
[expt_config] <ExperimentConfig 'comm202609' cycle='2026-3' run_keys=0>
[startup_ophyd] hklpy2 diffractometer (psic) configured
[startup_ophyd] Importing plans ...
[startup_ophyd] Ready -- 59 device(s) connected, plans imported.
```

The `[expt_config]` line is worth reading every time — it tells you which
experiment the session thinks it is running. If that is wrong, everything
downstream writes to the wrong place. See [Configuration](configuration.md).

## When a device is offline

The session **still starts**. This is deliberate: one dead IOC must not stop you
from using the rest of the beamline. You get a red banner instead:

```
[startup_ophyd] 57 device(s) connected, 2 skipped

*** Devices not online: ['tetramm3', 'lambda2M'] ***
```

A skipped device is simply absent. Nothing crashes at startup; the failure
arrives later, at the point where a plan actually needs that device. A plan that
looks the device up per call names it:

```
KeyError: 'lambda2M' is not in the device registry -- skipped at startup
          (offline IOC or bad prefix), or the name is misspelled. See
          oregistry.device_names for what did load.
```

The older modules that bind at module scope instead
(`plans/align/scan_8id.py`, `utils/nexus_runtime.py`) hold `None` and fail with an
`AttributeError` on `None`, which names neither the device nor the registry.

That is the intended behaviour, not a bug. Read
[Devices → offline devices](devices.md#when-a-device-does-not-connect) for how
the check works and what it does *not* catch.

**⚠ `pv_registers` fails differently.** Nothing raises when it is missing. It
carries one value — `expt.measurement_num`, the `NNNN` in every folder and file
name, stored in `8ideSoft:Reg1` — and at first use the session prints a red
warning and falls back to the mirror of that number kept in
`state/run_state.yml`. The mirror is per-checkout, so the counter is no longer
shared with any other session or tool: check what is already on disk before
acquiring. (On a checkout with no mirrored value either, that first read raises
instead and tells you to set `expt.measurement_num` yourself.) See
[Configuration → the measurement counter](configuration.md#the-measurement-counter-stays-in-epics).

## Checking things by hand

```python
"eiger4M" in oregistry          # did it load?
len(oregistry)                  # how many devices came up
oregistry.device_names          # what is actually here
eiger4M.connected               # is it connected right now
eiger4M.component_names         # what parts does it have (motors, plugins, signals)
eiger4M.summary()               # every signal on it, with its PV and read/write mode
expt                            # which experiment am I in
```

Every device is also a plain name at the prompt — `eiger4M.cam.acquire.get()`
works without going through the registry.

**⚠ Do not conclude a device is dead from an Ophyd error alone.** Check the real
PV from a host that can see the beamline:

```bash
caget 8idEiger4m:cam1:Acquire_RBV
```

An Ophyd-side failure can equally mean a wrong prefix or a network/subnet
problem.

## There is no RunEngine

This session runs Ophyd directly. There is no `RE`, no `%wa` / `%ct` magics, no
databroker. Plans are ordinary Python functions — call them:

```python
run_measurement_info()             # not RE(run_measurement_info())
dscan_ophyd(huber.delta, -0.5, 0.5, 41, 1.0)
```

Scans write a `.csv` next to the detector `.h5`, flushed after every point, so
you can watch one while it runs.

**⚠ Scans consume measurement numbers too.** They take their `A####` prefix from
the same `expt.measurement_num` that acquisitions use, but write under
`<cycle>/<experiment>/data/bluesky/`, not `data/`. So the highest number on disk
may well be in `data/bluesky/`; look in both before assuming where the counter
has got to.
