# Devices

[← index](../README.md)

How devices get into the session, and how to add or remove one.

## How devices are imported

Nothing in the Python code creates a device. Devices are **data**: every one is
an entry in a YAML file, instantiated at startup by a loader.

```
  configs/devices.yml            motors, slits, shutters, filters, QNW, transfocators, psic …
  configs/ad_devices.yml         area detectors (eiger4M, lambda2M, rigaku3M)
  configs/devices_aps_only.yml   things that only answer on the APS subnet (aps machine status)
            │
            ▼
  registry.py : safe_make_devices()
            │   for each YAML entry, in isolation:
            │     1. import the class named by the dotted path
            │     2. build it            → skip on failure
            │     3. wait_for_connection → skip if offline
            │     4. register it
            ▼
  oregistry["eiger4M"]        and also a bare name at the prompt:  eiger4M
```

Each top-level key in the YAML is a **dotted path to a class or factory**; each
list item under it is the keyword arguments for one instance:

```yaml
id8_common.devices.fast_shutter.FastShutter:
- name: shutter_8ide            # becomes the registry key AND the prompt name
  prefix: "8ideSoft:fastshutter:"

id8_common.devices.hv_motors.HV_Motors:
- name: bd6a                    # the kwargs are the class's own signature:
  prefix: "8ideSoft:CR8-E2:"    #   prefix + PV *suffixes*, not whole PVs
  pv_h: m9
  pv_v: m10
```

The loader does `import_module(...)` → `getattr(...)` → `target(**kwargs)`, so a
plain device class and a factory function are handled identically. That is why
`ad_devices.yml` can use `apstools.devices.area_detector_factory.ad_creator`
without the loader knowing anything about area detectors.

**Consequence worth internalising:** a device you never see mentioned in
`startup_ophyd.py` can still be fully live. `filter_8ide` is not named anywhere
in the startup script — it is loaded because it is in `devices.yml`. The three
area detectors and `psic` are the exceptions, and the reason is below.

`registry.py` is also the *lookup* side, with one definition of each helper:
`oregistry` itself, `get_connected_device("eiger4M")`, and
`get_ophyd_object("huber.x")`. The last one resolves its first segment through
`get_connected_device` and then walks plain attributes, which is how a device
named as a string in YAML — `inner_motor: huber.x` in `sample_info.yaml`, a
dual protocol's `motors:` block — becomes an object. `plans/set/select_device.py`
is the holdout: it resolves `device_position.yaml`'s dotted paths with its own
private `_resolve()`, which does not check `.connected`.

## Adding a device

1. Make sure the device class exists in [`src/id8_common/devices/`](../src/id8_common/devices/).
   Device classes are shared by 8-ID-E and 8-ID-I and belong there.
2. Add an entry under that class's dotted path in `configs/devices.yml` (or
   `devices_aps_only.yml` if it only exists on the APS subnet, or
   `ad_devices.yml` for an area detector). Reuse an existing top-level key if
   the class is already listed.
3. Restart the session. **No Python change is needed** — the loader is
   data-driven.
4. If a plan needs it, resolve it **inside the function** with
   `get_connected_device("my_device")` from `id8_common.registry`: a device that
   was skipped at startup gives a `KeyError` naming the banner, and one that
   registered but has since dropped gives a `RuntimeError`. Use
   `get_ophyd_object("my_device.x")` for the dotted form.
   Older modules (`plans/align/scan_8id.py`, `utils/nexus_runtime.py`) still bind
   devices at module scope with `oregistry.get("my_device")` — `.get()`, not
   `[...]`, so a missing device leaves `None` rather than aborting the import.
   That defers the failure to an `AttributeError` on `None` far from its cause;
   prefer the per-call form in new code. The exception is a device present in
   every session — `softglue` is still bound at module scope on purpose, because
   a per-call check buys nothing there.

### Area detectors need two extra lines

This is one of only two exceptions to "YAML only" — the other is `psic`, whose
hklpy2 wiring (`configure_hklpy2`) is a post-construction step for the same
reason. `ad_creator` builds the detector object, but two further steps cannot be
expressed in YAML and must be added to `startup_ophyd.py` (and `startup.py`):

```python
if "myDet" in oregistry:
    ad_setup(oregistry["myDet"], iconfig)          # plugin config + HDF5 priming
if "myDet" in oregistry:
    stream_rois(oregistry["myDet"], stats_nums=(1,))
```

**⚠ `stream_rois` is not copy-paste.** It defaults to `stats_nums=(1, 2, 3)`, but
`ad_creator` only builds the plugins listed in the YAML. eiger4M and lambda2M
declare `stats1`–`stats4`; **rigaku3M declares `stats1` only**, so it is called
with `stats_nums=(1,)`. A verbatim copy raises `AttributeError` on `stats2` and
aborts startup.

**⚠ Priming fires a real exposure.** `ad_setup` runs `AD_prime_plugin2` when
`ALLOW_AREA_DETECTOR_WARMUP` is set: `image_mode=Single`, `trigger_mode=0`,
`acquire=1`, 2 s, restore. rigaku3M is deliberately exempted (it is passed an
iconfig copy with the flag off) because its `cam.data_type` and `hdf1.data_type`
differ permanently, so the "already primed" check never passes and it would fire
an exposure on *every* startup.

Also update `read_path_template` / `write_path_template` in `ad_devices.yml`
each run cycle — they are real paths under `/gdata/dm/8ID/8IDE/<cycle>/`.

## Removing a device

Comment out or delete its list entry, then restart. If that was the last entry
under a class key, remove the key too.

There is no runtime unregister — `DeviceRegistry` has `register` and `clear` but
no `remove`. Taking a device out of a running session is not supported.

**Before removing, grep for it — and not only for `oregistry`.** A device name
reaches the code in four shapes, and only the first shows up in an `oregistry`
grep:

* bound at import time in the older plan modules (`scan_8id.py`,
  `nexus_runtime.py`), which then hold `None` and fail at the call site rather
  than at startup;
* passed as a bare string to `get_connected_device()` / `get_ophyd_object()`;
* as a string in a Python table — a mode's `required_devices` /
  `hardware_device` in the `*_MODES` dicts in `plans/acquire/*_modes.py`;
* as a string in YAML — a sample's `inner_motor` / `outer_motor` in
  `sample_info.yaml`, and the `device:` / `registers:` / `valve:` dotted paths
  in `plans/set/device_position.yaml`, which `select_device()` moves and writes.

So grep for the name itself:

```bash
grep -rn 'my_device' src/id8_common src/user_plans
```

A device that is simply *not connected* does not need removing at all — leave
the YAML entry and let the loader skip it. That is the design.

## When a device does not connect

`safe_make_devices` isolates every entry:

* **Import the class** — a bad dotted path skips that whole key, with a warning.
* **Build the instance** in its own `try/except` — a bad prefix or unknown kwarg
  skips only that device. (Stock `guarneri` builds in a plain loop with no
  per-entry guard, so one bad entry starves every device after it in the file.)
* **Check connectivity** with `wait_for_connection(timeout=20)` through a bounded
  4-worker thread pool, so many detectors do not all open Channel Access at once.

Only devices that pass are registered. The rest are named in the red banner.

**What this does *not* catch:** a device whose top-level signals connect but
which has a dead **lazy** PV — an uncommon area-detector plugin field, say. That
registers normally and only fails when a plan touches that specific PV. This is
intentional; an earlier version force-connected every lazy component and flooded
Channel Access.

Note also that connecting is not the same as being healthy. A device can connect
and still report stale or wrong values.
