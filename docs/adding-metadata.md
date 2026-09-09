[← index](../README.md)

# Adding a field to the NeXus metadata file

Every measurement writes `<name>_metadata.hdf` next to the data. This page is how
to get a new value into it — typically an EPICS device you have just added.

**Since 2026-09-08 the writer is Miaoqi Chu's `nexus_xpcs_aps`.** Our own
`nexus_utils.py`, `xpcs_schema.py` and `default_metadata.py` are retired to
`utils/Archive/` as `.txt` — kept for reference, not importable. There is no
env-var switch any more; there is one writer.

## The 30-second version

Three files, and knowing which one to edit is most of the job:

```
   configs/devices.yml  ──►  oregistry  ──►  live EPICS values
                                                    │
   ┌────────────────────────────────────────────────▼──────────────────────────┐
   │  utils/nexus_runtime.py            WHERE THE NUMBERS COME FROM      OURS  │
   │                                                                           │
   │    runtime_updates = {                                                    │
   │      "/entry/sample/huber_nu": huber.nu.position,                         │
   │      ...                                                                  │
   │    }                          one line per field, read at write time      │
   └────────────────────────────────────────────────┬──────────────────────────┘
                                                    │  {path: value}
   ┌────────────────────────────────────────────────▼──────────────────────────┐
   │  utils/xpcs_schema.py           WHAT FIELDS EXIST                OURS, │
   │                                                              CALLING       │
   │    "sl4": _tag(make_slits(4), ...)      ← upstream factory, our placement      │
   │    instrument = {...}  sample = {...}   ← 8-ID's composition             │
   └────────────────────────────────────────────────┬──────────────────────────┘
                                                    │  schema + values
   ┌────────────────────────────────────────────────▼──────────────────────────┐
   │  utils/nexus_writer.py             THE FRONT DOOR                   OURS  │
   │    create_nexus_format_metadata(filename, det)                            │
   │    + workarounds for two upstream bugs, + the missing unit categories     │
   └────────────────────────────────────────────────┬──────────────────────────┘
                                                    ▼
                      nexus_xpcs_aps.core.utils     WRITES THE HDF5   UPSTREAM
```

**Both a schema entry and a runtime line are required.** A runtime path with no
schema node raises `KeyError` at write time — and `det_acq_series()` swallows it,
so the measurement finishes and writes *no metadata file at all*. Read the
session output after adding a field; silence is not success.

## Which case are you in?

**Ask: does an upstream factory already produce the group you want?**

| factory | signature | gives you |
|---|---|---|
| `make_slits` | `(index, description=None)` | 4 fields, `NXslit` |
| `make_attenuator` | `(index)` | 2 fields, `NXattenuator` |
| `make_undulator` | `(index)` | 3 fields, `NXinsertion_device` |
| `make_detector` | `(index, name=)` | 19 fields, `NXdetector` |
| `make_diffractometer` | `(name)` | 6 fields, `NXpositioner` |
| `make_sample` | `(qnw=, rheometer=, huber_stage=, lakeshore=, keithley=, bk_pid=)` | 17–31 fields, `NXsample` |
| `make_entry` | `(beamline, instrument, sample, user)` | the top-level wrapper |

The first five take an **index** and can be called again for a new instance.
`make_sample` takes **flags** and returns a fixed menu. There is no `make_stage`,
no `make_motor`, and no arbitrary-leaf helper.

---

## Case 1 — a device the upstream factories already model

One line in `xpcs_schema.py`, inside the `instrument` dict. Upstream untouched:

```python
"sl9": _tag(make_slits(9, description="Slits 9"), "/entry/instrument/sl9", "upstream:make_slits"),
```

That is a complete `NXslit` group. The index need not be numeric — the schema
already does `make_slits("wb", …)` and `make_slits("mono", …)`.

Then the runtime lines in `nexus_runtime.py`:

```python
"/entry/instrument/sl9/horizontal_gap":    sl9.h_size.get(),
"/entry/instrument/sl9/horizontal_center": sl9.h_center.get(),
"/entry/instrument/sl9/vertical_gap":      sl9.v_size.get(),
"/entry/instrument/sl9/vertical_center":   sl9.v_center.get(),
```

…with the device bound at module scope near the top, **always `.get()`**:

```python
sl9 = oregistry.get("sl9")
```

**⚠ Never `oregistry["sl9"]` here.** `nexus_runtime` is imported by `nexus_writer`
→ `ad_acq` → startup. A bracket lookup raises `KeyError` at import and the session
will not start — one dead IOC costing you the whole beamline instead of one field.
`.get()` binds `None`, and the failure lands at write time inside
`det_acq_series()`'s except-block.

---

## Case 2 — a field he does not model (e.g. a pressure controller)

`make_sample()` has six fixed flags and none is pressure, so you patch the
composed dict. Same pattern as the existing `_rename()` and `_redescribe()` calls
in `xpcs_schema.py`:

```python
_sample["pressure"] = {
    "type": "NX_FLOAT", "required": False, "units": "NX_PRESSURE",
    "description": "Sample pressure readback", "data": 0.0,
}
_sample["pressure_set"] = {
    "type": "NX_FLOAT", "required": False, "units": "NX_PRESSURE",
    "description": "Sample pressure setpoint", "data": 0.0,
}
```

Then in `nexus_runtime.py`, with `pcd1 = oregistry.get("pcd1")` at module scope:

```python
"/entry/sample/pressure":     pcd1.pressure.get(),
"/entry/sample/pressure_set": pcd1.setpoint_rbv.get(),
```

`NX_PRESSURE` already resolves to `Pa` — see the units note below.

## Removing a field

Delete it from **both** places: the schema node in `xpcs_schema.py` and the
line in `nexus_runtime.py`. Schema only → `KeyError` at write time (swallowed).
Runtime only → the schema's literal placeholder is written as if it were a
measurement. Check that `boost_corr` and any plotting scripts do not read the
field before removing it.

## Units

The category string in a schema node (`"units": "NX_LENGTH"`) is looked up in
`nexus_xpcs_aps.core.utils.default_units_keymap` to get what is written. **The
map has 10 entries and an unknown category silently becomes `"any"`** — no
warning anywhere.

`nexus_writer.EXTRA_UNITS` patches in the three he lacks, via `setdefault` so an
upstream fix wins automatically:

```
NX_VOLTAGE → V     NX_FREQUENCY → Hz     NX_PRESSURE → Pa
```

**⚠ That fixes the map, not the schema.** Several upstream leaves *declare*
`NX_ANY` — `keysight_freq` and `keysight_amp` among them — so they write `"any"`
no matter what the map contains. Fixing those means changing the declared
category in `xpcs_schema.py`, or upstream. Low priority: a field's unit is
fixed by its device and recoverable by looking at the device.

## Per-measurement overrides

`create_nexus_format_metadata(..., additional_metadata={path: value})` merges
last. The dual path uses it for a leg with a `geometry:` block. The path must
already exist in the schema — this overrides a value, it does not create a field.

## Local patches currently carried

**None.** All five landed upstream in
[PR #1](https://github.com/AZjk/nexus_xpcs_aps/pull/1) (merged 2026-09-08) and
were deleted from this repo in `8b83c67`: the `make_sample()` leaf aliasing, the
`id()`-keyed plan cache, the invalid-JSON `workflow_kwargs` default, the three
missing unit categories, and the `beam_center_position_x/y` description.

The `deepcopy` in `nexus_writer.py` is **not** a leftover patch and must stay:
`update_schema_at_runtime()` assigns `node["data"] = value` in place by design,
so any caller has to protect its own module-level template. See
[How the NeXus file is written](reference/nexus-writer.md#local-patches-we-carry).

## When this needs Miaoqi rather than you

Short version: **if a value is wrong, it is ours; if a field cannot exist, it is
upstream's.** The full decision table is in
[How the NeXus file is written](reference/nexus-writer.md#-when-to-contact-miaoqi).

The one case that comes up in practice: a new **sample-environment** device.
`make_sample()` takes six fixed flags, so unlike `make_slits(9)` you cannot add a
seventh from outside the upstream package. Patch the composed dict as in Case 2 above to
keep running, and open a PR against
[AZjk/nexus_xpcs_aps](https://github.com/AZjk/nexus_xpcs_aps) for the real fix.

## Related

* [How the NeXus file is written](reference/nexus-writer.md) — how we got here, and the comparison
* [Devices](devices.md) — adding the EPICS device in the first place
