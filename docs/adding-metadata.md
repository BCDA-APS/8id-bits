# Adding a field to the NeXus metadata file

Every measurement writes `<name>_metadata.hdf` next to the data. This page is
how to get a new value into that file — typically a new EPICS device you have
just added to the beamline.

## ⚠ First, which writer is actually running

**The default writer is this repository's `utils/nexus_utils.py`.** Miaoqi Chu's
`nexus_xpcs_aps` package is wired up but **off** unless you deliberately set an
environment variable:

```bash
ID8_NEXUS_WRITER=mc    # not set in any launcher -- his writer only runs if you export this
```

So unless you exported that yourself, the instructions below are the ones you
want. Note also that even under `ID8_NEXUS_WRITER=mc`, his module imports
`create_runtime_metadata_dict` from ours — the device readings all still come
from `nexus_utils.py`, and only the final file-writing step differs.

## The 30-second version

```
  1. THE DEVICE EXISTS IN OPHYD
     configs/devices.yml   ──►  oregistry["my_device"]
     one YAML stanza                  can you read it at the prompt?

                 │
                 ▼
  2. THE FIELD EXISTS IN THE SCHEMA          utils/xpcs_schema.py
     "my_field": {                           the template: what the file
         "type": "NX_FLOAT",                 CAN contain, and its units
         "units": "NX_TEMPERATURE",          and description
         "description": "...",
         "data": None,          ← placeholder, overwritten at write time
     }
                 │
                 ▼
  3. THE VALUE IS READ AT WRITE TIME         utils/nexus_utils.py
     my_device = oregistry.get("my_device")  module level, with the others
     ...
     runtime_updates = {                     create_runtime_metadata_dict()
         "/entry/.../my_field":              path ──► live value
             my_device.readback_ch1.get(),
     }
                 │
                 ▼
  4. THE FILE                                <name>_metadata.hdf
     /entry/.../my_field   = 297.3
        attrs: NX_Class="NX_FLOAT"  unit="K"  description="..."
```

**Both step 2 and step 3 are required.**

A path in step 3 with no node in step 2 raises `KeyError` at write time. **⚠ That
error is swallowed** — `det_acq_series()` wraps the whole measurement in a
blanket `except Exception`, so a typo'd path prints one line and the run carries
on without metadata. Read the session output after adding a field; silence is
not success.

A node in step 2 with no entry in step 3 is written with **whatever literal the
schema declares** — not an empty dataset. All 127 nodes in the real schema carry
a literal placeholder (`0`, `""`, `"N/A"` …), so a field you add and never
populate reports that placeholder as if it were a measurement. Only
`"data": None` gives a zero-length dataset, and nothing in the schema uses it.

There is a **third** source the diagram leaves out: `utils/default_metadata.py`.
`create_runtime_metadata_dict()` starts from a copy of that dict and
`runtime_updates` overrides it, so a path may already have a value from there.

## Worked example: a new temperature sensor

Say you have added an EPICS temperature readout and want it in every file.

### Step 1 — make Ophyd see it

Add a stanza to `src/id8_common/configs/devices.yml`:

```yaml
# There is already a block for this class in devices.yml -- append to it rather
# than starting a second one with the same key, which would silently win.
id8_common.devices.lakeshore.Lakeshore:
- name: lakeshore1
  prefix: "8ideSoft:LS336:1:"
- name: lakeshore2
  prefix: "8ideSoft:LS336:2:"
- name: my_sensor                 # <-- the new one
  prefix: "8ideSoft:LS336:3:"
```

The top-level key is the dotted path to the class; each list item is the kwargs
it is built with. The `prefix` is the EPICS PV prefix, and the class supplies the
suffixes (`Lakeshore` declares `readback_ch1 = Component(EpicsSignalRO, "IN1")`,
so `my_sensor.readback_ch1` is `8ideSoft:LS336:3:IN1`).

Restart the session and check it loaded:

```python
"my_sensor" in oregistry            # True
my_sensor.readback_ch1.get()        # 297.34
```

If it is missing, look for it in the `*** Devices not online: … ***` banner —
see [Devices](devices.md).

### Step 2 — add the field to the schema

In `src/id8_common/utils/xpcs_schema.py`, find the group you want it under and
add a node. To put it at `/entry/sample/temperature_set`:

```python
"sample": {
    "type": "NXsample",
    "required": False,        # as it really is in the file -- do not "fix" this
    ...
    "temperature_set": {
        "type": "NX_FLOAT",
        "required": True,
        "units": "NX_TEMPERATURE",
        "description": "Sample temperature setpoint",
        "data": None,
    },
},
```

The keys mean:

| key | effect |
|---|---|
| `"data"` present | written as a **dataset**. Absent → written as a **group** |
| `"data": None` | the field exists in the file with no value (zero-length dataset) |
| `"type"` | goes in as the `NX_Class` attribute |
| `"units"` | a *category*, not a unit string — see the table below |
| `"description"` | goes in as the `description` attribute |
| `"required": False` | would be skipped if the writer were called with `ignore=True` — nothing in the repo ever does, so this has no effect today |

`"units"` is looked up in `default_units_keymap` in `nexus_utils.py`, which maps
the category to the string actually written:

```
NX_LENGTH → m      NX_TIME → s        NX_ENERGY → keV     NX_ANGLE → degree
NX_TEMPERATURE → K NX_CURRENT → mA    NX_PER_LENGTH → 1/Å NX_VOLTAGE → V
NX_COUNT → one     NX_DIMENSIONLESS → dim.less            NX_ANY → any
```

An unrecognised category silently becomes `"any"`. If your unit is not in that
map, add it there rather than inventing a category.

### Step 3 — read the live value

In `create_runtime_metadata_dict()` in `src/id8_common/utils/nexus_utils.py`,
add a line to the `runtime_updates` dict:

This file resolves its devices at **module level**, alongside `filter_8ide`,
`lakeshore1`, `mono` and the rest near the top:

```python
my_sensor = oregistry.get("my_sensor")
```

then reads it in `runtime_updates`:

```python
runtime_updates = {
    ...
    "/entry/sample/temperature_set": my_sensor.readback_ch1.get(),
}
```

**⚠ Follow that convention here, even though it is not the one used elsewhere.**
Plans use `get_connected_device("…")` *inside* the function so a missing device
names itself. `nexus_utils.py` does not, and does not even import it — writing
that call here raises `NameError` the first time a measurement writes metadata.
The cost of the module-level style is that a device which failed to connect at
startup binds `None`, and you find out as `AttributeError: 'NoneType' object has
no attribute 'readback_ch1'` while writing the file, after the data is on disk.

The value must be something h5py can store: a number, a string, or a numpy
array. An Ophyd `Signal` object is not — always call `.get()`.

### Step 4 — check it

Run one short measurement, then:

With `use_subfolder: no` (the current setting) the file is at
`…/data/<name>/<name>_metadata.hdf`; with `yes` there is an extra
`<file_header>/` level above it.

```bash
h5dump -n /gdata/dm/8ID/8IDE/<cycle>/<expt>/data/<name>/<name>_metadata.hdf | grep temperature
python -c "import h5py; f=h5py.File('<name>_metadata.hdf'); print(f['/entry/sample/temperature_set'][()])"
```

## Removing a field

Delete it from **both** places: the node in `xpcs_schema.py` and the line in
`runtime_updates`. Removing only the schema node leaves a `KeyError` at write
time; removing only the runtime line leaves an empty dataset in every file.

Removing a field changes what downstream analysis sees. Check that `boost_corr`
and any plotting scripts do not read it before deleting.

## Per-measurement values, without touching the schema

Both writers accept `additional_metadata`, a dict of the same
`path → value` shape, merged after everything derived from
`device_position.yaml`. The dual-detector path uses it for a leg that declares a
`geometry:` block — the exception, for a detector remounted somewhere
`device_position.yaml` does not describe, not something every leg does:

```python
create_nexus_format_metadata(fname, det=det,
                             additional_metadata={"/entry/instrument/detector_1/distance": 8.0})
```

The path still has to exist in the schema — this overrides a value, it does not
create a field.

## ⚠ Traps

* **Never hand `xpcs_schema` itself to the writer.** `create_nexus_entry()`
  `pop()`s five keys — `required`/`data`/`type`/`units`/`description` — out of
  every node as it writes,
  so a shallow copy guts the module-level template. `create_nexus_format_metadata()`
  takes a `deepcopy` for exactly this reason. Before that fix, the first file in
  a session had all 143 objects carrying their attributes and every later one had
  **none** of them, with six datasets degraded to empty groups — and the loss
  propagated into the analysed `*_results.hdf`.
* **`/entry/start_time` and `/entry/end_time` are both written at file-write
  time**, i.e. after the acquisition. Neither is the acquisition start.
* The class attribute is written as `NX_Class`; the NeXus standard spells it
  `NX_class`. Left as-is deliberately — every file ever written at 8-ID uses this
  spelling and nothing reads it back.

## If you switch to Miaoqi's writer

`ID8_NEXUS_WRITER=mc` swaps only the final step, and **only at one call site** —
the normal path in `det_acq_series()`. The abort path (`cleanup_acquisition()`)
and the whole dual-detector path ignore the variable and always use
`nexus_utils.py`, so with it set a run can produce files from both writers.
Steps 1 and 3 are unchanged;
step 2 moves to `utils/xpcs_schema_mc.py`, and the clone of `nexus_xpcs_aps`
must be on `PYTHONPATH`. If you add a field to `xpcs_schema.py` and not to
`xpcs_schema_mc.py`, the two writers produce different files.

## Related

* [Using `expt`](using-expt.md) — reading the values you may want to record
* [Devices](devices.md) — adding the EPICS device in the first place
* [Data Management](reference/data-management.md) — what happens to the file afterwards
