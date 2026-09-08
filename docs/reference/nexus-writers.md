[← index](../../README.md) · [Adding metadata fields](../adding-metadata.md)

# Two NeXus writers: ours and Miaoqi's

Why the beamline still writes metadata with `utils/nexus_utils.py` when
`nexus_xpcs_aps` exists, and what switching would actually involve.

## The 30-second version

```
  TODAY (default, nothing exported)

    devices.yml ─► oregistry ─► create_runtime_metadata_dict()   ← OURS
                                        │  path ──► live value
                                        ▼
                              xpcs_schema.py                     ← OURS
                              993 lines of literal dicts
                                        │
                                        ▼
                              create_nexus_entry()               ← OURS
                                        │
                                        ▼
                              <name>_metadata.hdf

  WITH  ID8_NEXUS_WRITER=mc  (test path, one call site only)

    devices.yml ─► oregistry ─► create_runtime_metadata_dict()   ← STILL OURS
                                        │
                                        ▼
                              xpcs_schema_mc.py                  ← 244 lines,
                              composed from his factories          calling HIS code
                                        │
                                        ▼
                       nexus_xpcs_aps.core.utils                 ← HIS
                       .create_nexus_format_metadata()
                                        │
                                        ▼
                              <name>_metadata.hdf
```

**Only the bottom two boxes change.** The beamline layer — which device is read,
what value goes at which path — is ours in both cases and is not something his
package models.

## Side by side

| | ours (`nexus_utils.py`) | his (`nexus_xpcs_aps.core`) |
|---|---|---|
| **schema leaves** | 125, as 993 lines of literal dicts | all 125 reachable, from 244 lines of factory calls |
| **class attribute** | `NX_Class` | `NX_class` ✅ **correct per the NeXus standard** |
| **units attribute** | `unit` | `units` ✅ **correct per the standard** |
| **units keymap** | 11 categories, incl. `NX_VOLTAGE` | 10 categories, **no** `NX_VOLTAGE` |
| **unit strings lost as `any`** | 1 (`keysight_freq`, needs `NX_FREQUENCY`) | 20 |
| **field descriptions** | ours on 41 leaves | his on 41; 2 are typo fixes, ~13 are losses |
| **new slit/attenuator/undulator/detector** | ~30 lines of literal dict | **1 line** — `make_slits(9)` |
| **new sample-environment device** | ~10 lines | ~10 lines; `make_sample` is 6 fixed flags, not extensible |
| **generic / catch-all class** | n/a | **none** — no `make_stage`, no arbitrary-leaf helper |
| **runtime layer** (signal → path) | ours, 269 of 420 lines | **not modelled at all**; ours is carried over unchanged |
| **who maintains the schema** | us | him; fixes arrive on a `git pull` |
| **installed?** | yes, it is the repo | **no** — clone + `PYTHONPATH` only |
| **known silent bugs** | none found | 2: aliased `make_sample` leaves, `id()`-keyed plan cache |

**On the attribute spelling, the standard is unambiguous.** The NeXus manual
gives `@NX_class` and `@units`. Ours is wrong on both counts and has been since
the `legacy/` code — which is why every file in the 8-ID archive carries the
wrong spelling, and why fixing it is an archive-consistency decision rather than
a typo fix.

## Status right now

| | |
|---|---|
| installed in `8id_bits`? | **No.** `import nexus_xpcs_aps` → `ModuleNotFoundError` |
| where it lives | `~/Documents/Miaoqi/nexus_xpcs_aps_95ab368/src`, reached by `PYTHONPATH` |
| declared dependencies | `h5py`, `numpy` — both already in the environment |
| how it is selected | `export ID8_NEXUS_WRITER=mc`, set by no launcher |
| which call sites honour it | 1 of 3 — the normal path in `det_acq_series()` only |
| verified equivalent? | Yes: 125/125 leaves path-for-path identical, 0 hand-written |

**⚠ Two clones sit side by side.** `~/Documents/Miaoqi/nexus_xpcs_aps` is `main`
at `73d0be4` (Aug 2025) and stale; `nexus_xpcs_aps_95ab368` is the `mc_refact`
commit the adapter was written against. Point `PYTHONPATH` at the second.

## Adding a device: which case are you in?

Ask one question — **does one of his factories already produce the group you
want?** The complete list, with what each returns:

| factory | signature | gives you |
|---|---|---|
| `make_slits` | `(index, description=None)` | 4 fields, `NXslit` |
| `make_attenuator` | `(index)` | 2 fields, `NXattenuator` |
| `make_undulator` | `(index)` | 3 fields, `NXinsertion_device` |
| `make_detector` | `(index, name="Eiger4m")` | 19 fields, `NXdetector` |
| `make_diffractometer` | `(name)` | 6 fields, `NXpositioner` |
| `make_sample` | `(qnw=, rheometer=, huber_stage=, lakeshore=, keithley=, bk_pid=)` | 17–31 fields, `NXsample` |
| `make_entry` | `(beamline, instrument, sample, user)` | the top-level wrapper |

**There is no `make_stage` and no `make_motor`**, and no arbitrary-leaf helper.

Note the shape difference: the first five take an **index** and can be called
again for a new instance. `make_sample` takes **flags** and returns a fixed
menu — you cannot ask it for a device it does not already know about.

---

## Worked example 1 — another slit (the easy case)

**With his package.** One line in *our* `xpcs_schema_mc.py`, inside the
`instrument` dict. His code is untouched:

```python
"sl9": _tag(make_slits(9, description="Slits 9"), "/entry/instrument/sl9", "his:make_slits"),
```

That is a complete `NXslit` group with `horizontal_gap`, `horizontal_center`,
`vertical_gap`, `vertical_center`, correctly typed and described. The index need
not be numeric — the schema already does `make_slits("wb", …)` and
`make_slits("mono", …)`.

**With ours.** The same group, written out by hand in `xpcs_schema.py`:

```python
"sl9": {
    "type": "NXslit",
    "required": False,
    "description": "Slits 9",
    "horizontal_gap":    {"type": "NX_FLOAT", "units": "NX_LENGTH", "required": False,
                          "description": "Horizontal size of the slits",   "data": 1.0},
    "horizontal_center": {"type": "NX_FLOAT", "units": "NX_LENGTH", "required": False,
                          "description": "Horizontal center of the slits", "data": 0.0},
    "vertical_gap":      {"type": "NX_FLOAT", "units": "NX_LENGTH", "required": False,
                          "description": "Vertical size of the slits",     "data": 1.0},
    "vertical_center":   {"type": "NX_FLOAT", "units": "NX_LENGTH", "required": False,
                          "description": "Vertical center of the slits",   "data": 0.0},
},
```

**One line versus about thirty**, and his version cannot drift from the NeXus
class definition because he owns it.

Both then need the same four runtime lines in **our**
`create_runtime_metadata_dict()` — this half never changes hands:

```python
"/entry/instrument/sl9/horizontal_gap":    sl9.h_size.get(),
"/entry/instrument/sl9/horizontal_center": sl9.h_center.get(),
"/entry/instrument/sl9/vertical_gap":      sl9.v_size.get(),
"/entry/instrument/sl9/vertical_center":   sl9.v_center.get(),
```

---

## Worked example 2 — a pressure controller (the case you actually have)

A set pressure and a readback pressure: structurally a temperature controller.

**⚠ Neither writer has this today.** His `make_sample()` can produce 31 fields
across all six flags and **none of them is a pressure**. Our schema has no
pressure field either. The closest existing shapes are his
`qnw1_temperature` / `qnw1_temperature_set` and `bk_pid_VAL` / `bk_pid_RDBK` —
both setpoint+readback pairs, both hard-wired to their device.

And because `make_sample` is flags-based, **there is no `make_sample(pressure=True)`
to call.** This is the case where his package gives you nothing extra.

**With ours** — add the node to `xpcs_schema.py` under `"sample"`:

```python
"pressure": {
    "type": "NX_FLOAT", "required": False, "units": "NX_PRESSURE",
    "description": "Sample pressure readback", "data": 0.0,
},
"pressure_set": {
    "type": "NX_FLOAT", "required": False, "units": "NX_PRESSURE",
    "description": "Sample pressure setpoint", "data": 0.0,
},
```

**With his** — call `make_sample()` as now, then patch the returned dict in
`xpcs_schema_mc.py`, the same way `_rename()` already patches
`flightpath_swing_horizontal`:

```python
_sample = make_sample(qnw=False, rheometer=False, lakeshore=True)
_sample["pressure"]     = {...}    # the same two literal blocks as above
_sample["pressure_set"] = {...}
```

**Same work either way — about ten lines.** The difference is what happens next:
with his, `pressure` is a standard `NXsample` field, so it is a reasonable thing
to ask him to add as a `pressure=True` flag. Once he does, your local patch
deletes and everyone at the APS gets it. With ours, it stays yours forever.

Then, identically, the runtime lines in our `create_runtime_metadata_dict()`:

```python
"/entry/sample/pressure":     my_pressure.readback.get(),
"/entry/sample/pressure_set": my_pressure.setpoint.get(),
```

**⚠ One gotcha that bites either way.** `default_units_keymap` in
`nexus_utils.py` has no `NX_PRESSURE` entry, so the unit attribute silently
falls through to `"any"` instead of `"Pa"`. Add it:

```python
"NX_PRESSURE": "Pa",
```

Unrecognised unit categories do not raise — they just write `"any"` — so nothing
tells you this happened except reading the file.

---

## Summary

| your device | with ours | with his |
|---|---|---|
| another slit / attenuator / undulator / detector | ~30 lines of literal dict | **one line** |
| one more field in a group he models | 2 steps | **1 step**, and upstream fixes flow in |
| a pressure controller, a bare motor stage | ~10 lines | ~10 lines, **but can become his** |

## Field-by-field: what he covers, and what neither writer has

Measured by building both schemas and diffing the leaf paths, not by reading
prose. **Ours: 125 leaves. His, with every factory and every `make_sample` flag
on: 119.**

### He covers all 125 of our leaves

18 leaf *names* differ, but not one is a missing capability — all three groups
are renames or re-placements, and `xpcs_schema_mc.py` already resolves them,
which is how it reaches 125/125 with nothing hand-written:

| ours | his | difference |
|---|---|---|
| `detector_1/flightpath_swing` | `flightpath_swing_horizontal` | rename (already handled by `_rename()`) |
| `/entry/sample/huber_{chi,delta,eta,mu,nu,phi}` | `make_diffractometer("huber")` → `chi, delta, eta, mu, nu, phi` | same six angles, no `huber_` prefix, and he puts them under `instrument` as `NXpositioner` |
| `/entry/sample/keysight_*` (11) | `instrument/keysight_waveform_generator` → `amp, freq, func, phase, output, pulse_width, burst_{count,mode,state}, trigg_{edge,source}` | same eleven, no `keysight_` prefix, under `instrument` |

**⚠ The last two are re-placements, not renames.** Migrating means repointing 17
paths in our `create_runtime_metadata_dict()`, and any downstream reader that
looks for `/entry/sample/keysight_freq` will not find it at the new address.

### The detector's degrees of freedom are fully covered

`make_detector(1)` returns 19 fields, a superset of what we use:

```
position_x  position_y                      translations
rotation_x  rotation_y  rotation_z          angles (he has three; we use none)
flightpath_swing_horizontal / _vertical     the two swing angles
x_pixel_size  y_pixel_size                  pixel size
distance  beam_center_x/y  beam_center_position_x/y
count_time  frame_time  compression  detector_name  qmap_file
```

Nothing to add here for a detector — including pixel size, which you asked about.

### The devices you asked about

| | ours | his | verdict |
|---|---|---|---|
| `keithley` | 8 fields | 8 fields | **both** — commented out in our runtime only |
| `bk_pid` | 2 | 2 | **both** — same |
| `keysight` | 11 | 11 | **both**, at different paths (see above) |
| `pixel_size` | 2 | 2 | **both** |
| **`pressure`** | 0 | 0 | **neither** |
| **`fofb`** | 0 | 0 | **neither** — and see the trap below |
| **`xbpm`** | 0 | 0 | **neither**; the device is also commented out in `devices_aps_only.yml` |

So `keithley` and `bk_pid` are not gaps at all — the schema nodes exist on both
sides and only the runtime lines are commented out. `pressure`, `fofb` and
`xbpm` are the three genuine holes, and they are holes in *both* writers.

### ⚠ Two commented-out lines are booby-trapped

Of the 30 commented-out paths in `create_runtime_metadata_dict()`, **28 are safe
to uncomment** — the schema node exists, so they need only a working device.
**Two are not:**

```
/entry/instrument/incident_beam/fofb_s09_horizontal
/entry/instrument/incident_beam/fofb_s09_vertical
```

`incident_beam` has no `fofb_*` node — it holds `extent`,
`incident_beam_intensity`, `incident_energy`, `incident_energy_spread`,
`incident_polarization_type`, `ring_current`, `transmitted_beam_intensity`.
Uncommenting either line raises `KeyError` at write time, and
`det_acq_series()` swallows it — so the measurement finishes, prints one line,
and writes **no metadata file at all**. Add the schema node first.

### There is no catch-all class

`make_entry(beamline, instrument, sample, user)` is the top-level wrapper: it
takes the three already-composed sub-dicts and adds the entry-level fields
(`definition`, `schema_version`, `start_time`, `end_time`). It is not a home for
a device that fits nowhere else. The complete inventory of `core/` is eleven
instrument modules, six sample-environment modules, `schema.py`, `user.py` and
`utils.py` — and **no generic module, no arbitrary-leaf helper, no
`make_stage`**.

## What each costs

| | `nexus_utils.py` (ours) | `nexus_xpcs_aps` (his) |
|---|---|---|
| schema | 993 lines of literals, all ours to maintain | 244 lines of composition; the leaves are his |
| who fixes a NeXus-standard bug | us | him, and it arrives on a pull |
| adding a device that reuses an existing factory | ~30 lines of literal dict | one line, no upstream change |
| adding something he does not model (e.g. pressure) | a literal node, ours forever | a local patch that can later become his |
| dependency risk | none | none in practice: `h5py` + `numpy` |
| coupling to `apsbits` | none | none, **provided only `core/` is used** — his `deployment/id8_{e,i}` needs apsbits and is obsolete |
| behaviour if it breaks mid-beamtime | ours to fix | ours to work around |

## Before 2026-09-09: what is actually blocked

The plan of "uncomment everything in `nexus_utils.py`, plus a sample pressure
field" is two very different jobs.

### Uncommenting the metadata lines — **do not do this before the run**

It is not one edit, it is a chain of four, and the third one can stop the
session from starting at all:

1. The **device** is commented out in `configs/devices.yml`.
2. The **module-level bind** is commented out in `nexus_utils.py`.
3. The bind uses `oregistry["name"]`, **not** `oregistry.get("name")` like the
   live ones above it. Bracket lookup **raises `KeyError` at import**, and
   `nexus_utils` is imported by `ad_acq`, which is imported by startup. A device
   whose IOC is down therefore does not degrade — **it stops the session
   booting.**
4. Only then the runtime line itself.

And the hardware is not there. Every PV behind these fields timed out from
`pearl`, which can see the beamline:

```
S08ID:USID:Gap.VAL                  not found      undulator_upstream
S08ID:DSID:Gap.VAL                  not found      undulator_downstream
8idKeithley2600:SMU:A:SrcLevelV_AO  not found      keithley_chA / chB
8idiSoft:BPMpid.VAL                 not found      bk_pid
```

So uncommenting the 28 lines would, today, either write nothing useful or break
startup. The work is real but belongs after the run, in this order: bring the
IOC up → uncomment the `devices.yml` entry → confirm it appears in `oregistry`
→ change the bind to `oregistry.get()` → uncomment the runtime line.

**⚠ Whatever else happens, change those binds from `oregistry["x"]` to
`oregistry.get("x")`** so a dead IOC costs you a field rather than the session.

### The sample pressure field — safe, and worth doing

Additive and self-contained: two literal nodes in `xpcs_schema.py`, two runtime
lines, and `"NX_PRESSURE": "Pa"` in `default_units_keymap`. It touches no
existing field, so nothing already working can regress. The only prerequisite is
the device itself in `devices.yml` and its IOC up. Needed 09/10, so there is
time to add it on the 9th once the controller is on.

### Switching writers — not before the run

Changing the attribute spelling on every object of every file, three days before
a beamtime, with no check yet on what reads those files, is the definition of a
bad time. The switch is right; the date is wrong.

## Can we migrate with the code that exists today?

**Yes, technically — but it is not a flip-a-switch job, and the schema is the
easy half.** He covers 125/125 of our leaves. What makes it work is deciding
what happens to the *contents* of every file.

### What stays ours no matter what

His package models the schema and the file writer. It does **not** model the
beamline layer: which EPICS signal feeds which path. That is
`create_runtime_metadata_dict()` plus its helpers and device binds — **269 of
`nexus_utils.py`'s 420 lines** — and it is carried over unchanged. Migrating
replaces the schema and the writer, not the part that knows about your hardware.

### ⚠ Every file changes, on every object

This is the decision to make first, because these files are copied into
`*_results.hdf`:

| change | scope |
|---|---|
| `NX_Class` → `NX_class` | **143 of 143 objects** (his spelling is the correct NeXus one) |
| `unit` → `units` | **97 of 97** unit-bearing datasets |
| unit string lost, becomes `any` | **20 datasets** — his keymap has 10 entries to our 11 and lacks `NX_VOLTAGE`, so the four keithley `*V` fields and `keysight_amp` degrade on top of 16 others |
| description text | 41 leaves |
| actual data value | 1 — `/entry/instrument/datamanagement/workflow_kwargs`, whose default is invalid JSON on his side; the other 11 diffs are masked by `default_metadata.py` |

**⚠ Do not take his descriptions wholesale.** Only 2 of the 41 are typo fixes.
At least 13 are losses: his `beam_center_position_x/y` says "position of beam
center" but we store the *detector translation preset* there, so his text would
mislead; `bk_pid_RDBK/VAL` lose the word "temperature" at the same moment their
unit degrades from `K` to `any`; and all 8 keithley strings become bare labels
that drop the current/voltage semantics.

### ⚠ Two upstream bugs our adapter is already working around

Both are in his code, and both are silent:

1. **`make_sample()` returns leaves aliased to module-level singletons** — all 17
   shared between two calls, and shared with his own `core.schema.xpcs_schema`.
   His documented `xpcs_schema.copy()` pattern therefore corrupts the template
   permanently. `nexus_utils_mc.py` deepcopies instead.
2. **`_compiled_plans` is cached on `id(schema)`.** Handing it a fresh deepcopy
   each call means a freed address can be reused and return a stale plan — it
   collided 30 times in 200 cycles under test and wrote a file with a leaf
   silently missing. `nexus_utils_mc.py` clears the cache every call, which
   negates the cache's only purpose.

Neither is a reason not to migrate; both are reasons the adapter must stay, or
be fixed upstream.

### Environment hazards

* Not installed in any beamline environment. `import nexus_xpcs_aps` fails in
  `8id_bits`.
* **An editable install already exists in `e2507_timepix`, and it points at the
  wrong clone** — `.../Miaoqi/nexus_xpcs_aps/src`, which is `main` at `73d0be4`
  (Aug 2025) *and has an uncommitted modification*. So "just pip install -e it"
  has already been done once against the stale tree. Install
  `nexus_xpcs_aps_95ab368`, and pin it.
* The source of truth would move to a directory in a user's home with no pinned
  version, where a `git pull` silently changes what the beamline writes.

### Call sites

Three write a metadata file, and **only one honours `ID8_NEXUS_WRITER`**:

| site | honours the flag? | note |
|---|---|---|
| `ad_acq.det_acq_series()` normal path | yes | the only one |
| `ad_acq.cleanup_acquisition()` abort path | no | **this is the safety net** — it rewrites a failed file with our known-good writer. Switching all three removes it |
| `dual_acq_eiger4m_rigaku3m` | no | per-leg overrides make the deepcopy discipline mandatory here |

### One thing that is *not* nearly free

`make_detector(2)` really does return a second 19-field detector group. But
`detector_1` is hardcoded **16 times in `nexus_utils.py`, 18 in
`default_metadata.py`, and 13 in dual_acq's `OVERRIDE_PATHS`** — all in the
beamline layer he does not model. The factory call is the cheapest part.

### The order of work

1. Grep every consumer of `*_metadata.hdf` and `*_results.hdf` — boost_corr,
   Miaoqi's readers, user scripts — for `NX_Class` and `unit`. **Nothing else
   should start until this is answered.**
2. `pip install -e` the **95ab368** clone into `8id_bits`, pinned. Delete the
   stale clone and the `e2507_timepix` editable install that points at it.
3. Add `NX_VOLTAGE` (and `NX_PRESSURE`) to whichever keymap ends up in use;
   reconcile the other 15 unit categories his schema leaves as `NX_ANY`.
4. Decide the 41 descriptions field by field. Do not bulk-accept.
5. Collapse the three call sites behind one `write_nexus_metadata()` shim, so
   the writer is chosen in one place rather than three.
6. Keep the deepcopy and the `_compiled_plans.clear()`, or fix them upstream.
7. Re-run the comparison on one measurement per detector mode.

Steps 1 and 4 are the ones that need a person, not a script.

## Why we have not switched

Honestly: **it is a test that was left switched off, not a considered rejection.**
The comparison was run, came out clean at the path level, and then the
DM/analysis investigation took over.

The original two objections have not aged equally. *"Installing risks disturbing
`8id_bits`"* is largely gone — the declared dependencies are `h5py` and `numpy`,
both already present. *"`deployment/` is unusable"* is still true and still
irrelevant: it needs `apsbits`, which this instrument is leaving, and we import
only `core/`.

What replaced them is the list above: the switch changes an attribute on every
object of every file, and nobody has yet checked what reads those files.

## Loose end

`ad_acq.py` names a launcher `run_mc_writer_test.sh` that does not exist
anywhere. To run the test path today, set `PYTHONPATH` yourself:

```bash
PYTHONPATH=~/Documents/Miaoqi/nexus_xpcs_aps_95ab368/src \
ID8_NEXUS_WRITER=mc ~/bin/start_ophyd.sh
```

## Related

* [Adding metadata fields](../adding-metadata.md) — the step-by-step for the writer in use today
* [Data Management](data-management.md) — what happens to the file afterwards
