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

## Adding something: three cases, and they differ a lot

### Case 1 — a whole new device that reuses one of his factories

**His wins outright, and needs nothing from him.** The factories are plain
functions that take an index and return a fresh dict — no registration, no
central list. The index is only interpolated into the description, so it can even
be a string:

```python
make_slits(9)              ->  4 fields, NXslit
make_slits("my_new_slit")  ->  4 fields, NXslit
make_attenuator(9)         ->  2 fields, NXattenuator
make_undulator(9)          ->  3 fields, NXinsertion_device
make_detector(9)           -> 19 fields, NXdetector
make_diffractometer("kappa") -> 6 fields, NXpositioner
```

You place the result under whatever key you like, in **our**
`xpcs_schema_mc.py` — his code is not touched. The current schema already does
this four times for slits, twice with non-numeric indices:

```python
"sl4":       _tag(make_slits(4, description="Slits 4"),          "/entry/instrument/sl4",       "his:make_slits"),
"wb_slit":   _tag(make_slits("wb", description="white beam slit"), "/entry/instrument/wb_slit",   "his:make_slits"),
"mono_slit": _tag(make_slits("mono", description="mono beam slit"), "/entry/instrument/mono_slit", "his:make_slits"),
```

One line for a complete, correctly-classed NeXus group. Under our writer the
same group is ~30 lines of literal dict, hand-written and hand-maintained.

**⚠ There is no `make_stage` and no `make_motor`.** The complete list is
`make_undulator`, `make_detector`, `make_slits`, `make_attenuator`,
`make_diffractometer`, `make_sample`, `make_entry`. A motor stage is Case 3.

### Case 2 — one more field inside a group he already models

Free with his: the leaf arrives in the factory output, and an upstream fix or
addition reaches you on a `git pull`. Two steps with ours (schema node + runtime
line) versus one with his (runtime line only).

### Case 3 — a field or group he does not model at all

Neither writer helps. His factories return fixed leaf sets — `make_slits` gives
exactly four fields, and there is no arbitrary-leaf helper — so you either ask
him to add it upstream and wait for a release, or patch the composed dict
locally the way `_rename()` already patches `flightpath_swing_horizontal`. That
local patch is about the same work as adding the literal node to ours.

### Either way

The runtime half never changes hands: the line that reads the device and puts
its value at a path lives in **our** `create_runtime_metadata_dict()` in both
writers. His package does not model "which EPICS signal feeds this field".

## What each costs

| | `nexus_utils.py` (ours) | `nexus_xpcs_aps` (his) |
|---|---|---|
| schema | 993 lines of literals, all ours to maintain | 244 lines of composition; the leaves are his |
| who fixes a NeXus-standard bug | us | him, and it arrives on a pull |
| adding a device that reuses an existing factory | ~30 lines of literal dict | one line, no upstream change |
| adding something he does not model | a literal node | a local patch, or wait for his release |
| dependency risk | none | none in practice: `h5py` + `numpy` |
| coupling to `apsbits` | none | none, **provided only `core/` is used** — his `deployment/id8_{e,i}` needs apsbits and is obsolete |
| behaviour if it breaks mid-beamtime | ours to fix | ours to work around |

## Why we have not switched

Honestly: **it is a test that was left switched off, not a considered rejection.**
The comparison was run, came out clean, and then the DM/analysis investigation
took over. Two reasons were live at the time, and only one still is:

1. *"Installing risks disturbing the `8id_bits` environment before beamtime."*
   Largely gone. The declared dependencies are `h5py` and `numpy`, both already
   present. A `pip install -e` of the clone would add no third-party package.
2. *"`deployment/` is unusable."* Still true, and still fine — it needs
   `apsbits`, which this instrument is moving away from. `core/` is the only
   part we import, and it is self-contained.

The remaining reason to wait is timing, not technology: swapping the writer
changes the bytes of every metadata file, and `*_results.hdf` is a copy of that
file. That is not a change to make two days before a run.

**On the merits, though, Case 1 is the strongest argument for switching** — most
new hardware at this beamline is another slit, attenuator, undulator or
detector, and each of those is one line with his factories against ~30 lines of
literal dict with ours.

## What switching would take

1. `pip install -e ~/Documents/Miaoqi/nexus_xpcs_aps_95ab368` into `8id_bits`,
   so `PYTHONPATH` juggling stops.
2. Delete the stale `~/Documents/Miaoqi/nexus_xpcs_aps` clone so nobody points
   at Aug 2025 by accident.
3. Make `mc` the default in `det_acq_series()`, and extend it to the two call
   sites that ignore the variable today — `cleanup_acquisition()` (the abort
   path) and the dual-detector path — or a run can emit files from both writers.
4. Agree with Miaoqi which of our 125 leaves he owns, so a field we need is not
   silently dropped by an upstream refactor.
5. Fold `xpcs_schema.py` and `xpcs_schema_mc.py` into one, and delete
   `nexus_utils_mc.py`.
6. Re-run the comparison on one measurement per detector mode.

Steps 3 and 4 are the substantive ones. Step 4 especially: pulling his changes
means his refactor can change our files, which is the cost of him maintaining
the schema.

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
