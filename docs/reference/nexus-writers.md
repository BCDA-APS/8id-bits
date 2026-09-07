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

## Adding a field: the same job, both ways

### A field his factories already model

| | ours | his |
|---|---|---|
| get the field | write the literal node by hand | free — it arrives in the factory output |
| keep it current | nobody | his upstream change flows in on a `git pull` |
| steps | 2 (schema node + runtime line) | 1 (runtime line) |

His wins here, clearly.

### A new EPICS device you just added — the case that actually comes up

**Ours** — a literal node in `xpcs_schema.py`:

```python
"temperature_set": {
    "type": "NX_FLOAT",
    "required": True,
    "units": "NX_TEMPERATURE",
    "description": "Sample temperature setpoint",
    "data": 0,
},
```

**His** — there is no "add an arbitrary leaf" helper. `core/` exposes
per-group factories (`make_detector`, `make_undulator`, `make_slits`,
`make_attenuator`, `make_diffractometer`, `make_sample`, `make_entry`) that
return fixed dicts. So a field he does not model means either:

* ask him to add it upstream and wait for a release — he controls the schema,
  which is the flip side of him maintaining it; or
* patch the composed dict locally in `xpcs_schema_mc.py`, the same way the
  existing `_rename()` call patches `flightpath_swing_horizontal`.

Then the runtime line, identically, in **our** `create_runtime_metadata_dict()`.

**So for a beamline-specific field, his is not fewer steps — it is the same
number plus a layer of indirection**, or the same number plus a round trip
through someone else's release cycle. The saving is real only for fields that
are already in his model.

## What each costs

| | `nexus_utils.py` (ours) | `nexus_xpcs_aps` (his) |
|---|---|---|
| schema | 993 lines of literals, all ours to maintain | 244 lines of composition; the leaves are his |
| who fixes a NeXus-standard bug | us | him, and it arrives on a pull |
| who can add a beamline-specific field today | us, alone | us locally, or him upstream |
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
