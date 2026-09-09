[← index](../../README.md) · [Adding metadata fields](../adding-metadata.md)

# How the NeXus metadata file is written

Since 2026-09-08 the writer is **Miaoqi Chu's `nexus_xpcs_aps`**. This page is
the architecture, and — more usefully — when a problem is ours to fix and when
it needs Miaoqi.

## The 30-second version

```
   configs/devices.yml ─► oregistry ─► live EPICS values
                                              │
   ┌──────────────────────────────────────────▼───────────────────────────────┐
   │  utils/nexus_runtime.py       which signal feeds which path        OURS  │
   │                               ~180 lines of {path: value}                │
   └──────────────────────────────────────────┬───────────────────────────────┘
   ┌──────────────────────────────────────────▼───────────────────────────────┐
   │  utils/xpcs_schema.py      which fields exist            OURS,        │
   │                               244 lines CALLING his factories            │
   └──────────────────────────────────────────┬───────────────────────────────┘
   ┌──────────────────────────────────────────▼───────────────────────────────┐
   │  utils/nexus_writer.py        the entry point + workarounds        OURS  │
   │                               create_nexus_format_metadata(file, det)    │
   └──────────────────────────────────────────┬───────────────────────────────┘
                                              ▼
              nexus_xpcs_aps.core.utils        writes the HDF5           HIS
```

| file | answers | lines | whose |
|---|---|---|---|
| `utils/nexus_runtime.py` | *where the numbers come from* | ~300 | ours |
| `utils/xpcs_schema.py` | *what fields exist* | 244 | ours, calling his |
| `utils/nexus_writer.py` | *the front door* | ~90 | ours |
| `nexus_xpcs_aps.core.*` | *schema factories + HDF5 writer* | — | **his** |

**His package has no runtime layer and cannot have one.** Nothing upstream knows
that `/entry/instrument/detector_1/distance` comes from `device_position.yaml`.
That is why `nexus_runtime.py` is ours — it is 8-ID's wiring, not a gap in his
code.

## ⚠ When to contact Miaoqi

Work out which layer the problem is in *before* writing to him. Most things are
ours.

**Fix it yourself — do not contact him:**

| symptom | where |
|---|---|
| a field holds the wrong *value* | `nexus_runtime.py` — that is our device mapping |
| a field is missing and one of his factories makes it | `xpcs_schema.py` — call `make_slits(9)` etc. |
| a device reads `None` / `AttributeError` at write time | `configs/devices.yml`, or the IOC is down |
| `KeyError` naming a NeXus path | you added a runtime line with no schema node |
| a value is right but in the wrong place in the tree | `xpcs_schema.py` composition |

**Contact him — it is upstream:**

| symptom | why it is his |
|---|---|
| you need a device group **no factory produces** and it is not a slit / attenuator / undulator / detector / diffractometer | his factories are the only way to get a correctly-classed group. `make_sample()` in particular takes fixed flags, so a new *sample environment* cannot be added from outside |
| a leaf *declares* the wrong `units` category, e.g. `NX_ANY` where a real one exists | the category lives in his factory output |
| the same schema writes different files on successive calls | an aliasing or caching bug in `core/` — see the patch list below |
| a NeXus-standard question: class names, attribute spelling, required fields | he owns the schema's conformance |

**How to ask.** Open a pull request against
[AZjk/nexus_xpcs_aps](https://github.com/AZjk/nexus_xpcs_aps) rather than sending
an email — the first one is
[#1](https://github.com/AZjk/nexus_xpcs_aps/pull/1), branched from `mc_refact`.
Include the measured numbers; "two `make_sample()` calls share 17 of 20 leaf
objects" gets a fix faster than "the copy seems shallow".

**Do not do PR work in the installed clone.**
`~/Documents/Miaoqi/nexus_xpcs_aps_95ab368` is pip-installed editable into the
beamline environment, so a branch switch there changes what the beamline imports
*live*. Clone somewhere else.

## Local patches we carry

**None.** [PR #1](https://github.com/AZjk/nexus_xpcs_aps/pull/1) merged on
2026-09-08 and upstream now carries all five fixes we had been holding: the
`make_sample()` leaf aliasing, the `id()`-keyed plan cache, the invalid-JSON
`workflow_kwargs` default, the three missing unit categories
(`NX_VOLTAGE` / `NX_FREQUENCY` / `NX_PRESSURE`), and the
`beam_center_position_x/y` description. The patches were deleted from this repo
in `8b83c67`.

The one thing that looks like a leftover and is not: `nexus_writer.py`
deepcopies the schema on every call. That is **not** working around a bug. His
`update_schema_at_runtime()` assigns `node["data"] = value` in place — by
design — so handing it our module-level `xpcs_schema` would leave one
measurement's values sitting in the template for the next measurement to
inherit. Any caller of his writer needs to do this. Don't remove it.

## Upgrading his package

It is installed editable from a clone, so "upgrading" is a `git pull` in that
clone — which changes what the beamline writes, with no warning and no version
bump. Treat it as a change to this repository:

```bash
# NOT during beamtime
cd ~/Documents/Miaoqi/nexus_xpcs_aps_95ab368
git log --oneline HEAD..origin/mc_refact      # read what would land, first
```

Then re-acquire one measurement per detector mode and confirm each writes a
metadata file and completes analysis, before trusting it.

**⚠ There is a second, stale clone** at `~/Documents/Miaoqi/nexus_xpcs_aps`
(`main` @ 73d0be4, Aug 2025) with uncommitted edits. Nothing should point at it.

## Related

* [Adding metadata fields](../adding-metadata.md) — the step-by-step
