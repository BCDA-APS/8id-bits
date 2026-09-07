# Detector modes reference

[← index](../README.md) · [Running measurements](../running-measurements.md)

Thirteen `(detector, mode)` pairs, in nine rows below — each Rigaku ZDT triple
is one row. `detector:` in a protocol must be a key here and `mode:` one of its
modes; anything else is rejected at validation.

## The table

| `detector:` | `mode:` | needs `acq_period` | min `acq_time` | output |
|---|---|---|---|---|
| `eiger4M` | `Internal Series` | no | — | `.h5` |
| `eiger4M` | `Internal Enable` | yes | — | `.h5` |
| `eiger4M` | `External Series` | yes | — | `.h5` |
| `eiger4M` | `External Enable` | yes | 0.1 s | `.h5` |
| `lambda2M` | `Internal` | no | — | `.h5` |
| `lambda2M` | `External` | yes | 0.1 s | `.h5` |
| `rigaku3M` | `ZDT2bit` / `ZDT4bit` / `ZDT8bit` | no | 2/4/8 ×10⁻⁵ s | `.bin.000`…`.005` |
| `rigaku3M_ftf` | `ZDT2bit` / `ZDT4bit` / `ZDT8bit` | no | 2/4/8 ×10⁻⁵ s | `.h5.000`…`.005` \* |
| `rigaku3M_epics` | `EPICS` | no | 0.01 s | `.h5` |

\* The `.000`…`.005` split into six per-module files is confirmed, from runs made
while fast transfer still wrote `.bin`. That the IOC appends the same suffix to a
`.h5` name is **inferred, not yet observed** — `dm_util.py` builds the path on
that assumption and says so. Confirm it on the first fast-transfer run and
correct both places if it differs.

`rigaku3M`, `rigaku3M_ftf` and `rigaku3M_epics` are **the same physical
detector**. The name selects the output format, and it is recorded in the NeXus
`detector_name` field so downstream can tell which was used.

## eiger4M

**Internal Series** — the detector paces its own frames. `acq_period` is not
used: `validate_timing()` overwrites it with `acq_time`, so a value left in the
protocol is ignored rather than rejected (unlike `num_segments` and
`trigger_period`, which are). The simplest of the four: one arm, no external
trigger source, and the only Eiger mode with no timing floor.

**Internal Enable** — one **software** trigger per frame, fired in a loop paced by
`acq_period`.
**⚠ This cannot be pushed.** The loop sleeps a fixed `acq_period` between
triggers and then waits for the detector to leave Acquire, so a period the
detector cannot service leaves it waiting for frames that were never taken.
Use ~1 s / 2 s. It is the only mode whose frame rate is set by a Python loop
rather than by the detector or by softglue, which is what puts the floor so high.

**External Enable** — one softglue pulse per frame, and *that same pulse holds the
shutter open*. Both `acq_time` and `acq_period` must be **≥ 0.1 s** — the shutter
cannot follow faster. Validation enforces it.

**External Series** — the only mode taking `num_segments` and `trigger_period`.
`num_frames` is **per segment**: softglue sends `num_segments` pulses and the file
receives `num_frames × num_segments` frames. `num_segments` defaults to 1 when a
protocol omits it; `trigger_period` is required.

```yaml
    acq_time: 0.1
    acq_period: 0.2
    num_frames: 10        # per segment
    num_segments: 3       # → 30 frames total
    trigger_period: 3     # segment-to-segment spacing
```

The one softglue pulse both starts a segment and holds the shutter open, so each
segment is **two** shutter movements, and `validate_timing()` requires each to be
at least `SHUTTER_MIN_TIME` (0.1 s):

* open — `num_frames × acq_period`, the segment itself;
* closed — `trigger_period − num_frames × acq_period`, the gap to the next pulse.

The closed-gap rule is what keeps `trigger_period` longer than a segment. Too
short and the next pulse lands mid-segment, is dropped, and the acquisition
hangs waiting for a trigger already spent. Within a segment the Eiger paces
itself, so there is **no** 0.1 s floor on `acq_time` or on `acq_period` on its
own here — only on the two shutter times above.

## lambda2M

**Internal** / **External** — same shapes as the Eiger's `Internal Series` and
`External Enable`; `External` is one softglue pulse per frame and carries the
same shutter floor, `acq_time` and `acq_period` both **≥ 0.1 s**.

## rigaku3M — one detector, three names

**`rigaku3M` (ZDT)** — zero-deadtime sparsified output written by the detector's
own fast-file mechanism, six per-module `.bin.000`…`.bin.005` files.

**`rigaku3M_ftf`** — identical acquisition, `output_control = "FastTransfer -
HDF5"`, so the six files are `.h5.000`…`.h5.005`.

**`rigaku3M_epics`** — ordinary areaDetector acquisition through the HDF1 plugin,
one `.h5`. The only Rigaku mode that arms and drains HDF1. It is also the only
one whose output the beamline's boost_corr reads without special handling — see
the analysis table below, and the caveat under it about that being observed
behaviour rather than anything this repo can be checked against.

### ⚠ Analysis support differs by format

| format | boost_corr |
|---|---|
| `rigaku3M_epics` → `.h5` | works |
| `rigaku3M` → `.bin.000` | supported, but see below |
| `rigaku3M_ftf` → `.h5.000` | **not supported** — `create_dataset` has a `.bin.000` branch but no `.h5.000` one |

Also, `Rigaku3MDataset` loads each module file **entirely into RAM** before
correlating, and `-e <n>` does not limit that. A 100 000-frame ZDT dataset can be
hundreds of GB and is not correlatable in practice.

boost_corr is not part of this repo. Both points above are observed behaviour of
the version installed at the beamline, not anything this tree can be checked
against.

### ⚠ Sparsified size depends on the threshold

Sparsification only compresses when the detector sees sparse photon events. One
100 000-frame ZDT run, measured at three thresholds — observed sizes, not
specifications:

| lower threshold | ZDT | FTF |
|---|---|---|
| 3.5 keV | 1.5 TB | 73 GB |
| 4.0 keV | 702 GB | 60 GB |
| 4.5 keV | 118 GB | 28 GB |

A low threshold on a quiet beam records a noise flood and sparsifies nothing.

## Timing rules at a glance

| Rule | Applies to |
|---|---|
| `acq_time > 0` | all |
| `acq_period > 0` and `≥ acq_time` | modes needing `acq_period` |
| `acq_time ≥ min_acq_time` | rigaku modes only — the floor is the mode table's `min_acq_time` |
| `acq_time` and `acq_period` ≥ 0.1 s | eiger External Enable, lambda2M External |
| `num_frames × acq_period ≥ 0.1 s` | eiger External Series — shutter open, one segment |
| `trigger_period − num_frames × acq_period ≥ 0.1 s` | eiger External Series — shutter closed between segments |
| `num_segments`, `trigger_period` accepted | eiger External Series only — on any other mode they are **rejected**, not ignored |

The 0.1 s is `SHUTTER_MIN_TIME` in `master_plan.py`, and it is exported to the
session, so `SHUTTER_MIN_TIME` at the prompt is the number actually enforced.

All of the above are checked by `dry_run_measurement_info()` before anything
moves. The dry run stops there: whether the mode's devices are present and
connected is checked only on the real path, by `run_measurement_info()`.
