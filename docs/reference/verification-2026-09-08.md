[← index](../../README.md)

# Verification run, 2026-09-08

Six detector modes acquired, metadata checked, analysis run. Everything below is
reproducible from the file paths given — the point of this page is that you can
check it yourself rather than take the summary on trust.

## The 30-second version

```
   nexus_xpcs_aps 95ab368  ──installed──►  8id_bits AND 8ide_bits_test
                                           (editable, --no-deps, +1 pip line)
                                                  │
                                                  │  NOT exercised by this run
                                                  ▼
   NX_Class → NX_class                   ID8_NEXUS_WRITER unset, so every
   unit     → units          ──►         measurement used OUR writer:
   + NX_FREQUENCY, NX_PRESSURE           "[writer] id8_common (default)"
                                                  │
                                                  ▼
   6 acquisitions  A0120…A0125  ──►  6 metadata files  ──►  6 results.hdf
   eiger ×4, rigaku ZDT2bit, rigaku EPICS            boost_corr_bin on adamite
```

**⚠ Read this before interpreting anything below.** Two independent things
happened today and it is easy to conflate them:

1. Miaoqi's package was **installed**. It was **not used**. Every one of the six
   measurements went through `utils/nexus_utils.py`, our writer, as the session
   log line `[writer] id8_common (default)` records. His writer runs only under
   `ID8_NEXUS_WRITER=mc`, which nothing sets.
2. The `NX_class`/`units` correction was made to **our** writer. That is what
   these files demonstrate.

So this run validates our corrected writer. It says nothing about his.

## What to check, and where

All paths under `/gdata/dm/8ID/8IDE/2026-3/comm202609/`.

| | protocol | detector / mode | data dir | data file | size |
|---|---|---|---|---|---|
| A0120 | eiger_internal_series | eiger4M / Internal Series | `data/A0120_Test_a0002_f000100_r00001/` | `.h5` | 16 M |
| A0121 | eiger_internal_enable | eiger4M / Internal Enable | `data/A0121_Test_a0002_f000010_r00001/` | `.h5` | 1.8 M |
| A0122 | eiger_external_series | eiger4M / External Series | `data/A0122_Test_a0002_f000100_r00001/` | `.h5` | 31 M |
| A0123 | eiger_external_enable | eiger4M / External Enable | `data/A0123_Test_a0002_f000100_r00001/` | `.h5` | 16 M |
| A0124 | rigaku_zdt2bit | rigaku3M / ZDT2bit | `data/A0124_Test_a0002_f100000_r00001/` | `.bin.000`…`.005` | **118 G** |
| A0125 | rigaku_epics | rigaku3M_epics / EPICS | `data/A0125_Test_a0002_f000100_r00001/` | `.h5` | 4.7 M |

Results are flat in `analysis/Multitau/<name>_results.hdf`.

### Check the attribute spelling yourself

```bash
python3 - <<'PY'
import h5py, collections
f = "/gdata/dm/8ID/8IDE/2026-3/comm202609/data/A0120_Test_a0002_f000100_r00001/A0120_Test_a0002_f000100_r00001_metadata.hdf"
c = collections.Counter()
with h5py.File(f) as h:
    h.visititems(lambda n, o: c.update(o.attrs.keys()))
    c.update(h.attrs.keys())
print(c)          # expect NX_class 143, description 142, units 97 -- and NO NX_Class, NO unit
PY
```

Every one of the six gives the same counts:

```
{'NX_class': 143, 'description': 142, 'units': 97}      144 objects
```

Zero occurrences of the old `NX_Class` or `unit` in any of the six.

### Check the per-detector values

| | `detector_name` | `x_pixel_size` |
|---|---|---|
| A0120–A0123 | `eiger4M` | 7.5e-05 |
| A0124 | `rigaku3M` | 7.6e-05 |
| A0125 | `rigaku3M_epics` | 7.6e-05 |

The pixel size differing per detector is the fix from 2026-09-06; before that
every Rigaku file recorded the Eiger's 75 µm.

### Check the analysis

```bash
ls -la /gdata/dm/8ID/8IDE/2026-3/comm202609/analysis/Multitau/A012*_results.hdf
```

| | `boost_corr_bin` | wall time | result | g2 |
|---|---|---|---|---|
| A0120 | exit 0 | 6 s | 1.79 M | `(21, 36)` all NaN |
| A0121 | exit 0 | 4 s | 1.78 M | `(8, 36)` all NaN |
| A0122 | exit 0 | 6 s | 1.79 M | `(25, 36)` all NaN |
| A0123 | exit 0 | 6 s | 1.80 M | `(21, 36)` all NaN |
| A0124 | exit 0 | **407 s** | 9.18 M | `(61, 36)` **3.90 – 15.19** |
| A0125 | exit 0 | 5 s | 1.19 M | `(21, 36)` **0 – 14425.9** |

Command used, on adamite:

```bash
/home/beams/8IDIUSER/bin/boost_corr_bin \
    -r <data>/<name>.h5 \
    -q /gdata/dm/8ID/8IDE/2026-3/comm202609/data/eiger4m_qmap_default.hdf \
    -o /gdata/dm/8ID/8IDE/2026-3/comm202609/analysis/Multitau \
    -i 0 -t Multitau
```

`-o` is a **directory**. Passing a filename makes boost_corr `mkdir` that name
and write `<name>_results.hdf` inside it.

**⚠ The four Eiger g2 arrays are entirely NaN.** The ring was at 1.09 mA — no
beam — so the Eiger recorded zero counts and the correlation divides by zero.
This run proves the plumbing, not the physics. The two Rigaku datasets have real
values because the Rigaku registers noise at its threshold. **Repeat at least
one Eiger mode with beam before trusting any number from this pipeline.**

### Metadata reaches the results file

The corrected spelling propagates into `*_results.hdf`:

```
A0124_..._results.hdf   NX_class 143, units 97   (all under entry/, ours)
                        unit 4                   (under xpcs/qmap/, boost_corr's own)
```

The four `unit` attributes are boost_corr's own convention on its own datasets;
nothing of ours writes that spelling any more. boost_corr also logs that it read
our file: `metadata filename/type is …A0124…_metadata.hdf | nexus`.

## Cost worth knowing

**A0124 was 118 GB for 100 000 ZDT frames** — one module file alone is 35 GB, and
`rigaku_handler` loads each into CPU RAM because they exceed the A100's 40 GB.
407 s on adamite (1 TB RAM, 112 cores, 4 × A100). On a smaller analysis box this
would fail outright, not merely run slowly.

## Shared state

| | before | after |
|---|---|---|
| `8ideSoft:Reg1` | 120 | 126 | six measurements, by design |
| `8idiSoft:FLIGHT:bo1:8` | OUT | OUT | unchanged |
| `8idEiger4m:cam1:DetectorState_RBV` | **Aborted** | **Idle** | the run cleared a pre-existing fault |
| `configs/experiment.yml` | — | untouched | `analysis_machine` still `polaris` |

DM submission was off for the run (`analysis_machine` forced to `none` in the
test process only, never written to disk), so `boost_corr_bin` was invoked by
hand rather than through a DM workflow.

## One bug fixed along the way

`run_measurement_info(measurement_info_file="/some/path.yaml")` raised
`AttributeError: 'str' object has no attribute 'parent'` before running anything
— it only accepted a `Path`, although the dual module's own usage example shows
a string. Both it and `dry_run_dual_measurement_info()` now coerce with `Path()`.

## Related

* [Two NeXus writers](nexus-writers.md) — why we still use ours, and what switching costs
* [Adding metadata fields](../adding-metadata.md)
* [Data Management](data-management.md)
