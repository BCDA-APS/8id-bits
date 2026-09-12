# Using `expt` at the prompt

`expt` is the one object every plan reads. This page is about driving it by
hand: what you can read, when, and why a perfectly reasonable-looking read
raises instead of answering.

For *where each setting lives and how to change it*, see
[Configuration](configuration.md). This page is the interactive companion.

## The 30-second version

```
                     WHEN IS IT READABLE?

  ┌─ static ─────────────────── readable the moment the session starts ────┐
  │   expt.cycle_name  expt.experiment_name  expt.mount_point              │
  │   expt.workflow_name*  expt.analysis_machine  expt.use_subfolder       │
  │   from configs/experiment.yml, read once at import                     │
  │   * workflow_name is ALSO a run field -- see "Changing things" below   │
  └────────────────────────────────────────────────────────────────────────┘

  ┌─ persistent ─────────────── readable the moment the session starts ────┐
  │   expt.sample_index   expt.file_name   expt.sample_position(3)         │
  │   from state/run_state.yml, read back at startup                       │
  └────────────────────────────────────────────────────────────────────────┘

  ┌─ EPICS ──────────────────── readable if pv_registers connected ────────┐
  │   expt.measurement_num          ← 8ideSoft:Reg1, shared between        │
  │                                   sessions and machines                │
  │   if pv_registers is down you get a red warning and a per-checkout     │
  │   fallback value, and if there is no fallback either it raises         │
  └────────────────────────────────────────────────────────────────────────┘

  ┌─ run state ──────────── EMPTY until a measurement is loaded ───────────┐
  │   expt.det_name  expt.det_mode  expt.acq_time  expt.num_frames         │
  │   expt.sample_name  expt.header  expt.qmap_file  expt.analysis_type …  │
  │                                                                        │
  │   populated by:  run_measurement_info()      ← a real run              │
  │                  expt.set_measurement(...)   ← by hand                 │
  │                  expt.det_name = "eiger4M"   ← one field at a time     │
  │                                                                        │
  │   NOT populated by dry_run_measurement_info()                          │
  └────────────────────────────────────────────────────────────────────────┘
```

**Check which state you are in by typing `expt`:**

```python
In [1]: expt
Out[1]: <ExperimentConfig 'comm202609' cycle='2026-3' run_keys=0>
                                                      ^^^^^^^^^^^
                                       0 = no measurement loaded yet
```

`run_keys=0` is the answer to "why did `expt.det_name` raise?".

**⚠ Nonzero does not mean your field is set.** `run_keys` is just how many of the
21 run fields have values. A parallel (multi-detector) acquisition loads only the
sample half (`run_keys=10`) and leaves `det_name` unset on purpose, because it is
per leg.
Zero is proof that nothing is loaded; anything else proves nothing in
particular.

## The two errors everyone hits first

### 1. `'det_name' has not been set yet`

```python
In [9]: expt.det_name
AttributeError: 'det_name' has not been set yet -- run state is populated by
master_plan.run_measurement() ...
```

Nothing is broken. Run-state fields describe *the measurement currently being
run*, and no measurement has been run in this session. The field has no value
yet, so `expt` says so rather than handing back a stale one from last week.

**⚠ `dry_run_measurement_info()` does not fix this.** It validates and prints
the whole queue without loading any of it into `expt` — that is the point of a
dry run: it changes nothing. The two things that do call `set_measurement()` are
`run_single_measurement()` (everything) and `run_multi_measurement()` (the sample
half only) — both reached through `run_measurement()`.

To populate it without acquiring anything:

```python
# Option A -- load one measurement's worth of fields by hand.
expt.set_measurement(
    measurement={"detector": "eiger4M", "mode": "Internal Series",
                 "acq_time": 0.01, "acq_period": 0.01,
                 "num_frames": 100, "num_repeats": 1,
                 "sample_move": "no",          # omit this and det_acq_series()
                                               # fails partway, after burning a
                                               # measurement number
                 "qmap_file": "eiger4m_qmap_default.hdf"},
    sample={"sample_name": "Test", "header": "A"},
)

# Option B -- set the one field you care about.
expt.det_name = "eiger4M"
```

Both write `state/run_state.yml` as a side effect, so a failed run leaves a
record on disk of what it was trying to do, and a GUI can see the current
measurement without a live session. **They do not survive a restart** — only the
`persistent:` block is read back, so a new session starts at `run_keys=0` again.

### 2. `'PosixPath' object is not callable`

```python
In [6]: expt.measurement_info_file()
TypeError: 'PosixPath' object is not callable
```

Some of these take parentheses and some do not, and there is no rule that
predicts it — the table below is the list.

What *is* reliable: **none of the four field buckets ever takes parentheses.**
`expt.cycle_name`, `expt.acq_time`, `expt.sample_index` and
`expt.measurement_num` are all plain attribute reads (resolved by `__getattr__`,
not Python `@property` objects, though they behave the same way at the prompt).
Everything else — the path helpers and the verbs — is in the table.

| type it as | what you get |
|---|---|
| `expt.cycle_name` | property → `'2026-3'` |
| `expt.measurement_info_file` | property → `PosixPath('…/measurement_info.yaml')` |
| `expt.sample_info_file` | property → `PosixPath('…/sample_info.yaml')` |
| `expt.trio_measurement_info_file` | property → `PosixPath('…/trio_measurement_info.yaml')` |
| `expt.user_plan_dir` | property → the folder those three live in |

**⚠ The four path entries raise `FileNotFoundError`, not `None`,** when
`src/user_plans/<cycle_name>/<experiment_name>/` does not exist. The message
names both values so you can see which one in `experiment.yml` is wrong. They
are built fresh on every read, so fixing `experiment.yml` and calling
`expt.reload()` is enough — no restart needed.
| `expt.sample_position(3)` | **method** → mesh index for sample 3 |
| `expt.set_sample_position(3, -1)` | **method** → rewind sample 3's mesh |
| `expt.reload()` | **method** → re-read `experiment.yml` after editing it |
| `expt.as_dict()` | **method** → everything static + run, as one dict |

`expt.sample_position()` with no argument raises `TypeError: missing 1 required
positional argument` — it needs to know which sample. Sample indices start at
**1**, and `-1` means "this sample's mesh has not been started yet".

## Reading everything at once

```python
expt                      # one-line summary; run_keys tells you if run state is loaded
expt.as_dict()            # static + run state merged, as a plain dict
expt.sample_position(1)   # one sample's mesh index
[expt.sample_position(i) for i in range(1, 6)]     # the first five
```

To see the persistent block as it is on disk:

```bash
cat ~/bluesky/state/run_state.yml
```

## Changing things

| to change | do this | not this |
|---|---|---|
| cycle, experiment, mount point | edit `configs/experiment.yml`, then `expt.reload()` | `expt.cycle_name = …` (raises) |
| what a measurement does | edit `measurement_info.yaml` | — |
| one run field, for a manual `det_acq_series()` | `expt.acq_time = 0.05` | — |
| a sample's mesh position | `expt.set_sample_position(3, -1)` | hand-edit `run_state.yml` |
| the measurement counter | leave it alone; it is `8ideSoft:Reg1` | hand-edit `run_state.yml` |

**⚠ `workflow_name` is the one field in two buckets.** It is in
`STATIC_FIELDS` *and* `RUN_FIELDS`, and writes try run state first — so
`expt.workflow_name = "x"` does **not** raise like the other static fields. It
silently sets a per-measurement override that `expt.reload()` will not clear
(reload only rebuilds the static table), and `dm_util` and `nexus_utils` then
read it for every later measurement. To change it for the experiment, edit
`experiment.yml`; to clear an override, restart the session.

Assigning to any other static field raises on purpose, and names the file to
edit:

```python
In [2]: expt.cycle_name = "2026-4"
AttributeError: 'cycle_name' is a per-experiment setting: edit
…/configs/experiment.yml and call expt.reload(), rather than setting it here.
```

**⚠ A typo'd *read* raises; a typo'd *write* does not — and reading it back
will not tell you.** `expt.acq_tmie` on its own raises, because `__getattr__`
runs only when normal lookup fails. But `expt.acq_tmie = 0.01` falls through to
`object.__setattr__` and puts it in the instance dictionary, after which
`expt.acq_tmie` finds it by normal lookup and cheerfully returns `0.01`. The
value looks set, and no plan will ever read it.

Check the name itself, not the value:

```python
"acq_tmie" in vars(expt)     # True  -> you created a junk attribute
"acq_time" in expt.as_dict() # True  -> this is a real field
```

## Related

* [Configuration](configuration.md) — where each setting lives, and the four buckets in full
* [Running measurements](running-measurements.md) — protocol syntax, and what a run does
* [Adding metadata fields](adding-metadata.md) — getting a value out of `expt` and into the NeXus file
