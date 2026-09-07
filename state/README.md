# `state/`

Generated, mutable session state. **Not configuration** — nothing in here is
hand-edited, and nothing in here is committed except this file.

| file | written by | read back? |
|---|---|---|
| `run_state.yml` | [`id8_common.expt_config`](../src/id8_common/expt_config.py) | only its `persistent:` block |
| `.run_state-*` | the same, transiently | no — atomic-write temp files |

## Why it is a directory of its own

`run_state.yml` used to sit loose at the top of the repo. It does not belong
there, and it does not belong in either of the two places you might reach for
first:

- **`src/id8_common/configs/`** holds hand-edited, version-controlled settings
  (`experiment.yml`, `devices.yml`, …). Generated state mixed in there invites
  someone to edit it mid-run, or to commit it.
- **`src/id8_common/`** is the Python package. Writing runtime state into a
  package directory breaks under a non-editable install.

## `run_state.yml`

Three blocks, and only one of them matters:

- **`persistent:`** — read back at every startup. The current sample index, the
  name of the measurement in flight, and the per-sample mesh position table.
  These used to live in `8ideSoft:` `Reg6`, `StrReg8` and `Reg16`–`Reg41`.
  **Do not hand-edit this while a session is running** — the session holds it
  in memory and will overwrite you on its next write.
- **`persistent.measurement_num:`** — the exception. Its store is
  `8ideSoft:Reg1`, not this file; what is written here is only a mirror of the
  last value this checkout saw. See
  [docs/configuration.md](../docs/configuration.md#the-measurement-counter-stays-in-epics)
  for why the counter stayed in EPICS when everything else moved out.
- **`static:`** and **`run:`** — output only, never read back. A debugging
  mirror of what the session currently thinks the experiment settings and the
  in-flight measurement are. To change them, edit
  [`configs/experiment.yml`](../src/id8_common/configs/experiment.yml) or the
  experiment's `measurement_info.yaml`.

Deleting `run_state.yml` resets every sample's mesh position, so the next mesh
scan starts from the beginning. It does **not** reset the measurement counter —
that lives in `8ideSoft:Reg1` and survives anything done to this directory,
which is the whole reason it is kept there. Still, back the file up rather than
delete it.
