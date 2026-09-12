"""Experiment configuration and per-measurement run state.

Replaces most of what ``pv_registers`` was carrying. Four backends sit
behind one object, ``expt``, so a call site says what it wants and not where
the value happens to live:

``configs/experiment.yml``  -- static, per-experiment. Changes a few times a
    cycle and is edited by hand: cycle name, mount point, experiment name,
    analysis machine, workflow, subfolder policy. Detector pixel size is not
    one of them -- it is per detector, in plans/set/device_position.yaml.

run state (in memory)       -- the parameters of the measurement currently
    being run. Populated by ``master_plan.run_measurement()`` straight from
    ``measurement_info.yaml`` and read back by ``ad_acq``/``acq_helpers``/
    ``nexus_utils``.

    These used to be written into ~20 EPICS registers by
    ``write_measurement_registers()`` and read straight back out a few
    function calls later, in the same process -- a round trip through Channel
    Access purely to pass arguments. That is what this replaces.

    Mirrored to the ``run:`` block of ``state/run_state.yml`` after every
    change, so a failed run leaves behind what it was trying to do, and so a
    GUI or a shell script can still see the current measurement without a live
    session. Nothing reads that block back -- it is an output, not a store.

persistent state             -- the handful of values that must survive between
    measurements and between sessions: ``sample_index``, ``file_name``, and the
    per-sample mesh positions. Stored in the ``persistent:`` block of
    ``state/run_state.yml``, which unlike the other two blocks IS read back at
    startup. These were ``Reg6``/``StrReg8``/``Reg16-41``.

``pv_registers``            -- the values that genuinely want to be an EPICS
    register, because losing them costs data rather than convenience. Just
    ``measurement_num`` (``8ideSoft:Reg1``) -- see PV_FIELDS for why that one
    stayed behind when the rest moved out.

Usage::

    from id8_common.expt_config import expt

    expt.cycle_name             # static     -> str
    expt.acq_time               # run        -> float
    expt.sample_index           # persistent -> int
    expt.measurement_num        # EPICS Reg1 -> int
    expt.measurement_num += 1   # EPICS write, mirrored to run_state.yml

An unknown attribute raises ``AttributeError`` naming the four buckets,
rather than silently returning ``None``.
"""

import os
import tempfile
from pathlib import Path

import yaml


CONFIGS_DIR = Path(__file__).parent / "configs"
EXPERIMENT_FILE = CONFIGS_DIR / "experiment.yml"

# <repo>/state/run_state.yml. Deliberately NOT under /gdata (this is session
# run state, not data -- what matters for the record is copied into each
# measurement's NeXus metadata), and deliberately not under configs/ or inside
# the id8_common package either: configs/ holds hand-edited, version-controlled
# settings, and a package directory is code, which a non-editable install would
# make read-only. A separate gitignored state/ directory keeps generated,
# mutable, never-hand-edited files clearly apart from both. Moved here from the
# repo root on 2026-09-06.
STATE_DIR = Path(__file__).resolve().parents[2] / "state"
RUN_STATE_FILE = STATE_DIR / "run_state.yml"

# Per-experiment plan files live under <repo>/src/user_plans/<cycle>/<experiment>/.
# The cycle and experiment come from experiment.yml, so pointing the session at a
# different experiment is a two-line edit there -- nothing in the plans is pinned
# to a path. Before 2026-09-06 these were three hardcoded absolute paths in
# master_plan.py and select_sample.py, flat under user_plans/.
USER_PLANS_ROOT = Path(__file__).resolve().parents[1] / "user_plans"

# name -> coercion applied on load. Values come from experiment.yml.
STATIC_FIELDS = {
    "cycle_name": str,
    "mount_point": str,
    "experiment_name": str,
    "analysis_machine": str,
    "workflow_name": str,
    "use_subfolder": "yes_no",
}

# name -> coercion applied on write. Populated per measurement.
RUN_FIELDS = {
    # from the protocol block of measurement_info.yaml
    "det_name": str,
    "det_mode": str,
    "acq_time": float,
    "acq_period": float,
    "num_frames": int,
    "num_repeats": int,
    "num_segments": int,
    "trigger_period": float,
    "sample_move": "yes_no",
    "qmap_file": str,
    "analysis_type": str,
    # Also a STATIC_FIELDS key: experiment.yml gives the default, and a
    # protocol (or a trio-acquisition leg) may override it for one measurement.
    # Reads prefer the run state and fall back to the static value.
    "workflow_name": str,
    # from the sample block of sample_info.yaml
    "header": str,
    "sample_name": str,
    "inner_motor": str,
    "outer_motor": str,
    "inner_center": float,
    "outer_center": float,
    "inner_range": float,
    "outer_range": float,
    "inner_pts": int,
    "outer_pts": int,
}

# Persistent session state: the values that must survive between measurements
# and between sessions, and whose loss costs continuity rather than data. These
# used to be EPICS registers (8ideSoft: Reg6/StrReg8/Reg16-41); they are now
# stored in RUN_STATE_FILE, whose `persistent:` block IS read back at import --
# unlike the static/run mirrors in the same file, which are output only.
#
# Reg1 (measurement_num) deliberately did NOT come with them -- see PV_FIELDS.
PERSISTENT_FIELDS = {
    "sample_index": int,      # was Reg6
    "file_name": str,         # was StrReg8 -- name of the measurement in flight
}

# Fields whose store is an EPICS register on `pv_registers`, not run_state.yml.
# The attribute name here is also the Component name on
# devices/registers_device.py:EpicsPvStorageRegisters.
#
# measurement_num is the one value where losing the store means overwriting
# data. It is the NNNN in every folder and file name (A0061_Test_a0002_...),
# it only ever counts up, and nothing downstream checks whether a name is
# already taken -- so a counter that restarts at 0 silently writes over
# existing measurements. run_state.yml cannot carry that risk: it is
# gitignored, so a `git clean -x`, a fresh clone, or an accidental delete
# resets it, and it is per-checkout, so a second session on another machine
# would hand out numbers this one has already used.
#
# 8ideSoft:Reg1 has neither problem. It outlives the checkout, it is visible to
# `caget` and to non-Bluesky tools, and every process that opens it sees the
# same number. run_state.yml still keeps a MIRROR of it (see _pv_mirror) --
# read only to catch the register itself going backwards, never as the store.
#
# One counter, TWO naming streams -- worth knowing before checking what number
# is already taken. det_acq_series() and trio_acq_series() name their output
# <...>/data/<header>####... through acq_helpers.gen_folder_prefix(), while the
# align scans in plans/align/ name theirs <...>/data/bluesky/...  Both advance
# this counter once per scan, so the highest number on disk may be under
# data/bluesky/, not data/.
#
# ophyd_scan.py does NOT go through gen_folder_prefix(): it builds its own name
# from the scan's arguments and advances this counter itself -- see
# ophyd_scan.scan_file_name. Its files start "S" whatever sample is loaded;
# everything else keeps the sample's header.
PV_FIELDS = {
    "measurement_num": int,   # 8ideSoft:Reg1
}

# Fallback values for PERSISTENT_FIELDS, used when RUN_STATE_FILE has no entry
# for a field -- or no file at all, on a fresh checkout.
#
# The per-sample mesh positions (was Reg16..Reg41 = sample1_pos..sample26_pos)
# are deliberately NOT listed here: they live in a separate dict keyed by sample
# index (self._sample_positions) so they do not need 26 named fields, and an
# index that has never been visited defaults to -1 -- see sample_position().
PERSISTENT_DEFAULTS = {
    "sample_index": 1,
    "file_name": "",
}


#: snapshot_run() marker for a run-state key that had no value to save, so
#: restore_run() removes it again instead of writing one back. Used but never
#: defined until 2026-09-07 -- snapshot_run() raised NameError on any unset
#: field, which is precisely the trio-acquisition case it was written for.
_ABSENT = object()

#: `_pv_dev` before the first lookup -- distinct from None, which records that
#: the lookup already ran and pv_registers is not available this session.
_UNRESOLVED = object()


def _coerce(value, how):
    """Apply one of the coercions named in the field tables."""
    if value is None:
        return None

    if how == "yes_no":
        # YAML turns a bare `no` into False, so a hand-edited
        # `use_subfolder: no` or `sample_move: no` must still compare equal
        # to "no" downstream.
        if value is True:
            return "yes"
        if value is False:
            return "no"
        text = str(value).strip().lower()
        if text in ("yes", "no"):
            return text
        raise ValueError(f"expected yes or no, got {value!r}")

    if how is str:
        return str(value).strip()

    return how(value)


class ExperimentConfig:
    """One attribute namespace over experiment.yml, run state, and persistent state."""

    def __init__(self):
        # object.__setattr__ throughout __init__: our own __setattr__ routes
        # by field name and would not know what to do with these.
        object.__setattr__(self, "_static", {})
        object.__setattr__(self, "_run", {})
        object.__setattr__(self, "_persistent", dict(PERSISTENT_DEFAULTS))
        object.__setattr__(self, "_sample_positions", {})
        # PV-backed fields: the device is resolved lazily on first use, because
        # expt_config is imported before make_devices() has run.
        object.__setattr__(self, "_pv_dev", _UNRESOLVED)
        object.__setattr__(self, "_pv_mirror", {})
        object.__setattr__(self, "_pv_checked", set())
        self.reload()
        self.load_persistent()

    # ---------------------------------------------------------------- static

    def reload(self):
        """Re-read experiment.yml. Safe to call at the prompt after editing it."""
        if not EXPERIMENT_FILE.exists():
            raise FileNotFoundError(
                f"{EXPERIMENT_FILE} not found -- it holds the per-experiment "
                f"settings (cycle name, mount point, experiment name, ...) that "
                f"used to live in pv_registers StrReg1-7."
            )

        raw = yaml.safe_load(EXPERIMENT_FILE.read_text()) or {}
        static = {}

        for name, how in STATIC_FIELDS.items():
            if name not in raw:
                raise KeyError(f"{EXPERIMENT_FILE} is missing required key {name!r}")
            static[name] = _coerce(raw[name], how)

        object.__setattr__(self, "_static", static)
        return static

    # ------------------------------------------------------------ run state

    def set_measurement(self, measurement=None, sample=None):
        """Load one measurement's parameters into the run state.

        Called by master_plan.run_measurement() with the already-validated
        protocol and sample blocks. Keys absent from either block are left at
        whatever the previous measurement set, except the three that a mode may
        legitimately not declare -- num_segments, trigger_period and
        analysis_type -- which are reset to their defaults every time so an
        External Series run cannot leak a stale value into the mode that
        follows it.
        """
        if measurement is not None:
            # workflow_name defaults to the experiment.yml value: drop any
            # per-leg override left by a previous measurement rather than
            # letting it leak into this one.
            if "workflow_name" in measurement:
                self._run["workflow_name"] = _coerce(measurement["workflow_name"], str)
            else:
                self._run.pop("workflow_name", None)

            self._run["num_segments"] = int(measurement.get("num_segments", 1))
            self._run["trigger_period"] = float(measurement.get("trigger_period", 0))
            self._run["analysis_type"] = str(measurement.get("analysis_type", "Multitau"))

            # (key in measurement_info.yaml, field name in the run state).
            # Most pairs are the same word twice; only the first two are
            # renamed on the way in.
            for key, target in (
                ("detector", "det_name"),
                ("mode", "det_mode"),
                ("acq_time", "acq_time"),
                ("acq_period", "acq_period"),
                ("num_frames", "num_frames"),
                ("num_repeats", "num_repeats"),
                ("sample_move", "sample_move"),
                ("qmap_file", "qmap_file"),
            ):
                if key in measurement:
                    self._run[target] = _coerce(measurement[key], RUN_FIELDS[target])

        if sample is not None:
            for key in (
                "header",
                "sample_name",
                "inner_motor",
                "outer_motor",
                "inner_center",
                "outer_center",
                "inner_range",
                "outer_range",
                "inner_pts",
                "outer_pts",
            ):
                if key in sample:
                    self._run[key] = _coerce(sample[key], RUN_FIELDS[key])

        self.dump()
        return dict(self._run)

    # ------------------------------------------------------- per-experiment paths

    @property
    def user_plan_dir(self):
        """``<repo>/src/user_plans/<cycle_name>/<experiment_name>/``.

        Derived, not stored, so it follows experiment.yml automatically. Raises
        with both the missing path and the two settings it was built from --
        otherwise a wrong cycle or experiment name surfaces as a bare
        FileNotFoundError on a file the user never named.
        """
        folder = USER_PLANS_ROOT / self.cycle_name / self.experiment_name

        if not folder.is_dir():
            raise FileNotFoundError(
                f"{folder} does not exist. It is built from configs/experiment.yml: "
                f"cycle_name={self.cycle_name!r}, experiment_name={self.experiment_name!r}."
            )

        return folder

    @property
    def sample_info_file(self):
        return self.user_plan_dir / "sample_info.yaml"

    @property
    def measurement_info_file(self):
        return self.user_plan_dir / "measurement_info.yaml"

    @property
    def trio_measurement_info_file(self):
        return self.user_plan_dir / "trio_measurement_info.yaml"

    def snapshot_run(self, names):
        """Capture run-state values for `names`, for later restore_run().

        A name that has no run-state value yet is recorded as absent rather than
        raising, so a caller can temporarily set fields that were never
        populated -- which is exactly what a trio acquisition does: it sets
        det_name/qmap_file/analysis_type one leg at a time and they have no
        meaningful global value in between.
        """
        snap = {}
        for name in names:
            if name not in RUN_FIELDS:
                raise KeyError(f"{name!r} is not a run-state field")
            snap[name] = self._run[name] if name in self._run else _ABSENT
        return snap

    def restore_run(self, snapshot):
        """Undo a snapshot_run(), including removing keys that were absent."""
        for name, value in snapshot.items():
            if value is _ABSENT:
                self._run.pop(name, None)
            else:
                self._run[name] = value
        self.dump()

    def as_dict(self):
        """Everything the run state and static config currently hold."""
        merged = dict(self._static)
        merged.update(self._run)
        return merged

    def dump(self):
        """Mirror the run state to run_state.yml, atomically.

        Written to a temp file in the same directory and renamed over the
        target, so an interrupted write leaves the previous file intact
        rather than a half-written one.
        """
        # Build the persistent block in three plain steps rather than one
        # nested dict(mapping, **other, keyword=...) call: the persistent
        # fields, then a MIRROR of the PV-backed fields (not the store -- see
        # PV_FIELDS), then the mesh positions, sorted by sample index so the
        # file does not reshuffle itself between writes.
        persistent_block = dict(self._persistent)
        persistent_block.update(self._pv_mirror)
        persistent_block["sample_positions"] = dict(sorted(self._sample_positions.items()))

        payload = {
            "_note": "Written by id8_common.expt_config. The 'persistent' block IS read "
                     "back at startup -- it holds the sample index and mesh positions "
                     "that used to be EPICS registers, so do not hand-edit it while a "
                     "session is running. measurement_num appears there too but is only "
                     "a MIRROR: the store is 8ideSoft:Reg1, and editing it here does "
                     "nothing except (if you raise it) make the next session push the "
                     "register up to match. 'static' and 'run' are output only; edit "
                     "configs/experiment.yml or measurement_info.yaml instead.",
            "persistent": persistent_block,
            "static": dict(self._static),
            "run": dict(self._run),
        }

        try:
            # A fresh clone has no state/ yet, and mkstemp will not create it.
            RUN_STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp = tempfile.mkstemp(dir=str(RUN_STATE_FILE.parent), prefix=".run_state-")
            with os.fdopen(fd, "w") as handle:
                yaml.safe_dump(payload, handle, default_flow_style=False, sort_keys=False)
            os.replace(tmp, RUN_STATE_FILE)
        except OSError as exc:
            # Losing the mirror must never take down an acquisition.
            print(f"[expt_config] could not write {RUN_STATE_FILE}: {exc}")

    # ------------------------------------------------------------------ PVs

    def _pv_device(self):
        """The pv_registers device, or None if it is not usable this session.

        Resolved on first use rather than in __init__: this module is imported
        before make_devices() has run, and it must stay importable with no
        devices at all (offline dry runs, unit-testing a plan).
        """
        device = object.__getattribute__(self, "_pv_dev")

        if device is _UNRESOLVED:
            try:
                from id8_common.registry import get_connected_device

                device = get_connected_device("pv_registers")
            except Exception as exc:
                device = None
                print(
                    f"\033[91m[expt_config] pv_registers is unavailable ({exc}). "
                    f"{sorted(PV_FIELDS)} fall back to the {RUN_STATE_FILE.name} mirror "
                    f"for this session -- the counter is NOT shared with other "
                    f"processes, so check for existing files before acquiring.\033[0m"
                )
            object.__setattr__(self, "_pv_dev", device)

        return device

    def _read_pv(self, name):
        """Read a PV-backed field, keeping the run_state.yml mirror in step."""
        device = self._pv_device()
        mirror = self._pv_mirror.get(name)

        if device is None:
            if mirror is None:
                raise RuntimeError(
                    f"{name!r} lives in pv_registers, which is not available, and "
                    f"{RUN_STATE_FILE} has no mirrored value to fall back on. Set it "
                    f"explicitly (expt.{name} = ...) before acquiring, or start a "
                    f"session with pv_registers connected."
                )
            return mirror

        signal = getattr(device, name)

        # Uncached: ophyd serves EpicsSignal.get() from its monitor cache, which
        # lags a write made by another process by however long the CA callback
        # takes. For a counter whose whole purpose is that two sessions never
        # hand out the same number, a stale read is the one failure that
        # matters, so pay for a real Channel Access get -- this runs once per
        # acquisition, not in any loop. Signals without the kwarg (sim, plain
        # Signal) fall back to the ordinary read.
        try:
            raw = signal.get(use_monitor=False)
        except TypeError:
            raw = signal.get()

        value = _coerce(raw, PV_FIELDS[name])

        # Once per session, check the register has not gone backwards relative
        # to the last value this checkout saw. A soft-IOC restart that loses
        # Reg1, or a register cleared by hand, would otherwise hand out numbers
        # that already name files on disk. Only ever corrects upwards.
        if name not in self._pv_checked:
            self._pv_checked.add(name)

            if mirror is not None and mirror > value:
                print(
                    f"\033[91m[expt_config] pv_registers.{name} reads {value}, but "
                    f"{RUN_STATE_FILE.name} last saw {mirror}. The register has gone "
                    f"backwards (IOC restart?). Advancing it to {mirror} so existing "
                    f"files are not overwritten.\033[0m"
                )
                getattr(device, name).put(mirror)
                value = mirror

        if mirror != value:
            self._pv_mirror[name] = value
            self.dump()

        return value

    def _write_pv(self, name, value):
        """Write a PV-backed field, and mirror it into run_state.yml."""
        value = _coerce(value, PV_FIELDS[name])
        device = self._pv_device()

        if device is not None:
            getattr(device, name).put(value)

        self._pv_mirror[name] = value
        self.dump()

    # ------------------------------------------------- persistent session state

    def load_persistent(self):
        """Read the persistent block back out of RUN_STATE_FILE.

        Missing file or missing keys fall back to PERSISTENT_DEFAULTS rather
        than raising -- a fresh checkout has no run_state.yml yet.
        """
        data = {}
        try:
            if RUN_STATE_FILE.exists():
                data = (yaml.safe_load(RUN_STATE_FILE.read_text()) or {}).get("persistent") or {}
        except Exception as exc:
            print(f"[expt_config] could not read {RUN_STATE_FILE}: {exc}")

        persistent = dict(PERSISTENT_DEFAULTS)
        for name, how in PERSISTENT_FIELDS.items():
            if name in data:
                try:
                    persistent[name] = _coerce(data[name], how)
                except Exception:
                    pass

        positions = {}
        for key, value in (data.get("sample_positions") or {}).items():
            try:
                positions[int(key)] = int(value)
            except Exception:
                pass

        # Seed the PV mirror from the same block. This is NOT the store -- it
        # is only what this checkout last saw, used by _read_pv() to notice the
        # register going backwards. A missing key leaves the field unmirrored,
        # and the register's own value is then taken as correct.
        mirror = {}
        for name, how in PV_FIELDS.items():
            if name in data:
                try:
                    mirror[name] = _coerce(data[name], how)
                except Exception:
                    pass

        object.__setattr__(self, "_persistent", persistent)
        object.__setattr__(self, "_sample_positions", positions)
        object.__setattr__(self, "_pv_mirror", mirror)
        return persistent

    def sample_position(self, sample_index):
        """Last mesh index used for this sample. -1 means start from the beginning.

        Validates the index. The EPICS version of this store was 26 named registers
        (sample1_pos..sample26_pos), so an out-of-range index raised AttributeError
        and callers relied on that check; a dict has no such bound, so it is explicit.
        """
        index = int(sample_index)

        if index < 1:
            raise ValueError(f"sample_index must be >= 1, got {index}")

        return int(self._sample_positions.get(index, -1))

    def set_sample_position(self, sample_index, value):
        self._sample_positions[int(sample_index)] = int(value)
        self.dump()

    # ------------------------------------------------------- attribute routing

    def __getattr__(self, name):
        # Only reached when normal lookup fails, so methods and properties
        # defined above never come through here.
        #
        # Which table a name comes from, in the order they are consulted:
        #   1. RUN_FIELDS -- but only if the field has actually been set
        #   2. STATIC_FIELDS   -> experiment.yml
        #   3. RUN_FIELDS again -- a run field with nothing in it yet; this
        #      branch exists to raise a message naming who is meant to set it
        #   4. PERSISTENT_FIELDS -> the persistent block of run_state.yml
        #   5. PV_FIELDS       -> an EPICS register on pv_registers
        # Steps 1 and 3 sit either side of step 2 so that a name in BOTH tables
        # (workflow_name is the only one) uses the per-measurement override when
        # one has been set, and otherwise falls back to the experiment.yml value.
        #
        # object.__getattribute__ rather than self._run: plain attribute access
        # from inside __getattr__ would come straight back here and recurse
        # forever if the dict were ever missing (e.g. part-built instance).
        run = object.__getattribute__(self, "_run")

        if name in RUN_FIELDS and name in run:
            return run[name]

        if name in STATIC_FIELDS:
            return object.__getattribute__(self, "_static")[name]

        if name in RUN_FIELDS:
            if name not in run:
                raise AttributeError(
                    f"{name!r} has not been set yet -- run state is populated by "
                    f"master_plan.run_measurement() (or set it directly for a "
                    f"manual det_acq_series() run)."
                )
            return run[name]

        if name in PERSISTENT_FIELDS:
            return object.__getattribute__(self, "_persistent")[name]

        if name in PV_FIELDS:
            return self._read_pv(name)

        raise AttributeError(
            f"{name!r} is not an experiment setting. Static keys come from "
            f"experiment.yml, run keys from measurement_info.yaml, "
            f"{sorted(PERSISTENT_FIELDS)} are persistent session state, and "
            f"{sorted(PV_FIELDS)} live in pv_registers."
        )

    def __setattr__(self, name, value):
        # The same field tables as __getattr__, tried here in the order
        # RUN -> PERSISTENT -> PV -> STATIC; the first one that owns `name`
        # decides where the value is stored. So setting `expt.workflow_name`
        # writes a per-measurement override into the run state (RUN_FIELDS is
        # checked first) rather than hitting the STATIC_FIELDS refusal below.
        #
        # A name in none of the tables becomes an ordinary instance attribute,
        # as it would on any object.
        if name in RUN_FIELDS:
            self._run[name] = _coerce(value, RUN_FIELDS[name])
            self.dump()
            return

        if name in PERSISTENT_FIELDS:
            self._persistent[name] = _coerce(value, PERSISTENT_FIELDS[name])
            self.dump()
            return

        if name in PV_FIELDS:
            self._write_pv(name, value)
            return

        if name in STATIC_FIELDS:
            raise AttributeError(
                f"{name!r} is a per-experiment setting: edit {EXPERIMENT_FILE} "
                f"and call expt.reload(), rather than setting it here."
            )

        object.__setattr__(self, name, value)

    def __repr__(self):
        return (
            f"<ExperimentConfig {self._static.get('experiment_name')!r} "
            f"cycle={self._static.get('cycle_name')!r} "
            f"run_keys={len(self._run)}>"
        )


expt = ExperimentConfig()
