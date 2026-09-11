import copy
from pathlib import Path

from id8_common.plans.acquire.ad_acq import ACQ_MODES
from id8_common.plans.acquire.ad_acq import det_acq_series
from id8_common.plans.set.shutter_att import att
from id8_common.plans.set.select_device import AXIS_NAMES
from id8_common.plans.set.select_device import DETECTOR_ALIASES
from id8_common.plans.set.select_device import _detector_config
from id8_common.plans.set.select_device import _find_motor
from id8_common.plans.set.select_device import _load_config
from id8_common.plans.set.select_device import move_detector_axes
from id8_common.plans.set.select_device import select_device
from id8_common.plans.acquire.validators import VALID_ANALYSIS_TYPES
from id8_common.plans.acquire.validators import normalize_yes_no
from id8_common.plans.acquire.validators import read_yaml
from id8_common.plans.acquire.validators import require_fields
from id8_common.plans.acquire.validators import require_mode_devices
from id8_common.plans.acquire.validators import require_positive_int
from id8_common.plans.acquire.validators import reset_sample_position
from id8_common.plans.acquire.validators import validate_acq_time
from id8_common.plans.acquire.validators import validate_qmap_exists
from id8_common.plans.acquire.validators import yes_no
from id8_common.plans.acquire import validators
from id8_common.expt_config import expt

# Not used here. Both names reached the interactive prompt through this
# module's `import *` before they moved to registry.py; the __all__ at the
# bottom of this file now excludes them, so they no longer do. Each still gets
# to the prompt by its own route: `get_ophyd_object` is in acq_helpers.__all__
# and ad_acq star-imports acq_helpers, and startup.py binds `oregistry`
# explicitly (startup_ophyd.py imports it from registry.py directly).
from id8_common.registry import get_ophyd_object  # noqa: F401
from id8_common.registry import oregistry  # noqa: F401

#: Seconds the fast shutter needs to move. In eiger4M External Series the one
#: softglue pulse both starts a segment and holds the shutter open, so BOTH the
#: pulse width (shutter open = num_frames * acq_period) and the gap to the next
#: pulse (shutter closed = trigger_period - that) are shutter movements, and
#: neither may be shorter than this.
SHUTTER_MIN_TIME = 0.1

#: Slack so a value meant to be exactly SHUTTER_MIN_TIME is not rejected by
#: binary rounding -- e.g. 0.3 - 0.2 evaluates to 0.09999999999999998.
_TIME_EPS = 1e-9



# Plan files are resolved from configs/experiment.yml at CALL time, not import
# time -- see expt.user_plan_dir. They used to be three hardcoded absolute paths
# here, which broke the moment the files moved under <cycle>/<experiment>/ and
# would have broken again on any clone of this repo. Nothing is pinned now:
# change cycle_name/experiment_name in experiment.yml and the plans follow.


# read_yaml, yes_no and get_ophyd_object used to be defined here. They are now
# imported above -- from validators.py and registry.py -- so the trio path and
# this one cannot drift apart again. None of the three is in this module's
# __all__, so `import *` in startup.py no longer carries them to the prompt;
# import them by name from validators.py / registry.py instead.


# =============================================================================
# YAML expansion
# =============================================================================
# A `runs:` block in measurement_info.yaml is shorthand. expand_measurements()
# turns the whole file into one flat list of "measurement" dicts -- one dict per
# acquisition -- and run_measurement_info() then walks that list in order.
#
# Worked example. This block:
#
#     loop_order: sample_major
#     runs:
#       - name: mode_test
#         samples: [1, 2]
#         protocols: [eiger_internal_series, rigaku_epics]
#         repeats: 2
#
# expands to 2 samples x 2 protocols x 2 repeats = 8 measurements, in this order
# (sample_major means the sample index turns slower than the protocol name;
# protocol_major just swaps those two loops):
#
#     repeat 1:  (1, eiger_internal_series)  (1, rigaku_epics)
#                (2, eiger_internal_series)  (2, rigaku_epics)
#     repeat 2:  the same four again
#
# Each measurement is a deep copy of that protocol's own YAML body -- detector,
# mode, acq_time, num_frames and the rest -- plus four bookkeeping keys that
# add_measurement() stamps on it, e.g. for the last one above:
#
#     {..., "sample_index": 2, "protocol_name": "rigaku_epics",
#           "run_name": "mode_test", "run_repeat": 2}
#
# Call chain: expand_measurements -> expand_run_block (once per runs: entry)
#             -> expand_samples / expand_protocols -> add_measurement.


def get_sample(sample_info, sample_index):
    """One sample's settings from sample_info.yaml: the file's `defaults`, overridden by `sample_N`."""
    sample_key = f"sample_{sample_index}"

    defaults = sample_info.get("defaults", {})
    samples = sample_info["samples"]

    if sample_key not in samples:
        raise ValueError(f"{sample_key} is not defined in sample_info.yaml.")

    sample = {}
    sample.update(defaults)
    sample.update(samples[sample_key])

    return sample


def expand_samples(run_block):
    """Sample indices one run block covers: an explicit `samples:` list, or an inclusive `sample_range:`."""
    if "samples" in run_block:
        return [int(x) for x in run_block["samples"]]

    if "sample_range" in run_block:
        first = int(run_block["sample_range"][0])
        last = int(run_block["sample_range"][1])
        return list(range(first, last + 1))

    raise ValueError("Each run block needs either samples or sample_range.")


def expand_protocols(run_block):
    """Protocol names one run block covers. A bare string counts as a one-item list."""
    protocols = run_block["protocols"]

    if isinstance(protocols, str):
        return [protocols]

    return list(protocols)


def add_measurement(expanded, protocols, run_block, sample_index, protocol_name, repeat_index):
    """Append one measurement -- a copy of the protocol body plus its run bookkeeping -- to `expanded`.

    The copy is deep because every measurement is edited in place later
    (validate_timing fills in acq_period, validate_counts fills in
    num_segments), and the protocols dict is shared by every run block that
    names it.
    """
    measurement = copy.deepcopy(protocols[protocol_name])

    measurement["sample_index"] = int(sample_index)
    measurement["protocol_name"] = protocol_name
    measurement["run_name"] = run_block.get("name", "unnamed_run")
    measurement["run_repeat"] = int(repeat_index)

    # Only the first repeat may rewind the sample mesh. position_reset: yes
    # sends the mesh index back to -1 (see validators.reset_sample_position);
    # doing that on repeats 2..N would walk later repeats back over the same
    # spots the first one already exposed.
    if repeat_index > 1:
        measurement["position_reset"] = "no"

    expanded.append(measurement)


def check_duplicate_assignments(run_block, seen):
    """Reject a (sample, protocol) pair that an earlier run block in the same file already claimed.

    `seen` is carried across every run block, so the check is file-wide rather
    than per-block. Measuring the same pair twice is nearly always a copy-paste
    slip; asking for it on purpose is what runs[].repeats is for.
    """
    samples = expand_samples(run_block)
    protocol_names = expand_protocols(run_block)

    for sample_index in samples:
        for protocol_name in protocol_names:
            key = (sample_index, protocol_name)

            if key in seen:
                raise ValueError(
                    f"Duplicate measurement assignment: "
                    f"sample {sample_index}, protocol {protocol_name}. "
                    f"Use runs[].repeats instead of duplicating run blocks."
                )

            seen.add(key)


def expand_run_block(measurement_info, run_block):
    """Expand one `runs:` entry into its list of measurement dicts -- see the worked example above."""
    protocols = measurement_info["protocols"]
    global_loop_order = measurement_info.get("loop_order", "sample_major")

    samples = expand_samples(run_block)
    protocol_names = expand_protocols(run_block)

    block_repeats = int(run_block.get("repeats", 1))
    loop_order = run_block.get("loop_order", global_loop_order)

    if block_repeats < 1:
        raise ValueError("runs[].repeats must be >= 1.")

    for protocol_name in protocol_names:
        if protocol_name not in protocols:
            raise ValueError(f"Protocol '{protocol_name}' is not defined.")

    expanded = []

    if loop_order == "sample_major":
        for repeat_index in range(1, block_repeats + 1):
            for sample_index in samples:
                for protocol_name in protocol_names:
                    add_measurement(
                        expanded=expanded,
                        protocols=protocols,
                        run_block=run_block,
                        sample_index=sample_index,
                        protocol_name=protocol_name,
                        repeat_index=repeat_index,
                    )

    elif loop_order == "protocol_major":
        for repeat_index in range(1, block_repeats + 1):
            for protocol_name in protocol_names:
                for sample_index in samples:
                    add_measurement(
                        expanded=expanded,
                        protocols=protocols,
                        run_block=run_block,
                        sample_index=sample_index,
                        protocol_name=protocol_name,
                        repeat_index=repeat_index,
                    )

    else:
        raise ValueError("loop_order must be sample_major or protocol_major.")

    return expanded


def expand_measurements(measurement_info):
    """Flatten every `runs:` block in the file into one ordered list of measurement dicts."""
    runs = measurement_info["runs"]

    expanded = []
    seen = set()

    for run_block in runs:
        check_duplicate_assignments(run_block, seen)
        expanded.extend(expand_run_block(measurement_info, run_block))

    return expanded


# =============================================================================
# Validation
# =============================================================================

def normalize_measurement(measurement):
    normalize_yes_no(measurement)


def validate_detector_mode(measurement):
    validators.validate_detector_mode(measurement["detector"], measurement["mode"])


def validate_required_devices_connected(measurement):
    # Also covers hardware_device, which this check used to miss -- see
    # validators.require_mode_devices.
    require_mode_devices(measurement["detector"], measurement["mode"])


def validate_timing(measurement):
    """Check acq_time / acq_period / trigger_period against this detector mode's rules.

    Edits `measurement` in place: a mode that paces its own frames states no
    acq_period, and this fills one in equal to acq_time so everything
    downstream can read the key unconditionally.
    """
    detector = measurement["detector"]
    mode = measurement["mode"]
    mode_info = ACQ_MODES[detector][mode]

    acq_time = validate_acq_time(measurement["acq_time"], detector, mode)

    if mode_info["needs_acq_period"]:
        if "acq_period" not in measurement:
            raise ValueError(f"{detector} {mode} requires acq_period.")

        acq_period = float(measurement["acq_period"])

        if acq_period <= 0:
            raise ValueError("acq_period must be > 0.")

        if acq_period < acq_time:
            raise ValueError("acq_period must be >= acq_time.")

        if detector == "eiger4M" and mode == "External Enable":
            if acq_time < 0.1:
                raise ValueError("Eiger External Enable requires acq_time >= 0.1 s.")

            if acq_period < 0.1:
                raise ValueError("Eiger External Enable requires acq_period >= 0.1 s.")

        # External Series has no 0.1 s floor. That floor exists because the
        # softglue pulse drives the shutter, and a shutter cannot follow faster
        # than that. In External Series softglue only starts each segment --
        # the shutter stays open for the whole acquisition -- and acq_period is
        # the Eiger pacing its own frames internally, so nothing mechanical
        # limits it. Only trigger_period, which is still softglue-generated,
        # has a constraint (see below).

        if detector == "lambda2M" and mode == "External":
            if acq_time < 0.1:
                raise ValueError("Lambda External requires acq_time >= 0.1 s.")

            if acq_period < 0.1:
                raise ValueError("Lambda External requires acq_period >= 0.1 s.")

    else:
        measurement["acq_period"] = acq_time

    if mode_info.get("needs_trigger_period", False):
        if "trigger_period" not in measurement:
            raise ValueError(f"{detector} {mode} requires trigger_period.")

        trigger_period = float(measurement["trigger_period"])

        if trigger_period <= 0:
            raise ValueError("trigger_period must be > 0.")

        # The single softglue pulse does two jobs: its rising edge starts a
        # segment, and its width holds the shutter open. So the acquisition has
        # two shutter movements per segment -- open for num_frames * acq_period,
        # then closed for the remainder of trigger_period -- and each needs at
        # least SHUTTER_MIN_TIME. The closed-gap check also subsumes the older
        # "trigger_period must exceed one segment" rule: without it a pulse
        # would land while the detector was still busy, be dropped, and the
        # acquisition would hang on a trigger that was already spent.
        segment_time = float(measurement["acq_period"]) * int(measurement["num_frames"])
        closed_time = trigger_period - segment_time

        if segment_time < SHUTTER_MIN_TIME - _TIME_EPS:
            raise ValueError(
                f"{detector} {mode}: the shutter is held open for one segment, "
                f"num_frames * acq_period = {measurement['num_frames']} * "
                f"{measurement['acq_period']} = {segment_time:g} s, which is shorter "
                f"than the {SHUTTER_MIN_TIME} s the shutter needs to move."
            )

        if closed_time < SHUTTER_MIN_TIME - _TIME_EPS:
            raise ValueError(
                f"{detector} {mode}: the shutter is closed between segments for "
                f"trigger_period - num_frames * acq_period = {trigger_period:g} - "
                f"{segment_time:g} = {closed_time:g} s, which is shorter than the "
                f"{SHUTTER_MIN_TIME} s the shutter needs to move. Raise trigger_period "
                f"to at least {segment_time + SHUTTER_MIN_TIME:g} s."
            )

    elif "trigger_period" in measurement:
        raise ValueError(f"{detector} {mode} does not use trigger_period.")


def validate_detector_position_overrides(measurement):
    """Refuse a protocol that asks to move a detector axis this detector does not have.

    An override is any AXIS_NAMES key written straight into the protocol body --
    `horizontal`, `vertical`, `swing_angle_horizontal`, `swing_angle_vertical`.
    run_measurement() feeds those to move_detector_axes() after select_device(),
    so they override the position device_position.yaml would have parked at.
    """
    overrides = {axis: measurement[axis] for axis in AXIS_NAMES if axis in measurement}

    if not overrides:
        return

    config = _load_config()
    cfg = _detector_config(config, measurement["detector"])
    motors_cfg = cfg["motors"]

    for axis_name in overrides:
        motor = _find_motor(motors_cfg, axis_name)

        if motor.get("device") is None:
            raise ValueError(
                f"Protocol '{measurement.get('protocol_name', '')}': "
                f"'{measurement['detector']}' has no '{axis_name}' axis to move."
            )


def validate_analysis_type(measurement):
    validators.validate_analysis_type(measurement.get("analysis_type", "Multitau"))


def validate_counts(measurement):
    """Check the frame/repeat/segment counts. Writes num_segments back as an int for modes that use it."""
    detector = measurement["detector"]
    mode = measurement["mode"]
    mode_info = ACQ_MODES[detector][mode]

    require_positive_int(measurement["num_frames"], "num_frames")
    require_positive_int(measurement["num_repeats"], "num_repeats")

    if mode_info.get("needs_num_segments", False):
        measurement["num_segments"] = require_positive_int(
            measurement.get("num_segments", 1), "num_segments"
        )

    elif "num_segments" in measurement:
        # Reject rather than ignore: a num_segments on a mode that never reads
        # it would silently do nothing.
        raise ValueError(f"{detector} {mode} does not use num_segments.")


def validate_sample_motion(measurement, sample):
    # No forbidden_motors: the serial path may mesh on any axis. The trio path
    # passes FORBIDDEN_MOTORS to the same check.
    validators.validate_sample_motion(measurement, sample)


def validate_measurement(measurement, sample):
    """Run every check a measurement must pass before any hardware moves.

    Also NORMALISES `measurement` in place along the way -- yes/no fields to
    their string form, a missing acq_period filled in, num_segments coerced to
    int -- so the callers below can read those keys without re-deriving them.
    """
    required_measurement_fields = [
        "sample_index",
        "detector",
        "mode",
        "att_level",
        "acq_time",
        "num_frames",
        "num_repeats",
        "sample_move",
        "qmap_file",
    ]

    require_fields(measurement, required_measurement_fields, "measurement")

    normalize_measurement(measurement)

    required_sample_fields = [
        "sample_name",
        "header",
    ]

    require_fields(sample, required_sample_fields, "sample")

    validate_detector_mode(measurement)
    validate_detector_position_overrides(measurement)
    validate_analysis_type(measurement)
    validate_required_devices_connected(measurement)
    validate_timing(measurement)
    validate_counts(measurement)
    validate_qmap_exists(measurement["qmap_file"])
    validate_sample_motion(measurement, sample)


def reset_sample_position_register(measurement):
    reset_sample_position(measurement)


# =============================================================================
# Detector placeholder hooks
# =============================================================================
# Called from run_measurement() right after select_device(), for every detector.
# All three hooks are no-ops as they stand: nothing here moves the huber. The
# moves are still in the file, commented out -- eiger4M to delta 10 / nu 22.7,
# rigaku3M (and its rigaku3M_epics alias) to delta 10 / nu 0, lambda2M nothing.
# They were written on the assumption that every detector acquires at delta = 10
# -- at delta = 0 the lambda2M sits in the beam. Until one is uncommented, the
# diffractometer has to be positioned before the acquisition.

def placeholder_eiger4M():
    pass
    # huber = oregistry["huber"]
    # huber.delta.move(10, wait=True)
    # huber.nu.move(22.7, wait=True)


def placeholder_lambda2M():
    pass


def placeholder_rigaku3M():
    pass
    # huber = oregistry["huber"]
    # huber.delta.move(10, wait=True)
    # huber.nu.move(0, wait=True)


DETECTOR_PLACEHOLDERS = {
    "eiger4M": placeholder_eiger4M,
    "lambda2M": placeholder_lambda2M,
    "rigaku3M": placeholder_rigaku3M,
}

# rigaku3M_epics and rigaku3M_ftf are the same physical detector as rigaku3M
# (see DETECTOR_ALIASES in select_device.py) so they get the same placeholder
# hook. Any alias added there is picked up here automatically.
for _alias, _canonical in DETECTOR_ALIASES.items():
    if _canonical in DETECTOR_PLACEHOLDERS:
        DETECTOR_PLACEHOLDERS[_alias] = DETECTOR_PLACEHOLDERS[_canonical]


def run_detector_placeholder(name: str):
    """Run the placeholder hook for a detector (rigaku3M aliases share rigaku3M's hook)."""
    placeholder = DETECTOR_PLACEHOLDERS.get(name)
    if placeholder is not None:
        placeholder()


# =============================================================================
# Run functions
# =============================================================================

def run_measurement(measurement, sample_info):
    """Validate one expanded measurement, publish it to the run state, set up hardware, and acquire."""
    sample_index = int(measurement["sample_index"])
    sample = get_sample(sample_info, sample_index)

    validate_measurement(measurement, sample)

    # The run state is what the acquisition path actually reads.
    expt.sample_index = sample_index
    expt.set_measurement(measurement=measurement, sample=sample)

    reset_sample_position_register(measurement)

    att_level = int(measurement["att_level"])
    att(att_level)

    wait_time = float(measurement.get("wait_time", 0))
    position_overrides = {axis: measurement[axis] for axis in AXIS_NAMES if axis in measurement}

    print("")
    print("==============================================")
    print(f"Run name:       {measurement.get('run_name', '')}")
    print(f"Protocol:       {measurement.get('protocol_name', '')}")
    print(f"Run repeat:     {measurement.get('run_repeat', 1)}")
    print(f"Sample index:   {sample_index}")
    print(f"Sample name:    {expt.sample_name}")
    print(f"Detector:       {expt.det_name}")
    print(f"Mode:           {expt.det_mode}")
    print(f"Attenuation:    {att_level}")
    print(f"acq_time:       {expt.acq_time}")
    print(f"acq_period:     {expt.acq_period}")
    print(f"num_frames:     {expt.num_frames}")
    if "num_segments" in measurement:
        print(f"num_segments:   {expt.num_segments}  "
              f"({expt.num_frames * expt.num_segments} frames total)")
    if "trigger_period" in measurement:
        print(f"trigger_period: {expt.trigger_period} s")
    print(f"num_repeats:    {expt.num_repeats}")
    print(f"sample_move:    {expt.sample_move}")
    print(f"position_reset: {measurement.get('position_reset', 'No')}")
    print(f"qmap_file:      {expt.qmap_file}")
    print(f"Analysis type:  {expt.analysis_type}")
    if position_overrides:
        print(f"Position override: {position_overrides}")
    print("==============================================")
    print("")

    select_device(measurement["detector"])
    run_detector_placeholder(measurement["detector"])

    if position_overrides:
        move_detector_axes(measurement["detector"], position_overrides)

    det_acq_series(wait_time=wait_time, hooks=measurement.get("hooks"))


def reload_experiment_config():
    """Re-read configs/experiment.yml, and say so if anything changed.

    expt caches experiment.yml: it is read once at startup and then only by an
    explicit expt.reload(). Everything downstream -- which cycle and experiment
    the data lands under, which mount point, which analysis machine, which DM
    workflow -- comes from that cache, so editing the file mid-session changed
    nothing until someone remembered to reload, and a run could quietly write
    into the previous experiment's tree.

    Called at the top of the run_* and dry_run_* entry points, before any path is
    resolved from expt, so the plan files themselves also follow an edit. Only the
    static block is re-read; run state and persistent state are untouched (see
    ExperimentConfig.reload).

    A change is announced, and experiment_name/mount_point/cycle_name are called
    out as a group, because those three decide where the data goes.
    """
    before = dict(getattr(expt, "_static", {}) or {})
    after = expt.reload()

    changed = {k: (before.get(k), v) for k, v in after.items() if before.get(k) != v}
    if not changed:
        return changed

    print("experiment.yml changed since it was last read:")
    for key in sorted(changed):
        was, now = changed[key]
        print(f"    {key}: {was!r} -> {now!r}")
    if {"experiment_name", "mount_point", "cycle_name"} & set(changed):
        print(f"    -> data now goes to {expt.mount_point}{expt.cycle_name}/"
              f"{expt.experiment_name}/")
    return changed


def run_measurement_info(
    measurement_info_file=None,
    sample_info_file=None,
):
    """Read measurement_info.yaml + sample_info.yaml, expand them, and run every measurement in order."""
    # None, not a default argument: a default is bound once at import, so it
    # could not follow an experiment.yml edit or expt.reload().
    # Path(), because the caller may pass a plain string -- the usage examples
    # and the trio module's docstring both show one -- and .parent below is a
    # Path method. expt.measurement_info_file is already a Path; Path() on a
    # Path is a no-op.
    # experiment.yml is cached by expt; re-read it so an edit takes effect
    # without a restart, and BEFORE the paths below are resolved from it.
    reload_experiment_config()

    measurement_info_file = Path(measurement_info_file or expt.measurement_info_file)
    sample_info_file = Path(sample_info_file or expt.sample_info_file)

    print(f"Reading plans from {measurement_info_file.parent}")

    sample_info = read_yaml(sample_info_file)
    measurement_info = read_yaml(measurement_info_file)

    measurements = expand_measurements(measurement_info)

    print("")
    print(f"Loaded {len(measurements)} expanded measurement blocks.")
    print("")

    for measurement in measurements:
        run_measurement(
            measurement=measurement,
            sample_info=sample_info,
        )


# =============================================================================
# Dry-run preview (no acquisitions executed)
# =============================================================================

def dry_run_measurement_info(measurement_info_file=None, sample_info_file=None):
    """Validate and print what run_measurement_info() would do, with a time estimate. Moves nothing.

    Deliberately NOT the full set of checks: this repeats the detector/mode,
    position-override, analysis-type, timing and count checks, but skips the
    required-field check, the device-connected check and the sample-mesh check
    that validate_measurement() also makes. A protocol that passes here can
    still be rejected once run_measurement() gets to it.
    """
    # experiment.yml is cached by expt; re-read it so an edit takes effect
    # without a restart, and BEFORE the paths below are resolved from it.
    reload_experiment_config()
    measurement_info_file = measurement_info_file or expt.measurement_info_file
    sample_info_file = sample_info_file or expt.sample_info_file

    sample_info = read_yaml(sample_info_file)
    measurement_info = read_yaml(measurement_info_file)

    measurements = expand_measurements(measurement_info)

    print("")
    print(f"Total measurements planned: {len(measurements)}")
    print("")

    total_time = 0.0

    for measurement in measurements:
        sample_index = int(measurement["sample_index"])
        sample = get_sample(sample_info, sample_index)

        normalize_measurement(measurement)
        validate_detector_mode(measurement)
        validate_detector_position_overrides(measurement)
        validate_analysis_type(measurement)
        validate_timing(measurement)
        validate_counts(measurement)
        validate_qmap_exists(measurement["qmap_file"])

        acq_period = float(measurement["acq_period"])
        num_frames = int(measurement["num_frames"])
        num_repeats = int(measurement["num_repeats"])
        # num_frames is per segment where a mode uses segments, so the frames
        # actually collected in one acquisition is the product.
        num_segments = int(measurement.get("num_segments", 1))
        wait_time = float(measurement.get("wait_time", 0))
        # Where segments are triggered externally the wall-clock cost is set by
        # the pulse spacing, not by how long a segment is busy for.
        trigger_period = float(measurement.get("trigger_period", 0))
        if trigger_period > 0:
            est_time = (trigger_period * num_segments + wait_time) * num_repeats
        else:
            est_time = (acq_period * num_frames * num_segments + wait_time) * num_repeats
        total_time += est_time

        position_overrides = {axis: measurement[axis] for axis in AXIS_NAMES if axis in measurement}

        print("")
        print("==============================================")
        print(f"Run name:       {measurement.get('run_name', '')}")
        print(f"Protocol:       {measurement.get('protocol_name', '')}")
        print(f"Run repeat:     {measurement.get('run_repeat', 1)}")
        print(f"Sample index:   {sample_index}")
        print(f"Sample name:    {sample.get('sample_name', '')}")
        print(f"Detector:       {measurement['detector']}")
        print(f"Mode:           {measurement['mode']}")
        print(f"Attenuation:    {measurement['att_level']}")
        print(f"acq_time:       {measurement['acq_time']}")
        print(f"acq_period:     {measurement['acq_period']}")
        print(f"num_frames:     {num_frames}")
        if "num_segments" in measurement:
            print(f"num_segments:   {num_segments}  ({num_frames * num_segments} frames total)")
        if "trigger_period" in measurement:
            print(f"trigger_period: {trigger_period:g} s  (segment busy {acq_period * num_frames:g} s)")
        print(f"num_repeats:    {num_repeats}")
        print(f"sample_move:    {measurement['sample_move']}")
        print(f"position_reset: {measurement.get('position_reset', 'no')}")
        print(f"qmap_file:      {measurement['qmap_file']}")
        print(f"Analysis type:  {measurement.get('analysis_type', 'Multitau')}")
        if position_overrides:
            print(f"Position override: {position_overrides}")
        print(f"Est. acq time:  {est_time:.1f} s")
        print("==============================================")
        print("")

    print(f"Total estimated acquisition time: {total_time:.1f} s ({total_time / 60:.1f} min)")
    print("")


# =============================================================================
# Usage examples -- see docs/running-measurements.md for
# run_measurement_info()/dry_run_measurement_info() usage and the
# recommended edit-YAML-then-run workflow.
# =============================================================================


#: What ``from master_plan import *`` puts in the beamline session.
#:
#: startup.py and startup_ophyd.py both star-import this module, so every
#: public name here lands at the scientist's prompt. Without this list that
#: was all 48 of them: ``copy``, ``Path``, the ``validators`` module, the whole
#: YAML-expansion chain, and every ``validate_*`` helper -- so tab-completing
#: ``val`` offered eight internals and ``require_positive_int`` looked like
#: something you were meant to call. It also made the prompt fragile in the
#: other direction: ``oregistry`` was reaching the session ONLY because this
#: module happened to import it, so a tidy-up of the imports here would
#: silently have taken it away. Neither it nor ``att`` is exported from here
#: now -- startup.py binds ``oregistry`` itself, and ``att`` arrives with
#: tetramm_acq's star-import of shutter_att, which startup.py runs first.
#:
#: The rule: list what a scientist would type, keep everything else internal.
#: Nothing below is required for the module to work -- ``__all__`` affects
#: ``import *`` only. trio_master_plan_rigaku3m_eiger4m_lambda2m.py's explicit
#: ``from master_plan import expand_measurements`` is unaffected, and so is
#: anything else that imports a name from here by name.
__all__ = [
    # Entry points.
    "run_measurement_info",
    "dry_run_measurement_info",
    "run_measurement",
    # Callable against a hand-built measurement dict, for debugging a protocol
    # without running it.
    "validate_measurement",
    # Limits worth reading at the prompt when a protocol is rejected.
    "SHUTTER_MIN_TIME",
    "VALID_ANALYSIS_TYPES",
    # Per-detector hook run right after select_device(), for every detector.
    "DETECTOR_PLACEHOLDERS",
    "run_detector_placeholder",
]
