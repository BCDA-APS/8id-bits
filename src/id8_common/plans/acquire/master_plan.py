import copy
import yaml
from pathlib import Path

from apsbits.core.instrument_init import oregistry

from id8_common.plans.acquire.ad_acq import ACQ_MODES
from id8_common.plans.acquire.ad_acq import det_acq_series
from id8_common.plans.set.shutter_att import att


pv_registers = oregistry["pv_registers"]

USER_PLAN_DIR = Path("/home/beams10/8IDIUSER/bluesky/src/user_plans")
SAMPLE_INFO_FILE = USER_PLAN_DIR / "sample_info.yaml"
MEASUREMENT_INFO_FILE = USER_PLAN_DIR / "measurement_info.yaml"


# =============================================================================
# Basic file and object helpers
# =============================================================================

def read_yaml(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


def get_ophyd_object(name):
    parts = name.split(".")
    obj = oregistry[parts[0]]

    for part in parts[1:]:
        obj = getattr(obj, part)

    return obj


def get_sample_position_register(sample_index):
    register_name = f"sample{sample_index}_pos"
    return getattr(pv_registers, register_name)


def yes_no(value, field_name):
    if value is True:
        return "yes"

    if value is False:
        return "no"

    if isinstance(value, str):
        value = value.strip().lower()

        if value == "yes":
            return "yes"

        if value == "no":
            return "no"

    raise ValueError(f"{field_name} must be yes or no.")


# =============================================================================
# YAML expansion
# =============================================================================

def get_sample(sample_info, sample_index):
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
    if "samples" in run_block:
        return [int(x) for x in run_block["samples"]]

    if "sample_range" in run_block:
        first = int(run_block["sample_range"][0])
        last = int(run_block["sample_range"][1])
        return list(range(first, last + 1))

    raise ValueError("Each run block needs either samples or sample_range.")


def expand_protocols(run_block):
    protocols = run_block["protocols"]

    if isinstance(protocols, str):
        return [protocols]

    return list(protocols)


def add_measurement(expanded, protocols, run_block, sample_index, protocol_name, repeat_index):
    measurement = copy.deepcopy(protocols[protocol_name])

    measurement["sample_index"] = int(sample_index)
    measurement["protocol_name"] = protocol_name
    measurement["run_name"] = run_block.get("name", "unnamed_run")
    measurement["run_repeat"] = int(repeat_index)

    if repeat_index > 1:
        measurement["position_reset"] = "no"

    expanded.append(measurement)


def check_duplicate_assignments(run_block, seen):
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
    measurement["sample_move"] = yes_no(
        measurement["sample_move"],
        "sample_move",
    )

    measurement["position_reset"] = yes_no(
        measurement.get("position_reset", "no"),
        "position_reset",
    )


def validate_detector_mode(measurement):
    detector = measurement["detector"]
    mode = measurement["mode"]

    if detector not in ACQ_MODES:
        raise ValueError(f"Invalid detector: {detector}")

    if mode not in ACQ_MODES[detector]:
        raise ValueError(f"Invalid mode '{mode}' for detector '{detector}'.")


def validate_required_devices_connected(measurement):
    detector = measurement["detector"]
    mode = measurement["mode"]
    mode_info = ACQ_MODES[detector][mode]

    for device_name in mode_info["required_devices"]:
        device = oregistry[device_name]

        if not device.connected:
            raise RuntimeError(f"{device_name} is not connected.")


def validate_timing(measurement):
    detector = measurement["detector"]
    mode = measurement["mode"]
    mode_info = ACQ_MODES[detector][mode]

    acq_time = float(measurement["acq_time"])

    if acq_time <= 0:
        raise ValueError("acq_time must be > 0.")

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

    else:
        measurement["acq_period"] = acq_time


def validate_counts(measurement):
    num_frames = int(measurement["num_frames"])
    num_repeats = int(measurement["num_repeats"])

    if num_frames < 1:
        raise ValueError("num_frames must be >= 1.")

    if num_repeats < 1:
        raise ValueError("num_repeats must be >= 1.")


def validate_sample_motion(measurement, sample):
    sample_move = measurement["sample_move"]

    if sample_move == "no":
        return

    required_sample_fields = [
        "inner_motor",
        "outer_motor",
        "inner_center",
        "outer_center",
        "inner_range",
        "outer_range",
        "inner_pts",
        "outer_pts",
    ]

    for field in required_sample_fields:
        if field not in sample:
            raise ValueError(f"Missing sample field: {field}")

    get_ophyd_object(sample["inner_motor"])
    get_ophyd_object(sample["outer_motor"])

    inner_pts = int(sample["inner_pts"])
    outer_pts = int(sample["outer_pts"])

    if inner_pts < 1:
        raise ValueError("inner_pts must be >= 1.")

    if outer_pts < 1:
        raise ValueError("outer_pts must be >= 1.")

    sample_index = int(measurement["sample_index"])
    get_sample_position_register(sample_index)


def validate_measurement(measurement, sample):
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

    for field in required_measurement_fields:
        if field not in measurement:
            raise ValueError(f"Missing measurement field: {field}")

    normalize_measurement(measurement)

    required_sample_fields = [
        "sample_name",
        "header",
    ]

    for field in required_sample_fields:
        if field not in sample:
            raise ValueError(f"Missing sample field: {field}")

    validate_detector_mode(measurement)
    validate_required_devices_connected(measurement)
    validate_timing(measurement)
    validate_counts(measurement)
    validate_sample_motion(measurement, sample)


# =============================================================================
# Register writing
# =============================================================================

def write_sample_registers(sample_index, sample):
    pv_registers.sample_index.put(int(sample_index))
    pv_registers.header.put(sample["header"])
    pv_registers.sample_name.put(sample["sample_name"])

    if "inner_motor" in sample:
        pv_registers.inner_motor.put(sample["inner_motor"])

    if "outer_motor" in sample:
        pv_registers.outer_motor.put(sample["outer_motor"])

    if "inner_center" in sample:
        pv_registers.inner_center.put(float(sample["inner_center"]))

    if "outer_center" in sample:
        pv_registers.outer_center.put(float(sample["outer_center"]))

    if "inner_range" in sample:
        pv_registers.inner_range.put(float(sample["inner_range"]))

    if "outer_range" in sample:
        pv_registers.outer_range.put(float(sample["outer_range"]))

    if "inner_pts" in sample:
        pv_registers.inner_pts.put(int(sample["inner_pts"]))

    if "outer_pts" in sample:
        pv_registers.outer_pts.put(int(sample["outer_pts"]))


def write_measurement_registers(measurement):
    pv_registers.det_name.put(measurement["detector"])
    pv_registers.det_mode.put(measurement["mode"])

    pv_registers.acq_time.put(float(measurement["acq_time"]))
    pv_registers.acq_period.put(float(measurement["acq_period"]))

    pv_registers.num_frames.put(int(measurement["num_frames"]))
    pv_registers.num_repeats.put(int(measurement["num_repeats"]))

    pv_registers.sample_move.put(measurement["sample_move"])
    pv_registers.qmap_file.put(measurement["qmap_file"])


def reset_sample_position_register(measurement):
    if measurement["sample_move"] != "yes":
        return

    if measurement["position_reset"] != "yes":
        return

    sample_index = int(measurement["sample_index"])
    sample_position_register = get_sample_position_register(sample_index)

    sample_position_register.put(-1)


# =============================================================================
# Run functions
# =============================================================================

def run_measurement(measurement, sample_info):
    sample_index = int(measurement["sample_index"])
    sample = get_sample(sample_info, sample_index)

    validate_measurement(measurement, sample)

    write_sample_registers(sample_index, sample)
    write_measurement_registers(measurement)
    reset_sample_position_register(measurement)

    att_level = int(measurement["att_level"])
    att(att_level)

    wait_time = float(measurement.get("wait_time", 0))

    print("")
    print("==============================================")
    print(f"Run name:       {measurement.get('run_name', '')}")
    print(f"Protocol:       {measurement.get('protocol_name', '')}")
    print(f"Run repeat:     {measurement.get('run_repeat', 1)}")
    print(f"Sample index:   {sample_index}")
    print(f"Sample name:    {pv_registers.sample_name.get()}")
    print(f"Detector:       {pv_registers.det_name.get()}")
    print(f"Mode:           {pv_registers.det_mode.get()}")
    print(f"Attenuation:    {att_level}")
    print(f"acq_time:       {pv_registers.acq_time.get()}")
    print(f"acq_period:     {pv_registers.acq_period.get()}")
    print(f"num_frames:     {pv_registers.num_frames.get()}")
    print(f"num_repeats:    {pv_registers.num_repeats.get()}")
    print(f"sample_move:    {pv_registers.sample_move.get()}")
    print(f"position_reset: {measurement.get('position_reset', 'No')}")
    print(f"qmap_file:      {pv_registers.qmap_file.get()}")
    print("==============================================")
    print("")

    det_acq_series(wait_time=wait_time)


def run_measurement_info(
    measurement_info_file=MEASUREMENT_INFO_FILE,
    sample_info_file=SAMPLE_INFO_FILE,
):
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

def dry_run_measurement_info():
    sample_info = read_yaml(SAMPLE_INFO_FILE)
    measurement_info = read_yaml(MEASUREMENT_INFO_FILE)

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
        validate_timing(measurement)
        validate_counts(measurement)

        acq_period = float(measurement["acq_period"])
        num_frames = int(measurement["num_frames"])
        num_repeats = int(measurement["num_repeats"])
        wait_time = float(measurement.get("wait_time", 0))
        est_time = (acq_period * num_frames + wait_time) * num_repeats
        total_time += est_time

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
        print(f"num_repeats:    {num_repeats}")
        print(f"sample_move:    {measurement['sample_move']}")
        print(f"position_reset: {measurement.get('position_reset', 'no')}")
        print(f"qmap_file:      {measurement['qmap_file']}")
        print(f"Est. acq time:  {est_time:.1f} s")
        print("==============================================")
        print("")

    print(f"Total estimated acquisition time: {total_time:.1f} s ({total_time / 60:.1f} min)")
    print("")


# =============================================================================
# Usage examples
# =============================================================================

# Example 1:
# Run the default sample_info.yaml and measurement_info.yaml files:
#
# from id8_common.plans.acquire.master_plan import run_measurement_info
# run_measurement_info()


# Example 2:
# Run a specific measurement YAML file with the default sample YAML file:
#
# from pathlib import Path
# from id8_common.plans.acquire.master_plan import run_measurement_info
#
# measurement_file = Path("/home/beams10/8IDIUSER/bluesky/src/user_plans/measurement_info_test.yaml")
#
# run_measurement_info(
#     measurement_info_file=measurement_file,
# )


# Example 3:
# Run a specific sample YAML file and a specific measurement YAML file:
#
# from pathlib import Path
# from id8_common.plans.acquire.master_plan import run_measurement_info
#
# sample_file = Path("/home/beams10/8IDIUSER/bluesky/src/user_plans/sample_info_test.yaml")
# measurement_file = Path("/home/beams10/8IDIUSER/bluesky/src/user_plans/measurement_info_test.yaml")
#
# run_measurement_info(
#     measurement_info_file=measurement_file,
#     sample_info_file=sample_file,
# )


# Example 4:
# Recommended user workflow:
#
# 1. Edit sample_info.yaml.
# 2. Edit measurement_info.yaml.
# 3. Control repetition using runs[].repeats inside measurement_info.yaml.
# 4. Start IPython/Bluesky and run:
#
#       from id8_common.plans.acquire.master_plan import run_measurement_info
#       run_measurement_info()