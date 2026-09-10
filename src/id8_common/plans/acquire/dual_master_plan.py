"""
YAML front end for parallel two-detector acquisition.

The dual counterpart to master_plan.py. Reads user_plans/dual_measurement_info.yaml, validates
it hard enough that nothing can move before an error surfaces, and hands per-detector "legs" to
dual_acq.dual_acq_series().

Imported by startup.py, so these are already in the session:

    dry_run_dual_measurement_info(check_hardware=True)   # validate and preview, moves nothing
    run_dual_measurement_info()                          # go

Importing this alongside master_plan.py is safe: the two share no state and no names. A dual
run touches the same pv_registers only after its parallel window has closed, and restores them
afterwards, so single-detector acquisition behaves exactly as it did before.

Run expansion (runs/protocols/samples/loop_order) is reused wholesale from master_plan.py, so
the two files behave identically there. What differs is the protocol body: a dual protocol
carries a `detectors:` list, and the acquisition settings that are per-detector live inside it
rather than at protocol level.

sample_info.yaml is shared with the serial path and read unchanged.
"""

from apsbits.core.instrument_init import oregistry
from id8_common.plans.acquire.ad_acq import ACQ_MODES
from id8_common.plans.acquire.ad_acq import get_ophyd_object
from id8_common.plans.acquire.dual_acq import DUAL_LEGS
from id8_common.plans.acquire.dual_acq import FORBIDDEN_MOTORS
from id8_common.plans.acquire.dual_acq import dual_acq_series
from id8_common.plans.acquire.master_plan import SAMPLE_INFO_FILE
from id8_common.plans.acquire.master_plan import USER_PLAN_DIR
from id8_common.plans.acquire.master_plan import VALID_ANALYSIS_TYPES
from id8_common.plans.acquire.master_plan import expand_measurements
from id8_common.plans.acquire.master_plan import get_sample
from id8_common.plans.acquire.master_plan import get_sample_position_register
from id8_common.plans.acquire.master_plan import read_yaml
from id8_common.plans.acquire.master_plan import write_sample_registers
from id8_common.plans.acquire.master_plan import yes_no
from id8_common.plans.set.select_device import DETECTOR_ALIASES
from id8_common.plans.set.shutter_att import att

pv_registers = oregistry["pv_registers"]

DUAL_MEASUREMENT_INFO_FILE = USER_PLAN_DIR / "dual_measurement_info.yaml"

# Huber position a dual Eiger+Rigaku measurement acquires at.
DUAL_HUBER_DELTA = 10.0
DUAL_HUBER_NU = 0.0

REQUIRED_PROTOCOL_FIELDS = [
    "att_level",
    "num_repeats",
    "sample_move",
    "detectors",
]

REQUIRED_LEG_FIELDS = [
    "device",
    "mode",
    "acq_time",
    "num_frames",
    "qmap_file",
]

# geometry is optional, and so is every field in it. A detector at the mount
# device_position.yaml already describes needs no geometry block -- that file's values are
# correct and are used as-is. State a field here only to override it, which is what an Eiger
# remounted on the huber arm needs.
KNOWN_GEOMETRY_FIELDS = [
    "db_x",
    "db_y",
    "distance",
    "pixel_size",
    "position_x",
    "position_y",
    "beam_center_position_x",
    "beam_center_position_y",
    "swing_horizontal",
    "swing_vertical",
]


# =============================================================================
# Helpers
# =============================================================================


def as_bool(value, field_name):
    """YAML `yes`/`no` arrives as a bool; accept the string spellings too."""
    if isinstance(value, bool):
        return value

    return yes_no(value, field_name) == "yes"


def leg_label(leg):
    """Short detector name used in file paths: rigaku3M_epics -> rigaku3M."""
    if leg.get("label"):
        return str(leg["label"])

    return DETECTOR_ALIASES.get(leg["device"], leg["device"])


def leg_duration(leg):
    """Expected acquisition time for one repeat of one leg, in seconds."""
    return float(leg["acq_time"]) * int(leg["num_frames"])


# =============================================================================
# Validation
# =============================================================================


def validate_geometry(geometry, label):
    """Check an optional geometry block. Absent means 'use device_position.yaml as-is'."""
    if geometry is None:
        return

    if not isinstance(geometry, dict):
        raise ValueError(f"Leg '{label}': geometry must be a mapping.")

    for field, value in geometry.items():
        if field not in KNOWN_GEOMETRY_FIELDS:
            raise ValueError(f"Leg '{label}': unknown geometry field '{field}'. Known: {sorted(KNOWN_GEOMETRY_FIELDS)}")

        # Placeholders must never reach a data file. A dotted ophyd path is a legal value
        # (it is read live at metadata time), so only reject things that resolve to neither
        # a number nor a real device.
        if isinstance(value, str):
            if value.strip().upper() == "TBD":
                raise ValueError(
                    f"Leg '{label}': geometry field '{field}' is still TBD. "
                    f"Measure it for this detector mount before running."
                )
            try:
                get_ophyd_object(value)
            except Exception as exc:
                raise ValueError(
                    f"Leg '{label}': geometry field '{field}' = '{value}' is neither a number "
                    f"nor a resolvable ophyd path ({exc})."
                ) from exc
        else:
            float(value)


def validate_leg(leg):
    label = leg_label(leg)

    for field in REQUIRED_LEG_FIELDS:
        if field not in leg:
            raise ValueError(f"Leg '{label}': missing field '{field}'.")

    device = leg["device"]
    mode = leg["mode"]

    if device not in ACQ_MODES:
        raise ValueError(f"Leg '{label}': invalid detector '{device}'. Known: {sorted(ACQ_MODES)}")

    if mode not in ACQ_MODES[device]:
        raise ValueError(f"Leg '{label}': invalid mode '{mode}' for '{device}'. Known: {sorted(ACQ_MODES[device])}")

    if (device, mode) not in DUAL_LEGS:
        supported = ", ".join(f"{d}/{m}" for d, m in sorted(DUAL_LEGS))
        raise ValueError(
            f"Leg '{label}': {device}/{mode} is not supported for dual acquisition. Supported: {supported}"
        )

    acq_time = float(leg["acq_time"])

    if acq_time <= 0:
        raise ValueError(f"Leg '{label}': acq_time must be > 0.")

    min_acq_time = ACQ_MODES[device][mode].get("min_acq_time")

    if min_acq_time is not None and acq_time < min_acq_time:
        raise ValueError(
            f"Leg '{label}': {device} {mode} requires acq_time >= {min_acq_time:.2e} s (got {acq_time:.2e} s)."
        )

    if int(leg["num_frames"]) < 1:
        raise ValueError(f"Leg '{label}': num_frames must be >= 1.")

    if not str(leg["qmap_file"]).strip():
        raise ValueError(f"Leg '{label}': qmap_file must not be empty.")

    analysis_type = leg.get("analysis_type", "Multitau")

    if analysis_type not in VALID_ANALYSIS_TYPES:
        raise ValueError(f"Leg '{label}': analysis_type must be one of {VALID_ANALYSIS_TYPES} (got '{analysis_type}').")

    validate_geometry(leg.get("geometry"), label)

    for dotted in leg.get("motors") or {}:
        if dotted in FORBIDDEN_MOTORS:
            raise ValueError(
                f"Leg '{label}': '{dotted}' cannot appear in a dual protocol's motors block. "
                f"Both huber axes are positioned once before acquisition by "
                f"setup_huber_for_dual() (delta {DUAL_HUBER_DELTA}, nu {DUAL_HUBER_NU})."
            )
        get_ophyd_object(dotted)

    device_obj = oregistry[ACQ_MODES[device][mode].get("hardware_device", device)]

    if not device_obj.connected:
        raise RuntimeError(f"Leg '{label}': {device} is not connected.")


def validate_shutter_owner(legs):
    owners = [leg for leg in legs if as_bool(leg.get("shutter_owner", False), "shutter_owner")]

    if len(legs) == 1 and not owners:
        return

    if len(owners) != 1:
        labels = ", ".join(leg_label(leg) for leg in owners) or "none"
        raise ValueError(
            f"Exactly one leg must set shutter_owner: yes (got {len(owners)}: {labels}). "
            f"The owner is armed first and everyone else waits for it to confirm it is acquiring."
        )


def validate_sample_motion(measurement, sample):
    if measurement["sample_move"] != "yes":
        return

    required = [
        "inner_motor",
        "outer_motor",
        "inner_center",
        "outer_center",
        "inner_range",
        "outer_range",
        "inner_pts",
        "outer_pts",
    ]

    for field in required:
        if field not in sample:
            raise ValueError(f"Missing sample field: {field}")

    # The mesh is the other way huber.delta / huber.nu could be driven, via sample_info.yaml's
    # inner_motor / outer_motor. Refuse it for the same reason as a leg's motors block.
    for role in ("inner_motor", "outer_motor"):
        if sample[role] in FORBIDDEN_MOTORS:
            raise ValueError(
                f"sample_info.yaml sets {role} = '{sample[role]}', which a dual acquisition "
                f"must never move. Use a different sample axis, or set sample_move: no."
            )

    get_ophyd_object(sample["inner_motor"])
    get_ophyd_object(sample["outer_motor"])

    if int(sample["inner_pts"]) < 1:
        raise ValueError("inner_pts must be >= 1.")

    if int(sample["outer_pts"]) < 1:
        raise ValueError("outer_pts must be >= 1.")

    get_sample_position_register(int(measurement["sample_index"]))


def normalize_dual_measurement(measurement):
    measurement["sample_move"] = yes_no(measurement["sample_move"], "sample_move")
    measurement["position_reset"] = yes_no(measurement.get("position_reset", "no"), "position_reset")


def validate_dual_measurement(measurement, sample, check_hardware=True):
    for field in REQUIRED_PROTOCOL_FIELDS:
        if field not in measurement:
            raise ValueError(f"Missing protocol field: {field}")

    for field in ["sample_name", "header"]:
        if field not in sample:
            raise ValueError(f"Missing sample field: {field}")

    normalize_dual_measurement(measurement)

    legs = measurement["detectors"]

    if not isinstance(legs, list) or not legs:
        raise ValueError("protocol 'detectors' must be a non-empty list.")

    if int(measurement["num_repeats"]) < 1:
        raise ValueError("num_repeats must be >= 1.")

    labels = [leg_label(leg) for leg in legs]

    if len(set(labels)) != len(labels):
        raise ValueError(f"Duplicate detector labels in one protocol: {labels}. Give one of them an explicit label.")

    validate_shutter_owner(legs)

    if check_hardware:
        for leg in legs:
            validate_leg(leg)

    validate_sample_motion(measurement, sample)


# =============================================================================
# Leg construction
# =============================================================================


def build_leg_specs(measurement):
    """Turn validated YAML detector blocks into the dicts dual_acq_series() consumes."""
    specs = []

    for leg in measurement["detectors"]:
        spec = {
            "device": leg["device"],
            "mode": leg["mode"],
            "label": leg_label(leg),
            "acq_time": float(leg["acq_time"]),
            "num_frames": int(leg["num_frames"]),
            "qmap_file": str(leg["qmap_file"]),
            "analysis_type": leg.get("analysis_type", "Multitau"),
            "geometry": leg.get("geometry") or {},
            "motors": leg.get("motors") or {},
            "select_device": as_bool(leg.get("select_device", False), "select_device"),
            "shutter_owner": as_bool(leg.get("shutter_owner", False), "shutter_owner"),
        }

        if "start_timeout" in leg:
            spec["start_timeout"] = float(leg["start_timeout"])

        if "hdf_timeout" in leg:
            spec["hdf_timeout"] = float(leg["hdf_timeout"])

        if leg.get("workflow_name"):
            spec["workflow_name"] = str(leg["workflow_name"])

        specs.append(spec)

    return specs


def reset_sample_position_register(measurement):
    if measurement["sample_move"] != "yes":
        return

    if measurement["position_reset"] != "yes":
        return

    get_sample_position_register(int(measurement["sample_index"])).put(-1)


# =============================================================================
# Printing
# =============================================================================


def print_measurement_header(measurement, sample, sample_index, extra=None):
    legs = measurement["detectors"]

    print("")
    print("==============================================")
    print(f"Run name:       {measurement.get('run_name', '')}")
    print(f"Protocol:       {measurement.get('protocol_name', '')}")
    print(f"Run repeat:     {measurement.get('run_repeat', 1)}")
    print(f"Sample index:   {sample_index}")
    print(f"Sample name:    {sample.get('sample_name', '')}")
    print(f"Attenuation:    {measurement['att_level']}")
    print(f"num_repeats:    {measurement['num_repeats']}")
    print(f"sample_move:    {measurement['sample_move']}")
    print(f"position_reset: {measurement.get('position_reset', 'no')}")
    print(f"Detectors:      {len(legs)} in parallel")

    for leg in legs:
        owner = " [shutter owner]" if as_bool(leg.get("shutter_owner", False), "shutter_owner") else ""
        print(f"  - {leg_label(leg)}{owner}")
        print(f"      mode:          {leg['mode']}")
        print(f"      acq_time:      {leg['acq_time']}")
        print(f"      num_frames:    {leg['num_frames']}")
        print(f"      qmap_file:     {leg['qmap_file']}")
        print(f"      analysis_type: {leg.get('analysis_type', 'Multitau')}")
        print(f"      duration:      {leg_duration(leg):.1f} s per repeat")

        if leg.get("motors"):
            print(f"      motors:        {leg['motors']}")

    if extra:
        for line in extra:
            print(line)

    print("==============================================")
    print("")


# =============================================================================
# Huber setup
# =============================================================================
# The dual analogue of master_plan.py's DETECTOR_PLACEHOLDERS hooks, called from the same
# point in the sequence: once per measurement, before any acquisition starts.
#
# Deliberately does not reuse placeholder_rigaku3M(). That hook belongs to the serial path and
# should stay free to change for it -- if the two shared one function, editing it for a serial
# Rigaku run would silently move the dual geometry too.


def setup_huber_for_dual():
    """Move the huber to the dual-acquisition position: delta 10, nu 0.

    These two axes are the only huber motion in a dual run. Once this returns, nothing in the
    acquisition may touch them -- see FORBIDDEN_MOTORS in dual_acq.py, which refuses both a
    leg's motors block and the sample mesh.
    """
    huber = oregistry["huber"]

    print(f"Moving huber.delta to {DUAL_HUBER_DELTA}, huber.nu to {DUAL_HUBER_NU}")

    huber.delta.move(DUAL_HUBER_DELTA, wait=True)
    huber.nu.move(DUAL_HUBER_NU, wait=True)


# =============================================================================
# Run functions
# =============================================================================


def run_dual_measurement(measurement, sample_info):
    """Validate one expanded measurement, set up the sample, and run both detectors."""
    sample_index = int(measurement["sample_index"])
    sample = get_sample(sample_info, sample_index)

    validate_dual_measurement(measurement, sample)

    write_sample_registers(sample_index, sample)
    pv_registers.sample_move.put(measurement["sample_move"])
    reset_sample_position_register(measurement)

    att(int(measurement["att_level"]))

    print_measurement_header(measurement, sample, sample_index)

    setup_huber_for_dual()

    dual_acq_series(
        leg_specs=build_leg_specs(measurement),
        num_repeats=int(measurement["num_repeats"]),
        wait_time=float(measurement.get("wait_time", 0)),
        cam_timeout=measurement.get("cam_timeout"),
    )


def run_dual_measurement_info(
    measurement_info_file=DUAL_MEASUREMENT_INFO_FILE,
    sample_info_file=SAMPLE_INFO_FILE,
):
    """Expand dual_measurement_info.yaml and run every measurement it describes."""
    sample_info = read_yaml(sample_info_file)
    measurement_info = read_yaml(measurement_info_file)

    measurements = expand_measurements(measurement_info)

    print("")
    print(f"Loaded {len(measurements)} expanded dual-measurement blocks.")
    print("")

    for measurement in measurements:
        run_dual_measurement(measurement=measurement, sample_info=sample_info)


# =============================================================================
# Dry-run preview (no acquisitions executed)
# =============================================================================


def dry_run_dual_measurement_info(
    measurement_info_file=DUAL_MEASUREMENT_INFO_FILE,
    sample_info_file=SAMPLE_INFO_FILE,
    check_hardware=False,
):
    """Validate and preview without moving anything.

    Estimated time uses max() over the legs, not sum() -- that difference is the whole point
    of running them in parallel, so the preview shows the saving up front.

    check_hardware=False skips the device-connected and ophyd-path checks so this can be run
    off the beamline; pass True on a live session to validate those too.
    """
    sample_info = read_yaml(sample_info_file)
    measurement_info = read_yaml(measurement_info_file)

    measurements = expand_measurements(measurement_info)

    print("")
    print(f"Total dual measurements planned: {len(measurements)}")
    print("")

    total_parallel = 0.0
    total_serial = 0.0

    for measurement in measurements:
        sample_index = int(measurement["sample_index"])
        sample = get_sample(sample_info, sample_index)

        validate_dual_measurement(measurement, sample, check_hardware=check_hardware)

        legs = measurement["detectors"]
        num_repeats = int(measurement["num_repeats"])
        wait_time = float(measurement.get("wait_time", 0))

        durations = [leg_duration(leg) for leg in legs]

        parallel_time = (max(durations) + wait_time) * num_repeats
        serial_time = (sum(durations) + wait_time) * num_repeats

        total_parallel += parallel_time
        total_serial += serial_time

        extra = [
            f"Est. parallel:  {parallel_time:.1f} s",
            f"Est. if serial: {serial_time:.1f} s",
        ]

        print_measurement_header(measurement, sample, sample_index, extra=extra)

    saved = total_serial - total_parallel

    print(f"Total estimated time: {total_parallel:.1f} s ({total_parallel / 60:.1f} min)")
    print(f"Same measurements run serially: {total_serial:.1f} s ({total_serial / 60:.1f} min)")
    print(f"Saved by running in parallel: {saved:.1f} s ({saved / 60:.1f} min)")
    print("")


# =============================================================================
# Usage
# =============================================================================

# 1. Edit sample_info.yaml (shared with the serial path).
# 2. Edit dual_measurement_info.yaml.
# 3. In IPython/Bluesky:
#
#       from id8_common.plans.acquire.dual_master_plan import dry_run_dual_measurement_info
#       from id8_common.plans.acquire.dual_master_plan import run_dual_measurement_info
#
#       dry_run_dual_measurement_info(check_hardware=True)
#       run_dual_measurement_info()
#
# To point at a different file:
#
#       run_dual_measurement_info(
#           measurement_info_file="/home/beams10/8IDIUSER/bluesky/src/user_plans/my_dual.yaml",
#       )

__all__ = [
    "DUAL_HUBER_DELTA",
    "DUAL_HUBER_NU",
    "DUAL_MEASUREMENT_INFO_FILE",
    "build_leg_specs",
    "dry_run_dual_measurement_info",
    "run_dual_measurement",
    "run_dual_measurement_info",
    "setup_huber_for_dual",
]
