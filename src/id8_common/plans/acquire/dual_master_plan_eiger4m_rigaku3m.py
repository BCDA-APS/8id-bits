"""
YAML front end for parallel two-detector acquisition.

The dual counterpart to master_plan.py. Reads dual_measurement_info.yaml from expt.user_plan_dir, validates
it hard enough that nothing can move before an error surfaces, and hands per-detector "legs" to
dual_acq_eiger4m_rigaku3m.dual_acq_series().

Imported by startup.py, so these are already in the session (the __all__ at the bottom of
this file is the full list `import *` brings across):

    dry_run_dual_measurement_info(check_hardware=True)   # validate and preview, moves nothing
    run_dual_measurement_info()                          # go

Importing this alongside master_plan.py is safe: the two export no names in common, both
having an __all__. They do share one thing -- the `expt` run state -- and a dual run is
careful with it: it sets only the sample fields and sample_move globally, while the per-leg
fields (det_name, qmap_file, analysis_type, workflow_name) are swapped in one leg at a time
by dual_acq_eiger4m_rigaku3m.swapped_registers(), after the parallel window has closed, and restored
afterwards. So single-detector acquisition behaves exactly as it did before.

Run expansion (runs/protocols/samples/loop_order) is reused wholesale from master_plan.py, so
the two files behave identically there. What differs is the protocol body: a dual protocol
carries a `detectors:` list, and the acquisition settings that are per-detector live inside it
rather than at protocol level.

sample_info.yaml is shared with the serial path and read unchanged.
"""

from pathlib import Path

from id8_common.plans.acquire.dual_acq_eiger4m_rigaku3m import DUAL_LEGS
from id8_common.plans.acquire.dual_acq_eiger4m_rigaku3m import FORBIDDEN_MOTORS
from id8_common.plans.acquire.dual_acq_eiger4m_rigaku3m import dual_acq_series
from id8_common.plans.acquire.master_plan import expand_measurements
from id8_common.plans.acquire.master_plan import get_sample
from id8_common.plans.acquire.validators import as_bool
from id8_common.plans.acquire.validators import normalize_yes_no
from id8_common.plans.acquire.validators import read_yaml
from id8_common.plans.acquire.validators import require_fields
from id8_common.plans.acquire.validators import require_mode_devices
from id8_common.plans.acquire.validators import require_positive_int
from id8_common.plans.acquire.validators import reset_sample_position
from id8_common.plans.acquire.validators import validate_acq_time
from id8_common.plans.acquire import validators
from id8_common.plans.set.select_device import DETECTOR_ALIASES
from id8_common.plans.set.shutter_att import att
from id8_common.expt_config import expt
from id8_common.registry import get_ophyd_object


# The YAML plan files are resolved from configs/experiment.yml at CALL time, not
# import time -- see expt.user_plan_dir and the run functions at the bottom of
# this file. No plan path is hardcoded here. (This note describes the plan
# files; the constants immediately below are unrelated to it.)

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
            # Called for the exception, not the result -- float() raises on
            # anything that is not a number. A bad value has to be caught here,
            # while the run can still be fixed, rather than at metadata-writing
            # time with the data already on disk.
            float(value)


def validate_leg(leg):
    """Check one detector leg of a dual protocol.

    Resolves devices and motors, so this needs a live session -- it is the part
    the dry run skips unless asked for with check_hardware=True.
    """
    label = leg_label(leg)
    where = f"Leg '{label}'"

    require_fields(leg, REQUIRED_LEG_FIELDS, "leg", where=where)

    device = leg["device"]
    mode = leg["mode"]

    validators.validate_detector_mode(device, mode, where=where)

    if (device, mode) not in DUAL_LEGS:
        supported = ", ".join(f"{d}/{m}" for d, m in sorted(DUAL_LEGS))
        raise ValueError(
            f"Leg '{label}': {device}/{mode} is not supported for dual acquisition. Supported: {supported}"
        )

    validate_acq_time(leg["acq_time"], device, mode, where=where)
    require_positive_int(leg["num_frames"], "num_frames", where=where)

    if not str(leg["qmap_file"]).strip():
        raise ValueError(f"{where}: qmap_file must not be empty.")

    validators.validate_analysis_type(leg.get("analysis_type", "Multitau"), where=where)

    validate_geometry(leg.get("geometry"), label)

    # `or {}` covers both a leg with no motors block and one written as
    # `motors:` with nothing under it, which YAML reads as None.
    motors = leg.get("motors") or {}

    for dotted in motors:
        if dotted in FORBIDDEN_MOTORS:
            raise ValueError(
                f"Leg '{label}': '{dotted}' cannot appear in a dual protocol's motors block. "
                f"Both huber axes are positioned once before acquisition by "
                f"setup_huber_for_dual() (delta {DUAL_HUBER_DELTA}, nu {DUAL_HUBER_NU})."
            )
        get_ophyd_object(dotted)

    # Checks required_devices too -- softglue for an eiger4M External Series
    # leg used to go unchecked here and only fail once the run had started.
    require_mode_devices(device, mode, where=where)


def validate_shutter_owner(legs):
    """Exactly one leg of a dual protocol must set shutter_owner: yes.

    The owner is armed first and everyone else waits for it -- see the shutter
    contract at the top of dual_acq_eiger4m_rigaku3m.py. A single-leg protocol
    may leave it unset and owns the shutter by default; prepare_legs() marks it
    as the owner at run time.
    """
    owners = [leg for leg in legs if as_bool(leg.get("shutter_owner", False), "shutter_owner")]

    if len(legs) == 1 and not owners:
        return

    if len(owners) != 1:
        labels = ", ".join(leg_label(leg) for leg in owners)

        # No owner at all joins to the empty string, which reads as a truncated
        # message rather than as the real problem.
        if not labels:
            labels = "none"

        raise ValueError(
            f"Exactly one leg must set shutter_owner: yes (got {len(owners)}: {labels}). "
            f"The owner is armed first and everyone else waits for it to confirm it is acquiring."
        )


def validate_sample_motion(measurement, sample):
    # FORBIDDEN_MOTORS is the dual-only part: the mesh is the other way
    # huber.delta / huber.nu could be driven, via sample_info.yaml's
    # inner_motor / outer_motor. Refuse it for the same reason as a leg's
    # motors block -- setup_huber_for_dual() owns both axes.
    validators.validate_sample_motion(measurement, sample, forbidden_motors=FORBIDDEN_MOTORS)


def normalize_dual_measurement(measurement):
    """Coerce the protocol's yes/no fields to their string form, in place."""
    normalize_yes_no(measurement)


def validate_dual_measurement(measurement, sample, check_hardware=True):
    """Check one expanded dual measurement: protocol level first, then leg by leg.

    check_hardware=False skips validate_leg(), which is the part that resolves
    devices and so needs a live session. The sample-mesh check runs either way
    and resolves the mesh motors, so a sample_move: yes protocol still cannot be
    checked offline.
    """
    require_fields(measurement, REQUIRED_PROTOCOL_FIELDS, "protocol")
    require_fields(sample, ["sample_name", "header"], "sample")

    normalize_dual_measurement(measurement)

    legs = measurement["detectors"]

    if not isinstance(legs, list) or not legs:
        raise ValueError("protocol 'detectors' must be a non-empty list.")

    require_positive_int(measurement["num_repeats"], "num_repeats")

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
    reset_sample_position(measurement)


# =============================================================================
# Printing
# =============================================================================


def print_measurement_header(measurement, sample, sample_index, extra=None):
    """Print the per-measurement banner, shared by the real run and the dry run.

    `extra` is a list of already-formatted lines printed just before the closing
    rule -- how the dry run adds its parallel/serial time estimates.
    """
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
    """Huber positioning hook for a dual run. Motion is DISABLED: this moves nothing.

    As it stands the function only prints the delta 10 / nu 0 position it would have moved to
    -- see the comment below. Position the diffractometer yourself before the run.

    When the motion is re-enabled, these two axes are the only huber motion in a dual run, and
    once this returns nothing in the acquisition may touch them -- see FORBIDDEN_MOTORS in
    dual_acq_eiger4m_rigaku3m.py, which refuses both a leg's motors block and the sample mesh.
    """
    # DISABLED 2026-09-06 for testing: no motor motion. The dual geometry is
    # whatever the diffractometer is already at, so a leg's metadata may not
    # describe the true beam path -- see the geometry: block in
    # dual_measurement_info.yaml for how to override it per leg.
    # Re-enable by uncommenting the three lines below.
    #
    # huber = oregistry["huber"]
    # print(f"Moving huber.delta to {DUAL_HUBER_DELTA}, huber.nu to {DUAL_HUBER_NU}")
    # huber.delta.move(DUAL_HUBER_DELTA, wait=True)
    # huber.nu.move(DUAL_HUBER_NU, wait=True)
    print(
        f"setup_huber_for_dual: motion DISABLED -- leaving huber where it is "
        f"(would have moved delta to {DUAL_HUBER_DELTA}, nu to {DUAL_HUBER_NU})"
    )


# =============================================================================
# Run functions
# =============================================================================


def run_dual_measurement(measurement, sample_info):
    """Validate one expanded measurement, set up the sample, and run both detectors."""
    sample_index = int(measurement["sample_index"])
    sample = get_sample(sample_info, sample_index)

    validate_dual_measurement(measurement, sample)

    # Populate the SAMPLE half of the run state. Added 2026-09-06: without it
    # expt.header and expt.sample_name are never set, and the first call to
    # gen_folder_prefix() raises AttributeError before any hardware moves.
    # master_plan.run_measurement() has always done this; the dual path was
    # missed when the pv_registers -> expt migration landed.
    #
    # Only the sample half. The measurement half (det_name, mode, acq_time,
    # qmap_file, analysis_type) is PER LEG in a dual run and is applied one leg
    # at a time by swapped_registers() in dual_acq_eiger4m_rigaku3m -- setting it globally here
    # would stamp both legs with whichever leg happened to be written last.
    expt.sample_index = sample_index
    expt.set_measurement(sample=sample)

    expt.sample_move = measurement["sample_move"]
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
    measurement_info_file=None,
    sample_info_file=None,
):
    """Expand dual_measurement_info.yaml and run every measurement it describes."""
    # None, not a default argument: a default binds once at import and could not
    # follow an experiment.yml edit or expt.reload().
    measurement_info_file = measurement_info_file or expt.dual_measurement_info_file
    sample_info_file = sample_info_file or expt.sample_info_file

    print(f"Reading dual plans from {Path(measurement_info_file).parent}")

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
    measurement_info_file=None,
    sample_info_file=None,
    check_hardware=False,
):
    """Validate and preview without moving anything.

    Estimated time uses max() over the legs, not sum() -- that difference is the whole point
    of running them in parallel, so the preview shows the saving up front.

    check_hardware=False skips the per-leg device-connected checks; pass True on a live session
    to run them too. It does NOT skip the sample-mesh check, which resolves inner_motor and
    outer_motor either way, so a sample_move: yes protocol still needs a live session.
    """
    measurement_info_file = measurement_info_file or expt.dual_measurement_info_file
    sample_info_file = sample_info_file or expt.sample_info_file

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
#       from id8_common.plans.acquire.dual_master_plan_eiger4m_rigaku3m import dry_run_dual_measurement_info
#       from id8_common.plans.acquire.dual_master_plan_eiger4m_rigaku3m import run_dual_measurement_info
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
    "build_leg_specs",
    "dry_run_dual_measurement_info",
    "run_dual_measurement",
    "run_dual_measurement_info",
    "setup_huber_for_dual",
]
