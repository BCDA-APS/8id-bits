"""Beamline runtime layer for the NeXus metadata file.

This is the half that Miaoqi Chu's ``nexus_xpcs_aps`` does NOT provide, and
cannot: it maps 8-ID's live EPICS signals onto NeXus paths. His package supplies
the schema (``core.*`` factories) and the file writer
(``core.utils.create_nexus_format_metadata``); everything here is the wiring
between them and this instrument.

Ported verbatim from utils/nexus_utils.py on 2026-09-08 so the acquisition path
can run on his package alone. The two differences from the original:

* defaults come from ``nexus_xpcs_aps.core.utils.get_default_metadata()`` applied
  to the composed schema, instead of our hand-maintained default_metadata.py;
* nothing here writes a file -- see nexus_writer.py.

Device binds are ALWAYS ``oregistry.get()``, never ``oregistry["name"]``: this
module is imported by ad_acq, which is imported by startup, so a bracket lookup
raises KeyError at import and the session does not start. ``.get()`` binds None,
and the failure lands at write time inside det_acq_series()'s except-block --
one measurement's metadata instead of the whole beamline.
"""

import datetime
from typing import Any
from typing import Dict
from typing import Optional

from id8_common.expt_config import expt
from id8_common.registry import oregistry
from id8_common.devices.area_detector import frame_period
from id8_common.plans.set.select_device import DETECTOR_ALIASES
from id8_common.plans.set.select_device import _find_motor
from id8_common.plans.set.select_device import _load_config
from id8_common.plans.set.select_device import _resolve


# Devices this module may read when it writes metadata.
#
# ALWAYS oregistry.get(), NEVER oregistry["name"], even in a commented-out line
# that someone will one day uncomment. This module is imported by ad_acq, which
# is imported by startup, so a bracket lookup here raises KeyError at import and
# the SESSION DOES NOT START -- one dead IOC would cost you the whole beamline
# instead of one metadata field. `.get()` binds None instead; the failure then
# lands at write time as AttributeError on None, inside det_acq_series()'s
# except-block, so you lose that measurement's metadata and nothing else.
#
# Uncommenting a line below is not enough on its own: the device also has to be
# uncommented in configs/devices.yml and its IOC has to be up, or the value
# written is whatever None does. See docs/reference/nexus-writers.md.
# fofb_s09 = oregistry.get("fofb_s09")   # also needs a schema node -- none exists yet
filter_8ide = oregistry.get("filter_8ide")
lakeshore1 = oregistry.get("lakeshore1")
mono = oregistry.get("mono")
tetramm1 = oregistry.get("tetramm1")
# undulator_upstream = oregistry.get("undulator_upstream")
# undulator_downstream = oregistry.get("undulator_downstream")
huber = oregistry.get("huber")
sl4 = oregistry.get("sl4")
sl7 = oregistry.get("sl7")
wb_slit = oregistry.get("wb_slit")
mono_slit = oregistry.get("mono_slit")
# xbpm1 = oregistry.get("xbpm1")
# aps = oregistry.get("aps")
# keithley_chA = oregistry.get("keithley_chA")
# keithley_chB = oregistry.get("keithley_chB")
# bk_pid = oregistry.get("bk_pid")
# keysight = oregistry.get("keysight")
rheometer = oregistry.get("rheometer")
sample = oregistry.get("sample")
qnw_env1 = oregistry.get("qnw_env1")
qnw_env2 = oregistry.get("qnw_env2")
qnw_env3 = oregistry.get("qnw_env3")
pcd1 = oregistry.get("pcd1")   # Alicat PCD pressure controller, unit 1
pcd2 = oregistry.get("pcd2")   # Alicat PCD pressure controller, unit 2


#: Units already reported as unreadable, so a long run prints one note per unit
#: instead of two lines per measurement. Same convention as the scan CSV
#: template, which skips an unreadable column "with one note printed".
_env_warned = set()

#: Runtime paths already reported as unreadable or unknown, so a long run prints
#: one line per problem instead of one per measurement.
_metadata_warned = set()


def _resolve_runtime_values(deferred, schema):
    """Evaluate the deferred runtime values, dropping any that cannot be read.

    NOTE the name: this module imports ``_resolve`` from plans.set.select_device
    to turn a dotted device path into an object. A second ``_resolve`` defined
    here shadowed that import and broke every motor lookup -- caught on a trio
    acquisition, 2026-09-10. Nothing defined here may collide with an import.

    Each entry in ``deferred`` is a zero-argument callable rather than a value,
    which is the whole point: the dict used to be a literal, so the FIRST device
    that was missing or disconnected raised while the literal was still being
    built and took every other field with it. ``det_acq_series()`` then swallowed
    the exception, and the measurement wrote its data with no metadata file and
    no analysis job -- silently. That cost A0013-A0016 on 2026-09-09.

    Now a field that cannot be read simply contributes no key, so the value the
    schema declares stands instead, and the file keeps its full shape. One dead
    IOC costs its own fields and nothing else.

    Paths the schema does not declare are dropped here too. Passing one through
    would raise KeyError out of update_schema_at_runtime() and lose the whole
    file again -- the very failure this is meant to end -- so an unknown path is
    reported and skipped. That also makes the schema/runtime duplication safe:
    the two lists no longer have to agree for a measurement to survive, and where
    they disagree it is said out loud instead of discovered days later.

    Both kinds of problem are reported once per process, not once per repeat.
    """
    resolved = {}
    for path, getter in deferred.items():
        if not _schema_has(schema, path):
            if path not in _metadata_warned:
                _metadata_warned.add(path)
                print(f"[nexus] {path} is not in xpcs_schema -- not written. "
                      f"Add the leaf there, or drop it from nexus_runtime.")
            continue
        try:
            resolved[path] = getter()
        except Exception as exc:
            if path not in _metadata_warned:
                _metadata_warned.add(path)
                print(f"[nexus] {path} unreadable ({type(exc).__name__}: {exc}) "
                      f"-- left at the schema default.")
    return resolved


def _schema_has(schema, path):
    """True if `path` names a leaf the schema declares."""
    node = schema
    for part in path.strip("/").split("/"):
        if not isinstance(node, dict) or part not in node:
            return False
        node = node[part]
    return True


def _env_readings(device, label, fields):
    """Return the /entry/sample entries for one sample-environment unit, or {}.

    ``fields`` is ((leaf, attribute), ...) -- the leaf names under /entry/sample
    this unit owns, and the ophyd attribute each one reads from.

    An omitted path is NOT an omitted field. create_runtime_metadata_dict seeds
    itself from ``get_default_metadata(xpcs_schema)``, so a path left out here
    keeps the default declared in the schema and the dataset is still written.
    That is the point: the file keeps the same shape whether or not the unit's
    IOC is up, and the default is upstream's to define -- make_pressure and
    make_sample both set 0.0 -- rather than a number duplicated here.

    Two distinct failures are guarded. ``device is None`` means the unit was
    skipped at startup, because ``oregistry.get()`` binds None for an IOC that
    was down then. An exception from ``.get()`` means it connected at startup
    and has gone away since, or is timing out now.

    Before this guard either one raised out of the dict literal below and failed
    the WHOLE metadata write -- all 50 sample leaves and every other field, not
    just this unit's -- inside ``det_acq_series()``'s except-block, where it is
    printed and swallowed. One dead sample-environment IOC cost the entire file.

    Only the units in _SAMPLE_ENV_UNITS are guarded. huber, sample and rheometer
    are bound the same way and are deliberately left unguarded: they are
    positioners read through ``.position``, and a measurement whose sample stage
    is missing has bigger problems than its metadata.
    """
    if device is None:
        _warn_env_once(label, "not in the registry")
        return {}

    # Resolve the attributes BEFORE the try. Components are class attributes, so
    # getattr succeeds whatever the connection state -- an AttributeError here
    # can only mean a name in _SAMPLE_ENV_UNITS does not match the device class,
    # i.e. a bug in this file. Raise it. Falling back to the schema default for a
    # typo would write 0.0 over a perfectly healthy reading, for every file, with
    # nothing in the log to say so -- silently losing data is worse than the
    # crash this guard exists to prevent.
    try:
        signals = [(leaf, getattr(device, attr)) for leaf, attr in fields]
    except AttributeError as exc:
        raise AttributeError(f"_SAMPLE_ENV_UNITS names {label}.{exc}") from exc

    # Reading, by contrast, is a hardware condition: fall back quietly.
    try:
        return {f"/entry/sample/{leaf}": float(signal.get()) for leaf, signal in signals}
    except Exception as exc:
        _warn_env_once(label, f"unreadable ({exc})")
        return {}


def _warn_env_once(label, reason):
    """Report an unreadable sample-environment unit the first time only."""
    if label not in _env_warned:
        _env_warned.add(label)
        print(f"[nexus] {label} {reason} -- its fields left at the schema default")


#: Sample-environment units guarded by _env_readings, each independent of the
#: others: a dead 8idAlicat does not cost the QNW fields, or the reverse.
#: qnw_env1 is the Air QNW. qnwN_temperature_set is the demand, while
#: pcdN_pressure_set is the controller's own readback of the demand rather than
#: the value we last wrote, so it records what the unit is actually aiming at.
_SAMPLE_ENV_UNITS = (
    ("qnw_env1", (("qnw1_temperature", "readback"), ("qnw1_temperature_set", "setpoint"))),
    ("qnw_env2", (("qnw2_temperature", "readback"), ("qnw2_temperature_set", "setpoint"))),
    ("qnw_env3", (("qnw3_temperature", "readback"), ("qnw3_temperature_set", "setpoint"))),
    ("pcd1", (("pcd1_pressure", "pressure"), ("pcd1_pressure_set", "setpoint_rbv"))),
    ("pcd2", (("pcd2_pressure", "pressure"), ("pcd2_pressure_set", "setpoint_rbv"))),
)


def _get_ring_current():
    """Return the APS storage-ring current in mA.

    Raises IndexError if no device named "aps" is in the registry (i.e. the machine
    status IOC was not loaded at startup).
    """
    matching_devices = list(oregistry.findall(name="aps"))
    aps_device = matching_devices[0]

    return float(aps_device.current.get())


def _get_detector_config(det_name):
    """Look up det_name's config block from device_position.yaml.

    Resolves DETECTOR_ALIASES first (e.g. rigaku3M_epics -> rigaku3M, which share the
    same physical detector and config).
    """
    config = _load_config()
    key = DETECTOR_ALIASES.get(det_name, det_name)

    if key in config.get("detectors", {}):
        return config["detectors"][key]

    raise KeyError(f"No configuration found for detector '{det_name}' in device_position.yaml.")


def _get_motor_value(motors_cfg, name):
    """Return the current value of a named axis.

    An entry with no `device` is a fixed number, so its literal `position` is returned;
    otherwise the device named by `device` is resolved and its live `.position` is read.
    """
    for m in motors_cfg:
        if m["name"] == name:
            device_path = m.get("device")
            if device_path is None:
                return m.get("position", 0)
            return _resolve(device_path).position

    raise KeyError(f"No motor named '{name}' in device_position.yaml config.")


def create_runtime_metadata_dict(
    det: Optional[Any] = None,
    additional_metadata: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Create a dictionary with runtime metadata.

    A full list of possible metadata is given in the default_metadata dictionary.
    This function should be maintained by beamline staff to include all relevant
    metadata as needed for the experiment.

    Args:
        det: The detector object (default: None)
        additional_metadata: Additional metadata to include (default: None)

    Returns:
        The runtime metadata dictionary
    """
    # Defaults from HIS schema, via his own extractor, instead of our
    # hand-maintained default_metadata.py. get_default_metadata() walks the
    # composed schema and returns {path: declared default}, which is exactly
    # what default_metadata.py was maintained by hand to be.
    from nexus_xpcs_aps.core.utils import get_default_metadata

    from id8_common.utils.xpcs_schema import xpcs_schema

    runtime_metadata = dict(get_default_metadata(xpcs_schema))

    # Resolve the current detector's swing/translation axes and distance from device_position.yaml
    det_name = expt.det_name
    det_cfg = _get_detector_config(det_name)
    motors_cfg = det_cfg["motors"]

    # Per DETECTOR, not per experiment: eiger4M 75 um, rigaku3M 76 um,
    # lambda2M 55 um. Until 2026-09-06 this came from a single global in
    # experiment.yml, so every rigaku and lambda file recorded the eiger value.
    det_pixel_size = float(det_cfg["pixel_size"])
    horizontal = _get_motor_value(motors_cfg, "horizontal")
    vertical = _get_motor_value(motors_cfg, "vertical")
    swing_angle_horizontal = _get_motor_value(motors_cfg, "swing_angle_horizontal")
    swing_angle_vertical = _get_motor_value(motors_cfg, "swing_angle_vertical")
    sample_detector_distance = det_cfg["distance"]

    # Beam-center metadata is read from device_position.yaml's presets, not from
    # pv_registers, so it always matches the calibrated direct-beam values on file
    # (even if this measurement moved off that preset via a position override).
    # Detectors with no translation stage (e.g. lambda2M) have no `position` on their
    # horizontal/vertical entries at all -- default to 0 (there is no motion to report).
    beam_center_position_x = _find_motor(motors_cfg, "horizontal").get("position", 0)
    beam_center_position_y = _find_motor(motors_cfg, "vertical").get("position", 0)

    # Update the metadata with runtime values
    deferred = {
        # Entry level metadata
        "/entry/entry_identifier": lambda: "xpcs_20240214_120000",
        "/entry/entry_identifier_uuid": lambda: "550e8400-e29b-41d4-a716-446655440000",
        "/entry/scan_number": lambda: 1,
        "/entry/user/cycle": lambda: expt.cycle_name,
        "/entry/start_time": lambda: str(datetime.datetime.now()),
        "/entry/end_time": lambda: str(datetime.datetime.now()),  # fixme later
        "/entry/instrument/datamanagement/workflow_name": lambda: expt.workflow_name,

        # TODO: Change the detector direct beam position and detector position to real numbers

        # Read detector name and use that name to decide what fields to use to populate the rest
        "/entry/instrument/detector_1/detector_name": lambda: expt.det_name,

        # Define all degrees of freedom of the detector
        "/entry/instrument/detector_1/position_x": lambda: horizontal / 1000.0,
        "/entry/instrument/detector_1/position_y": lambda: vertical / 1000.0,

        "/entry/instrument/detector_1/beam_center_x": lambda: det_cfg["db_x"],
        "/entry/instrument/detector_1/beam_center_y": lambda: det_cfg["db_y"],
        "/entry/instrument/detector_1/beam_center_position_x": lambda: beam_center_position_x / 1000.0,
        "/entry/instrument/detector_1/beam_center_position_y": lambda: beam_center_position_y / 1000.0,

        # These below are shared by all detectors 
        "/entry/instrument/detector_1/count_time": lambda: det.cam.acquire_time.get(),
        # frame_period(), not cam.acquire_period: on the Rigaku that PV is the
        # gap AFTER each exposure, so the period is acquire_time + acquire_period
        # (see Rigaku3MCam). Reading the PV directly recorded the wrong number
        # twice over -- the intended period while the detector ran at half that
        # rate, and, once the setups started writing 0, a frame_time of 0.
        "/entry/instrument/detector_1/frame_time": lambda: frame_period(det.cam),
        "/entry/instrument/detector_1/qmap_file": lambda: expt.qmap_file,
        "/entry/instrument/detector_1/distance": lambda: sample_detector_distance,
        "/entry/instrument/detector_1/x_pixel_size": lambda: det_pixel_size,
        "/entry/instrument/detector_1/y_pixel_size": lambda: det_pixel_size,

        "/entry/instrument/detector_1/flightpath_swing": lambda: swing_angle_horizontal,
        "/entry/instrument/detector_1/flightpath_swing_vertical": lambda: swing_angle_vertical,

        "/entry/sample/lakeshore1": lambda: lakeshore1.readback_ch1.get(),
        # "/entry/sample/keithley_chA_SrcLevelV": keithley_chA.SrcLevelV_AO.value,
        # "/entry/sample/keithley_chA_SrcLevelI": keithley_chA.SrcLevelI_AO.value,
        # "/entry/sample/keithley_chB_SrcLevelV": keithley_chB.SrcLevelV_AO.value,
        # "/entry/sample/keithley_chB_SrcLevelI": keithley_chB.SrcLevelI_AO.value,
        # "/entry/sample/keithley_chA_MeasRangeV": keithley_chA.MeasRangeV_AO.value,
        # "/entry/sample/keithley_chA_MeasRangeI": keithley_chA.MeasRangeI_AO.value,
        # "/entry/sample/keithley_chB_MeasRangeV": keithley_chB.MeasRangeV_AO.value,
        # "/entry/sample/keithley_chB_MeasRangeI": keithley_chB.MeasRangeI_AO.value,
        # "/entry/sample/bk_pid_VAL": bk_pid.VAL.value,
        # "/entry/sample/bk_pid_RDBK": bk_pid.RDBK.value,
        # "/entry/sample/keysight_func": keysight.func_rbv.value,
        # "/entry/sample/keysight_freq": keysight.frequency_rbv.value,
        # "/entry/sample/keysight_amp": keysight.amplitude_rbv.value,
        # "/entry/sample/keysight_phase": keysight.phase_rbv.value,
        # "/entry/sample/keysight_pulse_width": keysight.pulse_width_rbv.value,
        # "/entry/sample/keysight_trigg_source": keysight.trigger_source_rbv.value,
        # "/entry/sample/keysight_trigg_edge": keysight.trigger_edge_rbv.value,
        # "/entry/sample/keysight_burst_count": keysight.burst_count_rbv.value,
        # "/entry/sample/keysight_burst_mode": keysight.burst_mode_rbv.value,
        # "/entry/sample/keysight_burst_state": keysight.burst_state_rbv.value,
        # "/entry/sample/keysight_output": keysight.output_rbv.value,
        
        "/entry/instrument/wb_slit/vertical_gap": lambda: wb_slit.vgap.position,
        "/entry/instrument/wb_slit/vertical_center": lambda: wb_slit.vcen.position,
        "/entry/instrument/wb_slit/horizontal_gap": lambda: wb_slit.hgap.position,
        "/entry/instrument/wb_slit/horizontal_center": lambda: wb_slit.hcen.position,
        "/entry/instrument/mono_slit/vertical_gap": lambda: mono_slit.vgap.position,
        "/entry/instrument/mono_slit/vertical_center": lambda: mono_slit.vcen.position,
        "/entry/instrument/mono_slit/horizontal_gap": lambda: mono_slit.hgap.position,
        "/entry/instrument/mono_slit/horizontal_center": lambda: mono_slit.hcen.position,
        "/entry/instrument/sl4/vertical_gap": lambda: sl4.v.size.position,
        "/entry/instrument/sl4/vertical_center": lambda: sl4.v.center.position,
        "/entry/instrument/sl4/horizontal_gap": lambda: sl4.h.size.position,
        "/entry/instrument/sl4/horizontal_center": lambda: sl4.h.center.position,
        "/entry/instrument/sl7/vertical_gap": lambda: sl7.v.size.position,
        "/entry/instrument/sl7/vertical_center": lambda: sl7.v.center.position,
        "/entry/instrument/sl7/horizontal_gap": lambda: sl7.h.size.position,
        "/entry/instrument/sl7/horizontal_center": lambda: sl7.h.center.position,
        "/entry/instrument/monochromator/energy": lambda: mono.energy.user_readback.value,
        "/entry/instrument/monochromator/wavelength": lambda: mono.wavelength.user_readback.value,
        "/entry/instrument/incident_beam/incident_energy": lambda: mono.energy.user_readback.value,
        "/entry/instrument/incident_beam/incident_energy_spread": lambda: 0.0001,

        # "/entry/instrument/incident_beam/fofb_s09_horizontal": fofb_s09.h_loop.get(),
        # "/entry/instrument/incident_beam/fofb_s09_vertical": fofb_s09.v_loop.get(),
        
        # "/entry/instrument/incident_beam/incident_beam_intensity": (
        #     xbpm1.current1.mean_value.get() / xbpm1.current_scales.ch1.get()
        #     + xbpm1.current2.mean_value.get() / xbpm1.current_scales.ch2.get()
        #     + xbpm1.current3.mean_value.get() / xbpm1.current_scales.ch3.get()
        #     + xbpm1.current4.mean_value.get() / xbpm1.current_scales.ch4.get()
        # ), 
        "/entry/instrument/incident_beam/ring_current": lambda: _get_ring_current(),
        # "/entry/instrument/undulator_1/gap": undulator_upstream.gap.position,
        # "/entry/instrument/undulator_1/energy": undulator_upstream.energy.position,
        # "/entry/instrument/undulator_1/taper": undulator_upstream.gap_taper.position,
        # "/entry/instrument/undulator_2/gap": undulator_downstream.gap.position,
        # "/entry/instrument/undulator_2/energy": undulator_downstream.energy.position,
        # "/entry/instrument/undulator_2/taper": undulator_downstream.gap_taper.position,
        "/entry/instrument/attenuator_1/attenuator_transmission": lambda: (filter_8ide.transmission.readback.get()),
        "/entry/instrument/attenuator_1/attenuator_index": lambda: (filter_8ide.index.readback.get()),
        "/entry/instrument/attenuator_2/attenuator_transmission": lambda: (0),
        "/entry/instrument/attenuator_2/attenuator_index": lambda: (0),
        
        "/entry/sample/position_x": lambda: sample.x.position,
        "/entry/sample/position_y": lambda: sample.y.position,
        "/entry/sample/position_z": lambda: sample.z.position,
        "/entry/sample/position_rheo_x": lambda: rheometer.x.position,
        "/entry/sample/position_rheo_y": lambda: rheometer.y.position,
        "/entry/sample/position_rheo_z": lambda: rheometer.z.position,

        "/entry/sample/huber_nu": lambda: huber.nu.position,
        "/entry/sample/huber_delta": lambda: huber.delta.position,
        "/entry/sample/huber_mu": lambda: huber.mu.position,
        "/entry/sample/huber_eta": lambda: huber.eta.position,
        "/entry/sample/huber_chi": lambda: huber.chi.position,
        "/entry/sample/huber_phi": lambda: huber.phi.position,
        "/entry/sample/huber_y": lambda: huber.y.position,
        "/entry/sample/huber_z": lambda: huber.z.position,
        "/entry/sample/huber_x": lambda: huber.x.position,
        "/entry/instrument/bluesky/parent_folder": lambda: (
            f"{expt.mount_point}/{expt.cycle_name}/"
            f"{expt.experiment_name}/data/"
        ),
    }
    # Sample environments -- three QNW cells and two Alicat PCD controllers.
    # Kept as its own table because those five share one shape and one guard;
    # everything else is guarded per field by _resolve_runtime_values() below.
    for _label, _fields in _SAMPLE_ENV_UNITS:
        deferred.update(
            {path: (lambda v=value: v)
             for path, value in _env_readings(globals().get(_label), _label, _fields).items()}
        )

    # Read every field, skipping the ones that cannot be read or that the schema
    # does not declare. See _resolve_runtime_values().
    runtime_updates = _resolve_runtime_values(deferred, xpcs_schema)

    runtime_metadata.update(runtime_updates)
    if additional_metadata is not None:
        runtime_metadata.update(additional_metadata)
    return runtime_metadata
