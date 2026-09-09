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
    runtime_updates = {
        # Entry level metadata
        "/entry/entry_identifier": "xpcs_20240214_120000",
        "/entry/entry_identifier_uuid": "550e8400-e29b-41d4-a716-446655440000",
        "/entry/scan_number": 1,
        "/entry/user/cycle": expt.cycle_name,
        "/entry/start_time": str(datetime.datetime.now()),
        "/entry/end_time": str(datetime.datetime.now()),  # fixme later
        "/entry/instrument/datamanagement/workflow_name": expt.workflow_name,

        # TODO: Change the detector direct beam position and detector position to real numbers

        # Read detector name and use that name to decide what fields to use to populate the rest
        "/entry/instrument/detector_1/detector_name": expt.det_name,

        # Define all degrees of freedom of the detector
        "/entry/instrument/detector_1/position_x": horizontal / 1000.0,
        "/entry/instrument/detector_1/position_y": vertical / 1000.0,

        "/entry/instrument/detector_1/beam_center_x": det_cfg["db_x"],
        "/entry/instrument/detector_1/beam_center_y": det_cfg["db_y"],
        "/entry/instrument/detector_1/beam_center_position_x": beam_center_position_x / 1000.0,
        "/entry/instrument/detector_1/beam_center_position_y": beam_center_position_y / 1000.0,

        # These below are shared by all detectors 
        "/entry/instrument/detector_1/count_time": det.cam.acquire_time.get(),
        "/entry/instrument/detector_1/frame_time": det.cam.acquire_period.get(),
        "/entry/instrument/detector_1/qmap_file": expt.qmap_file,
        "/entry/instrument/detector_1/distance": sample_detector_distance,
        "/entry/instrument/detector_1/x_pixel_size": det_pixel_size,
        "/entry/instrument/detector_1/y_pixel_size": det_pixel_size,

        "/entry/instrument/detector_1/flightpath_swing": swing_angle_horizontal,
        "/entry/instrument/detector_1/flightpath_swing_vertical": swing_angle_vertical,

        "/entry/sample/lakeshore1": lakeshore1.readback_ch1.get(),
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
        
        "/entry/instrument/wb_slit/vertical_gap": wb_slit.vgap.position,
        "/entry/instrument/wb_slit/vertical_center": wb_slit.vcen.position,
        "/entry/instrument/wb_slit/horizontal_gap": wb_slit.hgap.position,
        "/entry/instrument/wb_slit/horizontal_center": wb_slit.hcen.position,
        "/entry/instrument/mono_slit/vertical_gap": mono_slit.vgap.position,
        "/entry/instrument/mono_slit/vertical_center": mono_slit.vcen.position,
        "/entry/instrument/mono_slit/horizontal_gap": mono_slit.hgap.position,
        "/entry/instrument/mono_slit/horizontal_center": mono_slit.hcen.position,
        "/entry/instrument/sl4/vertical_gap": sl4.v.size.position,
        "/entry/instrument/sl4/vertical_center": sl4.v.center.position,
        "/entry/instrument/sl4/horizontal_gap": sl4.h.size.position,
        "/entry/instrument/sl4/horizontal_center": sl4.h.center.position,
        "/entry/instrument/sl7/vertical_gap": sl7.v.size.position,
        "/entry/instrument/sl7/vertical_center": sl7.v.center.position,
        "/entry/instrument/sl7/horizontal_gap": sl7.h.size.position,
        "/entry/instrument/sl7/horizontal_center": sl7.h.center.position,
        "/entry/instrument/monochromator/energy": mono.energy.user_readback.value,
        "/entry/instrument/monochromator/wavelength": mono.wavelength.user_readback.value,
        "/entry/instrument/incident_beam/incident_energy": mono.energy.user_readback.value,
        "/entry/instrument/incident_beam/incident_energy_spread": 0.0001,

        # "/entry/instrument/incident_beam/fofb_s09_horizontal": fofb_s09.h_loop.get(),
        # "/entry/instrument/incident_beam/fofb_s09_vertical": fofb_s09.v_loop.get(),
        
        # "/entry/instrument/incident_beam/incident_beam_intensity": (
        #     xbpm1.current1.mean_value.get() / xbpm1.current_scales.ch1.get()
        #     + xbpm1.current2.mean_value.get() / xbpm1.current_scales.ch2.get()
        #     + xbpm1.current3.mean_value.get() / xbpm1.current_scales.ch3.get()
        #     + xbpm1.current4.mean_value.get() / xbpm1.current_scales.ch4.get()
        # ), 
        "/entry/instrument/incident_beam/ring_current": _get_ring_current(),
        # "/entry/instrument/undulator_1/gap": undulator_upstream.gap.position,
        # "/entry/instrument/undulator_1/energy": undulator_upstream.energy.position,
        # "/entry/instrument/undulator_1/taper": undulator_upstream.gap_taper.position,
        # "/entry/instrument/undulator_2/gap": undulator_downstream.gap.position,
        # "/entry/instrument/undulator_2/energy": undulator_downstream.energy.position,
        # "/entry/instrument/undulator_2/taper": undulator_downstream.gap_taper.position,
        "/entry/instrument/attenuator_1/attenuator_transmission": (filter_8ide.transmission.readback.get()),
        "/entry/instrument/attenuator_1/attenuator_index": (filter_8ide.index.readback.get()),
        "/entry/instrument/attenuator_2/attenuator_transmission": (0),
        "/entry/instrument/attenuator_2/attenuator_index": (0),
        
        "/entry/sample/position_x": sample.x.position,
        "/entry/sample/position_y": sample.y.position,
        "/entry/sample/position_z": sample.z.position,
        "/entry/sample/position_rheo_x": rheometer.x.position,
        "/entry/sample/position_rheo_y": rheometer.y.position,
        "/entry/sample/position_rheo_z": rheometer.z.position,
        "/entry/sample/qnw1_temperature": qnw_env1.readback.get(),  # Air QNW
        "/entry/sample/qnw1_temperature_set": qnw_env1.setpoint.get(),
        "/entry/sample/qnw2_temperature": qnw_env2.readback.get(),
        "/entry/sample/qnw2_temperature_set": qnw_env2.setpoint.get(),
        "/entry/sample/qnw3_temperature": qnw_env3.readback.get(),
        "/entry/sample/qnw3_temperature_set": qnw_env3.setpoint.get(),

        # Alicat PCD, two units. pcdN_pressure_set is the controller's own
        # readback of the demand, not the value we last wrote, so it records
        # what the unit is actually aiming at.
        "/entry/sample/pcd1_pressure": pcd1.pressure.get(),
        "/entry/sample/pcd1_pressure_set": pcd1.setpoint_rbv.get(),
        "/entry/sample/pcd2_pressure": pcd2.pressure.get(),
        "/entry/sample/pcd2_pressure_set": pcd2.setpoint_rbv.get(),

        "/entry/sample/huber_nu": huber.nu.position,
        "/entry/sample/huber_delta": huber.delta.position,
        "/entry/sample/huber_mu": huber.mu.position,
        "/entry/sample/huber_eta": huber.eta.position,
        "/entry/sample/huber_chi": huber.chi.position,
        "/entry/sample/huber_phi": huber.phi.position,
        "/entry/sample/huber_y": huber.y.position,
        "/entry/sample/huber_z": huber.z.position,
        "/entry/sample/huber_x": huber.x.position,
        "/entry/instrument/bluesky/parent_folder": (
            f"{expt.mount_point}/{expt.cycle_name}/"
            f"{expt.experiment_name}/data/"
        ),
    }
    # update the runtime metadata with the runtime updates
    runtime_metadata.update(runtime_updates)
    if additional_metadata is not None:
        runtime_metadata.update(additional_metadata)
    return runtime_metadata
