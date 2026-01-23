"""
NeXus file utilities for the 8ID-I beamline.

This module provides utilities for creating and manipulating NeXus files,
including metadata handling and schema validation. It supports the XPCS
experiment format used at the APS 8-ID-I beamline.
"""

import datetime
from typing import Any
from typing import Dict
from typing import Optional
from typing import Union

import h5py
from apsbits.core.instrument_init import oregistry

from .default_metadata import default_metadata
from .xpcs_schema import xpcs_schema

def _get_ring_current():
   
    aps_list = list(oregistry.findall(name="aps"))
    machine = aps_list[0]

    return float(machine.current.get())
    


# detector = oregistry["detector"]
# sample = oregistry["sample"]
filter_8ide = oregistry["filter_8ide"]
lakeshore1 = oregistry["lakeshore1"]
flight_path_8idi = oregistry["flight_path_8idi"]
mono_8id = oregistry["mono_8id"]
# qnw_env1 = oregistry["qnw_env1"]
# qnw_env2 = oregistry["qnw_env2"]
# qnw_env3 = oregistry["qnw_env3"]
tetramm1 = oregistry["tetramm1"]
pv_registers = oregistry["pv_registers"]
# undulator_upstream = oregistry["undulator_upstream"]
# undulator_downstream = oregistry["undulator_downstream"]
huber = oregistry["huber"]
sl4 = oregistry["sl4"]
xbpm1 = oregistry["xbpm1"]
# aps = oregistry["aps"]

default_units_keymap = {
    "NX_COUNT": "one",  # Used for frame_sum, frame_average, delay_difference
    "NX_DIMENSIONLESS": "dim.less",  # Used for g2, g2_derr, dynamic_roi_map
    "NX_LENGTH": "m",  # Used for beam_center_x, beam_center_y, distance
    "NX_TIME": "s",  # Used for count_time, frame_time
    "NX_ENERGY": "keV",  # Used for incident_energy, incident_energy_spread
    "NX_PER_LENGTH": "1/Å",  # Used for dynamic_phi_list, dynamic_q_list
    "NX_TEMPERATURE": "K",  # Used for temperature, temperature_set
    "NX_CURRENT": "mA",  # Used for current, milliamper
    "NX_ANY": "any",  # Used for G2_unnormalized, two_time_corr_func
    "NX_ANGLE": "degree",  # Used for rotation_x, rotation_y, rotation_z
}


default_storage_dtype = {
    "NX_CHAR": "S1",
    "NX_NUMBER": "f8",
    "NX_INT": "i8",
}


def create_nexus_entry(
    group_or_fhdl: Union[h5py.Group, h5py.File],
    runtime_schema: Dict[str, Any],
    ignore: bool = False,
):
    """Create a NeXus entry from a dictionary that represents a Nexus format.

    Args:
        group_or_fhdl: The group or file handle to create the entry in
        runtime_schema: The dictionary that represents the Nexus format
        ignore: If True, ignore optional fields (default: False)
    """
    for key, val in runtime_schema.items():
        if key == "attributes":
            for attr_key, attr_val in val.items():
                group_or_fhdl.attrs[attr_key] = attr_val["data"]
        elif key not in ("type", "required", "deprecated", "units", "description"):
            required = val.pop("required", False)
            if ignore and not required:
                continue

            if "data" in val:
                data_item = val.pop("data")
                if data_item is not None:
                    handle = group_or_fhdl.create_dataset(key, data=data_item)
                else:
                    handle = group_or_fhdl.create_dataset(key, shape=(0,))
            else:
                handle = group_or_fhdl.create_group(key)

            dtype = val.pop("type", None)
            if dtype is not None:
                handle.attrs["NX_Class"] = dtype
            units = val.pop("units", None)
            if units is not None:
                handle.attrs["unit"] = default_units_keymap.get(units, "any")
            description = val.pop("description", None)
            if description is not None:
                handle.attrs["description"] = description

            create_nexus_entry(handle, val, ignore=ignore)


def update_schema_at_runtime(
    schema: Dict[str, Any],
    runtime_metadata: Dict[str, Any],
) -> Dict[str, Any]:
    """Update the metadata dictionary with runtime data.

    Args:
        schema: The metadata dictionary
        runtime_metadata: The runtime metadata dictionary

    Returns:
        The updated schema with runtime metadata
    """
    for path, value in runtime_metadata.items():
        components = path.lstrip("/").split("/")
        current = schema
        for comp in components:
            current = current[comp]
        current["data"] = value
    return schema


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
    # Create a copy of the default metadata dictionary
    runtime_metadata = default_metadata.copy()

    # Update the metadata with runtime values
    runtime_updates = {
        # Entry level metadata
        "/entry/entry_identifier": "xpcs_20240214_120000",
        "/entry/entry_identifier_uuid": "550e8400-e29b-41d4-a716-446655440000",
        "/entry/scan_number": 1,
        "/entry/user/cycle": pv_registers.cycle_name.get(),
        "/entry/start_time": str(datetime.datetime.now()),
        "/entry/end_time": str(datetime.datetime.now()),  # fixme later
        "/entry/instrument/datamanagement/workflow_name": pv_registers.workflow_name.get(),
        "/entry/instrument/detector_1/beam_center_x": pv_registers.current_db_x0.get(),
        "/entry/instrument/detector_1/beam_center_y": pv_registers.current_db_y0.get(),
        # All lengths changed to unit of meters
        "/entry/instrument/detector_1/beam_center_position_x": pv_registers.current_det_x0.get() / 1000.0,
        "/entry/instrument/detector_1/beam_center_position_y": pv_registers.current_det_y0.get() / 1000.0,
        # "/entry/instrument/detector_1/position_x": detector.x.position / 1000.0,
        # "/entry/instrument/detector_1/position_y": detector.y.position / 1000.0,
        "/entry/instrument/detector_1/count_time": det.cam.acquire_time.get(),
        "/entry/instrument/detector_1/frame_time": det.cam.acquire_period.get(),
        "/entry/instrument/detector_1/detector_name": pv_registers.det_name.get(),
        "/entry/instrument/detector_1/qmap_file": pv_registers.qmap_file.get(),
        "/entry/instrument/detector_1/distance": 2.02,
        # "/entry/instrument/detector_1/distance": flight_path_8idi.length.position / 1000.0,  # Not calibrated
        # "/entry/instrument/detector_1/flightpath_swing": flight_path_8idi.swing.position,
        "/entry/instrument/detector_1/flightpath_swing": -huber.delta.position,
        "/entry/instrument/detector_1/x_pixel_size": pv_registers.det_pixel_size.get(),
        "/entry/instrument/detector_1/y_pixel_size": pv_registers.det_pixel_size.get(),
        "/entry/instrument/beam_stop/x_position": flight_path_8idi.ds_x.position,
        "/entry/instrument/beam_stop/y_position": flight_path_8idi.ds_y.position,
        "/entry/instrument/slits_4/vertical_gap": sl4.v.size.position,
        "/entry/instrument/slits_4/vertical_center": sl4.v.center.position,
        "/entry/instrument/slits_4/horizontal_gap": sl4.h.size.position,
        "/entry/instrument/slits_4/horizontal_center": sl4.h.center.position,
        "/entry/instrument/monochromator/energy": mono_8id.energy_readback.get(),
        "/entry/instrument/monochromator/wavelength": mono_8id.energy_readback.get() / 12.4,
        "/entry/instrument/incident_beam/incident_energy": mono_8id.energy_readback.get(),
        "/entry/instrument/incident_beam/incident_energy_spread": 0.0001,
        "/entry/instrument/incident_beam/incident_beam_intensity": (
            xbpm1.current1.mean_value.get() / xbpm1.current_scales.ch1.get()
            + xbpm1.current2.mean_value.get() / xbpm1.current_scales.ch2.get()
            + xbpm1.current3.mean_value.get() / xbpm1.current_scales.ch3.get()
            + xbpm1.current4.mean_value.get() / xbpm1.current_scales.ch4.get()
        ),
        "/entry/instrument/incident_beam/ring_current": _get_ring_current(),        
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
        "/entry/sample/huber_nu": huber.nu.position,
        "/entry/sample/huber_delta": huber.delta.position,
        "/entry/sample/huber_mu": huber.mu.position,
        "/entry/sample/huber_eta": huber.eta.position,
        "/entry/sample/huber_chi": huber.chi.position,
        "/entry/sample/huber_phi": huber.phi.position,
        "/entry/sample/huber_y": huber.sample_y.position,
        "/entry/sample/huber_z": huber.sample_z.position,
        "/entry/sample/huber_x": huber.sample_x.position,
        "/entry/sample/qnw_lakeshore": lakeshore1.readback_ch3.get(),
        # "/entry/sample/qnw1_temperature": qnw_env1.readback.get(),  # Air QNW
        # "/entry/sample/qnw1_temperature_set": qnw_env1.setpoint.get(),
        # "/entry/sample/qnw2_temperature": qnw_env2.readback.get(),
        # "/entry/sample/qnw2_temperature_set": qnw_env2.setpoint.get(),
        # "/entry/sample/qnw3_temperature": qnw_env3.readback.get(),
        # "/entry/sample/qnw3_temperature_set": qnw_env3.setpoint.get(),
        "/entry/instrument/bluesky/parent_folder": (
            f"{pv_registers.mount_point.get()}/{pv_registers.cycle_name.get()}/"
            f"{pv_registers.experiment_name.get()}/data/"
        ),
        "/entry/instrument/bluesky/spec_file": pv_registers.spec_file.get(),
    }
    # update the runtime metadata with the runtime updates
    runtime_metadata.update(runtime_updates)
    if additional_metadata is not None:
        runtime_metadata.update(additional_metadata)
    return runtime_metadata


def create_nexus_format_metadata(
    filename: str,
    det: Any,
    additional_metadata: Optional[Dict[str, Any]] = None,
):
    """Create a nexus file with the metadata for the xpcs experiment.

    Args:
        filename: Name of the nexus file to create
        det: Detector object to get the metadata from
        additional_metadata: Additional metadata to add (optional)
    """
    # create a copy of schema from the template, tree-structure of the nexus file
    runtime_schema = xpcs_schema.copy()

    # create a dictionary of the runtime metadata
    runtime_metadata = create_runtime_metadata_dict(det, additional_metadata)
    # update the schema with the metadata
    runtime_schema = update_schema_at_runtime(runtime_schema, runtime_metadata)

    # save the schema to a nexus file
    with h5py.File(filename, "w") as f:
        create_nexus_entry(f, runtime_schema)
    return


# if __name__ == "__main__":
#     # create_nexus_template(filename="template_metadata.hdf")
#     create_nexus_format_metadata("test02.hdf")
