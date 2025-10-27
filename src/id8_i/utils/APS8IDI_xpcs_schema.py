"""
NeXus schema definition for APS 8-ID-I XPCS.

This module defines the customized NeXus schema used for storing XPCS data at the
APS 8-ID-I beamline. The schema follows the NXxpcs application definition and
includes metadata about the experiment, instrument, and data collection.
"""

# customized NeXus schema for APS 8-ID-I XPCS

xpcs_schema = {
    "entry": {
        "type": "NXentry",
        "required": True,
        "definition": {
            "type": "NX_CHAR",
            "required": True,
            "description": ("Official NeXus NXDL schema to which this file conforms"),
            "data": "NXxpcs",
        },
        "schema_version": {
            "type": "NX_CHAR",
            "required": True,
            "description": (
                "Version of the XPCS-Nexus schema to which this file conforms"
            ),
            "data": "0.1.0",
        },
        "entry_identifier": {
            "type": "NX_CHAR",
            "required": True,
            "description": (
                "Locally unique identifier for the experiment " "(a.k.a. run or scan)"
            ),
            "data": "entry_identifier",
        },
        "entry_identifier_uuid": {
            "type": "NX_CHAR",
            "required": False,
            "description": "UUID identifier for this entry",
            "data": "entry_identifier_uuid",
        },
        "beamline": {
            "type": "NX_CHAR",
            "required": True,
            "description": "Beamline identifier, e.g., 8-ID-I",
            "data": "APS-8-ID-I",
        },
        "scan_number": {
            "type": "NX_INT",
            "required": True,
            "deprecated": True,
            "description": (
                "DEPRECATED: Use the entry_identifier field. "
                "Scan number (must be an integer)"
            ),
            "data": 1,
        },
        "start_time": {
            "type": "NX_DATE_TIME",
            "required": True,
            "description": (
                "Starting time of experiment, such as " '"2021-02-11 11:22:33.445566Z"'
            ),
            "data": "start_time",
        },
        "end_time": {
            "type": "NX_DATE_TIME",
            "required": False,
            "description": (
                "Ending time of experiment, such as " '"2021-02-11 11:23:45Z"'
            ),
            "data": "start_time",
        },
        "instrument": {
            "type": "NXinstrument",
            "required": True,
            "description": "XPCS instrument Metadata",
            "detector_1": {
                "type": "NXdetector",
                "required": True,
                "description": (
                    "XPCS data is typically produced by area detector "
                    "(likely EPICS AreaDetector)"
                ),
                "beam_center_x": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": True,
                    "description": (
                        "Position of beam center, x axis, in detector's " "coordinates"
                    ),
                    "data": 1.0,
                },
                "beam_center_y": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": True,
                    "description": (
                        "Position of beam center, y axis, in detector's " "coordinates"
                    ),
                    "data": 1.0,
                },
                "beam_center_position_x": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": True,
                    "description": (
                        "Position of the detector, x axis, during data " "collection"
                    ),
                    "data": 1.0,
                },
                "beam_center_position_y": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": True,
                    "description": (
                        "Position of the detector, y axis, during data " "collection"
                    ),
                    "data": 1.0,
                },
                "position_x": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": True,
                    "description": (
                        "Position of the detector, x axis, during data " "collection"
                    ),
                    "data": 1.0,
                },
                "position_y": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": True,
                    "description": (
                        "Position of the detector, y axis, during data " "collection"
                    ),
                    "data": 1.0,
                },
                "rotation_x": {
                    "type": "NX_NUMBER",
                    "units": "NX_ANGLE",
                    "required": True,
                    "description": (
                        "Rotation of the detector, x axis, during data " "collection"
                    ),
                    "data": 0.0,
                },
                "rotation_y": {
                    "type": "NX_NUMBER",
                    "units": "NX_ANGLE",
                    "required": True,
                    "description": (
                        "Rotation of the detector, y axis, during data " "collection"
                    ),
                    "data": 0.0,
                },
                "rotation_z": {
                    "type": "NX_NUMBER",
                    "units": "NX_ANGLE",
                    "required": True,
                    "description": (
                        "Rotation of the detector, z axis, during data " "collection"
                    ),
                    "data": 0.0,
                },
                "count_time": {
                    "type": "NX_NUMBER",
                    "units": "NX_TIME",
                    "required": True,
                    "description": "Exposure time of frames, s",
                    "data": 1.0,
                },
                "detector_name": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "Detector name",
                    "data": "Eiger4m",
                },
                "qmap_file": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "Qmap name",
                    "data": "eiger4M_qmap_default.hdf",
                },
                "distance": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Distance between sample and detector",
                    "data": 12.0,
                },
                "flightpath_swing": {
                    "type": "NX_NUMBER",
                    "units": "NX_ANGLE",
                    "required": False,
                    "description": "Swing angle of the flight path",
                    "data": 0,
                },
                "frame_time": {
                    "type": "NX_NUMBER",
                    "units": "NX_TIME",
                    "required": True,
                    "description": (
                        "Exposure period (time between frame starts) of frames, s"
                    ),
                    "data": 1.0,
                },
                "x_pixel_size": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Length of pixel in x direction",
                    "data": 75e-6,
                },
                "y_pixel_size": {
                    "type": "NX_NUMBER",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Length of pixel in y direction",
                    "data": 75e-6,
                },
                "compression": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "Compression algorithm used for data",
                    "data": "bslz4",
                },
            },
            "incident_beam": {
                "type": "NXbeam",
                "description": "Incident beam Metadata",
                "required": True,
                "extent": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Size (2-D) of the beam at this position",
                    "data": 1.0e-6,
                },
                "incident_energy": {
                    "type": "NX_FLOAT",
                    "units": "NX_ENERGY",
                    "required": True,
                    "description": "Incident beam line energy (either keV or eV)",
                    "data": 12.0,
                },
                "incident_energy_spread": {
                    "type": "NX_FLOAT",
                    "units": "NX_ENERGY",
                    "required": False,
                    "description": (
                        "Spread of incident beam line energy (either keV or eV)"
                    ),
                    "data": 0.0001,
                },
                "incident_polarization_type": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": (
                        "Terse description of the incident beam polarization"
                    ),
                    "data": "linear_horizontal",
                },
                "incident_beam_intensity": {
                    "type": "NX_FLOAT",
                    "units": "NX_ENERGY",
                    "required": False,
                    "description": "Incident beam intensity, aka I0",
                    "data": 0.0001,
                },
                "transmitted_beam_intensity": {
                    "type": "NX_FLOAT",
                    "units": "NX_ENERGY",
                    "required": False,
                    "description": "Transmitted beam intensity, aka I1",
                    "data": 0.0001,
                },
                "ring_current": {
                    "type": "NX_FLOAT",
                    "units": "NX_CURRENT",
                    "required": True,
                    "description": "Storage ring current in mA",
                    "data": 0.0,
                },
            },
            "undulator_1": {
                "type": "NXinsertion_device",
                "required": False,
                "description": "Undulator 1 Metadata",
                "gap": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Undulator gap",
                    "data": 1.0,
                },
                "energy": {
                    "type": "NX_FLOAT",
                    "units": "NX_ENERGY",
                    "required": False,
                    "description": "Undulator energy",
                    "data": 1.0,
                },
                "taper": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Undulator taper",
                    "data": 1.0,
                },
            },
            "undulator_2": {
                "type": "NXinsertion_device",
                "required": False,
                "description": "Undulator 2 Metadata",
                "gap": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Undulator gap",
                    "data": 1.0,
                },
                "energy": {
                    "type": "NX_FLOAT",
                    "units": "NX_ENERGY",
                    "required": False,
                    "description": "Undulator energy",
                    "data": 1.0,
                },
                "taper": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Undulator taper",
                    "data": 1.0,
                },
            },
            "monochromator": {
                "type": "NXmonochromator",
                "required": False,
                "description": "Monochromator Metadata",
                "energy": {
                    "type": "NX_FLOAT",
                    "units": "NX_ENERGY",
                    "required": False,
                    "description": "Monochromator energy",
                    "data": 1.0,
                },
                "wavelength": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Monochromator wavelength",
                    "data": 1.0,
                },
            },
            "attenuator_1": {
                "type": "NXattenuator",
                "required": False,
                "description": "Attenuator 1 in 8IDE Metadata",
                "attenuator_transmission": {
                    "type": "NX_FLOAT",
                    "units": "NX_DIMENSIONLESS",
                    "required": False,
                    "description": "Attenuator transmission",
                    "data": 1.0,
                },
                "attenuator_index": {
                    "type": "NX_INT",
                    "units": "NX_DIMENSIONLESS",
                    "required": False,
                    "description": "Attenuator index",
                    "data": 1,
                },
            },
            "attenuator_2": {
                "type": "NXattenuator",
                "required": False,
                "description": "Attenuator 2 in 8IDI Metadata",
                "attenuator_transmission": {
                    "type": "NX_FLOAT",
                    "units": "NX_DIMENSIONLESS",
                    "required": False,
                    "description": "Attenuator transmission",
                    "data": 1.0,
                },
                "attenuator_index": {
                    "type": "NX_INT",
                    "units": "NX_DIMENSIONLESS",
                    "required": False,
                    "description": "Attenuator index",
                    "data": 1,
                },
            },
            "beam_stop": {
                "type": "NXbeam_stop",
                "required": False,
                "description": "Beam stop Metadata",
                "x_position": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Beam stop x position",
                    "data": 1.0,
                },
                "y_position": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Beam stop y position",
                    "data": 1.0,
                },
                "size": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Beam stop size",
                    "data": 1.0,
                },
            },
            "datamanagement": {
                "type": "NXnote",
                "required": False,
                "description": "Data management Metadata",
                "workflow_name": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "Data management workflow name",
                    "data": "boost_corr workflow",
                },
                "workflow_version": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "Data management workflow version",
                    "data": "0.0.1",
                },
                "workflow_kwargs": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "input arguments for the data management workflow",
                    "data": """{"qmap": "sample_name", "mask": "mask_name"}""",
                },
            },
            "bluesky": {
                "type": "NXnote",
                "required": False,
                "description": "Bluesky Metadata",
                "scan_id": {
                    "type": "NX_INT",
                    "required": False,
                    "description": "Scan ID",
                    "data": 1,
                },
                "bluesky_version": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "Bluesky version",
                    "data": "1.0.0",
                },
                "bluesky_plan": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "Bluesky Plan for this dataset",
                    "data": "mesh_scan",
                },
                "bluesky_plan_kwargs": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": (
                        "Configration for the Bluesky Plan for this dataset"
                    ),
                    "data": "mesh_scan_sample_kwargs",
                },
                "spec_file": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "spec file name for this dataset",
                    "data": "/path/to/this/spec_file_name",
                },
                "parent_folder": {
                    "type": "NX_CHAR",
                    "required": False,
                    "description": "parent folder of the data file",
                    "data": "/path/to/this/spec_file_name",
                },
            },
            "slits_4": {
                "type": "NXslit",
                "required": False,
                "description": "Slits 4",
                "horizontal_gap": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Horizontal size of the slits",
                    "data": 1.0,
                },
                "horizontal_center": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "Horizontal center of the slits",
                    "data": 0.0,
                },
                "vertical_gap": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "vertical size of the slits",
                    "data": 1.0,
                },
                "vertical_center": {
                    "type": "NX_FLOAT",
                    "units": "NX_LENGTH",
                    "required": False,
                    "description": "vertical center of the slits",
                    "data": 0.0,
                },
            },
        },
        "sample": {
            "type": "NXsample",
            "required": False,
            "description": "Sample environment and position metadata",
            "short_description": {
                "type": "NX_CHAR",
                "required": False,
                "description": "Short description of the sample",
                "data": "Short sample description",
            },
            "full_description": {
                "type": "NX_CHAR",
                "required": False,
                "description": "Full description of the sample",
                "data": "Full sample description",
            },
            "position_x": {
                "type": "NXpositioner",
                "required": False,
                "description": "Sample position, x",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "position_y": {
                "type": "NXpositioner",
                "required": False,
                "description": "Sample position, y",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "position_z": {
                "type": "NXpositioner",
                "required": False,
                "description": "Sample position, z",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "position_rheo_x": {
                "type": "NXpositioner",
                "required": False,
                "description": "Sample position, x",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "position_rheo_y": {
                "type": "NXpositioner",
                "required": False,
                "description": "Sample position, y",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "position_rheo_z": {
                "type": "NXpositioner",
                "required": False,
                "description": "Sample position, z",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "huber_x": {
                "type": "NXpositioner",
                "required": False,
                "description": "Huber sample position, x",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "huber_y": {
                "type": "NXpositioner",
                "required": False,
                "description": "Huber sample position, y",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "huber_z": {
                "type": "NXpositioner",
                "required": False,
                "description": "Huber sample position, z",
                "units": "NX_LENGTH",
                "data": 1.0,
            },
            "huber_nu": {
                "type": "NXpositioner",
                "units": "NX_ANGLE",
                "required": True,
                "description": "Huber nu angle",
                "data": 0.0,
            },
            "huber_delta": {
                "type": "NXpositioner",
                "units": "NX_ANGLE",
                "required": True,
                "description": "Huber delta angle",
                "data": 0.0,
            },
            "huber_mu": {
                "type": "NXpositioner",
                "units": "NX_ANGLE",
                "required": True,
                "description": "Huber mu angle",
                "data": 0.0,
            },
            "huber_eta": {
                "type": "NXpositioner",
                "units": "NX_ANGLE",
                "required": True,
                "description": "Huber eta angle",
                "data": 0.0,
            },
            "huber_chi": {
                "type": "NXpositioner",
                "units": "NX_ANGLE",
                "required": True,
                "description": "Huber chi angle",
                "data": 0.0,
            },
            "huber_phi": {
                "type": "NXpositioner",
                "units": "NX_ANGLE",
                "required": True,
                "description": "Huber phi angle",
                "data": 0.0,
            },
            "qnw_lakeshore": {
                "type": "NX_NUMBER",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "backup temperature reading from lakeshore for qnw",
                "data": 1.0,
            },
            "qnw1_temperature": {
                "type": "NX_NUMBER",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "Sample temperature actual in celsius",
                "data": 1.0,
            },
            "qnw1_temperature_set": {
                "type": "NX_NUMBER",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "Sample temperature setpoint in celsius",
                "data": 1.0,
            },
            "qnw2_temperature": {
                "type": "NX_NUMBER",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "Sample temperature actual in celsius",
                "data": 1.0,
            },
            "qnw2_temperature_set": {
                "type": "NX_NUMBER",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "Sample temperature setpoint in celsius",
                "data": 1.0,
            },
            "qnw3_temperature": {
                "type": "NX_NUMBER",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "Sample temperature actual in celsius",
                "data": 1.0,
            },
            "qnw3_temperature_set": {
                "type": "NX_NUMBER",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "Sample temperature setpoint in celsius",
                "data": 1.0,
            },
            "rheometer_shear_rate": {
                "type": "NX_FLOAT",
                "units": "NX_ANY",
                "required": False,
                "description": "The shear rate of the rheometer",
                "data": 1.0,
            },
            "rheometer_temperature": {
                "type": "NX_FLOAT",
                "units": "NX_TEMPERATURE",
                "required": False,
                "description": "The temperature of the rheometer",
                "data": 1.0,
            },
        },
        "user": {
            "type": "NXuser",
            "required": False,
            "description": "User Metadata",
            "name": {
                "type": "NX_CHAR",
                "required": True,
                "description": "User name",
                "data": "John Doe",
            },
            "email": {
                "type": "NX_CHAR",
                "required": True,
                "description": "User email",
                "data": "JohnDoe@mail.edu",
            },
            "institution": {
                "type": "NX_CHAR",
                "required": True,
                "description": "User institution",
                "data": "Inititution Name",
            },
            "cycle": {
                "type": "NX_CHAR",
                "required": True,
                "description": (
                    "Cycle during which the experiment was performed, e.g., 2025-2"
                ),
                "data": "cycle",
            },
            "proposal_id": {
                "type": "NX_CHAR",
                "required": True,
                "description": "proposal_id",
                "data": "none",
            },
        },
    },
}
