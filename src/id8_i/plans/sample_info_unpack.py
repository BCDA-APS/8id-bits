"""
Module for selecting samples and reading sample information from a YAML
configuration file. Supports both rheometer and regular sample stages.
"""

import yaml
from pathlib import Path
from typing import Dict
from typing import Union

import numpy as np
from apsbits.core.instrument_init import oregistry
# from bluesky import plan_stubs as bps

sample = oregistry["sample"]
rheometer = oregistry["rheometer"]
filter = oregistry["filter_8ide"]
pv_registers = oregistry["pv_registers"]
huber = oregistry["huber"]

SAMPLE_INFO_PATH = Path("/home/beams/8IDIUSER/bluesky/src/user_plans/sample_info.yaml")


def select_sample(env: int):
    """Select and move to a sample position.

    This plan reads sample position information from a YAML file and moves
    either the rheometer (env=0) or sample stage (env=1-27) to the specified
    position.

    Args:
        env: Sample environment index (0 for rheometer, 1-27 for samples)

    Yields:
        Generator: Bluesky plan messages
    """
    with open(SAMPLE_INFO_PATH, "r") as f:
        loaded_dict = yaml.safe_load(f)

    sample_key = f"sample_{env}"
    x_cen = loaded_dict[sample_key]["x_cen"]
    y_cen = loaded_dict[sample_key]["y_cen"]
    # sample_name = loaded_dict[sample_key]["sample_name"] # Variable unused in original, kept commented or remove

    print(f"Moving {sample_key} x to {x_cen} and y to {y_cen}")

    if env == 0:
        rheometer.x.move(x_cen)
        rheometer.y.move(y_cen)
    elif 1 <= env <= 27:
        sample.x.move(x_cen)
        sample.y.move(y_cen)
    elif env == 31:
        huber.sample_x.move(x_cen)
        huber.sample_y.move(y_cen)
    else:
        pass

    pv_registers.qnw_index.put(env)


def sort_qnw() -> Dict[str, Union[int, float, str]]:
    """Read and organize sample information from the configuration file.

    This function reads the sample information YAML file and returns a dictionary
    containing the current sample's metadata, including position, dimensions,
    and measurement parameters.

    Returns:
        Dictionary containing sample metadata with the following keys:
        - qnw_index: Sample environment index
        - meas_num: Measurement number
        - sample_name: Name of the sample
        - header: Measurement header prefix
        - x_cen, y_cen: Center position coordinates
        - x_radius, y_radius: Scan range in each direction
        - x_pts, y_pts: Number of points in each direction
        - temp_zone: Temperature zone information
    """
    qnw_index = int(pv_registers.qnw_index.get())
    sample_key = f"sample_{qnw_index}"
    with open(SAMPLE_INFO_PATH, "r") as f:
        loaded_dict = yaml.safe_load(f)

    sam_dict = {
        "qnw_index": int(pv_registers.qnw_index.get()),
        "meas_num": int(pv_registers.measurement_num.get()),
        "sample_name": loaded_dict[sample_key]["sample_name"],
        "header": loaded_dict[sample_key]["header"],
        "x_cen": float(loaded_dict[sample_key]["x_cen"]),
        "y_cen": float(loaded_dict[sample_key]["y_cen"]),
        "x_radius": float(loaded_dict[sample_key]["x_radius"]),
        "y_radius": float(loaded_dict[sample_key]["y_radius"]),
        "x_pts": int(loaded_dict[sample_key]["x_pts"]),
        "y_pts": int(loaded_dict[sample_key]["y_pts"]),
        "temp_zone": loaded_dict[sample_key]["temp_zone"],
    }

    return sam_dict


def gen_folder_prefix() -> str:
    """Generate a folder name prefix for the current measurement.

    This function creates a standardized folder name using the sample metadata
    and current attenuation level.

    Returns:
        Formatted folder name string
    """
    sam_dict = sort_qnw()
    pv_registers.measurement_num.put(int(sam_dict["meas_num"]) + 1)

    pv_registers.sample_name.put(sam_dict["sample_name"])
    sample_name = sam_dict["sample_name"]
 

    att_level = int(filter.attenuation.readback.get())

    header_name = f'{sam_dict["header"]}{sam_dict["meas_num"]:04d}'
    folder_name = f"{header_name}_{sample_name}_a{att_level:04}"
    # print(folder_name)

    return folder_name


def mesh_grid_move():
    """Move to the next position in a mesh grid scan.

    This plan calculates and moves to the next position in a mesh grid pattern,
    using either the rheometer or sample stage depending on the current
    environment.

    Yields:
        Generator: Bluesky plan messages
    """
    sam_dict = sort_qnw()

    sample_pos_register = pv_registers.sample_position_register(sam_dict["qnw_index"])
    sam_pos = int(sample_pos_register.get())

    samx_list = np.linspace(
        sam_dict["x_cen"] - sam_dict["x_radius"],
        sam_dict["x_cen"] + sam_dict["x_radius"],
        num=sam_dict["x_pts"],
    )
    samy_list = np.linspace(
        sam_dict["y_cen"] - sam_dict["y_radius"],
        sam_dict["y_cen"] + sam_dict["y_radius"],
        num=sam_dict["y_pts"],
    )

    pos_index = np.mod(sam_pos + 1, sam_dict["x_pts"] * sam_dict["y_pts"])
    x_pos = samx_list[np.mod(pos_index, sam_dict["x_pts"])]
    y_pos = samy_list[int(np.floor(pos_index / sam_dict["x_pts"]))]

    if sam_dict["qnw_index"] == 0:
        rheometer.x.move(x_pos)
        rheometer.y.move(y_pos)
    elif sam_dict["qnw_index"] >= 1 and sam_dict["qnw_index"] <= 27:
        sample.x.move(x_pos)
        sample.y.move(y_pos)
    elif sam_dict["qnw_index"] == 31:
        huber.sample_x.move(x_pos)
        huber.sample_y.move(y_pos)
    else:
        pass

    sample_pos_register.put(pos_index)