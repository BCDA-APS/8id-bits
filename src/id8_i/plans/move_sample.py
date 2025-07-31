"""
Sample movement plans for the 8ID-I beamline.

This module provides plans for moving samples in a grid pattern, supporting both
rheometer and regular sample stages.
"""

import numpy as np
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps

rheometer = oregistry["rheometer"]
sample = oregistry["sample"]
pv_registers = oregistry["pv_registers"]


def mesh_grid_move(
    sam_index: int,
    x_cen: float,
    x_radius: float,
    x_pts: int,
    y_cen: float,
    y_radius: float,
    y_pts: int,
):
    """Move sample in a grid pattern.

    Args:
        sam_index: Index of the sample (0 for rheometer, 1-27 for regular samples)
        x_cen: Center x position
        x_radius: Radius of x movement
        x_pts: Number of points in x direction
        y_cen: Center y position
        y_radius: Radius of y movement
        y_pts: Number of points in y direction
    """
    sample_pos_register = pv_registers.sample_position_register(sam_index)
    sam_pos = int(sample_pos_register.get())

    samx_list = np.linspace(x_cen - x_radius, x_cen + x_radius, num=x_pts)
    samy_list = np.linspace(y_cen - y_radius, y_cen + y_radius, num=y_pts)

    pos_index = np.mod(sam_pos + 1, x_pts * y_pts)
    x_pos = samx_list[np.mod(pos_index, x_pts)]
    y_pos = samy_list[int(np.floor(pos_index / x_pts))]

    if sam_index == 0:
        yield from bps.mv(rheometer.x, x_pos, rheometer.y, y_pos)
    elif sam_index >= 1 and sam_index <= 27:
        yield from bps.mv(sample.x, x_pos, sample.y, y_pos)
    elif sam_index >= 31 and sam_index <= 38:
        yield from bps.mv(rheometer.x, x_pos, rheometer.y, y_pos)
    pass

    yield from bps.mv(sample_pos_register, pos_index)
