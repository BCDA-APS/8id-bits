"""Utility functions for 8ID-I beamline operations.

This module provides utility functions for common operations at the 8ID-I beamline,
including machine name selection and temperature string formatting.
"""

import random


def get_machine_name() -> str:
    """Returns the name of a randomly selected machine (A or B).

    Returns:
        str: Either 'adamite' or 'amazonite'
    """
    # machines = ["adamite", "califone", "amazonite"]
    machines = ["adamite", "califone"]
    return random.choice(machines)


def temp2str(temp: float) -> str:
    """Convert a temperature value to a formatted string.

    Args:
        temp (float): Temperature value to convert

    Returns:
        str: Formatted temperature string in the format 'XXXpY' where XXX is the
             zero-padded integer part and Y is the first decimal place
    """
    temp_str = str(temp)
    parts = temp_str.split(".")
    integer_part = parts[0]
    padded_integer_part = integer_part.zfill(3)
    decimal_part = parts[1]
    decimal_part_cropped = decimal_part[0]

    temp_name = f"{padded_integer_part}p{decimal_part_cropped}"
    return temp_name
