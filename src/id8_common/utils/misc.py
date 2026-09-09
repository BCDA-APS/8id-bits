"""Utility functions for 8ID beamline operations."""

import random

def ioc_alive(pv: str, timeout: float=0.5, retries: int=2) -> bool:
    """
    Check if an IOC is alive, return True if a PV responds, False otherwise.
    """
    import time
    from epics import caget

    for _ in range(retries):
        try:
            val = caget(pv, timeout=timeout)
            if val is not None:
                return True
        except Exception:
            pass
        time.sleep(0.1)
    return False


def stream_rois(det, stats_nums=(1, 2, 3), fields=("total",), hinted=("total",)):
    """Read out only the wanted stats-plugin fields for this detector.

    For each plugin number in stats_nums the plugin is marked "hinted" so the
    BestEffortCallback picks it up, and its read is narrowed to `fields`. Fields also
    listed in `hinted` get plotted; the rest are recorded but not plotted.
    """
    # Drop the dotted sub-attributes (e.g. "stats1.total") that the detector class ships
    # with -- the loop below re-adds each stats plugin by its plain name and then lets
    # the plugin's own read_attrs decide which sub-signals come along.
    det.read_attrs = [a for a in det.read_attrs if "." not in a]
    for n in stats_nums:
        stats_attr = f"stats{n}"

        if stats_attr not in det.read_attrs:
            det.read_attrs.append(stats_attr)

        stats = getattr(det, stats_attr)

        stats.kind = "hinted"

        stats.read_attrs = list(fields)

        for f in fields:
            sig = getattr(stats, f)
            sig.kind = "hinted" if f in hinted else "normal"

def get_machine_name() -> str:
    """Returns the name of a randomly selected analysis machine.

    Returns:
        str: one of the names in the `machines` list below -- currently
             'adamite' or 'amazonite', with equal probability.

    Both are accepted by the xpcs8-boost-corr workflow, whose choices are
    ['polaris', 'pioneer', 'adamite', 'amazonite']. califone was dropped from
    this list on 2026-09-09: it is NOT one of those choices, and being weighted
    two-thirds it was handing DM a machine it does not take on most calls.
    Anything added here must appear in the workflow's choices too.
    """
    # machines = ["adamite", "califone", "amazonite"]
    machines = ["adamite", "amazonite"]
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

