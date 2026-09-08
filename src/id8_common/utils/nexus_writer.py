"""NeXus metadata writer, built on Miaoqi Chu's nexus_xpcs_aps.

The single entry point the acquisition path calls. Replaces
utils/nexus_utils.py, which is retained as nexus_utils.txt for reference.

Division of labour:

    nexus_xpcs_aps.core.*           the schema factories        HIS
    utils/xpcs_schema_mc.py         8-ID's composition of them  ours (244 lines)
    utils/nexus_runtime.py          EPICS signal -> NeXus path  ours (~180 lines)
    nexus_xpcs_aps.core.utils       the HDF5 writer             HIS

His package has no runtime layer and cannot have one -- nothing upstream knows
that /entry/instrument/detector_1/distance comes from device_position.yaml. That
is why nexus_runtime.py stays ours.

Two upstream bugs are worked around here rather than in his tree, so that a
`git pull` of his repo cannot silently reintroduce them:

1. make_sample() returns leaves aliased to his module-level singletons -- 17 of
   17 shared between two calls, and shared with core.schema.xpcs_schema. His own
   documented `xpcs_schema.copy()` pattern therefore corrupts the template for
   the rest of the process. We deepcopy every call.
2. _compiled_plans is cached on id(schema). Handing it a fresh deepcopy each
   time means a freed address can be reused and return a stale plan -- observed
   colliding 30 times in 200 cycles, once writing a file with a leaf missing. We
   clear the cache every call, which costs nothing at one file per measurement.
"""

from copy import deepcopy
from typing import Any
from typing import Dict
from typing import Optional

from id8_common.utils.nexus_runtime import create_runtime_metadata_dict

#: Unit categories his keymap does not define. LOCAL PATCH, reported to Miaoqi
#: Chu 2026-09-08; delete each entry as it lands upstream.
#:
#: nexus_xpcs_aps.core.utils.default_units_keymap has 10 entries. A category it
#: does not know silently becomes the string "any" -- no warning, nothing in the
#: session output. Before the move to his writer our own keymap carried these
#: three, so without this patch the move would be a REGRESSION:
#:
#:   NX_VOLTAGE    4 keithley *V leaves and keysight_amp  ->  "any"
#:   NX_FREQUENCY  keysight_freq                          ->  "any"
#:   NX_PRESSURE   the Alicat pressure fields, when added ->  "any"
EXTRA_UNITS = {
    "NX_VOLTAGE": "V",
    "NX_FREQUENCY": "Hz",
    "NX_PRESSURE": "Pa",
}


def create_nexus_format_metadata(
    filename: str,
    det: Any,
    additional_metadata: Optional[Dict[str, Any]] = None,
):
    """Write one measurement's NeXus metadata file.

    Same signature and contract as the retired utils.nexus_utils version, so
    call sites did not have to change.

    Args:
        filename: full path of the .hdf to write
        det: the detector object the measurement used
        additional_metadata: {nexus_path: value} merged last, used by the dual
            path to give each leg its own geometry
    """
    from nexus_xpcs_aps.core import utils as mc_utils

    from id8_common.utils.xpcs_schema_mc import xpcs_schema as mc_schema

    # Teach his keymap the categories it is missing -- see EXTRA_UNITS. Applied
    # here rather than at import so it survives a reload of his module, and
    # setdefault so an upstream fix wins over our patch automatically.
    for category, unit in EXTRA_UNITS.items():
        mc_utils.default_units_keymap.setdefault(category, unit)

    # See docstring, bug 2.
    mc_utils._compiled_plans.clear()

    # See docstring, bug 1; also his update_schema_at_runtime mutates in place,
    # so without this one measurement's values persist into the next.
    runtime_schema = deepcopy(mc_schema)

    runtime_metadata = create_runtime_metadata_dict(det, additional_metadata)
    runtime_schema = mc_utils.update_schema_at_runtime(runtime_schema, runtime_metadata)

    mc_utils.create_nexus_format_metadata(filename, runtime_schema)
    return filename
