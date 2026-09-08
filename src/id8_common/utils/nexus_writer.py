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

The local patches this module used to carry are gone: AZjk/nexus_xpcs_aps#1
merged on 2026-09-08 and upstream now fixes the make_sample() aliasing, the
id()-keyed plan cache, the invalid-JSON workflow_kwargs default, and the three
missing unit categories.

The deepcopy below is NOT one of those and must stay. It is not working around
a bug in his code -- his writer no longer pops keys out of the schema at all.
It protects OUR module-level `xpcs_schema_mc.xpcs_schema`: his
`update_schema_at_runtime()` assigns `node["data"] = value` in place, so handing
it the module-level dict would leave each measurement's values sitting in the
template, to be inherited by the next measurement that does not overwrite them.
"""

from copy import deepcopy
from typing import Any
from typing import Dict
from typing import Optional

from id8_common.utils.nexus_runtime import create_runtime_metadata_dict



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



    # See the module docstring: his update_schema_at_runtime() mutates in place,
    # so the module-level template must never be handed to it directly.
    runtime_schema = deepcopy(mc_schema)

    runtime_metadata = create_runtime_metadata_dict(det, additional_metadata)
    runtime_schema = mc_utils.update_schema_at_runtime(runtime_schema, runtime_metadata)

    mc_utils.create_nexus_format_metadata(filename, runtime_schema)
    return filename
