"""TEST-ONLY NeXus writer built on Miaoqi Chu's nexus_xpcs_aps (mc_refact, 95ab368).

Parallel to utils/nexus_utils.py -- that module is untouched and remains the
production writer. This one exists so a single acquisition can be run against
the upstream writer and the result compared, without disturbing anything.

How it differs from utils/nexus_utils.py:

* the schema comes from utils/xpcs_schema_mc.py, which composes our 125 leaves
  from ``nexus_xpcs_aps.core.*`` factories instead of a 977-line literal.
  Path-for-path identical to ours (verified: 0 missing, 0 extra).
* the entry writer is his ``core.utils.create_nexus_format_metadata``.
* the beamline layer is UNCHANGED: it reuses our own
  ``create_runtime_metadata_dict``, which reads ``expt`` and our oregistry.
  His ``deployment/id8_{e,i}`` packages are obsolete and need apsbits; they are
  never imported.

``nexus_xpcs_aps`` is deliberately NOT installed. It is reachable only when
PYTHONPATH includes the clone at
``~/Documents/Miaoqi/nexus_xpcs_aps_95ab368/src``. The imports below are
therefore done inside the function, so a normal session that never calls this
starts perfectly well without the package present.

To revert: restore the original call in plans/acquire/ad_acq.py and delete this
file plus utils/xpcs_schema_mc.py. Nothing else references either.
"""

import copy

from id8_common.utils.nexus_utils import create_runtime_metadata_dict


def create_nexus_format_metadata_mc(filename, det, additional_metadata=None):
    """Same signature and contract as utils.nexus_utils.create_nexus_format_metadata."""
    from nexus_xpcs_aps.core import utils as mc_utils

    from id8_common.utils.xpcs_schema_mc import xpcs_schema as mc_schema

    # His writer caches a compiled plan keyed by id(schema) (core/utils.py:36).
    # We hand it a fresh deepcopy every call, so a freed address could be reused
    # and return a stale plan. Clearing the cache costs nothing at our rate and
    # removes the hazard entirely.
    mc_utils._compiled_plans.clear()

    # deepcopy because his update_schema_at_runtime() mutates in place; without
    # it, one measurement's values would persist into the next.
    runtime_schema = copy.deepcopy(mc_schema)

    runtime_metadata = create_runtime_metadata_dict(det, additional_metadata)
    runtime_schema = mc_utils.update_schema_at_runtime(runtime_schema, runtime_metadata)

    mc_utils.create_nexus_format_metadata(filename, runtime_schema)
    return filename
