"""The 8-ID NeXus schema: which fields exist in the metadata file.

Composed by *calling* the factories in ``nexus_xpcs_aps.core.*`` rather than
spelling out 977 lines of literal dicts. This file is OURS -- it decides what
8-ID records; the factories it calls are Miaoqi Chu's. It replaced a
hand-written dict of the same name, kept for reference as
``Archive/xpcs_schema.txt``, and reproduces it path for path.

Rules obeyed here
-----------------
* Only ``nexus_xpcs_aps.core.*`` is imported.  ``nexus_xpcs_aps.deployment.*`` is
  obsolete and needs apsbits -- it is never touched.
* OUR path names win.  Where an upstream factory names a field differently we
  rename the leaf (e.g. ``flightpath_swing_horizontal`` -> our
  ``flightpath_swing``) so
  that dual_acq's OVERRIDE_PATHS keep resolving.
* Every bare upstream module-level dict used more than once is deep-copied
  first, so no two instances are aliased to the same object.  In fact we deepcopy
  *everything* taken from upstream, because ``update_schema_at_runtime`` mutates the
  schema in place and aliasing would smear one detector's data onto another.
* ``xpcs_schema`` is exported with exactly the same shape as ours.

Provenance of every leaf is recorded in ``LEAF_PROVENANCE`` (``upstream:<factory>``)
so it is possible to report how much came from upstream vs. from us. Nothing
reads it today -- the one-off comparison script it served is gone -- but it is
cheap to maintain and answers "where did this leaf come from?" at the prompt.
"""

from copy import deepcopy

# --- upstream core factories ------------------------------------------------
# "upstream" throughout this file means the installed nexus_xpcs_aps package
# (Miaoqi Chu's -- see the README for the checkout path). Two conventions use it:
#
#   _upstream_<name>   an object imported from that package and used AS IS.
#                      The alias is not decoration: several of these names --
#                      user, bluesky, beam_stop, incident_beam -- are also names
#                      this module exports, so importing them bare would make
#                      `user = _tag(deepcopy(user), ...)` a self-reference.
#                      The prefix keeps "theirs" and "ours" apart in one file.
#
#   "upstream:<thing>" the LEAF_PROVENANCE label recording which upstream
#                      factory produced a leaf, e.g. "upstream:make_slits".
#
# Names WITHOUT the prefix (make_slits, make_detector, ...) are upstream too --
# they are factories we call rather than objects we reuse, so there is no
# ours/theirs pair to keep apart.
from nexus_xpcs_aps.core.schema import make_entry
from nexus_xpcs_aps.core.instrument.attenuator import make_attenuator
from nexus_xpcs_aps.core.instrument.beam_stop import beam_stop as _upstream_beam_stop
from nexus_xpcs_aps.core.instrument.bluesky import bluesky as _upstream_bluesky
from nexus_xpcs_aps.core.instrument.datamanagement import (
    datamanagement as _upstream_datamanagement,
)
from nexus_xpcs_aps.core.instrument.detector import make_detector
from nexus_xpcs_aps.core.instrument.diffractometer import make_diffractometer
from nexus_xpcs_aps.core.instrument.incident_beam import (
    incident_beam as _upstream_incident_beam,
)
from nexus_xpcs_aps.core.instrument.keysight_waveform_generator import (
    keysight_waveform_generator as _upstream_keysight,
)
from nexus_xpcs_aps.core.instrument.monochromator import (
    monochromator as _upstream_monochromator,
)
from nexus_xpcs_aps.core.instrument.slits import make_slits
from nexus_xpcs_aps.core.instrument.undulator import make_undulator
from nexus_xpcs_aps.core.sample import make_sample
from nexus_xpcs_aps.core.user import user as _upstream_user

BEAMLINE = "APS-8-ID-E,I"

_META = frozenset({"required", "deprecated", "type", "units", "description"})

# leaf NeXus path -> short provenance label
LEAF_PROVENANCE = {}


# --- small composition helpers ----------------------------------------------
def _leaves(node):
    """Yield (key, value) for direct children that are leaf datasets."""
    for k, v in node.items():
        if k in _META or not isinstance(v, dict):
            continue
        if "data" in v:
            yield k, v


def _tag(node, prefix, label):
    """Record ``label`` as the provenance of every leaf under ``node``."""
    for k, v in node.items():
        if k in _META or not isinstance(v, dict):
            continue
        path = f"{prefix}/{k}"
        if "data" in v:
            LEAF_PROVENANCE[path] = label
        else:
            _tag(v, path, label)
    return node


def _pick(src, names):
    """Deep-copied subset of a flat leaf dict, in the order given."""
    return {n: deepcopy(src[n]) for n in names}


def _drop(node, names):
    """Deep copy of ``node`` with the named leaves removed."""
    out = deepcopy(node)
    for n in names:
        out.pop(n)
    return out


def _rename(node, old, new, description=None):
    """Deep copy of ``node`` with leaf ``old`` re-keyed to ``new``.

    Key order is preserved (the renamed leaf stays where it was), so the HDF5
    write order is unchanged.
    """
    out = {}
    for k, v in node.items():
        if k == old:
            leaf = deepcopy(v)
            if description is not None:
                leaf["description"] = description
            out[new] = leaf
        else:
            out[k] = deepcopy(v)
    return out


def _flatten(group, prefix):
    """Flatten a nested upstream group into ``prefix``-named sibling leaves.

    Upstream's ``make_diffractometer`` and ``keysight_waveform_generator`` model their
    axes as children of an NXpositioner / NXnote group; our schema carries the
    same leaves flat under /entry/sample with a name prefix.  The leaf dicts
    themselves are upstream's, verbatim (deep-copied); only the key changes.
    """
    return {f"{prefix}{k}": deepcopy(v) for k, v in _leaves(group)}


# =============================================================================
# instrument
# =============================================================================

# detector_1 -- upstream make_detector, with the one field we name differently.
_detector_1 = _rename(
    make_detector(1, name="Eiger4m"),
    "flightpath_swing_horizontal",
    "flightpath_swing",
    description="Swing angle of the flight path",
)


instrument = {
    "type": "NXinstrument",
    "required": True,
    "description": "XPCS instrument Metadata",
    "detector_1": _tag(_detector_1, "/entry/instrument/detector_1", "upstream:make_detector"),
    "incident_beam": _tag(
        deepcopy(_upstream_incident_beam),
        "/entry/instrument/incident_beam",
        "upstream:incident_beam",
    ),
    "undulator_1": _tag(
        make_undulator(1), "/entry/instrument/undulator_1", "upstream:make_undulator"
    ),
    "undulator_2": _tag(
        make_undulator(2), "/entry/instrument/undulator_2", "upstream:make_undulator"
    ),
    "monochromator": _tag(
        deepcopy(_upstream_monochromator),
        "/entry/instrument/monochromator",
        "upstream:monochromator",
    ),
    "attenuator_1": _tag(
        make_attenuator(1), "/entry/instrument/attenuator_1", "upstream:make_attenuator"
    ),
    "attenuator_2": _tag(
        make_attenuator(2), "/entry/instrument/attenuator_2", "upstream:make_attenuator"
    ),
    "beam_stop": _tag(
        deepcopy(_upstream_beam_stop), "/entry/instrument/beam_stop", "upstream:beam_stop"
    ),
    "datamanagement": _tag(
        deepcopy(_upstream_datamanagement),
        "/entry/instrument/datamanagement",
        "upstream:datamanagement",
    ),
    "bluesky": _tag(
        deepcopy(_upstream_bluesky), "/entry/instrument/bluesky", "upstream:bluesky"
    ),
    # Four physically distinct slit assemblies, each an independent call so the
    # instances are never aliased.
    "sl4": _tag(
        make_slits(4, description="Slits 4"), "/entry/instrument/sl4", "upstream:make_slits"
    ),
    "sl7": _tag(
        make_slits(7, description="Slits 7"), "/entry/instrument/sl7", "upstream:make_slits"
    ),
    "wb_slit": _tag(
        make_slits("wb", description="white beam slit"),
        "/entry/instrument/wb_slit",
        "upstream:make_slits",
    ),
    "mono_slit": _tag(
        make_slits("mono", description="mono beam slit"),
        "/entry/instrument/mono_slit",
        "upstream:make_slits",
    ),
}

# =============================================================================
# sample
# =============================================================================

# Upstream make_sample gives us base + qnw + rheometer + huber_stage + lakeshore +
# keithley + bk_pid in one call.  We drop the two rheometer readings our schema
# does not carry, then bolt on the flattened Huber angles and Keysight fields.
_sample_full = make_sample(
    qnw=True,
    rheometer=True,
    huber_stage=True,
    lakeshore=True,
    keithley=True,
    bk_pid=True,
    # Alicat PCD pressure controllers. This is a COUNT, not a flag: unit i
    # contributes pcd{i}_pressure and pcd{i}_pressure_set. Two units are
    # installed at 8-ID (8idAlicat:PCD1: and PCD2:), so 2.
    #
    # Upstream added pressure support in AZjk/nexus_xpcs_aps#1 at our request
    # -- as a bool, one unnumbered pair -- then generalised it to N units in
    # 7e2b8b7. Adopting that renamed our leaves: pressure -> pcd1_pressure and
    # pressure_set -> pcd1_pressure_set. Metadata files written before
    # 2026-09-09 carry the old unnumbered paths.
    pressure=2,
)
_sample_core = _drop(_sample_full, ["rheometer_shear_rate", "rheometer_temperature"])

# Huber goniometer angles: upstream make_diffractometer, flattened to huber_<axis>.
_huber_angles = _flatten(make_diffractometer("Huber"), "huber_")

# Keysight waveform generator: upstream module dict, flattened to keysight_<field>.
_keysight_flat = _flatten(_upstream_keysight, "keysight_")

sample = {}
sample.update(_sample_core)
sample.update(_huber_angles)
sample.update(_keysight_flat)

_tag(_sample_core, "/entry/sample", "upstream:make_sample")
_tag(_huber_angles, "/entry/sample", "upstream:make_diffractometer")
_tag(_keysight_flat, "/entry/sample", "upstream:keysight_waveform_generator")

# =============================================================================
# user
# =============================================================================

user = _tag(deepcopy(_upstream_user), "/entry/user", "upstream:user")

# =============================================================================
# entry
# =============================================================================

xpcs_schema = {"entry": make_entry(BEAMLINE, instrument, sample, user)}

# entry-level scalars come from make_entry itself; tag only the ones it owns
for _k, _v in _leaves(xpcs_schema["entry"]):
    LEAF_PROVENANCE.setdefault(f"/entry/{_k}", "upstream:make_entry")

# Nothing in this schema is hand-written: every leaf below is sourced from a
# nexus_xpcs_aps.core factory or module-level dict.
HANDWRITTEN_LEAVES = {}
