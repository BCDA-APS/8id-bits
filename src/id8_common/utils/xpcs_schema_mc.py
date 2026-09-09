"""Candidate 8-ID NeXus schema composed from Miaoqi Chu's mc_refact core factories.

This module reproduces, path-for-path, the schema currently defined in
    /home/beams10/8IDIUSER/bluesky/src/id8_common/utils/xpcs_schema.py
but builds it by *calling* the factories in ``nexus_xpcs_aps.core.*`` instead of
spelling out 977 lines of literal dicts.

Rules obeyed here
-----------------
* Only ``nexus_xpcs_aps.core.*`` is imported.  ``nexus_xpcs_aps.deployment.*`` is
  obsolete and needs apsbits -- it is never touched.
* OUR path names win.  Where his factory names a field differently we rename the
  leaf (e.g. his ``flightpath_swing_horizontal`` -> our ``flightpath_swing``) so
  that dual_acq's OVERRIDE_PATHS keep resolving.
* Every bare module-level dict of his that is used more than once is deep-copied
  first, so no two instances are aliased to the same object.  In fact we deepcopy
  *everything* we take from him, because ``update_schema_at_runtime`` mutates the
  schema in place and aliasing would smear one detector's data onto another.
* ``xpcs_schema`` is exported with exactly the same shape as ours.

Provenance of every leaf is recorded in ``LEAF_PROVENANCE`` so the comparison
script can report how much came from his code vs. ours.
"""

from copy import deepcopy

# --- his core factories -----------------------------------------------------
from nexus_xpcs_aps.core.schema import make_entry
from nexus_xpcs_aps.core.instrument.attenuator import make_attenuator
from nexus_xpcs_aps.core.instrument.beam_stop import beam_stop as _mc_beam_stop
from nexus_xpcs_aps.core.instrument.bluesky import bluesky as _mc_bluesky
from nexus_xpcs_aps.core.instrument.datamanagement import (
    datamanagement as _mc_datamanagement,
)
from nexus_xpcs_aps.core.instrument.detector import make_detector
from nexus_xpcs_aps.core.instrument.diffractometer import make_diffractometer
from nexus_xpcs_aps.core.instrument.incident_beam import (
    incident_beam as _mc_incident_beam,
)
from nexus_xpcs_aps.core.instrument.keysight_waveform_generator import (
    keysight_waveform_generator as _mc_keysight,
)
from nexus_xpcs_aps.core.instrument.monochromator import (
    monochromator as _mc_monochromator,
)
from nexus_xpcs_aps.core.instrument.slits import make_slits
from nexus_xpcs_aps.core.instrument.undulator import make_undulator
from nexus_xpcs_aps.core.sample import make_sample
from nexus_xpcs_aps.core.user import user as _mc_user

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
    """Flatten a nested group of his into ``prefix``-named sibling leaves.

    His ``make_diffractometer`` and ``keysight_waveform_generator`` model their
    axes as children of an NXpositioner / NXnote group; our schema carries the
    same leaves flat under /entry/sample with a name prefix.  The leaf dicts
    themselves are his, verbatim (deep-copied); only the key changes.
    """
    return {f"{prefix}{k}": deepcopy(v) for k, v in _leaves(group)}


# =============================================================================
# instrument
# =============================================================================

# detector_1 -- his make_detector, with the one field we name differently.
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
    "detector_1": _tag(_detector_1, "/entry/instrument/detector_1", "his:make_detector"),
    "incident_beam": _tag(
        deepcopy(_mc_incident_beam),
        "/entry/instrument/incident_beam",
        "his:incident_beam",
    ),
    "undulator_1": _tag(
        make_undulator(1), "/entry/instrument/undulator_1", "his:make_undulator"
    ),
    "undulator_2": _tag(
        make_undulator(2), "/entry/instrument/undulator_2", "his:make_undulator"
    ),
    "monochromator": _tag(
        deepcopy(_mc_monochromator),
        "/entry/instrument/monochromator",
        "his:monochromator",
    ),
    "attenuator_1": _tag(
        make_attenuator(1), "/entry/instrument/attenuator_1", "his:make_attenuator"
    ),
    "attenuator_2": _tag(
        make_attenuator(2), "/entry/instrument/attenuator_2", "his:make_attenuator"
    ),
    "beam_stop": _tag(
        deepcopy(_mc_beam_stop), "/entry/instrument/beam_stop", "his:beam_stop"
    ),
    "datamanagement": _tag(
        deepcopy(_mc_datamanagement),
        "/entry/instrument/datamanagement",
        "his:datamanagement",
    ),
    "bluesky": _tag(
        deepcopy(_mc_bluesky), "/entry/instrument/bluesky", "his:bluesky"
    ),
    # Four physically distinct slit assemblies, each an independent call so the
    # instances are never aliased.
    "sl4": _tag(
        make_slits(4, description="Slits 4"), "/entry/instrument/sl4", "his:make_slits"
    ),
    "sl7": _tag(
        make_slits(7, description="Slits 7"), "/entry/instrument/sl7", "his:make_slits"
    ),
    "wb_slit": _tag(
        make_slits("wb", description="white beam slit"),
        "/entry/instrument/wb_slit",
        "his:make_slits",
    ),
    "mono_slit": _tag(
        make_slits("mono", description="mono beam slit"),
        "/entry/instrument/mono_slit",
        "his:make_slits",
    ),
}

# =============================================================================
# sample
# =============================================================================

# His make_sample gives us base + qnw + rheometer + huber_stage + lakeshore +
# keithley + bk_pid in one call.  We drop the two rheometer readings our schema
# does not carry, then bolt on the flattened Huber angles and Keysight fields.
_mc_sample = make_sample(
    qnw=True,
    rheometer=True,
    huber_stage=True,
    lakeshore=True,
    keithley=True,
    bk_pid=True,
    # Alicat PCD pressure controller -- pressure and pressure_set. Upstream
    # added this flag in AZjk/nexus_xpcs_aps#1 at our request; before that a
    # sample-environment device could not be added from outside the package.
    pressure=True,
)
_sample_core = _drop(_mc_sample, ["rheometer_shear_rate", "rheometer_temperature"])

# Huber goniometer angles: his make_diffractometer, flattened to huber_<axis>.
_huber_angles = _flatten(make_diffractometer("Huber"), "huber_")

# Keysight waveform generator: his module dict, flattened to keysight_<field>.
_keysight_flat = _flatten(_mc_keysight, "keysight_")

sample = {}
sample.update(_sample_core)
sample.update(_huber_angles)
sample.update(_keysight_flat)

_tag(_sample_core, "/entry/sample", "his:make_sample")
_tag(_huber_angles, "/entry/sample", "his:make_diffractometer")
_tag(_keysight_flat, "/entry/sample", "his:keysight_waveform_generator")

# =============================================================================
# user
# =============================================================================

user = _tag(deepcopy(_mc_user), "/entry/user", "his:user")

# =============================================================================
# entry
# =============================================================================

xpcs_schema = {"entry": make_entry(BEAMLINE, instrument, sample, user)}

# entry-level scalars come from make_entry itself; tag only the ones it owns
for _k, _v in _leaves(xpcs_schema["entry"]):
    LEAF_PROVENANCE.setdefault(f"/entry/{_k}", "his:make_entry")

# Nothing in this schema is hand-written: every leaf below is sourced from a
# nexus_xpcs_aps.core factory or module-level dict.
HANDWRITTEN_LEAVES = {}
