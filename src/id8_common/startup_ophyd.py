"""
Start an Ophyd-only session for 8-ID -- no RunEngine, and no Bluesky or
apstools beyond the two deliberate exceptions listed below.

This is the Ophyd-only counterpart to ``startup.py``. It stays alongside
that file rather than replacing it: either script can be imported for a
given session, and a colleague who wants to roll back simply imports
``startup.py`` again. See ``docs/starting-a-session.md`` for how to run it.

The script runs top to bottom in these sections, and the order matters --
each one needs the section above it to have happened already:

1. configuration         -- read ``configs/iconfig.yml``
2. EPICS signal defaults -- must be set before any EPICS signal object exists
3. devices               -- build everything in the three device YAML files
4. area detectors        -- plugin wiring that needs the detectors to exist
5. experiment settings   -- ``expt``, the static-config + run-state object
6. hklpy2                -- diffractometer support (see the exception below)
7. detector ROI streams
8. plans                 -- last, because several plan modules look their
                            devices up in ``oregistry`` at import time

Differences from ``startup.py``:

* Device loading uses ``id8_common.registry.safe_make_devices`` (hand-written
  -- no ``guarneri``/``ophyd-registry``) instead of
  ``apsbits.core.instrument_init``.
* There is no ``RunEngine``, no ``BestEffortCallback``, no databroker
  catalog subscription, and no queueserver branch -- none of them have a
  job to do without a ``RunEngine``. Plans call devices directly
  (``.put()``/``.get()``/``.move()``), synchronously.
* ``plans/align/scan_8id.py`` is not imported here -- it is a Bluesky
  generator plan and cannot run without a ``RunEngine``. Use
  ``plans/align/ophyd_scan.py`` instead.
* ``hklpy2`` (diffractometer UB-matrix support) is the one deliberate
  exception: it is built on Bluesky's ``Movable``/``Readable`` protocols,
  so importing it makes Bluesky importable in this process. No
  ``RunEngine`` is created and no Bluesky plan is ever run because of it --
  ``hklpy2``'s diffractometer math is plain synchronous code. Comment out
  the hklpy2 block below if you don't need diffractometer support.
  Nothing else anywhere in this Ophyd-only path needs ``apsbits`` or
  ``bluesky`` -- if you write a new plan/util module and it needs hklpy2,
  import it directly in that module (``from hklpy2.user import ...``)
  rather than relying on this file having already imported it; don't add
  a new top-level ``apsbits``/``bluesky`` import here for anything else.
* The two real area detectors (``eiger4M``, ``lambda2M``) are still built
  by ``apstools.devices.area_detector_factory.ad_creator`` from
  ``ad_devices.yml`` -- unlike hklpy2, this never imports Bluesky, and it
  is the same battle-tested detector-triggering/HDF5-writing code that
  ``id8_common/devices/area_detector.py`` already depends on directly (that
  file imports ``apstools.devices`` unconditionally, regardless of which
  startup script is used). Porting it out of apstools would mean
  re-implementing ~1000 lines of correctness-critical detector code for no
  behavior change; this is treated the same as that existing exception.
"""

import logging
import sys
from pathlib import Path

import yaml

from id8_common.registry import oregistry
from id8_common.registry import safe_make_devices

MAIN_NAMESPACE = "__main__"


# ---------------------------------------------------------- 1. configuration

# Read the same iconfig.yml that startup.py uses, but only the keys this
# script needs (OPHYD.*, TILED_*, DM_SETUP_FILE, ALLOW_AREA_DETECTOR_WARMUP).
# RUN_ENGINE/BEC/SPEC_DATA_FILES etc. are for startup.py only.
instrument_path = Path(__file__).parent
iconfig_path = instrument_path / "configs" / "iconfig.yml"
iconfig = yaml.safe_load(iconfig_path.read_text())

logger = logging.getLogger(__name__)
logger.info("Starting Ophyd-only session with iconfig: %s", iconfig_path)
print(f"[startup_ophyd] Starting Ophyd-only session ({iconfig_path})")

# -------------------------------------------------- 2. EPICS signal defaults

# Match apsbits.core.run_engine_init.init_RE()'s call to
# apsbits.utils.controls_setup.set_timeouts(): raise EpicsSignalBase's
# per-signal connection_timeout default from ophyd's own hardcoded 1.0s to
# this repo's configured OPHYD.TIMEOUTS.PV_CONNECTION (5s), so an
# individual signal doesn't give up before the overall `timeout=` deadline
# passed to safe_make_devices() has elapsed. Must run before any
# EpicsSignalBase is constructed anywhere in this process (ophyd raises if
# called later), so this has to happen here, before safe_make_devices().
from ophyd.signal import EpicsSignalBase

_timeouts = iconfig.get("OPHYD", {}).get("TIMEOUTS", {})
EpicsSignalBase.set_defaults(
    auto_monitor=True,
    timeout=_timeouts.get("PV_READ", 5),
    write_timeout=_timeouts.get("PV_WRITE", 5),
    connection_timeout=_timeouts.get("PV_CONNECTION", 5),
)


# ---------------------------------------------------------------- 3. devices


def _register(device):
    # Match apsbits.core.instrument_init.guarneri_namespace_loader (and
    # id8_common/utils/safe_devices.py's own _register for the Bluesky
    # path): a device is available as a bare name at the interactive
    # prompt (e.g. `eiger4M.connected`), not just via oregistry["eiger4M"].
    # Lost otherwise, since this bypasses apsbits' own loader entirely.
    oregistry.register(device)
    setattr(sys.modules[MAIN_NAMESPACE], device.name, device)


# Create the devices. Same devices.yml / ad_devices.yml / devices_aps_only.yml
# as startup.py, read by the hand-written loader instead of apsbits.
# safe_make_devices() is safe by default: a device is only skipped if it
# fails to build, or fails its basic wait_for_connection() check (genuinely
# offline) -- instead of taking down the rest of this script. It doesn't
# try to prove every lazily-declared PV on a device connects; a missing
# lazy PV only surfaces when a plan actually reads or writes it. See
# id8_common/registry.py:safe_make_devices.
print("[startup_ophyd] Connecting devices (devices.yml, ad_devices.yml, devices_aps_only.yml) ...")
offline_devices = []
offline_devices += safe_make_devices(file="devices.yml", register=_register)
offline_devices += safe_make_devices(file="ad_devices.yml", register=_register)
offline_devices += safe_make_devices(file="devices_aps_only.yml", register=_register)
_skip_note = f", {len(offline_devices)} skipped" if offline_devices else ""
print(f"[startup_ophyd] {len(oregistry)} device(s) connected{_skip_note}")
if offline_devices:
    print(f"\033[91m\n*** Devices not online: {offline_devices} ***\n\033[0m")

# --------------------------------------------------------- 4. area detectors

from id8_common.devices.area_detector import ad_setup

# Guard, not assume: safe_make_devices() above may have skipped any of
# these detectors as offline. Without this check, a skipped detector would
# still crash the whole session right here via oregistry["eiger4M"] --
# exactly the failure mode this file exists to avoid.
if "eiger4M" in oregistry:
    ad_setup(oregistry["eiger4M"], iconfig)
    print("[startup_ophyd] eiger4M area-detector plugins configured")
if "lambda2M" in oregistry:
    ad_setup(oregistry["lambda2M"], iconfig)
    print("[startup_ophyd] lambda2M area-detector plugins configured")
# rigaku3M: plugin config yes, warmup/priming no.
#
# AD_plugin_primed() compares cam.data_type with hdf1.data_type; on this
# detector they differ permanently (cam Int32, HDF1 UInt8 -- the ZDT
# sparsified output path), so it reports "not primed" on EVERY startup and
# AD_prime_plugin2() would fire a real exposure each time: image_mode ->
# Single, trigger_mode -> 0, acquire -> 1, 2 s wait, then restore. That
# would disturb a detector that is often mid-acquisition when a session
# starts, and it is unnecessary here -- hdf1 runs LazyOpen=Yes in Stream
# mode, which per apstools' own AD_plugin_primed docstring removes the need
# to prime at all. So hand ad_setup an iconfig with the warmup flag off:
# everything else (wait_for_plugins, blocking_callbacks, stage_sigs
# cleanup, hdf1.kind) still applies. Verified against live PVs 2026-09-03.
_iconfig_no_warmup = dict(iconfig, ALLOW_AREA_DETECTOR_WARMUP=False)
if "rigaku3M" in oregistry:
    ad_setup(oregistry["rigaku3M"], _iconfig_no_warmup)
    print("[startup_ophyd] rigaku3M area-detector plugins configured (warmup skipped)")

# ---------------------------------------------------- 5. experiment settings

# pv_registers is down to one live field: expt.measurement_num is 8ideSoft:Reg1
# (see PV_FIELDS in expt_config.py -- a counter that restarts overwrites data,
# so it must outlive the checkout). Everything else it used to carry moved to
# configs/experiment.yml and state/run_state.yml on 2026-09-06.

# Experiment settings (configs/experiment.yml) and the current measurement's
# run state. `expt` is the single source for everything the acquisition path
# reads: static settings from configs/experiment.yml, per-measurement values
# from measurement_info.yaml, persistent session state (sample index, mesh
# positions) in state/run_state.yml, and the measurement counter in
# 8ideSoft:Reg1. See id8_common/expt_config.py.
from id8_common.expt_config import expt  # noqa: E402

print(f"[expt_config] {expt}")

# ------------------------------------------------ 6. hklpy2 / diffractometer

# hklpy2 setup -- see module docstring for why this is the one place
# Bluesky becomes importable in this process.
from hklpy2.user import *  # noqa: F401, F403
from .utils.hklpy2_setup import configure_hklpy2

# Guard: configure_hklpy2() immediately calls set_diffractometer(psic) and
# psic.add_reflection(...) -- real use of the device, not just a lookup --
# so it cannot be made to degrade gracefully the way a plain oregistry.get()
# does elsewhere. If psic is offline, skip diffractometer setup entirely
# rather than crash the session over it.
if "psic" in oregistry:
    configure_hklpy2(oregistry)
    print("[startup_ophyd] hklpy2 diffractometer (psic) configured")
else:
    print("\033[91m*** psic not online: skipping hklpy2/diffractometer setup ***\033[0m")

# --------------------------------------------------- 7. detector ROI streams

from .utils.misc import stream_rois

if "eiger4M" in oregistry:
    stream_rois(oregistry["eiger4M"])
if "lambda2M" in oregistry:
    stream_rois(oregistry["lambda2M"])
# stats_nums=(1,): ad_creator only builds the plugins listed in
# ad_devices.yml, and rigaku3M declares stats1 only (eiger4M/lambda2M
# declare stats1-4). The default stats_nums=(1, 2, 3) would raise
# AttributeError on stats2 here and abort startup.
if "rigaku3M" in oregistry:
    stream_rois(oregistry["rigaku3M"], stats_nums=(1,))

# ------------------------------------------------------------------ 8. plans

print("[startup_ophyd] Importing plans ...")

# import acquire plans
from .plans.acquire.ad_acq import *  # noqa: F401, F403
from .plans.acquire.tetramm_acq import *  # noqa: F401, F403
from .plans.acquire.master_plan import *  # noqa: F401, F403
from .plans.acquire.dual_master_plan_eiger4m_rigaku3m import *  # noqa: F401, F403

# import align plans -- ophyd_scan, not scan_8id (see module docstring).
#
# The PLAIN names here, deliberately: this session never imports scan_8id, so
# there is nothing to collide with and `dscan(...)` is what people type. The
# suffixed aliases exist for the Bluesky session, which does star-import
# scan_8id -- see the alias block at the bottom of ophyd_scan.py.
from .plans.align.ophyd_scan import (  # noqa: F401
    a2scan,
    ascan,
    auto_att,
    d2scan,
    dmesh,
    dscan,
    dscan_ophyd,
    huber_x_lup,
    huber_y_lup,
    mesh,
    rheo_x_lup,
    rheo_y_lup,
    save_images,
    x_lup,
    y_lup,
)

# import set plans
from .plans.set.select_sample import select_sample  # noqa: F401
from .plans.set.select_device import *  # noqa: F401, F403
from .plans.set.qnw_plans import *  # noqa: F401, F403

print(f"[startup_ophyd] Ready -- {len(oregistry)} device(s) connected, plans imported.")
