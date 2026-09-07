"""
Start Bluesky Data Acquisition session for 8-ID-E and 8-ID-I.

Includes:

* Python script
* IPython console
* Jupyter notebook
* Bluesky queueserver
"""

# Standard Library Imports
import logging
from pathlib import Path

# Not used by this module: `yaml` reached the interactive prompt only as a
# side effect of `from master_plan import *`, and master_plan stopped needing
# it when read_yaml moved to validators.py. startup_ophyd.py imports it the
# same way, explicitly.
import yaml  # noqa: F401

# Core Functions
from apsbits.core.best_effort_init import init_bec_peaks
from apsbits.core.catalog_init import init_catalog
from apsbits.core.instrument_init import init_instrument
from apsbits.core.run_engine_init import init_RE

# Utility functions
from apsbits.utils.baseline_setup import setup_baseline_stream

# Configuration functions
from apsbits.utils.config_loaders import load_config
from apsbits.utils.helper_functions import register_bluesky_magics
from apsbits.utils.helper_functions import running_in_queueserver
from apsbits.utils.logging_setup import configure_logging

# Core Functions
from tiled.client import from_profile

# from apstools.devices import load_devices_from_yaml
# from id8_common.utils.misc import ioc_alive
from id8_common.utils.safe_devices import safe_make_devices
from id8_common.utils.plot_mesh import plot_mesh

# Configuration block
# Get the path to the instrument package
# Load configuration to be used by the instrument.
instrument_path = Path(__file__).parent
iconfig_path = instrument_path / "configs" / "iconfig.yml"
iconfig = load_config(iconfig_path)

# Additional logging configuration, only needed if using different logging setup from the one in the apsbits package
extra_logging_configs_path = instrument_path / "configs" / "extra_logging.yml"
configure_logging(extra_logging_configs_path=extra_logging_configs_path)

logger = logging.getLogger(__name__)
logger.info("Starting Instrument with iconfig: %s", iconfig_path)

# initialize instrument
instrument, oregistry = init_instrument("guarneri")

# Discard oregistry items loaded above.
oregistry.clear()

# Configure the session with callbacks, devices, and plans.
# aps_dm_setup(iconfig.get("DM_SETUP_FILE"))

# Command-line tools, such as %wa, %ct, ...
register_bluesky_magics()

# Bluesky initialization block
if iconfig.get("TILED_PROFILE_NAME", {}):
    profile_name = iconfig.get("TILED_PROFILE_NAME")
    tiled_client = from_profile(profile_name)

bec, peaks = init_bec_peaks(iconfig)
cat = init_catalog(iconfig)
RE, sd = init_RE(iconfig, subscribers=[bec, cat])

# These imports must come after the above setup.
# Queue server block
if running_in_queueserver():
    ### To make all the standard plans available in QS, import by '*', otherwise import
    ### plan by plan.
    from apstools.plans import lineup2  # noqa: F401
    from bluesky.plans import *  # noqa: F403
else:
    # Import bluesky plans and stubs with prefixes set by common conventions.
    # The apstools plans and utils are imported by '*'.
    from apstools.plans import *  # noqa: F403
    from apstools.utils import *  # noqa: F403
    from bluesky import plan_stubs as bps  # noqa: F401
    from bluesky import plans as bp  # noqa: F401

# Experiment specific logic, device and plan loading. # Create the devices.
# safe_make_devices skips (with a warning) a device that fails to build, or
# fails its basic wait_for_connection() check (genuinely offline) -- instead
# of letting it take down the rest of startup. It doesn't try to prove every
# lazily-declared PV on a device connects; a missing lazy PV only surfaces
# when a plan actually reads or writes it. See
# id8_common/registry.py:safe_make_devices for how.
offline_devices = []
offline_devices += safe_make_devices(file="devices.yml", device_manager=instrument)
offline_devices += safe_make_devices(file="ad_devices.yml", device_manager=instrument)
offline_devices += safe_make_devices(file="devices_aps_only.yml", device_manager=instrument)
if offline_devices:
    print(f"\033[91m\n*** Devices not online: {offline_devices} ***\n\033[0m")

# Bridging step: copy the devices just loaded into id8_common's own registry,
# so plan modules written against `from id8_common.registry import oregistry`
# resolve devices the same way whether this Bluesky startup or
# startup_ophyd.py launched the session.
from id8_common.registry import oregistry as shared_oregistry

for _device in oregistry.root_devices:
    shared_oregistry.register(_device)

from id8_common.devices.area_detector import ad_setup

# Guard, not assume: safe_make_devices() above may have skipped any of
# these detectors as offline. Checked against shared_oregistry, not the bare
# guarneri `oregistry` -- that object defines __getitem__ but neither
# __contains__ nor __iter__, so `in` on it falls back to the legacy
# __getitem__(0), __getitem__(1), ... protocol and raises
# ComponentNotFound instead of doing a membership test.
if "eiger4M" in shared_oregistry:
    ad_setup(shared_oregistry["eiger4M"], iconfig)
if "lambda2M" in shared_oregistry:
    ad_setup(shared_oregistry["lambda2M"], iconfig)
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
if "rigaku3M" in shared_oregistry:
    ad_setup(shared_oregistry["rigaku3M"], _iconfig_no_warmup)

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

# Setup baseline stream with connect=False is default
# Devices with the label 'baseline' will be added to the baseline stream.
setup_baseline_stream(sd, oregistry, connect=False)

# Import useful tools
from .utils.check_file_dim import check_h5_shape
# from .utils.peak import rock_and_move, center_x, center_y, center_delta

# from .plans.sim_plan import sim_count_plan  # noqa: E402, F401
# from .plans.sim_plan import sim_print_plan  # noqa: E402, F401
# from .plans.sim_plan import sim_rel_scan_plan  # noqa: E402, F401

# from .plans.shutter_logic import *

# hklpy2 setup - only for 8ide
from hklpy2.user import *
from .utils.hklpy2_setup import configure_hklpy2

# Guard, checked against shared_oregistry for the same reason as the
# ad_setup guard above: configure_hklpy2() immediately calls
# set_diffractometer(psic) and psic.add_reflection(...) -- real use of the
# device -- so if psic is offline, skip diffractometer setup entirely
# rather than crash the session over it.
if "psic" in shared_oregistry:
    configure_hklpy2(oregistry)
else:
    print("\033[91m*** psic not online: skipping hklpy2/diffractometer setup ***\033[0m")

from .utils.misc import stream_rois
if "eiger4M" in shared_oregistry:
    stream_rois(shared_oregistry["eiger4M"])
if "lambda2M" in shared_oregistry:
    stream_rois(shared_oregistry["lambda2M"])
# stats_nums=(1,): ad_creator only builds the plugins listed in
# ad_devices.yml, and rigaku3M declares stats1 only (eiger4M/lambda2M
# declare stats1-4). The default stats_nums=(1, 2, 3) would raise
# AttributeError on stats2 here and abort startup.
if "rigaku3M" in shared_oregistry:
    stream_rois(shared_oregistry["rigaku3M"], stats_nums=(1,))

# import acquire plans

from .plans.acquire.ad_acq import *
from .plans.acquire.tetramm_acq import *
from .plans.acquire.master_plan import *

# Parallel two-detector acquisition (Eiger + Rigaku in one beam window). Separate from the
# serial path above and shares no state with it -- single-detector runs are unaffected.
from .plans.acquire.dual_master_plan_eiger4m_rigaku3m import *

# The prompt's `oregistry` is the id8_common registry, not the guarneri object
# bound at the top of this file: it supports `in`, len() and iteration, which
# the guarneri one does not (see the comment above the ad_setup guards). Until
# 2026-09-06 that happened only as a side effect of `from master_plan import *`
# re-exporting the name after guarneri had bound it. master_plan now has an
# __all__, so state it here rather than depend on import order.
oregistry = shared_oregistry

# import align plans
from .plans.align.scan_8id import *

# QZ added on 08/14:
# The Ophyd-only scans, under their suffixed names ONLY. This session also does
# `from .plans.align.scan_8id import *` above, and both modules define dscan,
# ascan, dmesh, mesh, d2scan and a2scan. Importing the plain names here would
# shadow Sam's generators with functions that run immediately -- and because
# Python evaluates arguments first, `RE(dscan(...))` would then perform the
# whole scan before raising on the non-generator. See the alias block at the
# bottom of ophyd_scan.py.
from .plans.align.ophyd_scan import (  # noqa: F401
    a2scan_ophyd,
    ascan_ophyd,
    auto_att_ophyd,
    d2scan_ophyd,
    dmesh_ophyd,
    dscan_ophyd,
    mesh_ophyd,
)

# import set plans
from .plans.set.select_sample import select_sample
from .plans.set.select_device import *
from .plans.set.qnw_plans import *

# import calibrate plans

# 8ide plan import (legacy)
# from .plans.master_plan import run_measurement_info #, set_temp_lakeshore2
# from .plans.sample_info_unpack import *

# 8idi plan import (legacy)
# from .plans.select_sample_env import select_sample_env
# from .plans.select_diagnostics import *
# from .plans.sample_info_unpack import select_sample
# from .plans.select_detector import *
# # from .plans.scan_8idi import *
# from .plans.qnw_plans import *



