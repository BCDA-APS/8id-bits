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

# Core Functions
from apsbits.core.best_effort_init import init_bec_peaks
from apsbits.core.catalog_init import init_catalog
from apsbits.core.instrument_init import init_instrument
from apsbits.core.instrument_init import make_devices
from apsbits.core.run_engine_init import init_RE

# Utility functions
from apsbits.utils.aps_functions import host_on_aps_subnet
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
offline_devices = []
offline_devices += safe_make_devices(file="devices.yml", device_manager=instrument)
offline_devices += safe_make_devices(file="ad_devices.yml", device_manager=instrument)
if offline_devices:
    print(f"\033[91m\n*** Devices not online: {offline_devices} ***\n\033[0m")
if host_on_aps_subnet(): # test this 
    make_devices(clear=False, file="devices_aps_only.yml", device_manager=instrument)
# make_devices(clear=False, file="devices.yml", device_manager=instrument)
# make_devices(clear=False, file="ad_devices.yml", device_manager=instrument)
# make_devices(clear=False, file="devices_aps_only.yml", device_manager=instrument)

from id8_common.devices.area_detector import ad_setup
ad_setup(oregistry["eiger4M"], iconfig)
ad_setup(oregistry["lambda2M"], iconfig)

pv_registers = oregistry["pv_registers"]

# Setup baseline stream with connect=False is default
# Devices with the label 'baseline' will be added to the baseline stream.
setup_baseline_stream(sd, oregistry, connect=False)

# Import useful tools
from .utils.check_file_dim import check_h5_shape
# from .utils.peak import rock_and_move, center_x, center_y, center_delta

# from .plans.sim_plan import sim_count_plan  # noqa: E402, F401
# from .plans.sim_plan import sim_print_plan  # noqa: E402, F401
# from .plans.sim_plan import sim_rel_scan_plan  # noqa: E402, F401

from .plans.shutter_logic import *

# hklpy2 setup - only for 8ide
from hklpy2.user import *  
from .utils.hklpy2_setup import configure_hklpy2
configure_hklpy2(oregistry)

from .utils.misc import stream_rois
stream_rois(oregistry["eiger4M"])
stream_rois(oregistry["lambda2M"])

# import acquire plans

from .plans.acquire.ad_acq import *

from .plans.nexus_acq_eiger_int import *
# from .plans.nexus_acq_eiger_ext import *
# from .plans.nexus_acq_lambda_int import *
from .plans.nexus_acq_lambda_ext import *
# from .plans.nexus_acq_rigaku_zdt import *
# from .plans.tetramm_acq import *
# from .plans.nexus_acq_eiger_int_wei import *

# import align plans
from .plans.align.scan_8id import *

# import calibrate plans

# 8ide plan import (legacy)
# from .plans.master_plan import run_measurement_info #, set_temp_lakeshore2
# from .plans.sample_info_unpack import *

# 8idi plan import (legacy)
from .plans.select_sample_env import select_sample_env
from .plans.select_diagnostics import *
from .plans.sample_info_unpack import select_sample
from .plans.select_detector import *
# from .plans.scan_8idi import *
from .plans.qnw_plans import *



