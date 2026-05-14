"""
Start Bluesky Data Acquisition sessions of all kinds.

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

from tiled.client import from_profile



# from apstools.devices import load_devices_from_yaml
# from id8_common.utils.misc import ioc_alive

# Configuration block
# Get the path to the instrument package
# Load configuration to be used by the instrument.
instrument_path = Path(__file__).parent
iconfig_path = instrument_path / "configs" / "iconfig.yml"
iconfig = load_config(iconfig_path)

# Additional logging configuration
# only needed if using different logging setup
# from the one in the apsbits package
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
bec, peaks = init_bec_peaks(iconfig)
cat = init_catalog(iconfig)
RE, sd = init_RE(iconfig, subscribers=[bec, cat])

# Optional SPEC callback block
# delete this block if not using SPEC
# if iconfig.get("SPEC_DATA_FILES", {}).get("ENABLE", False):
#     from .callbacks.demo_spec_callback import init_specwriter_with_RE
#     from .callbacks.demo_spec_callback import newSpecFile  # noqa: F401
#     from .callbacks.demo_spec_callback import spec_comment  # noqa: F401
#     from .callbacks.demo_spec_callback import specwriter  # noqa: F401

    # init_specwriter_with_RE(RE)

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
make_devices(clear=False, file="devices.yml", device_manager=instrument)
# make_devices(clear=False, file="devices_aps_only.yml", device_manager=instrument)
make_devices(clear=False, file="ad_devices.yml", device_manager=instrument)

from id8_common.devices.area_detector import ad_setup
ad_setup(oregistry["eiger4M"], iconfig)
ad_setup(oregistry["lambda2M"], iconfig)

# RIGAKU_TEST_PV = "8idRigaku3m:cam1:Manufacturer_RBV"
# EIGER_TEST_PV = "8idEiger4m:cam1:Manufacturer_RBV"

# if ioc_alive(RIGAKU_TEST_PV):
#     print("Rigaku3M IOC is up — loading detector")
#     load_devices_from_yaml("ad_devices.yml", oregistry)

# else:
#     print("Rigaku3M IOC not reachable — skipping detector load")


if host_on_aps_subnet():
    make_devices(clear=False, file="devices_aps_only.yml", device_manager=instrument)

pv_registers = oregistry["pv_registers"]

# Setup baseline stream with connect=False is default
# Devices with the label 'baseline' will be added to the baseline stream.
setup_baseline_stream(sd, oregistry, connect=False)

from .utils.check_file_dim import check_h5_shape
from id8_e.utils.peak import rock_and_move, center_x, center_y, center_delta
from id8_e.utils.check_file_dim import check_h5_shape

from .plans.sim_plan import sim_count_plan  # noqa: E402, F401
from .plans.sim_plan import sim_print_plan  # noqa: E402, F401
from .plans.sim_plan import sim_rel_scan_plan  # noqa: E402, F401

from id8_common.utils.misc import stream_rois
stream_rois(oregistry["eiger4M"])
stream_rois(oregistry["lambda2M"])

from id8_e.plans.lakeshore import *
from .plans.master_plan import run_measurement_info #, set_temp_lakeshore2
from .plans.sample_info_unpack import *
from .plans.scan_8ide import *
from .plans.nexus_acq_eiger_int import *
from .plans.nexus_acq_eiger_ext import *

from .plans.nexus_acq_lambda_int import *
# from .plans.nexus_acq_lambda_ext import *

from .plans.sample_info_unpack import *

from .plans.tetramm_acq import *
# from .plans.nexus_acq_rigaku_zdt import *

from hklpy2.user import *  
from .utils.hklpy2_setup import configure_hklpy2
configure_hklpy2(oregistry)

