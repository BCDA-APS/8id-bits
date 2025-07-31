"""
Start Bluesky Data Acquisition sessions of all kinds.

Includes:

* Python script
* IPython console
* Jupyter notebook
* Bluesky queueserver
"""

import logging
from pathlib import Path

from apsbits.core.best_effort_init import init_bec_peaks
from apsbits.core.catalog_init import init_catalog
from apsbits.core.instrument_init import make_devices
from apsbits.core.instrument_init import oregistry
from apsbits.core.run_engine_init import init_RE
from apsbits.utils.aps_functions import aps_dm_setup
from apsbits.utils.aps_functions import host_on_aps_subnet
from apsbits.utils.config_loaders import get_config
from apsbits.utils.config_loaders import load_config
from apsbits.utils.helper_functions import register_bluesky_magics
from apsbits.utils.helper_functions import running_in_queueserver
from id8_i.plans.ad_setup_plans import ad_initial_setup

logger = logging.getLogger(__name__)
logger.bsdev(__file__)

# Get the path to the instrument package
instrument_path = Path(__file__).parent

# Load configuration to be used by the instrument.
iconfig_path = instrument_path / "configs" / "iconfig.yml"
load_config(iconfig_path)

# Get the configuration
iconfig = get_config()

logger.info("Starting Instrument with iconfig: %s", iconfig_path)

# Discard oregistry items loaded above.
oregistry.clear()

# Configure the session with callbacks, devices, and plans.
# aps_dm_setup(iconfig.get("DM_SETUP_FILE")) #TODO: This line is broken

# Command-line tools, such as %wa, %ct, ...
register_bluesky_magics()

# Initialize core bluesky components
bec, peaks = init_bec_peaks(iconfig)
cat = init_catalog(iconfig)
RE, sd = init_RE(iconfig, bec_instance=bec, cat_instance=cat)

# Import optional components based on configuration
if iconfig.get("NEXUS_DATA_FILES", {}).get("ENABLE", False):
    from .callbacks.nexus_data_file_writer import nxwriter_init

    nxwriter = nxwriter_init(RE)

if iconfig.get("SPEC_DATA_FILES", {}).get("ENABLE", False):
    from .callbacks.spec_data_file_writer import init_specwriter_with_RE
    from .callbacks.spec_data_file_writer import newSpecFile  # noqa: F401
    from .callbacks.spec_data_file_writer import spec_comment  # noqa: F401
    from .callbacks.spec_data_file_writer import specwriter

    specwriter.write_new_scan_header = False  # issue #1032
    init_specwriter_with_RE(RE)

# These imports must come after the above setup.
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


def _startup_create_devices_plan():
    """Create ALL devices using a single call to RE on startup."""

    ############################
    # These device files MUST load or startup will stop.
    yield from make_devices(clear=False, file="devices.yml")
    yield from make_devices(clear=False, file="transfocator.yml")

    if host_on_aps_subnet():
        yield from make_devices(clear=False, file="devices_aps_only.yml")
        yield from make_devices(clear=False, file="ad_devices.yml")
        yield from ad_initial_setup()

    ############################
    # These device files can fail gracefully.  Startup will continue.
    device_files = [
        "flight_tube_devices.yml",
        "aerotech_stages_devices.yml",
    ]
    for device_file in device_files:
        try:
            yield from make_devices(clear=False, file=device_file)
        except Exception as excuse:
            print(f"Could not import {device_file!r}: {excuse}")


RE(_startup_create_devices_plan())

pv_registers = oregistry["pv_registers"]

spec_file_name = pv_registers.spec_file
spec_file_name.wait_for_connection()
_fname = spec_file_name.get()

if len(_fname) > 4 and _fname.endswith(".dat"):
    # PV should contain a valid file name
    specwriter.newfile(_fname)
else:
    logger.warning(
        f"SPEC file name {_fname!r} from EPICS PV"
        f" {spec_file_name.pvname!r} is unacceptable."
        "  File name must be of form 'NAME.dat' where NAME"
        " is at least 1 character."
        f"  Using {specwriter.spec_filename}."
            )

# flag1ad = oregistry["flag1ad"]
# flag1ad.stats1.kind = "hinted"
# flag1ad.stats1.total.kind = "hinted"

# flag2ad = oregistry["flag2ad"]
# flag2ad.stats1.kind = "hinted"
# flag2ad.stats1.total.kind = "hinted"

# flag3ad = oregistry["flag3ad"]
# flag3ad.stats1.kind = "hinted"
# flag3ad.stats1.total.kind = "hinted"

from .plans.master_plan import run_measurement_info
from .plans.nexus_acq_eiger_int import setup_eiger_int_series, eiger_acq_int_series
from .plans.nexus_acq_eiger_ext import setup_eiger_ext_trig, eiger_acq_ext_trig
from .plans.nexus_acq_rigaku_zdt import setup_rigaku_ZDT_series, rigaku_acq_ZDT_series, rigaku_zdt_acquire
from .plans.pv_break_test import break_pv
from .plans.rheometer_wait import wait_for_mcr
from .plans.sample_info_unpack import select_sample, gen_folder_prefix
from .plans.scan_8idi import att, x_lup, y_lup, rheo_x_lup, rheo_y_lup, rheo_set_x_lup
from .plans.select_detector import select_detector
from .plans.select_sample_env import select_sample_env
from .plans.shutter_logic import showbeam, blockbeam, shutteron, shutteroff, pre_align, post_align 
from .plans.spec_8IDE import submit_Nexus_DM

from .plans.qnw_plans import set_qnw, find_qnw_index, set_temperature 
from .plans.qnw_plans import set_ramp_rate, set_temperature_with_ramp
from .plans.qnw_plans import set_temperature_env, set_ramp_rate_env

from .plans.nexus_acq_tempus import setup_tempus_int_series, tempus_acquire, tempus_acq_int_series

from id8_user_plans.write_measurement_info import write_measurement_info
