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
from apsbits.utils.aps_functions import host_on_aps_subnet
from apsbits.utils.config_loaders import get_config
from apsbits.utils.config_loaders import load_config
from apsbits.utils.helper_functions import register_bluesky_magics
from apsbits.utils.helper_functions import running_in_queueserver

# from id8_i.plans.ad_setup_plans import ad_initial_setup

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
    from .callbacks.demo_spec_callback import init_specwriter_with_RE
    from .callbacks.demo_spec_callback import newSpecFile  # noqa: F401
    from .callbacks.demo_spec_callback import spec_comment  # noqa: F401
    from .callbacks.demo_spec_callback import specwriter

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
        # yield from ad_initial_setup()

    ############################
    # These device files can fail gracefully.  Startup will continue.
    device_files = [
        # "flight_tube_devices.yml",
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

# from .plans.nexus_acq_tempus import setup_tempus_int_series, tempus_acquire, tempus_acq_int_series


# from .plans.nexus_acq_rigaku_zdt import setup_rigaku_ZDT_series, rigaku_acq_ZDT_series, rigaku_zdt_acquire
