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
if iconfig.get("SPEC_DATA_FILES", {}).get("ENABLE", False):
    from .callbacks.demo_spec_callback import init_specwriter_with_RE
    from .callbacks.demo_spec_callback import newSpecFile  # noqa: F401
    from .callbacks.demo_spec_callback import spec_comment  # noqa: F401
    from .callbacks.demo_spec_callback import specwriter  # noqa: F401

    init_specwriter_with_RE(RE)

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
make_devices(clear=False, file="devices_aps_only.yml", device_manager=instrument)
make_devices(clear=False, file="ad_devices.yml", device_manager=instrument)

# from apstools.devices import load_devices_from_yaml
# from id8_common.utils.misc import ioc_alive

# RIGAKU_TEST_PV = "8idRigaku3m:cam1:Manufacturer_RBV"
# EIGER_TEST_PV = "8idEiger4m:cam1:Manufacturer_RBV"

# if ioc_alive(RIGAKU_TEST_PV):
#     print("Rigaku3M IOC is up — loading detector")
#     load_devices_from_yaml("ad_devices.yml", oregistry)

# else:
#     print("Rigaku3M IOC not reachable — skipping detector load")


# LivePlot for area-detector ROI sum
try:
    from bluesky.callbacks import LivePlot
    from apsbits.utils.helper_functions import running_in_queueserver

    if not running_in_queueserver():
        # Choose the detector name here; 'eiger4M', 'rigaku3M', etc.
        det_name = "eiger4M"  # change to the detector you want to plot
        det = oregistry.get(det_name)
        if det is None:
            print(f"LivePlot: detector {det_name!r} not in oregistry; skipping LivePlot.")
        else:
            plugin = getattr(det, "roi1", None) or getattr(det, "stats1", None) or getattr(det, "image", None)
            if plugin is None:
                print(f"LivePlot: detector {det_name} has no roi1/stats1/image plugin attribute; available: {det.component_names}")
            else:
                # Find a numeric field to plot. Common choices: 'sum', 'total', 'value', 'mean_value'
                candidates = ["sum", "total", "value", "mean_value", "total_value"]
                field = None
                for c in candidates:
                    if c in plugin.component_names or hasattr(plugin, c):
                        field = c
                        break
                # fallback: look for any readable attribute that returns a scalar
                if field is None:
                    for name in plugin.component_names:
                        try:
                            val = getattr(plugin, name).get()
                            # scalar numeric test
                            if isinstance(val, (int, float)):
                                field = name
                                break
                        except Exception:
                            continue

                if field is None:
                    print(f"LivePlot: couldn't find numeric field in {det_name}.roi1; plugin components: {plugin.component_names}")
                else:
                    # Use 'seq_num' as x axis; for position-based scans you can use the motor signal instead
                    lp = LivePlot(field, plugin, x="seq_num")
                    RE.subscribe(lp)
                    print(f"LivePlot subscribed: {det_name}.roi1.{field} vs seq_num")
    else:
        print("Running in queueserver mode; skipping LivePlot setup.")
except Exception as e:
    print("LivePlot setup failed:", e)

if host_on_aps_subnet():
    make_devices(clear=False, file="devices_aps_only.yml", device_manager=instrument)

pv_registers = oregistry["pv_registers"]

# Setup baseline stream with connect=False is default
# Devices with the label 'baseline' will be added to the baseline stream.
setup_baseline_stream(sd, oregistry, connect=False)

from .plans.sim_plan import sim_count_plan  # noqa: E402, F401
from .plans.sim_plan import sim_print_plan  # noqa: E402, F401
from .plans.sim_plan import sim_rel_scan_plan  # noqa: E402, F401

# from .plans.master_plan import run_measurement_info
# from .plans.select_sample_env import select_sample_env
from .plans.sample_info_unpack import *
# from .plans.select_detector import *
from .plans.scan_8ide import *
# from .plans.qnw_plans import *
from .plans.nexus_acq_eiger_int import *
# from .plans.nexus_acq_rigaku_zdt import *

psic = oregistry["psic"]
import hklpy2
from hklpy2.user import *
from hklpy2.user import set_diffractometer
set_diffractometer(psic)

# from id8_common.utils.misc import stream_rois
# stream_rois(eiger4M, stats_nums=(1,2,3), fields=("total",))
