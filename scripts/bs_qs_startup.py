"""
Configure for data collection using bluesky-queueserver.
"""

import getpass
import logging
import os
import socket

import apstools
import bluesky
import bluesky_queueserver
import databroker
import epics
import h5py
import matplotlib
import numpy
import ophyd
import pyRestTable
import spec2nexus
from aps_8id_bs_instrument.callbacks import *  # noqa
from aps_8id_bs_instrument.devices import *  # noqa
from aps_8id_bs_instrument.plans import *  # noqa
from aps_8id_bs_instrument.utils.catalog import load_catalog
from aps_8id_bs_instrument.utils.epics_tools import set_control_layer
from aps_8id_bs_instrument.utils.iconfig_loader import iconfig
from aps_8id_bs_instrument.utils.run_engine import run_engine

# guides choice of module to import cat

logger = logging.getLogger(__name__)

logger.info(__file__)
print(__file__)

HOSTNAME = socket.gethostname() or "localhost"
USERNAME = getpass.getuser() or "queueserver user"

# useful diagnostic to record with all data
versions = dict(
    apstools=apstools.__version__,
    bluesky=bluesky.__version__,
    bluesky_queueserver=bluesky_queueserver.__version__,
    databroker=databroker.__version__,
    epics=epics.__version__,
    h5py=h5py.__version__,
    matplotlib=matplotlib.__version__,
    numpy=numpy.__version__,
    ophyd=ophyd.__version__,
    pyRestTable=pyRestTable.__version__,
    spec2nexus=spec2nexus.__version__,
)

#################
# Connect with our mongodb database
catalog_name = iconfig.get("DATABROKER_CATALOG", "training")
try:
    # We don't actually run this part because it is unable to find a yaml
    # corresponding to the catalog name. For this to execute the yaml has to be
    # store in:
    # import databroker; print(databroker.catalog_search_path())
    # https://blueskyproject.io/databroker/reference/configuration.html?highlight=search%20path
    cat = load_catalog(catalog_name)
    logger.info("using databroker catalog '%s'", cat.name)  # FIXME: now showing on console
except KeyError:
    cat = databroker.temp().v2
    logger.info("using TEMPORARY databroker catalog '%s'", cat.name)  # FIXME: now showing on console

print(f"cat = {cat!r}")  # TODO: remove this diagnostic after FIXME items resolved

# Set up a RunEngine.
RE = run_engine(cat=cat)  # TODO: No idea why when I put cat it no longer works can you please check?
# RE = run_engine()

set_control_layer("PyEpics")
###################

RE.md["databroker_catalog"] = cat.name
RE.md["login_id"] = USERNAME + "@" + HOSTNAME
RE.md.update(iconfig.get("RUNENGINE_METADATA", {}))
RE.md["versions"] = versions
RE.md["pid"] = os.getpid()
# if scan_id_epics is not None:
#     RE.md["scan_id"] = scan_id_epics.get()

# Set up SupplementalData.
sd = bluesky.SupplementalData()
RE.preprocessors.append(sd)
