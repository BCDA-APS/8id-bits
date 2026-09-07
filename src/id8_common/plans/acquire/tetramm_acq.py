"""
Simple, modular Ophyd scripts for users.

One plan: tetramm_acq_series(), which records a series of TetrAMM
picoammeter readings straight to an HDF file. Plain Ophyd puts and gets, no
RunEngine -- callable from the prompt or from a user script.
"""

from datetime import datetime
import time

# Load-bearing star import, not a stray one: this module has no __all__, so
# startup.py's "from .plans.acquire.tetramm_acq import *" is what carries
# att(), showbeam(), blockbeam() and the rest of shutter_att to the
# scientist's prompt. See the __all__ note in master_plan.py. Removing it
# takes those names off the prompt even though nothing here calls them.
from id8_common.plans.set.shutter_att import *
from id8_common.expt_config import expt
from id8_common.registry import get_connected_device

# The TetrAMMs are resolved per call by get_connected_device() -- see
# id8_common/registry.py. The tetramm1/2/4 and huber module-level captures
# that used to sit here were unused by this module; every device is still
# available as a bare name at the prompt, put there by the device loader.
# Path settings come from expt (configs/experiment.yml).

DEFAULT_TETRAMM = "tetramm3"


def tetramm_acq_series(
    det=None,
    filename=None,
    num_capture=1,
):
    """Record num_capture TetrAMM samples to an HDF file, then wait for it to close.

    Args:
        det: the TetrAMM to read. Defaults to DEFAULT_TETRAMM, looked up when
            the function runs (see below).
        filename: base name of the output file, no path and no extension.
            The rest of the path is built from the experiment config.
        num_capture: how many samples the HDF plugin should collect.

    The beam is NOT opened: the showbeam()/blockbeam() calls are commented out,
    so whatever the shutter was doing when you called this is what it keeps
    doing.
    """
    # det=None rather than det=tetramm3: a default argument is evaluated once,
    # at import time, so the old form captured whatever tetramm3 was then --
    # None if it had been skipped as offline. Resolved here instead.
    if det is None:
        det = get_connected_device(DEFAULT_TETRAMM)

    cycle_name = expt.cycle_name
    exp_name = expt.experiment_name
    mount_point = expt.mount_point
    
    # No slash after mount_point on purpose -- it already ends with one. This
    # is the same layout get_common_file_path() builds in acq_helpers.py:
    #   <mount_point><cycle_name>/<experiment_name>/data/<filename>
    # The whole path goes into hdf1.file_name; unlike the area detectors,
    # nothing sets hdf1.file_path separately here.
    full_filename = f"{mount_point}{cycle_name}/{exp_name}/data/{filename}"

    det.hdf1.file_name.put(full_filename)
    det.hdf1.num_capture.put(num_capture)

    # showbeam()
    time.sleep(0.1)
    det.hdf1.capture.put(1)

    # Wait for the plugin to collect its samples and close the file: capture
    # drops to 0 by itself once num_capture have been written. Unbounded on
    # purpose for now -- if the TetrAMM stops sending, this waits forever.
    while True:
        time.sleep(0.1)
        det_status = det.hdf1.capture.get()
        if det_status == 1:
            time.sleep(0.1)
        if det_status == 0:
            break
    # blockbeam()
