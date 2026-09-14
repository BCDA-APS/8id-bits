#!/bin/bash
#
# Start a Bluesky session against the MIGRATION checkout, ~/bluesky_trio.
#
# This is a copy of ~/bin/start_bluesky.sh with one difference: the import path
# points at ~/bluesky_trio/src instead of ~/bluesky/src. ~/bluesky is untouched,
# so `start_bluesky.sh` still starts exactly the session it always did.
#
# What this session has that ~/bluesky does not:
#   * startup does not fire an Eiger exposure (ad_setup no longer primes), and a
#     guard reports any EPICS write attempted during startup
#   * run_trio_measurement_info() / dry_run_trio_measurement_info()
#   * measurement plans may be written as CSV as well as YAML
#
# Both sessions can be run, but NOT at the same time against the same detector.

source /home/dm_id/etc/dm.setup.sh
conda deactivate

CONDA_BASE=$(conda info --base)

if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    . "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "Error: conda.sh not found at $CONDA_BASE/etc/profile.d/conda.sh"
    exit 1
fi

# Same environment as start_bluesky.sh -- 8ide_bits_test, NOT 8id_bits. See the
# comment in that script for why the two are not interchangeable.
conda activate 8ide_bits_test || { echo "Error: Failed to activate '8ide_bits_test'."; exit 1; }

cd ~/bluesky_trio || { echo "Error: ~/bluesky_trio not found."; exit 1; }

# Pin the import path to THIS checkout.
#
# `cd` alone is NOT enough: the editable install __editable__.id8_bits-0.0.1.pth
# in 8ide_bits_test hardcodes /home/beams10/8IDIUSER/bluesky/src, so
# `import id8_common` resolves THERE whatever directory you start from.
# PYTHONPATH lands earlier in sys.path than site-packages, so this wins.
#
# id8_common is a NAMESPACE package (no __init__.py), so two src roots MERGE
# rather than one shadowing the other. Do NOT add ~/bluesky/src or ~/ophyd/src
# here as well: a module missing from the first would silently fall back to the
# other checkout. One src root only.
export PYTHONPATH="$HOME/bluesky_trio/src${PYTHONPATH:+:$PYTHONPATH}"

# Add robocart EPICS access to whichever machine is running Bluesky
export EPICS_CA_ADDR_LIST="$EPICS_CA_ADDR_LIST 164.54.116.40"

ipython -i -c "from id8_common.startup import *" || { echo "Error: Failed to start ipython."; exit 1; }
