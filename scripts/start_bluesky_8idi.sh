#!/bin/bash

# Find the base conda installation path (adjust if necessary)
# This command tries to locate your conda.sh script reliably

source /home/dm_id/etc/dm.setup.sh
conda deactivate

CONDA_BASE=$(conda info --base)

# Source the conda.sh script directly from the base installation
# This is crucial for making 'conda activate' available in non-interactive shells
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    . "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "Error: conda.sh not found at $CONDA_BASE/etc/profile.d/conda.sh"
    echo "Please verify your conda installation path and adjust CONDA_BASE in the script."
    exit 1
fi

# Change to the bluesky directory
cd ~/bluesky || { echo "Error: ~/bluesky directory not found."; exit 1; }

# Activate the conda environment
conda activate 8idi_bits_test || { echo "Error: Failed to activate conda environment '8idi_bits'."; exit 1; }

# Add robocart EPICS Access to whichever machine running Bluesky
export EPICS_CA_ADDR_LIST="$EPICS_CA_ADDR_LIST 164.54.116.40"

# Start ipython with the specified startup command
ipython -i -c "from id8_i.startup import *" || { echo "Error: Failed to start ipython."; exit 1; }

