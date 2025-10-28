#!/bin/bash

# # Start the bluesky queueserver.

# Standard BCDA setup at APS defines this pointing to APSshare.
# We choose here to get PyEpics from the conda environment.
unset PYEPICS_LIBCA # TODO: Need to fix this more permenantly
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # This allows for command to be run from anywhere

#--------------------
# change the program defaults here
SCRIPT_DIR=$(dirname $(readlink -f "${0}"))
DATABROKER_CATALOG=$(grep DATABROKER_CATALOG ${SCRIPT_DIR}/../src/instrument/configs/iconfig.yml  | awk '{print $NF}') #TODO: De-hardcode
echo "Using catalog ${DATABROKER_CATALOG}"


export QS_SERVER_HOST=$(hostname)  # or host (that passes $(hostname) test below)
export QS_UPDATE_PLANS_DEVICES=ENVIRONMENT_OPEN
export QS_USER_GROUP_PERMISSIONS_FILE="${SCRIPT_DIR}/../src/instrument/configs/user_group_permissions.yaml" #TODO: De-hardcode
export QS_EXISTING_PLANS_DEVICES_FILE="${SCRIPT_DIR}/../src/instrument/configs/" #TODO: De-hardcode
export QS_USER_GROUP_PERMISSIONS_RELOAD=ON_STARTUP
# export STARTUP_DIR="${SCRIPT_DIR}/../src/instrument/"
export STARTUP_SCRIPT="${SCRIPT_DIR}/bs_qs_startup.py"

# REDIS_ADDR is __always__ localhost.
# Override if it is not, but you may encounter access issues.
export REDIS_ADDR=localhost
#--------------------

# QS and redis must be on the same workstation
if [ "$(hostname)" != "${QS_SERVER_HOST}" ]; then
    echo "Must run queueserver on ${QS_SERVER_HOST}.  This is $(hostname)"
    exit 1
fi

# echo "CONDA_EXE = '${CONDA_EXE}'"
if [ ! -f "${CONDA_EXE}" ]; then
    echo "No 'conda' command available."
    exit 1
fi

################### TODO: See if needs to be kept
# In GitHub Actions workflow,
# # $ENV_NAME is an environment variable naming the conda environment to be used
# if [ -z "${ENV_NAME}" ] ; then
#     ENV_NAME="${CONDA_ENVIRONMENT}"
# fi

# if [ "${CONDA_DEFAULT_ENV}" != "${ENV_NAME}" ]; then
#     echo "Activating conda environment ${ENV_NAME}"
#     CONDA_BASE="$(dirname $(dirname $(readlink -f ${CONDA_EXE})))"
#     source "${CONDA_BASE}/etc/profile.d/conda.sh"
#     conda activate "${ENV_NAME}"
# fi
# echo "conda env list = $(conda env list)"
###################

# #--------------------
# echo "Environment: $(env | sort)"
# echo "------"
# echo "CONDA_ENVIRONMENT=${CONDA_ENVIRONMENT}"
# echo "CONDA=${CONDA}"
# echo "DATABROKER_CATALOG=${DATABROKER_CATALOG}"
# echo "QS_SERVER_HOST=${QS_SERVER_HOST}"
# echo "QS_UPDATE_PLANS_DEVICES=${QS_EXISTING_PLANS_DEVICES_FILE}"
# echo "QS_USER_GROUP_PERMISSIONS_FILE=${QS_USER_GROUP_PERMISSIONS_FILE}"
# echo "QS_USER_GROUP_PERMISSIONS_RELOAD=${QS_USER_GROUP_PERMISSIONS_RELOAD}"
# echo "REDIS_ADDR=${REDIS_ADDR}"
# echo "SHELL_SCRIPT_NAME=${SHELL_SCRIPT_NAME}"
# echo "STARTUP_DIR=${STARTUP_DIR}"
# #--------------------

# Start the bluesky queueserver (QS)
start-re-manager \
    --redis-addr "${REDIS_ADDR}" \
    --startup-script "${STARTUP_SCRIPT}" \
    --update-existing-plans-devices "${QS_UPDATE_PLANS_DEVICES}" \
    --user-group-permissions "${QS_USER_GROUP_PERMISSIONS_FILE}" \
    --user-group-permissions-reload "${QS_USER_GROUP_PERMISSIONS_RELOAD}" \
    --existing-plans-devices "${QS_EXISTING_PLANS_DEVICES_FILE}" \
    --zmq-publish-console ON \
    --keep-re
