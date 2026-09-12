#!/bin/bash
#
# Start the BLUETELLA scan viewer (scanviewer.py) from anywhere.
#
#   start_scanviewer.sh                  # 1D: this experiment's bluesky data folder
#   start_scanviewer.sh --mesh           # 2D: same folder, X/Y coloured by Z
#   start_scanviewer.sh --dir <folder>   # some other folder
#   start_scanviewer.sh --scan-prefix A  # older scans, named before 2026-09-12
#   start_scanviewer.sh --help           # the viewer's own options
#
# 8-ID writes BLUETELLA's own extended CSV format AND, since 2026-09-12, its
# "S" filename prefix, so the stock 9-ID viewers read our scans with no
# 8-ID-specific code in them and no flags on this command line.
#
# Scans written BEFORE that date took the sample's header letter instead, so a
# single folder can hold A, C and D files. The viewer browses one prefix at a
# time, so reach those with --scan-prefix A (or C, or D). Nothing was renamed.
# See ophyd_scan.SCAN_FILE_HEADER for why the letter is fixed now.
#
# --mesh selects meshviewer.py, for the 2-motor rasters dmesh/mesh. A scan whose
# header declares a raster shape is drawn as a filled grid, each cell binned by
# nearest commanded position, unmeasured cells left blank; anything without a
# declared shape stays a coloured scatter of motor1 vs motor2. Without --mesh you
# get scanviewer.py, the 1D viewer, which is also the right one for d2scan/a2scan
# (two motors along ONE line, so there is no 2D field to grid -- just pick either
# motor as X).
#
# Under heavy attenuation the lambda2M stats columns are all zero -- pick
# tetramm1_sum_all as Y (or Z) to see real signal.
#
# With no --dir, the folder is built the same way the scans themselves build it,
# by reading src/id8_common/configs/experiment.yml:
#
#     <mount_point><cycle_name>/<experiment_name>/data/bluesky
#
# so after a cycle or experiment rolls over this script follows along with no
# edit. --mesh is the only flag this script consumes; everything else you pass
# is handed straight to the viewer it launches.
#
# Tick "Follow latest" in the GUI to track the running scan. The viewer is
# read-only: no EPICS, no ophyd, no Data Management -- it only ever reads files.
#
# Like start_specr_py.sh, this deliberately does NOT source dm.setup.sh:
# that overwrites LD_LIBRARY_PATH, which the viewer's Qt/Tk stack does not want.
#
# VIEWER_ENV picks the conda environment. The default, "bluetella_viewer", is a
# dedicated read-only-viewer environment -- deliberately NOT one of the Bluesky
# environments (8id_bits, 8ide_bits_test, ...). The viewer must never be a
# reason to touch an environment that runs the instrument: upgrading pandas or
# matplotlib for a plotting tweak has no business affecting ophyd or the
# RunEngine, and vice versa. Keep them separate.
#
# To rebuild it (needs the internet, so from a 164.* host such as kouga; conda
# envs live under the shared home, so pearl picks it up with no further work):
#
#     conda create -y -n bluetella_viewer python=3.11 \
#         pandas numpy scipy matplotlib-base tk

# --mesh picks the 2D viewer. It must be stripped before the exec: meshviewer.py's
# argparse does not accept it and would exit "unrecognized arguments: --mesh".
#
# No --scan-prefix is added here any more. 8-ID and BLUETELLA now agree on S, so
# the viewer's own default is right; --scan-prefix is passed straight through for
# the pre-2026-09-12 files that still carry a sample header.
VIEWER_PY="scanviewer.py"
ARGS=()
for a in "$@"; do
    case "$a" in
        --mesh)          VIEWER_PY="meshviewer.py" ;;
        *)               ARGS+=("$a") ;;
    esac
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

VIEWER_DIR=${VIEWER_DIR:-$HOME/Documents/BLUETELLA_9ID}
VIEWER_ENV=${VIEWER_ENV:-bluetella_viewer}
# ~/ophyd, not ~/bluesky. Those split into two checkouts on 2026-09-11 -- dev
# and main -- and only ~/ophyd carries the align scans and their experiment.yml,
# which is also the tree start_ophyd.sh puts on PYTHONPATH. Pointing here at
# ~/bluesky made every launch fail with "experiment.yml not found".
EXPT_YML=${EXPT_YML:-$HOME/ophyd/src/id8_common/configs/experiment.yml}

# --- somewhere to draw ------------------------------------------------------
# Only a note, never a stop: if there is no display Tk says so itself, just less
# clearly ("no display name and no $DISPLAY environment variable").
if [ -z "$DISPLAY" ] && [ -z "$WAYLAND_DISPLAY" ] && [ -z "$QT_QPA_PLATFORM" ]; then
    echo "Note: no DISPLAY set -- reconnect with 'ssh -Y' if the window does not appear."
fi

# --- the viewer -------------------------------------------------------------
if [ ! -f "$VIEWER_DIR/$VIEWER_PY" ]; then
    echo "Error: $VIEWER_DIR/$VIEWER_PY not found."
    echo "       Set VIEWER_DIR if the checkout lives somewhere else, or clone it:"
    echo
    echo "           git clone -b 8id_test \\"
    echo "               github-qzhang234:AdvancedPhotonSource/BLUETELLA_9ID.git \\"
    echo "               $VIEWER_DIR"
    exit 1
fi

# --- where the scan files are ----------------------------------------------
# Skipped entirely when the caller passes their own --dir.
case " $* " in
    *" --dir "*|*" --dir="*|*" -h "*|*" --help "*) DIR_ARGS=() ;;
    *)
        if [ ! -f "$EXPT_YML" ]; then
            echo "Error: $EXPT_YML not found, so the data folder cannot be worked out."
            echo "       Pass one yourself:  $(basename "$0") --dir <folder>"
            exit 1
        fi
        # Flat YAML, one "key: value" per line, values optionally quoted.
        yval() { sed -n "s/^$1:[[:space:]]*[\"']\?\([^\"'#]*\)[\"']\?[[:space:]]*$/\1/p" "$EXPT_YML" | tail -1; }
        MOUNT=$(yval mount_point)
        CYCLE=$(yval cycle_name)
        EXPERIMENT=$(yval experiment_name)
        if [ -z "$MOUNT" ] || [ -z "$CYCLE" ] || [ -z "$EXPERIMENT" ]; then
            echo "Error: could not read mount_point/cycle_name/experiment_name from"
            echo "       $EXPT_YML -- got '$MOUNT' '$CYCLE' '$EXPERIMENT'."
            echo "       Pass a folder yourself:  $(basename "$0") --dir <folder>"
            exit 1
        fi
        DATA_DIR="${MOUNT}${CYCLE}/${EXPERIMENT}/data/bluesky"
        DIR_ARGS=(--dir "$DATA_DIR")

        # A stale mount_point names another station's tree (8IDI vs 8IDE) and is
        # the usual thing left pointing at the previous experiment. Say so here
        # rather than open an empty browser on a folder that does not exist.
        if [ ! -d "$DATA_DIR" ]; then
            echo "Note: $DATA_DIR does not exist."
            echo "      It is built from $EXPT_YML -- check mount_point"
            echo "      ('$MOUNT'), cycle_name ('$CYCLE') and experiment_name"
            echo "      ('$EXPERIMENT'). Opening anyway; use Browse... in the GUI, or --dir."
        else
            echo "Scan folder: $DATA_DIR"
        fi
        ;;
esac

# Scan files live on GPFS. Some analysis workstations (kouga) do not mount it.
if [ ! -d /gdata ]; then
    echo "Note: /gdata is not mounted here, so beamline scan files are not visible."
    echo "      Run this on pearl, or copy a .csv over and use --dir."
fi

# --- conda ------------------------------------------------------------------
# Last, so the checks above still report usefully if the environment is missing.
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    echo "Error: conda.sh not found at '$CONDA_BASE/etc/profile.d/conda.sh'."
    echo "       Check the conda installation, or set CONDA_BASE by hand."
    exit 1
fi
. "$CONDA_BASE/etc/profile.d/conda.sh"
conda deactivate 2>/dev/null

conda activate "$VIEWER_ENV" || {
    echo "Error: conda environment '$VIEWER_ENV' not found."
    echo
    echo "       Rebuild it from a 164.* host (conda needs the internet, and"
    echo "       10.* machines do not have it -- 'ssh 8idiuser@kouga'):"
    echo
    echo "           conda create -y -n $VIEWER_ENV python=3.11 \\"
    echo "               pandas numpy scipy matplotlib-base tk"
    echo
    echo "       Conda environments live under the shared home, so pearl sees it"
    echo "       straight away. Do NOT point VIEWER_ENV at a Bluesky environment."
    echo "       In a pinch, pearl's system python3 has all five dependencies:"
    echo
    echo "           python3 $VIEWER_DIR/$VIEWER_PY $*"
    exit 1
}

# exec: no leftover shell, and ^C goes straight to the viewer.
exec python "$VIEWER_DIR/$VIEWER_PY" "${DIR_ARGS[@]}" "$@"
