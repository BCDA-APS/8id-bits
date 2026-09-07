"""
Consolidated acquisition script for 8-ID detectors.

Per-detector setup/acquire functions and their mode tables are just data and
callables until this module assembles them -- they live in sibling
*_modes.py files (see each file's docstring for what each mode does and when
to use it):
    eiger4m_modes.py  : eiger4M       -- "Internal Series", "Internal Enable",
                                          "External Series", "External Enable"
    lambda2m_modes.py : lambda2M      -- "Internal", "External"
    rigaku3m_modes.py : rigaku3M      -- "ZDT2bit", "ZDT4bit", "ZDT8bit"
                        rigaku3M_ftf   -- "ZDT2bit", "ZDT4bit", "ZDT8bit"
                        rigaku3M_epics -- "EPICS"

This module assembles those into ACQ_MODES and provides the shared
acquisition loop (det_acq_series) and its cleanup helper. det_acq_series() is
the usual entry point but not the only caller: dual_acq_series() pulls each
leg's setup callable straight out of ACQ_MODES and runs it, and reuses
cleanup_acquisition().

scan_8id.py and scan_8id_dev.py do `from id8_common.plans.acquire.ad_acq
import *` and call setup_eiger_internal(...) etc. directly, so this module
re-exports everything from the per-detector files the same way -- no __all__
needed: unlike named imports, ruff's --fix leaves `import *` alone (it can't
tell what's "unused" behind a star import), so nothing here is at risk of
being silently stripped on the next pre-commit run. The per-detector files
have no __all__ either, so all of their public names come through; acq_helpers
does have one, and it is what decides which of ITS names (gen_folder_prefix,
get_ophyd_object, sample_mesh_move, ...) reach here and, through startup.py,
the prompt.
"""

import os
import time as ttime
from datetime import datetime

from id8_common.expt_config import expt
from id8_common.plans.acquire import acq_helpers
from id8_common.plans.acquire.acq_helpers import *
from id8_common.plans.acquire.eiger4m_modes import *
from id8_common.plans.acquire.lambda2m_modes import *
from id8_common.plans.acquire.rigaku3m_modes import *
from id8_common.plans.set.shutter_att import blockbeam
from id8_common.plans.set.shutter_att import post_align
from id8_common.plans.set.shutter_att import shutteroff
from id8_common.utils.dm_util import dm_run_job
from id8_common.utils.dm_util import dm_setup
from id8_common.utils.nexus_utils import create_nexus_format_metadata

# =============================================================================
# Mode table
# =============================================================================

ACQ_MODES = {
    "eiger4M": EIGER4M_MODES,
    "lambda2M": LAMBDA2M_MODES,
    "rigaku3M": RIGAKU3M_BIN_MODES,
    "rigaku3M_ftf": RIGAKU3M_FTF_MODES,
    "rigaku3M_epics": RIGAKU3M_EPICS_MODES,
}


# =============================================================================
# Cleanup helper
# =============================================================================

def cleanup_acquisition(det=None, mode_info=None, metadata_fname=None):
    """
    Close the shutter, stop softglue and the detector, write metadata.

    metadata_fname:
        NeXus metadata path of the measurement that was still in flight. The partial
        dataset left behind by the abort keeps its metadata this way. Pass None (the
        default) once the metadata has already been written, so an abort later in the
        same repeat does not overwrite a complete file.
    """
    # First, and in its own try, so no later failure can leave the beam on the sample.
    try:
        blockbeam()
    except Exception as e:
        print(f"Could not block beam during cleanup: {e}")

    if mode_info is not None and "softglue" in mode_info.get("required_devices", []):
        try:
            softglue = get_connected_device("softglue")
            softglue.stop_pulses.put("1!")
        except Exception:
            pass

    if det is not None:
        try:
            det.cam.acquire.put(0)
        except Exception:
            pass
        if hasattr(det, "hdf1"):
            try:
                det.hdf1.capture.put(0)
            except Exception:
                pass

    # Last, so the detector is stopped and the data file is closed before it is read.
    if metadata_fname is not None and det is not None:
        try:
            print(f"Writing metadata for interrupted measurement: {metadata_fname}")
            create_nexus_format_metadata(metadata_fname, det=det)
        except Exception as e:
            print(f"Could not write metadata {metadata_fname}: {e}")


# =============================================================================
# Main user-facing acquisition function
# =============================================================================

def det_acq_series(wait_time=0, hooks=None):
    """Run repeated detector acquisitions from the current run state.

    Parameters come from ``expt`` (id8_common/expt_config.py): the protocol
    values are loaded from measurement_info.yaml by
    master_plan.run_measurement(), the per-experiment ones from
    configs/experiment.yml, and measurement_num from 8ideSoft:Reg1.
    ``expt.file_name`` is an output, not an input -- this function builds it
    per repeat and writes it back. See docs/running-measurements.md for
    sample-motion behavior, hooks, and file-naming rules.

    To drive this directly, without master_plan, set the run state first::

        expt.set_measurement({"detector": "eiger4M", "mode": "Internal Series",
                              "acq_time": 0.01, "acq_period": 0.01,
                              "num_frames": 100, "num_repeats": 1,
                              "sample_move": "no", "qmap_file": "..."},
                             {"header": "A", "sample_name": "Test"})
        det_acq_series()

    Error handling, which matters to anything calling this in a loop: Ctrl+C is
    turned into a RuntimeError (so an upstream ``except KeyboardInterrupt``
    will NOT catch it), while any other exception is cleaned up after, printed,
    and swallowed -- this function then returns normally. A caller therefore
    cannot tell a failed measurement from a successful one by exception alone.
    Either way cleanup_acquisition() runs first and attempts to close the
    shutter before anything else.
    """
    det = None
    mode_info = None
    metadata_fname = None
    try:
        # The [:] overwrites the contents of the existing list instead of
        # rebinding the name. acq_helpers is star-imported here and, through
        # startup.py, into the interactive namespace, so several names are
        # already bound to this one list object; a plain assignment would
        # replace only acq_helpers' own attribute and leave the rest pointing
        # at the old, empty list.
        acq_helpers.active_hooks[:] = load_hooks(hooks)

        post_align()
        shutteroff()

        workflowProcApi, dmuser = dm_setup()

        detector = expt.det_name
        mode = expt.det_mode

        acq_time = expt.acq_time
        acq_period = expt.acq_period
        num_frames = expt.num_frames
        num_reps = expt.num_repeats
        num_segments = expt.num_segments
        trigger_period = expt.trigger_period

        mode_info = ACQ_MODES[detector][mode]

        for device_name in mode_info["required_devices"]:
            get_connected_device(device_name)

        # Most ACQ_MODES keys are themselves device names, but a few are aliases
        # for different ways of driving one detector -- rigaku3M_ftf and
        # rigaku3M_epics both drive the rigaku3M device -- and those table
        # entries carry a "hardware_device" naming the real one.
        det = get_connected_device(mode_info.get("hardware_device", detector))
        setup_func = mode_info["setup"]
        acquire_func = mode_info["acquire"]

        folder_prefix = gen_folder_prefix()
        file_header = f"{folder_prefix}_f{num_frames:06d}"

        for rep in range(num_reps):
            ttime.sleep(wait_time)

            sample_mesh_move()

            file_name = f"{file_header}_r{rep + 1:05d}"

            # Persistent session state: written to the `persistent:` block of
            # state/run_state.yml, so a GUI or a shell script can see which
            # measurement is running without a live session. It used to be
            # 8ideSoft:StrReg8 -- that Component still exists on pv_registers,
            # but nothing writes it any more, so don't caget it.
            # Bare name, no extension and no path -- the folder is rebuildable
            # from mount_point/cycle_name/experiment_name/use_subfolder, and
            # each detector appends its own suffix (.h5, .bin, _metadata.hdf).
            expt.file_name = file_name

            # Only the parameters a mode declares are passed, so a setup
            # function takes exactly the arguments it uses. Add a
            # "needs_<name>" flag to a mode table entry to extend this.
            setup_kwargs = {
                "acq_time": acq_time,
                "num_frames": num_frames,
                "file_header": file_header,
                "file_name": file_name,
            }

            # needs_acq_period is spelled with plain [...] because every mode
            # entry must define it; the two below use .get(..., False) because
            # they are optional and most entries leave them out.
            if mode_info["needs_acq_period"]:
                setup_kwargs["acq_period"] = acq_period

            if mode_info.get("needs_num_segments", False):
                setup_kwargs["num_segments"] = num_segments

            if mode_info.get("needs_trigger_period", False):
                setup_kwargs["trigger_period"] = trigger_period

            metadata_fname = setup_func(**setup_kwargs)

            time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"\n{time_now}, Starting measurement {file_name}")

            acquire_func()

            time_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"{time_now}, Complete measurement {file_name}")

            print(f"{time_now}, Writing metadata, {file_name}")

            # ---------------------------------------------------------------
            # TEST 2026-09-04: the NeXus writer is selectable, so the DEFAULT
            # behaviour is unchanged -- a normal session uses our writer and
            # needs nothing extra on PYTHONPATH.
            #
            # Set ID8_NEXUS_WRITER=mc to use Miaoqi Chu's mc_refact writer
            # instead (requires PYTHONPATH to include the nexus_xpcs_aps clone;
            # the launcher run_mc_writer_test.sh does that).
            #
            # TO REMOVE THE TEST ENTIRELY: replace this whole block with the
            # single line `create_nexus_format_metadata(metadata_fname, det=det)`
            # and delete utils/nexus_utils_mc.py + utils/xpcs_schema_mc.py.
            #
            # cleanup_acquisition() (the abort path, above) always uses OUR
            # writer, so an aborted run falls back to known-good code.
            if os.environ.get("ID8_NEXUS_WRITER") == "mc":
                from id8_common.utils.nexus_utils_mc import (
                    create_nexus_format_metadata_mc,
                )

                print("    [writer] mc_refact (nexus_xpcs_aps)")
                create_nexus_format_metadata_mc(metadata_fname, det=det)
            else:
                print("    [writer] id8_common (default)")
                create_nexus_format_metadata(metadata_fname, det=det)
            # ---------------------------------------------------------------

            # Complete, so an abort during DM submission must not rewrite it.
            metadata_fname = None

            print(f"{time_now}, Submitting to DM, {file_name}")
            dm_run_job(workflowProcApi, dmuser, file_name)

    except KeyboardInterrupt:
        cleanup_acquisition(det, mode_info, metadata_fname)
        raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).")
    except Exception as e:
        # DELIBERATE, do not "fix" this into a raise. A measurement_info.yaml run
        # is often an overnight queue of dozens of measurements; one that fails --
        # a detector hiccup, a full disk, a DM submission that bounces -- must not
        # take the rest of the night with it. So the error is reported and this
        # measurement is abandoned, and master_plan.run_measurement() goes on to
        # the next one.
        #
        # The cost, which is real: master_plan cannot tell a failed measurement
        # from a successful one, so a run can "complete" with gaps in it. Read the
        # printed errors before trusting a night's output. Decided by the owner
        # 2026-09-07 with that tradeoff understood.
        cleanup_acquisition(det, mode_info, metadata_fname)
        print(f"Error occurred during measurement: {e}")
    finally:
        pass

# =============================================================================
# Example use cases -- see docs/running-measurements.md for
# run-state-driven det_acq_series() examples (Eiger/Lambda/Rigaku, with and
# without sample motion).
# =============================================================================
