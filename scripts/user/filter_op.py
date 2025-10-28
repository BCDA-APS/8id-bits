"""
This one is to operate the filters (FL2 in the 8-ID-E hutch through
  the EPICS interface pydev_filter_12
The filter PV prefix can be changed to use for a different set of
  hardware.

The porper way to do this, of course, is to use opyd layer but this
  is a quick and dirty hack.

This one can set attnenuation or transmission, though the EPICS layer
  takes only attenuation, which should be >= 1 always.
For transmission, it would be always <= 1.

Attnuation and transmission can also be read back.

Another function is to enable/disable the blades, most likely the whole
  bank altogether, when the next hutch is running.

Note diable a blade means pulling out the blade and locking it in the
  out position.

============ Zhan Zhang, 10/21/2024, =============

To do list:
        - check the current setting of the blade, if the same as the target,
          do not operate.  The EPICS operates the blades every time it gets
          a input, even if it is the same value. And it first insert all
          blades, then remove those not needed.
          Or maybe this should be in the EPICS layer.

"""

import epics as pe

EMPTY_LIST = []

# Define some global variables for PV access and default values
# If no value given, use this default value for attnuation
ATTN_DEFAULT = 1e6
# PV prefix for the filter unit
FILTER_PV_PREFIX = "8idPyFilter:FL2:"
# contro PV
FILTER_CONTROL = "attenuation"
# readback PV
FILTER_RDBK = "attenuation_actual"


def filter_attn(attn=ATTN_DEFAULT):
    """
    Set the filter based on the given attenuation factors.
    Usage:
            filter_attn()
            filter_attn(attn)
    Parameters:
            attn (optional): attenuation factor (>=1)
                    default value: defined above with ATTN_DEFAULT
    Example:
            filter_attn(200)
            filter_attn(1e10)

    """
    if attn < 1:
        print(" Attenuation factor should be >= 1.")
        print(f" Default attenuation factor {ATTN_DEFAULT:5.3e} is used.")
        attn = ATTN_DEFAULT
    else:
        print(f" Set attenuation factor to {attn:5.3e}.")
    pe.caput((FILTER_PV_PREFIX + FILTER_CONTROL), attn)


def filter_read():
    """
    Read back the current filter setting as attenuation factor
    Usage:
            filter_read()
    Parameters:
            none
    Example:
            filter_read()

    """
    curr_attn = pe.caget(FILTER_PV_PREFIX + FILTER_RDBK)
    print(f" Current attenuation is {curr_attn:5.3e}.")
    print(f"         transmission is {1/curr_attn:5.3e}")
    return pe.caget(FILTER_PV_PREFIX + FILTER_RDBK)


def filter_trans(trans=1 / ATTN_DEFAULT):
    """
    Set the filter based on the given transmission factors.
    NOTE:
            transmission = 1/attenuation
    Usage:
            filter_trans()
            filter_trans(trans)
    Parameters:
            trans (optional): transmission factor (<=1)
                    default value: defined with 1/ATTN_DEFAULT
    Example:
            filter_trans(0.01)
            filter_trans(1e-10)

    """
    if (trans > 1) | (trans <= 0):
        print(" Transmission factor should be <= 1. ")
        print(" Default transmission factor is used.")
        attn = ATTN_DEFAULT
    else:
        attn = 1 / trans
        print(f" Set transmission factor to {1/attn:5.3e}.")
    pe.caput((FILTER_PV_PREFIX + FILTER_CONTROL), attn)


def attn(attn):
    """
    Alias for filter_attn()
    """
    filter_attn(attn)


def FL2_enable(enable_flag="Enabled", blade_list=EMPTY_LIST):
    """
    Enable or disable filter blade(s)
    Usage:
            FL2_enable()
            FL2_enable(enable_flag = "Enabled", blade_list=[])
    Parameters:
            enable_flag (optional): "Enabled" or "Disabled"
                    default value: "Enabled"
            blade_list (optional): list of blades.
                    default value: empty list [], meaning all blades
    Example:
            FL2_enable()
            FL2_enable("Disabled", [1,2,3])
            FL2_enable("Eabled")
    """
    if not blade_list:
        blade_list = range(1, 13)
    if enable_flag != "Disabled":
        enable_flag = "Enabled"
    for blade_num in blade_list:
        if blade_num in range(1, 13):
            pe.caput((FILTER_PV_PREFIX + f"filter{blade_num:02d}" + "_Enable"), enable_flag)


def FL2_disable():
    """
    Disable the filter blades. A shortcut.
    """
    FL2_enable("Disabled")
