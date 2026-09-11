"""
12-bank Filters from A-V-S

Device uses PyDevice for attenuation calculation and filter configuration

    Parameters
    ==========
    prefix:
      EPICS prefix required to communicate with filter IOC, ex: "100idPyFilter:FL2:"
    translation_motor:
      The motor record PV controlling the lateral translation of the filter system

"""

from ophyd import Component as Cpt
from ophyd import Device
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO
from ophyd import FormattedComponent as FCpt
from ophyd import PVPositioner

#: Seconds before a filter move gives up. See the note on done_value below.
FILTER_MOVE_TIMEOUT = 60.0


# filterBusy is an enum: [0] Done, [1] Changing. Ophyd's PVPositioner defaults
# to done_value=1, which inverts the sense of all three positioners below:
# at rest (0) the positioner believes it is moving, and the first transition to
# Changing (1) satisfies done and completes the move. `.move()` therefore
# returns the instant the blades start, not when they arrive. Diagnosed
# 2026-09-11, after att(10) reported success at 1.78x and a measurement ran at
# the wrong attenuation.
#
# A timeout is set alongside it because done_value=0 introduces the opposite
# failure: _done_moving() only fires on the Changing -> Done transition, so a
# request that changes nothing never raises filterBusy and `.move()` would
# otherwise wait forever. shutter_att.att() polls filterBusy directly and does
# not rely on either.
class filter_index(PVPositioner):
    """
    filter index; increasing index, increasing attenuation
    """

    readback = Cpt(EpicsSignalRO, "sortedIndex_RBV")
    setpoint = Cpt(EpicsSignal, "sortedIndex")
    done = Cpt(EpicsSignalRO, "filterBusy")
    done_value = 0


class filter_atten(PVPositioner):
    """
    filter attenuation positioner
    """

    readback = Cpt(EpicsSignalRO, "attenuation_actual")
    setpoint = Cpt(EpicsSignal, "attenuation")
    done = Cpt(EpicsSignalRO, "filterBusy")
    done_value = 0


class filter_trans(PVPositioner):
    """
    filter transmission positioner
    """

    readback = Cpt(EpicsSignalRO, "transmission_RBV")
    setpoint = Cpt(EpicsSignal, "transmission")
    done = Cpt(EpicsSignalRO, "filterBusy")
    done_value = 0


class AVSfilters(Device):
    """
    Ophyd device for avs filters
    """

    def __init__(
        self,
        prefix: str,
        # translation_motor: str,
        *args,
        **kwargs,
    ):
        """
        Initialize the AVS filter device
        """
        # self._translation = translation_motor

        super().__init__(prefix, *args, **kwargs)

    index = FCpt(filter_index, "{prefix}", timeout=FILTER_MOVE_TIMEOUT)
    attenuation = FCpt(filter_atten, "{prefix}", timeout=FILTER_MOVE_TIMEOUT)
    transmission = FCpt(filter_trans, "{prefix}", timeout=FILTER_MOVE_TIMEOUT)
    # translation = FCpt(EpicsMotor, "{_translation}", labels={'motors'})

    binary_crl1_config = Cpt(EpicsSignalRO, "filterConfig", kind="hinted")
    bw_crl1_config = Cpt(EpicsSignalRO, "filterConfig_BW")
    rbv_crl1_config = Cpt(EpicsSignalRO, "filterConfig_RBV", kind="hinted")
    inMask_config = Cpt(EpicsSignalRO, "inMask_RBV", kind="config")
    outMask_config = Cpt(EpicsSignalRO, "outMask_RBV", kind="config")
