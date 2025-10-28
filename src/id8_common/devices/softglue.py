"""
The BCDA SoftGlue FPGA generates pulses to trigger an area detector.

Example acquisition sequence::

    fpga = SoftGlue("", name="fpga")
    ...
    # configure fpga & area_detector
    ...
    yield from bps.mv(fpga.start_pulses, "1!")
    yield from bp.count([area_detector])
    yield from bps.mv(fpga.stop_pulses, "1!")

"""

from ophyd import Device
from ophyd import EpicsSignal
from ophyd import FormattedComponent as FCpt


class SoftGlue(Device):
    """
    BCDA SoftGlue FPGA to provide external trigger signal for detectors.

    The value to write to start & stop is the same text: 1!
    """

    def __init__(
        self,
        prefix: str,
        *args,
        pv_acq_period: str = "",
        pv_acq_time: str = "",
        pv_num_triggers: str = "",
        pv_start_pulses: str = "",
        pv_stop_pulses: str = "",
        **kwargs,
    ):
        """Initialize the SoftGlue device.

        Args:
            prefix: The EPICS PV prefix for the device
            pv_acq_period: The PV name for acquisition period
            pv_acq_time: The PV name for acquisition time
            pv_num_triggers: The PV name for number of triggers
            pv_start_pulses: The PV name for start pulses
            pv_stop_pulses: The PV name for stop pulses
            **kwargs: Additional keyword arguments passed to the Device constructor
        """
        self._pv_acq_period = pv_acq_period
        self._pv_acq_time = pv_acq_time
        self._pv_num_triggers = pv_num_triggers
        self._pv_start_pulses = pv_start_pulses
        self._pv_stop_pulses = pv_stop_pulses
        super().__init__(prefix, *args, **kwargs)

    acq_period = FCpt(EpicsSignal, "{prefix}{_pv_acq_period}", kind="config")
    acq_time = FCpt(EpicsSignal, "{prefix}{_pv_acq_time}", kind="config")
    num_triggers = FCpt(EpicsSignal, "{prefix}{_pv_num_triggers}", kind="config")
    start_pulses = FCpt(
        EpicsSignal,
        "{prefix}{_pv_start_pulses}",
        kind="omitted",
        string=True,
    )
    stop_pulses = FCpt(
        EpicsSignal,
        "{prefix}{_pv_stop_pulses}",
        kind="omitted",
        string=True,
    )
