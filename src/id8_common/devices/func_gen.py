"""
Function generator at 8-ID-E, first used on 12/04/2024
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal

# epics_put("dpKeysight:KEY1:1:FUNC",func)
# epics_put("dpKeysight:KEY1:1:BURST:NCYCLES",ncycle) ##number of waves
# epics_put("dpKeysight:KEY1:1:FREQ", freq) ##freq in Hz
# epics_put("dpKeysight:KEY1:1:AMP", amp) #peak to peak amplitude in Volts
# epics_put("dpKeysight:KEY1:TRIG.PROC",send_trigger) #send the waves


class Function_Generator(Device):
    """A device class for controlling function generators in the beamline.

    This class provides control over function generators used for signal generation
    and control. It includes functionality for setting frequency, amplitude, and
    waveform parameters.
    """

    func = Component(EpicsSignal, "1:FUNC")
    ncycle = Component(EpicsSignal, "1:BURST:NCYCLES")
    freq = Component(EpicsSignal, "1:FREQ")
    amp = Component(EpicsSignal, "1:AMP")
    send_trigger = Component(EpicsSignal, "TRIG.PROC")
