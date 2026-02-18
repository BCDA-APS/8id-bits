"""
Function generator at 8-ID-E, first used on 12/04/2024
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal, EpicsSignalRO

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
    """Control PVs"""
    func = Component(EpicsSignal, "1:FUNC")
    frequency = Component(EpicsSignal, "1:FREQ")
    amplitude = Component(EpicsSignal, "1:AMP")
    phase = Component(EpicsSignal, "1:PHASE")
    pulse_width = Component(EpicsSignal, "1:PULSEW")
    trigger_source = Component(EpicsSignal, "1:TRIG:SOURCE")
    trigger_edge = Component(EpicsSignal, "1:TRIG:SLOPE")
    burst_count = Component(EpicsSignal, "1:BURST:NCYCLES")
    burst_mode = Component(EpicsSignal, "1:BURST:MODE")
    burst_state = Component(EpicsSignal, "1:BURST:STATE")
    send_trigger = Component(EpicsSignal, "TRIG.PROC")
    output = Component(EpicsSignal, "1:OUT")

    """Readback PVs"""
    func_rbv = Component(EpicsSignalRO, "1:FUNC:RBV")
    frequency_rbv = Component(EpicsSignalRO, "1:FREQ:RBV")
    amplitude_rbv = Component(EpicsSignalRO, "1:AMP:RBV")
    phase_rbv = Component(EpicsSignalRO, "1:PHASE:RBV")
    pulse_width_rbv = Component(EpicsSignalRO, "1:PULSEW:RBV")
    trigger_source_rbv = Component(EpicsSignalRO, "1:TRIG:SOURCE:RBV")
    trigger_edge_rbv = Component(EpicsSignalRO, "1:TRIG:SLOPE:RBV")
    burst_count_rbv = Component(EpicsSignalRO, "1:BURST:NCYCLES:RBV")
    burst_mode_rbv = Component(EpicsSignalRO, "1:BURST:MODE:RBV")
    burst_state_rbv = Component(EpicsSignalRO, "1:BURST:STATE:RBV")
    output_rbv = Component(EpicsSignalRO, "1:OUT:RBV")

