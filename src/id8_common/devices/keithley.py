"""
Keithley 2612 readout and control
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal

class Keithley2400(Device):
    """Keithley 2400 SourceMeter -- one channel that both sources and measures.

    Used to bias a sample during an acquisition (see plans/set/volt_seq.py).
    Its configs/devices.yml entry -- currently commented out -- uses prefix
    "8idKeithley2400:K1:". Nothing reaches the sample until `output` is set to 1.
    """

    output = Component(EpicsSignal, "enableBO")

    # --- source: what the instrument drives onto the sample ---
    set_volt = Component(EpicsSignal, "setVoltAO")
    set_curr = Component(EpicsSignal, "setCurrAO")

    set_compl_volt = Component(EpicsSignal, "setComplVoltAO")
    set_compl_curr = Component(EpicsSignal, "setComplCurrAO")

    # --- measure: what the instrument reads back ---
    meas_volt = Component(EpicsSignalRO, 'measVoltAI')
    meas_curr = Component(EpicsSignalRO, 'measCurrAI')
    output_status = Component(EpicsSignalRO, "enabledBI")

class Keithley(Device):
    """One channel (SMU A or B) of a Keithley 2600-series SourceMeter.

    Its configs/devices.yml entries -- currently commented out -- instantiate
    one per channel, prefixes "8idKeithley2600:SMU:A:" and ":B:". Component
    names below mirror the EPICS
    record names one for one: the ``_AO``/``_BO`` suffixes are the settable
    records and ``_AI``/``_BI`` the matching readbacks.
    """

    output = Component(EpicsSignal, "SourceOutputBO")
    
    # --- source: what the instrument drives onto the sample ---
    
    SrcAutorangeV_BO = Component(EpicsSignal, "SrcAutorangeV_BO")
    SrcLowrangeV_AO = Component(EpicsSignal, "SrcLowrangeV_AO")
    SrcRangeV_AO = Component(EpicsSignal, "SrcRangeV_AO")
    SrcLevelV_AO = Component(EpicsSignal, "SrcLevelV_AO")
    SrcLimitV_AO = Component(EpicsSignal, "SrcLimitV_AO")

    SrcAutorangeV_BI = Component(EpicsSignalRO, "SrcAutorangeV_BI")
    SrcLowrangeV_AI = Component(EpicsSignalRO, "SrcLowrangeV_AI")
    SrcRangeV_AI = Component(EpicsSignalRO, "SrcRangeV_AI")
    SrcLevelV_AI = Component(EpicsSignalRO, "SrcLevelV_AI")
    SrcLimitV_AI = Component(EpicsSignalRO, "SrcLimitV_AI")

    SrcAutorangeI_BO = Component(EpicsSignal, "SrcAutorangeI_BO")
    SrcLowrangeI_AO = Component(EpicsSignal, "SrcLowrangeI_AO")
    SrcRangeI_AO = Component(EpicsSignal, "SrcRangeI_AO")
    SrcLevelI_AO = Component(EpicsSignal, "SrcLevelI_AO")
    SrcLimitI_AO = Component(EpicsSignal, "SrcLimitI_AO")

    SrcAutorangeI_BI = Component(EpicsSignalRO, "SrcAutorangeI_BI")
    SrcLowrangeI_AI = Component(EpicsSignalRO, "SrcLowrangeI_AI")
    SrcRangeI_AI = Component(EpicsSignalRO, "SrcRangeI_AI")
    SrcLevelI_AI = Component(EpicsSignalRO, "SrcLevelI_AI")
    SrcLimitI_AI = Component(EpicsSignalRO, "SrcLimitI_AI")

    SettlingDelayMO = Component(EpicsSignal, "SettlingDelayMO")
    OfflimitI_AO = Component(EpicsSignal, "OfflimitI_AO")
    OffmodeMO = Component(EpicsSignal, "OffmodeMO")

    # --- measure: what the instrument reads back ---

    MeasAutorangeI_BO = Component(EpicsSignal, "MeasAutorangeI_BO")
    MeasLowrangeI_AO = Component(EpicsSignal, "MeasLowrangeI_AO")
    MeasRangeI_AO = Component(EpicsSignal, "MeasRangeI_AO")
    MeasRelEnableI_BO = Component(EpicsSignal, "MeasRelEnableI_BO")
    MeasRelLevelI_AO = Component(EpicsSignal, "MeasRelLevelI_AO")

    MeasAutorangeV_BO = Component(EpicsSignal, "MeasAutorangeV_BO")
    MeasLowrangeV_AO = Component(EpicsSignal, "MeasLowrangeV_AO")
    MeasRangeV_AO = Component(EpicsSignal, "MeasRangeV_AO")
    MeasRelEnableV_BO = Component(EpicsSignal, "MeasRelEnableV_BO")
    MeasRelLevelV_AO = Component(EpicsSignal, "MeasRelLevelV_AO")

    MeasAutorangeI_BO = Component(EpicsSignalRO, "MeasAutorangeI_BO")
    MeasLowrangeI_AO = Component(EpicsSignalRO, "MeasLowrangeI_AO")
    MeasRangeI_AO = Component(EpicsSignalRO, "MeasRangeI_AO")
    MeasRelEnableI_BO = Component(EpicsSignalRO, "MeasRelEnableI_BO")
    MeasRelLevelI_AO = Component(EpicsSignalRO, "MeasRelLevelI_AO")

    MeasAutorangeV_BI = Component(EpicsSignalRO, "MeasAutorangeV_BI")
    MeasLowrangeV_AI = Component(EpicsSignalRO, "MeasLowrangeV_AI")
    MeasRangeV_AI = Component(EpicsSignalRO, "MeasRangeV_AI")
    MeasRelEnableV_BI = Component(EpicsSignalRO, "MeasRelEnableV_BI")
    MeasRelLevelV_AI = Component(EpicsSignalRO, "MeasRelLevelV_AI")



