"""
Keithley 2612 readout and control
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal

class Keithley(Device):
    """Device representing a Keithley 2600 class device.
    Component(EpicsSignal, "")
    """

    output = Component(EpicsSignal, "SourceOutputBO")
    
    """Source"""
    
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

    """Measure"""

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



