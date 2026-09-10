"""
Lakeshore 336 (temperature readout and control)
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal


class Lakeshore(Device):
    """Device representing a Lakeshore 336 temperature controller.

    This device provides read-only access to temperature measurements from up to four
    input channels of a Lakeshore 336 temperature controller.
    """

    readback_ch1 = Component(EpicsSignalRO, "IN1")
    readback_ch2 = Component(EpicsSignalRO, "IN2")
    readback_ch3 = Component(EpicsSignalRO, "IN3")
    readback_ch4 = Component(EpicsSignalRO, "IN4")

    setpoint_out1 = Component(EpicsSignal, "OUT1:SP")
    setpoint_out2 = Component(EpicsSignal, "OUT2:SP")

def set_temp_lakeshore(temp, wait):

### Set lakeshore2 setpoint in loop 1 and wait for (wait) seconds. 
###

    lakeshore2.setpoint_out1.put(temp)
    time.sleep(wait)