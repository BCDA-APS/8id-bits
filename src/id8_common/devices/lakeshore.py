"""
Lakeshore 336 (temperature readout and control)
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignalRO, EpicsSignal


class Lakeshore(Device):
    """One Lakeshore 336 temperature controller.

    readback_ch1..4 are the four sensor inputs (read-only). setpoint_out1 and
    setpoint_out2 are the two control-loop setpoints, and are writable.

    Built from configs/devices.yml as `lakeshore1` and `lakeshore2`, prefixes
    "8ideSoft:LS336:1:" and ":2:".
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