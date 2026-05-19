"""Module for controlling Sydor TP4U detectors in the beamline.

This module provides classes and utilities for interfacing with Sydor TP4U
detectors, including configuration, data acquisition, and monitoring
capabilities.
"""

from ophyd import Component
from ophyd import EpicsSignalRO
from ophyd import QuadEM
from ophyd.areadetector.plugins import ImagePlugin_V34
from ophyd.areadetector.plugins import StatsPlugin_V34
from ophyd.quadem import QuadEMPort


class SydorTP4U(QuadEM):
    """
    Sydor TP4U quad electrometer.

    Example::

        xbpm2 = SydorTP4U("8idiSoft:T4U_BPM:", name="xbpm2")
    """

    conf = Component(QuadEMPort, port_name="T4U_BPM")

    # latest versions of the AD plugins
    current1 = Component(StatsPlugin_V34, "Current1:")
    current2 = Component(StatsPlugin_V34, "Current2:")
    current3 = Component(StatsPlugin_V34, "Current3:")
    current4 = Component(StatsPlugin_V34, "Current4:")
    image = Component(ImagePlugin_V34, "image1:")
    sum_all = Component(StatsPlugin_V34, "SumAll:")

    # better as text than numbers
    firmware = Component(EpicsSignalRO, "Firmware", string=True, kind="config")
    model = Component(EpicsSignalRO, "Model", string=True, kind="config")
