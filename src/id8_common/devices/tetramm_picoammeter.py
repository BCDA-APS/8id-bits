"""
Caen picoammeter - TetraAMM

GitHub apstools issue #878 has some useful documentation in the comments.

.. see:: https://github.com/BCDA-APS/apstools/issues/878
"""

import logging
import time as ttime

from ophyd import Component
from ophyd import TetrAMM
from ophyd.areadetector.plugins import ImagePlugin_V34
from ophyd.areadetector.plugins import StatsPlugin_V34
from ophyd.device import Staged
from ophyd.quadem import QuadEMPort

logger = logging.getLogger(__name__)
logger.info(__file__)


class MyTetrAMM(TetrAMM):
    """Caen picoammeter - TetraAMM."""

    conf = Component(QuadEMPort, port_name="QUAD_PORT")

    current1 = Component(StatsPlugin_V34, "Current1:")
    current2 = Component(StatsPlugin_V34, "Current2:")
    current3 = Component(StatsPlugin_V34, "Current3:")
    current4 = Component(StatsPlugin_V34, "Current4:")
    image = Component(ImagePlugin_V34, "image1:")
    sum_all = Component(StatsPlugin_V34, "SumAll:")

    def __init__(self, *args, port_name="TetrAMM", **kwargs):
        """custom port name"""
        super().__init__(*args, **kwargs)
        self.conf.port_name.put(port_name)  # fix the port name here
        self.stage_sigs = {}

        # Mark some components as "config" so they do not appear on data rows.
        for attr_name in self.component_names:
            attr = getattr(self, attr_name)
            if attr_name.startswith("current_"):
                for ch_name in attr.component_names:
                    getattr(attr, ch_name).kind = "config"
            elif attr_name.startswith("position_"):
                attr.kind = "config"

        self.sum_all.mean_value.kind = "hinted"  # Show as a data column in SPEC file.
        # tetramm1.sum_all.mean_value.kind = "hinted"
        self.current1.mean_value.kind = "hinted"
        self.current2.mean_value.kind = "hinted"
        self.current3.mean_value.kind = "hinted"
        self.current4.mean_value.kind = "hinted"

    def trigger(self):
        """
        Tetramms are operated in continuous mode. Will just accept the current value.
        """
        if self._staged != Staged.yes:
            raise RuntimeError(
                "This detector is not ready to trigger."
                "Call the stage() method before triggering."
            )

        self._status = self._status_type(self)
        self._acquisition_signal.put(1, wait=False)
        self.generate_datum(self._image_name, ttime.time(), {})
        self._status.set_finished()
        return self._status

