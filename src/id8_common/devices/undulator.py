"""
8ID undulator readout.
"""

from apstools.devices.aps_undulator import Revolver_Undulator


class RevolverUndulator_8ID(Revolver_Undulator):

    # Set these Components to None to avoid missing PV error
    version_hdmu = None
