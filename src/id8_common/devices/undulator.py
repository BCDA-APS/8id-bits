"""
8ID undulator readout.
"""

from apstools.devices.aps_undulator import Revolver_Undulator


class RevolverUndulator_8ID(Revolver_Undulator):
    """One of the two 8-ID revolver undulators (upstream S08ID:USID:, downstream :DSID:).

    apstools' Revolver_Undulator with the one Component the 8-ID IOCs do not
    serve removed -- see below.
    """

    # Set these Components to None to avoid missing PV error
    version_hdmu = None
