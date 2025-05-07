
from ophyd import EpicsMotor
from ophyd import FormattedComponent as FCpt
from ophyd import MotorBundle


class HV_Motors(MotorBundle):
    """
    Horizontal/Vertical Motor Stage with h & v translations.

    Parameters

    prefix
        (str) : EPICS PV prefix common to both motors.
        Includes trailing delimiter, such as ':', as needed.
        (Same convention as used with GUI screens, IOC startup scripts,
        and other EPICS configuration.)
    pv_h
        (str) : EPICS PV suffix for the horizontal motor PV.
    pv_v
        (str) : EPICS PV suffix for the vertical motor PV.

    Example::

        base = HV_Motors("IOC:DB:", name="base", pv_h="m9", pv_v="m10")
        # which creates:
        # base.h.pvname: IOC:DB:m9
        # base.v.pvname: IOC:DB:m10
    """

    h = FCpt(EpicsMotor, "{prefix}{_pv_h}", labels={"motors"})
    v = FCpt(EpicsMotor, "{prefix}{_pv_v}", labels={"motors"})

    def __init__(
        self,
        prefix: str,
        *args,
        pv_h: str = "",
        pv_v: str = "",
        **kwargs,
    ):
        self._pv_h = pv_h
        self._pv_v = pv_v
        super().__init__(prefix, *args, **kwargs)
