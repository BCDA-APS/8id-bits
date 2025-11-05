"""Write-only EpicsSignal variant.

Workaround a problem when random failures occur when writing repeatedly
to an EpicsSignal.  During a series of identical data acquisition events,
RuntimeError exceptions have been raised when writing to an EpicsSignal
object, despite many successful previous writes in the same series.

Common to all these failures, is a failed status object where the target
value is not received.  (done=True, success=False)

Create a "write-only" EpicsSignal subclass that restores the simpler
method of EpicsSignalBase.set().
"""

from ophyd import EpicsSignal
from ophyd.signal import DEFAULT_WRITE_TIMEOUT
from ophyd.status import Status


class EpicsSignalWO(EpicsSignal):
    """Write-only EpicsSignal variant.

    Here, the .set() method is always "successful".

    Bypass put_complete and readback verification.
    Hard-codes the status(done-True, success=True).
    """

    def set(self, value, *, timeout=DEFAULT_WRITE_TIMEOUT, settle_time=None):
        """Must replace set() method to avoid default status object."""
        if timeout is DEFAULT_WRITE_TIMEOUT:
            timeout = self.write_timeout

        # bypass put_complete and readback verification
        self.put(value)  # <== changed
        st = Status(self, timeout=timeout, settle_time=settle_time)
        st.set_finished()  # <== changed
        st.wait()  # <== added
        return st
