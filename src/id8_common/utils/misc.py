"""Utility functions for 8ID beamline operations."""


def ioc_alive(pv: str, timeout: float=0.5, retries: int=2) -> bool:
    """
    Check if an IOC is alive, return True if a PV responds, False otherwise.
    """
    import time
    from epics import caget

    for _ in range(retries):
        try:
            val = caget(pv, timeout=timeout)
            if val is not None:
                return True
        except Exception:
            pass
        time.sleep(0.1)
    return False


def stream_rois(det, stats_nums=(1, 2, 3), fields=("total",), hinted=("total",)):
    """Configure 'kind' attribute for area detector signals."""
    det.read_attrs = [a for a in det.read_attrs if "." not in a]
    for n in stats_nums:
        stats_attr = f"stats{n}"

        if stats_attr not in det.read_attrs:
            det.read_attrs.append(stats_attr)

        stats = getattr(det, stats_attr)

        stats.kind = "hinted"

        stats.read_attrs = list(fields)

        for f in fields:
            sig = getattr(stats, f)
            sig.kind = "hinted" if f in hinted else "normal"
