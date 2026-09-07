"""Safe device loading for the Bluesky (guarneri) startup path.

Thin adapter over ``id8_common.registry.safe_make_devices`` -- see that
function for how a device that fails to build, or fails its basic
``wait_for_connection()`` check, gets skipped instead of taking down the
rest of startup.
"""

import sys

from id8_common.registry import safe_make_devices as _safe_make_devices

MAIN_NAMESPACE = "__main__"


def safe_make_devices(file: str, device_manager, clear: bool = False, timeout: float = 20.0):
    if clear:
        device_manager.devices.clear()

    def _register(device):
        # Match apsbits.core.instrument_init.guarneri_namespace_loader: a
        # device is available as a bare name at the interactive prompt
        # (e.g. `sample_x.move(5)`), not just via oregistry["sample_x"].
        # Lost otherwise, since this bypasses apsbits' own loader entirely.
        device_manager.devices.register(device)
        setattr(sys.modules[MAIN_NAMESPACE], device.name, device)

    return _safe_make_devices(file, _register, timeout=timeout)
