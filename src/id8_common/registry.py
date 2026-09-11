"""Hand-written device registry and YAML device loader.

Replaces ``apsbits.core.instrument_init`` (``oregistry`` / ``make_devices``)
so plan code has no Bluesky/apsbits dependency. Both ``startup.py``
(Bluesky) and ``startup_ophyd.py`` (Ophyd-only) populate this same shared
registry after they load devices, so any module that does
``from id8_common.registry import oregistry`` and looks devices up by
``oregistry["name"]`` works unchanged under either startup script. See
``docs/devices.md`` for how devices are loaded and looked up.
"""

import importlib
import inspect
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

CONFIGS_DIR = Path(__file__).parent / "configs"

# Shared, bounded across every safe_make_devices() call in a session, so
# devices.yml/ad_devices.yml/devices_aps_only.yml checking each other's
# concurrent EPICS Channel Access traffic in startup.py stays throttled
# instead of each file opening its own unbounded burst of connections.
_CONNECTIVITY_POOL = ThreadPoolExecutor(max_workers=4)


class DeviceRegistry:
    """Dict-backed device registry: name -> ophyd device."""

    def __init__(self):
        self._devices = {}

    def register(self, device):
        self._devices[device.name] = device
        return device

    def __getitem__(self, name):
        return self._devices[name]

    def get(self, name, default=None):
        """Dict-style lookup that returns `default` instead of raising.

        For module-level `x = oregistry.get("x")` in plan files: a device
        skipped by safe_make_devices() (offline IOC, disconnected PV) then
        resolves to None at import time instead of raising KeyError and
        aborting the whole startup script. The plan still fails if that
        device is actually used while None -- just later, at the call site
        that needs it, instead of at import time for every other plan too.
        """
        return self._devices.get(name, default)

    def __contains__(self, name):
        return name in self._devices

    def __iter__(self):
        return iter(self._devices.values())

    def __len__(self):
        return len(self._devices)

    @property
    def device_names(self):
        return list(self._devices)

    def findall(self, name=None):
        """Every registered device, or just the one(s) with this name.

        Names are unique here (they are the dict keys), so filtering by name
        gives either a one-item list or an empty one. Returning a list is the
        point: ``utils/nexus_utils.py`` calls ``findall(name="aps")`` to test
        whether the storage-ring device loaded, without having to catch a
        KeyError. The name and signature come from the ophyd-registry API this
        class stands in for.
        """
        if name is None:
            return list(self._devices.values())

        matching = []
        for device in self._devices.values():
            if device.name == name:
                matching.append(device)
        return matching

    def clear(self):
        self._devices.clear()


oregistry = DeviceRegistry()


def get_connected_device(name, registry=oregistry):
    """Resolve a device by name, raising if it is missing or disconnected.

    Use this inside a plan function instead of caching the lookup at module
    scope (``x = oregistry.get("x")``). A device that ``safe_make_devices``
    skipped is absent from the registry entirely, so the module-level form
    binds ``None`` at import time and the real failure surfaces much later,
    as an ``AttributeError`` on ``None`` somewhere unrelated to the cause.

    ``softglue`` is deliberately still resolved at module scope
    (``plans/set/shutter_att.py``, ``plans/align/ophyd_scan.py``) -- it is
    present in every session, so a per-call check buys nothing.
    ``pv_registers`` was in that category too, but nothing binds it at module
    scope any more: ``expt_config.py`` resolves it through this function,
    lazily, on the first read or write of ``expt.measurement_num``.

    Raises
    ------
    KeyError
        Not in the registry: it failed to build or connect at startup (see the
        "*** Devices not online: ... ***" banner), or the name is wrong.
        ``oregistry.device_names`` lists what did load.
    RuntimeError
        Registered at startup, but not connected now -- e.g. its IOC has since
        gone away.
    """
    try:
        device = registry[name]
    except KeyError:
        raise KeyError(
            f"{name!r} is not in the device registry -- skipped at startup "
            f"(offline IOC or bad prefix), or the name is misspelled. "
            f"See oregistry.device_names for what did load."
        ) from None

    if not device.connected:
        raise RuntimeError(f"{name!r} is in the registry but is not connected now.")

    return device


def get_ophyd_object(dotted, registry=oregistry):
    """Resolve a dotted device path like ``"huber.x"`` to the ophyd object.

    The dotted extension of :func:`get_connected_device`: the first segment is
    a registry key and must be a connected device; the rest are plain
    attribute lookups on it. ``"softglue"`` and ``"huber.x.user_setpoint"``
    are both valid.

    Used wherever a YAML file names a motor or a signal as a string --
    sample_info.yaml's ``inner_motor``, a trio protocol's ``motors:`` block,
    a geometry field given as a live PV rather than a number.

    Deliberately strict: the root must be CONNECTED, not merely registered.
    Callers use this to prove a motor named in YAML is usable *before* anything
    moves, so resolving to a device whose IOC is down would defeat the check.
    Earlier copies of this helper did a bare ``oregistry[name]`` and skipped
    that test; if you are tempted to reintroduce a lax variant, don't -- add a
    keyword here instead.

    This is the only definition. ``acq_helpers`` re-exports it, and both
    ``master_plan`` and ``trio_master_plan_rigaku3m_eiger4m_lambda2m`` reach it through
    the same ``validators.validate_sample_motion``, so the serial and trio
    paths cannot drift apart again.

    One holdout, if you are auditing: ``plans/set/select_device._resolve``
    walks dotted paths from device_position.yaml its own way.

    Raises
    ------
    KeyError, RuntimeError
        From ``get_connected_device`` -- see there.
    AttributeError
        A segment after the first does not exist on the object so far.
    """
    parts = str(dotted).split(".")
    obj = get_connected_device(parts[0], registry=registry)

    for part in parts[1:]:
        try:
            obj = getattr(obj, part)
        except AttributeError:
            raise AttributeError(
                f"{dotted!r}: {type(obj).__name__} has no attribute {part!r}."
            ) from None

    return obj


def _load_yaml_spec(file, path=None):
    """Read one device YAML file: `path` if given, otherwise configs/<file>.

    An empty file parses as None, so `or {}` turns it into an empty spec that
    the caller can loop over instead of crashing.
    """
    if path is not None:
        yaml_path = Path(path)
    else:
        yaml_path = CONFIGS_DIR / file

    return yaml.safe_load(yaml_path.read_text()) or {}


def make_devices(file, registry=oregistry, path=None, clear=False, timeout=20.0):
    """Instantiate devices from a Guarneri-style YAML file and register them.

    Each top-level key in the YAML is a dotted class or factory path; each
    list item under it is the kwargs used to build one device --
    ``target(**kwargs)``. This works whether ``target`` is a plain device
    class (``id8_common.devices.xy_motors.XY_Motors``) or a factory function
    (``apstools.devices.area_detector_factory.ad_creator``) -- both are just
    ``callable(**kwargs) -> device``.

    A factory that returns a generator instead of a device (e.g.
    ``apsbits.utils.sim_creator.predefined_device``) is iterated, so one
    YAML entry can register more than one device.

    Safe by default: delegates to ``safe_make_devices`` so a device that
    fails to build or connect is skipped (with a warning) instead of
    raising -- see that function for why. Returns the list of skipped
    device names (previously returned ``registry``; no caller in this repo
    used that return value).
    """
    if clear:
        registry.clear()

    return safe_make_devices(file, registry.register, path=path, timeout=timeout)


def safe_make_devices(file, register, *, path=None, timeout=20.0):
    """Instantiate devices from YAML, skipping any that don't build or connect.

    Same YAML shape as ``make_devices``. Unlike
    ``guarneri.Instrument.make_devices`` (which calls ``Klass(**kwargs)`` in
    a plain loop with no per-entry ``try/except`` -- one bad entry raises
    and starves every device defined *after* it in the same file), each
    instance here is isolated:

    - construct it in its own ``try/except``, so a bad prefix, unknown
      class, or bad kwarg only skips that one device.
    - check connectivity with a plain ``device.wait_for_connection(timeout=
      timeout)`` -- ophyd's own default, which only waits on signals the
      device already eagerly created during ``__init__`` (lazy components
      are left alone entirely, never forced into existence). A device is
      skipped only if that fails, i.e. it genuinely looks offline (IOC
      down, wrong prefix).

    This intentionally does *not* try to prove every lazily-declared PV on
    a device is reachable (e.g. an area-detector plugin's less-common
    fields) -- an earlier version of this function did that by walking and
    force-connecting every lazy component up front
    (``include_lazy=True``), which caught a single dropped PV in an
    otherwise-fine device but also flooded EPICS Channel Access when many
    heavy area-detector devices were checked at once, and produced a
    "connected with N PV(s) not found" log line long enough to be
    unusable in practice. A device with a missing lazy PV now registers
    normally at startup, and that PV's own connection is only attempted
    -- and only then reported as missing -- when an align or acquisition
    plan actually reads or writes it, matching stock ophyd behavior.

    Only devices that connect are handed to ``register``. Returns the list
    of skipped device names, so the caller can report them (e.g. the
    ``offline_devices`` banner in ``startup.py`` / ``startup_ophyd.py``).

    ``register`` is a plain callable so this same implementation works for
    both the Ophyd-only registry (``id8_common.registry.oregistry.register``)
    and the Bluesky guarneri registry (``instrument.devices.register``).

    Connectivity checks run through a small, bounded thread pool
    (``_CONNECTIVITY_POOL``, module-level and shared across every call --
    ``devices.yml``, then ``ad_devices.yml``, then ``devices_aps_only.yml``
    in ``startup.py``) rather than one unbounded pool per call, so that
    many devices' EPICS Channel Access connection attempts don't all fire
    at the same instant.
    """
    spec = _load_yaml_spec(file, path=path)

    candidates = []  # (device, fallback_name) pairs awaiting a connectivity check
    skipped = []

    for class_path, instances in spec.items():
        try:
            module_name, _, attr_name = class_path.rpartition(".")
            module = importlib.import_module(module_name)
            target = getattr(module, attr_name)
        except (ImportError, AttributeError) as exc:
            logger.warning("Skipping class %r -- cannot import: %s", class_path, exc)
            skipped.append(class_path)
            continue

        for kwargs in instances or []:
            fallback_name = kwargs.get("name", f"{class_path}({kwargs})")
            try:
                obj = target(**kwargs)
            except Exception as exc:
                logger.warning("Skipping device %r -- failed to construct: %s", fallback_name, exc)
                skipped.append(fallback_name)
                continue

            devices = obj if inspect.isgenerator(obj) else [obj]
            for device in devices:
                candidates.append((device, fallback_name))

    def _check(device):
        """Raise if the device's already-instantiated signals never connect."""
        if hasattr(device, "wait_for_connection"):
            device.wait_for_connection(timeout=timeout)
        elif hasattr(device, "connected") and not device.connected:
            raise TimeoutError(f"{getattr(device, 'name', device)} not connected")

    # Submit one connectivity check per device, keyed by the future it returns.
    # as_completed() below hands the futures back in whatever order they finish,
    # not in candidate order, so this mapping is the only way to say afterwards
    # which device a failure belonged to.
    futures = {}
    for device, fallback_name in candidates:
        future = _CONNECTIVITY_POOL.submit(_check, device)
        futures[future] = (device, fallback_name)

    for future in as_completed(futures):
        device, fallback_name = futures[future]
        name = getattr(device, "name", fallback_name)
        try:
            future.result()
        except Exception as exc:
            logger.warning("Skipping device %r -- not connected: %s", name, exc)
            skipped.append(name)
            continue
        register(device)

    return skipped
