from concurrent.futures import ThreadPoolExecutor, as_completed
import logging, tempfile
from pathlib import Path
import yaml
from apsbits.core.instrument_init import make_devices
from id8_common.utils.misc import ioc_alive

logger = logging.getLogger(__name__)

def safe_make_devices(file: str, device_manager, clear: bool = False, timeout: float = 1.0):
    config_path = Path(__file__).parent.parent / "configs" / file
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Collect all (device_class, idx, pv) to check in parallel
    checks = []
    for cls, entries in config.items():
        if isinstance(entries, list):
            for i, entry in enumerate(entries):
                pv = entry.get("heartbeat_pv")
                if pv:
                    checks.append((cls, i, pv))

    # Run all PV checks concurrently
    alive = {}
    with ThreadPoolExecutor() as pool:
        futures = {pool.submit(ioc_alive, pv, timeout): (cls, i, pv) for cls, i, pv in checks}
        for fut in as_completed(futures):
            cls, i, pv = futures[fut]
            alive[(cls, i)] = fut.result()

    # Build filtered config, stripping heartbeat_pv
    filtered = {}
    skipped, loaded = [], []
    for cls, entries in config.items():
        if not isinstance(entries, list):
            filtered[cls] = entries
            continue
        kept = []
        for i, entry in enumerate(entries):
            entry = dict(entry)  # copy so we don't mutate original
            entry.pop("heartbeat_pv", None)  # strip before passing to ophyd
            name = entry.get("name", f"{cls}[{i}]")
            if (cls, i) in alive:
                if alive[(cls, i)]:
                    kept.append(entry)
                    loaded.append(name)
                else:
                    skipped.append(name)
                    logger.warning("Skipping %r — IOC not responding", name)
            else:
                kept.append(entry)  # no heartbeat_pv = always load
        if kept:
            filtered[cls] = kept

    if loaded:
        logger.info("Loaded devices after IOC check: %s", loaded)

    # Write temp filtered YAML, call make_devices, clean up
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yml", dir=config_path.parent,
        delete=False, prefix="_hb_filtered_"
    ) as tf:
        yaml.dump(filtered, tf)
        temp_name = Path(tf.name).name

    try:
        make_devices(clear=clear, file=temp_name, device_manager=device_manager)
    finally:
        (config_path.parent / temp_name).unlink(missing_ok=True)

    return skipped
