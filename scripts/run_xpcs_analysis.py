#!/usr/bin/env python3
"""Batch XPCS correlation over a range of 8-ID dataset folders.

Wraps ``boost_corr_bin``.  Everything is configured from a YAML file; the only
command-line argument is the path to that file::

    ./run_xpcs_analysis.py                      # reads ./xpcs_analysis.yaml
    ./run_xpcs_analysis.py my_config.yaml

Run it on a machine that can see both the data (``/gdata/...``) and
``boost_corr_bin`` -- e.g. a GPU node such as adamite.  The script itself does
no ssh and makes no assumption about the host.

Dataset folders are named ``<header>_<sample>_a####_f######_<detector>_<repeat>``,
for example ``A0056_HEA-0GPa_a0002_f003000_eiger4M_r00001``.  The header index
(``A0056``) selects the measurement, the repeat index (``r00001``) counts
repeats of the same condition.  The header letter defaults to ``A`` and is
changed with the ``header`` key.

For every selected folder the script runs one correlation per peak:

* ``rigaku3M`` -- a single run, no peak in the result name.
* ``eiger4M``  -- two runs, one per peak (``BCC`` and ``FCC``), passed to
  ``boost_corr`` as ``-u BCC`` / ``-u FCC`` so the results land side by side as
  ``..._BCC_results.hdf`` and ``..._FCC_results.hdf``.

Results are overwritten in place, so re-running a range replaces the previous
``..._BCC_results.hdf`` rather than piling up ``..._BCC_results_01.hdf``.  Set
``overwrite: false`` to keep the old files and auto-number instead.

Qmap selection, highest priority first:

1. A *dedicated* qmap sitting in the data directory, named strictly
   ``<detector>_qmap[_<peak>]_<header>_<repeat>.hdf``, e.g.
   ``eiger4m_qmap_BCC_A0057_r00001.hdf`` or ``rigaku3m_qmap_A0057_r00001.hdf``.
   Matched case-insensitively.  Set ``ignore_dedicated: true`` to bypass.
2. The default for that detector/peak from the ``qmaps`` block of the config.
3. The built-in ``DEFAULT_QMAPS`` names below.

Relative paths in the config are resolved against the directory holding the
config file, so a config in ``analysis/`` can say ``data_dir: ../data``.

Run ``./run_xpcs_analysis.py --write-example`` to print a fully commented
starter config.
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from typing import Dict
from typing import List
from typing import NamedTuple
from typing import Optional
from typing import Sequence
from typing import Set
from typing import Tuple

try:
    import yaml
except ImportError:  # pragma: no cover - depends on the interpreter used
    yaml = None

BOOST_CORR_BIN = "boost_corr_bin"

DEFAULT_CONFIG_NAME = "xpcs_analysis.yaml"

DEFAULT_HEADER_LETTER = "A"

# Peaks analysed per detector.  ``None`` means a single run with no peak token.
DETECTOR_PEAKS: Dict[str, Tuple[Optional[str], ...]] = {
    "eiger4M": ("BCC", "FCC"),
    "rigaku3M": (None,),
    "lambda2M": (None,),
}

# Fallback qmap names, looked up in the data directory.  Same naming scheme as
# the dedicated qmaps, with "default" in place of "<header>_<repeat>".
DEFAULT_QMAPS: Dict[Tuple[str, Optional[str]], str] = {
    ("eiger4M", "BCC"): "eiger4m_qmap_BCC_default.hdf",
    ("eiger4M", "FCC"): "eiger4m_qmap_FCC_default.hdf",
    ("rigaku3M", None): "rigaku3m_qmap_default.hdf",
    ("lambda2M", None): "lambda2m_qmap_default.hdf",
}

QMAP_EXTENSIONS = (".hdf", ".h5")

# Tried in order inside a dataset folder; first hit wins.
RAW_EXTENSIONS = (".h5", ".hdf", ".bin.000", ".bin", ".imm")

ANALYSIS_TYPES = ("Multitau", "Twotime", "Both")

# Every key the config file may contain, with its default.
CONFIG_DEFAULTS: Dict[str, Any] = {
    "range": None,  # required
    "header": DEFAULT_HEADER_LETTER,
    "data_dir": None,  # cwd, or ../data when run from an analysis directory
    "output_dir": None,  # <data_dir>/../analysis/<type>
    "detectors": ["eiger4M", "rigaku3M"],
    "peaks": None,  # all peaks known for each detector
    "type": "Both",
    "gpu_id": -2,
    "overwrite": True,
    "ignore_dedicated": False,
    "qmaps": {},
    "extra_args": [],
    "boost_corr": BOOST_CORR_BIN,
    "dry_run": False,
}

EXAMPLE_CONFIG = """\
# Configuration for run_xpcs_analysis.py
#
# Relative paths are resolved against the directory holding THIS file.

# ---------------------------------------------------------------- what to run
# Header index range.  Accepts "56", "53-57", "53-57,60,62-64", or a list.
# A leading header letter is tolerated, so "A0053-A0057" works too.
range: 53-57

# Letter in front of the index in folder names (A0056_... -> "A").
header: A

# Where the dataset folders and qmap files live.
# Omit to use the current directory, or ../data when run from analysis/.
data_dir: ../data

# Where results are written.  Omit for <data_dir>/../analysis/<type>.
# output_dir: ../analysis/Both

# Detectors to process.  eiger4M is run once per peak, rigaku3M once total.
detectors:
  - eiger4M
  - rigaku3M

# Restrict which peaks run.  Omit for all peaks of each detector (BCC + FCC).
# peaks:
#   - BCC

# ------------------------------------------------------------------ boost_corr
type: Both          # Multitau, Twotime or Both
gpu_id: -2          # -1 CPU, -2 auto-schedule, >=0 a specific GPU
overwrite: true     # false keeps old results and auto-numbers (_01, _02, ...)

# Extra flags appended verbatim to every boost_corr command.
extra_args: []
# extra_args: ["-b", "100", "-v"]

# Executable; an absolute path works if it is not on PATH.
boost_corr: boost_corr_bin

# ----------------------------------------------------------------------- qmaps
# A per-dataset qmap named <detector>_qmap[_<peak>]_<header>_<repeat>.hdf in
# data_dir always wins over the defaults below, e.g.
#   eiger4m_qmap_BCC_A0057_r00001.hdf
#   rigaku3m_qmap_A0057_r00001.hdf
# Set ignore_dedicated: true to always use the defaults instead.
ignore_dedicated: false

# Default qmaps, used when no dedicated qmap exists for a dataset.
# Bare names are looked up in data_dir; paths are used as given.
qmaps:
  rigaku3M: rigaku3m_qmap_default.hdf
  eiger4M:
    BCC: eiger4m_qmap_BCC_default.hdf
    FCC: eiger4m_qmap_FCC_default.hdf

# ----------------------------------------------------------------------- other
# true prints the commands and exits without running anything.
dry_run: false
"""


def folder_pattern(header_letter: str) -> "re.Pattern[str]":
    """Regex matching ``<letter>####_<middle>_r#####`` folder names."""
    return re.compile(
        r"^(?P<header>{letter}\d+)_(?P<middle>.+)_(?P<repeat>r\d+)$".format(letter=re.escape(header_letter)),
        re.IGNORECASE,
    )


class ConfigError(Exception):
    """The config file is missing something or says something impossible."""


class Dataset(NamedTuple):
    """One acquisition folder on disk."""

    folder: Path
    header: str  # e.g. "A0056"
    index: int  # e.g. 56
    repeat: str  # e.g. "r00001"
    detector: str  # canonical spelling, e.g. "eiger4M"


class Job(NamedTuple):
    """One boost_corr invocation."""

    dataset: Dataset
    peak: Optional[str]
    raw: Path
    qmap: Path
    qmap_origin: str  # "dedicated" or "default"

    @property
    def label(self) -> str:
        tail = f" {self.peak}" if self.peak else ""
        return f"{self.dataset.folder.name}{tail}"


def parse_index_range(spec: Any, header_letter: str) -> Set[int]:
    """Expand ``"53-57,60,62-64"`` into ``{53..57, 60, 62..64}``.

    Accepts a string, a bare integer, or a list of either.  A leading header
    letter is tolerated, so ``A0053-A0057`` works too.
    """
    if isinstance(spec, (list, tuple)):
        text = ",".join(str(item) for item in spec)
    else:
        text = str(spec)

    strip = header_letter + header_letter.swapcase()
    indices: Set[int] = set()
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        bounds = chunk.split("-")
        if len(bounds) > 2:
            raise ConfigError(f"cannot parse index range {chunk!r}")
        try:
            numbers = [int(b.strip().lstrip(strip)) for b in bounds]
        except ValueError:
            raise ConfigError(f"cannot parse index range {chunk!r}") from None
        first, last = numbers[0], numbers[-1]
        if last < first:
            raise ConfigError(f"index range {chunk!r} runs backwards")
        indices.update(range(first, last + 1))
    if not indices:
        raise ConfigError(f"index range {spec!r} selects nothing")
    return indices


def default_data_dir() -> Path:
    """Guess the data directory: cwd, or its ``../data`` sibling from ``analysis/``."""
    cwd = Path.cwd()
    if cwd.name.lower() == "analysis" and (cwd.parent / "data").is_dir():
        return cwd.parent / "data"
    return cwd


def resolve_path(value: Any, base: Path) -> Path:
    """Resolve a config path against the config file's directory."""
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def as_string_list(value: Any, key: str) -> List[str]:
    """Accept either a YAML list or a comma-separated string."""
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    if isinstance(value, (list, tuple)):
        return [str(item).strip() for item in value if str(item).strip()]
    raise ConfigError(f"{key!r} must be a list or a comma-separated string, got {value!r}")


def build_default_qmaps(raw: Any, detectors: Sequence[str]) -> Dict[Tuple[str, Optional[str]], str]:
    """Merge the config ``qmaps`` block over the built-in defaults."""
    defaults = dict(DEFAULT_QMAPS)
    if raw is None:
        return defaults
    if not isinstance(raw, dict):
        raise ConfigError(f"'qmaps' must be a mapping of detector -> qmap, got {raw!r}")

    canonical = {name.lower(): name for name in DETECTOR_PEAKS}
    for detector, value in raw.items():
        key = canonical.get(str(detector).lower())
        if key is None:
            known = ", ".join(DETECTOR_PEAKS)
            raise ConfigError(f"unknown detector {detector!r} in 'qmaps'; known: {known}")
        peaks = DETECTOR_PEAKS[key]
        if isinstance(value, dict):
            for peak, name in value.items():
                if str(peak) not in [p for p in peaks if p]:
                    listed = ", ".join(p for p in peaks if p) or "none"
                    raise ConfigError(f"unknown peak {peak!r} for {key} in 'qmaps'; known: {listed}")
                defaults[(key, str(peak))] = str(name)
        else:
            # A single name applies to every peak of that detector.
            for peak in peaks:
                defaults[(key, peak)] = str(value)
    return defaults


def load_config(path: Path) -> SimpleNamespace:
    """Read, validate and normalise the YAML config."""
    if yaml is None:
        raise ConfigError(
            "PyYAML is not available to this interpreter. Install it "
            "(pip install pyyaml) or run the script with a python that has it."
        )
    if not path.is_file():
        raise ConfigError(f"config file not found: {path}\nWrite a starter one with --write-example")

    with open(path, "r") as handle:
        loaded = yaml.safe_load(handle)
    if loaded is None:
        loaded = {}
    if not isinstance(loaded, dict):
        raise ConfigError(f"config must be a mapping of key: value, got {type(loaded).__name__}")

    unknown = sorted(set(loaded) - set(CONFIG_DEFAULTS))
    if unknown:
        known = ", ".join(sorted(CONFIG_DEFAULTS))
        raise ConfigError(f"unknown config key(s): {', '.join(unknown)}\nknown keys: {known}")

    settings = dict(CONFIG_DEFAULTS)
    settings.update(loaded)
    base = path.parent.resolve()

    if settings["range"] is None:
        raise ConfigError("'range' is required, e.g. range: 53-57")

    header = str(settings["header"]).strip()
    if not header.isalpha():
        raise ConfigError(f"'header' must be alphabetic, got {settings['header']!r}")

    indices = parse_index_range(settings["range"], header)

    analysis_type = str(settings["type"])
    if analysis_type not in ANALYSIS_TYPES:
        raise ConfigError(f"'type' must be one of {', '.join(ANALYSIS_TYPES)}, got {analysis_type!r}")

    canonical = {name.lower(): name for name in DETECTOR_PEAKS}
    detectors = as_string_list(settings["detectors"], "detectors")
    if not detectors:
        raise ConfigError("'detectors' selects nothing")
    unknown_detectors = [d for d in detectors if d.lower() not in canonical]
    if unknown_detectors:
        known = ", ".join(DETECTOR_PEAKS)
        raise ConfigError(f"unknown detector(s) {', '.join(unknown_detectors)}; known: {known}")
    detectors = [canonical[d.lower()] for d in detectors]

    peaks = as_string_list(settings["peaks"], "peaks") if settings["peaks"] is not None else None

    data_dir = resolve_path(settings["data_dir"], base) if settings["data_dir"] else default_data_dir().resolve()

    if settings["output_dir"]:
        output_dir = resolve_path(settings["output_dir"], base)
    else:
        output_dir = data_dir.parent / "analysis" / analysis_type

    try:
        gpu_id = int(settings["gpu_id"])
    except (TypeError, ValueError):
        raise ConfigError(f"'gpu_id' must be an integer, got {settings['gpu_id']!r}") from None

    return SimpleNamespace(
        indices=indices,
        header=header,
        data_dir=data_dir,
        output_dir=output_dir,
        detectors=detectors,
        peaks=peaks,
        type=analysis_type,
        gpu_id=gpu_id,
        overwrite=bool(settings["overwrite"]),
        ignore_dedicated=bool(settings["ignore_dedicated"]),
        qmaps=build_default_qmaps(settings["qmaps"], detectors),
        extra_args=as_string_list(settings["extra_args"], "extra_args"),
        boost_corr=str(settings["boost_corr"]),
        dry_run=bool(settings["dry_run"]),
    )


def find_datasets(
    data_dir: Path,
    indices: Set[int],
    detectors: Sequence[str],
    header_letter: str,
) -> List[Dataset]:
    """Collect the acquisition folders whose header index is in ``indices``."""
    pattern = folder_pattern(header_letter)
    strip = header_letter + header_letter.swapcase()
    by_lower = {det.lower(): det for det in detectors}
    found: List[Dataset] = []
    for entry in sorted(data_dir.iterdir()):
        if not entry.is_dir():
            continue
        match = pattern.match(entry.name)
        if match is None:
            continue
        index = int(match.group("header").lstrip(strip))
        if index not in indices:
            continue
        # The detector is the last underscore-separated token before the repeat.
        token = match.group("middle").rsplit("_", 1)[-1]
        detector = by_lower.get(token.lower())
        if detector is None:
            continue
        found.append(
            Dataset(
                folder=entry,
                header=match.group("header"),
                index=index,
                repeat=match.group("repeat"),
                detector=detector,
            )
        )
    return found


def find_raw_file(dataset: Dataset) -> Optional[Path]:
    """Locate the raw data file inside a dataset folder."""
    for extension in RAW_EXTENSIONS:
        candidate = dataset.folder / f"{dataset.folder.name}{extension}"
        if candidate.is_file():
            return candidate
    # Fall back to any non-metadata file with a known extension.
    for entry in sorted(dataset.folder.iterdir()):
        if not entry.is_file() or "_metadata" in entry.name:
            continue
        if entry.name.endswith(RAW_EXTENSIONS):
            return entry
    return None


def dedicated_qmap_stem(detector: str, peak: Optional[str], header: str, repeat: str) -> str:
    """Build the strict dedicated-qmap name, minus the extension."""
    peak_token = f"{peak}_" if peak else ""
    return f"{detector.lower()}_qmap_{peak_token}{header}_{repeat}"


def resolve_qmap(
    dataset: Dataset,
    peak: Optional[str],
    data_dir: Path,
    files_by_lower: Dict[str, Path],
    defaults: Dict[Tuple[str, Optional[str]], str],
    ignore_dedicated: bool,
) -> Tuple[Optional[Path], str]:
    """Pick the qmap for one job.

    Returns ``(path, origin)``; ``path`` is ``None`` when nothing was found, in
    which case ``origin`` describes what was looked for.
    """
    if not ignore_dedicated:
        stem = dedicated_qmap_stem(dataset.detector, peak, dataset.header, dataset.repeat)
        for extension in QMAP_EXTENSIONS:
            hit = files_by_lower.get(f"{stem}{extension}".lower())
            if hit is not None:
                return hit, "dedicated"

    default = defaults.get((dataset.detector, peak))
    if default is None:
        return None, f"no default qmap configured for {dataset.detector}/{peak or 'no peak'}"

    # An explicit path is used as given; a bare name is looked up in the data dir.
    if os.sep in default or default.startswith("~"):
        candidate = Path(default).expanduser()
        if not candidate.is_absolute():
            candidate = (data_dir / candidate).resolve()
        if candidate.is_file():
            return candidate, "default"
        return None, f"default qmap not found: {candidate}"

    hit = files_by_lower.get(default.lower())
    if hit is not None:
        return hit, "default"
    return None, f"default qmap not found in {data_dir}: {default}"


def build_jobs(
    datasets: Sequence[Dataset],
    peaks_filter: Optional[Sequence[str]],
    data_dir: Path,
    files_by_lower: Dict[str, Path],
    defaults: Dict[Tuple[str, Optional[str]], str],
    ignore_dedicated: bool,
) -> Tuple[List[Job], List[str]]:
    """Turn datasets into jobs, collecting a message for everything skipped."""
    jobs: List[Job] = []
    skipped: List[str] = []

    for dataset in datasets:
        raw = find_raw_file(dataset)
        if raw is None:
            skipped.append(f"{dataset.folder.name}: no raw data file in folder")
            continue

        for peak in DETECTOR_PEAKS.get(dataset.detector, (None,)):
            if peaks_filter is not None and peak is not None and peak not in peaks_filter:
                continue
            qmap, origin = resolve_qmap(dataset, peak, data_dir, files_by_lower, defaults, ignore_dedicated)
            if qmap is None:
                tail = f" {peak}" if peak else ""
                skipped.append(f"{dataset.folder.name}{tail}: {origin}")
                continue
            jobs.append(Job(dataset=dataset, peak=peak, raw=raw, qmap=qmap, qmap_origin=origin))

    return jobs, skipped


def build_command(job: Job, config: SimpleNamespace) -> List[str]:
    """Assemble the boost_corr command line for one job."""
    command = [
        config.boost_corr,
        "-r",
        str(job.raw),
        "-q",
        str(job.qmap),
        "-t",
        config.type,
        "-o",
        str(config.output_dir),
        "-i",
        str(config.gpu_id),
    ]
    if job.peak:
        command += ["-u", job.peak]
    if config.overwrite:
        command.append("-w")
    command += list(config.extra_args)
    return command


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if any(flag in argv for flag in ("-h", "--help")):
        print(__doc__)
        return 0
    if "--write-example" in argv:
        print(EXAMPLE_CONFIG, end="")
        return 0
    if len(argv) > 1:
        print(f"error: expected at most one argument (the config file), got {len(argv)}", file=sys.stderr)
        return 2

    config_path = Path(argv[0]).expanduser() if argv else Path(DEFAULT_CONFIG_NAME)

    try:
        config = load_config(config_path)
    except ConfigError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if not config.data_dir.is_dir():
        print(f"error: data directory does not exist: {config.data_dir}", file=sys.stderr)
        return 2

    files_by_lower = {p.name.lower(): p for p in config.data_dir.iterdir() if p.is_file()}

    datasets = find_datasets(config.data_dir, config.indices, config.detectors, config.header)
    if not datasets:
        low, high = min(config.indices), max(config.indices)
        span = f"{config.header}{low:04d}-{config.header}{high:04d}"
        print(f"No dataset folders for {span} in {config.data_dir}", file=sys.stderr)
        return 1

    jobs, skipped = build_jobs(
        datasets,
        config.peaks,
        config.data_dir,
        files_by_lower,
        config.qmaps,
        config.ignore_dedicated,
    )

    low, high = min(config.indices), max(config.indices)
    print(f"config:     {config_path.resolve()}")
    print(f"data dir:   {config.data_dir}")
    print(f"output dir: {config.output_dir}")
    print(f"range:      {config.header}{low:04d}-{config.header}{high:04d}")
    print(f"type:       {config.type}   gpu-id: {config.gpu_id}   overwrite: {config.overwrite}")
    print(f"{len(datasets)} dataset(s) -> {len(jobs)} job(s)")
    for job in jobs:
        print(f"  {job.label:<58s} qmap={job.qmap.name} ({job.qmap_origin})")
    if skipped:
        print(f"skipping {len(skipped)}:")
        for message in skipped:
            print(f"  {message}")
    print()

    if config.dry_run:
        for job in jobs:
            print(shlex.join(build_command(job, config)))
        return 0

    if not jobs:
        return 1

    config.output_dir.mkdir(parents=True, exist_ok=True)

    failures: List[Tuple[str, int]] = []
    for number, job in enumerate(jobs, start=1):
        command = build_command(job, config)
        print(f"===== [{number}/{len(jobs)}] {job.label} =====")
        print(shlex.join(command), flush=True)
        result = subprocess.run(command)
        if result.returncode != 0:
            print(f"FAILED ({result.returncode}): {job.label}", file=sys.stderr)
            failures.append((job.label, result.returncode))
        print(flush=True)

    print(f"===== done: {len(jobs) - len(failures)}/{len(jobs)} succeeded =====")
    for label, code in failures:
        print(f"  FAILED ({code}): {label}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
