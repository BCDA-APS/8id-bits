"""Read a measurement plan written as CSV, and return the dict the YAML would have given.

Why this exists: the people who set up a run are beamline scientists, not
programmers, and a batch table of attenuation x repeats genuinely wants to be a
spreadsheet. YAML indentation is a real barrier -- and worse, a *silent* one: a
key indented two spaces too far is still valid YAML, it just lands on the wrong
parent, and nothing notices until a required field turns up missing at run time.

The contract of this module is one line:

    read_plan_csv(path)  ==  yaml.safe_load(open(equivalent_yaml))

Nothing downstream is allowed to know which format was on disk. This returns the
same nested dict ``yaml.safe_load`` returns, so ``expand_measurements()``,
``REQUIRED_LEG_FIELDS``, ``validate_acq_time`` and ``validate_qmap_exists`` all
run exactly as they do today and do all the semantic checking. This module's only
job is shape and type: fail loudly on a malformed file, and never silently invent
a value. See ``scripts/check_plan_csv.py`` for the equivalence check.

No file is ever written. There is no generated YAML on disk to go stale, and
nothing here touches hardware, EPICS or the registry.

----------------------------------------------------------------- file formats

A row whose first cell starts with ``#`` opens a section. A wholly blank row is
ignored, so the trailing empty rows Excel leaves behind are harmless.

``sample_info.csv``::

    #DEFAULTS
    inner_motor,huber.x
    outer_motor,huber.y

    #SAMPLES
    index,header,inner_center,outer_center,inner_pts,outer_pts
    1,A,-0.1,0.1,21,31
    2,B,-0.1,0.1,21,31

  -> {"defaults": {...}, "samples": {"sample_1": {...}, "sample_2": {...}}}

  A blank cell is OMITTED from the sample rather than written as an empty
  string, so the ``defaults:`` block applies -- which is what a blank cell in a
  spreadsheet means to the person who left it blank.

``measurement_info.csv`` -- which samples run which protocols, in order::

    #REPEAT,1
    #BATCH
    sample,protocol
    1,test_protocol
    2,
    3,

  -> {"loop_order": ..., "runs": [...], "protocols": {...}}

  A blank cell inherits from the row above (forward fill), again matching
  spreadsheet intuition. Every fill is printed, so an accidental blank shows up
  in the terminal instead of quietly running the wrong protocol.

  Each row becomes ONE run block, in file order, so what you read top-to-bottom
  in the spreadsheet is the order the measurements happen. ``protocol`` names a
  sibling ``<name>.csv`` holding the protocol itself.

``<protocol>.csv`` -- one protocol family::

    #DETECTOR
    device,"rigaku3M_epics, eiger4M"
    mode,"EPICS, Internal Series"
    qmap_file,"rigaku3m_qmap_default.hdf, eiger4m_qmap_default.hdf"
    analysis_type,"Multitau, Both"

    #MEASUREMENT
    wait_time,0
    sample_move,1
    position_reset,1
    acq_time,1
    num_frames,3000

    #BATCH
    att_level,num_repeat
    1,1
    2,1

  ``#DETECTOR`` values are comma-separated, one entry per leg. A single value is
  broadcast to every leg; any other length mismatch is an error rather than a
  guess.

  WHICH SHAPE IS EMITTED. This tree has two protocol schemas, and the number of
  detectors picks between them:

    one detector   -> the SERIAL shape master_plan.py consumes: `detector:` and
                      `mode:` singular, at protocol level.
    two or more    -> the PARALLEL shape dual_master_plan.py consumes: a
                      `detectors:` list, one entry per leg.

  Add ``parallel,yes`` to ``#MEASUREMENT`` to force the list shape for a single
  detector -- a one-leg parallel run is legal (it becomes its own shutter owner)
  and is the right way to smoke-test one leg of a trio.

  ``#MEASUREMENT`` keys are split by ``LEG_FIELDS`` below: a per-leg field
  (``acq_time``, ``num_frames``, ...) is copied into EVERY leg, because
  ``master_plan.REQUIRED_LEG_FIELDS`` requires it per leg. Everything else stays
  at protocol level.

  ``#BATCH`` produces one protocol per row: ``test_protocol`` with rows
  att_level 1 and 2 yields protocols ``test_protocol_att1`` and
  ``test_protocol_att2``, both listed in the run block that named
  ``test_protocol``.

Always emits the parallel ``detectors:`` shape, even for one detector.
master_plan accepts both shapes, and one shape means one code path to reason
about.
"""

import csv
from pathlib import Path

__all__ = ["read_plan_csv", "PlanCSVError"]


class PlanCSVError(ValueError):
    """A CSV plan file could not be read. Always names the file and the row."""


#: Fields that belong to an individual detector leg, not to the protocol.
#: Mirrors master_plan.REQUIRED_LEG_FIELDS plus the optional per-leg fields it
#: lifts into the spec. Anything not listed here stays at protocol level.
LEG_FIELDS = {
    "device",
    "mode",
    "acq_time",
    "acq_period",
    "num_frames",
    "num_segments",
    "trigger_period",
    "qmap_file",
    "analysis_type",
    "label",
    "select_device",
    "start_timeout",
    "stop_timeout",
    "hdf_timeout",
    "workflow_name",
}

#: CSV spellings that differ from the YAML key. Kept deliberately short: a long
#: alias table is a sign the CSV header should have been renamed instead.
KEY_ALIASES = {
    "num_repeat": "num_repeats",
    "repeat": "num_repeats",
}

#: Booleans. YAML 1.1 loads a bare `yes` as Python True, and validators.yes_no()
#: is built around that, so the CSV must produce real bools -- not "1", not "yes"
#: -- or a dict comparison against the YAML will not match.
_TRUE = {"1", "yes", "y", "true", "on"}
_FALSE = {"0", "no", "n", "false", "off"}

#: Protocol-level fields that are genuinely boolean in the YAML.
BOOL_FIELDS = {"sample_move", "position_reset", "select_device", "parallel"}


def _coerce(value, key=None):
    """Turn one CSV cell into the Python type ``yaml.safe_load`` would produce."""
    text = str(value).strip()
    if text == "":
        return None
    if key in BOOL_FIELDS:
        low = text.lower()
        if low in _TRUE:
            return True
        if low in _FALSE:
            return False
        raise PlanCSVError(f"{key}: expected yes/no (or 1/0), got {text!r}")
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        pass
    return text


def _key(name):
    return KEY_ALIASES.get(str(name).strip(), str(name).strip())


def _rows(path):
    """Every row of the file, blank rows dropped, cells stripped."""
    path = Path(path)
    try:
        with open(path, newline="") as handle:
            raw = list(csv.reader(handle))
    except OSError as exc:
        raise PlanCSVError(f"cannot read {path}: {exc}") from exc
    out = []
    for lineno, row in enumerate(raw, start=1):
        cells = [c.strip() for c in row]
        if not any(cells):
            continue  # Excel leaves trailing empty rows; ignore them
        out.append((lineno, cells))
    return out


#: The only ``#`` rows that open a section. A fixed vocabulary, because the
#: drafted files also write DATA rows with a leading '#' (``#device,"a, b"``),
#: and a parser that treated every '#' row as a section would silently read
#: those as empty sections and drop the detectors. Anything else starting with
#: '#' is either such a data row or a free-text comment -- see _sections().
KNOWN_SECTIONS = {
    "DEFAULTS",
    "SAMPLES",
    "DETECTOR",
    "MEASUREMENT",
    "BATCH",
    "REPEAT",
    "LOOP_ORDER",
}


def _sections(path):
    """Split rows into ``{SECTION: [(lineno, cells), ...]}`` on ``#NAME`` rows.

    Three things can start with '#', and they are told apart by name and shape:

    * ``#BATCH`` / ``#DETECTOR`` / ... -- a name in KNOWN_SECTIONS opens a
      section. ``#REPEAT,1`` is a section name that also carries its value.
    * ``#device,"rigaku3M_epics, eiger4M"`` -- not a known section, but it has a
      value, so it is a data row whose key happens to be written with a '#'.
      The '#' is stripped and the row is kept.
    * ``# anything else`` with no value -- a free-text comment, ignored. Put
      notes in the spreadsheet this way and they cost nothing.

    Rows before the first marker go under ``""``, so a file may open with bare
    ``key,value`` lines.
    """
    sections = {"": []}
    current = ""
    for lineno, cells in _rows(path):
        head = cells[0]
        if head.startswith("#"):
            name = head.lstrip("#").strip().upper()
            rest = [c for c in cells[1:] if c]
            if name in KNOWN_SECTIONS:
                sections.setdefault(name, [])
                if rest:
                    # `#REPEAT,1` is a marker AND a value: keep the value with it.
                    sections[name].append((lineno, cells))
                current = name
                continue
            if rest:
                # A data row written with a leading '#'. Strip it and keep going.
                sections[current].append((lineno, [head.lstrip("#").strip()] + cells[1:]))
            # else: a comment row. Dropped.
            continue
        sections[current].append((lineno, cells))
    return sections


def _pairs(rows, path):
    """``key,value`` rows -> dict, with types coerced and blanks dropped."""
    out = {}
    for lineno, cells in rows:
        if len(cells) < 2:
            raise PlanCSVError(f"{path}:{lineno}: expected 'key,value', got {cells!r}")
        key = _key(cells[0])
        value = _coerce(cells[1], key)
        if value is not None:
            out[key] = value
    return out


def _table(rows, path, section):
    """A header row plus data rows -> list of dicts. Blank cells become None."""
    if not rows:
        return []
    header_lineno, header = rows[0]
    keys = [_key(c) for c in header]
    if not keys or not keys[0]:
        raise PlanCSVError(f"{path}:{header_lineno}: #{section} needs a header row")
    table = []
    for lineno, cells in rows[1:]:
        padded = cells + [""] * (len(keys) - len(cells))
        record = {}
        for key, cell in zip(keys, padded):
            if not key:
                continue
            record[key] = _coerce(cell, key)
        table.append((lineno, record))
    return table


# --------------------------------------------------------------- sample_info


def read_sample_csv(path):
    """``sample_info.csv`` -> the dict ``sample_info.yaml`` would have given."""
    path = Path(path)
    sections = _sections(path)
    defaults = _pairs(sections.get("DEFAULTS", []), path)

    rows = _table(sections.get("SAMPLES", []), path, "SAMPLES")
    if not rows:
        raise PlanCSVError(f"{path}: no #SAMPLES rows found")

    samples = {}
    for lineno, record in rows:
        index = record.get("index")
        if index is None:
            raise PlanCSVError(f"{path}:{lineno}: sample row has no index")
        name = f"sample_{index}"
        if name in samples:
            raise PlanCSVError(f"{path}:{lineno}: sample index {index} appears twice")
        # A blank cell is omitted, not written as "", so defaults: still applies.
        samples[name] = {k: v for k, v in record.items() if k != "index" and v is not None}

    out = {"samples": samples}
    if defaults:
        out["defaults"] = defaults
    return out


# ---------------------------------------------------------------- protocols


def _split_list(value, count, key, path):
    """One ``#DETECTOR`` cell -> a list of ``count`` per-leg values."""
    parts = [p.strip() for p in str(value).split(",")]
    parts = [p for p in parts if p != ""]
    if len(parts) == count:
        return parts
    if len(parts) == 1:
        return parts * count  # one value, broadcast to every leg
    raise PlanCSVError(
        f"{path}: '{key}' has {len(parts)} value(s) but there are {count} detector(s). "
        f"Give one value per detector, or a single value to use for all."
    )


def read_protocol_csv(path, stem=None):
    """``<protocol>.csv`` -> ``{protocol_name: protocol_dict, ...}``.

    One entry per ``#BATCH`` row. With no ``#BATCH`` section the file defines a
    single protocol under its own file stem.
    """
    path = Path(path)
    stem = stem or path.stem
    sections = _sections(path)

    detector = _pairs(sections.get("DETECTOR", []), path)
    if "device" not in detector:
        raise PlanCSVError(f"{path}: #DETECTOR has no 'device' row")
    devices = [d.strip() for d in str(detector["device"]).split(",") if d.strip()]
    n_legs = len(devices)

    per_leg_lists = {"device": devices}
    for key, value in detector.items():
        if key == "device":
            continue
        per_leg_lists[key] = _split_list(value, n_legs, key, path)

    measurement = _pairs(sections.get("MEASUREMENT", []), path)
    leg_shared = {k: v for k, v in measurement.items() if k in LEG_FIELDS}
    protocol_level = {k: v for k, v in measurement.items() if k not in LEG_FIELDS}

    legs = []
    for i in range(n_legs):
        leg = dict(leg_shared)
        for key, values in per_leg_lists.items():
            leg[key] = _coerce(values[i], key)
        legs.append(leg)

    # Default: one detector means the serial path, more than one means parallel.
    parallel = protocol_level.pop("parallel", n_legs > 1)

    batch = _table(sections.get("BATCH", []), path, "BATCH")
    if not batch:
        return {stem: {**protocol_level, **_shape_legs(legs, parallel)}}

    protocols = {}
    for lineno, record in batch:
        row = {k: v for k, v in record.items() if v is not None}
        if not row:
            continue
        name = _batch_name(stem, row)
        if name in protocols:
            raise PlanCSVError(f"{path}:{lineno}: two #BATCH rows both produce protocol '{name}'")
        protocols[name] = {
            **protocol_level,
            **row,
            **_shape_legs(legs, parallel),
        }
    if not protocols:
        raise PlanCSVError(f"{path}: #BATCH section has a header but no rows")
    return protocols


def _shape_legs(legs, parallel):
    """Return the protocol fragment carrying the detector(s), in the right schema.

    master_plan.py (serial) reads `detector`/`mode` singular at protocol level;
    dual_master_plan.py (dual and trio) reads a `detectors:` list. Emitting the
    wrong one is not a silent failure -- the serial path raises KeyError on
    'detector' and the parallel path raises on 'detectors' -- but picking by leg
    count means the common cases need no marker at all.
    """
    if parallel:
        return {"detectors": [dict(leg) for leg in legs]}

    if len(legs) != 1:
        raise PlanCSVError(
            f"the serial shape needs exactly one detector, got {len(legs)}. "
            f"Remove the extra detectors, or set parallel,yes in #MEASUREMENT."
        )

    leg = dict(legs[0])
    # The serial schema spells it `detector`; a leg spells it `device`.
    leg["detector"] = leg.pop("device")
    return leg


def _batch_name(stem, row):
    """``test_protocol`` + ``{att_level: 20, num_repeats: 1}`` -> ``test_protocol_att20``.

    num_repeats is excluded: it is how many times the protocol runs, not which
    protocol it is. If nothing else varies, the stem is used unchanged.
    """
    parts = []
    for key, value in row.items():
        if key == "num_repeats":
            continue
        short = key[:-6] if key.endswith("_level") else key
        parts.append(f"{short}{value}")
    return f"{stem}_{'_'.join(parts)}" if parts else stem


# ------------------------------------------------------------ measurement_info


def read_measurement_csv(path):
    """``measurement_info.csv`` -> the dict ``measurement_info.yaml`` would have given.

    Resolves each ``protocol`` cell to a sibling ``<name>.csv`` and folds every
    protocol it defines into one ``protocols:`` block.
    """
    path = Path(path)
    sections = _sections(path)

    settings = {}
    for marker, key in (("REPEAT", "repeats"), ("LOOP_ORDER", "loop_order")):
        rows = sections.get(marker, [])
        if rows:
            _, cells = rows[0]
            settings[key] = _coerce(cells[1] if len(cells) > 1 else "", key)
    repeats = settings.get("repeats", 1)

    rows = _table(sections.get("BATCH", []), path, "BATCH")
    if not rows:
        raise PlanCSVError(f"{path}: no #BATCH rows found (expected a sample,protocol table)")

    protocols = {}
    runs = []
    last = {"sample": None, "protocol": None}

    for lineno, record in rows:
        for field in ("sample", "protocol"):
            if record.get(field) is None:
                if last[field] is None:
                    raise PlanCSVError(
                        f"{path}:{lineno}: '{field}' is blank and there is no earlier row to copy."
                    )
                record[field] = last[field]
                # Announced, not silent: a stray blank cell is a real risk here.
                print(f"[plan_csv] {path.name}:{lineno}: blank '{field}' -> {record[field]!r} (from the row above)")
            last[field] = record[field]

        protocol_name = str(record["protocol"])
        protocol_path = path.parent / f"{protocol_name}.csv"
        if not protocol_path.exists():
            raise PlanCSVError(
                f"{path}:{lineno}: protocol '{protocol_name}' names no file -- "
                f"expected {protocol_path}"
            )
        defined = read_protocol_csv(protocol_path, stem=protocol_name)
        for name, body in defined.items():
            if name in protocols and protocols[name] != body:
                raise PlanCSVError(f"{path}:{lineno}: protocol '{name}' defined twice, differently")
            protocols[name] = body

        run = {
            "name": f"{protocol_name}_sample{record['sample']}",
            "samples": [record["sample"]],
            "protocols": sorted(defined),
            "repeats": repeats,
        }
        runs.append(run)

    out = {"runs": runs, "protocols": protocols}
    if "loop_order" in settings:
        out["loop_order"] = settings["loop_order"]
    return out


# ------------------------------------------------------------------ dispatch


def read_plan_csv(path):
    """Read any plan CSV. Picks the reader from the file name.

    ``sample_info*.csv`` is a sample table; everything else is a measurement
    plan. Called only from ``validators.read_yaml`` when it is handed a ``.csv``.
    """
    path = Path(path)
    if path.stem.lower().startswith("sample_info"):
        return read_sample_csv(path)
    return read_measurement_csv(path)
