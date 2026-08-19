"""
Turn ``configs/spec_template.yml`` into the header lines and columns of a scan.

This is the layer that lets users change what a ``.spec`` file contains without
editing Python: the template names the header entries, the ``#O``/``#P``
positioner snapshot, and the data columns, and this module resolves those names
against the live ``oregistry``.

Used by the Ophyd scans in :mod:`ophyd_scan`; the actual file writing lives in
:mod:`ophyd_spec_writer`, which has no hardware imports at all.

Note on loading: this deliberately uses ``yaml.safe_load`` rather than
``apsbits.utils.config_loaders.load_config``. That helper merges whatever it
reads into a process-global ``_iconfig`` and overwrites ``ICONFIG_PATH`` /
``INSTRUMENT_PATH``, so pointing it at the SPEC template would corrupt the
session's real instrument config.
"""

import os
from pathlib import Path

import yaml
from apsbits.core.instrument_init import oregistry

#: Shipped default, edited in place by users (the repo is an editable install).
DEFAULT_TEMPLATE = Path(__file__).resolve().parents[2] / "configs" / "spec_template.yml"

#: Sources computed by the scan rather than read from a signal.
_SPECIAL_SOURCES = ("motor", "motor_setpoint", "epoch", "epoch_float")


# =============================================================================
# Loading
# =============================================================================


def template_path(explicit=None):
    """Resolve which template to use: argument, then ``$ID8_SPEC_TEMPLATE``, then default."""
    if explicit:
        return Path(os.path.expanduser(str(explicit)))
    from_env = os.environ.get("ID8_SPEC_TEMPLATE")
    if from_env:
        return Path(os.path.expanduser(from_env))
    return DEFAULT_TEMPLATE


def load_template(path=None):
    """Read a template file. Returns ``(config_dict, resolved_path)``."""
    resolved = template_path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"SPEC template not found: {resolved}")
    with open(resolved, "r") as handle:
        config = yaml.safe_load(handle) or {}
    if not isinstance(config, dict):
        raise ValueError(f"SPEC template {resolved} must be a mapping at the top level")
    return config, resolved


# =============================================================================
# Source resolution
# =============================================================================


def _resolve_object(dotted, motor, det):
    """``"det.stats1.total"`` / ``"tetramm1.current1.mean_value"`` -> ophyd object.

    ``det.`` and ``motor.`` refer to the objects passed to *this* scan, which is
    what lets a single template serve both eiger4M and lambda2M.
    """
    head, _, rest = str(dotted).partition(".")
    if head == "det":
        if det is None:
            raise KeyError(f"{dotted!r} needs a detector but this scan has none")
        obj = det
    elif head == "motor":
        obj = motor
    else:
        try:
            obj = oregistry[head]
        except Exception as exc:
            raise KeyError(f"no device {head!r} in the oregistry ({dotted!r})") from exc
    for part in filter(None, rest.split(".")):
        try:
            obj = getattr(obj, part)
        except AttributeError as exc:
            raise KeyError(f"{dotted!r}: {obj.name!r} has no attribute {part!r}") from exc
    return obj


def _substitute(text, context):
    """Apply ``{motor}``/``{det}``/... substitutions with a helpful error."""
    try:
        return str(text).format(**context)
    except KeyError as exc:
        raise KeyError(
            f"unknown substitution {exc} in {text!r}; "
            f"available: {', '.join(sorted(context))}"
        ) from exc


def _read_position(obj):
    """Positioners expose ``.position``; plain signals only ``.get()``."""
    if hasattr(obj, "position"):
        return obj.position
    return obj.get()


# =============================================================================
# Rendered result
# =============================================================================


class _Column:
    """One data column: either a computed special or a signal to read."""

    __slots__ = ("label", "special", "signal")

    def __init__(self, label, special=None, signal=None):
        self.label = label
        self.special = special
        self.signal = signal

    def value(self, position, setpoint, elapsed):
        if self.special == "motor":
            return position
        if self.special == "motor_setpoint":
            return setpoint
        if self.special == "epoch":
            return round(elapsed)
        if self.special == "epoch_float":
            return elapsed
        return self.signal.get()


class RenderedSpec:
    """Everything a scan needs to write its SPEC block, resolved against hardware."""

    def __init__(self, columns, metadata, comments, positioner_names,
                 positioner_positions, source):
        self._columns = columns
        self.metadata = metadata
        self.comments = comments
        self.positioner_names = positioner_names
        self.positioner_positions = positioner_positions
        self.source = source

    @property
    def labels(self):
        return [c.label for c in self._columns]

    @property
    def counter_indices(self):
        """Indices of columns read from hardware, i.e. not motor/setpoint/epoch.

        The live table already shows position and elapsed time in their own
        columns, so it prints only these.
        """
        return [i for i, c in enumerate(self._columns) if c.special is None]

    @property
    def counter_labels(self):
        return [self._columns[i].label for i in self.counter_indices]

    def read(self, position, setpoint, elapsed):
        """Values for one scan point, in column order."""
        return [c.value(position, setpoint, elapsed) for c in self._columns]

    def __repr__(self):
        return (f"<RenderedSpec {len(self._columns)} columns, "
                f"{len(self.positioner_names)} positioners, from {self.source.name}>")


# =============================================================================
# Rendering
# =============================================================================


def render(motor, det=None, scan_type="", num_points=0, count_time=0.0,
           image_file="", comment="", template=None, extra=None):
    """Resolve a template against the current hardware.

    args:
        motor: the scanned positioner
        det: detector for ``det.*`` sources, or None
        scan_type, num_points, count_time, image_file: substitution values
        comment: free-text note for this scan. Written as the FIRST ``#MD``
            line so it is easy to spot and easy to Ctrl+F for. Also available
            to the template as ``{comment}``.
        template: path override; otherwise ``$ID8_SPEC_TEMPLATE`` or the default
        extra: additional substitutions (a future d2scan passes ``motor2`` here)

    returns:
        :class:`RenderedSpec`

    Entries marked ``optional: true`` are skipped with a warning when their
    device is missing; anything else raises here -- before the scan starts --
    rather than failing partway through.
    """
    config, source = load_template(template)

    context = {
        "motor": motor.name,
        "det": det.name if det is not None else "",
        "scan_type": scan_type,
        "num_points": num_points,
        "count_time": count_time,
        "image_file": image_file,
        "comment": "" if comment is None else str(comment),
    }
    context.update(extra or {})

    positioner_names, positioner_positions = _render_positioners(config, motor, det)

    # The scan comment is an ordinary template entry (shipped first in the
    # default template), so it can be moved, renamed or removed like any other.
    # The fallback below only fires when the template did not write the note
    # anywhere -- a note the user bothered to type should never vanish silently.
    # Testing the rendered values rather than the key name matters: a template
    # that renames the key to 'note' has still recorded it, and prepending a
    # second copy under 'comment' would duplicate it.
    metadata = _render_metadata(config, motor, det, context)
    comment_text = context["comment"].strip()
    if comment_text and comment_text not in metadata.values():
        metadata = {"comment": comment_text, **metadata}
    comments = [_substitute(c, context) for c in config.get("comments", []) or []]
    columns = _render_columns(config, motor, det, context)

    if not columns:
        raise ValueError(f"SPEC template {source} defines no usable columns")
    if columns[-1].special in ("motor", "motor_setpoint"):
        print(
            f"WARNING: the last SPEC column is {columns[-1].label!r}, which the "
            "specr_py viewer will use as the default Y axis. Put your primary "
            "counter last in the template's 'columns:' list."
        )

    return RenderedSpec(columns, metadata, comments, positioner_names,
                        positioner_positions, source)


def _render_positioners(config, motor, det):
    """``#O`` names and ``#P`` positions, with the scanned motor guaranteed present."""
    names, positions = [], []
    for dotted in config.get("positioners", []) or []:
        try:
            obj = _resolve_object(dotted, motor, det)
            positions.append(_read_position(obj))
            names.append(obj.name)
        except Exception as exc:
            # A missing stage must never stop a scan -- it only costs us context.
            print(f"WARNING: SPEC positioner {dotted!r} unavailable ({exc}); skipped.")
    if motor.name not in names:
        try:
            positions.append(motor.position)
            names.append(motor.name)
        except Exception as exc:
            print(f"WARNING: could not read scanned motor position ({exc}).")
    return names, positions


def _join_sources(spec, motor, det):
    """Read several PVs and join them into one header value.

    A mapping keeps the names, which is usually what you want in a header::

        {min_x: det.roi1.min_xyz.min_x, size_x: det.roi1.size.x}
        -> "min_x=730 size_x=100"

    A list gives bare values in order::

        [det.roi1.min_xyz.min_x, det.roi1.size.x]  ->  "730 100"

    Readers treat a ``#MD`` value as free text and split only on the *first*
    ``=``, so the ``name=value`` pairs survive intact.
    """
    if isinstance(spec, dict):
        return " ".join(
            f"{name}={_resolve_object(path, motor, det).get()}"
            for name, path in spec.items()
        )
    if isinstance(spec, (list, tuple)):
        return " ".join(
            str(_resolve_object(path, motor, det).get()) for path in spec
        )
    raise TypeError(f"'sources' must be a mapping or a list, got {type(spec).__name__}")


def _render_metadata(config, motor, det, context):
    """``#MD key = value`` entries."""
    metadata = {}
    for entry in config.get("metadata", []) or []:
        if not isinstance(entry, dict) or "key" not in entry:
            print(f"WARNING: skipping malformed SPEC metadata entry {entry!r}.")
            continue
        key = _substitute(entry["key"], context)
        optional = bool(entry.get("optional"))
        try:
            if "sources" in entry:
                # Several PVs collapsed onto one #MD line. Stays generic -- the
                # template names the PVs, so nothing here knows what an ROI is.
                value = _join_sources(entry["sources"], motor, det)
            elif "source" in entry:
                # Name the individual PV, not a parent device: a bare device's
                # .get() returns a namedtuple of every component, which is not
                # a value anyone wants in a header line.
                value = _resolve_object(entry["source"], motor, det).get()
            else:
                value = _substitute(entry.get("value", ""), context)
        except Exception as exc:
            if optional:
                print(f"WARNING: SPEC metadata {key!r} unavailable ({exc}); skipped.")
                continue
            raise
        text = str(value).strip()
        if entry.get("skip_if_empty") and not text:
            continue
        metadata[key] = text
    return metadata


def _render_columns(config, motor, det, context):
    """Data columns, in template order."""
    columns = []
    for entry in config.get("columns", []) or []:
        if isinstance(entry, dict) and "sources" in entry:
            # Would emit a cell containing spaces and shift every later column.
            raise ValueError(
                f"SPEC column {entry.get('label', entry['sources'])!r} uses "
                "'sources'; that is for metadata only, because a data cell must "
                "be a single number. Use one 'source' per column instead."
            )
        if not isinstance(entry, dict) or "source" not in entry:
            print(f"WARNING: skipping malformed SPEC column entry {entry!r}.")
            continue
        source = str(entry["source"])
        label = _substitute(entry.get("label", source), context)
        if source in _SPECIAL_SOURCES:
            columns.append(_Column(label, special=source))
            continue
        try:
            signal = _resolve_object(source, motor, det)
        except Exception as exc:
            if entry.get("optional"):
                print(f"WARNING: SPEC column {label!r} unavailable ({exc}); skipped.")
                continue
            raise
        columns.append(_Column(label, signal=signal))
    return columns
