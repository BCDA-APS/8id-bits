#!/usr/bin/env python
"""Prove a plan CSV means the same thing as the YAML it replaces.

The converter's whole contract is::

    read_plan_csv(x.csv)  ==  yaml.safe_load(x.yaml)

so the test is a deep dict comparison, not an eyeball of the output. Reading 200
lines of nested YAML is how a wrong att_level survives; a dict diff is how it
does not. Runs anywhere -- no EPICS, no /gdata, no beam.

    python scripts/check_plan_csv.py                       # the bundled examples
    python scripts/check_plan_csv.py a.csv a.yaml [...]    # your own pair(s)

Exit status is 0 only if every pair matches, so it can gate a commit.

Pair the CSV with a YAML a HUMAN wrote. A YAML dumped from the converter proves
only that the converter equals itself.
"""

import difflib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import yaml  # noqa: E402

from id8_common.plans.acquire.validators import read_yaml  # noqa: E402

EXAMPLES = REPO / "Developer_Notes" / "csv_examples"
GREEN, RED, RESET = "\033[92m", "\033[91m", "\033[0m"


def normalise(obj):
    """Sorted-key YAML text, so a diff reads as a diff instead of `False`."""
    return yaml.dump(obj, sort_keys=True, default_flow_style=False).splitlines()


def compare(csv_path, yaml_path):
    csv_path, yaml_path = Path(csv_path), Path(yaml_path)
    print(f"\n=== {csv_path.name}  vs  {yaml_path.name} ===")
    try:
        # Through read_yaml, not read_plan_csv: this also proves the dispatch
        # in validators.py routes by suffix the way the real run will.
        from_csv = read_yaml(csv_path)
    except Exception as exc:  # noqa: BLE001 -- a converter failure IS the result
        print(f"{RED}CSV FAILED TO PARSE: {type(exc).__name__}: {exc}{RESET}")
        return False
    from_yaml = read_yaml(yaml_path)

    if from_csv == from_yaml:
        print(f"{GREEN}MATCH{RESET} — the CSV and the YAML produce identical dicts")
        return True

    print(f"{RED}MISMATCH{RESET} — unified diff, '-' is the YAML, '+' is the CSV:")
    diff = difflib.unified_diff(
        normalise(from_yaml), normalise(from_csv),
        fromfile=str(yaml_path), tofile=str(csv_path), lineterm="",
    )
    for line in diff:
        colour = GREEN if line.startswith("+") else RED if line.startswith("-") else ""
        print(f"  {colour}{line}{RESET}" if colour else f"  {line}")
    return False


def main(argv):
    if len(argv) >= 2:
        if len(argv) % 2:
            print("Give CSV/YAML pairs: check_plan_csv.py a.csv a.yaml [b.csv b.yaml ...]")
            return 2
        pairs = list(zip(argv[0::2], argv[1::2]))
    else:
        pairs = [
            (EXAMPLES / "sample_info.csv", EXAMPLES / "sample_info.expected.yaml"),
            (EXAMPLES / "measurement_info.csv", EXAMPLES / "measurement_info.expected.yaml"),
        ]

    results = [compare(c, y) for c, y in pairs]
    ok = all(results)
    print(f"\n{GREEN if ok else RED}{sum(results)}/{len(results)} pair(s) matched{RESET}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
