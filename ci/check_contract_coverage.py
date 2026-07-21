#!/usr/bin/env python3
"""Fail closed if the frozen contract or private test denominator shrinks."""

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = (ROOT / "contract_coverage.toml").read_text(encoding="utf-8")


def manifest_path(field: str) -> Path:
    match = re.search(rf'(?m)^{field}\s*=\s*"([^"]+)"\s*$', MANIFEST)
    if not match:
        raise SystemExit(f"coverage manifest is missing {field!r}")
    return ROOT / match.group(1)


contract_text = manifest_path("contract").read_text(encoding="utf-8")
visible_text = manifest_path("tests").read_text(encoding="utf-8")
private_text = "\n".join(
    path.read_text(encoding="utf-8") for path in sorted((ROOT / "test").rglob("*.jl"))
)

contract_pattern = r"(?:API|CRY|MOD|BAS|SCF|OBS|RSP|DYN|PHN|POL|QEI)-\d{3}"
test_pattern = r"(?:UT|PT|RT|IT)-\d{3}"

contract_ids = set(re.findall(contract_pattern, contract_text))
visible_ids = set(re.findall(test_pattern, visible_text))
private_ids = set(re.findall(test_pattern, private_text))

coverage = {}
for contract_id, body in re.findall(
    rf'(?m)^"({contract_pattern})"\s*=\s*\[([^\]]*)\]\s*$', MANIFEST
):
    coverage[contract_id] = set(re.findall(test_pattern, body))

errors = []
missing_rows = contract_ids - coverage.keys()
extra_rows = coverage.keys() - contract_ids
if missing_rows:
    errors.append("missing contract rows: " + ", ".join(sorted(missing_rows)))
if extra_rows:
    errors.append("unknown contract rows: " + ", ".join(sorted(extra_rows)))

for contract_id, mapped_tests in sorted(coverage.items()):
    if not mapped_tests:
        errors.append(f"{contract_id} has an empty coverage row")
    unknown = mapped_tests - visible_ids
    if unknown:
        errors.append(f"{contract_id} maps unknown visible tests: {sorted(unknown)}")
    absent = mapped_tests - private_ids
    if absent:
        errors.append(f"{contract_id} maps tests absent from private suite: {sorted(absent)}")

unimplemented = visible_ids - private_ids
if unimplemented:
    errors.append("visible test IDs missing privately: " + ", ".join(sorted(unimplemented)))

spec_ids = set(re.findall(r"quantumdft-v[^`\s\"]+", contract_text + visible_text + MANIFEST))
if len(spec_ids) != 1:
    errors.append("specification ID is absent or inconsistent: " + repr(sorted(spec_ids)))

if errors:
    print("contract coverage FAILED", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    raise SystemExit(1)

print(
    f"contract coverage OK: {len(contract_ids)} contract IDs, "
    f"{len(visible_ids)} test IDs, spec {next(iter(spec_ids))}"
)
