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
motivation_text = (manifest_path("contract").parent / "MOTIVATION.md").read_text(
    encoding="utf-8"
)
private_files = {
    path: path.read_text(encoding="utf-8")
    for path in sorted((ROOT / "test").rglob("*.jl"))
    if "sentinels" not in path.parts
}
private_text = "\n".join(private_files.values())

contract_pattern = r"(?:API|CRY|MOD|BAS|SCF|OBS|RSP|DYN|PHN|POL|QEI)-\d{3}"
test_pattern = r"(?:UT|PT|RT|IT)-\d{3}"
verifier_test_pattern = r"VT-\d{3}"

contract_ids = set(re.findall(contract_pattern, contract_text))
visible_ids = set(re.findall(test_pattern, visible_text))

definition_pattern = re.compile(
    rf'@testset\s+"({test_pattern})(?:\s+[^"]*)?"\s+begin'
)
definitions = []
testset_bodies = {}
for path, text in private_files.items():
    matches = list(definition_pattern.finditer(text))
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        test_id = match.group(1)
        definitions.append((test_id, path))
        testset_bodies[test_id] = text[match.end():end]
private_ids = {test_id for test_id, _ in definitions}

visible_verifier_ids = set(re.findall(verifier_test_pattern, visible_text))
verifier_definitions = re.findall(
    rf'@testset\s+"({verifier_test_pattern})(?:\s+[^"]*)?"\s+begin', private_text
)

coverage = {}
for contract_id, body in re.findall(
    rf'(?m)^"({contract_pattern})"\s*=\s*\[([^\]]*)\]\s*$', MANIFEST
):
    coverage[contract_id] = set(re.findall(test_pattern, body))

errors = []
definition_counts = {test_id: 0 for test_id in private_ids}
for test_id, _ in definitions:
    definition_counts[test_id] += 1
duplicates = sorted(test_id for test_id, count in definition_counts.items() if count != 1)
if duplicates:
    errors.append("test IDs must be defined exactly once: " + ", ".join(duplicates))

for test_id, body in sorted(testset_bodies.items()):
    if not re.search(r"@test(?:_throws)?\b", body):
        errors.append(f"{test_id} has no executable assertion")
    if re.search(r"@test_(?:broken|skip)\b", body):
        errors.append(f"{test_id} contains a skipped or broken assertion")

if set(verifier_definitions) != visible_verifier_ids:
    errors.append(
        "verifier self-test definitions differ from TESTS.md: "
        f"defined={sorted(set(verifier_definitions))}, "
        f"declared={sorted(visible_verifier_ids)}"
    )
if len(verifier_definitions) != len(set(verifier_definitions)):
    errors.append("verifier self-test IDs must be defined exactly once")

missing_rows = contract_ids - coverage.keys()
extra_rows = coverage.keys() - contract_ids
if missing_rows:
    errors.append("missing contract rows: " + ", ".join(sorted(missing_rows)))
if extra_rows:
    errors.append("unknown contract rows: " + ", ".join(sorted(extra_rows)))

for contract_id, mapped_tests in sorted(coverage.items()):
    if not mapped_tests:
        errors.append(f"{contract_id} has an empty coverage row")
    if mapped_tests <= {"IT-005"}:
        errors.append(f"{contract_id} relies only on the umbrella IT-005 test")
    unknown = mapped_tests - visible_ids
    if unknown:
        errors.append(f"{contract_id} maps unknown visible tests: {sorted(unknown)}")
    absent = mapped_tests - private_ids
    if absent:
        errors.append(f"{contract_id} maps tests absent from private suite: {sorted(absent)}")

unimplemented = visible_ids - private_ids
if unimplemented:
    errors.append("visible test IDs missing privately: " + ", ".join(sorted(unimplemented)))
undeclared = private_ids - visible_ids
if undeclared:
    errors.append("private test IDs absent from TESTS.md: " + ", ".join(sorted(undeclared)))

spec_ids = set(re.findall(
    r"espressodft-v[^`\s\"]+", motivation_text + contract_text + visible_text + MANIFEST
))
if len(spec_ids) != 1:
    errors.append("specification ID is absent or inconsistent: " + repr(sorted(spec_ids)))

heldout_body = testset_bodies.get("IT-005", "")
for required_call in (
    "candidate_state", "energy", "forces", "stress", "assert_density_matches",
    "assert_gamma_bands_match", "dynamical_matrix", "phonon_modes",
    "born_effective_charges", "dielectric_tensor", "run_qe_input",
):
    if required_call not in heldout_body:
        errors.append(f"IT-005 does not exercise {required_call}")

if errors:
    print("contract coverage FAILED", file=sys.stderr)
    for error in errors:
        print(f"- {error}", file=sys.stderr)
    raise SystemExit(1)

print(
    f"contract coverage OK: {len(contract_ids)} contract IDs, "
    f"{len(visible_ids)} test IDs, spec {next(iter(spec_ids))}"
)
