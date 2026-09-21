"""Run main's identity-bound morphology diagnostic on frozen ST052-M.

This is delivery/diagnostic plumbing only.  It joins two already-governed assets:
main's generic morphology wrapper and the constrained branch's authenticated
ST052-M whole-child runtime.  It does not change the velocity field, fit public
images, or promote PDE/visual truth.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import re
from typing import Any

SCHEMA = "st052-identity-bound-morphology-execution/v1"
TASK_ID = "CR-A9-102"
CANDIDATE_ID = "ST052-M-linear-temporal-child-v1"
MAIN_MORPHOLOGY_BASE = "502936d119a0b582e6f271813804581abf9cc9e5"
CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE = "6a293b3c870d70b8d8aece1b12d68685c9b8b010"
SOURCE_HEAD = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
SOURCE_TREE = "cbaf188a7079e17984b0f78d7c8c38df8e11401b"
SOURCE_RUNTIME_IDENTITY_SHA256 = "a0a4682d62677911758e6f9240439c133ee1e3787d3801b87c25de6dcd32f89c"
DEPENDENCY_CONTRACT_REL = Path("configs/st052m_runtime_dependency_acceptance.json")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")

TRUTH_BOUNDARY = {
    "candidate_changed": False,
    "velocity_coefficients_changed": False,
    "pressure_changed": False,
    "forcing_changed": False,
    "scientific_threshold_changed": False,
    "public_source_numeric_targets_used": False,
    "renderer_or_camera_used": False,
    "pixel_loss_used": False,
    "velocity_export_ready": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 64-hex SHA-256 digest")
    return value


def build_velocity_identity_payload(
    candidate_manifest: dict[str, Any],
    dependency_contract: dict[str, Any],
    *,
    dependency_contract_sha256: str,
) -> dict[str, Any]:
    """Bind callable semantics to the accepted exact-source runtime contract."""
    if candidate_manifest.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("unexpected ST052 whole-child candidate id")
    whole_id = _require_sha(candidate_manifest.get("whole_candidate_identity_sha256"), "whole_candidate_identity_sha256")
    runtime = candidate_manifest.get("runtime")
    truth = candidate_manifest.get("truth_boundary")
    if not isinstance(runtime, dict) or not isinstance(truth, dict):
        raise ValueError("malformed ST052 whole-child manifest")
    if runtime.get("exact_source_runtime_identity_sha256") != SOURCE_RUNTIME_IDENTITY_SHA256:
        raise ValueError("whole-child exact source runtime identity drift")
    for key in ("visualization_ready", "pde_validated", "source_correspondence_verified"):
        if truth.get(key) is not False:
            raise ValueError(f"whole-child truth boundary promoted: {key}")

    if dependency_contract.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("runtime-dependency contract candidate mismatch")
    if dependency_contract.get("decision") != "accept_authenticated_exact_source_checkout_as_external_visualization_runtime_dependency":
        raise ValueError("runtime-dependency acceptance decision drift")
    required = dependency_contract.get("required_source_runtime")
    policy = dependency_contract.get("runtime_policy")
    dep_truth = dependency_contract.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (required, policy, dep_truth)):
        raise ValueError("malformed runtime-dependency acceptance contract")
    for key, expected in {
        "source_head": SOURCE_HEAD,
        "source_tree": SOURCE_TREE,
        "source_runtime_identity_sha256": SOURCE_RUNTIME_IDENTITY_SHA256,
    }.items():
        if required.get(key) != expected:
            raise ValueError(f"runtime-dependency source identity drift: {key}")
    if policy.get("exact_source_checkout_dependency_explicitly_accepted") is not True:
        raise ValueError("exact source checkout is not explicitly accepted")
    if policy.get("resolved_dependency_versions_are_part_of_frozen_candidate_identity") is not False:
        raise ValueError("resolved dependency version identity policy drift")
    for key in ("velocity_export_ready", "visualization_ready", "pde_validated"):
        if dep_truth.get(key) is not False:
            raise ValueError(f"runtime-dependency truth boundary promoted: {key}")

    return {
        "schema": "st052-exact-source-callable-velocity-identity/v1",
        "candidate_id": CANDIDATE_ID,
        "whole_candidate_identity_sha256": whole_id,
        "exact_source_runtime": {
            "source_head": SOURCE_HEAD,
            "source_tree": SOURCE_TREE,
            "source_runtime_identity_sha256": SOURCE_RUNTIME_IDENTITY_SHA256,
        },
        "runtime_dependency_acceptance": {
            "constrained_merge_commit": CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE,
            "contract_sha256": _require_sha(dependency_contract_sha256, "dependency_contract_sha256"),
        },
        "callable": "St052LinearTemporalWholeCandidate.velocity(x,y,z,t)->[...,3]",
        "resolved_runtime_versions_identity_bound": False,
    }


def validate_execution_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected ST052 morphology execution schema/task")
    provenance = receipt.get("provenance")
    identity = receipt.get("candidate_identity")
    morphology = receipt.get("morphology_receipt")
    runtime = receipt.get("execution_runtime")
    truth = receipt.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (provenance, identity, morphology, runtime, truth)):
        raise ValueError("malformed ST052 morphology execution receipt")
    if provenance.get("main_morphology_base") != MAIN_MORPHOLOGY_BASE:
        raise ValueError("generic morphology base drift")
    if provenance.get("constrained_runtime_acceptance_merge") != CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE:
        raise ValueError("constrained runtime-acceptance merge drift")
    if provenance.get("source_head") != SOURCE_HEAD or provenance.get("source_tree") != SOURCE_TREE:
        raise ValueError("ST052 source provenance drift")

    if identity.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("candidate identity drift")
    candidate_sha = _require_sha(identity.get("candidate_sha256"), "candidate_sha256")
    payload = identity.get("velocity_identity_payload")
    velocity_sha = _require_sha(identity.get("velocity_identity_sha256"), "velocity_identity_sha256")
    if not isinstance(payload, dict) or _canonical_sha256(payload) != velocity_sha:
        raise ValueError("velocity identity digest mismatch")
    if payload.get("whole_candidate_identity_sha256") != candidate_sha:
        raise ValueError("candidate/velocity identity linkage mismatch")
    if payload.get("resolved_runtime_versions_identity_bound") is not False:
        raise ValueError("resolved runtime versions must remain execution-only")

    if morphology.get("candidate_identity_bound") is not True or morphology.get("cylindrical_morphology_diagnostic_ready") is not True:
        raise ValueError("morphology receipt is not an executed identity-bound diagnostic")
    expected_identity = {
        "candidate_id": CANDIDATE_ID,
        "candidate_sha256": candidate_sha,
        "velocity_identity_sha256": velocity_sha,
    }
    if morphology.get("candidate_identity") != expected_identity:
        raise ValueError("morphology receipt consumed a different candidate identity")
    protocol = morphology.get("protocol")
    if not isinstance(protocol, dict):
        raise ValueError("missing morphology protocol")
    for key in ("source_numeric_targets_used", "renderer_or_camera_used", "pixel_loss_used"):
        if protocol.get(key) is not False:
            raise ValueError(f"morphology protocol promoted forbidden source/render target: {key}")
    for key in ("visualization_ready", "visual_correspondence_verified", "source_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        if morphology.get(key) is not False:
            raise ValueError(f"morphology receipt promoted truth state: {key}")
    if runtime.get("resolved_versions_identity_bound") is not False:
        raise ValueError("execution runtime versions must not silently change candidate identity")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("ST052 morphology execution truth boundary drift")

    if receipt.get("receipt_sha256") is not None:
        unsigned = dict(receipt)
        recorded = unsigned.pop("receipt_sha256")
        if recorded != _canonical_sha256(unsigned):
            raise ValueError("execution receipt checksum mismatch")


def _install_constrained_package_path(constrained_root: Path) -> None:
    import openai_ns_reconstruction as package

    constrained_pkg = (constrained_root / "src" / "openai_ns_reconstruction").resolve()
    if not constrained_pkg.is_dir():
        raise FileNotFoundError(constrained_pkg)
    existing = [str(Path(p).resolve()) for p in package.__path__]
    package.__path__ = [str(constrained_pkg), *[p for p in existing if p != str(constrained_pkg)]]


def _version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def execute(*, constrained_root: str | Path, source_root: str | Path, bundle_dir: str | Path) -> dict[str, Any]:
    constrained_root = Path(constrained_root).resolve()
    source_root = Path(source_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()

    # Pin the diagnostic implementation to main before exposing constrained
    # package modules.  This prevents an older constrained sibling from silently
    # shadowing the already-merged #930/#834 morphology code.
    morphology_mod = importlib.import_module("openai_ns_reconstruction.candidate_cylindrical_morphology")
    _install_constrained_package_path(constrained_root)
    acceptance = importlib.import_module("openai_ns_reconstruction.constrained_st052_runtime_dependency_acceptance")
    whole = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")

    acceptance_receipt = acceptance.audit(constrained_root)
    contract_path = constrained_root / DEPENDENCY_CONTRACT_REL
    dependency_contract = json.loads(contract_path.read_text(encoding="utf-8"))
    dependency_contract_sha = _sha256_file(contract_path)
    candidate = whole.load_bundle_runtime(bundle_dir, exact_source_root=source_root)

    velocity_payload = build_velocity_identity_payload(
        candidate.manifest,
        dependency_contract,
        dependency_contract_sha256=dependency_contract_sha,
    )
    velocity_sha = _canonical_sha256(velocity_payload)
    candidate_sha = str(candidate.identity_sha256)
    morphology = morphology_mod.fingerprint_identified_velocity_field(
        candidate,
        candidate_id=CANDIDATE_ID,
        candidate_sha256=candidate_sha,
        velocity_identity_sha256=velocity_sha,
    )
    morphology_mod.validate_candidate_morphology_receipt(morphology)

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "provenance": {
            "main_morphology_base": MAIN_MORPHOLOGY_BASE,
            "constrained_runtime_acceptance_merge": CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE,
            "source_head": SOURCE_HEAD,
            "source_tree": SOURCE_TREE,
            "source_runtime_identity_sha256": SOURCE_RUNTIME_IDENTITY_SHA256,
            "runtime_dependency_contract_sha256": dependency_contract_sha,
            "runtime_dependency_acceptance_receipt_sha256": _canonical_sha256(acceptance_receipt),
        },
        "candidate_identity": {
            "candidate_id": CANDIDATE_ID,
            "candidate_sha256": candidate_sha,
            "velocity_identity_payload": velocity_payload,
            "velocity_identity_sha256": velocity_sha,
        },
        "execution_runtime": {
            "python": platform.python_version(),
            "numpy": _version("numpy"),
            "scipy": _version("scipy"),
            "sympy": _version("sympy"),
            "resolved_versions_identity_bound": False,
            "scope": "execution receipt only; A8-DELIVERY-01 keeps resolved dependency versions outside frozen candidate identity",
        },
        "morphology_receipt": morphology,
        "direct_contribution": "executes the merged renderer-independent morphology fingerprint on the same authenticated frozen ST052-M save/load velocity identity",
        "remaining_limits": [
            "no official-public numerical morphology target is used",
            "no visual/source correspondence verdict is made",
            "ST052 velocity_export_ready remains false pending the broader same-identity delivery smoke/promotion",
            "no compatible pressure or restricted forcing is rebuilt",
            "no fresh 4096-point complete-NS validation is performed",
        ],
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    receipt["receipt_sha256"] = _canonical_sha256(receipt)
    validate_execution_receipt(receipt)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constrained-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    receipt = execute(constrained_root=args.constrained_root, source_root=args.source_root, bundle_dir=args.bundle_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "candidate_id": receipt["candidate_identity"]["candidate_id"],
        "candidate_sha256": receipt["candidate_identity"]["candidate_sha256"],
        "velocity_identity_sha256": receipt["candidate_identity"]["velocity_identity_sha256"],
        "morphology_measurement_sha256": receipt["morphology_receipt"]["measurement_sha256"],
        "receipt_sha256": receipt["receipt_sha256"],
        "visualization_ready": receipt["truth_boundary"]["visualization_ready"],
        "pde_validated": receipt["truth_boundary"]["pde_validated"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
