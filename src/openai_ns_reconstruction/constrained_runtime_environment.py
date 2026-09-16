"""Capture resolved runtime dependencies for a frozen public velocity artifact.

This receipt is reproducibility metadata only.  It binds the exact public velocity
candidate to the Python/runtime package versions that consumed it.  A successful
receipt does not establish visualization correspondence, Navier--Stokes acceptance,
paper exactness, or a blow-up result.
"""
from __future__ import annotations

import argparse
import hashlib
from importlib import metadata as importlib_metadata
import json
from pathlib import Path
import platform
import sys
from typing import Any, Iterable


SCHEMA = "constrained_runtime_environment_v1"
PROJECT_DISTRIBUTION = "openai-ns-velocity-field-reconstruction"
DEFAULT_REQUIRED = (PROJECT_DISTRIBUTION, "numpy", "scipy")
DEFAULT_OPTIONAL = ("matplotlib", "pytest", "sympy")


def _canonical_sha256(payload: dict[str, Any]) -> str:
    raw = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _resolved_versions(names: Iterable[str], *, required: bool) -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for raw_name in names:
        name = str(raw_name).strip()
        if not name:
            raise ValueError("distribution names must be nonempty")
        try:
            version = importlib_metadata.version(name)
        except importlib_metadata.PackageNotFoundError as exc:
            if required:
                raise ValueError(f"required distribution is not installed: {name}") from exc
            version = None
        result[name] = version
    return dict(sorted(result.items()))


def _require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError(f"{label} must be a 64-character SHA256 digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ValueError(f"{label} must be hexadecimal") from exc
    return value.lower()


def build_runtime_environment_receipt(
    field=None,
    *,
    required_distributions: Iterable[str] = DEFAULT_REQUIRED,
    optional_distributions: Iterable[str] = DEFAULT_OPTIONAL,
) -> dict[str, Any]:
    """Return a deterministic environment receipt for one public velocity field.

    ``field`` must expose the same narrow contract used by the final delivery:
    ``sha256`` plus ``metadata()``.  When omitted, the checked-in ``VelocityField``
    is loaded, so the receipt is bound to the candidate users actually evaluate.
    """
    if field is None:
        from .velocity_components import VelocityField

        field = VelocityField()

    candidate_sha256 = _require_sha256(getattr(field, "sha256", None), "field.sha256")
    field_metadata = field.metadata()
    if not isinstance(field_metadata, dict):
        raise ValueError("public velocity metadata must be a mapping")
    if field_metadata.get("candidate_sha256") != candidate_sha256:
        raise ValueError("public velocity metadata candidate SHA256 mismatch")
    family = field_metadata.get("family")
    if not isinstance(family, str) or not family:
        raise ValueError("public velocity metadata requires a nonempty family")

    required = _resolved_versions(required_distributions, required=True)
    optional = _resolved_versions(optional_distributions, required=False)
    if PROJECT_DISTRIBUTION not in required:
        raise ValueError("required distributions must include the project distribution")

    declared = importlib_metadata.requires(PROJECT_DISTRIBUTION) or []
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "scope": "runtime_reproducibility_metadata_only",
        "candidate": {
            "family": family,
            "sha256": candidate_sha256,
            "public_evaluator": f"{type(field).__module__}.{type(field).__qualname__}.at_points",
        },
        "runtime": {
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "python_cache_tag": getattr(sys.implementation, "cache_tag", None),
            "sys_platform": sys.platform,
            "platform": platform.platform(),
        },
        "dependencies": {
            "required_resolved_versions": required,
            "optional_resolved_versions": optional,
            "project_declared_requirements": sorted(str(item) for item in declared),
        },
        "delivery_state": {
            "velocity_export_ready": "not_assessed_here",
            "visualization_ready": "not_assessed_here",
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
        "limitations": [
            "environment replayability is not Navier-Stokes validation",
            "resolved dependency versions do not verify visual correspondence",
            "no exact OpenAI field, singularity, or blow-up theorem is claimed",
        ],
    }
    payload["receipt_sha256"] = _canonical_sha256(payload)
    return validate_runtime_environment_receipt(payload)


def validate_runtime_environment_receipt(receipt: dict[str, Any]) -> dict[str, Any]:
    """Fail-close validate a runtime receipt without promoting scientific claims."""
    if not isinstance(receipt, dict) or receipt.get("schema") != SCHEMA:
        raise ValueError(f"receipt schema must be {SCHEMA!r}")

    candidate = receipt.get("candidate")
    if not isinstance(candidate, dict):
        raise ValueError("candidate receipt section must be an object")
    _require_sha256(candidate.get("sha256"), "candidate.sha256")
    if not isinstance(candidate.get("family"), str) or not candidate["family"]:
        raise ValueError("candidate family must be nonempty")
    if not isinstance(candidate.get("public_evaluator"), str) or not candidate["public_evaluator"]:
        raise ValueError("public evaluator identity must be nonempty")

    runtime = receipt.get("runtime")
    if not isinstance(runtime, dict) or not runtime.get("python_version"):
        raise ValueError("runtime metadata is incomplete")

    dependencies = receipt.get("dependencies")
    if not isinstance(dependencies, dict):
        raise ValueError("dependencies section must be an object")
    required = dependencies.get("required_resolved_versions")
    if not isinstance(required, dict) or PROJECT_DISTRIBUTION not in required:
        raise ValueError("required resolved dependencies must include the project")
    if not all(isinstance(value, str) and value for value in required.values()):
        raise ValueError("required dependency versions must be resolved nonempty strings")
    optional = dependencies.get("optional_resolved_versions")
    if not isinstance(optional, dict) or not all(
        value is None or (isinstance(value, str) and value) for value in optional.values()
    ):
        raise ValueError("optional dependency versions must be strings or null")
    declared = dependencies.get("project_declared_requirements")
    if not isinstance(declared, list) or not all(isinstance(item, str) for item in declared):
        raise ValueError("declared project requirements must be a string list")

    state = receipt.get("delivery_state")
    expected_state = {
        "velocity_export_ready": "not_assessed_here",
        "visualization_ready": "not_assessed_here",
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    if state != expected_state:
        raise ValueError("runtime receipt cannot promote independent delivery/scientific states")

    claimed_digest = _require_sha256(receipt.get("receipt_sha256"), "receipt_sha256")
    payload = dict(receipt)
    payload.pop("receipt_sha256", None)
    actual_digest = _canonical_sha256(payload)
    if claimed_digest != actual_digest:
        raise ValueError("runtime receipt digest mismatch")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = build_runtime_environment_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
