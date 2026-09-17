"""Truth-preserving consistency checks for constrained reconstruction artifacts.

A passing manifest check means the recorded files and metadata agree. It does not
mean the candidate satisfies the Navier--Stokes acceptance thresholds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


_ALLOWED_CONFIG_BINDINGS = {"exact", "inherited_not_exact_experiment_config"}
_ALLOWED_CONSTRAINT_STATUS = {"pass", "fail", "pending", "not_applicable"}


def _json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return data


def git_blob_sha1(path: str | Path) -> str:
    path = Path(path)
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _resolve_inside(root: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise ValueError("artifact path must be a nonempty string")
    root = root.resolve()
    path = (root / relative).resolve()
    if path != root and root not in path.parents:
        raise ValueError(f"artifact path escapes repository root: {relative}")
    return path


def _close(a: Any, b: Any) -> bool:
    try:
        return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-15)
    except (TypeError, ValueError):
        return False


def validate_manifest(manifest_path: str | Path, repo_root: str | Path = ".") -> list[str]:
    manifest_path = Path(manifest_path)
    root = Path(repo_root)
    errors: list[str] = []
    try:
        manifest = _json(manifest_path)
    except Exception as exc:
        return [f"manifest unreadable: {exc}"]

    if manifest.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    if manifest.get("artifact_kind") != "constrained_reconstruction_reproducibility_manifest":
        errors.append("unexpected artifact_kind")
    truth = manifest.get("truth_boundary", {})
    if truth.get("paper_exact") is not False:
        errors.append("truth_boundary.paper_exact must be false")
    if truth.get("full_blowup_proof") is not False:
        errors.append("truth_boundary.full_blowup_proof must be false")

    binding = manifest.get("config_binding", {})
    if binding.get("status") not in _ALLOWED_CONFIG_BINDINGS:
        errors.append("unsupported config_binding.status")

    files = manifest.get("files")
    if not isinstance(files, dict):
        return errors + ["files must be an object"]

    loaded: dict[str, dict[str, Any]] = {}
    for role in ("candidate", "training", "validation", "config"):
        spec = files.get(role)
        if not isinstance(spec, dict):
            errors.append(f"files.{role} missing")
            continue
        try:
            path = _resolve_inside(root, spec.get("path"))
        except Exception as exc:
            errors.append(str(exc))
            continue
        if not path.is_file():
            errors.append(f"files.{role} does not exist: {spec.get('path')}")
            continue
        expected = spec.get("git_blob_sha1")
        actual = git_blob_sha1(path)
        if expected != actual:
            errors.append(f"files.{role} git blob hash mismatch: expected {expected}, actual {actual}")
        try:
            loaded[role] = _json(path)
        except Exception as exc:
            errors.append(f"files.{role} invalid JSON: {exc}")

    if set(("candidate", "training", "validation", "config")) - loaded.keys():
        return errors

    candidate = loaded["candidate"]
    training = loaded["training"]
    validation = loaded["validation"]
    config = loaded["config"]

    if candidate.get("status") != "candidate":
        errors.append("candidate status must remain 'candidate'")
    if manifest.get("candidate_family") != candidate.get("family"):
        errors.append("candidate_family does not match candidate artifact")
    if manifest.get("validation_status") != validation.get("status"):
        errors.append("validation_status does not match validation artifact")

    seeds = manifest.get("seeds", {})
    train_seed = training.get("seed")
    val_seed = validation.get("seed")
    if seeds.get("training") != train_seed:
        errors.append("training seed mismatch")
    if seeds.get("validation") != val_seed:
        errors.append("validation seed mismatch")
    if train_seed == val_seed:
        errors.append("training and validation seeds must differ")

    cfg_opt_seed = config.get("optimization", {}).get("seed")
    cfg_val_seed = config.get("validation", {}).get("seed")
    if cfg_opt_seed != train_seed:
        errors.append("config optimization seed does not match training seed")
    if cfg_val_seed != val_seed:
        errors.append("config validation seed does not match validation seed")

    validation_candidate = validation.get("candidate")
    candidate_path = files["candidate"].get("path")
    if validation_candidate != candidate_path:
        errors.append("validation candidate path does not match manifest candidate path")

    rows = validation.get("rows")
    if not isinstance(rows, list) or not rows:
        errors.append("validation rows missing")
    else:
        try:
            finest = min(float(row["step"]) for row in rows)
            finest_rows = [row for row in rows if _close(row.get("step"), finest)]
            pde_max = max(float(row["residual_sampled_max"]) for row in finest_rows)
            div_max = max(float(row["divergence_sampled_max"]) for row in finest_rows)
            observed = manifest.get("observed_finest_validation", {})
            if not _close(observed.get("step"), finest):
                errors.append("recorded finest validation step mismatch")
            if not _close(observed.get("pde_residual_sampled_max"), pde_max):
                errors.append("recorded PDE sampled maximum mismatch")
            if not _close(observed.get("divergence_sampled_max"), div_max):
                errors.append("recorded divergence sampled maximum mismatch")

            thresholds = config.get("validation", {}).get("thresholds", {})
            expected_status = {
                "pde_residual_max": "pass" if pde_max <= float(thresholds["pde_residual_max"]) else "fail",
                "divergence_max": "pass" if div_max <= float(thresholds["divergence_max"]) else "fail",
            }
            snapshot = manifest.get("constraint_snapshot", {})
            for name, status in expected_status.items():
                entry = snapshot.get(name, {})
                if entry.get("status") != status:
                    errors.append(f"constraint_snapshot.{name}.status should be {status}")
        except (KeyError, TypeError, ValueError) as exc:
            errors.append(f"could not derive finest validation metrics: {exc}")

    snapshot = manifest.get("constraint_snapshot", {})
    statuses = []
    if not isinstance(snapshot, dict) or not snapshot:
        errors.append("constraint_snapshot must be a nonempty object")
    else:
        for name, entry in snapshot.items():
            if not isinstance(entry, dict) or entry.get("status") not in _ALLOWED_CONSTRAINT_STATUS:
                errors.append(f"invalid constraint status for {name}")
            else:
                statuses.append(entry["status"])

    acceptance_ready = manifest.get("acceptance_ready")
    if not isinstance(acceptance_ready, bool):
        errors.append("acceptance_ready must be boolean")
    elif acceptance_ready:
        if binding.get("status") != "exact":
            errors.append("acceptance_ready requires an exact experiment config binding")
        if validation.get("status") != "passed_validation":
            errors.append("acceptance_ready requires passed_validation")
        if any(status != "pass" for status in statuses):
            errors.append("acceptance_ready requires every recorded constraint to pass")

    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest")
    parser.add_argument("--repo-root", default=".")
    args = parser.parse_args(argv)
    errors = validate_manifest(args.manifest, args.repo_root)
    result = {
        "status": "manifest_consistency_pass" if not errors else "manifest_consistency_fail",
        "acceptance_claimed": False,
        "errors": errors,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
