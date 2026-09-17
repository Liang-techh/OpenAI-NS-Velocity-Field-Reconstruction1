"""Checked delivery inventory for the support-connected Eq. (4.5) field.

The capsule is deliberately about *what can be consumed reproducibly* from the
current support-connected velocity ancestry.  It does not import scientific
results from sibling PRs and it never upgrades visualization or PDE claims from
callability, serialization, or green CI.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .eq45_supported_delivery import default_field


SCHEMA = "eq45_supported_delivery_capsule_v1"
EXPECTED_CHILD_SHA256 = "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
EXPECTED_PARENT_SHA256 = "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"


def _load_constraints(repo_root: Path) -> dict[str, Any]:
    path = repo_root / "configs" / "constraints.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise ValueError("unsupported constrained configuration schema")
    return data


def _integrated(repo_root: Path, relative_path: str) -> bool:
    return (repo_root / relative_path).is_file()


def build_eq45_supported_delivery_capsule(repo_root: str | Path) -> dict[str, Any]:
    """Build the current support-connected delivery bill of materials.

    Only evidence present in this ancestry is consumed.  Open/sibling PR
    results remain explicitly unconsumed so that a delivery snapshot cannot
    launder them into accepted scientific state.
    """

    root = Path(repo_root)
    constraints = _load_constraints(root)
    field = default_field()
    metadata = field.metadata()

    if field.sha256 != EXPECTED_CHILD_SHA256:
        raise ValueError("default supported child identity drifted")
    if metadata["parent_sha256"] != EXPECTED_PARENT_SHA256:
        raise ValueError("default supported parent identity drifted")

    truth = dict(metadata["truth_boundary"])
    required_truth = {
        "callable_serializable": True,
        "velocity_export_ready": True,
        "physical_support_connection_implemented": True,
        "physical_support_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    for key, value in required_truth.items():
        if truth.get(key) is not value:
            raise ValueError(f"supported delivery truth boundary drifted at {key}")

    domain = constraints["domain"]
    validation = constraints["validation"]
    t0, t1 = (float(value) for value in domain["time_interval"])
    reference_times = [t0, 0.5 * (t0 + t1), t1]

    return {
        "schema": SCHEMA,
        "candidate": {
            "family": metadata["family"],
            "child_sha256": field.sha256,
            "parent_sha256": metadata["parent_sha256"],
            "support_connected": True,
        },
        "public_interface": {
            "velocity": metadata["entrypoint"],
            "point_batch": "openai_ns_reconstruction.eq45_supported_delivery:default_field().at_points",
            "grid": "openai_ns_reconstruction.eq45_supported_delivery:default_field().grid",
            "save_candidate": "openai_ns_reconstruction.eq45_supported_delivery:default_field().save_candidate",
            "load_candidate": "openai_ns_reconstruction.eq45_supported_delivery:Eq45SupportedDeliveryField.load_candidate",
            "components": list(metadata["components"]),
            "grid_layout": list(metadata["grid_layout"]),
        },
        "recommended_window": {
            "evaluation_box": domain["evaluation_box"],
            "physical_support": domain["support"],
            "time_interval": [t0, t1],
            "reference_times": reference_times,
            "reference_times_classification": "autonomous_delivery_choice",
        },
        "preregistered_validation_contract": {
            "experiment_id": constraints["experiment_id"],
            "validation_seed": int(validation["seed"]),
            "held_out_points": int(validation["held_out_points"]),
            "derivative_steps": [float(value) for value in validation["derivative_steps"]],
            "thresholds": dict(validation["thresholds"]),
        },
        "integrated_assets": {
            "named_python_velocity_api": True,
            "direct_grid_evaluator": True,
            "candidate_save_load_api": True,
            "standalone_supported_candidate_json": _integrated(
                root, "artifacts/constrained/eq45_supported_velocity_candidate.json"
            ),
            "matlab_mat_exporter": _integrated(
                root, "src/openai_ns_reconstruction/constrained_matlab_export.py"
            ),
            "vtk_exporter": _integrated(
                root, "src/openai_ns_reconstruction/constrained_vtk_export.py"
            ),
            "python_slice_renderer": _integrated(
                root, "src/openai_ns_reconstruction/constrained_velocity_slice_smoke.py"
            ),
        },
        "scientific_evidence_in_this_ancestry": {
            "physical_support_connection": "implemented_not_independently_validated",
            "energy": "pending_unconsumed_sibling_evidence",
            "divergence": "pending_unconsumed_sibling_evidence",
            "pde": "pending_unconsumed_sibling_evidence",
            "public_visual_correspondence": "pending",
        },
        "states": required_truth,
        "scope": {
            "velocity_changed": False,
            "candidate_promoted": False,
            "claim": "reproducible_delivery_inventory_only",
        },
    }


def validate_eq45_supported_delivery_capsule(
    payload: Mapping[str, Any], repo_root: str | Path
) -> dict[str, Any]:
    """Fail closed unless ``payload`` exactly matches the current governed build."""

    if not isinstance(payload, Mapping):
        raise ValueError("delivery capsule must be an object")
    expected = build_eq45_supported_delivery_capsule(repo_root)
    if dict(payload) != expected:
        raise ValueError("supported Eq45 delivery capsule does not match governed state")
    return expected


def load_checked_eq45_supported_delivery_capsule(
    path: str | Path, repo_root: str | Path
) -> dict[str, Any]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return validate_eq45_supported_delivery_capsule(data, repo_root)


def write_eq45_supported_delivery_capsule(path: str | Path, repo_root: str | Path) -> None:
    payload = build_eq45_supported_delivery_capsule(repo_root)
    Path(path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def public_velocity_probe() -> np.ndarray:
    """Small public-interface smoke used by reproducibility tests."""

    values = default_field().velocity(0.37, -0.29, 0.41, 0.5)
    values = np.asarray(values, dtype=float)
    if values.shape != (3,) or not np.all(np.isfinite(values)):
        raise RuntimeError("supported public velocity probe is malformed")
    return values
