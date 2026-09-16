"""Build a truth-bounded bill of materials for the constrained velocity artifact.

This module inventories the files needed to rerun and inspect the currently frozen
velocity candidate.  It intentionally does not recompute acceptance thresholds or
promote numerical/visual claims.  In particular, a reproducible callable field may
be export-ready while PDE validation and visual correspondence remain false.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "constrained_delivery_capsule_v1"
DEFAULT_ROOT = Path(__file__).resolve().parents[2]

CORE_PATHS = {
    "selected_candidate": "artifacts/constrained/coupled_joint/candidate.json",
    "packaged_candidate": "src/openai_ns_reconstruction/data/velocity_candidate.json",
    "velocity_evaluator": "src/openai_ns_reconstruction/velocity_components.py",
    "velocity_api_doc": "docs/VELOCITY_API.md",
    "constraint_config": "configs/constraints_coupled.json",
    "training_record": "artifacts/constrained/coupled_joint/training.json",
    "validation_record": "artifacts/constrained/coupled_joint/validation.json",
    "structure_report": "artifacts/constrained/structure_identities.json",
    "structure_doc": "docs/STRUCTURE_IDENTITIES.md",
    "reference_grid": "artifacts/visual/velocity_api/grid.npz",
    "reference_grid_metadata": "artifacts/visual/velocity_api/metadata.json",
    "project_status": "project_status.json",
}

OPTIONAL_DELIVERY_PATHS = {
    "matlab_export_entry": "src/openai_ns_reconstruction/constrained_matlab_export.py",
    "vtk_export_entry": "src/openai_ns_reconstruction/constrained_vtk_export.py",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON artifact: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON artifact must contain an object: {path}")
    return value


def _inventory_entry(root: Path, relative: str, *, required: bool) -> dict[str, Any]:
    path = root / relative
    if not path.is_file():
        return {
            "path": relative,
            "required": bool(required),
            "status": "pending",
            "sha256": None,
        }
    return {
        "path": relative,
        "required": bool(required),
        "status": "present",
        "sha256": _sha256(path),
    }


def build_delivery_capsule(root: str | Path = DEFAULT_ROOT) -> dict[str, Any]:
    """Build a deterministic inventory for the currently selected velocity artifact.

    The selected candidate and packaged candidate must be byte-identical.  Training
    and validation seeds must remain distinct, and the validation record must point
    at the selected candidate path.  These are artifact-consistency checks only; a
    passing capsule is not a PDE or visual-correspondence acceptance result.
    """
    root = Path(root).resolve()
    missing_core = [name for name, rel in CORE_PATHS.items() if not (root / rel).is_file()]
    if missing_core:
        raise ValueError(f"missing required delivery inputs: {missing_core}")

    candidate = _read_json(root / CORE_PATHS["selected_candidate"])
    packaged = _read_json(root / CORE_PATHS["packaged_candidate"])
    config = _read_json(root / CORE_PATHS["constraint_config"])
    training = _read_json(root / CORE_PATHS["training_record"])
    validation = _read_json(root / CORE_PATHS["validation_record"])
    structure = _read_json(root / CORE_PATHS["structure_report"])
    project_status = _read_json(root / CORE_PATHS["project_status"])

    selected_sha = _sha256(root / CORE_PATHS["selected_candidate"])
    packaged_sha = _sha256(root / CORE_PATHS["packaged_candidate"])
    if selected_sha != packaged_sha:
        raise ValueError("selected candidate bytes differ from packaged public candidate")

    selected_family = candidate.get("family")
    if not isinstance(selected_family, str) or not selected_family:
        raise ValueError("selected candidate requires a nonempty family")
    if packaged.get("family") != selected_family:
        raise ValueError("packaged candidate family differs from selected candidate family")

    training_seed = training.get("seed")
    validation_seed = validation.get("seed")
    if not isinstance(training_seed, int) or not isinstance(validation_seed, int):
        raise ValueError("training and validation records require integer seeds")
    if training_seed == validation_seed:
        raise ValueError("training and validation seeds must be distinct")

    validation_candidate = validation.get("candidate")
    if validation_candidate != CORE_PATHS["selected_candidate"]:
        raise ValueError("validation record does not point at the selected candidate")

    domain = config.get("domain")
    if not isinstance(domain, dict):
        raise ValueError("constraint config requires a domain object")
    evaluation_box = domain.get("evaluation_box")
    time_interval = domain.get("time_interval")
    if not (
        isinstance(evaluation_box, list)
        and len(evaluation_box) == 3
        and all(isinstance(axis, list) and len(axis) == 2 for axis in evaluation_box)
    ):
        raise ValueError("constraint config has malformed evaluation_box")
    if not (isinstance(time_interval, list) and len(time_interval) == 2):
        raise ValueError("constraint config has malformed time_interval")

    inventory = {
        name: _inventory_entry(root, rel, required=True)
        for name, rel in CORE_PATHS.items()
    }
    inventory.update(
        {
            name: _inventory_entry(root, rel, required=False)
            for name, rel in OPTIONAL_DELIVERY_PATHS.items()
        }
    )

    optional_present = {
        name: entry["status"] == "present"
        for name, entry in inventory.items()
        if name in OPTIONAL_DELIVERY_PATHS
    }
    all_core_present = all(
        entry["status"] == "present"
        for name, entry in inventory.items()
        if name in CORE_PATHS
    )

    reference_times = [
        float(time_interval[0]),
        0.5 * (float(time_interval[0]) + float(time_interval[1])),
        float(time_interval[1]),
    ]

    validation_status = validation.get("status")
    visual_correspondence = project_status.get("visual_correspondence")

    capsule = {
        "schema": SCHEMA,
        "scope": "reproducible_velocity_artifact_inventory_only",
        "candidate": {
            "family": selected_family,
            "selected_path": CORE_PATHS["selected_candidate"],
            "packaged_path": CORE_PATHS["packaged_candidate"],
            "sha256": selected_sha,
            "status": candidate.get("status", "candidate"),
        },
        "reproducibility": {
            "config_path": CORE_PATHS["constraint_config"],
            "training_path": CORE_PATHS["training_record"],
            "validation_path": CORE_PATHS["validation_record"],
            "training_seed": training_seed,
            "validation_seed": validation_seed,
            "training_calls": training.get("calls"),
            "training_budget": training.get("max_nfev"),
            "validation_points": validation.get("points"),
            "validation_status": validation_status,
            "structure_status": structure.get("status"),
        },
        "recommended_evaluation": {
            "units": config.get("units"),
            "evaluation_box": evaluation_box,
            "time_interval": time_interval,
            "reference_times": reference_times,
            "reference_times_origin": "autonomous_design",
        },
        "entry_points": {
            "python_components": "from openai_ns_reconstruction.velocity_components import velocity",
            "python_grid_export": "python -m openai_ns_reconstruction.velocity_components --export <output-dir>",
            "matlab_export_module": OPTIONAL_DELIVERY_PATHS["matlab_export_entry"]
            if optional_present["matlab_export_entry"]
            else None,
            "vtk_export_module": OPTIONAL_DELIVERY_PATHS["vtk_export_entry"]
            if optional_present["vtk_export_entry"]
            else None,
        },
        "inventory": inventory,
        "claim_states": {
            "velocity_export_ready": bool(all_core_present),
            "visualization_ready": False,
            "visual_correspondence_verified": visual_correspondence == "verified",
            "pde_validated": False,
            "paper_exact": False,
        },
        "current_validation_label": validation_status,
        "final_artifact_complete": bool(
            all_core_present
            and optional_present["matlab_export_entry"]
            and optional_present["vtk_export_entry"]
            and visual_correspondence == "verified"
        ),
        "limitations": [
            "inventory/hash consistency is not Navier-Stokes validation",
            "saved training convergence is not independent validation",
            "visualization/export availability is independent of PDE acceptance",
            "no exact OpenAI field or blow-up theorem is claimed",
        ],
    }
    return validate_delivery_capsule(capsule)


def validate_delivery_capsule(capsule: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(capsule, dict) or capsule.get("schema") != SCHEMA:
        raise ValueError(f"capsule schema must be {SCHEMA!r}")
    states = capsule.get("claim_states")
    required_states = {
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
    }
    if not isinstance(states, dict) or set(states) != required_states:
        raise ValueError("capsule must contain exactly the independent delivery claim states")
    if not all(isinstance(states[name], bool) for name in required_states):
        raise ValueError("delivery claim states must be booleans")
    if states["paper_exact"]:
        raise ValueError("paper_exact cannot be promoted by a reproducibility capsule")
    if states["pde_validated"]:
        raise ValueError("pde_validated cannot be promoted by a reproducibility capsule")

    inventory = capsule.get("inventory")
    if not isinstance(inventory, dict) or not inventory:
        raise ValueError("capsule inventory must be nonempty")
    for name, entry in inventory.items():
        if not isinstance(entry, dict):
            raise ValueError(f"inventory entry {name!r} must be an object")
        if entry.get("status") not in {"present", "pending"}:
            raise ValueError(f"inventory entry {name!r} has invalid status")
        digest = entry.get("sha256")
        if entry["status"] == "present":
            if not isinstance(digest, str) or len(digest) != 64:
                raise ValueError(f"inventory entry {name!r} lacks a SHA256 digest")
        elif digest is not None:
            raise ValueError(f"pending inventory entry {name!r} must not carry a digest")

    if capsule.get("final_artifact_complete"):
        pending = [
            name for name, entry in inventory.items()
            if entry.get("required") and entry.get("status") != "present"
        ]
        if pending:
            raise ValueError(f"complete artifact cannot have pending required entries: {pending}")
        if not states["velocity_export_ready"] or not states["visualization_ready"]:
            raise ValueError("complete artifact requires export and visualization readiness")

    return capsule


def write_delivery_capsule(output: str | Path, root: str | Path = DEFAULT_ROOT) -> Path:
    output = Path(output)
    capsule = build_delivery_capsule(root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(capsule, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/constrained/delivery_capsule.json"),
    )
    args = parser.parse_args()
    path = write_delivery_capsule(args.output, args.root)
    print(path)


if __name__ == "__main__":
    main()
