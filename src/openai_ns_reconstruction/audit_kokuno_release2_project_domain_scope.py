"""Fail-closed CR002 audit for Kokuno release2 project-domain delivery scope.

This module does not alter the velocity field.  It binds the exact A1 #1204
stage representation to the canonical constrained-domain contract and prevents
stage-local callability/save-load evidence from being promoted into global
Kokuno delivery, visual correspondence, PDE validation, or exact-field claims.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

from .kokuno_pa16_current_cartesian_postswirl_release2 import (
    KokunoPA16CurrentCartesianPostSwirlRelease2,
)

SCHEMA = "cr002-kokuno-release2-project-domain-scope-v1"
TASK_ID = "CR002-KOKUNO-RELEASE2-PROJECT-DOMAIN-SCOPE-146"
CONTRACT_REL = Path("configs/kokuno_release2_project_domain_scope.json")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _canonical_constraints_checks(cfg: Mapping[str, Any], expected: Mapping[str, Any]) -> None:
    _require(cfg["nu"] == expected["nu"] == 0.01, "canonical nu drift")
    domain = cfg["domain"]
    _require(domain["physical"] == expected["physical_domain"] == "R^3", "physical domain drift")
    _require(domain["evaluation_box"] == expected["evaluation_box"] == [[-2, 2], [-2, 2], [-2, 2]], "evaluation box drift")
    _require(domain["support"] == expected["support"] == "r < 2 and abs(z) < 2", "support drift")
    _require(domain["time_interval"] == expected["time_interval"] == [0.25, 0.75], "time interval drift")

    forcing = cfg["forcing"]
    _require(forcing["mode"] == expected["forcing_mode"] == "restricted_two_parameter_family", "forcing family drift")
    _require("No residual-dependent basis or pointwise free force" in forcing["restriction"], "free-force firewall drift")

    nontrivial = cfg["nontriviality"]
    _require(nontrivial["reference_energy"] == expected["reference_energy"] == 1.0, "reference energy drift")
    _require(
        nontrivial["reference_energy_abs_tolerance"]
        == expected["reference_energy_abs_tolerance"]
        == 0.001,
        "reference energy tolerance drift",
    )
    _require("reject collapsed candidates" in nontrivial["enforcement"], "collapse firewall drift")

    val = cfg["validation"]
    _require(val["seed"] == expected["validation_seed"] == 914027, "validation seed drift")
    _require(val["held_out_points"] == expected["held_out_points"] == 4096, "held-out size drift")
    _require(val["derivative_steps"] == expected["derivative_steps"] == [0.02, 0.01, 0.005], "FD ladder drift")
    _require(
        val["quadrature_orders_per_axis"]
        == expected["quadrature_orders_per_axis"]
        == [24, 48, 96],
        "quadrature ladder drift",
    )
    th = val["thresholds"]
    _require(th["pde_residual_max"] == expected["momentum_max"] == 0.001, "momentum max gate drift")
    _require(th["pde_residual_L2"] == expected["momentum_l2"] == 0.001, "momentum L2 gate drift")
    _require(th["divergence_max"] == expected["divergence_max"] == 1e-5, "divergence max gate drift")
    _require(th["divergence_L2"] == expected["divergence_l2"] == 1e-5, "divergence L2 gate drift")


def audit_contract(
    payload: Mapping[str, Any] | None = None,
    *,
    root: Path | None = None,
    run_runtime_probe: bool = True,
) -> dict[str, Any]:
    """Replay the governance contract and return a scoped receipt.

    The runtime probe is deliberately a canonical CR001 point, not a source
    claim.  Its rejection proves only that the exact #1204 object is a
    stage-local evaluator rather than a total evaluator for the project box.
    """

    repo = Path(root) if root is not None else _repo_root()
    contract = dict(payload) if payload is not None else _load_json(repo / CONTRACT_REL)

    _require(contract["schema"] == SCHEMA, "schema drift")
    _require(contract["task_id"] == TASK_ID, "task id drift")
    _require(
        set(contract["provenance_classes"])
        == {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"},
        "four-way provenance partition drift",
    )

    lineage = contract["audited_lineage"]
    _require(lineage["a1_pr"] == 1204, "A1 PR drift")
    _require(
        lineage["exact_head"] == "18e75e6e2df43206147db34788c68ee6a7fa0f37",
        "A1 exact-head drift",
    )
    impl = repo / lineage["implementation_path"]
    _require(
        _git_blob_sha1(impl) == lineage["implementation_blob_sha1"] == "101d8da6b3deff2322e5c2efaea278fefa7d19e3",
        "release2 implementation identity drift",
    )
    source_text = impl.read_text()
    _require(
        "post-swirl release2 requires 0<=log(X/X_release2_start)<=1" in source_text,
        "release2 source-chart hard guard disappeared",
    )
    _require('"unified_global_cartesian_velocity_export_ready": False' in source_text, "global-export truth guard drift")
    _require('"outer_global_leading_velocity_materialized": False' in source_text, "global-leading truth guard drift")

    canonical = contract["canonical_contract"]
    constraints_path = repo / canonical["constraints_path"]
    status_path = repo / canonical["project_status_path"]
    _require(
        _git_blob_sha1(constraints_path)
        == canonical["constraints_blob_sha1"]
        == "6c559e42895a606e2ef025ade4cb448966d75814",
        "canonical constraints blob drift",
    )
    _require(
        _git_blob_sha1(status_path)
        == canonical["project_status_blob_sha1"]
        == "f8e05e9d25ab6cc83a7221fc431ab42319d38ff0",
        "local canonical status snapshot drift",
    )
    constraints = _load_json(constraints_path)
    _canonical_constraints_checks(constraints, canonical)

    status = _load_json(status_path)
    _require(status["states"]["velocity_export_ready"] is True, "canonical Eq45 export unexpectedly downgraded")
    for key in ("visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        _require(status["states"][key] is False, f"canonical {key} unexpectedly promoted")

    stage = contract["audited_stage_contract"]
    _require(stage["stage_callable_save_load_materialized"] is True, "release2 stage-callable fact drift")
    for key in (
        "source_terminal_multiplier_materialized",
        "source_exterior_heat_replacement_materialized",
        "outer_global_leading_velocity_materialized",
        "unified_global_cartesian_velocity_export_ready",
    ):
        _require(stage[key] is False, f"release2 stage boundary promoted: {key}")

    firewall = contract["transfer_firewall"]
    _require(all(value is False for value in firewall.values()), "stage/global evidence-transfer firewall opened")

    states = contract["independent_states"]
    _require(states["canonical_eq45_velocity_export_ready"] is True, "canonical export state drift")
    _require(states["kokuno_release2_stage_callable"] is True, "Kokuno stage-callable state drift")
    for key in (
        "kokuno_velocity_export_ready",
        "kokuno_visual_correspondence_verified",
        "kokuno_pde_validated",
        "kokuno_paper_exact",
        "kokuno_openai_field_identified",
    ):
        _require(states[key] is False, f"Kokuno state unexpectedly promoted: {key}")

    witness = contract["canonical_probe_witness"]
    x, y, z, t = (float(v) for v in witness["point"])
    box = canonical["evaluation_box"]
    _require(all(lo <= q <= hi for q, (lo, hi) in zip((x, y, z), box)), "witness left canonical box")
    _require(canonical["time_interval"][0] <= t <= canonical["time_interval"][1], "witness left time interval")
    _require(math.hypot(x, y) < 2.0 and abs(z) < 2.0, "witness left declared support")
    _require(witness["evidence_role"] == "autonomous_representation_mechanics_only", "witness role drift")
    _require(
        witness["expected_release2_stage_result"] == "reject_outside_release2_source_chart",
        "witness expectation drift",
    )

    rejected = None
    error_text = None
    semantic_sha = None
    if run_runtime_probe:
        field = KokunoPA16CurrentCartesianPostSwirlRelease2()
        semantic_sha = field.semantic_sha256
        truth = field.truth_boundary
        _require(truth["unified_global_cartesian_velocity_export_ready"] is False, "runtime global-export truth promoted")
        _require(truth["outer_global_leading_velocity_materialized"] is False, "runtime global-leading truth promoted")
        try:
            field.velocity(x, y, z, t)
        except ValueError as exc:
            error_text = str(exc)
            rejected = "post-swirl release2 requires" in error_text
        else:
            rejected = False
        _require(rejected is True, "canonical in-support probe was accepted as release2-stage total velocity")

    return {
        "schema": "cr002-kokuno-release2-project-domain-scope-receipt-v1",
        "task_id": TASK_ID,
        "a1_exact_head": lineage["exact_head"],
        "release2_semantic_sha256": semantic_sha,
        "canonical_probe": witness["point"],
        "canonical_probe_rejected_by_release2_stage": rejected,
        "runtime_error": error_text,
        "canonical_eq45_velocity_export_ready": True,
        "kokuno_release2_stage_callable": True,
        "kokuno_velocity_export_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }


def main() -> None:
    print(json.dumps(audit_contract(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
