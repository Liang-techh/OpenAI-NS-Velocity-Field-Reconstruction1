"""Fail-closed CR002 audit for the Kokuno strict-inner radial-stress lower-limit scope."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

CONTRACT_REL = Path("configs/kokuno_strict_inner_radial_stress_lower_limit_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
EXPECTED_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class AuditError(RuntimeError):
    """Raised when the governed lower-limit contract drifts."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot load JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise AuditError(f"{path} must contain a JSON object")
    return value


def _git_blob_sha(text: str) -> str:
    payload = text.encode("utf-8")
    return hashlib.sha1(
        b"blob " + str(len(payload)).encode("ascii") + b"\0" + payload
    ).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def _audit_classification(contract: dict[str, Any]) -> None:
    classification = contract.get("classification")
    _require(isinstance(classification, dict), "classification must be an object")
    _require(
        set(classification) == EXPECTED_CLASSES,
        "classification must use exactly four provenance classes",
    )
    seen: set[str] = set()
    for name in sorted(EXPECTED_CLASSES):
        values = classification[name]
        _require(
            isinstance(values, list) and values,
            f"classification.{name} must be a nonempty list",
        )
        for value in values:
            _require(
                isinstance(value, str) and value,
                f"classification.{name} contains invalid item",
            )
            _require(
                value not in seen,
                f"provenance item appears in more than one class: {value}",
            )
            seen.add(value)

    public = set(classification["public_source_fact"])
    autonomous = set(classification["autonomous_design"])
    _require(
        "registered_compact_inverse_structure_Me_Pe_sigmae_with_formal_axis_lower_limit_zero"
        in public,
        "registered formal inverse structure must remain a public-source fact",
    )
    for item in (
        "strict_inner_positive_radius_sampling",
        "finite_sample_trapezoid_moment",
        "finite_sample_moment_complement",
        "cumulative_trapezoid_initialized_at_first_supplied_radius",
        "manufactured_lower_limit_witness",
        "cr001_validation_protocol",
    ):
        _require(
            item in autonomous,
            f"autonomous discretization/protocol item drifted: {item}",
        )


def _audit_target_sources(repo_root: Path, contract: dict[str, Any]) -> None:
    target = contract["target"]
    _require(target["pr"] == 858, "target PR must remain #858")
    _require(
        target["head"] == "9140b0bab1bfc34d496460675e92e63fc6dd721b",
        "target exact head drifted",
    )

    target_path = repo_root / target["path"]
    helper_path = repo_root / target["helper_path"]
    try:
        target_text = target_path.read_text(encoding="utf-8")
        helper_text = helper_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise AuditError(f"cannot read governed source: {exc}") from exc

    _require(
        _git_blob_sha(target_text) == target["blob_sha"],
        "target source blob drifted",
    )
    _require(
        _git_blob_sha(helper_text) == target["helper_blob_sha"],
        "radial-stress helper blob drifted",
    )

    for needle in (
        "sigma_e(r) = -r^(-e) integral_0^r s^e P_e(F)(s) ds",
        "radii must be finite and strictly positive",
        "_compact_radial_stress(",
        '"complete_ns_defect": False',
        '"pde_validated": False',
        '"paper_exact": False',
        '"openai_field_identified": False',
    ):
        _require(needle in target_text, f"target source invariant missing: {needle}")

    for needle in (
        "out = np.zeros_like(values)",
        "out[1:] = np.cumsum(increments)",
        "moment = float(_cumulative_trapezoid(weight * f, r)[-1])",
        "primitive = _cumulative_trapezoid(weight * complement, r)",
        "stress = -primitive / weight",
        '"stress_inner_edge": float(stress[0])',
    ):
        _require(needle in helper_text, f"helper lower-limit invariant missing: {needle}")


def _audit_lower_limit_semantics(contract: dict[str, Any]) -> None:
    scope = contract["radial_inverse_scope"]
    _require(
        scope["formal_lower_limit"] == 0.0,
        "formal radial inverse lower limit must remain zero",
    )
    _require(
        scope["realized_discrete_lower_limit"]
        == "r_min = first supplied strictly-positive radial sample",
        "realized discrete lower limit must remain the first positive sample",
    )
    _require(
        scope["realized_primitive_initialization"]
        == "primitive[0]=0 by cumulative trapezoid construction",
        "inner-edge zero must remain attributed to cumulative integration initialization",
    )
    _require(
        scope["current_equivalence_to_formal_axis_inverse_proven"] is False,
        "formal/truncated equivalence is not proven",
    )
    _require(
        scope["stress_inner_edge_zero_is_axis_regularity_evidence"] is False,
        "inner-edge zero is not axis-regularity evidence",
    )
    _require(
        scope["sampled_weighted_moment_is_formal_full_domain_moment"] is False,
        "sampled moment is not a formal full-domain moment",
    )
    _require(
        scope["strict_inner_radii_include_axis"] is False,
        "strict-inner positive-radius grid must not be relabeled as including the axis",
    )
    required = set(scope["equivalence_requires"])
    _require(
        required
        == {
            "missing integral_0^rmin s^e P_e(F)(s) ds equals zero",
            "finite sampled weighted moment/complement matches the intended formal radial domain semantics",
        },
        "formal/truncated equivalence requirements drifted",
    )


def _audit_mechanics_witness(contract: dict[str, Any]) -> dict[str, float]:
    witness = contract["manufactured_mechanics_witness"]
    _require(
        witness["status"] == "autonomous_mechanics_only_not_candidate_evidence",
        "manufactured witness must remain autonomous mechanics only",
    )
    _require(witness["profile"] == "P_e(s)=1", "manufactured profile drifted")
    r_min = float(witness["r_min"])
    _require(
        math.isclose(r_min, 0.2, rel_tol=0.0, abs_tol=1e-15),
        "manufactured r_min drifted",
    )
    _require(witness["exponents"] == [1, 2], "manufactured exponent set drifted")
    computed: dict[str, float] = {}
    for exponent in (1, 2):
        formal = -r_min / (exponent + 1)
        key = f"e{exponent}"
        expected = float(witness["formal_sigma_at_r_min"][key])
        _require(
            math.isclose(formal, expected, rel_tol=0.0, abs_tol=1e-15),
            f"formal witness drifted for e={exponent}",
        )
        computed[key] = formal
    _require(
        float(witness["truncated_discrete_sigma_at_r_min"]) == 0.0,
        "truncated primitive must start at zero",
    )
    _require(
        witness["missing_term_formula_at_r"] == "-r_min^(e+1)/((e+1)*r^e)",
        "missing lower-limit term formula drifted",
    )
    return computed


def _audit_cr001(repo_root: Path, contract: dict[str, Any]) -> None:
    live = _load_json(repo_root / CONSTRAINTS_REL)
    snap = contract["cr001_snapshot"]

    _require(live["nu"] == snap["nu"] == 0.01, "CR001 viscosity drifted")
    domain = live["domain"]
    _require(
        domain["physical"] == snap["physical_domain"] == "R^3",
        "CR001 physical domain drifted",
    )
    _require(
        domain["evaluation_box"] == snap["evaluation_box"],
        "CR001 evaluation box drifted",
    )
    _require(domain["support"] == snap["support"], "CR001 support drifted")
    _require(
        domain["time_interval"] == snap["time_interval"],
        "CR001 time interval drifted",
    )

    forcing = live["forcing"]
    _require(
        forcing["mode"] == snap["forcing_mode"],
        "CR001 forcing mode drifted",
    )
    _require(
        forcing["parameters"] == snap["forcing_parameters"],
        "CR001 forcing bounds drifted",
    )
    _require(
        snap["free_residual_force_allowed"] is False,
        "free residual force must remain forbidden",
    )
    restriction = forcing["restriction"].lower()
    _require(
        "no residual-dependent basis" in restriction
        and "pointwise free force" in restriction,
        "CR001 free-force prohibition drifted",
    )

    nontrivial = live["nontriviality"]
    _require(
        nontrivial["reference_energy"] == snap["reference_energy"],
        "CR001 reference energy drifted",
    )
    _require(
        nontrivial["reference_energy_abs_tolerance"]
        == snap["reference_energy_abs_tolerance"],
        "CR001 energy tolerance drifted",
    )
    _require(
        snap["amplitude_collapse_allowed"] is False,
        "amplitude collapse must remain forbidden",
    )
    _require(
        "reject collapsed candidates" in nontrivial["enforcement"],
        "CR001 collapsed-candidate rejection drifted",
    )

    validation = live["validation"]
    _require(
        validation["seed"] == snap["validation_seed"],
        "CR001 validation seed drifted",
    )
    _require(
        validation["held_out_points"] == snap["held_out_points"],
        "CR001 held-out point count drifted",
    )
    _require(
        validation["times"] == snap["validation_times"],
        "CR001 validation times drifted",
    )
    _require(
        validation["derivative_steps"] == snap["derivative_steps"],
        "CR001 derivative ladder drifted",
    )
    _require(
        validation["quadrature_orders_per_axis"]
        == snap["quadrature_orders_per_axis"],
        "CR001 quadrature ladder drifted",
    )
    thresholds = validation["thresholds"]
    for key in (
        "divergence_max",
        "divergence_L2",
        "pde_residual_max",
        "pde_residual_L2",
    ):
        _require(thresholds[key] == snap[key], f"CR001 threshold drifted: {key}")
    _require(
        "volume-weighted L2 spatial norm at each time" in validation["norms"],
        "CR001 volume-weighted L2 norm semantics drifted",
    )
    _require(
        "changing thresholds requires a new experiment version"
        in validation["failure_policy"],
        "CR001 no-posthoc-threshold-relaxation policy drifted",
    )


def _audit_truth_boundary(contract: dict[str, Any]) -> None:
    truth = contract["truth_boundary"]
    _require(
        truth["scoped_truncated_radial_stress_diagnostic_exists"] is True,
        "scoped diagnostic truth drifted",
    )
    for key in (
        "formal_axis_based_inverse_verified",
        "formal_full_domain_moment_verified",
        "axis_regularity_verified",
        "complete_ns_defect",
        "correction_target_admitted",
        "velocity_export_ready_promoted_by_this_audit",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(truth[key] is False, f"unauthorized truth promotion: {key}")
    _require(
        truth["pde_pending_blocks_callable_velocity_delivery"] is False,
        "PDE pending must not block an independently callable velocity delivery",
    )

    forbidden = set(contract["forbidden_promotions"])
    expected = {
        "truncated_positive_radius_primitive_to_formal_axis_inverse",
        "stress_inner_edge_zero_to_axis_regularity",
        "sampled_weighted_moment_to_formal_full_domain_moment",
        "scoped_transport_stress_to_global_correction_target",
        "scoped_transport_stress_to_complete_ns_or_pde_validation",
        "autonomous_discretization_to_public_source_fact",
        "posthoc_threshold_relaxation",
        "residual_defined_or_pointwise_free_forcing",
        "amplitude_collapse",
        "pde_pending_to_callable_velocity_delivery_blocker",
    }
    _require(forbidden == expected, "forbidden-promotion set drifted")


def audit(repo_root: Path) -> dict[str, Any]:
    repo_root = Path(repo_root).resolve()
    contract = _load_json(repo_root / CONTRACT_REL)
    _require(contract.get("schema_version") == 1, "unsupported contract schema")
    _require(
        contract.get("task")
        == "CR002-KOKUNO-STRICT-INNER-RADIAL-STRESS-LOWER-LIMIT-SCOPE",
        "task identity drifted",
    )
    _require(
        contract.get("status") == "governance_only_not_scientific_admission",
        "governance status drifted",
    )
    _audit_classification(contract)
    _audit_target_sources(repo_root, contract)
    _audit_lower_limit_semantics(contract)
    witness = _audit_mechanics_witness(contract)
    _audit_cr001(repo_root, contract)
    _audit_truth_boundary(contract)
    siblings = contract["sibling_governance"]
    _require(
        siblings["duplicate_of_sibling"] is False,
        "this increment must remain non-duplicate",
    )
    return {
        "ok": True,
        "task": contract["task"],
        "target_pr": contract["target"]["pr"],
        "target_head": contract["target"]["head"],
        "formal_lower_limit": contract["radial_inverse_scope"]["formal_lower_limit"],
        "realized_lower_limit": contract["radial_inverse_scope"][
            "realized_discrete_lower_limit"
        ],
        "manufactured_formal_sigma_at_rmin": witness,
        "formal_axis_based_inverse_verified": False,
        "pde_validated": False,
        "pde_pending_blocks_callable_velocity_delivery": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args()
    try:
        result = audit(args.repo_root)
    except AuditError as exc:
        if args.as_json:
            print(json.dumps({"ok": False, "error": str(exc)}, sort_keys=True))
        else:
            print(f"FAIL: {exc}")
        raise SystemExit(1)
    if args.as_json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print("PASS: CR002 Kokuno strict-inner radial-stress lower-limit scope")


if __name__ == "__main__":
    main()
