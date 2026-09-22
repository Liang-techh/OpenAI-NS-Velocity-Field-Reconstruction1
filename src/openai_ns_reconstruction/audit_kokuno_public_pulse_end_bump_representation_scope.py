"""Fail-closed CR002 audit for the public pulse-end bump representation seam.

The public two-row end-compensation algebra and public bump centers/width are
source-form facts. The pointwise C-infinity bump used by exact A1 #1116 is a
repository-autonomous numerical realization. Algebraic closure under that
realization does not recover a source-exact bump or promote Cartesian/global,
visual, PDE, paper-exact, or OpenAI-field states.
"""
from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

CONTRACT_PATH = Path("configs/kokuno_public_pulse_end_bump_representation_scope.json")
SOURCE_PATH = Path("src/openai_ns_reconstruction/kokuno_public_pulse_end_compensator.py")
CONSTRAINTS_PATH = Path("configs/constraints.json")
STATUS_PATH = Path("project_status.json")

EXPECTED_BASE_HEAD = "c144d00fd81a938e497ca5d841fed6d1fe448a8b"
EXPECTED_PARENT_HEAD = "45da043dd2b4cd067f005a72c7e21fd0f2bcf309"
EXPECTED_SOURCE_BLOB = "3c4d7a2d51cbccad7322a9bf6825a69e2b45e486"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_SCHEMA = "kokuno-public-pulse-end-compensator-v1"

_REQUIRED_FALSE = {
    "source_exact_bump_shape_recovered",
    "source_exact_amplitude_root_materialized",
    "current_lineage_J_entry_materialized",
    "current_cartesian_end_compensation_composed",
    "source_hidden_parameters_recovered",
    "matched_global_pressure_materialized",
    "restricted_forcing_materialized",
    "heldout_ns_residual_assessed",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
}
_REQUIRED_TRUE = {
    "public_reconstruction_source",
    "public_end_compensation_linear_system_materialized",
    "public_bump_centers_and_width_materialized",
    "repository_autonomous_bump_shape_materialized",
}
_REQUIRED_SOLVE_ARGS = ["self", "amplitude", "m_entry_ratio", "j_entry_ratio"]


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def _git_blob_sha1(data: bytes) -> str:
    header = b"blob " + str(len(data)).encode("ascii") + b"\0"
    return hashlib.sha1(header + data).hexdigest()


def _module_literals_and_methods(path: Path) -> tuple[dict[str, Any], dict[str, list[str]]]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    literals: dict[str, Any] = {}
    methods: dict[str, list[str]] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                literals[node.targets[0].id] = ast.literal_eval(node.value)
            except Exception:
                pass
        if isinstance(node, ast.ClassDef) and node.name == "KokunoPublicPulseEndCompensator":
            for item in node.body:
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    methods[item.name] = [arg.arg for arg in item.args.args]
    return literals, methods


def _standard_bump(y: float, center: float, width: float, power: int) -> float:
    """Synthetic mechanics witness only, not a source reconstruction."""
    z = 2.0 * (y - center) / width
    if abs(z) >= 1.0:
        return 0.0
    base = math.exp(1.0 - 1.0 / (1.0 - z * z))
    return base**power


def _midpoint_integral(center: float, width: float, slope: float, y1: float, power: int) -> float:
    n = 2400
    lo = center - 0.5 * width
    step = width / n
    total = 0.0
    for k in range(n):
        y = lo + (k + 0.5) * step
        total += math.exp(slope * (y - y1)) * _standard_bump(y, center, width, power)
    return total * step


def _solve_2x2(
    matrix: tuple[tuple[float, float], tuple[float, float]], rhs: tuple[float, float]
) -> tuple[float, float]:
    (a, b), (c, d) = matrix
    det = a * d - b * c
    if abs(det) <= 1.0e-14:
        raise RuntimeError("synthetic witness matrix unexpectedly singular")
    x = (rhs[0] * d - b * rhs[1]) / det
    y = (a * rhs[1] - rhs[0] * c) / det
    return x, y


def _residual(
    matrix: tuple[tuple[float, float], tuple[float, float]],
    coeff: tuple[float, float],
    rhs: tuple[float, float],
) -> tuple[float, float]:
    return (
        matrix[0][0] * coeff[0] + matrix[0][1] * coeff[1] - rhs[0],
        matrix[1][0] * coeff[0] + matrix[1][1] * coeff[1] - rhs[1],
    )


def autonomous_bump_identity_witness() -> dict[str, float]:
    """Show that source-compatible pointwise bump choices change c1/c2 receipts.

    Both synthetic shapes are nonnegative C-infinity functions with the same
    centers, full support width, and peak one. They are mechanics examples only.
    """
    lam = 0.05
    width = 0.3
    y1 = 13.0 / lam - 3.0
    y2 = 13.0 / lam - 1.0
    slopes = (0.5 - lam, 0.5 - 2.0 * lam)
    centers = (y1, y2)

    matrices = []
    coeffs = []
    max_residuals = []
    rhs = (-0.1, -0.2)
    for power in (1, 2):
        matrix = tuple(
            tuple(_midpoint_integral(center, width, slope, y1, power) for center in centers)
            for slope in slopes
        )
        coeff = _solve_2x2(matrix, rhs)
        residual = _residual(matrix, coeff, rhs)
        matrices.append(matrix)
        coeffs.append(coeff)
        max_residuals.append(max(abs(residual[0]), abs(residual[1])))

    matrix_delta = max(
        abs(matrices[0][i][j] - matrices[1][i][j])
        for i in range(2)
        for j in range(2)
    )
    coefficient_delta = max(abs(coeffs[0][i] - coeffs[1][i]) for i in range(2))
    return {
        "matrix_linf_delta": matrix_delta,
        "coefficient_linf_delta": coefficient_delta,
        "shape1_closure_max": max_residuals[0],
        "shape2_closure_max": max_residuals[1],
    }


def _expect_equal(violations: list[str], label: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        violations.append(f"{label}: expected {expected!r}, got {actual!r}")


def audit_scope(root: Path = Path("."), contract: Mapping[str, Any] | None = None) -> list[str]:
    violations: list[str] = []
    payload = dict(contract) if contract is not None else _read_json(root / CONTRACT_PATH)

    _expect_equal(
        violations,
        "schema",
        payload.get("schema"),
        "cr002-kokuno-pulse-end-bump-representation-v1",
    )
    _expect_equal(
        violations,
        "scope",
        payload.get("scope"),
        "public_end_compensator_pointwise_bump_and_evidence_transfer_only",
    )

    stack = payload.get("exact_stack", {})
    _expect_equal(violations, "base head", stack.get("base_head"), EXPECTED_BASE_HEAD)
    _expect_equal(violations, "parent head", stack.get("parent_a1_head"), EXPECTED_PARENT_HEAD)
    _expect_equal(violations, "source blob", stack.get("source_module_git_blob_sha1"), EXPECTED_SOURCE_BLOB)
    _expect_equal(violations, "A2 consumes #1116", stack.get("current_a2_consumes_base_pr1116"), False)

    source_path = root / SOURCE_PATH
    source_bytes = source_path.read_bytes()
    _expect_equal(violations, "actual source blob", _git_blob_sha1(source_bytes), EXPECTED_SOURCE_BLOB)
    literals, methods = _module_literals_and_methods(source_path)
    _expect_equal(violations, "source schema", literals.get("SCHEMA"), EXPECTED_SCHEMA)
    _expect_equal(violations, "source parent", literals.get("PARENT_EXACT_HEAD"), EXPECTED_PARENT_HEAD)
    _expect_equal(
        violations,
        "source repository",
        literals.get("SOURCE_REPOSITORY"),
        "KokunoYumeto/yang-mills-interacting-workbench",
    )
    _expect_equal(
        violations,
        "source commit",
        literals.get("SOURCE_COMMIT"),
        "143f6773feb424ad9ed3a8d116653200f20346b7",
    )
    _expect_equal(violations, "lambda", literals.get("CURRENT_LAMBDA"), 0.05)
    _expect_equal(violations, "public bump width", literals.get("PUBLIC_BUMP_WIDTH"), 0.3)
    _expect_equal(violations, "quadrature order", literals.get("QUADRATURE_ORDER"), 96)
    _expect_equal(violations, "solve signature", methods.get("solve"), _REQUIRED_SOLVE_ARGS)

    source_truth = literals.get("_TRUTH_BOUNDARY", {})
    for key in _REQUIRED_TRUE:
        _expect_equal(violations, f"source truth {key}", source_truth.get(key), True)
    for key in _REQUIRED_FALSE:
        _expect_equal(violations, f"source truth {key}", source_truth.get(key), False)

    numerical = literals.get("_NUMERICAL_REALIZATION", {})
    bump_text = str(numerical.get("bump_shape", ""))
    if "repository-autonomous" not in bump_text or "peak=1" not in bump_text:
        violations.append(
            "source numerical realization must explicitly identify autonomous peak-one pointwise bump"
        )
    if "96-point Gauss-Legendre" not in str(numerical.get("quadrature", "")):
        violations.append("source numerical realization must retain the 96-point quadrature label")

    source_formulas = literals.get("_SOURCE_FORMULAS", {})
    _expect_equal(
        violations,
        "row slopes formula",
        source_formulas.get("row_slopes"),
        "s1=1/2-lambda, s2=1/2-2lambda",
    )
    if "width 0.3" not in str(source_formulas.get("bump_width", "")):
        violations.append("public bump-width source formula missing")

    sr = payload.get("source_vs_realization", {})
    expected_sr = {
        "public_bump_centers_and_width_materialized": True,
        "public_two_row_end_compensation_algebra_materialized": True,
        "repository_autonomous_pointwise_bump_materialized": True,
        "source_exact_pointwise_bump_recovered": False,
        "source_exact_amplitude_root_materialized": False,
        "current_lineage_J_entry_materialized": False,
        "current_cartesian_end_compensation_composed": False,
        "caller_supplies_amplitude": True,
        "caller_supplies_normalized_M_entry": True,
        "caller_supplies_normalized_J_entry": True,
    }
    for key, value in expected_sr.items():
        _expect_equal(violations, f"source_vs_realization {key}", sr.get(key), value)

    transfer = payload.get("evidence_transfer", {})
    required_transfer = {
        "same_public_centers_width_imply_same_pointwise_bump_identity": False,
        "same_public_algebra_implies_same_numerical_c1_c2_across_bump_realizations": False,
        "row_closure_under_autonomous_bump_proves_source_exact_bump": False,
        "row_closure_under_autonomous_bump_proves_source_exact_cartesian_end_compensation": False,
        "coefficient_receipts_transfer_across_pointwise_bump_replacement": False,
        "fresh_semantic_identity_required_if_pointwise_bump_changes": True,
        "fresh_save_load_and_numerical_replay_required_if_pointwise_bump_changes": True,
        "a2_pr1117_receipts_cover_pr1116_end_compensation": False,
    }
    for key, value in required_transfer.items():
        _expect_equal(violations, f"evidence_transfer {key}", transfer.get(key), value)

    truth = payload.get("truth_states", {})
    _expect_equal(
        violations, "truth public algebra", truth.get("public_end_compensator_algebra_ready"), True
    )
    for key in (
        "source_exact_bump_shape_recovered",
        "source_exact_amplitude_root_materialized",
        "current_lineage_J_entry_materialized",
        "current_cartesian_end_compensation_composed",
        "kokuno_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _expect_equal(violations, f"contract truth {key}", truth.get(key), False)

    cr001 = payload.get("canonical_cr001", {})
    canonical_expected = {
        "constraints_git_blob_sha1": EXPECTED_CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_max": 0.001,
        "momentum_L2": 0.001,
        "divergence_max": 1.0e-5,
        "divergence_L2": 1.0e-5,
        "residual_defined_free_forcing_forbidden": True,
        "candidate_collapse_forbidden": True,
        "post_hoc_threshold_relaxation_forbidden": True,
    }
    for key, value in canonical_expected.items():
        _expect_equal(violations, f"canonical_cr001 {key}", cr001.get(key), value)

    constraints_path = root / CONSTRAINTS_PATH
    _expect_equal(
        violations,
        "actual constraints blob",
        _git_blob_sha1(constraints_path.read_bytes()),
        EXPECTED_CONSTRAINTS_BLOB,
    )
    constraints = _read_json(constraints_path)
    _expect_equal(violations, "constraints nu", constraints.get("nu"), 0.01)
    domain = constraints.get("domain", {})
    _expect_equal(violations, "constraints physical domain", domain.get("physical"), "R^3")
    _expect_equal(
        violations,
        "constraints evaluation box",
        domain.get("evaluation_box"),
        [[-2, 2], [-2, 2], [-2, 2]],
    )
    _expect_equal(violations, "constraints support", domain.get("support"), "r < 2 and abs(z) < 2")
    _expect_equal(violations, "constraints time", domain.get("time_interval"), [0.25, 0.75])
    _expect_equal(
        violations,
        "constraints forcing",
        constraints.get("forcing", {}).get("mode"),
        "restricted_two_parameter_family",
    )
    _expect_equal(
        violations,
        "constraints energy",
        constraints.get("nontriviality", {}).get("reference_energy"),
        1.0,
    )
    validation = constraints.get("validation", {})
    _expect_equal(violations, "constraints validation seed", validation.get("seed"), 914027)
    _expect_equal(violations, "constraints held out", validation.get("held_out_points"), 4096)
    _expect_equal(
        violations,
        "constraints derivative steps",
        validation.get("derivative_steps"),
        [0.02, 0.01, 0.005],
    )
    _expect_equal(
        violations,
        "constraints quadrature",
        validation.get("quadrature_orders_per_axis"),
        [24, 48, 96],
    )
    thresholds = validation.get("thresholds", {})
    _expect_equal(violations, "constraints pde max", thresholds.get("pde_residual_max"), 0.001)
    _expect_equal(violations, "constraints pde L2", thresholds.get("pde_residual_L2"), 0.001)
    _expect_equal(violations, "constraints div max", thresholds.get("divergence_max"), 1.0e-5)
    _expect_equal(violations, "constraints div L2", thresholds.get("divergence_L2"), 1.0e-5)

    delivery = payload.get("canonical_delivery", {})
    delivery_expected = {
        "candidate_family": "eq45_supported_velocity_candidate_v1",
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, value in delivery_expected.items():
        _expect_equal(violations, f"canonical_delivery {key}", delivery.get(key), value)

    status = _read_json(root / STATUS_PATH)
    _expect_equal(
        violations,
        "status candidate family",
        status.get("candidate_family"),
        "eq45_supported_velocity_candidate_v1",
    )
    status_states = status.get("states", {})
    for key, value in delivery_expected.items():
        if key == "candidate_family":
            continue
        _expect_equal(violations, f"status states {key}", status_states.get(key), value)

    provenance = payload.get("four_way_provenance", {})
    for key in ("user_requirements", "public_source_facts", "autonomous_design", "pending_unknown"):
        values = provenance.get(key)
        if not isinstance(values, list) or not values or not all(
            isinstance(v, str) and v.strip() for v in values
        ):
            violations.append(f"four_way_provenance {key} must be a nonempty string list")

    forbidden = payload.get("promotion_rule", {}).get("forbidden_implications", [])
    required_forbidden = {
        "two_row_algebra_closes -> source_exact_pointwise_bump",
        "two_row_algebra_closes -> current_cartesian_end_compensation_composed",
        "two_row_algebra_closes -> pde_validated",
        "same_centers_width -> transferable_c1_c2_receipt",
    }
    if not required_forbidden.issubset(set(forbidden) if isinstance(forbidden, list) else set()):
        violations.append("promotion_rule is missing required forbidden implications")

    witness = autonomous_bump_identity_witness()
    if not witness["matrix_linf_delta"] > 1.0e-4:
        violations.append("autonomous bump witness matrix delta is vacuous")
    if not witness["coefficient_linf_delta"] > 1.0e-3:
        violations.append("autonomous bump witness coefficient delta is vacuous")
    if witness["shape1_closure_max"] > 1.0e-12 or witness["shape2_closure_max"] > 1.0e-12:
        violations.append("autonomous bump witness systems do not each close numerically")

    return violations


def main() -> int:
    violations = audit_scope()
    if violations:
        for item in violations:
            print(f"CR002 VIOLATION: {item}")
        return 1
    witness = autonomous_bump_identity_witness()
    print(json.dumps({"status": "PASS", "witness": witness}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
