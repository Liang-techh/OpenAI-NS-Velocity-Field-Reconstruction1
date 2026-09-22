from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_main_pulse_float64_materialization_scope import (
    CANDIDATE_PATH,
    CONSTRAINTS_PATH,
    CONTRACT_PATH,
    PROJECT_STATUS_PATH,
    audit_contract,
    float64_coordinate_witness,
    load_and_audit,
)

ROOT = Path(__file__).resolve().parents[1]


def _contract() -> dict:
    return json.loads((ROOT / CONTRACT_PATH).read_text(encoding="utf-8"))


def _candidate_text() -> str:
    return "\n".join(
        (
            'SCHEMA = "kokuno-pa16-current-cartesian-main-pulse-v1"',
            "float_end = math.log(np.finfo(float).max / _FLOAT_X_SAFETY_DIVISOR)",
            "return min(source_end, float_end)",
            "return self.lambda_value * (self.log_X_materializable_end - self.log_X_p)",
            "self.xi_materializable_max < MAIN_XI_MAX",
            '"full_source_xi_11_current_cartesian_materialized": False',
            '"unified_global_cartesian_velocity_export_ready": False',
            'raise ValueError("Cartesian point lies beyond the finite-X current main-pulse domain")',
            "def velocity(self, x: Any, y: Any, z: Any, t: Any)",
        )
    )


def _constraints() -> dict:
    return {
        "nu": 0.01,
        "domain": {
            "physical": "R^3",
            "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
            "support": "r < 2 and abs(z) < 2",
            "time_interval": [0.25, 0.75],
        },
        "forcing": {
            "mode": "restricted_two_parameter_family",
            "restriction": "Only a,c may be fitted. No residual-dependent basis or pointwise free force. Freeze family before optimization.",
        },
        "nontriviality": {
            "reference_energy": 1.0,
            "reference_energy_abs_tolerance": 0.001,
        },
        "validation": {
            "seed": 914027,
            "held_out_points": 4096,
            "derivative_steps": [0.02, 0.01, 0.005],
            "quadrature_orders_per_axis": [24, 48, 96],
            "thresholds": {
                "pde_residual_max": 0.001,
                "pde_residual_L2": 0.001,
                "divergence_max": 1e-5,
                "divergence_L2": 1e-5,
            },
        },
    }


def _project_status() -> dict:
    return {
        "candidate_family": "eq45_supported_velocity_candidate_v1",
        "velocity_api": "openai_ns_reconstruction.eq45_supported_delivery:velocity",
        "states": {
            "velocity_export_ready": True,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def _audit(cfg: dict) -> list[str]:
    return audit_contract(cfg, _candidate_text(), _constraints(), _project_status())


def test_baseline_contract_is_fail_closed_and_consistent() -> None:
    assert _audit(_contract()) == []


def test_float64_coordinate_witness_is_nonvacuous_and_scoped() -> None:
    witness = float64_coordinate_witness()
    assert witness["classification"] == "autonomous_mechanics_only"
    assert witness["source_coordinate_endpoint_is_finite"] is True
    assert witness["source_endpoint_requires_X_beyond_float64"] is True
    assert witness["float64_prefix_ends_before_source_xi_end"] is True
    assert witness["not_source_numeric_evidence"] is True
    assert witness["not_openai_numeric_evidence"] is True
    assert witness["not_pde_evidence"] is True


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("scope", "base_pr"), 1094),
        (("scope", "base_exact_head"), "0" * 40),
        (("scope", "audited_candidate_blob_sha"), "f" * 40),
        (("scope", "candidate_bytes_changed"), True),
        (("scope", "source_formula_changed"), True),
        (("scope", "pressure_changed"), True),
        (("scope", "forcing_changed"), True),
        (("scope", "scientific_threshold_changed"), True),
        (("current_truth", "source_coordinate_main_pulse_kernel_defined_through_xi_11"), False),
        (("current_truth", "current_cartesian_main_pulse_finite_X_prefix_materialized"), False),
        (("current_truth", "float64_materialization_ceiling_is_public_source_support_endpoint"), True),
        (("current_truth", "float64_materialization_ceiling_is_physical_domain_boundary"), True),
        (("current_truth", "full_source_xi_11_current_cartesian_materialized"), True),
        (("current_truth", "terminal_global_velocity_materialized"), True),
        (("current_truth", "velocity_export_ready_for_kokuno_route"), True),
        (("current_truth", "visual_correspondence_verified"), True),
        (("current_truth", "pde_validated"), True),
        (("current_truth", "paper_exact"), True),
        (("current_truth", "openai_field_identified"), True),
        (("representation_contract", "ceiling_classification"), "public_source_support"),
        (("representation_contract", "source_endpoint_classification"), "float64_runtime_limit"),
        (("cr001_lock", "nu"), 0.02),
        (("cr001_lock", "forcing_mode"), "free_residual_force"),
        (("cr001_lock", "momentum_max"), 0.01),
        (("cr001_lock", "divergence_L2"), 1e-4),
        (("cr001_lock", "residual_defined_free_forcing_forbidden"), False),
        (("cr001_lock", "candidate_collapse_forbidden"), False),
        (("cr001_lock", "post_hoc_threshold_relaxation_forbidden"), False),
        (("canonical_eq45_independent_state", "velocity_export_ready"), False),
        (("canonical_eq45_independent_state", "pde_validated"), True),
    ],
)
def test_forbidden_contract_mutations_fail_closed(path: tuple[str, ...], value: object) -> None:
    cfg = copy.deepcopy(_contract())
    target = cfg
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    assert _audit(cfg), f"mutation unexpectedly admitted: {path}={value!r}"


@pytest.mark.parametrize(
    "token",
    [
        "float_end = math.log(np.finfo(float).max / _FLOAT_X_SAFETY_DIVISOR)",
        "return min(source_end, float_end)",
        "self.xi_materializable_max < MAIN_XI_MAX",
        '"full_source_xi_11_current_cartesian_materialized": False',
        '"unified_global_cartesian_velocity_export_ready": False',
        'raise ValueError("Cartesian point lies beyond the finite-X current main-pulse domain")',
        "def velocity(self, x: Any, y: Any, z: Any, t: Any)",
    ],
)
def test_candidate_representation_guards_are_required(token: str) -> None:
    text = _candidate_text().replace(token, "")
    assert audit_contract(_contract(), text, _constraints(), _project_status())


def test_source_fact_cannot_reclassify_float64_ceiling_as_source_endpoint() -> None:
    cfg = copy.deepcopy(_contract())
    cfg["provenance_classes"]["public_source_fact"] = [
        "pinned source proves the float64 Cartesian materialization ceiling is the source support endpoint",
        "xi=11",
    ]
    assert _audit(cfg)


def test_future_representation_must_get_new_identity_and_fresh_audits() -> None:
    cfg = copy.deepcopy(_contract())
    cfg["representation_contract"]["future_representation_identity_rule"] = "reuse all finite-prefix receipts unchanged"
    assert _audit(cfg)


def test_cr001_live_mutations_are_rejected() -> None:
    constraints = _constraints()
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.1
    assert audit_contract(_contract(), _candidate_text(), constraints, _project_status())


def test_eq45_live_state_is_not_downgraded_by_kokuno_representation_gap() -> None:
    status = _project_status()
    status["states"]["velocity_export_ready"] = False
    assert audit_contract(_contract(), _candidate_text(), _constraints(), status)


def test_repository_live_replay_on_exact_stacked_parent() -> None:
    assert (ROOT / CANDIDATE_PATH).exists(), "exact #1100 candidate must exist on stacked branch"
    assert (ROOT / CONSTRAINTS_PATH).exists()
    assert (ROOT / PROJECT_STATUS_PATH).exists()
    assert load_and_audit(ROOT) == []
