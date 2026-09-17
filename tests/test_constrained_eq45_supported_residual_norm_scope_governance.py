import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_supported_residual_norm_scope_governance import (
    audit_supported_residual_norm_scope,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = json.loads((ROOT / "configs/eq45_supported_residual_norm_scope.json").read_text())


def _constraints():
    return {
        "experiment_id": "compact_axisymmetric_window_v1",
        "forcing": {"mode": "restricted_two_parameter_family"},
        "validation": {
            "held_out_points": 4096,
            "derivative_steps": [0.02, 0.01, 0.005],
            "norms": [
                "max Euclidean vector norm",
                "volume-weighted L2 spatial norm at each time",
            ],
            "thresholds": {
                "pde_residual_max": 0.001,
                "pde_residual_L2": 0.001,
            },
        },
    }


def _diagnostic():
    return {
        "schema": "eq45_supported_vorticity_revalidation_v1",
        "claim_scope": "independent pressure_free_vorticity_equation_on_supported_child",
        "parent_sha256": "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7",
        "supported_child_sha256": "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d",
        "spatial_steps": [0.02, 0.01, 0.005],
        "probe_count": 16,
        "finest_supported_zero_force_rms": 3.58728819303,
        "finest_supported_zero_force_max": 13.0812827018,
        "finest_supported_term_normalized_rms": 0.9468023011,
        "forcing_family": "restricted_two_parameter_family",
        "supported_force_coefficients_status": "pending_refit_and_independent_revalidation",
        "untapered_fitted_force_transferred": False,
        "pde_thresholds": {"max": 0.001, "L2": 0.001},
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "pressure and supported-child restricted-force coefficients are not bound here; "
            "probe RMS is not the preregistered volume-weighted spatial L2"
        ),
    }


def _truth():
    return {
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def test_current_scope_keeps_probe_metrics_out_of_formal_pde_gate():
    result = audit_supported_residual_norm_scope(CONTRACT, _constraints(), _diagnostic(), _truth())
    assert result["probe_count"] == 16
    assert result["registered_held_out_points"] == 4096
    assert result["fixed_probe_metrics_are_formal_pde_acceptance"] is False
    assert result["formal_pde_gate_assessed"] is False
    assert result["velocity_export_ready"] is True
    assert result["pde_validated"] is False


@pytest.mark.parametrize(
    "mutation, match",
    [
        (lambda c, k, d, t: d.__setitem__("formal_pde_gate_assessed", True), "cannot assess"),
        (lambda c, k, d, t: c["policy"].__setitem__("fixed_probe_rms_may_satisfy_registered_L2", True), "illegal policy promotion"),
        (lambda c, k, d, t: c["policy"].__setitem__("fixed_probe_max_may_satisfy_registered_max", True), "illegal policy promotion"),
        (lambda c, k, d, t: c["policy"].__setitem__("green_ci_may_set_pde_validated", True), "illegal policy promotion"),
        (lambda c, k, d, t: k["validation"]["thresholds"].__setitem__("pde_residual_L2", 10.0), "threshold drifted"),
        (lambda c, k, d, t: d.__setitem__("forcing_family", "residual_defined"), "forcing family drifted"),
        (lambda c, k, d, t: d.__setitem__("untapered_fitted_force_transferred", True), "cannot be inherited"),
        (lambda c, k, d, t: t.__setitem__("pde_validated", True), "promoted without formal validation"),
        (lambda c, k, d, t: t.__setitem__("velocity_export_ready", False), "lost export readiness"),
        (lambda c, k, d, t: c["classification"].__setitem__("openai_visual_correspondence", "public_source_fact"), "cannot be promoted"),
    ],
)
def test_norm_scope_mutations_fail_closed(mutation, match):
    contract = copy.deepcopy(CONTRACT)
    constraints = _constraints()
    diagnostic = _diagnostic()
    truth = _truth()
    mutation(contract, constraints, diagnostic, truth)
    with pytest.raises(ValueError, match=match):
        audit_supported_residual_norm_scope(contract, constraints, diagnostic, truth)
