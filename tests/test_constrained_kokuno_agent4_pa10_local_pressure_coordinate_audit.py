from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_agent4_pa10_local_pressure_coordinate_audit import (
    AGENT1_CANDIDATE_SHA256,
    FINAL_PROJECT_GATES,
    _fd8_s,
    run_audit,
)
from openai_ns_reconstruction.kokuno_pa10_local_pressure_coordinate import (
    KokunoPA10LocalPressureCoordinate,
)


def test_fd8_public_value_derivative_matches_local_pressure_reference() -> None:
    interface = KokunoPA10LocalPressureCoordinate(radial_quadrature_points=64)
    y = np.asarray([0.61, 1.37, 2.44, 3.51], dtype=float)
    s = np.asarray([-0.49, -0.17, 0.21, 0.53], dtype=float)
    estimate = _fd8_s(interface, y, s, 0.02)
    reference = np.asarray(interface.evaluate_local(y, s)["p_s"], dtype=float)
    denominator = max(float(np.sqrt(np.mean(reference * reference))), 1.0e-300)
    relative_rms = float(np.sqrt(np.mean((estimate - reference) ** 2))) / denominator
    assert relative_rms < 5.0e-6


def test_independent_local_coordinate_receipt_passes_without_pde_promotion() -> None:
    receipt = run_audit()
    assert receipt["candidate_sha256"] == AGENT1_CANDIDATE_SHA256
    assert receipt["local_pressure_coordinate_independent_preflight_passed"] is True
    assert receipt["failed_guards"] == []
    assert receipt["p_s_sign_flip_relative_rms"] >= 0.5
    assert receipt["truth_boundary"]["prior_eta_only_preflight_passed"] is False
    assert receipt["truth_boundary"]["global_pressure_matched"] is False
    assert receipt["truth_boundary"]["global_leading_profile_reconstructed"] is False
    assert receipt["truth_boundary"]["heldout_ns_momentum_residual_assessed"] is False
    assert receipt["truth_boundary"]["formal_full_domain_pde_gate_assessed"] is False
    assert receipt["truth_boundary"]["pde_validated"] is False
    assert FINAL_PROJECT_GATES["normalized_momentum_max"] == 1.0e-3
    assert FINAL_PROJECT_GATES["normalized_momentum_L2"] == 1.0e-3
    assert FINAL_PROJECT_GATES["divergence_max"] == 1.0e-5
    assert FINAL_PROJECT_GATES["divergence_L2"] == 1.0e-5
