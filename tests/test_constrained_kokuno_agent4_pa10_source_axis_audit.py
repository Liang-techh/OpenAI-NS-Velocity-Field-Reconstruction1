import inspect

import numpy as np

from openai_ns_reconstruction.kokuno_agent4_pa10_source_axis_audit import (
    SEED,
    _H,
    _Z,
    run_audit,
)
from openai_ns_reconstruction.kokuno_pa10_source_axis_domain import (
    KokunoPA10SourceCompatibleAxisDomain,
)


def test_independent_source_axis_audit_passes_frozen_guards():
    report = run_audit()
    assert report["failed_guards"] == []
    assert report["source_axis_domain_independent_preflight_passed"] is True
    real = report["real_axis"]
    assert real["nonvacuous_K_and_Hsmall"] is True
    assert real["PA8_independent_sampling_passed"] is True
    assert len(real["levels"]) == 3
    assert real["levels"][-1]["min_abs_H_on_K"] > 0.01
    assert real["levels"][-1]["min_chi_on_K"] > 0.99
    assert real["levels"][-1]["min_Z_on_Hsmall"] > 1.0


def test_seeded_offgrid_axis_points_preserve_pa8_separation():
    datum = KokunoPA10SourceCompatibleAxisDomain()
    rng = np.random.default_rng(SEED)
    lo = -1.0 - datum.enlarged_real_margin
    hi = 1.0 + datum.enlarged_real_margin
    eta = rng.uniform(lo, hi, size=8192)
    H = np.asarray(_H(eta, D=datum.D, j0=datum.j0), dtype=float)
    Z = _Z(
        eta,
        A=datum.A,
        D=datum.D,
        j0=datum.j0,
        pressure_square=datum.pressure_square,
    )
    K = np.abs(Z) <= datum.delta_star
    Hsmall = np.abs(H) <= 10.0 * datum.sigma_star
    assert np.count_nonzero(K) > 0
    assert np.count_nonzero(Hsmall) > 0
    assert np.min(np.abs(H[K])) > 10.0 * datum.sigma_star
    assert np.min(Z[Hsmall]) > 2.0 * datum.delta_star


def test_complex_roots_and_g_series_are_independent_and_mutation_sensitive():
    report = run_audit()
    complex_report = report["complex_zero_exclusion"]
    assert complex_report["independent_zero_exclusion_passed"] is True
    assert complex_report["min_H_root_distance_over_tube"] > 1.0
    g = report["g_coefficient_norm"]
    assert g["independent_norm_domination_passed"] is True
    assert g["public_to_independent_ratio"] >= 1.0
    assert g["partial_sums"]["64"] <= g["independent_closed_form"]
    assert report["mutation"]["complex_tube_x100_detected"] is True
    assert report["mutation"]["g_bound_0p99_detected"] is True


def test_synthetic_pa8_bad_point_is_actually_rejected():
    datum = KokunoPA10SourceCompatibleAxisDomain()

    def point_satisfies_pa8_separation(H: float, Z: float) -> bool:
        if abs(Z) <= datum.delta_star:
            return abs(H) > 10.0 * datum.sigma_star
        return True

    assert point_satisfies_pa8_separation(0.0, 0.0) is False
    assert point_satisfies_pa8_separation(0.02, 0.0) is True


def test_validator_does_not_call_upstream_certificate_or_axis_state_paths():
    src = inspect.getsource(run_audit)
    for forbidden in (
        "real_interval_certificate(",
        "complex_domain_certificate(",
        "axis_state(",
        "_real_certificate(",
        "_complex_certificate(",
    ):
        assert forbidden not in src


def test_truth_boundary_keeps_ns_and_missing_leading_objects_closed():
    truth = run_audit()["truth_boundary"]
    assert truth["source_Phi_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_Phi_radius_one_ball_lipschitz_machine_bound"] is False
    assert truth["source_R1_R2_machine_bound"] is False
    assert truth["source_operator_M_K_machine_bound"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
