import numpy as np

from openai_ns_reconstruction.constrained_eq45_bipolar_f20_theta_overlap import (
    audit_f20_theta_overlap,
)


def test_f20_theta_overlap_is_finite_pure_and_more_flank_selective():
    report = audit_f20_theta_overlap()

    assert report["sample_count"] == 216
    assert report["baseline_theta_error"]["rms_all"] > 0.0
    assert np.isfinite(report["baseline_theta_error"]["radial_flank_to_interior_rms"])

    for name in ("F10", "F20"):
        mode = report["modes"][name]
        assert mode["swirl_fraction"] > 1.0 - 1e-10
        assert mode["regions"]["interior"]["swirl_response_rms"] > 0.0
        assert mode["regions"]["radial_flank"]["swirl_response_rms"] > 0.0
        assert np.isfinite(mode["global_overlap"]["error_enrichment"])
        assert 0.0 <= mode["global_overlap"]["abs_response_abs_error_cosine"] <= 1.0
        assert 0.0 <= mode["global_overlap"]["top_quartile_error_response_energy_share"] <= 1.0

    assert report["F20_vs_F10_radial_flank_selectivity_ratio"] > 1.5
    assert report["truth_boundary"]["coefficient_fitted"] is False
    assert report["truth_boundary"]["perturbed_candidate_pde_residual_evaluated"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False
