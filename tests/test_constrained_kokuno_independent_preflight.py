import numpy as np

from openai_ns_reconstruction.kokuno_complete_curl import KokunoCompleteCurlCorrection
from openai_ns_reconstruction.kokuno_independent_preflight import (
    MUTATION_DIVERGENCE,
    STEPS,
    TIMES,
    run_independent_preflight,
)


def test_kokuno_oscillatory_preflight_is_fourth_order_and_mutation_sensitive():
    report = run_independent_preflight(KokunoCompleteCurlCorrection())

    assert report["structural_preflight_passed"] is True
    assert report["pde_validated"] is False
    assert report["registered_full_ns_gate"]["assessed"] is False
    assert report["registered_full_ns_gate"]["normalized_residual_threshold"] == 1.0e-3
    assert report["truth_boundary"]["training_loss_used"] is False
    assert report["truth_boundary"]["candidate_internal_derivative_helpers_used"] is False
    assert report["outside_support_max"] == 0.0

    rows = report["rows"]
    assert len(rows) == len(TIMES) * len(STEPS)
    fine = [row for row in rows if row["step"] == min(STEPS)]
    assert len(fine) == len(TIMES)
    assert max(row["curl_error_rms"] for row in fine) < 1.0e-7
    assert max(row["divergence_rms"] for row in fine) < 1.0e-7
    assert max(abs(row["mutation_divergence_rms"] - MUTATION_DIVERGENCE) for row in fine) < 5.0e-5

    orders = [
        value
        for row in rows
        for key, value in row.items()
        if key in {"curl_rms_order_to_next", "divergence_rms_order_to_next"}
    ]
    assert min(orders) > 3.5

    activity = report["velocity_activity"]
    assert min(value["rms"] for value in activity.values()) > 1.0e-3


def test_kokuno_oscillatory_preflight_covers_axis_and_support_strata():
    report = run_independent_preflight()
    regions = report["finest_region_metrics"]

    for time in map(str, TIMES):
        assert set(regions[time]) == {"interior", "collar", "axis_near", "corner"}
        for region in regions[time].values():
            values = np.array(
                [
                    region["curl_error"]["max"],
                    region["curl_error"]["rms"],
                    region["divergence"]["max"],
                    region["divergence"]["rms"],
                ]
            )
            assert np.all(np.isfinite(values))
            assert np.all(values >= 0.0)
        assert regions[time]["axis_near"]["curl_error"]["max"] < 1.0e-7
        assert regions[time]["axis_near"]["divergence"]["max"] < 1.0e-7
