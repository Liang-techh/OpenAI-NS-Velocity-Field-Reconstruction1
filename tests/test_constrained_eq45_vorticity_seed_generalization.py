import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_vorticity_seed_generalization import (
    DEFAULT_SEEDS,
    audit_eq45_vorticity_seed_generalization,
)


EXPECTED_ZERO_RMS = np.array(
    [
        26.93980218230254,
        17.688755891104446,
        21.8224785687113,
        21.025893787310384,
        23.67579761046505,
    ]
)
EXPECTED_RESTRICTED_RMS = np.array(
    [
        26.44838938851789,
        17.14741307923976,
        21.278429756208784,
        20.57403239434794,
        23.363677692690082,
    ]
)
EXPECTED_REDUCTION = np.array(
    [
        0.018241143363237938,
        0.030603781023227403,
        0.024930660868310573,
        0.021490710337134445,
        0.013183079316280625,
    ]
)


def test_checked_eq45_vorticity_obstruction_generalizes_across_probe_seeds():
    report = audit_eq45_vorticity_seed_generalization()

    assert report["schema"] == "eq45_vorticity_seed_generalization_v1"
    assert report["candidate_sha256"] == (
        "48f1845fcfd71ec95495c98bb7aac3fca4653e748a8856a5421234e9601525a7"
    )
    assert report["validation_seeds"] == list(DEFAULT_SEEDS)
    assert report["spatial_steps"] == [0.02, 0.01, 0.005]
    assert report["reported_spatial_step"] == 0.005
    assert report["sample_count_per_seed"] == 12

    rows = report["rows"]
    zero_rms = np.array([row["zero_force_rms"] for row in rows])
    restricted_rms = np.array([row["restricted_force_rms"] for row in rows])
    reduction = np.array([row["force_rms_reduction_fraction"] for row in rows])
    np.testing.assert_allclose(zero_rms, EXPECTED_ZERO_RMS, rtol=0.0, atol=5e-6)
    np.testing.assert_allclose(
        restricted_rms, EXPECTED_RESTRICTED_RMS, rtol=0.0, atol=5e-6
    )
    np.testing.assert_allclose(reduction, EXPECTED_REDUCTION, rtol=0.0, atol=5e-7)

    summary = report["cross_seed_summary"]
    assert summary["zero_force_rms"]["min"] == pytest.approx(17.688755891104446)
    assert summary["zero_force_rms"]["median"] == pytest.approx(21.8224785687113)
    assert summary["zero_force_rms"]["mean"] == pytest.approx(22.230545607978744)
    assert summary["zero_force_rms"]["max"] == pytest.approx(26.93980218230254)
    assert summary["zero_force_rms"]["coefficient_of_variation"] == pytest.approx(
        0.13719700342279872
    )
    assert summary["force_rms_reduction_fraction"]["min"] == pytest.approx(
        0.013183079316280625
    )
    assert summary["force_rms_reduction_fraction"]["max"] == pytest.approx(
        0.030603781023227403
    )

    assert report["velocity_changed"] is False
    assert report["pressure_fitted"] is False
    assert report["forcing_fitted_in_this_audit"] is False
    assert report["training_holdout_separate"] is True
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
    assert report["blowup_proved"] is False


@pytest.mark.parametrize(
    "seeds",
    [
        (1, 2),
        (1, 1, 2),
        (1, 2, -3),
        (1, 2, 3.5),
    ],
)
def test_seed_contract_fails_closed(seeds):
    with pytest.raises(ValueError):
        audit_eq45_vorticity_seed_generalization(seeds=seeds)
