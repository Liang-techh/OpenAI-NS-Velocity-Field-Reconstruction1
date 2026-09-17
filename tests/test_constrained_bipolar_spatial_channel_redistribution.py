import numpy as np

from openai_ns_reconstruction.bipolar_spatial_channel_redistribution import REGION_NAMES, run


def test_bipolar_capped_spatial_channel_audit():
    report = run()
    assert report["resolutions"] == [24, 32, 40]
    assert report["times"] == [0.3125, 0.5, 0.75]
    assert report["odd3_sha256"] != report["capped_sha256"]

    # The capped search changed only the spatial swirl profile to numerical
    # precision; the poloidal Phi profile is frozen.
    assert report["phi_coefficient_max_abs_delta"] < 1e-12
    np.testing.assert_allclose(report["swirl_coefficient_max_abs_delta"], 2.609296672782034, rtol=2e-8)

    for row in report["rows"]:
        for candidate in ("odd3", "capped"):
            metrics = row[candidate]
            assert np.isfinite(metrics["energy"]["total"])
            assert metrics["energy"]["total"] > 0.0
            assert abs(sum(metrics["channel_fraction"].values()) - 1.0) < 1e-12
            assert abs(sum(metrics["region_fraction_of_total"].values()) - 1.0) < 1e-12
            assert set(metrics["regions"]) == set(REGION_NAMES)

    fine = report["finest_comparison"]
    np.testing.assert_allclose(
        [row["total_energy_relative_change"] for row in fine],
        [-4.5521620600467703e-05, -1.6814888876143792e-04, -2.0248241328281452e-04],
        rtol=2e-6,
        atol=2e-10,
    )
    np.testing.assert_allclose(
        [row["swirl_energy_relative_change"] for row in fine],
        [-8.559617264990799e-05, -3.1406857022264474e-04, -3.769233122414952e-04],
        rtol=2e-6,
        atol=2e-10,
    )
    # No compensating migration into the poloidal radial/axial energy channels.
    assert max(abs(row["radial_energy_relative_change"]) for row in fine) < 1e-12
    assert max(abs(row["axial_energy_relative_change"]) for row in fine) < 1e-12
    # All support-collar totals fall rather than grow; the late radial collar is
    # the clearest redistribution signal.
    for row in fine:
        assert row["radial_collar_total_relative_change"] <= 0.0
        assert row["axial_collar_total_relative_change"] <= 0.0
        assert row["corner_total_relative_change"] <= 0.0
    assert fine[-1]["radial_collar_total_relative_change"] < -0.10
    np.testing.assert_allclose(
        [row["angular_momentum_rms_relative_change"] for row in fine],
        [-4.078677757658826e-04, -9.517801100867158e-04, -1.1577894632286228e-03],
        rtol=2e-6,
        atol=2e-10,
    )

    # Preserve the nonzero central rotation and the signed radial/axial flow
    # direction used by the constrained search.
    for row in report["central_signed_components"]:
        for candidate in ("odd3", "capped"):
            assert row[candidate]["radial"] < 0.0
            assert row[candidate]["swirl"] > 0.0
            assert row[candidate]["axial"] > 0.0
        assert abs(row["capped"]["radial"] - row["odd3"]["radial"]) < 1e-12
        assert abs(row["capped"]["axial"] - row["odd3"]["axial"]) < 1e-12
        assert abs((row["capped"]["swirl"] - row["odd3"]["swirl"]) / row["odd3"]["swirl"]) < 5e-5

    # Smooth global measures are already stable on the medium-to-fine step.
    for row in report["medium_to_fine_stability"]:
        for candidate in ("odd3", "capped"):
            assert row[candidate]["total_energy_relative_change"] < 0.003
            assert row[candidate]["angular_momentum_rms_relative_change"] < 0.001

    fine_rows = [row for row in report["rows"] if row["resolution"] == 40]
    for row in fine_rows:
        # At this grid resolution the q50/q90 shell locations do not move under
        # the capped swirl redistribution, consistent with a small bulk change.
        assert row["odd3"]["swirl_concentration"] == row["capped"]["swirl_concentration"]

    assert report["truth_boundary"] == {
        "pde_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
    }
