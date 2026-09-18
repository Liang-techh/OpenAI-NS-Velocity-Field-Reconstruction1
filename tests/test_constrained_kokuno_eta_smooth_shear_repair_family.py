import json
import math

import numpy as np

from openai_ns_reconstruction.kokuno_eta_smooth_shear_repair_family import (
    KokunoEtaSmoothShearRepairFamily,
)


def test_eta_smooth_family_closes_selected_modulation_on_heldout_eta():
    family = KokunoEtaSmoothShearRepairFamily()
    report = family.family_report()

    assert report["eta_even_by_construction"] is True
    assert report["analytic_eta_derivative"] is True
    assert report["source_admissible_loop_reconstructed"] is False
    assert report["selected_autonomous_modulation_consumed"] is True
    assert report["heldout_midpoint_max_abs_moment_residual"] <= report["closure_tolerance"]
    assert report["max_abs_coefficient"] < report["coefficient_limit"]

    etas = np.asarray([-0.93, -0.67, -0.41, -0.13, 0.19, 0.47, 0.73, 0.91])
    closure = family.closure_report(etas)
    assert closure["passed"] is True
    assert closure["max_abs_residual"] <= family.closure_tolerance


def test_even_polynomial_family_and_analytic_eta_derivative():
    family = KokunoEtaSmoothShearRepairFamily()
    etas = np.asarray([0.0, 0.17, 0.41, 0.78, 0.97])
    positive = family.coefficients(etas)
    negative = family.coefficients(-etas)
    assert np.allclose(positive, negative, rtol=0.0, atol=2.0e-14)

    d_positive = family.coefficient_eta(etas)
    d_negative = family.coefficient_eta(-etas)
    assert np.allclose(d_positive, -d_negative, rtol=0.0, atol=2.0e-13)
    assert np.array_equal(family.coefficient_eta(0.0), np.zeros(5))

    eta = 0.37
    step = 2.0e-6
    fd = (family.coefficients(eta + step) - family.coefficients(eta - step)) / (2.0 * step)
    analytic = family.coefficient_eta(eta)
    assert np.allclose(analytic, fd, rtol=2.0e-6, atol=2.0e-10)


def test_smooth_family_drives_nonzero_i1_profile_and_velocity_correction():
    family = KokunoEtaSmoothShearRepairFamily()
    repair = family.repair
    eta = 0.2
    xi = 1.18
    log_X = repair.log_X_1 + math.log(xi)
    profile = family.profile_correction_logX(log_X, eta)

    for key in (
        "delta_E",
        "delta_E_X",
        "delta_E_eta",
        "delta_F",
        "delta_U",
        "delta_U_X",
        "delta_U_eta",
        "delta_M",
        "delta_M_eta",
        "delta_v0",
    ):
        assert np.all(np.isfinite(profile[key]))
    assert abs(float(profile["delta_E"])) > 0.0

    # The existing native-coordinate contract has q=1, eta=0 at t=z=0, so
    # choosing r=sqrt(2X) evaluates the same I1 scale through the public xyz,t path.
    X = repair.X_1 * xi
    radius = math.sqrt(2.0 * X)
    velocity = family.velocity_correction(radius, 0.0, 0.0, 0.0)
    assert velocity.shape == (3,)
    assert np.all(np.isfinite(velocity))
    assert np.linalg.norm(velocity) > 0.0


def test_eta_smooth_family_serialization_preserves_truth_boundary(tmp_path):
    family = KokunoEtaSmoothShearRepairFamily(interpolation_nodes=9)
    path = family.save_json(tmp_path / "eta_smooth_pa17.json")
    loaded = KokunoEtaSmoothShearRepairFamily.load_json(path)
    assert loaded.sha256 == family.sha256
    assert loaded.to_payload() == family.to_payload()

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["truth_boundary"]["eta_smooth_repair_coefficient_family_reconstructed"] is True
    assert payload["truth_boundary"]["I1_repair_applied_to_selected_autonomous_modulation"] is True
    assert payload["truth_boundary"]["source_admissible_loop_reconstructed"] is False
    assert payload["truth_boundary"]["actual_source_admissible_loop_discrepancy_supplied"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False

    payload["truth_boundary"]["source_admissible_loop_reconstructed"] = True
    path.write_text(json.dumps(payload), encoding="utf-8")
    try:
        KokunoEtaSmoothShearRepairFamily.load_json(path)
    except ValueError as exc:
        assert "hash or content mismatch" in str(exc)
    else:
        raise AssertionError("tampered truth metadata must fail closed")
