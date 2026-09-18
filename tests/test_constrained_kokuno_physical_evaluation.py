import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_physical_evaluation import (
    SOURCE_B_G,
    SOURCE_J_G,
    SOURCE_KAPPA_S,
    SOURCE_LAMBDA_G,
    SOURCE_T_G,
    SOURCE_V_R,
    SOURCE_V_T,
    KokunoPhysicalEvaluationMap,
    corrected_source_torus_constants,
)


def _fd4(fun, x, h):
    return (-fun(x + 2 * h) + 8 * fun(x + h) - 8 * fun(x - h) + fun(x - 2 * h)) / (12 * h)


def test_displayed_source_torus_constants_and_group_identities():
    source = corrected_source_torus_constants()
    sqrt2 = np.sqrt(2.0)
    assert SOURCE_B_G == pytest.approx(sqrt2 - 1.0)
    assert SOURCE_V_R == pytest.approx((1.0, 1.0 - sqrt2))
    assert SOURCE_V_T == pytest.approx((sqrt2 - 1.0, 1.0))
    assert SOURCE_LAMBDA_G == pytest.approx(4.0 - sqrt2)
    assert SOURCE_T_G == pytest.approx(4.0 + sqrt2)
    assert SOURCE_KAPPA_S == pytest.approx(1.0e-5)

    jg = np.asarray(SOURCE_J_G)
    vr = np.asarray(SOURCE_V_R)
    vt = np.asarray(SOURCE_V_T)
    np.testing.assert_allclose(jg @ vr, SOURCE_LAMBDA_G * vr, rtol=0, atol=2e-15)
    np.testing.assert_allclose(jg @ vt, SOURCE_T_G * vt, rtol=0, atol=2e-15)
    assert vr @ vt == pytest.approx(0.0, abs=2e-16)
    assert np.linalg.det(jg) == pytest.approx(14.0)
    diag = source["identity_diagnostics"]
    assert diag["Jg_vr_minus_Lambda_vr_max_abs"] < 2e-15
    assert diag["Jg_vt_minus_T_vt_max_abs"] < 2e-15
    assert abs(diag["vr_dot_vt"]) < 2e-16
    assert diag["det_J_g"] == pytest.approx(14.0)


def test_source_bound_constructor_uses_displayed_constants_and_radial_exponent():
    h = 0.004
    op = KokunoPhysicalEvaluationMap.from_corrected_source_constants(h=h)
    expected_rho = np.log(4.0 - np.sqrt(2.0)) / np.log(4.0 + np.sqrt(2.0))
    expected_dr = 2.0 * ((1.0 + h) * expected_rho - h * 1.0e-5)
    assert op.v_r == pytest.approx(SOURCE_V_R)
    assert op.v_t == pytest.approx(SOURCE_V_T)
    assert op.d_r == pytest.approx(expected_dr)
    assert op.binding == "corrected_source_displayed_constants"
    p = op.provenance()
    assert p["numeric_binding"]["uses_displayed_source_torus_constants"] is True
    assert p["source_displayed_torus_constants"]["Lambda_g"] == pytest.approx(SOURCE_LAMBDA_G)
    assert p["source_displayed_torus_constants"]["T_g"] == pytest.approx(SOURCE_T_G)


def test_torus_phase_and_source_radial_exponent():
    out = KokunoPhysicalEvaluationMap.source_radial_exponent(14.0, 9.0, 0.004)
    expected = 2 * ((1.004) * np.log(14.0) / np.log(9.0) - 0.004e-5)
    assert out["d_r"] == pytest.approx(expected)
    op = KokunoPhysicalEvaluationMap((0.7, -0.2), (0.4, 0.9), out["d_r"], 0.004)
    y = op.torus_phase(np.array([0.8, 1.1]), np.array([0.2, -0.1]))
    assert y.shape == (2, 2)
    assert np.all((y >= 0) & (y < 1))
    assert op.binding == "caller_supplied"


def test_chain_rule_matches_independent_fd4():
    op = KokunoPhysicalEvaluationMap((0.31, -0.27), (0.23, 0.41), 1.37, 0.004)
    r, z, t, q = 0.83, -0.22, 0.17, 0.61
    twopi = 2 * np.pi

    def lifted(rr, zz, tt):
        y = op.torus_phase(rr, tt)
        return rr**2 + zz * tt + np.sin(twopi * y[..., 0]) + 0.3 * np.cos(twopi * y[..., 1])

    y = op.torus_phase(r, t)
    gy = np.array(
        [twopi * np.cos(twopi * y[0]), -0.3 * twopi * np.sin(twopi * y[1])]
    )
    out = op.differentiate(
        r=r,
        q=q,
        partial_r=2 * r,
        partial_t=z,
        partial_z=t,
        grad_y=gy,
    )
    hr, ht, hz = 2e-4, 2e-4, 2e-4
    dr_fd = _fd4(lambda x: lifted(x, z, t), r, hr)
    dt_fd = _fd4(lambda x: lifted(r, z, x), t, ht)
    dz_fd = _fd4(lambda x: lifted(r, x, t), z, hz)
    assert out["mathsf_r"] == pytest.approx(dr_fd, rel=2e-9, abs=2e-9)
    assert out["mathsf_t"] == pytest.approx(dt_fd, rel=2e-9, abs=2e-9)
    assert out["partial_z_after_evaluation"] == pytest.approx(dz_fd, rel=2e-10, abs=2e-10)
    assert out["D_r"] == pytest.approx(np.sqrt(q) * dr_fd, rel=2e-9)
    assert out["D_z"] == pytest.approx(np.sqrt(q) * dz_fd, rel=2e-9)
    assert out["mathsf_t_star"] == pytest.approx(q ** (1 + op.h) * dt_fd, rel=2e-9)


def test_batch_and_fail_closed_guards():
    op = KokunoPhysicalEvaluationMap((0.2, 0.3), (0.5, -0.4), 1.2, 0.003)
    r = np.array([0.5, 0.8, 1.1])
    q = np.array([0.7, 0.8, 0.9])
    gy = np.stack([r, -2 * r], axis=-1)
    out = op.differentiate(
        r=r,
        q=q,
        partial_r=0,
        partial_t=1,
        partial_z=2,
        grad_y=gy,
    )
    assert out["D_r"].shape == (3,)
    assert np.all(np.isfinite(out["D_r"]))
    with pytest.raises(ValueError, match="r>0"):
        op.torus_phase(0.0, 0.0)
    with pytest.raises(ValueError, match="strictly positive"):
        op.differentiate(
            r=r,
            q=np.array([0.7, 0.0, 0.9]),
            partial_r=0,
            partial_t=0,
            partial_z=0,
            grad_y=gy,
        )
    with pytest.raises(ValueError, match="sample_shape"):
        op.differentiate(
            r=r,
            q=q,
            partial_r=0,
            partial_t=0,
            partial_z=0,
            grad_y=np.zeros((4, 2)),
        )
    with pytest.raises(ValueError, match="binding"):
        KokunoPhysicalEvaluationMap((0.2, 0.3), (0.5, -0.4), 1.2, 0.003, binding="paper_exact")


def test_provenance_keeps_truth_boundary():
    op = KokunoPhysicalEvaluationMap.from_corrected_source_constants(h=0.003)
    p = op.provenance()
    assert p["source"]["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert p["truth_boundary"]["auxiliary_torus_physical_evaluation_map_executable"] is True
    assert p["truth_boundary"]["source_torus_vectors_recovered"] is True
    assert p["truth_boundary"]["source_group_constants_recovered"] is True
    assert p["truth_boundary"]["actual_positive_order_background_bound"] is False
    assert p["truth_boundary"]["actual_auxiliary_torus_mode_family_bound"] is False
    assert p["truth_boundary"]["public_xyz_t_velocity_correction_materialized"] is False
    assert p["truth_boundary"]["paper_exact"] is False
    assert len(p["sha256"]) == 64
