import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_physical_evaluation import (
    KokunoPhysicalEvaluationMap,
)


def _fd4(fun, x, h):
    return (-fun(x + 2 * h) + 8 * fun(x + h) - 8 * fun(x - h) + fun(x - 2 * h)) / (12 * h)


def test_torus_phase_and_source_radial_exponent():
    out = KokunoPhysicalEvaluationMap.source_radial_exponent(14.0, 9.0, 0.004)
    expected = 2 * ((1.004) * np.log(14.0) / np.log(9.0) - 0.004e-5)
    assert out["d_r"] == pytest.approx(expected)
    op = KokunoPhysicalEvaluationMap((0.7, -0.2), (0.4, 0.9), out["d_r"], 0.004)
    y = op.torus_phase(np.array([0.8, 1.1]), np.array([0.2, -0.1]))
    assert y.shape == (2, 2)
    assert np.all((y >= 0) & (y < 1))


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


def test_provenance_keeps_truth_boundary():
    op = KokunoPhysicalEvaluationMap((0.2, 0.3), (0.5, -0.4), 1.2, 0.003)
    p = op.provenance()
    assert p["source"]["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert p["truth_boundary"]["auxiliary_torus_physical_evaluation_map_executable"] is True
    assert p["truth_boundary"]["source_torus_vectors_recovered"] is False
    assert p["truth_boundary"]["actual_positive_order_background_bound"] is False
    assert p["truth_boundary"]["public_xyz_t_velocity_correction_materialized"] is False
    assert p["truth_boundary"]["paper_exact"] is False
    assert len(p["sha256"]) == 64
