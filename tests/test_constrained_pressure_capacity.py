import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pressure_capacity import (
    diagnose_pressure_family_capacity,
)


class LinearPressureCandidate:
    def __init__(self, coefficient_count=3, coefficients=None, basis_kind="xyz"):
        self.pressure_coefficients = tuple(
            np.zeros(coefficient_count) if coefficients is None else coefficients
        )
        self.basis_kind = basis_kind

    def velocity(self, points, time):
        return np.zeros_like(np.asarray(points, dtype=float))

    def pressure_basis(self, points, time):
        p = np.asarray(points, dtype=float)
        if self.basis_kind == "xyz":
            return p[:, :len(self.pressure_coefficients)]
        if self.basis_kind == "x_only":
            return p[:, :1]
        if self.basis_kind == "constant":
            return np.ones((len(p), 1))
        raise AssertionError("bad test basis")

    def pressure(self, points, time):
        return self.pressure_basis(points, time) @ np.asarray(self.pressure_coefficients)


def constant_force(vector):
    vector = np.asarray(vector, dtype=float)
    def force(points, time):
        return np.broadcast_to(vector, (len(points), 3)).copy()
    return force


def test_pressure_capacity_exactly_recovers_representable_constant_gradient():
    c = LinearPressureCandidate()
    points = np.array([[0.1,0.2,0.3],[-0.2,0.1,-0.1],[0.4,-0.3,0.2]])
    report = diagnose_pressure_family_capacity(
        c, constant_force([0.25,-0.5,0.75]), points, 0.5, step=0.01
    )
    assert report.design_rank == 3
    assert report.design_nullity == 0
    assert report.capacity_coefficients == pytest.approx((0.25,-0.5,0.75), abs=1e-10)
    assert report.capacity_momentum_rms < 1e-10
    assert report.recoverable_fraction_rms > 1-1e-10
    assert report.same_sample_capacity is True
    assert report.pde_validated is False


def test_pressure_family_cannot_hide_orthogonal_residual():
    c = LinearPressureCandidate(1, basis_kind="x_only")
    points = np.array([[0.1,0.2,0.3],[-0.2,0.1,-0.1]])
    report = diagnose_pressure_family_capacity(
        c, constant_force([0.0,1.0,0.0]), points, 0.5, step=0.01
    )
    assert report.design_rank == 1
    assert report.capacity_momentum_rms == pytest.approx(1.0, abs=1e-10)
    assert report.recoverable_fraction_rms == pytest.approx(0.0, abs=1e-10)


def test_pressure_bounds_are_enforced_and_reported():
    c = LinearPressureCandidate(1, basis_kind="x_only")
    points = np.array([[0.1,0.2,0.3],[-0.2,0.1,-0.1]])
    report = diagnose_pressure_family_capacity(
        c, constant_force([2.0,0.0,0.0]), points, 0.5, step=0.01
    )
    assert report.capacity_coefficients[0] == pytest.approx(1.0, abs=1e-8)
    assert report.capacity_momentum_rms == pytest.approx(1.0, abs=1e-8)
    assert report.active_upper_bounds == (0,)


def test_gauge_like_constant_column_is_reported_as_nullity():
    c = LinearPressureCandidate(1, basis_kind="constant")
    points = np.array([[0.1,0.2,0.3],[-0.2,0.1,-0.1]])
    report = diagnose_pressure_family_capacity(
        c, constant_force([0.0,1.0,0.0]), points, 0.5, step=0.01
    )
    assert report.design_rank == 0
    assert report.design_nullity == 1
    assert report.design_condition is None
    assert report.capacity_momentum_rms == pytest.approx(1.0, abs=1e-10)


def test_fail_closed_on_out_of_contract_candidate_or_basis():
    points = np.array([[0.1,0.2,0.3]])
    bad = LinearPressureCandidate(1, coefficients=[1.1], basis_kind="x_only")
    with pytest.raises(ValueError, match=r"\[-1,1\]"):
        diagnose_pressure_family_capacity(
            bad, constant_force([0,0,0]), points, 0.5, step=0.01
        )

    class BadBasis(LinearPressureCandidate):
        def pressure_basis(self, points, time):
            return np.full((len(points), 1), np.nan)

    with pytest.raises(ValueError, match="pressure_basis"):
        diagnose_pressure_family_capacity(
            BadBasis(1), constant_force([0,0,0]), points, 0.5, step=0.01
        )
