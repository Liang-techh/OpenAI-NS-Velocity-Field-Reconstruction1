from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_partial_radial_stress_independent_audit import (
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    FROZEN_BUMP_CENTER,
    FROZEN_BUMP_HALFWIDTH,
    FROZEN_TIME,
    FROZEN_Z,
    _independent_compact_stress,
    _piecewise_cubic_weighted_cumulative,
    audit_serialized_receipts,
    enforce_audit,
    frozen_radial_grids,
    mutate_bump_geometry,
    mutate_fine_stress_tangent,
    truth_boundary,
)


def _piece_report(r: np.ndarray, source: np.ndarray, exponent: int) -> dict[str, object]:
    independent = _independent_compact_stress(
        r,
        source,
        exponent=exponent,
        bump_center=FROZEN_BUMP_CENTER,
        bump_halfwidth=FROZEN_BUMP_HALFWIDTH,
    )
    stress = np.asarray(independent["stress"], dtype=float)
    return {
        "weighted_moment": float(independent["weighted_moment"]),
        "stress": stress.tolist(),
        "stress_rms": float(np.sqrt(np.mean(stress * stress))),
        "stress_outer_edge": float(stress[-1]),
    }


def _manufactured_receipt(r: np.ndarray) -> dict[str, object]:
    # Smooth nontrivial source with exact piece attribution.  The same fixture is
    # only a validator mechanics test; the scientific workflow consumes #976.
    q = (0.7 + 0.4 * r + 0.2 * r * r) * np.exp(-4.0 * r)
    mixed = -0.35 * q + 0.03 * r**3
    aggregate = q + mixed

    def cyl(source: np.ndarray) -> np.ndarray:
        return np.column_stack((0.15 * source, 0.8 * source, -0.6 * source))

    q_cyl = cyl(q)
    m_cyl = cyl(mixed)
    a_cyl = cyl(aggregate)
    mean = {
        "mean_oscillatory_self_advection_cylindrical": q_cyl.tolist(),
        "mean_mixed_cross_cylindrical": m_cyl.tolist(),
        "mean_aggregate_nonlinear_cylindrical": a_cyl.tolist(),
    }

    receipt: dict[str, object] = {
        "schema": "kokuno-a3-current-partial-nonlinear-radial-stress-v1",
        "parent_agent3_971_source_blob": "9c62ae3c3b4b3ca2d9082f822d57feef700ea6cf",
        "geometry": {
            "time": FROZEN_TIME,
            "axial_z": FROZEN_Z,
            "radii": r.tolist(),
            "bump_center": FROZEN_BUMP_CENTER,
            "bump_halfwidth": FROZEN_BUMP_HALFWIDTH,
        },
        "current_partial_mean_witness": mean,
        "truth_boundary": {
            "complete_ns_defect": False,
            "scoped_current_partial_radial_stress_authorized_as_correction_target": False,
            "current_real_ns_correction_velocity_materialized": False,
            "heldout_normalized_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
        },
    }
    for channel, exponent, component in (("theta_e2", 2, 1), ("axial_e1", 1, 2)):
        receipt[channel] = {
            "quadratic": _piece_report(r, q_cyl[:, component], exponent),
            "mixed": _piece_report(r, m_cyl[:, component], exponent),
            "aggregate": _piece_report(r, a_cyl[:, component], exponent),
            "stress_piece_closure_relative_max": 0.0,
        }
    return receipt


def test_piecewise_cubic_gl_integral_is_exact_for_quadratic_source():
    r = np.linspace(0.013, 0.47, 37)
    f = 1.2 - 0.7 * r + 0.4 * r * r
    for exponent in (1, 2):
        got = _piecewise_cubic_weighted_cumulative(r, f, exponent=exponent)
        expected = (
            1.2 * r ** (exponent + 1) / (exponent + 1)
            - 0.7 * r ** (exponent + 2) / (exponent + 2)
            + 0.4 * r ** (exponent + 3) / (exponent + 3)
        )
        assert np.max(np.abs(got - expected)) <= 2.0e-12


def test_frozen_grids_are_nested_plus_deterministic_shifted_offgrid():
    grids = frozen_radial_grids()
    assert [len(v) for v in grids] == [43, 85, 169, 169]
    assert all(np.all(np.diff(v) > 0.0) for v in grids)
    assert grids[0][0] == grids[1][0] == grids[2][0] == grids[3][0]
    assert grids[0][-1] == grids[1][-1] == grids[2][-1] == grids[3][-1]
    assert not np.array_equal(grids[2][1:-1], grids[3][1:-1])
    assert np.array_equal(grids[3], frozen_radial_grids()[3])


def test_manufactured_serialized_receipts_pass_scoped_audit():
    receipts = [_manufactured_receipt(r) for r in frozen_radial_grids()]
    report = audit_serialized_receipts(copy.deepcopy(receipts))
    enforce_audit(report)
    assert report["candidate_residual_evidence"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["axis_evidence_scope"].startswith("positive-radius")


def test_outer_edge_preserving_stress_mutation_is_detected():
    receipts = [_manufactured_receipt(r) for r in frozen_radial_grids()]
    mutated = copy.deepcopy(receipts)
    mutated[2] = mutate_fine_stress_tangent(mutated[2])
    report = audit_serialized_receipts(mutated)
    with pytest.raises(AssertionError, match="fine stress"):
        enforce_audit(report)


def test_bump_parameter_drift_fails_closed_before_scientific_gate():
    receipts = [_manufactured_receipt(r) for r in frozen_radial_grids()]
    mutated = copy.deepcopy(receipts)
    mutated[3] = mutate_bump_geometry(mutated[3])
    with pytest.raises(ValueError, match="bump geometry drifted"):
        audit_serialized_receipts(mutated)


def test_mean_piece_mutation_is_detected_before_stress_comparison():
    receipts = [_manufactured_receipt(r) for r in frozen_radial_grids()]
    mutated = copy.deepcopy(receipts)
    arr = np.asarray(
        mutated[1]["current_partial_mean_witness"]["mean_aggregate_nonlinear_cylindrical"],
        dtype=float,
    )
    arr[:, 1] += 1.0e-4
    mutated[1]["current_partial_mean_witness"]["mean_aggregate_nonlinear_cylindrical"] = arr.tolist()
    with pytest.raises(ValueError, match="mean attribution"):
        audit_serialized_receipts(mutated)


def test_public_audit_has_no_threshold_or_physics_tuning_inputs():
    params = inspect.signature(audit_serialized_receipts).parameters
    assert list(params) == ["receipts"]
    for forbidden in (
        "threshold",
        "gate",
        "seed",
        "resolution",
        "forcing",
        "pressure",
        "viscosity",
        "nu",
        "residual",
        "correction",
    ):
        assert forbidden not in params


def test_truth_boundary_keeps_complete_ns_and_pde_claims_closed():
    boundary = truth_boundary()
    assert boundary["independent_reference_uses_parent_radial_inverse"] is False
    assert boundary["independent_reference_uses_public_serialized_mean_samples"] is True
    assert boundary["three_resolution_radial_stress_audit"] is True
    assert boundary["shifted_offgrid_radial_stress_audit"] is True
    assert boundary["exact_cartesian_axis_audit"] is False
    for key in (
        "complete_ns_defect",
        "scoped_radial_stress_authorized_as_correction_target",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "current_real_ns_correction_velocity_materialized",
        "after_correction_residual_assessed",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
    ):
        assert boundary[key] is False, key
    assert boundary["final_normalized_momentum_gate"] == FINAL_NORMALIZED_MOMENTUM_GATE == 1.0e-3
    assert boundary["final_normalized_divergence_gate"] == FINAL_NORMALIZED_DIVERGENCE_GATE == 1.0e-5
