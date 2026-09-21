from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_a4_current_rf40_axial_shutdown_radial_stress_independent_audit as audit


FORBIDDEN_API_NAMES = {
    "threshold",
    "gate",
    "seed",
    "resolution",
    "pressure",
    "forcing",
    "viscosity",
    "residual",
    "correction",
}


def _mean_arrays(radii: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.asarray(radii, dtype=float)
    quadratic = np.column_stack(
        (
            0.3 * r + 0.1 * r**2,
            0.7 * r**2 - 0.2 * r**3,
            -0.4 * r + 0.8 * r**2,
        )
    )
    mixed = np.column_stack(
        (
            -0.2 * r + 0.05 * r**3,
            0.15 * r + 0.4 * r**2,
            0.25 * r - 0.3 * r**3,
        )
    )
    return quadratic, mixed, quadratic + mixed


def _piece(radii: np.ndarray, source: np.ndarray, exponent: int) -> dict[str, object]:
    independent = audit._independent_compact_stress(
        radii,
        source,
        exponent=exponent,
        bump_center=audit.FROZEN_BUMP_CENTER,
        bump_halfwidth=audit.FROZEN_BUMP_HALFWIDTH,
    )
    stress = np.asarray(independent["stress"], dtype=float)
    return {
        "weighted_moment": float(independent["weighted_moment"]),
        "moment_complement_weighted_moment": 0.0,
        "stress": stress.tolist(),
        "stress_rms": float(independent["stress_rms"]),
        "stress_outer_edge": float(stress[-1]),
    }


def _receipt(radii: np.ndarray) -> dict[str, object]:
    q, m, n = _mean_arrays(radii)
    sources = {"quadratic": q, "mixed": m, "aggregate": n}
    mean = {
        "mean_oscillatory_self_advection_cylindrical": q.tolist(),
        "mean_mixed_cross_cylindrical": m.tolist(),
        "mean_aggregate_nonlinear_cylindrical": n.tolist(),
        "quadratic_mean_rms": float(np.sqrt(np.mean(q * q))),
    }
    theta = {
        piece: _piece(radii, values[:, 1], 2)
        for piece, values in sources.items()
    }
    axial = {
        piece: _piece(radii, values[:, 2], 1)
        for piece, values in sources.items()
    }
    return {
        "schema": audit.PARENT_SCHEMA,
        "geometry": {
            "time": audit.FROZEN_TIME,
            "axial_z": audit.FROZEN_Z,
            "radii": radii.tolist(),
            "bump_center": audit.FROZEN_BUMP_CENTER,
            "bump_halfwidth": audit.FROZEN_BUMP_HALFWIDTH,
        },
        "current_rf40_axial_shutdown_mean_witness": mean,
        "theta_e2": theta,
        "axial_e1": axial,
        "parent_agent3_1011_source_blob": audit.PARENT_AGENT3_1011_SOURCE_BLOB,
        "provenance": {
            "parent_agent3_head": "6aabf7b6dd28ae683a3776d518e9a977507911e9",
            "parent_agent3_source_blob": audit.PARENT_AGENT3_1011_SOURCE_BLOB,
            "radial_operator_source_blob": audit.RADIAL_OPERATOR_SOURCE_BLOB,
            "agent2_composite_head": audit.AGENT2_999_HEAD,
            "agent2_composite_source_blob": audit.AGENT2_999_SOURCE_BLOB,
            "agent1_leading_head": audit.AGENT1_993_HEAD,
            "agent1_leading_source_blob": audit.AGENT1_993_SOURCE_BLOB,
        },
        "truth_boundary": {
            "complete_ns_defect": False,
            "scoped_current_RF40_axial_shutdown_radial_stress_authorized_as_correction_target": False,
            "current_real_ns_correction_velocity_materialized": False,
            "real_candidate_finite_correction_cycle_run": False,
            "heldout_normalized_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
            "final_normalized_momentum_gate": 1e-3,
            "final_normalized_divergence_gate": 1e-5,
        },
    }


def _receipts() -> list[dict[str, object]]:
    return [_receipt(grid) for grid in audit.frozen_radial_grids()]


def test_frozen_protocol_and_public_api_have_no_tuning_knobs() -> None:
    assert audit.FROZEN_SEED == 9173751
    assert audit.FROZEN_GRID_COUNTS == (43, 85, 169)
    assert audit.FINE_STRESS_REL_RMS_GATE == 5e-2
    assert audit.FINE_STRESS_REL_MAX_GATE == 1.5e-1
    assert audit.FINAL_NORMALIZED_MOMENTUM_GATE == 1e-3
    assert audit.FINAL_NORMALIZED_DIVERGENCE_GATE == 1e-5
    for fn in (audit.audit_serialized_receipts, audit.enforce_scientific_gates):
        names = set(inspect.signature(fn).parameters)
        assert names.isdisjoint(FORBIDDEN_API_NAMES)


def test_offgrid_is_fresh_shifted_fine_grid() -> None:
    coarse, medium, fine, shifted = audit.frozen_radial_grids()
    assert (coarse.size, medium.size, fine.size, shifted.size) == (43, 85, 169, 169)
    assert np.all(np.diff(shifted) > 0.0)
    assert shifted[0] == pytest.approx(fine[0], abs=0.0)
    assert shifted[-1] == pytest.approx(fine[-1], abs=0.0)
    assert not np.array_equal(shifted[1:-1], fine[1:-1])


@pytest.mark.parametrize("exponent", [1, 2])
def test_cubic_gl8_weighted_integral_is_exact_on_cubic_manufactured_data(
    exponent: int,
) -> None:
    r = np.linspace(0.005, 0.44, 85)
    values = 1.2 - 0.7 * r + 0.3 * r**2 - 0.2 * r**3
    cumulative = audit._piecewise_cubic_weighted_cumulative(
        r, values, exponent=exponent
    )
    R = r[-1]
    expected = (
        1.2 * R ** (exponent + 1) / (exponent + 1)
        - 0.7 * R ** (exponent + 2) / (exponent + 2)
        + 0.3 * R ** (exponent + 3) / (exponent + 3)
        - 0.2 * R ** (exponent + 4) / (exponent + 4)
    )
    assert cumulative[-1] == pytest.approx(expected, rel=2e-12, abs=2e-14)


def test_manufactured_serialized_receipts_pass_scoped_contract() -> None:
    report = audit.audit_serialized_receipts(_receipts())
    assert report["audit_pass"] is True
    assert report["failures"] == []
    audit.enforce_scientific_gates(report)
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["after_correction_ns_residual_assessed"] is False


def test_endpoint_preserving_interior_stress_mutation_is_detected() -> None:
    receipts = _receipts()
    mutated = copy.deepcopy(receipts)
    radii = np.asarray(mutated[2]["geometry"]["radii"], dtype=float)
    stress = np.asarray(mutated[2]["theta_e2"]["aggregate"]["stress"], dtype=float)
    scale = max(float(np.max(np.abs(stress))), 1e-6)
    tangent = 0.35 * scale * np.sin(
        np.pi * (radii - radii[0]) / (radii[-1] - radii[0])
    )
    assert tangent[0] == pytest.approx(0.0, abs=1e-15)
    assert tangent[-1] == pytest.approx(0.0, abs=1e-15)
    mutated[2]["theta_e2"]["aggregate"]["stress"] = (stress + tangent).tolist()
    report = audit.audit_serialized_receipts(mutated)
    assert report["audit_pass"] is False
    with pytest.raises(RuntimeError):
        audit.enforce_scientific_gates(report)


def test_bump_geometry_drift_fails_before_comparison() -> None:
    receipts = _receipts()
    receipts[2]["geometry"]["bump_center"] = audit.FROZEN_BUMP_CENTER + 0.005
    with pytest.raises(ValueError, match="compact-bump geometry drifted"):
        audit.audit_serialized_receipts(receipts)


def test_mean_attribution_mutation_fails_before_stress_comparison() -> None:
    receipts = _receipts()
    mean = receipts[2]["current_rf40_axial_shutdown_mean_witness"]
    aggregate = np.asarray(mean["mean_aggregate_nonlinear_cylindrical"], dtype=float)
    aggregate[20, 1] += 1e-4
    mean["mean_aggregate_nonlinear_cylindrical"] = aggregate.tolist()
    with pytest.raises(ValueError, match="mean attribution no longer closes"):
        audit.audit_serialized_receipts(receipts)


@pytest.mark.parametrize(
    ("where", "key", "value"),
    [
        ("top", "parent_agent3_1011_source_blob", "0" * 40),
        ("provenance", "agent2_composite_head", "0" * 40),
        ("provenance", "agent1_leading_source_blob", "0" * 40),
        ("truth", "pde_validated", True),
    ],
)
def test_provenance_and_truth_mutations_fail_closed(
    where: str, key: str, value: object
) -> None:
    receipts = _receipts()
    target = receipts[2]
    if where == "provenance":
        target = target["provenance"]
    elif where == "truth":
        target = target["truth_boundary"]
    target[key] = value
    with pytest.raises(ValueError):
        audit.audit_serialized_receipts(receipts)
