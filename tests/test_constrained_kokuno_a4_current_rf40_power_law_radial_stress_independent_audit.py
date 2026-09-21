from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_current_rf40_power_law_radial_stress_independent_audit import (
    AGENT1_1005_HEAD,
    AGENT1_1005_SOURCE_BLOB,
    AGENT2_1010_HEAD,
    AGENT2_1010_SOURCE_BLOB,
    AGENT2_960_HEAD,
    AGENT2_960_SOURCE_BLOB,
    FROZEN_BUMP_CENTER,
    FROZEN_BUMP_HALFWIDTH,
    FROZEN_TIME,
    FROZEN_Z,
    FINAL_NORMALIZED_DIVERGENCE_GATE,
    FINAL_NORMALIZED_MOMENTUM_GATE,
    PARENT_AGENT3_1028_HEAD,
    PARENT_AGENT3_1028_SOURCE_BLOB,
    PARENT_SCHEMA,
    RADIAL_OPERATOR_SOURCE_BLOB,
    _independent_compact_stress,
    audit_serialized_receipts,
    enforce_scientific_gates,
    frozen_radial_grids,
)


def _sources(r: np.ndarray) -> dict[str, np.ndarray]:
    quadratic = np.column_stack(
        (
            0.08 * r + 0.03 * r**2,
            0.70 + 0.20 * r - 0.10 * r**2,
            -0.45 + 0.30 * r + 0.05 * r**3,
        )
    )
    mixed = np.column_stack(
        (
            -0.02 + 0.04 * r,
            -0.25 + 0.08 * r + 0.03 * r**3,
            0.16 - 0.07 * r + 0.02 * r**2,
        )
    )
    return {
        "quadratic": quadratic,
        "mixed": mixed,
        "aggregate": quadratic + mixed,
    }


def _receipt(r: np.ndarray) -> dict[str, object]:
    sources = _sources(r)
    mean = {
        "mean_oscillatory_self_advection_cylindrical": sources["quadratic"].tolist(),
        "mean_mixed_cross_cylindrical": sources["mixed"].tolist(),
        "mean_aggregate_nonlinear_cylindrical": sources["aggregate"].tolist(),
    }
    out: dict[str, object] = {
        "schema": PARENT_SCHEMA,
        "geometry": {
            "time": FROZEN_TIME,
            "axial_z": FROZEN_Z,
            "radii": r.tolist(),
            "bump_center": FROZEN_BUMP_CENTER,
            "bump_halfwidth": FROZEN_BUMP_HALFWIDTH,
        },
        "current_rf40_power_law_mean_witness": mean,
        "parent_agent3_1028_source_blob": PARENT_AGENT3_1028_SOURCE_BLOB,
        "provenance": {
            "parent_agent3_head": PARENT_AGENT3_1028_HEAD,
            "radial_operator_source_blob": RADIAL_OPERATOR_SOURCE_BLOB,
            "agent2_composite_head": AGENT2_1010_HEAD,
            "agent2_composite_source_blob": AGENT2_1010_SOURCE_BLOB,
            "agent2_differential_head": AGENT2_960_HEAD,
            "agent2_differential_source_blob": AGENT2_960_SOURCE_BLOB,
            "agent1_leading_head": AGENT1_1005_HEAD,
            "agent1_leading_source_blob": AGENT1_1005_SOURCE_BLOB,
        },
        "truth_boundary": {
            "complete_ns_defect": False,
            "scoped_current_RF40_power_law_radial_stress_authorized_as_correction_target": False,
            "current_real_ns_correction_velocity_materialized": False,
            "real_candidate_finite_correction_cycle_run": False,
            "heldout_normalized_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "same_protocol_comparable_to_st006": False,
            "pde_validated": False,
            "final_normalized_momentum_gate": FINAL_NORMALIZED_MOMENTUM_GATE,
            "final_normalized_divergence_gate": FINAL_NORMALIZED_DIVERGENCE_GATE,
        },
    }
    for channel, exponent, component in (("theta_e2", 2, 1), ("axial_e1", 1, 2)):
        pieces: dict[str, object] = {}
        for piece in ("quadratic", "mixed", "aggregate"):
            independent = _independent_compact_stress(
                r,
                sources[piece][:, component],
                exponent=exponent,
                bump_center=FROZEN_BUMP_CENTER,
                bump_halfwidth=FROZEN_BUMP_HALFWIDTH,
            )
            pieces[piece] = {
                "stress": independent["stress"].tolist(),
                "weighted_moment": independent["weighted_moment"],
            }
        out[channel] = pieces
    return out


def _all_receipts() -> list[dict[str, object]]:
    return [_receipt(grid) for grid in frozen_radial_grids()]


def test_frozen_grids_have_three_resolutions_plus_distinct_offgrid() -> None:
    coarse, medium, fine, offgrid = frozen_radial_grids()
    assert tuple(map(len, (coarse, medium, fine))) == (43, 85, 169)
    assert len(offgrid) == 169
    assert np.all(np.diff(offgrid) > 0.0)
    assert np.array_equal(offgrid[[0, -1]], fine[[0, -1]])
    assert not np.array_equal(offgrid[1:-1], fine[1:-1])


def test_synthetic_public_receipts_pass_exact_independent_reference() -> None:
    report = audit_serialized_receipts(_all_receipts())
    assert report["audit_pass"] is True
    assert report["failures"] == []
    assert report["fine_max_stress_relative_rms"] == pytest.approx(0.0)
    assert report["offgrid_max_stress_relative_rms"] == pytest.approx(0.0)
    enforce_scientific_gates(report)
    assert report["truth_boundary"]["pde_validated"] is False


def test_endpoint_preserving_interior_stress_mutation_is_detected() -> None:
    receipts = _all_receipts()
    r = np.asarray(receipts[2]["geometry"]["radii"], dtype=float)
    stress = np.asarray(receipts[2]["theta_e2"]["aggregate"]["stress"], dtype=float)
    scale = max(float(np.max(np.abs(stress))), 1e-12)
    tangent = 0.35 * scale * np.sin(np.pi * (r - r[0]) / (r[-1] - r[0]))
    receipts[2]["theta_e2"]["aggregate"]["stress"] = (stress + tangent).tolist()
    report = audit_serialized_receipts(receipts)
    assert report["audit_pass"] is False
    with pytest.raises(RuntimeError):
        enforce_scientific_gates(report)


def test_bump_geometry_drift_fails_closed() -> None:
    receipts = _all_receipts()
    receipts[2]["geometry"]["bump_center"] = FROZEN_BUMP_CENTER + 0.005
    with pytest.raises(ValueError, match="geometry drifted"):
        audit_serialized_receipts(receipts)


def test_mean_attribution_drift_fails_closed() -> None:
    receipts = _all_receipts()
    mean = receipts[2]["current_rf40_power_law_mean_witness"]
    aggregate = np.asarray(mean["mean_aggregate_nonlinear_cylindrical"], dtype=float)
    aggregate[len(aggregate) // 2, 1] += 1e-4
    mean["mean_aggregate_nonlinear_cylindrical"] = aggregate.tolist()
    with pytest.raises(ValueError, match="no longer closes"):
        audit_serialized_receipts(receipts)


@pytest.mark.parametrize("mutation", ["provenance", "truth"])
def test_identity_and_truth_mutations_fail_closed(mutation: str) -> None:
    receipts = _all_receipts()
    if mutation == "provenance":
        receipts[2]["provenance"]["agent2_composite_head"] = "0" * 40
        message = "provenance drifted"
    else:
        receipts[2]["truth_boundary"]["pde_validated"] = True
        message = "scientific boundary drifted"
    with pytest.raises(ValueError, match=message):
        audit_serialized_receipts(receipts)
