from __future__ import annotations

import copy
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction import kokuno_a4_current_i4_radial_stress_independent_audit as audit


def _manufactured_receipt(radii: np.ndarray) -> dict:
    r = np.asarray(radii, dtype=float)
    quadratic = np.column_stack(
        (0.1 * r, 0.4 + 0.2 * r + 0.1 * r**2, -0.3 + 0.15 * r)
    )
    mixed = np.column_stack(
        (-0.05 * r, -0.1 + 0.3 * r**2, 0.2 - 0.07 * r**2)
    )
    aggregate = quadratic + mixed
    sources = {"quadratic": quadratic, "mixed": mixed, "aggregate": aggregate}
    receipt = {
        "schema": audit.PARENT_SCHEMA,
        "parent_agent3_1101_source_blob": audit.PARENT_MEAN_SOURCE_BLOB,
        "geometry": {
            "time": audit.FROZEN_TIME,
            "axial_z": audit.FROZEN_Z,
            "radii": r.tolist(),
            "bump_center": audit.FROZEN_BUMP_CENTER,
            "bump_halfwidth": audit.FROZEN_BUMP_HALFWIDTH,
        },
        "current_i4_mean_witness": {
            "mean_oscillatory_self_advection_cylindrical": quadratic.tolist(),
            "mean_mixed_cross_cylindrical": mixed.tolist(),
            "mean_aggregate_nonlinear_cylindrical": aggregate.tolist(),
        },
        "provenance": {
            "parent_agent3_head": audit.PARENT_MEAN_HEAD,
            "radial_operator_source_blob": audit.RADIAL_OPERATOR_SOURCE_BLOB,
            "agent2_composite_head": audit.AGENT2_COMPOSITE_HEAD,
            "agent2_composite_source_blob": audit.AGENT2_COMPOSITE_SOURCE_BLOB,
            "agent2_differential_head": audit.AGENT2_DIFFERENTIAL_HEAD,
            "agent2_differential_source_blob": audit.AGENT2_DIFFERENTIAL_SOURCE_BLOB,
            "agent1_leading_head": audit.AGENT1_LEADING_HEAD,
            "agent1_leading_source_blob": audit.AGENT1_LEADING_SOURCE_BLOB,
        },
        "truth_boundary": {
            "source_I4_five_row_mean_correction_materialized": False,
            "current_I4_correction_velocity_materialized": False,
            "radial_force_materialized_in_this_increment": False,
            "complete_ns_defect": False,
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
    for channel, exponent, component in (("theta_e2", 2, 1), ("axial_e1", 1, 2)):
        receipt[channel] = {}
        for piece in ("quadratic", "mixed", "aggregate"):
            ref = audit._independent_compact_stress(
                r,
                sources[piece][:, component],
                exponent=exponent,
                bump_center=audit.FROZEN_BUMP_CENTER,
                bump_halfwidth=audit.FROZEN_BUMP_HALFWIDTH,
            )
            receipt[channel][piece] = {
                "stress": np.asarray(ref["stress"]).tolist(),
                "weighted_moment": ref["weighted_moment"],
            }
    return receipt


def _receipts() -> list[dict]:
    return [_manufactured_receipt(r) for r in audit.frozen_radial_grids()]


def test_manufactured_public_receipts_pass_frozen_protocol():
    report = audit.audit_serialized_receipts(_receipts())
    assert report["audit_pass"] is True
    assert report["failures"] == []
    audit.enforce_scientific_gates(report)
    assert report["truth_boundary"]["pde_validated"] is False


def test_endpoint_preserving_interior_stress_corruption_fails_gate():
    receipts = _receipts()
    r = np.asarray(receipts[2]["geometry"]["radii"])
    stress = np.asarray(receipts[2]["theta_e2"]["aggregate"]["stress"])
    scale = max(float(np.max(np.abs(stress))), 1e-12)
    stress = stress + 0.35 * scale * np.sin(
        np.pi * (r - r[0]) / (r[-1] - r[0])
    )
    receipts[2]["theta_e2"]["aggregate"]["stress"] = stress.tolist()
    report = audit.audit_serialized_receipts(receipts)
    assert report["audit_pass"] is False
    with pytest.raises(RuntimeError):
        audit.enforce_scientific_gates(report)


def test_mean_geometry_provenance_and_truth_mutations_fail_closed():
    receipts = _receipts()
    bad = copy.deepcopy(receipts)
    values = np.asarray(
        bad[2]["current_i4_mean_witness"]["mean_aggregate_nonlinear_cylindrical"]
    )
    values[len(values) // 2, 1] += 1e-4
    bad[2]["current_i4_mean_witness"]["mean_aggregate_nonlinear_cylindrical"] = (
        values.tolist()
    )
    with pytest.raises(ValueError):
        audit.audit_serialized_receipts(bad)

    bad = copy.deepcopy(receipts)
    bad[2]["geometry"]["bump_center"] += 0.005
    with pytest.raises(ValueError):
        audit.audit_serialized_receipts(bad)

    bad = copy.deepcopy(receipts)
    bad[2]["provenance"]["agent2_composite_head"] = "0" * 40
    with pytest.raises(ValueError):
        audit.audit_serialized_receipts(bad)

    bad = copy.deepcopy(receipts)
    bad[2]["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit.audit_serialized_receipts(bad)


def test_frozen_grids_nested_and_offgrid_distinct():
    coarse, medium, fine, offgrid = audit.frozen_radial_grids()
    assert np.allclose(coarse, fine[::4], rtol=0, atol=2e-15)
    assert np.allclose(medium, fine[::2], rtol=0, atol=2e-15)
    assert np.all(np.diff(offgrid) > 0)
    assert np.max(np.abs(offgrid[1:-1] - fine[1:-1])) > 0


def test_public_audit_surface_has_no_tuning_inputs():
    params = set(inspect.signature(audit.audit_serialized_receipts).parameters)
    assert params == {"receipts"}
    boundary = audit.truth_boundary()
    assert boundary["production_radial_operator_reused_by_reference"] is False
    assert boundary["serialized_public_receipts_only"] is True
    assert boundary["pde_validated"] is False
