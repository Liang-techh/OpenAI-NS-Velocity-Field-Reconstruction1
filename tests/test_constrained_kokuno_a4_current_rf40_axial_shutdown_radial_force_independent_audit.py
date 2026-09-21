from __future__ import annotations

import copy

import numpy as np
import pytest

import openai_ns_reconstruction.kokuno_a4_current_rf40_axial_shutdown_radial_force_independent_audit as mod


def _boundary() -> dict[str, object]:
    return {
        "complete_ns_defect": False,
        "current_real_ns_correction_velocity_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": 1e-3,
        "final_normalized_divergence_gate": 1e-5,
    }


def _stress_values(z: float, r: np.ndarray) -> dict[str, np.ndarray]:
    quadratic = (
        1.0 + 0.2 * r
        + (0.8 + 0.1 * r) * z
        + (0.2 + 0.05 * r) * z**2
        + 0.04 * z**3
    )
    mixed = (
        -0.3 + 0.1 * r
        + (-0.4 + 0.07 * r) * z
        + (0.03 - 0.02 * r) * z**2
        - 0.02 * z**3
    )
    return {
        "quadratic": quadratic,
        "mixed": mixed,
        "aggregate": quadratic + mixed,
    }


def _geometry(z: float, r: np.ndarray) -> dict[str, object]:
    return {
        "time": mod.FROZEN_TIME,
        "axial_z": float(z),
        "radii": [float(v) for v in r],
        "bump_center": mod.FROZEN_BUMP_CENTER,
        "bump_halfwidth": mod.FROZEN_BUMP_HALFWIDTH,
    }


def _stress_receipt(z: float, r: np.ndarray) -> dict[str, object]:
    values = _stress_values(z, r)
    return {
        "schema": mod.STRESS_SCHEMA,
        "geometry": _geometry(z, r),
        "parent_agent3_1011_source_blob": mod.STRESS_PARENT_AGENT3_1011_SOURCE_BLOB,
        "provenance": {
            "parent_agent3_head": "6aabf7b6dd28ae683a3776d518e9a977507911e9",
            "agent2_composite_head": mod.AGENT2_999_HEAD,
            "agent1_leading_head": mod.AGENT1_993_HEAD,
        },
        "axial_e1": {
            piece: {"stress": [float(v) for v in arr]}
            for piece, arr in values.items()
        },
        "truth_boundary": _boundary(),
    }


def _fd2(piece: str, center: float, h: float, r: np.ndarray) -> np.ndarray:
    return (
        _stress_values(center + h, r)[piece]
        - _stress_values(center - h, r)[piece]
    ) / (2.0 * h)


def _force_receipt(center: float, r: np.ndarray) -> dict[str, object]:
    out: dict[str, object] = {
        "schema": mod.PARENT_FORCE_SCHEMA,
        "geometry": _geometry(center, r),
        "parent_agent3_1015_source_blob": mod.STRESS_AGENT3_SOURCE_BLOB,
        "source_radial_force_formula": "(div T)_r = partial_z sigma_1",
        "z_derivative_step_ladder": list(mod.FROZEN_Z_STEPS),
        "piece_closure_relative_max_levels": [0.0, 0.0, 0.0],
        "aggregate_derivative_stability_preflight_passed": True,
        "provenance": {
            "parent_agent3_head": mod.STRESS_AGENT3_HEAD,
            "agent2_composite_head": mod.AGENT2_999_HEAD,
            "agent1_leading_head": mod.AGENT1_993_HEAD,
        },
        "truth_boundary": _boundary(),
    }
    for piece in ("quadratic", "mixed", "aggregate"):
        out[piece] = {
            "levels": [
                {
                    "z_step": h,
                    "radial_force": [float(v) for v in _fd2(piece, center, h, r)],
                }
                for h in mod.FROZEN_Z_STEPS
            ]
        }
    return out


def _case(center: float) -> dict[str, object]:
    r = mod.frozen_radial_grid()
    stresses = [
        _stress_receipt(center + off, r)
        for off in mod.required_z_offsets()
    ]
    return {
        "force_receipt": _force_receipt(center, r),
        "stress_receipts": stresses,
    }


def _cases() -> list[dict[str, object]]:
    return [_case(center) for center in mod.FROZEN_Z_CENTERS]


def test_frozen_protocol_is_deterministic_and_offgrid() -> None:
    r1 = mod.frozen_radial_grid()
    r2 = mod.frozen_radial_grid()
    np.testing.assert_array_equal(r1, r2)
    regular = np.linspace(mod.FROZEN_R_MIN, mod.FROZEN_R_MAX, mod.FROZEN_R_COUNT)
    assert np.any(np.abs(r1[1:-1] - regular[1:-1]) > 1e-12)
    assert np.all(np.diff(r1) > 0.0)
    assert np.count_nonzero(r1 <= mod.FROZEN_AXIS_NEAR_MAX_R) >= 2
    assert mod.required_z_offsets() == (
        -0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04
    )


def test_fd4_recovers_cubic_stress_derivative() -> None:
    r = mod.frozen_radial_grid()
    center = mod.FROZEN_Z_CENTERS[0]
    indexed = {
        off: _stress_values(center + off, r)
        for off in mod.required_z_offsets()
    }
    h = mod.FROZEN_Z_STEPS[-1]
    got = mod._fd4(indexed, "quadratic", h)
    expected = (
        (0.8 + 0.1 * r)
        + 2.0 * (0.2 + 0.05 * r) * center
        + 3.0 * 0.04 * center**2
    )
    np.testing.assert_allclose(got, expected, rtol=0.0, atol=2e-12)


def test_clean_synthetic_public_receipts_pass() -> None:
    report = mod.audit_serialized_cases(_cases())
    assert report["audit_pass"] is True
    assert report["failures"] == []
    assert report["fine_max_force_relative_rms"] < mod.FINE_FORCE_REL_RMS_GATE
    assert report["fine_max_force_relative_integral_l2"] < (
        mod.FINE_FORCE_REL_INTEGRAL_L2_GATE
    )
    assert report["truth_boundary"]["pde_validated"] is False
    mod.enforce_scientific_gates(report)


def test_case_and_shift_order_do_not_change_aggregate_metrics() -> None:
    cases = _cases()
    report_a = mod.audit_serialized_cases(cases)
    shuffled = list(reversed(copy.deepcopy(cases)))
    for case in shuffled:
        case["stress_receipts"] = list(reversed(case["stress_receipts"]))
    report_b = mod.audit_serialized_cases(shuffled)
    assert report_a["fine_max_force_relative_rms"] == pytest.approx(
        report_b["fine_max_force_relative_rms"], abs=1e-15
    )
    assert report_b["audit_pass"] is True


def test_force_sign_mutation_is_rejected_by_scoped_gates() -> None:
    cases = _cases()
    levels = cases[0]["force_receipt"]["aggregate"]["levels"]
    for level in levels:
        level["radial_force"] = [-float(v) for v in level["radial_force"]]
    report = mod.audit_serialized_cases(cases)
    assert report["audit_pass"] is False
    with pytest.raises(RuntimeError, match="radial-force audit failed"):
        mod.enforce_scientific_gates(report)


def test_shifted_stress_mutation_breaks_public_piece_attribution() -> None:
    cases = _cases()
    target = cases[0]["stress_receipts"][0]["axial_e1"]["aggregate"]["stress"]
    target[len(target) // 2] += 0.1
    with pytest.raises(ValueError, match="piece attribution"):
        mod.audit_serialized_cases(cases)


def test_geometry_parameter_drift_fails_closed() -> None:
    cases = _cases()
    cases[1]["force_receipt"]["geometry"]["bump_center"] += 0.005
    with pytest.raises(ValueError, match="bump center drifted"):
        mod.audit_serialized_cases(cases)


def test_truth_boundary_promotion_fails_closed() -> None:
    cases = _cases()
    cases[2]["force_receipt"]["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        mod.audit_serialized_cases(cases)


def test_missing_shift_fails_closed() -> None:
    cases = _cases()
    cases[0]["stress_receipts"] = cases[0]["stress_receipts"][:-1]
    with pytest.raises(ValueError, match="exactly eight"):
        mod.audit_serialized_cases(cases)


def test_public_audit_signature_has_no_scientific_tuning_knobs() -> None:
    boundary = mod.truth_boundary()
    assert boundary["implementation_distinct_from_parent_fd2"] is True
    assert boundary["independent_derivative_rule"] == "centered_fd4_five_point"
    assert boundary["caller_scientific_tuning_absent"] is True
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == pytest.approx(1e-3)
    assert boundary["final_normalized_divergence_gate"] == pytest.approx(1e-5)
