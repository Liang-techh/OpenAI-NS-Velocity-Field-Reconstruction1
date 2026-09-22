from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_bridge_common_scalar_refutation import (
    KokunoCurrentBridgeCommonScalarRefutation,
    _strict_decrease_margin,
)


def test_strict_decrease_margin_matches_two_analytic_sign_cases() -> None:
    h = 0.005
    q = np.asarray([3.0, 3.0], dtype=float)
    p = np.asarray([1.0, -0.75], dtype=float)
    got = _strict_decrease_margin(q, p, h)
    expected = np.asarray(
        [
            (1.0 - h) * (q[0] - p[0]),
            (1.0 - h) * q[1] + h * p[1],
        ],
        dtype=float,
    )
    assert np.array_equal(got, expected)
    assert np.all(got > 0.0)


def test_truth_boundary_materializes_only_one_way_witness_mechanics() -> None:
    witness = KokunoCurrentBridgeCommonScalarRefutation()
    truth = witness.truth_boundary
    assert truth["current_l_minus_h_q_s_transport_materialized"] is True
    assert truth["current_common_scalar_bridge_witness_check_materialized"] is True
    assert truth["full_eta_interval_entry_above_qp_established"] is False
    assert truth["full_eta_interval_target_totality_established"] is False
    assert truth["full_eta_interval_root_uniqueness_established"] is False
    assert truth["full_eta_interval_transversality_established"] is False
    assert truth["smooth_eta_target_time_map_established"] is False
    assert truth["full_eta_interval_common_scalar_bridge_length_established"] is False
    assert truth["current_l_minus_h_matching_bridge_materialized"] is False
    assert truth["eta_dependent_cartesian_matching_boundary_materialized"] is False
    assert truth["current_cartesian_terminal_multiplier_composed"] is False
    assert truth["outer_global_leading_velocity_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_fixed_witness_grid_is_deterministic_sorted_and_covers_interval() -> None:
    witness = KokunoCurrentBridgeCommonScalarRefutation()
    eta1 = witness.witness_eta()
    eta2 = witness.witness_eta()
    assert np.array_equal(eta1, eta2)
    assert eta1.shape == (65,)
    assert np.all(np.diff(eta1) > 0.0)
    lo, hi = witness.eta_interval
    assert eta1[0] == pytest.approx(lo, rel=0.0, abs=4e-16)
    assert eta1[-1] == pytest.approx(hi, rel=0.0, abs=4e-16)


def test_unique_root_certificate_is_finite_and_one_way() -> None:
    witness = KokunoCurrentBridgeCommonScalarRefutation()
    eta = witness.witness_eta()
    cert = witness.unique_root_certificate(eta)
    for key in (
        "Q_s_entry",
        "P_entry",
        "entry_minus_Q_p",
        "strict_decrease_margin",
        "asymptotic_positive_coefficient",
        "numerical_margin_floor",
    ):
        arr = np.asarray(cert[key], dtype=float)
        assert arr.shape == eta.shape
        assert np.all(np.isfinite(arr))
    flags = np.asarray(cert["unique_positive_target_root_certified"], dtype=bool)
    assert flags.shape == eta.shape
    assert np.all(np.asarray(cert["numerical_margin_floor"], dtype=float) > 0.0)


def test_certified_actual_roots_have_sign_oriented_brackets_when_available() -> None:
    witness = KokunoCurrentBridgeCommonScalarRefutation()
    eta = witness.witness_eta()
    cert = witness.unique_root_certificate(eta)
    flags = np.asarray(cert["unique_positive_target_root_certified"], dtype=bool)
    selected = eta[flags][:4]
    if selected.size == 0:
        pytest.skip("bounded sufficient uniqueness certificate did not fire on witness grid")
    bracket = witness.unique_root_bracket(selected)
    lower = np.asarray(bracket["root_lower"], dtype=float)
    upper = np.asarray(bracket["root_upper"], dtype=float)
    assert np.all(np.isfinite(lower))
    assert np.all(np.isfinite(upper))
    assert np.all(lower >= 0.0)
    assert np.all(upper >= lower)
    scale = max(1.0, abs(witness.q_p))
    assert np.all(np.asarray(bracket["lower_residual"], dtype=float) >= -8e-15 * scale)
    assert np.all(np.asarray(bracket["upper_residual"], dtype=float) <= 8e-15 * scale)


def test_refutation_report_is_internally_consistent_without_truth_promotion() -> None:
    witness = KokunoCurrentBridgeCommonScalarRefutation()
    report = witness.refutation_report()
    assert report["witness_count"] == 65
    assert 0 <= report["unique_root_certified_count"] <= 65
    assert report["finite_witness_is_sufficient_only_for_refutation"] is True
    assert report["absence_of_refutation_does_not_establish_common_scalar_bridge"] is True
    assert report["full_eta_interval_entry_above_qp_established"] is False
    assert report["full_eta_interval_target_totality_established"] is False
    assert report["full_eta_interval_root_uniqueness_established"] is False
    assert report["full_eta_interval_transversality_established"] is False
    assert report["smooth_eta_target_time_map_established"] is False
    assert report["full_eta_interval_common_scalar_bridge_length_established"] is False
    assert report["current_l_minus_h_matching_bridge_materialized"] is False
    assert report["eta_dependent_cartesian_matching_boundary_materialized"] is False

    refuted = report["common_scalar_bridge_refuted_by_disjoint_unique_root_brackets"]
    if refuted:
        pair = report["disjoint_root_pair"]
        assert pair is not None
        e0, e1 = pair["earlier_root_bracket"]
        l0, l1 = pair["later_root_bracket"]
        assert e0 <= e1 < l0 <= l1
        assert pair["disjoint_gap_lower_bound"] > 0.0
        assert report["disjoint_gap_lower_bound"] == pair["disjoint_gap_lower_bound"]
    else:
        assert report["disjoint_root_pair"] is None
        assert report["disjoint_gap_lower_bound"] == 0.0


def test_save_load_replays_semantic_identity_and_report(tmp_path) -> None:
    witness = KokunoCurrentBridgeCommonScalarRefutation()
    path = witness.save_configuration(tmp_path / "scalar_refutation.json")
    loaded = KokunoCurrentBridgeCommonScalarRefutation.load_configuration(path)
    assert loaded.semantic_sha256 == witness.semantic_sha256
    assert loaded.refutation_report() == witness.refutation_report()


def test_configuration_truth_mutation_fails_closed(tmp_path) -> None:
    witness = KokunoCurrentBridgeCommonScalarRefutation()
    path = witness.save_configuration(tmp_path / "scalar_refutation.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["truth_boundary"]["full_eta_interval_common_scalar_bridge_length_established"] = True
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="configuration/provenance mismatch"):
        KokunoCurrentBridgeCommonScalarRefutation.load_configuration(path)


def test_public_api_exposes_no_scientific_or_bridge_selection_tuning_knobs() -> None:
    witness_sig = inspect.signature(KokunoCurrentBridgeCommonScalarRefutation.witness_eta)
    cert_sig = inspect.signature(KokunoCurrentBridgeCommonScalarRefutation.unique_root_certificate)
    bracket_sig = inspect.signature(KokunoCurrentBridgeCommonScalarRefutation.unique_root_bracket)
    report_sig = inspect.signature(KokunoCurrentBridgeCommonScalarRefutation.refutation_report)
    assert tuple(witness_sig.parameters) == ("self",)
    assert tuple(cert_sig.parameters) == ("self", "eta")
    assert tuple(bracket_sig.parameters) == ("self", "eta")
    assert tuple(report_sig.parameters) == ("self",)
    forbidden = {
        "residual",
        "forcing",
        "pressure",
        "optimizer",
        "tolerance",
        "threshold",
        "gain",
        "h",
        "precision",
        "iterations",
        "bridge_length",
        "eta_representative",
        "witness_count",
    }
    for sig in (witness_sig, cert_sig, bracket_sig, report_sig):
        assert forbidden.isdisjoint(sig.parameters)
    assert not hasattr(KokunoCurrentBridgeCommonScalarRefutation, "velocity")
