from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import st052_stable_temporal_no_growth_decision as decision


def _receipt() -> dict:
    return {
        "task_id": decision.UPSTREAM_TASK,
        "issue": 1220,
        "measurement_sha256": decision.UPSTREAM_MEASUREMENT_SHA256,
        "stable_candidate_identity": {
            "candidate_semantic_identity_sha256": decision.CANDIDATE_SEMANTIC_SHA256,
            "velocity_semantic_identity_sha256": decision.VELOCITY_SEMANTIC_SHA256,
        },
        "public_observable": {
            "id": "shrinking_central_region_while_speed_increases",
            "numerical_target": None,
        },
        "measurement": {
            "derived": {
                "speed_weighted_rms_radius_by_time": list(decision.EXPECTED_WIDTHS),
                "peak_ring_mean_speed_by_time": list(decision.EXPECTED_SPEEDS),
                "endpoint_width_fraction": decision.EXPECTED_WIDTH_FRACTION,
                "endpoint_peak_speed_fraction": decision.EXPECTED_SPEED_FRACTION,
                "endpoint_shrink_and_speedup_proxy": True,
                "three_time_monotone_proxy": True,
                "three_time_width_monotone_decrease": True,
                "three_time_peak_speed_monotone_increase": True,
            },
            "routing": {
                "sixth_basis_authorized": False,
                "candidate_mutation_authorized": False,
            },
        },
        "truth_boundary": {
            "candidate_changed": False,
            "basis_dimension_changed": False,
            "coefficient_selected": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "source_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "visualization_fingerprint_direct_improvement": 0.0,
        },
    }


def test_exact_admitted_receipt_closes_only_coarse_temporal_absence_trigger() -> None:
    report = decision.decide(_receipt())
    route = report["decision"]
    assert route["coarse_missing_shrink_speedup_trigger_closed"] is True
    assert route["new_temporal_basis_authorized_from_this_observable"] is False
    assert route["sixth_basis_authorized"] is False
    assert route["coefficient_selection_authorized"] is False
    assert route["candidate_mutation_authorized"] is False
    assert route["existing_offgrid_temporal_capacity_pr"] == 1178
    assert route["offgrid_capacity_promotion_authorized"] is False
    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["basis_dimension_changed"] is False
    assert truth["visualization_fingerprint_direct_improvement"] == 0.0
    assert truth["pde_validated"] is False


@pytest.mark.parametrize(
    ("mutator", "match"),
    [
        (lambda r: r.__setitem__("measurement_sha256", "0" * 64), "measurement identity"),
        (
            lambda r: r["stable_candidate_identity"].__setitem__(
                "velocity_semantic_identity_sha256", "1" * 64
            ),
            "velocity semantic identity",
        ),
        (
            lambda r: r["measurement"]["derived"]["speed_weighted_rms_radius_by_time"].__setitem__(1, 0.99),
            "core widths numerical drift",
        ),
        (
            lambda r: r["public_observable"].__setitem__("numerical_target", 0.12),
            "invented public numerical target",
        ),
        (
            lambda r: r["truth_boundary"].__setitem__("pde_validated", True),
            "truth promotion",
        ),
    ],
)
def test_evidence_drift_fails_closed(mutator, match: str) -> None:
    receipt = copy.deepcopy(_receipt())
    mutator(receipt)
    with pytest.raises(ValueError, match=match):
        decision.decide(receipt)


def test_nonmonotone_temporal_evidence_cannot_be_laundered_into_no_growth() -> None:
    receipt = copy.deepcopy(_receipt())
    receipt["measurement"]["derived"]["speed_weighted_rms_radius_by_time"] = [
        decision.EXPECTED_WIDTHS[0],
        decision.EXPECTED_WIDTHS[0] + 0.01,
        decision.EXPECTED_WIDTHS[2],
    ]
    with pytest.raises(ValueError):
        decision.decide(receipt)
