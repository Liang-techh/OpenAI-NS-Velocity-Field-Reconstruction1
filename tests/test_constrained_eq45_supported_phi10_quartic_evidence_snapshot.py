from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_evidence_snapshot import (
    EXPECTED_CANDIDATE_SHA256,
    audit_snapshot,
    load_snapshot,
)


def test_quartic_evidence_snapshot_binds_delivery_and_completed_receipts() -> None:
    snapshot = load_snapshot()
    receipt = audit_snapshot(snapshot)

    assert receipt["candidate_sha256"] == EXPECTED_CANDIDATE_SHA256
    assert receipt["velocity_export_ready"] is True
    assert receipt["visual_baseline_for_next_public_comparison"] == (
        "quartic_derivative_balanced_phi10"
    )
    assert receipt["candidate_selection_resolved"] is False
    assert receipt["quartic_vs_cubic_pde_ordering"] == "unresolved_seed_sensitive"
    assert receipt["visualization_ready"] is False
    assert receipt["pde_validated"] is False
    assert receipt["acceptance_ready"] is False

    morphology = snapshot["morphology_evidence"]
    assert morphology["pr"] == 170
    assert morphology["grid_sizes"] == [49, 65, 81]
    assert morphology["finest_81"]["core_superlevel_support_collar_occupancy"] == [
        0.0,
        0.0,
        0.0,
    ]

    fresh = snapshot["fresh_seed_pde_evidence"]
    assert fresh["pr"] == 171
    assert fresh["quartic_vs_static_wins"] == 3
    assert fresh["quartic_vs_cubic_wins"] == 1
    assert fresh["ordering"] == "unresolved_seed_sensitive"

    siblings = snapshot["unconsumed_sibling_evidence"]
    assert siblings["c2_compact_family_prs"] == [175, 178, 180]
    assert siblings["c2_compact_candidate_selection_consumed"] is False


def test_quartic_evidence_snapshot_keeps_export_visual_and_pde_states_independent() -> None:
    status = load_snapshot()["status"]

    assert status["velocity_export_ready"] is True
    assert status["visualization_candidate_only"] is True
    assert status["visualization_ready"] is False
    assert status["visual_correspondence_verified"] is False
    assert status["pde_validated"] is False
    assert status["paper_exact"] is False
    assert status["openai_field_identified"] is False
    assert status["blowup_proved"] is False


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("candidate_sha256",), "0" * 64),
        (("morphology_evidence", "public_openai_reference_used"), True),
        (("fresh_seed_pde_evidence", "quartic_vs_cubic_wins"), 3),
        (("fresh_seed_pde_evidence", "ordering"), "quartic_superior"),
        (("routing", "candidate_selection_resolved"), True),
        (("routing", "quartic_pde_superiority_over_cubic"), True),
        (("routing", "acceptance_ready"), True),
        (("status", "visualization_ready"), True),
        (("status", "pde_validated"), True),
    ],
)
def test_quartic_evidence_snapshot_fails_closed_on_promotion_or_provenance_drift(
    path: tuple[str, ...], value: object
) -> None:
    payload = deepcopy(load_snapshot())
    target = payload
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value

    with pytest.raises(ValueError):
        audit_snapshot(payload)
