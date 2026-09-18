from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st006_visual_evidence_card import (
    CANDIDATE_SHA256,
    OBSERVABLE_IDS,
    SOURCE_BINDINGS,
    audit_st006_visual_evidence_card,
    make_st006_visual_evidence_card,
)


def test_card_binds_actual_st006_receipts_and_public_observables():
    card = make_st006_visual_evidence_card()
    audit_st006_visual_evidence_card(card)
    assert card["candidate_sha256"] == CANDIDATE_SHA256
    assert tuple(card["observables"]) == OBSERVABLE_IDS
    assert card["evidence_sources"] == SOURCE_BINDINGS
    assert card["observables"]["inward_spiraling"]["measurements"]["inward_path_count"] == 48
    assert card["observables"]["increasing_elongation"]["measurements"]["vorticity_aspect_ratio_final_over_initial"] == pytest.approx(1.1215380885092212)
    assert card["observables"]["central_region_shrinks"]["measurements"]["covariance_volume_proxy_final_over_initial"] == pytest.approx(0.774185498561495)
    assert card["observables"]["speed_increases"]["measurements"]["velocity_rms_final_over_initial"] == pytest.approx(1.0728912572743456)


def test_card_is_deterministic_and_does_not_define_a_visual_score():
    first = make_st006_visual_evidence_card()
    second = make_st006_visual_evidence_card()
    assert first == second
    assert first["card_sha256"] == second["card_sha256"]
    assert first["routing_summary"]["visual_score_defined"] is False
    assert first["routing_summary"]["candidate_selected_by_this_card"] is False
    assert first["truth_boundary"]["visualization_ready"] is False
    assert first["truth_boundary"]["pde_validated"] is False


@pytest.mark.parametrize(
    "mutation",
    [
        lambda c: c.__setitem__("candidate_sha256", "0" * 64),
        lambda c: c["public_reference"].__setitem__("hidden_frame_time_used", True),
        lambda c: c["public_reference"].__setitem__("quantitative_visual_target_used", True),
        lambda c: c["evidence_sources"]["material_paths"].__setitem__("workflow_run", 1),
        lambda c: c["observables"]["inward_spiraling"].__setitem__("target", 1.0),
        lambda c: c["observables"]["relative_angular_rotation_color_encoding"].__setitem__("status", "measured_candidate_side_only"),
        lambda c: c["routing_summary"].__setitem__("visual_score_defined", True),
        lambda c: c["routing_summary"].__setitem__("candidate_selected_by_this_card", True),
        lambda c: c["truth_boundary"].__setitem__("visual_correspondence_verified", True),
        lambda c: c["truth_boundary"].__setitem__("pde_validated", True),
    ],
)
def test_fail_closed_on_provenance_target_or_truth_laundering(mutation):
    card = make_st006_visual_evidence_card()
    mutation(card)
    with pytest.raises(ValueError):
        audit_st006_visual_evidence_card(card)


def test_fail_closed_on_measurement_or_digest_tampering():
    card = make_st006_visual_evidence_card()
    altered = deepcopy(card)
    altered["observables"]["inward_spiraling"]["measurements"]["mean_radius_change"] = 0.1
    with pytest.raises(ValueError):
        audit_st006_visual_evidence_card(altered)

    altered = deepcopy(card)
    altered["card_sha256"] = "f" * 64
    with pytest.raises(ValueError):
        audit_st006_visual_evidence_card(altered)
