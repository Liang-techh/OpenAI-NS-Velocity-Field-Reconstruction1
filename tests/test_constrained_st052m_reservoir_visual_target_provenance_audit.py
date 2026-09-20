from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import constrained_st052m_reservoir_visual_target_provenance_audit as audit


def _inputs():
    return audit.load_default_inputs()


def test_live_reservoir_visual_target_provenance_audit_passes():
    receipt = audit.audit()
    assert receipt["audited_pr"] == 775
    assert receipt["numeric_probe_provenance"] == "autonomous_design"
    assert receipt["public_numeric_target"] is None
    assert receipt["targeted_development_evidence_only"] is True
    assert receipt["independent_visual_validation"] is False
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["pde_validated"] is False
    assert receipt["openai_field_identified"] is False
    assert receipt["cr001_unchanged"] is True


def test_numeric_visible_tip_geometry_cannot_be_relabelled_public_source_fact():
    contract, source, constraints = _inputs()
    mutated = copy.deepcopy(contract)
    moved = [item for item in mutated["classification"]["autonomous_design"] if "abs(z)=1.55 and 1.75" in item]
    assert len(moved) == 1
    mutated["classification"]["autonomous_design"].remove(moved[0])
    mutated["classification"]["public_source_fact"].append(moved[0])
    with pytest.raises(AssertionError, match="visible z probes must remain autonomous design"):
        audit.audit(mutated, source, constraints)


def test_public_numeric_target_invention_fails_closed():
    contract, source, constraints = _inputs()
    mutated = copy.deepcopy(contract)
    mutated["frozen_numeric_screen"]["public_numeric_target"] = -0.01
    with pytest.raises(AssertionError, match="public numerical target must remain null"):
        audit.audit(mutated, source, constraints)


def test_targeted_development_screen_cannot_be_called_independent_visual_validation():
    contract, source, constraints = _inputs()
    mutated = copy.deepcopy(contract)
    mutated["frozen_numeric_screen"]["independent_visual_validation"] = True
    with pytest.raises(AssertionError, match="relabelled independent"):
        audit.audit(mutated, source, constraints)


def test_targeted_development_screen_cannot_promote_visual_correspondence():
    contract, source, constraints = _inputs()
    mutated = copy.deepcopy(contract)
    mutated["truth_boundary"]["visual_correspondence_verified"] = True
    with pytest.raises(AssertionError, match="improperly promoted visual_correspondence_verified"):
        audit.audit(mutated, source, constraints)


def test_targeted_development_screen_cannot_promote_pde_or_exact_identity():
    contract, source, constraints = _inputs()
    for key in ("pde_validated", "paper_exact", "openai_field_identified"):
        mutated = copy.deepcopy(contract)
        mutated["truth_boundary"][key] = True
        with pytest.raises(AssertionError, match=f"improperly promoted {key}"):
            audit.audit(mutated, source, constraints)


def test_agent7_probe_coordinate_drift_is_detected():
    contract, source, constraints = _inputs()
    mutated_source = source.replace("VISIBLE_ABS_Z = (1.55, 1.75)", "VISIBLE_ABS_Z = (1.50, 1.75)", 1)
    assert mutated_source != source
    with pytest.raises(AssertionError, match="#775 visible z probes drifted"):
        audit.audit(contract, mutated_source, constraints)


def test_agent7_public_numeric_target_disclaimer_is_required():
    contract, source, constraints = _inputs()
    original = "outer return-flow location is a project representation choice; no OpenAI/public numeric threshold exists"
    mutated_source = source.replace(original, "outer return-flow location passed a public OpenAI threshold", 1)
    assert mutated_source != source
    with pytest.raises(AssertionError, match="public-target disclaimer drifted"):
        audit.audit(contract, mutated_source, constraints)


def test_agent7_truth_promotion_is_detected():
    contract, source, constraints = _inputs()
    mutated_source = source.replace('"visual_correspondence_verified": False', '"visual_correspondence_verified": True', 1)
    assert mutated_source != source
    with pytest.raises(AssertionError, match="#775 truth boundary promoted visual_correspondence_verified"):
        audit.audit(contract, mutated_source, constraints)


def test_cr001_threshold_relaxation_is_detected():
    contract, source, constraints = _inputs()
    mutated = copy.deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(AssertionError, match="CR001 momentum gates drifted"):
        audit.audit(contract, source, mutated)


def test_cr001_free_force_escape_is_detected():
    contract, source, constraints = _inputs()
    mutated = copy.deepcopy(constraints)
    mutated["forcing"]["restriction"] = "Pointwise force may be fitted from the residual."
    with pytest.raises(AssertionError, match="free-force shortcut appeared"):
        audit.audit(contract, source, mutated)


def test_open_external_public_observable_pr_cannot_be_laundered_live():
    contract, source, constraints = _inputs()
    mutated = copy.deepcopy(contract)
    mutated["public_observable_authority"]["integration_state_at_freeze"] = "live_integrated_and_green"
    with pytest.raises(AssertionError, match="may not be laundered"):
        audit.audit(mutated, source, constraints)
