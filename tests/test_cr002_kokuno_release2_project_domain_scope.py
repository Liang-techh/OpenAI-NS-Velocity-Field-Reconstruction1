from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_release2_project_domain_scope import (
    CONTRACT_REL,
    audit_contract,
)
from openai_ns_reconstruction.kokuno_pa16_current_cartesian_postswirl_release2 import (
    KokunoPA16CurrentCartesianPostSwirlRelease2,
)


ROOT = Path(__file__).resolve().parents[1]


def _contract():
    return json.loads((ROOT / CONTRACT_REL).read_text())


def test_release2_contract_replays_and_canonical_probe_is_rejected():
    receipt = audit_contract(root=ROOT)
    assert receipt["canonical_probe"] == [0.5, 0.0, 0.0, 0.5]
    assert receipt["canonical_probe_rejected_by_release2_stage"] is True
    assert "post-swirl release2 requires" in receipt["runtime_error"]
    assert receipt["canonical_eq45_velocity_export_ready"] is True
    assert receipt["kokuno_velocity_export_ready"] is False
    assert receipt["pde_validated"] is False
    assert receipt["paper_exact"] is False
    assert receipt["openai_field_identified"] is False


def test_canonical_probe_is_genuinely_inside_registered_box_support_and_time():
    c = _contract()
    x, y, z, t = c["canonical_probe_witness"]["point"]
    assert -2 <= x <= 2 and -2 <= y <= 2 and -2 <= z <= 2
    assert x * x + y * y < 4
    assert abs(z) < 2
    assert 0.25 <= t <= 0.75

    field = KokunoPA16CurrentCartesianPostSwirlRelease2()
    with pytest.raises(ValueError, match="post-swirl release2 requires"):
        field.velocity(x, y, z, t)


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("transfer_firewall", "stage_callable_implies_project_domain_coverage", True),
        ("transfer_firewall", "stage_receipt_transferable_to_global_export", True),
        ("transfer_firewall", "stage_receipt_transferable_to_visual_correspondence", True),
        ("transfer_firewall", "stage_receipt_transferable_to_pde_validation", True),
        ("transfer_firewall", "stage_receipt_transferable_to_paper_exactness", True),
        ("transfer_firewall", "stage_receipt_transferable_to_openai_field_identity", True),
        ("independent_states", "kokuno_velocity_export_ready", True),
        ("independent_states", "kokuno_visual_correspondence_verified", True),
        ("independent_states", "kokuno_pde_validated", True),
        ("independent_states", "kokuno_paper_exact", True),
        ("independent_states", "kokuno_openai_field_identified", True),
    ],
)
def test_stage_to_global_or_scientific_promotion_mutations_fail_closed(section, key, value):
    c = copy.deepcopy(_contract())
    c[section][key] = value
    with pytest.raises(ValueError):
        audit_contract(c, root=ROOT, run_runtime_probe=False)


def test_canonical_threshold_mutation_cannot_be_hidden_in_scope_contract():
    c = copy.deepcopy(_contract())
    c["canonical_contract"]["momentum_max"] = 0.01
    with pytest.raises(ValueError, match="momentum max gate drift"):
        audit_contract(c, root=ROOT, run_runtime_probe=False)


def test_canonical_constraints_and_eq45_delivery_are_not_downgraded_by_stage_gap():
    c = copy.deepcopy(_contract())
    c["independent_states"]["canonical_eq45_velocity_export_ready"] = False
    with pytest.raises(ValueError, match="canonical export state drift"):
        audit_contract(c, root=ROOT, run_runtime_probe=False)


def test_four_provenance_classes_remain_exactly_separate():
    c = copy.deepcopy(_contract())
    assert set(c["provenance_classes"]) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    c["provenance_classes"]["source_exact_hidden_parameter"] = []
    with pytest.raises(ValueError, match="four-way provenance partition drift"):
        audit_contract(c, root=ROOT, run_runtime_probe=False)
