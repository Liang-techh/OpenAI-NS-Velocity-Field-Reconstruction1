import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_heat_frontier_checkpoint import (
    AGENT2_UPSTREAM,
    AGENT3_UPSTREAM,
    AGENT4_UPSTREAM,
    REFERENCE_METRICS,
    STATES,
    _sha,
    build_checkpoint,
    summarize_heat_audit,
    validate_checkpoint,
)
from openai_ns_reconstruction.kokuno_independent_heat_repair import (
    DIVERGENCE_REFERENCE,
    LEVELS,
    RESIDUAL_REFERENCE,
    SEED,
    TASK,
    UPSTREAM_HEAD,
)


def _heat_fixture():
    return {
        "task": TASK,
        "upstream_agent1_head": UPSTREAM_HEAD,
        "seed": SEED,
        "levels": list(LEVELS),
        "sample_count": 6,
        "summary": {**REFERENCE_METRICS, "structural_preflight_passed": True},
        "fixed_references": {
            "normalized_ns_residual": RESIDUAL_REFERENCE,
            "divergence": DIVERGENCE_REFERENCE,
            "changed": False,
        },
        "truth_boundary": {
            "actual_heat_discrepancy_supplied": False,
            "global_leading_profile_reconstructed": False,
            "complete_kokuno_composite_velocity": False,
            "formal_full_domain_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e-3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def _resign(payload):
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)
    return payload


def test_heat_audit_summary_binds_exact_agent4_reference():
    summary = summarize_heat_audit(_heat_fixture())
    assert summary["structural_preflight_passed"] is True
    assert summary["seed"] == SEED
    assert summary["levels"] == list(LEVELS)
    for key, expected in REFERENCE_METRICS.items():
        assert summary[key] == pytest.approx(expected)


def test_checkpoint_keeps_latest_frontier_fail_closed():
    checkpoint = build_checkpoint(_heat_fixture())
    assert checkpoint["upstream"]["agent2"] == AGENT2_UPSTREAM
    assert checkpoint["upstream"]["agent3"] == AGENT3_UPSTREAM
    assert checkpoint["upstream"]["agent4"] == AGENT4_UPSTREAM
    assert checkpoint["routing"]["instantiate_heat_repair_into_global_velocity"] is False
    assert checkpoint["routing"]["compose_agent2_leading_bridge_into_complete_curl"] is False
    assert checkpoint["routing"]["apply_agent3_signed_curl_correction"] is False
    assert checkpoint["routing"]["run_formal_full_domain_gate"] is False
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["independent_heat_repair_preflight_ready"] is True
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert validate_checkpoint(checkpoint)["checkpoint_sha256"] == checkpoint["checkpoint_sha256"]


def test_checkpoint_rejects_laundering_or_unavailable_composition(tmp_path):
    base = build_checkpoint(_heat_fixture())
    mutations = []

    relaxed = copy.deepcopy(base)
    relaxed["fixed_gates"]["held_out_normalized_full_momentum_residual"] = 1.0e-2
    mutations.append(_resign(relaxed))

    pde = copy.deepcopy(base)
    pde["states"]["pde_validated"] = True
    mutations.append(_resign(pde))

    heat = copy.deepcopy(base)
    heat["routing"]["instantiate_heat_repair_into_global_velocity"] = True
    mutations.append(_resign(heat))

    phase = copy.deepcopy(base)
    phase["routing"]["compose_agent2_leading_bridge_into_complete_curl"] = True
    mutations.append(_resign(phase))

    signed = copy.deepcopy(base)
    signed["routing"]["apply_agent3_signed_curl_correction"] = True
    mutations.append(_resign(signed))

    ancestry = copy.deepcopy(base)
    ancestry["upstream"]["agent2"]["consumed_in_executable_ancestry"] = True
    mutations.append(_resign(ancestry))

    for index, payload in enumerate(mutations):
        path = tmp_path / f"mutation_{index}.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(ValueError):
            validate_checkpoint(json.loads(path.read_text(encoding="utf-8")))
