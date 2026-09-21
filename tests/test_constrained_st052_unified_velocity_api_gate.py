from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

from openai_ns_reconstruction.constrained_st052_unified_velocity_api_gate import (
    CONTRACT_PATH,
    audit_contract,
    audit_repository,
)


ROOT = Path(__file__).resolve().parents[1]


def _json(path: str | Path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _inputs():
    contract = _json(CONTRACT_PATH)
    pinned = contract["pinned_repository_evidence"]
    return {
        "contract": contract,
        "acceptance": _json(pinned["runtime_acceptance_contract"]),
        "project_status": _json(pinned["project_status"]),
        "readiness_contract": _json(pinned["delivery_readiness_contract"]),
        "velocity_delivery_contract": _json(pinned["velocity_delivery_contract"]),
        "capsule_source": (ROOT / pinned["whole_child_module"]).read_text(encoding="utf-8"),
    }


def _audit(data):
    return audit_contract(
        data["contract"],
        acceptance=data["acceptance"],
        project_status=data["project_status"],
        readiness_contract=data["readiness_contract"],
        velocity_delivery_contract=data["velocity_delivery_contract"],
        capsule_source=data["capsule_source"],
    )


def test_repository_live_replay_passes_exact_pinned_gate():
    receipt = audit_repository(ROOT)
    assert receipt["passed"] is True, receipt["violations"]
    assert receipt["violations"] == []
    assert receipt["truth_boundary"] == {
        "st052_velocity_export_ready": False,
        "canonical_eq45_velocity_export_ready": True,
        "pde_pending_blocks_callable_delivery": False,
        "visual_correspondence_pending_blocks_callable_delivery": False,
    }


def test_baseline_contract_has_no_semantic_violation():
    assert _audit(_inputs()) == []


def test_forbidden_state_and_scope_promotions_fail_closed():
    baseline = _inputs()

    mutations = []

    d = deepcopy(baseline)
    d["contract"]["current_st052_delivery_scope"]["repository_unified_st052_velocity_entrypoint_registered"] = True
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["current_st052_delivery_scope"]["velocity_export_ready"] = True
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["readiness_semantics"]["accepted_external_runtime_dependency_implies_velocity_export_ready"] = True
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["readiness_semantics"]["pde_failure_or_pending_blocks_velocity_delivery"] = True
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["promotion_gate"]["state"] = "ready"
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["canonical_eq45_independence"]["velocity_export_ready"] = False
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["cr001_nonmutation"]["momentum_max"] = 0.002
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["cr001_nonmutation"]["thresholds_changed"] = True
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["cr001_nonmutation"]["residual_defined_free_forcing_forbidden"] = False
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["source_classification"]["public_source_fact"] = ["repository runtime decision"]
    mutations.append(d)

    d = deepcopy(baseline)
    d["contract"]["mutation_scope"]["velocity_changed"] = True
    mutations.append(d)

    for mutated in mutations:
        assert _audit(mutated), "forbidden contract mutation unexpectedly passed"


def test_upstream_evidence_laundering_fails_closed():
    baseline = _inputs()

    d = deepcopy(baseline)
    d["acceptance"]["truth_boundary"]["standalone_package_parent_runtime_ready"] = True
    assert _audit(d)

    d = deepcopy(baseline)
    d["acceptance"]["truth_boundary"]["velocity_export_ready"] = True
    assert _audit(d)

    d = deepcopy(baseline)
    d["project_status"]["latest_st052m_visual_candidate_evidence"]["velocity_export_ready"] = True
    assert _audit(d)

    d = deepcopy(baseline)
    d["readiness_contract"]["readiness_semantics"]["unqualified_velocity_export_ready"]["external_source_checkout_reconstruction_alone_is_sufficient"] = True
    assert _audit(d)

    d = deepcopy(baseline)
    d["velocity_delivery_contract"]["claim_gates"]["pde_validation"]["blocking_for_velocity_delivery"] = True
    assert _audit(d)

    d = deepcopy(baseline)
    d["velocity_delivery_contract"]["primary_deliverable"]["api"] = "openai_ns_reconstruction.st052_linear_temporal_capsule:load_bundle_runtime"
    assert _audit(d)


def test_capsule_api_shape_mutations_require_explicit_contract_update():
    baseline = _inputs()

    d = deepcopy(baseline)
    needle = "    exact_source_root: str | Path,\n"
    assert needle in d["capsule_source"]
    d["capsule_source"] = d["capsule_source"].replace(
        needle,
        "    exact_source_root: str | Path = Path('.'),\n",
        1,
    )
    assert _audit(d)

    d = deepcopy(baseline)
    d["capsule_source"] += "\n\ndef velocity(x, y, z, t):\n    raise RuntimeError('contract-update-required')\n"
    assert _audit(d)
