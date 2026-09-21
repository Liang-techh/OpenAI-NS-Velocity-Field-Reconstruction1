from __future__ import annotations

import copy
import json
from pathlib import Path

from openai_ns_reconstruction.audit_kokuno_post_xr_leading_independent_audit_scope import audit_scope


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/kokuno_post_xr_leading_independent_audit_scope.json"


def _payload() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def _must_fail(mutator) -> None:
    payload = copy.deepcopy(_payload())
    mutator(payload)
    errors = audit_scope(ROOT, payload)
    assert errors, "governance mutation must fail closed"


def test_baseline_contract_passes() -> None:
    assert audit_scope(ROOT) == []


def test_cannot_promote_post_xr_independent_cartesian_audit() -> None:
    _must_fail(
        lambda p: p["machine_locked_distinctions"].__setitem__(
            "post_xr_leading_independent_cartesian_public_velocity_audit_available", True
        )
    )


def test_cannot_relabel_through_xr_a4_as_post_xr_coverage() -> None:
    _must_fail(
        lambda p: p["machine_locked_distinctions"].__setitem__(
            "through_xr_a4_audit_may_be_relabelled_as_post_xr_audit", True
        )
    )


def test_cannot_relabel_source_coordinate_self_check_as_independent_cartesian_validation() -> None:
    _must_fail(
        lambda p: p["machine_locked_distinctions"].__setitem__(
            "source_coordinate_handoff_or_derivative_replay_may_be_relabelled_as_independent_cartesian_post_xr_validation",
            True,
        )
    )


def test_post_xr_leading_does_not_imply_post_xr_composite() -> None:
    _must_fail(
        lambda p: p["machine_locked_distinctions"].__setitem__(
            "post_xr_leading_materialization_implies_post_xr_composite_materialization", True
        )
    )


def test_post_xr_leading_does_not_imply_global_totality_or_pde() -> None:
    def mutate(payload: dict) -> None:
        locked = payload["machine_locked_distinctions"]
        locked["post_xr_leading_materialization_implies_global_project_domain_totality"] = True
        locked["post_xr_leading_materialization_implies_pde_validation"] = True

    _must_fail(mutate)


def test_exact_a1_986_identity_is_bound() -> None:
    _must_fail(
        lambda p: p["upstream_scope"]["agent1_986"].__setitem__(
            "head", "0" * 40
        )
    )


def test_exact_a4_989_identity_is_bound() -> None:
    _must_fail(
        lambda p: p["upstream_scope"]["agent4_989"].__setitem__(
            "head", "f" * 40
        )
    )


def test_exact_a4_983_identity_is_bound() -> None:
    _must_fail(
        lambda p: p["upstream_scope"]["agent4_983"].__setitem__(
            "head", "e" * 40
        )
    )


def test_momentum_gate_cannot_be_relaxed() -> None:
    _must_fail(lambda p: p["cr001_snapshot"].__setitem__("momentum_max", 1.0e-2))


def test_divergence_gate_cannot_be_relaxed() -> None:
    _must_fail(lambda p: p["cr001_snapshot"].__setitem__("divergence_l2", 1.0e-4))


def test_residual_defined_free_forcing_cannot_be_enabled() -> None:
    _must_fail(
        lambda p: p["cr001_snapshot"].__setitem__(
            "residual_defined_free_forcing_forbidden", False
        )
    )


def test_candidate_collapse_cannot_be_enabled() -> None:
    _must_fail(
        lambda p: p["cr001_snapshot"].__setitem__("candidate_collapse_forbidden", False)
    )


def test_post_hoc_threshold_relaxation_cannot_be_enabled() -> None:
    _must_fail(
        lambda p: p["cr001_snapshot"].__setitem__(
            "post_hoc_threshold_relaxation_forbidden", False
        )
    )


def test_canonical_eq45_delivery_cannot_be_downgraded_by_kokuno_gap() -> None:
    _must_fail(
        lambda p: p["canonical_delivery_independence"].__setitem__(
            "velocity_export_ready", False
        )
    )


def test_canonical_eq45_pde_state_cannot_be_promoted() -> None:
    _must_fail(
        lambda p: p["canonical_delivery_independence"].__setitem__(
            "pde_validated", True
        )
    )


def test_four_way_provenance_bucket_names_are_fixed() -> None:
    def mutate(payload: dict) -> None:
        payload["four_way_provenance"]["paper_fact"] = payload["four_way_provenance"].pop(
            "public_source_facts"
        )

    _must_fail(mutate)


def test_public_source_bucket_cannot_claim_exact_openai_field() -> None:
    def mutate(payload: dict) -> None:
        payload["four_way_provenance"]["public_source_facts"].append(
            "exact OpenAI field"
        )

    _must_fail(mutate)
