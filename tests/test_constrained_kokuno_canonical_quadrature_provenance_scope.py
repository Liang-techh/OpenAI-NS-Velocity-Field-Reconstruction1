from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from openai_ns_reconstruction import kokuno_a4_blackbox_full_ns_validator as parent
from openai_ns_reconstruction import kokuno_a4_canonical_quadrature_evidence as gate


ROOT = Path(__file__).resolve().parents[1]
SCOPE_PATH = ROOT / "configs/kokuno_canonical_quadrature_provenance_scope.json"
AUDITOR_PATH = ROOT / "scripts/audit_kokuno_canonical_quadrature_provenance_scope.py"


def _load_auditor():
    spec = importlib.util.spec_from_file_location("cr002_quadrature_provenance_auditor", AUDITOR_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fabricated_mechanics_levels():
    # Autonomous provenance witness only: these rows are deliberately authored
    # in this test and are NOT the output of any canonical quadrature run.
    levels = []
    for order, momentum, divergence in (
        (24, 6.0e-4, 6.0e-6),
        (48, 5.0e-4, 5.0e-6),
        (96, 4.9e-4, 4.9e-6),
    ):
        rows = tuple(
            gate.QuadratureTimeMetrics(
                time=float(t),
                momentum_volume_l2=momentum,
                divergence_volume_l2=divergence,
                kinetic_energy=1.0,
            )
            for t in parent.VALIDATION_TIMES
        )
        levels.append(gate.QuadratureOrderMetrics(order_per_axis=order, by_time=rows))
    return tuple(levels)


def test_auditor_fail_closed_scope_matches_exact_903_source_and_cr001_contract():
    report = _load_auditor().audit()
    assert report["audit_passed"] is True
    assert report["audited_head"] == "065810d94a37ade8d81acd2c7f622c1ff8cc6bdd"
    assert report["audited_source_blob"] == "4a14a2e2bd6c06280d992802414109444d41e84b"
    assert report["canonical_constraints_blob"] == "6c559e42895a606e2ef025ade4cb448966d75814"
    assert report["external_execution_provenance_input_present"] is False
    assert report["builder_self_asserts_convergence"] is True
    assert report["builder_self_asserts_independence"] is True
    assert report["builder_blocks_residual_defined_forcing"] is True
    assert report["receipt_integrity_is_execution_provenance"] is False
    assert report["checksum_valid_self_built_receipt_proves_independent_run"] is False
    assert report["pde_validated"] is False


def test_checksum_valid_self_built_receipt_is_only_an_integrity_witness():
    # The current #903 builder can produce a checksum-valid, gate-clean receipt
    # from rows supplied directly by this unit test.  That is useful receipt
    # integrity behavior, but it cannot prove where the rows came from.
    evidence = gate.build_canonical_quadrature_evidence(
        candidate_sha256="autonomous-mechanics-candidate",
        stage="leading_only",
        physical_contract_sha256="autonomous-mechanics-contract",
        levels=_fabricated_mechanics_levels(),
    )
    assert evidence.receipt_sha256 == evidence.recomputed_sha256()
    assert evidence.convergence_assessed is True
    assert evidence.independent_of_training_and_held_in_selection is True
    assert gate.canonical_quadrature_failures(
        evidence,
        candidate_sha256="autonomous-mechanics-candidate",
        stage="leading_only",
        physical_contract_sha256="autonomous-mechanics-contract",
    ) == []

    scope = json.loads(SCOPE_PATH.read_text(encoding="utf-8"))
    assert scope["scope"]["answer"] is False
    forbidden = scope["required_claim_separation"][
        "checksum_valid_receipt_may_not_claim_without_separate_execution_provenance"
    ]
    assert any("actually produced" in item for item in forbidden)
    assert any("independent of training" in item for item in forbidden)
    assert any("pde_validated" in item for item in forbidden)


def test_four_way_provenance_keeps_mechanics_witness_out_of_public_source_facts():
    scope = json.loads(SCOPE_PATH.read_text(encoding="utf-8"))
    provenance = scope["four_way_provenance"]
    assert set(provenance) == {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    assert provenance["public_source_fact"] == []
    assert any("mechanics/provenance witness" in item for item in provenance["autonomous_design"])
    assert any("No real canonical 24/48/96" in item for item in provenance["pending_unknown"])


def test_cr001_gates_and_anti_shortcuts_are_not_relaxed_by_this_audit():
    scope = json.loads(SCOPE_PATH.read_text(encoding="utf-8"))
    invariants = scope["canonical_cr001_invariants"]
    assert invariants["nu"] == 0.01
    assert invariants["time_interval"] == [0.25, 0.75]
    assert invariants["forcing_mode"] == "restricted_two_parameter_family"
    assert invariants["reference_energy"] == 1.0
    assert invariants["reference_energy_abs_tolerance"] == 0.001
    assert invariants["validation_seed"] == 914027
    assert invariants["held_out_points"] == 4096
    assert invariants["derivative_steps"] == [0.02, 0.01, 0.005]
    assert invariants["quadrature_orders_per_axis"] == [24, 48, 96]
    assert invariants["divergence_max"] == 1e-5
    assert invariants["divergence_L2"] == 1e-5
    assert invariants["pde_residual_max"] == 1e-3
    assert invariants["pde_residual_L2"] == 1e-3
    assert invariants["residual_defined_free_forcing_allowed"] is False
    assert invariants["amplitude_collapse_allowed"] is False
    assert invariants["threshold_relaxation_allowed"] is False


def test_delivery_visual_pde_and_exactness_states_remain_independent():
    scope = json.loads(SCOPE_PATH.read_text(encoding="utf-8"))
    states = scope["independent_truth_states"]
    assert states["canonical_eq45_velocity_export_ready"] is True
    assert states["kokuno_real_canonical_quadrature_execution_verified"] is False
    assert states["kokuno_pde_validated"] is False
    assert states["visual_correspondence_verified"] is False
    assert states["paper_exact"] is False
    assert states["openai_field_identified"] is False
