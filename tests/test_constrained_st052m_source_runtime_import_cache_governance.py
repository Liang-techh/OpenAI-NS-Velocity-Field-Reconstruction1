from __future__ import annotations

import pytest

from openai_ns_reconstruction.constrained_st052m_source_runtime_import_cache_governance import (
    audit_contract,
    load_contract,
    mutated,
    observe_loader,
)


def test_current_import_cache_contract_passes() -> None:
    report = audit_contract(load_contract())
    assert report["passed"], report["errors"]
    assert report["truth_boundary"]["clean_exact_ci_receipt_retained"] is True
    assert report["truth_boundary"]["runtime_identity_closed"] is False
    assert report["truth_boundary"]["velocity_export_ready"] is False


def test_live_loader_observation_records_unisolated_bare_import() -> None:
    observed = observe_loader()
    assert observed == {
        "uses_sys_path_prepend": True,
        "uses_bare_replay_import": True,
        "references_sys_modules": False,
        "checks_imported_module_file_origin": False,
    }


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("source_runtime", "preexisting_direct_module_cache_rejected"), True),
        (("source_runtime", "post_import_direct_module_file_origin_checked"), True),
        (("source_runtime", "transitive_module_cache_isolation_verified"), True),
        (("allowed_states", "direct_replay_import_cache_isolation_verified"), True),
        (("allowed_states", "transitive_source_runtime_import_cache_isolation_verified"), True),
        (("allowed_states", "runtime_identity_closed"), True),
        (("allowed_states", "standalone_reproducible_velocity_identity_ready"), True),
        (("allowed_states", "velocity_export_ready"), True),
        (("allowed_states", "visual_correspondence_verified"), True),
        (("allowed_states", "pde_validated"), True),
        (("allowed_states", "paper_exact"), True),
        (("allowed_states", "openai_field_identified"), True),
    ],
)
def test_premature_import_runtime_or_scientific_promotions_fail_closed(
    path: tuple[str, str], value: bool
) -> None:
    report = audit_contract(mutated(load_contract(), *path, value=value))
    assert report["passed"] is False
    assert report["errors"]


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("integration_base", "head"), "0" * 40),
        (("source_runtime", "source_head"), "0" * 40),
        (("evidence_scope", "exact_source_parity_gate"), 1e-9),
        (("cr001_nonmutation", "validation_seed"), 1),
        (("cr001_nonmutation", "held_out_points"), 1024),
        (("cr001_nonmutation", "pde_residual_max"), 0.01),
        (("cr001_nonmutation", "divergence_max"), 0.001),
        (("cr001_nonmutation", "reference_energy_abs_tolerance"), 0.01),
    ],
)
def test_identity_receipt_and_cr001_drift_fail_closed(
    path: tuple[str, str], value: object
) -> None:
    report = audit_contract(mutated(load_contract(), *path, value=value))
    assert report["passed"] is False
    assert report["errors"]


def test_runtime_bridge_remains_autonomous_and_cache_closure_pending() -> None:
    contract = load_contract()
    autonomous = contract["classification"]["autonomous_design"]
    pending = contract["classification"]["pending_unknown"]
    assert any("sys.path" in item for item in autonomous)
    assert any("sys.modules" in item for item in pending)

    wrong = mutated(
        contract,
        "classification",
        "autonomous_design",
        value=["ST052-M transform only"],
    )
    assert audit_contract(wrong)["passed"] is False
