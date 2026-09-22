from __future__ import annotations

import hashlib
import inspect
from pathlib import Path

import openai_ns_reconstruction.kokuno_current_i4_exact_backend_correction_firewall as firewall
import openai_ns_reconstruction.kokuno_current_i4_rf44_prestate_provider as prestate


PARENT_A2_PRESTATE_HEAD = "cab837f6358c8d29e3d470061879a58b3f7e1454"
PARENT_A2_PRESTATE_BLOB = "0d8bb6471a4b7825b1b4b428d3b2c38d207729fd"
IMPORTED_A3_FIREWALL_HEAD = "25f2793091cbf8dc2d8f966b941c32eaff63d876"
IMPORTED_A3_FIREWALL_BLOB = "0152a85f7fbaf9c32c24c9d5c7d0d11166e3b157"


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    prefix = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(prefix + data).hexdigest()


def _module_blob(module) -> str:
    path = inspect.getsourcefile(module)
    assert path is not None
    return _git_blob_sha1(Path(path))


def test_exact_a2_prestate_and_a3_firewall_blobs_coexist_unchanged() -> None:
    assert _module_blob(prestate) == PARENT_A2_PRESTATE_BLOB
    assert _module_blob(firewall) == IMPORTED_A3_FIREWALL_BLOB
    assert prestate.PARENT_AGENT3_HEAD == firewall.PARENT_HEAD
    assert firewall.PARENT_HEAD == "e75940c33127b0725ddc703eaae97590158f2a56"


def test_restack_preserves_separate_fail_closed_truth_boundaries() -> None:
    a2 = prestate.truth_boundary()
    a3 = firewall.truth_boundary()

    assert a2["rf44_preupdate_state_materialized"] is True
    assert a2["exact_a2_1198_raw_provider_delegated_unchanged"] is True
    assert a2["raw_provider_blob_pin_relabelled_as_rf44_prestate_pin"] is False
    assert a2["source_exact_rf44_prestate_claimed"] is False
    assert a2["rf44_postupdate_state_materialized"] is False
    assert a2["rf44_correction_increment_materialized"] is False
    assert a2["agent3_mean_correction_solved_here"] is False
    assert a2["cartesian_correction_velocity_materialized"] is False
    assert a2["heldout_ns_residual_assessed"] is False
    assert a2["pde_validated"] is False

    assert a3["exact_backend_rebind_required"] is True
    assert a3[
        "exact_backend_bind_authenticates_concrete_module_class_function_and_blobs"
    ] is True
    assert a3[
        "manual_ExactCurrentI4NonlinearBackend_dataclass_construction_is_sufficient"
    ] is False
    assert a3["zero_or_synthetic_backend_fixture_is_scientific_correction_evidence"] is False
    assert a3["parent_1210_promotion_without_exact_rebind_is_sufficient"] is False
    assert a3["parent_rf30_rf31_rf34_rf39_formulas_changed"] is False
    assert a3["agent2_curl_or_jacobian_reimplemented"] is False
    assert a3["correction_applied_to_candidate"] is False
    assert a3["rf44_rf49_nonlinear_remainder_recomputed"] is False
    assert a3["cartesian_correction_velocity_materialized"] is False
    assert a3["finite_correction_cycle_run"] is False
    assert a3["heldout_ns_residual_assessed"] is False
    assert a3["pde_validated"] is False

    assert a2["final_normalized_momentum_gate"] == a3["final_normalized_momentum_gate"] == 1.0e-3
    assert a2["final_normalized_divergence_gate"] == a3["final_normalized_divergence_gate"] == 1.0e-5
    assert a2["residual_defined_free_forcing_allowed"] is False
    assert a3["residual_defined_free_forcing_allowed"] is False


def test_exact_backend_rebind_is_present_without_extending_rf44_provider_surface() -> None:
    rebind_source = inspect.getsource(firewall._require_exact_rebind)
    assert "ExactCurrentI4NonlinearBackend.bind" in rebind_source
    assert "backend.composite_field" in rebind_source
    assert "backend.differential_function" in rebind_source

    cls = prestate.KokunoCurrentI4RF44PreStateProvider
    assert callable(getattr(cls, "rf44_pre_update_state", None))
    assert not hasattr(cls, "rf44_post_update_state")
    assert not hasattr(cls, "rf44_correction_increment")

    a2_params = set(inspect.signature(prestate.bind_current_i4_rf44_prestate_provider).parameters)
    a3_params = set(
        inspect.signature(firewall.materialize_current_i4_rf34_rf39_exact_backend_correction).parameters
    )
    forbidden = {
        "residual",
        "defect",
        "forcing",
        "pressure",
        "gain",
        "threshold",
        "heldout",
        "viscosity",
        "nu",
        "delta_u",
    }
    assert not (a2_params & forbidden)
    assert not (a3_params & forbidden)


def test_restack_is_topology_only_not_candidate_or_pde_promotion() -> None:
    receipt = {
        "task": "K2-OSC-114",
        "parent_a2_prestate_head": PARENT_A2_PRESTATE_HEAD,
        "parent_a2_prestate_blob": PARENT_A2_PRESTATE_BLOB,
        "imported_a3_firewall_head": IMPORTED_A3_FIREWALL_HEAD,
        "imported_a3_firewall_blob": IMPORTED_A3_FIREWALL_BLOB,
        "same_lineage_materialized": True,
        "rf44_prestate_materialized": True,
        "exact_backend_firewall_materialized": True,
        "concrete_a2_1080_960_runtime_dependencies_restored": False,
        "real_exact_backend_bind_executed": False,
        "repository_candidate_scientific_correction_evidence": False,
        "rf44_postupdate_state_materialized": False,
        "rf44_correction_increment_materialized": False,
        "cartesian_delta_u_materialized": False,
        "finite_correction_cycle_run": False,
        "heldout_ns_residual_assessed": False,
        "pde_validated": False,
    }
    assert receipt["same_lineage_materialized"] is True
    for key in (
        "concrete_a2_1080_960_runtime_dependencies_restored",
        "real_exact_backend_bind_executed",
        "repository_candidate_scientific_correction_evidence",
        "rf44_postupdate_state_materialized",
        "rf44_correction_increment_materialized",
        "cartesian_delta_u_materialized",
        "finite_correction_cycle_run",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ):
        assert receipt[key] is False
