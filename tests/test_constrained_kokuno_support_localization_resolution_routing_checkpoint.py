from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_support_localized_axis_refinement import (
    build_report as build_agent4_report,
)
from openai_ns_reconstruction.kokuno_support_localization_resolution_routing_checkpoint import (
    FIXED_GATES,
    STATES,
    _sha,
    _validate_agent4_report,
    build_checkpoint,
    validate_checkpoint,
    write_bundle,
)


@pytest.fixture(scope="module")
def agent4_report() -> dict:
    return build_agent4_report()


def test_checkpoint_records_resolution_of_parent_local_rejection(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    assert validate_checkpoint(checkpoint) == checkpoint
    audit = checkpoint["agent4_axis_refinement"]
    finest = audit["finest"]
    assert audit["manufactured_structural_preflight_passed"] is True
    assert audit["source_wave_shell_preflight_assessed"] is False
    assert finest["step"] == pytest.approx(2.5e-4)
    assert finest["curl_relative_rms"] == pytest.approx(2.1468147348041597e-10, rel=1.0e-3)
    # Agent-4's independent FD4 receipt is reproducible at the scientific-guard
    # scale but can move a few 1e-4 relatively across current NumPy/SciPy builds.
    # Keep the actual 2e-4 scientific guard frozen; do not turn a historical
    # floating value into a stricter accidental CI threshold.
    assert finest["divergence_max_abs"] == pytest.approx(
        1.3652247786009714e-08, rel=1.0e-3
    )
    assert finest["regions"]["small_radius"]["divergence_max_abs"] == pytest.approx(
        1.3652247786009714e-08, rel=1.0e-3
    )
    assert min(audit["refinement"]["curl_relative_rms_coarse_to_fine_ratios"]) > 15.9
    assert min(audit["refinement"]["divergence_rms_coarse_to_fine_ratios"]) > 15.9
    assert audit["mutation"]["mutated_divergence_max_abs"] == pytest.approx(
        0.06275470589885247, rel=1.0e-3
    )


def test_checkpoint_routes_new_source_and_correction_blockers(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["manufactured_support_localization_axis_refinement_passed"] is True
    assert checkpoint["states"]["source_wave_shell_independent_cartesian_preflight_passed"] is False
    assert checkpoint["upstream"]["agent2"]["pr"] == 358
    assert checkpoint["upstream"]["agent2"]["source_wave_shell_wrapper_ready"] is True
    assert checkpoint["upstream"]["agent2"]["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert checkpoint["upstream"]["agent3"]["pr"] == 357
    assert checkpoint["upstream"]["agent3"]["radial_force_derivative_stability_preflight_passed"] is True
    assert checkpoint["upstream"]["agent3"]["radial_to_tangential_force_rms_ratio"] == pytest.approx(
        0.4638871624620116
    )
    assert checkpoint["upstream"]["agent3"]["finite_correction_cycle_rerun_allowed"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False


def test_checkpoint_rejects_false_source_domain_promotion(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["agent4_axis_refinement"]["source_wave_shell_preflight_assessed"] = True
    unsigned = dict(mutated)
    unsigned.pop("checkpoint_sha256")
    mutated["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="source-wave-shell audit was falsely promoted"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_dropped_radial_force(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["truth_boundary"]["radial_force_side_effect_dropped"] = True
    unsigned = dict(mutated)
    unsigned.pop("checkpoint_sha256")
    mutated["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="truth boundary promoted"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_threshold_change(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["fixed_gates"]["held_out_normalized_full_momentum_max"] = 2.0e-3
    unsigned = dict(mutated)
    unsigned.pop("checkpoint_sha256")
    mutated["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="fixed gates changed"):
        validate_checkpoint(mutated)


def test_agent4_same_guard_pass_is_required(agent4_report: dict) -> None:
    rejected = copy.deepcopy(agent4_report)
    rejected["local_guards"]["structural_preflight_passed"] = False
    with pytest.raises(ValueError, match="manufactured local pass is absent"):
        _validate_agent4_report(rejected)


def test_bundle_round_trip_retains_global_truth_boundary(tmp_path) -> None:
    paths = write_bundle(tmp_path / "bundle")
    checkpoint = json.loads(paths["checkpoint"].read_text(encoding="utf-8"))
    validate_checkpoint(checkpoint)
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["truth_boundary"]["pde_validated"] is False
