from __future__ import annotations

import json

from openai_ns_reconstruction.kokuno_agent4_autonomous_signed_rectangle_independent_audit import (
    CASES,
    FROZEN_GUARDS,
    run_audit,
)


def test_independent_autonomous_signed_rectangle_audit_passes_frozen_guards() -> None:
    report = run_audit()
    metrics = report["metrics"]

    assert report["local_structural_preflight_passed"] is True
    assert report["formal_full_domain_pde_gate_assessed"] is False
    assert report["heldout_ns_residual_assessed"] is False
    assert report["pde_validated"] is False

    assert metrics["schedule_relative_max"] <= FROZEN_GUARDS["schedule_relative_max"]
    assert metrics["center_float_relative_max"] <= FROZEN_GUARDS["center_float_relative_max"]
    assert metrics["permutation_relative_max"] <= FROZEN_GUARDS["permutation_relative_max"]
    assert metrics["covering_level_mismatches"] == 0
    assert metrics["center_index_mismatches"] == 0
    assert metrics["center_exact_mismatches"] == 0
    assert metrics["delta_max_mismatches"] == 0
    assert metrics["d_c_squared_exact_mismatches"] == 0
    assert metrics["strict_guard_violations"] == 0
    assert metrics["permutation_exact_mismatches"] == 0
    assert metrics["forged_partition_rejections"] == len(CASES)

    assert (
        metrics["order_dependent_mutation_mismatch_fraction"]
        >= FROZEN_GUARDS["order_dependent_mutation_mismatch_fraction_min"]
    )
    assert metrics["swapped_sign_mutation_mismatch_fraction"] == 1.0
    assert (
        metrics["wrong_pulse_formula_relative_rms"]
        >= FROZEN_GUARDS["wrong_pulse_formula_relative_rms_min"]
    )


def test_independent_audit_receipt_is_json_serializable_and_truthful() -> None:
    report = run_audit()
    payload = json.loads(json.dumps(report, allow_nan=False))

    assert payload["oracle"]["public_values_only"] is True
    assert payload["oracle"]["agent2_rational_witness_helper_reused"] is False
    assert payload["oracle"]["agent2_torus_distance_helper_reused"] is False
    assert payload["oracle"]["agent2_band_covering_helper_reused"] is False
    assert payload["oracle"]["agent2_center_assignment_helper_reused"] is False
    assert payload["physical_coordinate_probe_applicable"] is False
    assert payload["truth_boundary"]["autonomous_geometry_only"] is True
    assert payload["truth_boundary"]["actual_source_rectangle_centers_recovered"] is False
    assert payload["truth_boundary"]["actual_source_partition_labels_bound"] is False
    assert payload["truth_boundary"]["actual_source_h_sigma_bound"] is False
    assert payload["truth_boundary"]["public_source_bound_velocity_osc_materialized"] is False
    assert payload["truth_boundary"]["full_ns_momentum_gate_value"] is None
