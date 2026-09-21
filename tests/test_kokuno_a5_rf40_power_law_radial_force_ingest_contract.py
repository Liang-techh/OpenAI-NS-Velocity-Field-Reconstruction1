from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import kokuno_a5_rf40_power_law_radial_force_ingest_contract as a5


def _mutated(registration, mutate):
    out = copy.deepcopy(registration)
    mutate(out)
    payload = {k: copy.deepcopy(v) for k, v in out.items() if k not in {"schema", "task_id", "digest"}}
    out["digest"] = a5._sha256(payload)
    return out


def test_exact_x4_radial_force_and_a4_audit_are_registered_without_pde_promotion():
    reg = a5.build_registration()
    a5.validate_registration(reg)

    assert reg["parent_a5"]["pr"] == 1039
    assert reg["parent_a5"]["head"] == "6d88e8a83b51704f65aa125006b7671546c77c43"

    force = reg["agent3_power_law_force"]
    assert force["pr"] == 1045
    assert force["head"] == "3ec20a77b551be819a71e3308e9ecab2a2226406"
    assert force["source_blob"] == "39c98d3f4df05846848b897c240b3b56fb6fbe2c"
    assert force["parent_stress_head"] == "9a5cdcdf78b7862d5ebe171efb9ca56d53664612"
    assert force["agent2_composite_head"] == "e36d9da4b7f037e998f5b1658f8c0ea291a76b80"
    assert force["agent1_leading_head"] == "2c76ebdc41d6c566f43a1305034ba2ff9dce410b"
    assert force["formula"] == "(div T)_r = partial_z sigma_1"
    assert force["production_derivative"] == "centered_fd2"
    assert force["z_step_ladder"] == [0.02, 0.01, 0.005]

    audit = reg["agent4_power_law_force_audit"]
    assert audit["pr"] == 1046
    assert audit["head"] == "df5af9c503d69d665cea93d66fa30bf6f00863a8"
    assert audit["audited_agent3_head"] == force["head"]
    assert audit["audited_agent3_source_blob"] == force["source_blob"]
    assert audit["production_operator"] == "centered_fd2"
    assert audit["independent_operator"] == "five_point_centered_fd4"
    assert audit["implementation_distinct"] is True
    assert audit["frozen_seed"] == 9173791
    assert audit["heldout_z_centers"] == [-0.143, -0.067, 0.049]
    assert audit["radial_count"] == 169

    truth = reg["truth_boundary"]
    assert truth["current_nonlinear_radial_force_through_power_law_materialized"] is True
    assert truth["agent4_dedicated_power_law_radial_force_audit_present"] is True
    assert truth["agent4_independent_power_law_radial_force_audit_registered"] is True
    assert truth["agent4_independent_power_law_radial_force_audit_admitted"] is False
    assert truth["scoped_power_law_radial_force_authorized_as_ns_correction_target"] is False
    assert truth["newer_modulated_a1_lineage_consumed"] is False
    assert truth["evidence_transferred_from_newer_modulated_lineage"] is False
    assert truth["matched_cartesian_pressure_materialized"] is False
    assert truth["restricted_forcing_materialized"] is False
    assert truth["complete_ns_defect_materialized"] is False
    assert truth["real_agent3_ns_correction_velocity_materialized"] is False
    assert truth["real_candidate_finite_correction_cycle_run"] is False
    assert truth["complete_candidate_api_ready"] is False
    assert truth["pde_validated"] is False

    assert reg["readiness"] == a5.parent.READINESS
    assert reg["final_normalized_momentum_gate"] == 1.0e-3
    assert reg["final_normalized_divergence_gate"] == 1.0e-5
    assert reg["canonical_quadrature"] == [24, 48, 96]


@pytest.mark.parametrize(
    "mutate",
    [
        lambda r: r["truth_boundary"].__setitem__("complete_ns_defect_materialized", True),
        lambda r: r["truth_boundary"].__setitem__("pde_validated", True),
        lambda r: r["truth_boundary"].__setitem__("newer_modulated_a1_lineage_consumed", True),
        lambda r: r["truth_boundary"].__setitem__("agent4_independent_power_law_radial_force_audit_admitted", True),
        lambda r: r.__setitem__("final_normalized_momentum_gate", 2.0e-3),
        lambda r: r["agent4_power_law_force_audit"].__setitem__("independent_operator", "centered_fd2"),
        lambda r: r["agent4_power_law_force_audit"].__setitem__("authorized_as_ns_correction_target", True),
    ],
)
def test_forbidden_promotion_or_identity_drift_fails_closed(mutate):
    reg = a5.build_registration()
    bad = _mutated(reg, mutate)
    with pytest.raises(ValueError):
        a5.validate_registration(bad)


def test_registration_json_is_deterministic():
    first = a5.registration_json()
    second = a5.registration_json()
    assert first == second
    assert '"pde_validated": false' in first
    assert '"current_nonlinear_radial_force_through_power_law_materialized": true' in first
