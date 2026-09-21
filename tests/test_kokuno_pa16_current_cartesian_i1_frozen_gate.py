import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa16_current_cartesian_i1_frozen_gate import (
    I1_CLOSURE_TOLERANCE,
    KokunoPA16CurrentCartesianI1FrozenGate,
)


@pytest.fixture(scope="module")
def candidate():
    return KokunoPA16CurrentCartesianI1FrozenGate()


def test_frozen_gate_matches_admitted_parent_without_changing_profiles(candidate):
    assert candidate.closure_tolerance == I1_CLOSURE_TOLERANCE == 5.0e-7
    assert candidate.candidate.closure_tolerance == I1_CLOSURE_TOLERANCE
    assert candidate.i1_family.closure_tolerance == I1_CLOSURE_TOLERANCE

    X = candidate.i1_family.repair.X_1
    eta = np.asarray([-0.45, 0.0, 0.37, 0.8])
    wrapped = candidate.similarity_profile_values(X, eta)
    parent = candidate.candidate.similarity_profile_values(X, eta)
    for key in (
        "F_current_i1_repair",
        "U_current_i1_repair",
        "E_current_i1_repair",
        "M_over_X_current_i1_repair",
        "M_eta_over_X_current_i1_repair",
        "v0_current_i1_repair",
    ):
        np.testing.assert_array_equal(wrapped[key], parent[key])

    np.testing.assert_array_equal(
        candidate.velocity(0.0, 0.0, 0.0, 0.5),
        candidate.candidate.velocity(0.0, 0.0, 0.0, 0.5),
    )


def test_configuration_roundtrip_binds_frozen_gate(candidate, tmp_path):
    path = tmp_path / "current_i1_frozen_gate.json"
    payload = candidate.save_configuration(path)
    assert payload["frozen_i1_closure_gate"] == {
        "closure_tolerance": I1_CLOSURE_TOLERANCE,
        "mutable": False,
    }
    assert (
        payload["parent_i1_candidate"]["i1_repair"]["closure_tolerance"]
        == I1_CLOSURE_TOLERANCE
    )

    restored = KokunoPA16CurrentCartesianI1FrozenGate.load_configuration(path)
    assert restored.configuration() == candidate.configuration()
    assert restored.semantic_sha256 == candidate.semantic_sha256


def test_serialized_parent_gate_mutation_5e7_to_1e5_fails_closed(candidate):
    payload = copy.deepcopy(candidate.configuration())
    payload["parent_i1_candidate"]["i1_repair"]["closure_tolerance"] = 1.0e-5
    with pytest.raises(ValueError, match="frozen I1 closure gate"):
        KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(payload)


def test_serialized_wrapper_gate_mutation_5e7_to_1e5_fails_closed(candidate):
    payload = copy.deepcopy(candidate.configuration())
    payload["frozen_i1_closure_gate"]["closure_tolerance"] = 1.0e-5
    with pytest.raises(ValueError, match="frozen I1 closure gate"):
        KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(payload)


def test_mutability_and_schema_mutations_fail_closed(candidate):
    mutable = copy.deepcopy(candidate.configuration())
    mutable["frozen_i1_closure_gate"]["mutable"] = True
    with pytest.raises(ValueError, match="immutable"):
        KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(mutable)

    wrong_schema = copy.deepcopy(candidate.configuration())
    wrong_schema["schema"] = "wrong"
    with pytest.raises(ValueError, match="schema"):
        KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(wrong_schema)


def test_report_and_truth_boundary_do_not_promote_science(candidate):
    report = candidate.report()
    assert report["velocity_equivalence_to_parent"] is True
    assert report["frozen_i1_closure_gate"]["closure_tolerance"] == 5.0e-7
    assert report["frozen_i1_closure_gate"]["mutable"] is False

    truth = candidate.truth_boundary
    assert truth["current_I1_closure_gate_frozen"] is True
    assert truth["configuration_mutation_of_I1_closure_gate_fails_closed"] is True
    assert truth["velocity_formula_changed_by_this_increment"] is False
    assert truth["PA17_coefficients_changed_by_this_increment"] is False
    for key in (
        "current_I2_overlay_applied",
        "current_I3_overlay_applied",
        "current_I4_overlay_applied",
        "outer_global_leading_velocity_materialized",
        "matched_global_pressure_materialized",
        "restricted_forcing_materialized",
        "heldout_ns_residual_assessed",
        "same_protocol_comparable_to_st006",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        assert truth[key] is False
