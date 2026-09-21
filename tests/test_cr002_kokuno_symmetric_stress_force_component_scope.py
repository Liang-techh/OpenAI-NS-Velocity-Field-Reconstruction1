from __future__ import annotations

import pytest

from openai_ns_reconstruction.audit_kokuno_symmetric_stress_force_component_scope import (
    assert_contract,
    audit_contract,
    load_contract,
    mechanics_force,
    mutated,
)


def test_canonical_contract_passes() -> None:
    assert_contract()


def test_mechanics_witness_has_three_nonzero_force_components_from_two_scalars() -> None:
    force = mechanics_force(2.0, 3.0)
    assert force == {"radial": 2.0, "tangential": 8.0, "axial": 6.0}
    assert all(abs(value) > 0.0 for value in force.values())


@pytest.mark.parametrize(
    ("path", "value", "needle"),
    [
        (
            ("representation_resolution", "source_specific_independent_third_radial_stress_channel_required"),
            True,
            "redundant third radial stress",
        ),
        (
            ("representation_resolution", "agent3_982_current_full_three_component_force_vector_materialized"),
            True,
            "full force vector",
        ),
        (
            ("representation_resolution", "agent3_982_complete_ns_defect_materialized"),
            True,
            "complete defect",
        ),
        (
            ("representation_resolution", "agent3_982_authorized_correction_target"),
            True,
            "prematurely authorized",
        ),
        (
            ("representation_resolution", "agent3_982_cartesian_correction_velocity_materialized"),
            True,
            "correction velocity",
        ),
        (
            ("provenance", "public_source_fact", "independent_third_scalar_radial_stress_present"),
            True,
            "invented third source scalar stress",
        ),
        (("truth_states", "paper_exact"), True, "paper_exact"),
        (("truth_states", "openai_field_identified"), True, "openai_field_identified"),
        (("truth_states", "kokuno_pde_validated"), True, "kokuno_pde_validated"),
        (("cr001_snapshot", "momentum_max_gate"), 0.01, "momentum max gate changed"),
        (
            ("cr001_snapshot", "residual_defined_pointwise_free_force_allowed"),
            True,
            "free residual forcing enabled",
        ),
        (("cr001_snapshot", "candidate_collapse_allowed"), True, "candidate collapse enabled"),
        (
            ("cr001_snapshot", "post_hoc_threshold_relaxation_allowed"),
            True,
            "post-hoc threshold relaxation enabled",
        ),
    ],
)
def test_fail_closed_mutations(path: tuple[str, ...], value: object, needle: str) -> None:
    contract = load_contract()
    errors = audit_contract(mutated(contract, path, value))
    assert any(needle in error for error in errors), errors


def test_mechanics_witness_cannot_be_laundered_into_source_fact() -> None:
    contract = load_contract()
    changed = mutated(contract, ("mechanics_only_witness", "classification"), "public_source_fact")
    errors = audit_contract(changed)
    assert any("mechanics witness provenance laundering" in error for error in errors)


def test_generic_cr002_979_witness_cannot_override_source_specific_map() -> None:
    contract = load_contract()
    changed = mutated(
        contract,
        ("representation_resolution", "generic_vector_witness_is_not_the_registered_source_tensor_representation"),
        False,
    )
    errors = audit_contract(changed)
    assert any("generic witness laundered into source tensor" in error for error in errors)


def test_corrected_reader_cannot_be_promoted_to_independently_verified_original_paper() -> None:
    contract = load_contract()
    changed = mutated(
        contract,
        ("provenance", "public_source_fact", "scope"),
        "verified exact original OpenAI paper formula",
    )
    errors = audit_contract(changed)
    assert any("paper-exact fact" in error for error in errors)


def test_complete_force_vector_must_remain_pending() -> None:
    contract = load_contract()
    changed = dict(contract)
    changed["provenance"] = dict(contract["provenance"])
    changed["provenance"]["pending_unknown"] = [
        item
        for item in contract["provenance"]["pending_unknown"]
        if "complete three-component div(T) force vector" not in item
    ]
    errors = audit_contract(changed)
    assert any("full force vector no longer pending" in error for error in errors)


def test_forbidden_promotion_firewall_cannot_be_weakened() -> None:
    contract = load_contract()
    changed = dict(contract)
    changed["forbidden_promotions"] = [
        item
        for item in contract["forbidden_promotions"]
        if item != "current radial force -> complete NS defect"
    ]
    errors = audit_contract(changed)
    assert any("forbidden-promotion firewall weakened" in error for error in errors)


def test_eq45_delivery_cannot_be_downgraded_by_kokuno_scope() -> None:
    contract = load_contract()
    changed = mutated(contract, ("truth_states", "canonical_eq45_velocity_export_ready"), False)
    errors = audit_contract(changed)
    assert any("incorrectly downgraded Eq45 delivery" in error for error in errors)
