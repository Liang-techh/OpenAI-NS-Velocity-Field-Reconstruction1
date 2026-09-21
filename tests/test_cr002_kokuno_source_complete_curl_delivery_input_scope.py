from __future__ import annotations

from copy import deepcopy

import numpy as np
import pytest

from openai_ns_reconstruction.audit_kokuno_source_complete_curl_delivery_input_scope import (
    audit_scope,
    load_contract,
)
from openai_ns_reconstruction.kokuno_source_normalized_complete_curl_reference import (
    complete_harmonic_amplitude_from_coefficient_jet,
    source_reference_contract,
)


def test_baseline_scope_passes_against_live_repository() -> None:
    assert audit_scope() == []


def test_exact_upstream_contract_keeps_directional_jets_external_and_promotions_false() -> None:
    upstream = source_reference_contract()
    assert upstream["source_normalized_curl_formula_executable"] is True
    assert upstream["source_complete_harmonic_amplitude_algebra_executable"] is True
    assert upstream["source_annular_domain_requires_positive_radius"] is True
    assert upstream["source_directional_derivatives_supplied_by_caller"] is True
    assert upstream["current_runtime_phase_source_exact"] is False
    assert upstream["current_runtime_support_source_exact"] is False
    assert upstream["current_runtime_complete_curl_source_equivalence_verified"] is False
    assert upstream["source_to_runtime_parameter_map_complete"] is False
    assert upstream["paper_exact"] is False
    assert upstream["pde_validated"] is False


def test_mechanics_witness_shows_missing_jet_realization_is_output_relevant() -> None:
    # Autonomous mechanics only: these numbers are not Kokuno/OpenAI parameters.
    common = dict(
        radius=1.0,
        n_phi=np.array([1.0, 0.0, 0.0]),
        t_m=np.array([0.0, 1.0, 0.0], dtype=np.complex128),
        k=2.0,
        m=1,
    )
    first = complete_harmonic_amplitude_from_coefficient_jet(
        **common,
        coefficient_dr=np.zeros(3, dtype=np.complex128),
        coefficient_dz=np.zeros(3, dtype=np.complex128),
    )
    second = complete_harmonic_amplitude_from_coefficient_jet(
        **common,
        coefficient_dr=np.zeros(3, dtype=np.complex128),
        coefficient_dz=np.array([0.0, 1.0, 0.0], dtype=np.complex128),
    )

    assert np.allclose(first.amplitude, np.array([0.0, 1.0, 0.0]))
    assert np.allclose(second.amplitude, np.array([-1.0, 1.0, 0.0]))
    assert not np.allclose(first.amplitude, second.amplitude)


@pytest.mark.parametrize(
    ("section", "key", "bad_value"),
    [
        (
            "delivery_input_boundary",
            "source_reference_self_contained_cartesian_velocity_provider",
            True,
        ),
        ("delivery_input_boundary", "source_reference_accepts_x_y_z_t_only", True),
        ("delivery_input_boundary", "source_reference_defines_global_axis_extension", True),
        (
            "delivery_input_boundary",
            "source_reference_establishes_axis_safe_global_cartesian_runtime",
            True,
        ),
        (
            "delivery_input_boundary",
            "source_reference_establishes_project_domain_totality",
            True,
        ),
        (
            "delivery_input_boundary",
            "source_reference_establishes_identity_preserving_velocity_save_load",
            True,
        ),
        ("delivery_input_boundary", "source_reference_velocity_export_ready", True),
        (
            "upstream_source_reference",
            "source_directional_derivatives_supplied_by_caller",
            False,
        ),
        (
            "upstream_source_reference",
            "current_runtime_complete_curl_source_equivalence_verified",
            True,
        ),
        ("upstream_source_reference", "source_to_runtime_parameter_map_complete", True),
        (
            "machine_locked_nonimplications",
            "manufactured_curl_or_divergence_identity_implies_runtime_source_equivalence",
            True,
        ),
        (
            "machine_locked_nonimplications",
            "source_formula_reference_implies_visual_correspondence",
            True,
        ),
        (
            "machine_locked_nonimplications",
            "source_formula_reference_implies_pde_validation",
            True,
        ),
        (
            "machine_locked_nonimplications",
            "source_formula_reference_implies_paper_exact",
            True,
        ),
        (
            "machine_locked_nonimplications",
            "source_formula_reference_implies_openai_field_identity",
            True,
        ),
        ("kokuno_source_reference_state", "velocity_export_ready", True),
        ("kokuno_source_reference_state", "visual_correspondence_verified", True),
        ("kokuno_source_reference_state", "pde_validated", True),
        ("kokuno_source_reference_state", "paper_exact", True),
        ("kokuno_source_reference_state", "openai_field_identified", True),
        ("canonical_delivery_independence", "velocity_export_ready", False),
    ],
)
def test_forbidden_scope_promotions_fail_closed(
    section: str, key: str, bad_value: object
) -> None:
    contract = load_contract()
    contract[section][key] = bad_value
    errors = audit_scope(contract, verify_live_files=False)
    assert errors, (section, key, bad_value)


@pytest.mark.parametrize(
    ("key", "bad_value"),
    [
        ("momentum_max", 0.01),
        ("momentum_l2", 0.01),
        ("divergence_max", 1e-4),
        ("divergence_l2", 1e-4),
        ("residual_defined_free_forcing_forbidden", False),
        ("candidate_collapse_forbidden", False),
        ("post_hoc_threshold_relaxation_forbidden", False),
    ],
)
def test_cr001_mutations_fail_closed(key: str, bad_value: object) -> None:
    contract = load_contract()
    contract["cr001_snapshot"][key] = bad_value
    assert audit_scope(contract, verify_live_files=False)


def test_four_way_provenance_cannot_be_collapsed_or_renamed() -> None:
    contract = load_contract()
    moved = contract["four_way_provenance"].pop("pending_or_unknown")
    contract["four_way_provenance"]["public_source_facts"].extend(moved)
    assert audit_scope(contract, verify_live_files=False)


def test_autonomous_mechanics_witness_cannot_be_promoted_to_source_or_candidate_evidence() -> None:
    contract = load_contract()
    contract["scope_logic_witness"]["classification"] = "public_source_fact"
    contract["scope_logic_witness"]["not_a_public_source_parameter_claim"] = False
    contract["scope_logic_witness"]["not_candidate_numerical_evidence"] = False
    errors = audit_scope(contract, verify_live_files=False)
    assert len(errors) >= 3


def test_formula_execution_and_velocity_delivery_must_remain_independent() -> None:
    contract = load_contract()
    boundary = contract["delivery_input_boundary"]
    assert boundary["source_formula_executable_can_be_true_while_velocity_delivery_remains_false"]
    assert boundary["source_reference_self_contained_cartesian_velocity_provider"] is False
    assert boundary["source_reference_velocity_export_ready"] is False
    assert contract["canonical_delivery_independence"]["velocity_export_ready"] is True
    assert audit_scope(contract, verify_live_files=False) == []


def test_missing_forbidden_promotion_firewall_is_rejected() -> None:
    contract = load_contract()
    contract["forbidden_promotions"] = [
        value
        for value in contract["forbidden_promotions"]
        if "downgrade canonical Eq45 velocity_export_ready" not in value
    ]
    assert audit_scope(contract, verify_live_files=False)


def test_source_reference_axis_rejection_is_not_a_global_project_axis_claim() -> None:
    contract = load_contract()
    assert contract["upstream_source_reference"]["source_annular_domain_requires_positive_radius"]
    assert (
        contract["machine_locked_nonimplications"][
            "annular_fail_closed_reference_implies_project_velocity_cannot_be_axis_safe"
        ]
        is False
    )
    assert (
        contract["delivery_input_boundary"][
            "source_reference_establishes_axis_safe_global_cartesian_runtime"
        ]
        is False
    )
