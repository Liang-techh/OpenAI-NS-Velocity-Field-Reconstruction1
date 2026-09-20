from __future__ import annotations

from dataclasses import replace
import inspect

import pytest

from openai_ns_reconstruction.kokuno_finite_correction_stage import (
    FiniteCorrectionContractError,
    ValidationPartition,
)
from openai_ns_reconstruction.kokuno_finite_correction_fixed_physical_contract import (
    CR002_AUDIT_BLOB,
    CR002_CONFIG_BLOB,
    CR002_PHYSICAL_CONTRACT_HEAD,
    CR002_PHYSICAL_CONTRACT_PR,
    PARENT_AGENT3_HEAD,
    PARENT_AGENT3_PR,
    PARENT_AGENT3_SOURCE_BLOB,
    _FixedContractMechanicsBackend,
    build_mechanics_report,
    build_physical_contract_provenance,
    run_fixed_physical_contract_finite_correction_stage,
    truth_boundary,
    VelocityOnlyCorrectionBinding,
)


def _partitions() -> tuple[ValidationPartition, ValidationPartition]:
    return (
        ValidationPartition("held-in", ("train-0", "train-1", "train-2")),
        ValidationPartition("held-out", ("test-0", "test-1", "test-2")),
    )


def test_fixed_contract_mechanics_stage_preserves_identity_and_contracts():
    held_in, held_out = _partitions()
    result = run_fixed_physical_contract_finite_correction_stage(
        _FixedContractMechanicsBackend(),
        {"amplitude": 1.0},
        {"delta": -0.5},
        held_in,
        held_out,
    )
    receipt = result.receipt
    assert receipt.parent_stage.stage_accepted is True
    assert receipt.stage_accepted is True
    assert receipt.rejection_reasons == ()
    assert result.candidate_after_if_accepted == {"amplitude": 0.5}
    assert receipt.parent_stage.held_in_momentum_contraction_factor == pytest.approx(0.5)
    assert receipt.parent_stage.held_out_momentum_contraction_factor == pytest.approx(0.5)
    assert receipt.physical_contract_identity_preserved is True
    assert receipt.residual_evaluations_bound_to_physical_contract is True
    assert (
        receipt.physical_contract_before.physical_contract_sha256
        == receipt.physical_contract_after.physical_contract_sha256
        == receipt.held_in_before_contract_sha256
        == receipt.held_out_before_contract_sha256
        == receipt.held_in_after_contract_sha256
        == receipt.held_out_after_contract_sha256
    )
    assert receipt.mechanics_only is True
    assert receipt.fixed_physical_contract_residual_contraction_verified is False
    assert receipt.pde_validated is False


def test_secret_forcing_parameter_drift_rejects_parent_accepted_stage():
    class DriftBackend(_FixedContractMechanicsBackend):
        def physical_contract_provenance(self, candidate):
            base = super().physical_contract_provenance(candidate)
            amplitude = float(candidate["amplitude"])
            forcing_sha = "forcing-before" if amplitude >= 0.75 else "forcing-after"
            return build_physical_contract_provenance(
                contract_id=base.contract_id,
                viscosity=base.viscosity,
                physical_domain=base.physical_domain,
                evaluation_box=base.evaluation_box,
                time_interval=base.time_interval,
                support_condition=base.support_condition,
                nontriviality_normalization_id=base.nontriviality_normalization_id,
                residual_normalization_id=base.residual_normalization_id,
                matched_pressure_identity_sha256=base.matched_pressure_identity_sha256,
                restricted_forcing_family_id=base.restricted_forcing_family_id,
                restricted_forcing_parameters_sha256=forcing_sha,
                restricted_forcing_preregistration_sha256=base.restricted_forcing_preregistration_sha256,
            )

    held_in, held_out = _partitions()
    result = run_fixed_physical_contract_finite_correction_stage(
        DriftBackend(),
        {"amplitude": 1.0},
        {"delta": -0.5},
        held_in,
        held_out,
    )
    assert result.receipt.parent_stage.stage_accepted is True
    assert result.receipt.stage_accepted is False
    assert result.candidate_after_if_accepted is None
    assert result.receipt.physical_contract_identity_preserved is False
    assert "physical_contract_identity_changed" in result.receipt.rejection_reasons
    assert result.receipt.fixed_physical_contract_residual_contraction_verified is False


def test_residual_binding_mismatch_fails_closed_before_contraction_claim():
    class ResidualBindingDriftBackend(_FixedContractMechanicsBackend):
        def residual_physical_contract_sha256(self, candidate, partition):
            if partition.partition_id == "held-out":
                return "wrong-contract-sha"
            return super().residual_physical_contract_sha256(candidate, partition)

    held_in, held_out = _partitions()
    with pytest.raises(
        FiniteCorrectionContractError,
        match="residual evaluation is not bound",
    ):
        run_fixed_physical_contract_finite_correction_stage(
            ResidualBindingDriftBackend(),
            {"amplitude": 1.0},
            {"delta": -0.5},
            held_in,
            held_out,
        )


def test_joint_update_cannot_be_labeled_velocity_only():
    class JointUpdateBackend(_FixedContractMechanicsBackend):
        def correction_update_provenance(self, candidate, correction):
            source = self.candidate_provenance(candidate)
            return VelocityOnlyCorrectionBinding(
                correction_id="joint-update",
                source_candidate_sha256=source.candidate_sha256,
                update_kind="joint-physical-contract",
                pressure_changed=True,
                restricted_forcing_changed=True,
                physical_contract_changed=True,
            )

    held_in, held_out = _partitions()
    with pytest.raises(
        FiniteCorrectionContractError,
        match="not a velocity-only correction stage",
    ):
        run_fixed_physical_contract_finite_correction_stage(
            JointUpdateBackend(),
            {"amplitude": 1.0},
            {"delta": -0.5},
            held_in,
            held_out,
        )


def test_structured_contract_hash_rejects_laundered_identity():
    base = _FixedContractMechanicsBackend().physical_contract_provenance(
        {"amplitude": 1.0}
    )
    with pytest.raises(ValueError, match="does not match structured identity"):
        replace(base, physical_contract_sha256="laundered-contract-sha")


def test_pressure_identity_is_inside_contract_hash():
    base = _FixedContractMechanicsBackend().physical_contract_provenance(
        {"amplitude": 1.0}
    )
    changed = build_physical_contract_provenance(
        contract_id=base.contract_id,
        viscosity=base.viscosity,
        physical_domain=base.physical_domain,
        evaluation_box=base.evaluation_box,
        time_interval=base.time_interval,
        support_condition=base.support_condition,
        nontriviality_normalization_id=base.nontriviality_normalization_id,
        residual_normalization_id=base.residual_normalization_id,
        matched_pressure_identity_sha256="different-pressure",
        restricted_forcing_family_id=base.restricted_forcing_family_id,
        restricted_forcing_parameters_sha256=base.restricted_forcing_parameters_sha256,
        restricted_forcing_preregistration_sha256=base.restricted_forcing_preregistration_sha256,
    )
    assert changed.physical_contract_sha256 != base.physical_contract_sha256


def test_public_stage_api_has_no_contract_residual_forcing_or_threshold_knobs():
    parameters = inspect.signature(
        run_fixed_physical_contract_finite_correction_stage
    ).parameters
    assert list(parameters) == ["backend", "candidate", "correction", "held_in", "held_out"]
    forbidden = {
        "physical_contract",
        "physical_contract_sha256",
        "residual",
        "defect",
        "mean",
        "stress",
        "pressure",
        "forcing",
        "forcing_parameters",
        "target",
        "gain",
        "alpha",
        "damping",
        "viscosity",
        "nu",
        "delta_a",
        "normalized_score",
        "momentum_gate",
        "divergence_gate",
        "scientific_threshold",
    }
    assert forbidden.isdisjoint(parameters)


def test_truth_boundary_pins_cr002_and_does_not_promote_mechanics():
    report = build_mechanics_report()
    boundary = truth_boundary()
    assert PARENT_AGENT3_PR == 888
    assert PARENT_AGENT3_HEAD == "3fd33d3ca20d3ebed7dad872d3264036bdd9cb51"
    assert PARENT_AGENT3_SOURCE_BLOB == "4fc931fe7703b8b8e05efb41af957b54c9f7a4f9"
    assert CR002_PHYSICAL_CONTRACT_PR == 891
    assert CR002_PHYSICAL_CONTRACT_HEAD == "32c1c8f9779a1da077d0f21019e2658315e9eb62"
    assert CR002_CONFIG_BLOB == "a210483d413602f450cec3c769925c24a7a212da"
    assert CR002_AUDIT_BLOB == "c2c709e3f836975cd3e3d006cecf65b39591186a"
    assert report["mechanics_only"] is True
    assert report["candidate_residual_evidence"] is False
    assert report["mechanics_stage_receipt"]["stage_accepted"] is True
    assert boundary["structured_physical_contract_identity_defined"] is True
    assert boundary["nu_domain_support_time_normalization_bound"] is True
    assert boundary["matched_pressure_identity_bound"] is True
    assert boundary["restricted_forcing_family_and_parameter_identity_bound"] is True
    assert boundary["before_after_physical_contract_identity_required"] is True
    assert boundary["every_residual_partition_bound_to_physical_contract_sha"] is True
    assert boundary["velocity_only_update_typed"] is True
    assert boundary["joint_physical_contract_update_admitted_as_velocity_only"] is False
    assert boundary["current_scoped_transport_source_admitted"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["fixed_physical_contract_residual_contraction_verified"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed"] is False
    assert boundary["same_protocol_comparable_to_st006"] is False
    assert boundary["pde_validated"] is False
