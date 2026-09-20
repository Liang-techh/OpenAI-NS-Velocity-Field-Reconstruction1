from __future__ import annotations

"""Fixed-physical-contract guard for the Kokuno Agent-3 finite correction stage.

The parent stage (#888) already enforces complete-defect provenance, disjoint
held-in/held-out partitions, one residual-protocol SHA, and before/after
residual recomputation.  CR002 #891 identified a narrower attribution gap:
the same residual protocol does not prove that viscosity, domain/support,
normalization, matched pressure, or restricted-forcing parameters stayed
fixed.

This module adds exactly that missing guard.  It does not construct a new
pressure, forcing, defect, correction, or residual.  A correction accepted by
the parent stage is returned here only when the backend binds every residual
evaluation to the same explicit physical-contract SHA before and after and the
correction is typed as velocity-only.
"""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Protocol
import argparse
import hashlib
import json
import math

from openai_ns_reconstruction.kokuno_finite_correction_stage import (
    FiniteCorrectionBackend,
    FiniteCorrectionContractError,
    FiniteCorrectionStageReceipt,
    FiniteCorrectionStageResult,
    ValidationPartition,
    _MechanicsBackend,
    current_agent3_scoped_source_cycle_admission,
    run_finite_correction_stage,
)

PARENT_AGENT3_PR = 888
PARENT_AGENT3_HEAD = "3fd33d3ca20d3ebed7dad872d3264036bdd9cb51"
PARENT_AGENT3_SOURCE_BLOB = "4fc931fe7703b8b8e05efb41af957b54c9f7a4f9"
CR002_PHYSICAL_CONTRACT_PR = 891
CR002_PHYSICAL_CONTRACT_HEAD = "32c1c8f9779a1da077d0f21019e2658315e9eb62"
CR002_CONFIG_BLOB = "a210483d413602f450cec3c769925c24a7a212da"
CR002_AUDIT_BLOB = "c2c709e3f836975cd3e3d006cecf65b39591186a"
TASK = "KOKUNO-A3-FINITE-CORRECTION-PHYSICAL-CONTRACT-083"


def _sha256_payload(payload: dict[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class PhysicalContractProvenance:
    contract_id: str
    viscosity: float
    physical_domain: str
    evaluation_box: tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    time_interval: tuple[float, float]
    support_condition: str
    nontriviality_normalization_id: str
    residual_normalization_id: str
    matched_pressure_identity_sha256: str
    restricted_forcing_family_id: str
    restricted_forcing_parameters_sha256: str
    restricted_forcing_preregistration_sha256: str
    physical_contract_sha256: str

    def identity_payload(self) -> dict[str, Any]:
        return {
            "contract_id": self.contract_id,
            "viscosity": float(self.viscosity),
            "physical_domain": self.physical_domain,
            "evaluation_box": [list(v) for v in self.evaluation_box],
            "time_interval": list(self.time_interval),
            "support_condition": self.support_condition,
            "nontriviality_normalization_id": self.nontriviality_normalization_id,
            "residual_normalization_id": self.residual_normalization_id,
            "matched_pressure_identity_sha256": self.matched_pressure_identity_sha256,
            "restricted_forcing_family_id": self.restricted_forcing_family_id,
            "restricted_forcing_parameters_sha256": self.restricted_forcing_parameters_sha256,
            "restricted_forcing_preregistration_sha256": self.restricted_forcing_preregistration_sha256,
        }

    def __post_init__(self) -> None:
        text_fields = (
            "contract_id",
            "physical_domain",
            "support_condition",
            "nontriviality_normalization_id",
            "residual_normalization_id",
            "matched_pressure_identity_sha256",
            "restricted_forcing_family_id",
            "restricted_forcing_parameters_sha256",
            "restricted_forcing_preregistration_sha256",
            "physical_contract_sha256",
        )
        if any(not str(getattr(self, name)) for name in text_fields):
            raise ValueError("physical-contract identity fields must be nonempty")
        if not math.isfinite(float(self.viscosity)) or float(self.viscosity) <= 0.0:
            raise ValueError("viscosity must be finite and positive")
        if len(self.evaluation_box) != 3 or any(len(v) != 2 for v in self.evaluation_box):
            raise ValueError("evaluation_box must contain three finite intervals")
        for lo, hi in self.evaluation_box:
            if not (math.isfinite(float(lo)) and math.isfinite(float(hi)) and float(lo) < float(hi)):
                raise ValueError("evaluation_box intervals must be finite and increasing")
        if len(self.time_interval) != 2:
            raise ValueError("time_interval must contain two endpoints")
        t0, t1 = map(float, self.time_interval)
        if not (math.isfinite(t0) and math.isfinite(t1) and t0 < t1):
            raise ValueError("time_interval must be finite and increasing")
        expected = _sha256_payload(self.identity_payload())
        if self.physical_contract_sha256 != expected:
            raise ValueError("physical_contract_sha256 does not match structured identity")


def build_physical_contract_provenance(
    *,
    contract_id: str,
    viscosity: float,
    physical_domain: str,
    evaluation_box: tuple[tuple[float, float], tuple[float, float], tuple[float, float]],
    time_interval: tuple[float, float],
    support_condition: str,
    nontriviality_normalization_id: str,
    residual_normalization_id: str,
    matched_pressure_identity_sha256: str,
    restricted_forcing_family_id: str,
    restricted_forcing_parameters_sha256: str,
    restricted_forcing_preregistration_sha256: str,
) -> PhysicalContractProvenance:
    payload = {
        "contract_id": contract_id,
        "viscosity": float(viscosity),
        "physical_domain": physical_domain,
        "evaluation_box": [list(v) for v in evaluation_box],
        "time_interval": list(time_interval),
        "support_condition": support_condition,
        "nontriviality_normalization_id": nontriviality_normalization_id,
        "residual_normalization_id": residual_normalization_id,
        "matched_pressure_identity_sha256": matched_pressure_identity_sha256,
        "restricted_forcing_family_id": restricted_forcing_family_id,
        "restricted_forcing_parameters_sha256": restricted_forcing_parameters_sha256,
        "restricted_forcing_preregistration_sha256": restricted_forcing_preregistration_sha256,
    }
    return PhysicalContractProvenance(
        contract_id=contract_id,
        viscosity=float(viscosity),
        physical_domain=physical_domain,
        evaluation_box=evaluation_box,
        time_interval=time_interval,
        support_condition=support_condition,
        nontriviality_normalization_id=nontriviality_normalization_id,
        residual_normalization_id=residual_normalization_id,
        matched_pressure_identity_sha256=matched_pressure_identity_sha256,
        restricted_forcing_family_id=restricted_forcing_family_id,
        restricted_forcing_parameters_sha256=restricted_forcing_parameters_sha256,
        restricted_forcing_preregistration_sha256=restricted_forcing_preregistration_sha256,
        physical_contract_sha256=_sha256_payload(payload),
    )


@dataclass(frozen=True)
class VelocityOnlyCorrectionBinding:
    correction_id: str
    source_candidate_sha256: str
    update_kind: str
    pressure_changed: bool
    restricted_forcing_changed: bool
    physical_contract_changed: bool

    def __post_init__(self) -> None:
        if not self.correction_id or not self.source_candidate_sha256:
            raise ValueError("correction update identity must be nonempty")
        if self.update_kind not in {"velocity-only", "joint-physical-contract"}:
            raise ValueError("unsupported correction update kind")


class FixedPhysicalContractBackend(FiniteCorrectionBackend, Protocol):
    def physical_contract_provenance(self, candidate: Any) -> PhysicalContractProvenance: ...

    def correction_update_provenance(
        self, candidate: Any, correction: Any
    ) -> VelocityOnlyCorrectionBinding: ...

    def residual_physical_contract_sha256(
        self, candidate: Any, partition: ValidationPartition
    ) -> str: ...


@dataclass(frozen=True)
class FixedPhysicalContractStageReceipt:
    parent_stage: FiniteCorrectionStageReceipt
    physical_contract_before: PhysicalContractProvenance
    physical_contract_after: PhysicalContractProvenance | None
    correction_update: VelocityOnlyCorrectionBinding
    held_in_before_contract_sha256: str
    held_out_before_contract_sha256: str
    held_in_after_contract_sha256: str | None
    held_out_after_contract_sha256: str | None
    physical_contract_identity_preserved: bool
    residual_evaluations_bound_to_physical_contract: bool
    stage_accepted: bool
    rejection_reasons: tuple[str, ...]
    mechanics_only: bool
    fixed_physical_contract_residual_contraction_verified: bool
    pde_validated: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class FixedPhysicalContractStageResult:
    candidate_after_if_accepted: Any | None
    receipt: FixedPhysicalContractStageReceipt


def _require_extra_backend_contract(backend: Any) -> None:
    required = (
        "physical_contract_provenance",
        "correction_update_provenance",
        "residual_physical_contract_sha256",
    )
    missing = [name for name in required if not callable(getattr(backend, name, None))]
    if missing:
        raise FiniteCorrectionContractError(
            "backend is missing fixed-physical-contract methods: " + ", ".join(missing)
        )


def _get_contract(backend: Any, candidate: Any) -> PhysicalContractProvenance:
    contract = backend.physical_contract_provenance(candidate)
    if not isinstance(contract, PhysicalContractProvenance):
        raise FiniteCorrectionContractError(
            "backend must return PhysicalContractProvenance"
        )
    return contract


def _get_residual_contract_binding(
    backend: Any,
    candidate: Any,
    partition: ValidationPartition,
    expected_sha256: str,
) -> str:
    value = backend.residual_physical_contract_sha256(candidate, partition)
    if not isinstance(value, str) or not value:
        raise FiniteCorrectionContractError(
            "residual physical-contract binding must be a nonempty SHA identity"
        )
    if value != expected_sha256:
        raise FiniteCorrectionContractError(
            "residual evaluation is not bound to the candidate physical contract"
        )
    return value


def _validate_velocity_only_update(
    binding: VelocityOnlyCorrectionBinding,
    candidate_sha256: str,
) -> None:
    if not isinstance(binding, VelocityOnlyCorrectionBinding):
        raise FiniteCorrectionContractError(
            "backend must return VelocityOnlyCorrectionBinding"
        )
    if binding.source_candidate_sha256 != candidate_sha256:
        raise FiniteCorrectionContractError(
            "velocity-only correction source candidate hash mismatch"
        )
    if binding.update_kind != "velocity-only":
        raise FiniteCorrectionContractError(
            "joint physical-contract update is not a velocity-only correction stage"
        )
    if (
        binding.pressure_changed
        or binding.restricted_forcing_changed
        or binding.physical_contract_changed
    ):
        raise FiniteCorrectionContractError(
            "velocity-only correction provenance declares physical-contract drift"
        )


def run_fixed_physical_contract_finite_correction_stage(
    backend: FixedPhysicalContractBackend,
    candidate: Any,
    correction: Any,
    held_in: ValidationPartition,
    held_out: ValidationPartition,
) -> FixedPhysicalContractStageResult:
    """Run one parent finite stage, then fail closed on physical-contract drift.

    The caller still supplies only backend/candidate/correction and the two
    partition identities.  All physical-contract and residual bindings are
    obtained from the backend, not from caller-provided residual/forcing values.
    """

    _require_extra_backend_contract(backend)

    candidate_before = backend.candidate_provenance(candidate)
    contract_before = _get_contract(backend, candidate)
    update = backend.correction_update_provenance(candidate, correction)
    _validate_velocity_only_update(update, candidate_before.candidate_sha256)

    in_before_binding = _get_residual_contract_binding(
        backend, candidate, held_in, contract_before.physical_contract_sha256
    )
    out_before_binding = _get_residual_contract_binding(
        backend, candidate, held_out, contract_before.physical_contract_sha256
    )

    parent_result: FiniteCorrectionStageResult = run_finite_correction_stage(
        backend, candidate, correction, held_in, held_out
    )
    parent_receipt = parent_result.receipt

    if not parent_receipt.stage_accepted:
        receipt = FixedPhysicalContractStageReceipt(
            parent_stage=parent_receipt,
            physical_contract_before=contract_before,
            physical_contract_after=None,
            correction_update=update,
            held_in_before_contract_sha256=in_before_binding,
            held_out_before_contract_sha256=out_before_binding,
            held_in_after_contract_sha256=None,
            held_out_after_contract_sha256=None,
            physical_contract_identity_preserved=False,
            residual_evaluations_bound_to_physical_contract=False,
            stage_accepted=False,
            rejection_reasons=tuple(
                f"parent:{reason}" for reason in parent_receipt.rejection_reasons
            ),
            mechanics_only=parent_receipt.mechanics_only,
            fixed_physical_contract_residual_contraction_verified=False,
            pde_validated=False,
        )
        return FixedPhysicalContractStageResult(None, receipt)

    candidate_after = parent_result.candidate_after_if_accepted
    if candidate_after is None:
        raise FiniteCorrectionContractError(
            "parent accepted stage without returning candidate-after object"
        )

    contract_after = _get_contract(backend, candidate_after)
    in_after_binding = _get_residual_contract_binding(
        backend, candidate_after, held_in, contract_after.physical_contract_sha256
    )
    out_after_binding = _get_residual_contract_binding(
        backend, candidate_after, held_out, contract_after.physical_contract_sha256
    )

    rejection_reasons: list[str] = []
    identity_preserved = (
        contract_after.physical_contract_sha256
        == contract_before.physical_contract_sha256
    )
    if not identity_preserved:
        rejection_reasons.append("physical_contract_identity_changed")

    residuals_bound = bool(
        in_before_binding == contract_before.physical_contract_sha256
        and out_before_binding == contract_before.physical_contract_sha256
        and in_after_binding == contract_after.physical_contract_sha256
        and out_after_binding == contract_after.physical_contract_sha256
    )
    if not residuals_bound:
        rejection_reasons.append("residual_physical_contract_binding_changed")

    stage_accepted = bool(
        parent_receipt.stage_accepted
        and identity_preserved
        and residuals_bound
        and not rejection_reasons
    )
    scientific_fixed_contract_contraction = bool(
        stage_accepted and not parent_receipt.mechanics_only
    )

    receipt = FixedPhysicalContractStageReceipt(
        parent_stage=parent_receipt,
        physical_contract_before=contract_before,
        physical_contract_after=contract_after,
        correction_update=update,
        held_in_before_contract_sha256=in_before_binding,
        held_out_before_contract_sha256=out_before_binding,
        held_in_after_contract_sha256=in_after_binding,
        held_out_after_contract_sha256=out_after_binding,
        physical_contract_identity_preserved=identity_preserved,
        residual_evaluations_bound_to_physical_contract=residuals_bound,
        stage_accepted=stage_accepted,
        rejection_reasons=tuple(rejection_reasons),
        mechanics_only=parent_receipt.mechanics_only,
        fixed_physical_contract_residual_contraction_verified=scientific_fixed_contract_contraction,
        pde_validated=False,
    )
    return FixedPhysicalContractStageResult(
        candidate_after_if_accepted=candidate_after if stage_accepted else None,
        receipt=receipt,
    )


def _mechanics_contract() -> PhysicalContractProvenance:
    return build_physical_contract_provenance(
        contract_id="mechanics-fixed-contract-v1",
        viscosity=0.01,
        physical_domain="R^3",
        evaluation_box=((-2.0, 2.0), (-2.0, 2.0), (-2.0, 2.0)),
        time_interval=(0.25, 0.75),
        support_condition="mechanics-compact-support",
        nontriviality_normalization_id="mechanics-energy-normalization-v1",
        residual_normalization_id="mechanics-residual-normalization-v1",
        matched_pressure_identity_sha256="mechanics-pressure-sha",
        restricted_forcing_family_id="mechanics-fixed-forcing-family",
        restricted_forcing_parameters_sha256="mechanics-forcing-params-sha",
        restricted_forcing_preregistration_sha256="mechanics-forcing-preregistered-sha",
    )


class _FixedContractMechanicsBackend(_MechanicsBackend):
    def physical_contract_provenance(
        self, candidate: dict[str, float]
    ) -> PhysicalContractProvenance:
        return _mechanics_contract()

    def correction_update_provenance(
        self, candidate: dict[str, float], correction: dict[str, float]
    ) -> VelocityOnlyCorrectionBinding:
        source = self.candidate_provenance(candidate)
        return VelocityOnlyCorrectionBinding(
            correction_id="mechanics-half-amplitude",
            source_candidate_sha256=source.candidate_sha256,
            update_kind="velocity-only",
            pressure_changed=False,
            restricted_forcing_changed=False,
            physical_contract_changed=False,
        )

    def residual_physical_contract_sha256(
        self, candidate: dict[str, float], partition: ValidationPartition
    ) -> str:
        return self.physical_contract_provenance(candidate).physical_contract_sha256


def truth_boundary() -> dict[str, Any]:
    current = current_agent3_scoped_source_cycle_admission()
    return {
        "task": TASK,
        "parent_agent3_pr": PARENT_AGENT3_PR,
        "parent_agent3_head": PARENT_AGENT3_HEAD,
        "cr002_physical_contract_pr": CR002_PHYSICAL_CONTRACT_PR,
        "cr002_physical_contract_head": CR002_PHYSICAL_CONTRACT_HEAD,
        "structured_physical_contract_identity_defined": True,
        "physical_contract_sha_recomputed_from_structured_fields": True,
        "nu_domain_support_time_normalization_bound": True,
        "matched_pressure_identity_bound": True,
        "restricted_forcing_family_and_parameter_identity_bound": True,
        "restricted_forcing_preregistration_identity_bound": True,
        "before_after_physical_contract_identity_required": True,
        "every_residual_partition_bound_to_physical_contract_sha": True,
        "velocity_only_update_typed": True,
        "joint_physical_contract_update_admitted_as_velocity_only": False,
        "current_scoped_transport_source_admitted": current[
            "admitted_to_real_finite_correction_cycle"
        ],
        "real_candidate_finite_correction_cycle_run": False,
        "fixed_physical_contract_residual_contraction_verified": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "same_protocol_comparable_to_st006": False,
        "pde_validated": False,
    }


def build_mechanics_report() -> dict[str, Any]:
    held_in = ValidationPartition("held-in", ("train-0", "train-1", "train-2"))
    held_out = ValidationPartition("held-out", ("test-0", "test-1", "test-2"))
    result = run_fixed_physical_contract_finite_correction_stage(
        _FixedContractMechanicsBackend(),
        {"amplitude": 1.0},
        {"delta": -0.5},
        held_in,
        held_out,
    )
    return {
        "task": TASK,
        "mechanics_only": True,
        "candidate_residual_evidence": False,
        "current_scoped_source_admission": current_agent3_scoped_source_cycle_admission(),
        "mechanics_stage_receipt": result.receipt.to_dict(),
        "truth_boundary": truth_boundary(),
    }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_mechanics_report()
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(encoded, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")


if __name__ == "__main__":
    _main()
