"""Typed Agent-5 glue for assembling a same-cycle Kokuno composite.

This module is an integration seam, not a new mathematical lane.  It adapts
stage-level velocity / velocity-time-derivative providers into the exact
``VectorFieldProvider`` types consumed by Agent 3's same-cycle raw momentum
contract, while forwarding pressure and restricted forcing without fitting or
redefining either one.

The intended production wiring is

    Agent 1 leading bundle
      + Agent 2 frozen oscillatory bundle
      + optional Agent 3 correction bundle
      + matched pressure
      + independently restricted forcing
      -> one immutable same-cycle provider set
      -> Agent 3 raw-defect / finite-cycle machinery.

At this checkpoint the real Agent-1 global leading/matched-pressure handoff is
still unavailable, so the deterministic receipt below is a manufactured
integration regression only.  It must never be interpreted as a Kokuno
candidate residual or as the final independent PDE gate.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
    evaluate_disjoint_heldin_heldout,
)

TASK = "KOKUNO-A5-COMPOSITE-PROVIDER-GLUE-051"
SCHEMA = "kokuno-a5-composite-provider-glue-v1"
PARENT_A5_PR = 656
PARENT_A5_HEAD = "e95cdd5875c90824a55922cc27ee0d7c81da471f"
FINAL_MOMENTUM_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

VectorEvaluator = Callable[[float, float, float, float], Sequence[float]]
BundleEvaluator = Callable[[float, float, float, float], Mapping[str, Any]]


def _finite_vector(value: Any, *, label: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (3,) or not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must return one finite Cartesian 3-vector")
    return arr


@dataclass(frozen=True)
class StageVectorPair:
    """One immutable velocity / velocity_dt stage on a shared cycle identity."""

    identity: CycleIdentity
    stage: str
    source_ref: str
    velocity_evaluator: VectorEvaluator
    velocity_dt_evaluator: VectorEvaluator

    def __post_init__(self) -> None:
        if self.stage not in {"leading", "oscillatory", "correction"}:
            raise ValueError("stage must be leading, oscillatory, or correction")
        if not self.source_ref.strip():
            raise ValueError("source_ref must be non-empty")


def stage_from_bundle_evaluator(
    *,
    identity: CycleIdentity,
    stage: str,
    source_ref: str,
    bundle_evaluator: BundleEvaluator,
) -> StageVectorPair:
    """Adapt a public bundle exposing ``velocity`` and ``velocity_dt``.

    This is deliberately structural.  In particular it can consume Agent 2's
    public oscillatory bundle once that exact upstream module is present on the
    integration ancestry, without importing or copying Agent-2 mathematics here.
    """

    def velocity(x: float, y: float, z: float, t: float) -> np.ndarray:
        payload = bundle_evaluator(x, y, z, t)
        if "velocity" not in payload:
            raise ValueError("bundle evaluator omitted velocity")
        return _finite_vector(payload["velocity"], label=f"{source_ref}:velocity")

    def velocity_dt(x: float, y: float, z: float, t: float) -> np.ndarray:
        payload = bundle_evaluator(x, y, z, t)
        if "velocity_dt" not in payload:
            raise ValueError("bundle evaluator omitted velocity_dt")
        return _finite_vector(
            payload["velocity_dt"], label=f"{source_ref}:velocity_dt"
        )

    return StageVectorPair(
        identity=identity,
        stage=stage,
        source_ref=source_ref,
        velocity_evaluator=velocity,
        velocity_dt_evaluator=velocity_dt,
    )


def stage_from_vector_evaluators(
    *,
    identity: CycleIdentity,
    stage: str,
    source_ref: str,
    velocity_evaluator: VectorEvaluator,
    velocity_dt_evaluator: VectorEvaluator,
) -> StageVectorPair:
    """Typed adapter for stage providers that already expose separate callables."""
    return StageVectorPair(
        identity=identity,
        stage=stage,
        source_ref=source_ref,
        velocity_evaluator=velocity_evaluator,
        velocity_dt_evaluator=velocity_dt_evaluator,
    )


@dataclass(frozen=True)
class CompositeProviderSet:
    """Same-cycle provider set ready for the Agent-3 raw-defect evaluator."""

    identity: CycleIdentity
    velocity: VectorFieldProvider
    velocity_dt: VectorFieldProvider
    pressure: ScalarFieldProvider
    restricted_forcing: RestrictedForcingProvider
    stage_source_refs: tuple[str, ...]
    correction_included: bool

    def receipt(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "identity": {
                "cycle_id": self.identity.cycle_id,
                "cycle_index": self.identity.cycle_index,
                "state_token": self.identity.state_token,
            },
            "stage_source_refs": list(self.stage_source_refs),
            "pressure_source_ref": self.pressure.source_ref,
            "restricted_forcing_source_ref": self.restricted_forcing.source_ref,
            "restricted_forcing_receipt": self.restricted_forcing.restriction_receipt,
            "correction_included": self.correction_included,
        }


def assemble_composite_provider_set(
    *,
    identity: CycleIdentity,
    leading: StageVectorPair,
    oscillatory: StageVectorPair,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    correction: StageVectorPair | None = None,
) -> CompositeProviderSet:
    """Assemble stage velocity fields without accepting a residual or target.

    The function is intentionally strict about stage roles and one-cycle
    identity.  Pressure and forcing are forwarded as already-typed providers;
    this layer cannot fit either one and cannot replace forcing with a computed
    residual.
    """
    if leading.stage != "leading":
        raise ValueError("leading provider must declare stage='leading'")
    if oscillatory.stage != "oscillatory":
        raise ValueError("oscillatory provider must declare stage='oscillatory'")
    if correction is not None and correction.stage != "correction":
        raise ValueError("correction provider must declare stage='correction'")

    named_identities = {
        "leading": leading.identity,
        "oscillatory": oscillatory.identity,
        "pressure": pressure.identity,
        "restricted_forcing": restricted_forcing.identity,
    }
    if correction is not None:
        named_identities["correction"] = correction.identity
    for name, other in named_identities.items():
        if other != identity:
            raise ValueError(f"{name} does not share the requested cycle identity")

    stages = (leading, oscillatory) + (() if correction is None else (correction,))
    source_refs = tuple(stage.source_ref for stage in stages)
    if len(set(source_refs)) != len(source_refs):
        raise ValueError("stage source_ref values must be distinct")

    def total_velocity(x: float, y: float, z: float, t: float) -> np.ndarray:
        values = [
            _finite_vector(
                stage.velocity_evaluator(x, y, z, t),
                label=f"{stage.source_ref}:velocity",
            )
            for stage in stages
        ]
        return np.sum(np.vstack(values), axis=0)

    def total_velocity_dt(x: float, y: float, z: float, t: float) -> np.ndarray:
        values = [
            _finite_vector(
                stage.velocity_dt_evaluator(x, y, z, t),
                label=f"{stage.source_ref}:velocity_dt",
            )
            for stage in stages
        ]
        return np.sum(np.vstack(values), axis=0)

    joined = "+".join(source_refs)
    return CompositeProviderSet(
        identity=identity,
        velocity=VectorFieldProvider(
            identity=identity,
            source_ref=f"kokuno-composite-velocity[{joined}]",
            evaluator=total_velocity,
        ),
        velocity_dt=VectorFieldProvider(
            identity=identity,
            source_ref=f"kokuno-composite-velocity-dt[{joined}]",
            evaluator=total_velocity_dt,
        ),
        pressure=pressure,
        restricted_forcing=restricted_forcing,
        stage_source_refs=source_refs,
        correction_included=correction is not None,
    )


def evaluate_composite_raw_defect(
    providers: CompositeProviderSet,
    *,
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    viscosity: float,
    spatial_step: float,
) -> dict[str, Any]:
    """Route an assembled provider set through Agent 3's raw-defect contract."""
    result = evaluate_disjoint_heldin_heldout(
        providers.velocity,
        providers.velocity_dt,
        providers.pressure,
        providers.restricted_forcing,
        held_in_points,
        held_out_points,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    return {
        "assembly": providers.receipt(),
        "held_in": result["held_in"].to_receipt(),
        "held_out": result["held_out"].to_receipt(),
    }


def truth_boundary() -> dict[str, Any]:
    return {
        "typed_composite_velocity_provider_executable": True,
        "leading_plus_oscillatory_assembly_executable": True,
        "optional_correction_slot_executable": True,
        "bundle_velocity_velocity_dt_adapter_executable": True,
        "pressure_forwarded_without_fit": True,
        "restricted_forcing_forwarded_without_redefinition": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "real_agent1_global_leading_handoff_available_here": False,
        "real_agent1_matched_pressure_handoff_available_here": False,
        "restricted_forcing_semantics_independently_validated_here": False,
        "real_full_candidate_instantiated": False,
        "real_full_candidate_defect_consumed": False,
        "full_same_cycle_requested_stress_materialized": False,
        "correction_ready": False,
        "candidate_artifact_instantiated": False,
        "velocity_export_ready": False,
        "formal_full_domain_pde_gate_assessed": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "final_normalized_momentum_gate": FINAL_MOMENTUM_GATE,
        "final_normalized_divergence_gate": FINAL_DIVERGENCE_GATE,
    }


def deterministic_receipt() -> dict[str, Any]:
    """Manufactured integration regression; never a candidate science result."""
    identity = CycleIdentity("a5-composite-glue-regression", 0, "manufactured")

    leading = stage_from_vector_evaluators(
        identity=identity,
        stage="leading",
        source_ref="manufactured:leading-u=(y,0,0)",
        velocity_evaluator=lambda x, y, z, t: (y, 0.0, 0.0),
        velocity_dt_evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
    )

    def oscillatory_bundle(x: float, y: float, z: float, t: float) -> dict[str, Any]:
        return {"velocity": (0.0, x, 0.0), "velocity_dt": (0.0, 0.0, 0.0)}

    oscillatory = stage_from_bundle_evaluator(
        identity=identity,
        stage="oscillatory",
        source_ref="manufactured:oscillatory-u=(0,x,0)",
        bundle_evaluator=oscillatory_bundle,
    )
    pressure = ScalarFieldProvider(
        identity=identity,
        source_ref="manufactured:p=2x+3z",
        evaluator=lambda x, y, z, t: 2.0 * x + 3.0 * z,
    )
    forcing = RestrictedForcingProvider(
        identity=identity,
        source_ref="manufactured:f=0",
        restriction_receipt="manufactured-zero-force-fixed-before-residual",
        evaluator=lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    providers = assemble_composite_provider_set(
        identity=identity,
        leading=leading,
        oscillatory=oscillatory,
        pressure=pressure,
        restricted_forcing=forcing,
    )

    probe = (-0.31, 0.27, 0.14, 0.37)
    expected_velocity = np.asarray((probe[1], probe[0], 0.0), dtype=float)
    assembly_error = float(
        np.max(np.abs(np.asarray(providers.velocity.evaluator(*probe)) - expected_velocity))
    )

    held_in = np.asarray(
        [
            (-0.31, 0.27, 0.14, 0.37),
            (0.22, -0.19, -0.33, 0.63),
        ],
        dtype=float,
    )
    held_out = np.asarray(
        [
            (0.17, 0.41, -0.21, 0.44),
            (-0.26, -0.38, 0.29, 0.58),
        ],
        dtype=float,
    )
    defect = evaluate_composite_raw_defect(
        providers,
        held_in_points=held_in,
        held_out_points=held_out,
        viscosity=0.01,
        spatial_step=0.005,
    )

    expected_held_in = np.asarray(
        [(x + 2.0, y, 3.0) for x, y, z, t in held_in], dtype=float
    )
    actual_held_in = np.asarray(defect["held_in"]["residual_values"], dtype=float)
    defect_error = float(np.max(np.abs(actual_held_in - expected_held_in)))

    return {
        "task": TASK,
        "schema": SCHEMA,
        "provenance": {
            "parent_a5_pr": PARENT_A5_PR,
            "parent_a5_head": PARENT_A5_HEAD,
            "agent3_same_cycle_contract_present_on_parent": True,
            "agent3_finite_step_contract_present_on_parent": True,
        },
        "manufactured_regression": {
            "assembly_max_abs_error": assembly_error,
            "same_cycle_defect_max_abs_error": defect_error,
            "held_in_vector_rms": defect["held_in"]["vector_rms"],
            "held_out_vector_rms": defect["held_out"]["vector_rms"],
            "restricted_forcing_receipt": forcing.restriction_receipt,
        },
        "truth_boundary": truth_boundary(),
        "scientific_scope": (
            "typed integration regression only; no real Agent-1 leading/matched-pressure "
            "handoff, no real full-candidate defect, no correction readiness, and no "
            "independent normalized PDE validation"
        ),
    }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
