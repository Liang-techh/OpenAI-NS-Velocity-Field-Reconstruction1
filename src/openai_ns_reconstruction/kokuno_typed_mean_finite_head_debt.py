"""Apply the frozen finite-head factor directly to typed actual mean stress.

Agent-5 integration copy of the Agent-3 #686 seam.  The public API accepts only
raw same-cycle field providers plus physical radial geometry and derivative
settings.  It recomputes the actual momentum defect, projects the actual
mean channels, reconstructs compact radial stress, and applies the frozen
repository-autonomous finite-head factor

    Delta C_aut = -w_aut * requestedStress_actual.

The scalar factor is replayed through the pure Agent-5 adapter extracted from
Agent-3 #586, so this integration stack does not import the unrelated
oscillatory-candidate-specific Agent-2 graph.  No numerical formula or frozen
factor input is retuned.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .kokuno_autonomous_finite_head_factor_adapter import (
    AUTONOMOUS_FACTOR_ORIGIN_PR,
    AUTONOMOUS_FACTOR_REPLAY_PR,
    EXPECTED_AUTONOMOUS_MISSING_WEIGHT,
    replay_autonomous_finite_head_factor,
)
from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)
from .kokuno_typed_mean_radial_stress import (
    RadialMeanStressGeometry,
    TypedMeanRadialStressEvaluation,
    materialize_actual_mean_radial_stress,
)

TASK = "KOKUNO-A5-INTEGRATE-TYPED-MEAN-FINITE-HEAD-DEBT-054"
SCHEMA = "kokuno-a5-typed-mean-finite-head-debt-v1"
UPSTREAM_AGENT3_PR = 686
UPSTREAM_AGENT3_HEAD = "df05669aaee94cd8de2463a35c99a73c476118c1"
PARENT_A5_PR = 679
PARENT_A5_HEAD = "b97907632d5d9255045faa530eb00b8849cc25ac"


@dataclass(frozen=True)
class TypedMeanFiniteHeadDebtEvaluation:
    """Typed actual radial stress and its frozen autonomous finite-head debt."""

    identity: CycleIdentity
    radial_stress: TypedMeanRadialStressEvaluation
    autonomous_factor: dict[str, object]
    theta_debt: np.ndarray
    axial_debt: np.ndarray
    requested_stress_rms: float
    requested_stress_max_abs: float
    debt_rms: float
    debt_max_abs: float
    scalar_weight_rms_prediction: float
    identity_closure_max_abs: float

    def to_receipt(self) -> dict[str, object]:
        return {
            "identity": asdict(self.identity),
            "radial_stress": self.radial_stress.to_receipt(),
            "autonomous_factor": self.autonomous_factor,
            "finite_head_mean_debt": {
                "formula": "DeltaC_aut=-autonomous_missing_weight*requestedStress_actual",
                "ordering": ["theta_e2", "axial_e1"],
                "theta_e2_values": self.theta_debt.tolist(),
                "axial_e1_values": self.axial_debt.tolist(),
                "requested_stress_rms": self.requested_stress_rms,
                "requested_stress_max_abs": self.requested_stress_max_abs,
                "debt_rms": self.debt_rms,
                "debt_max_abs": self.debt_max_abs,
                "scalar_weight_rms_prediction": self.scalar_weight_rms_prediction,
                "identity_closure_max_abs": self.identity_closure_max_abs,
            },
        }


def materialize_typed_mean_finite_head_debt(
    velocity: VectorFieldProvider,
    velocity_dt: VectorFieldProvider,
    pressure: ScalarFieldProvider,
    restricted_forcing: RestrictedForcingProvider,
    geometry: RadialMeanStressGeometry,
    *,
    viscosity: float,
    spatial_step: float,
) -> TypedMeanFiniteHeadDebtEvaluation:
    """Map raw same-cycle fields to a frozen autonomous finite-head mean debt."""

    stress_eval = materialize_actual_mean_radial_stress(
        velocity,
        velocity_dt,
        pressure,
        restricted_forcing,
        geometry,
        viscosity=viscosity,
        spatial_step=spatial_step,
    )
    factor = replay_autonomous_finite_head_factor()
    weight = float(factor["autonomous_missing_weight"])
    if not math.isclose(
        weight,
        EXPECTED_AUTONOMOUS_MISSING_WEIGHT,
        rel_tol=0.0,
        abs_tol=5.0e-15,
    ):
        raise RuntimeError("frozen autonomous finite-head factor replay drifted")
    if factor["formal_theorem_machine_bump_identity_claimed"] is not False:
        raise RuntimeError("autonomous factor must not masquerade as theorem-machine bump")
    if factor["formal_missing_weight_equality_claimed"] is not False:
        raise RuntimeError("autonomous factor must not claim formal missingWeight equality")

    theta_stress = np.asarray(stress_eval.theta_stress["stress"], dtype=float)
    axial_stress = np.asarray(stress_eval.axial_stress["stress"], dtype=float)
    if theta_stress.ndim != 1 or axial_stress.shape != theta_stress.shape:
        raise RuntimeError("typed radial-stress channel shape contract changed")
    if not np.all(np.isfinite(theta_stress)) or not np.all(np.isfinite(axial_stress)):
        raise RuntimeError("typed radial stress contains non-finite values")

    stress_matrix = np.stack((theta_stress, axial_stress), axis=-1)
    theta_debt = -weight * theta_stress
    axial_debt = -weight * axial_stress
    debt_matrix = np.stack((theta_debt, axial_debt), axis=-1)
    closure = float(np.max(np.abs(debt_matrix + weight * stress_matrix)))
    requested_rms = float(np.sqrt(np.mean(np.sum(stress_matrix * stress_matrix, axis=-1))))
    debt_rms = float(np.sqrt(np.mean(np.sum(debt_matrix * debt_matrix, axis=-1))))
    requested_max = float(np.max(np.abs(stress_matrix)))
    debt_max = float(np.max(np.abs(debt_matrix)))
    predicted_rms = weight * requested_rms

    scale = max(1.0, debt_max)
    if closure > 128.0 * math.ulp(scale):
        raise ArithmeticError("finite-head scalar-weight identity lost floating-point closure")
    if not math.isclose(debt_rms, predicted_rms, rel_tol=2.0e-15, abs_tol=1.0e-15):
        raise ArithmeticError("finite-head debt RMS lost scalar-weight consistency")

    return TypedMeanFiniteHeadDebtEvaluation(
        identity=stress_eval.identity,
        radial_stress=stress_eval,
        autonomous_factor=factor,
        theta_debt=theta_debt,
        axial_debt=axial_debt,
        requested_stress_rms=requested_rms,
        requested_stress_max_abs=requested_max,
        debt_rms=debt_rms,
        debt_max_abs=debt_max,
        scalar_weight_rms_prediction=predicted_rms,
        identity_closure_max_abs=closure,
    )


def truth_boundary() -> dict[str, object]:
    signature = inspect.signature(materialize_typed_mean_finite_head_debt)
    forbidden = {
        "residual",
        "defect",
        "mean_source",
        "source",
        "stress",
        "requested_stress",
        "target",
        "factor",
        "missing_weight",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "prepared_n",
        "band",
        "coordinate_q",
    }
    return {
        "source_defect_formula": "u_t + (u dot grad)u + grad(p) - nu*Delta(u) - f",
        "actual_stress_provider": "kokuno_typed_mean_radial_stress:materialize_actual_mean_radial_stress",
        "upstream_agent3_pr": UPSTREAM_AGENT3_PR,
        "upstream_agent3_head": UPSTREAM_AGENT3_HEAD,
        "parent_a5_pr": PARENT_A5_PR,
        "parent_a5_head": PARENT_A5_HEAD,
        "autonomous_factor_origin_pr": AUTONOMOUS_FACTOR_ORIGIN_PR,
        "autonomous_factor_replay_pr": AUTONOMOUS_FACTOR_REPLAY_PR,
        "autonomous_factor_replay_provider": "kokuno_autonomous_finite_head_factor_adapter:replay_autonomous_finite_head_factor",
        "finite_head_debt_formula": "DeltaC_aut=-w_aut*requestedStress_actual",
        "raw_same_cycle_providers_required": True,
        "caller_supplied_precomputed_defect_allowed": False,
        "caller_supplied_mean_source_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_finite_head_factor_allowed": False,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "typed_actual_stress_fed_directly_to_frozen_finite_head_factor": True,
        "typed_actual_state_finite_head_debt_executable": True,
        "autonomous_factor_is_repository_engineering": True,
        "integration_adapter_only": True,
        "formal_theorem_machine_bump_identity_claimed": False,
        "formal_missing_weight_materialized": False,
        "formal_missing_weight_equality_claimed": False,
        "restricted_forcing_semantics_independently_validated_here": False,
        "real_full_candidate_defect_consumed": False,
        "full_same_cycle_composite_requested_stress_materialized": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed": False,
        "pde_validated": False,
    }


def _compact_profile(radius: float, *, center: float, halfwidth: float) -> float:
    s = (radius - center) / halfwidth
    if abs(s) >= 1.0:
        return 0.0
    return math.cos(0.5 * math.pi * s) ** 8


def deterministic_receipt() -> dict[str, object]:
    """Nonzero f=0 mechanics regression from raw fields through finite-head debt."""

    identity = CycleIdentity("analytic-a5-typed-mean-debt-v1", 0, "t0")
    center = 0.50
    halfwidth = 0.30

    def profile(x: float, y: float) -> float:
        return _compact_profile(math.hypot(x, y), center=center, halfwidth=halfwidth)

    velocity = VectorFieldProvider(
        identity,
        "analytic:u=t*psi(r)*(-y,x,1)",
        lambda x, y, z, t: (
            -t * y * profile(x, y),
            t * x * profile(x, y),
            t * profile(x, y),
        ),
    )
    velocity_dt = VectorFieldProvider(
        identity,
        "analytic:u_t=psi(r)*(-y,x,1)",
        lambda x, y, z, t: (
            -y * profile(x, y),
            x * profile(x, y),
            profile(x, y),
        ),
    )
    pressure = ScalarFieldProvider(identity, "analytic:p=0", lambda x, y, z, t: 0.0)
    forcing = RestrictedForcingProvider(
        identity,
        "analytic:f=0",
        "regression-zero-forcing; fixed independently of residual",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    geometry = RadialMeanStressGeometry(
        time=0.0,
        axial_z=0.13,
        radii=tuple(float(value) for value in np.linspace(0.20, 0.80, 49)),
        bump_center=center,
        bump_halfwidth=halfwidth,
        angular_count=32,
        phase=0.017,
    )
    evaluation = materialize_typed_mean_finite_head_debt(
        velocity,
        velocity_dt,
        pressure,
        forcing,
        geometry,
        viscosity=0.01,
        spatial_step=0.005,
    )
    weight = float(evaluation.autonomous_factor["autonomous_missing_weight"])
    ratio = evaluation.debt_rms / evaluation.requested_stress_rms
    return {
        "schema": SCHEMA,
        "task": TASK,
        "provenance": {
            "upstream_agent3_pr": UPSTREAM_AGENT3_PR,
            "upstream_agent3_head": UPSTREAM_AGENT3_HEAD,
            "parent_a5_pr": PARENT_A5_PR,
            "parent_a5_head": PARENT_A5_HEAD,
            "structural_source": "Kokuno corrected reconstruction finite-head mean-correction architecture",
            "repository_role": "A5 typed actual mean stress -> frozen autonomous finite-head debt integration",
            "autonomous_factor_origin_pr": AUTONOMOUS_FACTOR_ORIGIN_PR,
            "autonomous_factor_replay_pr": AUTONOMOUS_FACTOR_REPLAY_PR,
        },
        "analytic_regression": {
            "field": "u=t*psi(r)*(-y,x,1), u_t=psi(r)*(-y,x,1), p=0, f=0 at t=0",
            "expected_actual_defect": "R=psi(r)*(-y,x,1)",
            "autonomous_missing_weight": weight,
            "debt_to_requested_stress_rms_ratio": ratio,
            "evaluation": evaluation.to_receipt(),
        },
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
