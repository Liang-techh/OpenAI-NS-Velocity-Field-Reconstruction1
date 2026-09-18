"""Fail-closed Agent-3 bridge from formal signed mean debt to numeric correction input.

The pinned formal source now admits the finite-head identity

    meanBar(actualCross) - requestedStress
      = -missingWeight * requestedStress.

That is genuine source information, but :mod:`actual_signed_mean_defect` also
states explicitly that the requested stress, missing weight, cycle state, and
point remain opaque in Python.  A theorem identity is therefore not yet a
numeric covariance target and must not be substituted by a hand-filled
surrogate.

This module makes that boundary executable.  It consumes the theorem-grade
``ActualSignedMeanDefectAdmission`` directly, plus the downstream physical
signed-covariance rank receipt.  A numeric target can enter the next
``Delta C / epsilon`` inverse only when the formal admission itself reports
that the relevant actual values have been machine-materialized and the target
is bound to the same theorem application.  The current admission intentionally
fails that gate.

A PASS from this module would only mean that a real numeric signed-mean inverse
input is ready for the existing Agent-3 unit/budget/spacetime/radial/quadratic
guards.  It never, by itself, authorizes a finite correction-cycle rerun or a
Navier--Stokes residual claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from .actual_signed_mean_defect import ActualSignedMeanDefectAdmission


PHYSICAL_TARGET_UNITS = "physical phase-mean covariance/stress"
FORMAL_NUMERIC_EVIDENCE_KIND = "formal-machine-replay"


def _finite_vector(values: Sequence[float], name: str) -> tuple[float, ...]:
    array = np.asarray(tuple(values), dtype=float)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{name} must be a nonempty one-dimensional vector")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return tuple(float(value) for value in array)


@dataclass(frozen=True)
class ActualSignedMeanNumericTarget:
    """Numeric target that is cryptographically/semantically bound upstream.

    This object is deliberately stricter than a generic pair of arrays.  It is
    not accepted unless the accompanying formal admission says that the actual
    theorem values have been materialized.  Consequently constructing this
    dataclass by hand cannot turn the current opaque theorem witness into
    correction evidence.
    """

    application_id: str
    delta_c_theta: tuple[float, ...]
    delta_c_z: tuple[float, ...]
    epsilon: tuple[float, ...]
    provenance: str
    evidence_kind: str = FORMAL_NUMERIC_EVIDENCE_KIND
    units: str = PHYSICAL_TARGET_UNITS
    surrogate_defect_used: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.application_id, str) or not self.application_id.strip():
            raise ValueError("application_id must be nonempty")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("provenance must be nonempty")
        if self.evidence_kind != FORMAL_NUMERIC_EVIDENCE_KIND:
            raise ValueError("numeric target evidence_kind must be formal-machine-replay")
        if self.units != PHYSICAL_TARGET_UNITS:
            raise ValueError("numeric target units must be physical phase-mean covariance/stress")
        if self.surrogate_defect_used is not False:
            raise ValueError("surrogate mean-defect targets are forbidden")

        theta = _finite_vector(self.delta_c_theta, "delta_c_theta")
        axial = _finite_vector(self.delta_c_z, "delta_c_z")
        epsilon = _finite_vector(self.epsilon, "epsilon")
        if not (len(theta) == len(axial) == len(epsilon)):
            raise ValueError("delta_c_theta, delta_c_z, and epsilon must have equal length")
        if any(value <= 0.0 for value in epsilon):
            raise ValueError("epsilon must be strictly positive")
        object.__setattr__(self, "delta_c_theta", theta)
        object.__setattr__(self, "delta_c_z", axial)
        object.__setattr__(self, "epsilon", epsilon)


@dataclass(frozen=True)
class KokunoActualSignedMeanCorrectionGate:
    """Admit a real signed-mean inverse input, never a theorem-shaped surrogate."""

    def evaluate(
        self,
        admission: ActualSignedMeanDefectAdmission,
        signed_rank_receipt: Mapping[str, Any],
        *,
        numeric_target: ActualSignedMeanNumericTarget | None = None,
    ) -> dict[str, Any]:
        if not isinstance(admission, ActualSignedMeanDefectAdmission):
            raise TypeError("admission must be an ActualSignedMeanDefectAdmission")
        if admission.status != "formal-structure":
            raise ValueError("signed mean-defect admission must retain formal-structure status")
        if admission.requested_cross_defect_identity_admitted is not True:
            raise ValueError("requested_cross_defect theorem identity is not admitted")
        if admission.finite_head_band_admitted is not True:
            raise ValueError("finite-head theorem band is not admitted")

        required_rank_keys = {
            "local_supplied_signed_physical_covariance_rank_two",
            "physical_complete_curl_covariance_rank_two_assessed",
            "actual_positive_order_background_bound",
            "actual_source_h_sigma_pulse_integrals_bound",
            "actual_signed_auxiliary_rectangles_bound",
            "actual_auxiliary_torus_mode_family_bound",
            "source_actual_partition_labels_instantiated",
            "genuinely_independent_second_covariance_column_ready",
            "real_candidate_defect_consumed",
        }
        missing = sorted(required_rank_keys - set(signed_rank_receipt))
        if missing:
            raise ValueError(f"signed covariance rank receipt is missing required keys: {missing}")

        theorem_values_materialized = bool(
            admission.actual_mean_cross_values_materialized
            and admission.missing_weight_values_materialized
            and admission.finite_head_mean_debt_materialized
        )
        actual_source_family_bound = bool(
            signed_rank_receipt["actual_positive_order_background_bound"]
            and signed_rank_receipt["actual_source_h_sigma_pulse_integrals_bound"]
            and signed_rank_receipt["actual_signed_auxiliary_rectangles_bound"]
            and signed_rank_receipt["actual_auxiliary_torus_mode_family_bound"]
            and signed_rank_receipt["source_actual_partition_labels_instantiated"]
        )
        supplied_physical_rank_two = bool(
            signed_rank_receipt["physical_complete_curl_covariance_rank_two_assessed"]
            and signed_rank_receipt["local_supplied_signed_physical_covariance_rank_two"]
        )
        actual_physical_rank_two = bool(
            actual_source_family_bound
            and signed_rank_receipt["genuinely_independent_second_covariance_column_ready"]
        )

        if numeric_target is not None and not theorem_values_materialized:
            raise ValueError(
                "numeric target supplied while the formal admission still marks the actual "
                "mean-cross/missing-weight/debt values opaque"
            )

        reference_target_theta: tuple[float, ...] | None = None
        reference_target_z: tuple[float, ...] | None = None
        real_numeric_target_bound = False
        if numeric_target is not None:
            if not isinstance(numeric_target, ActualSignedMeanNumericTarget):
                raise TypeError("numeric_target must be an ActualSignedMeanNumericTarget")
            if numeric_target.application_id != admission.witness.application_id:
                raise ValueError("numeric target application_id does not match theorem admission")
            if numeric_target.provenance != admission.witness.provenance:
                raise ValueError("numeric target provenance does not match theorem admission")
            reference_target_theta = tuple(
                delta / epsilon
                for delta, epsilon in zip(
                    numeric_target.delta_c_theta, numeric_target.epsilon, strict=True
                )
            )
            reference_target_z = tuple(
                delta / epsilon
                for delta, epsilon in zip(
                    numeric_target.delta_c_z, numeric_target.epsilon, strict=True
                )
            )
            real_numeric_target_bound = True

        signed_mean_inverse_input_ready = bool(
            theorem_values_materialized
            and real_numeric_target_bound
            and actual_source_family_bound
            and actual_physical_rank_two
        )

        blockers: list[str] = []
        if not theorem_values_materialized:
            blockers.append("actual finite-head mean debt values remain opaque")
        if not actual_source_family_bound:
            blockers.append("actual signed complete-curl source family is not bound")
        if not actual_physical_rank_two:
            blockers.append("actual-source physical covariance rank two is not certified")
        if not real_numeric_target_bound:
            blockers.append("no theorem-bound numeric physical mean-defect target is available")

        return {
            "requested_cross_defect_identity_admitted": True,
            "finite_head_band_admitted": True,
            "requested_cross_defect_theorem_machine_replayed": bool(
                admission.requested_cross_defect_theorem_machine_replayed
            ),
            "actual_mean_cross_values_materialized": bool(
                admission.actual_mean_cross_values_materialized
            ),
            "missing_weight_values_materialized": bool(
                admission.missing_weight_values_materialized
            ),
            "finite_head_mean_debt_materialized": bool(
                admission.finite_head_mean_debt_materialized
            ),
            "theorem_values_materialized": theorem_values_materialized,
            "supplied_signed_physical_covariance_rank_two": supplied_physical_rank_two,
            "actual_signed_complete_curl_source_family_bound": actual_source_family_bound,
            "actual_source_physical_covariance_rank_two": actual_physical_rank_two,
            "real_numeric_target_bound": real_numeric_target_bound,
            "physical_to_reference_conversion": "H_ref y = Delta C / epsilon",
            "reference_target_theta": reference_target_theta,
            "reference_target_z": reference_target_z,
            "signed_mean_inverse_input_ready": signed_mean_inverse_input_ready,
            "blockers": tuple(blockers),
            "surrogate_defect_used": False,
            "real_candidate_defect_consumed": signed_mean_inverse_input_ready,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "finite_correction_cycle_run": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "scope": (
                "theorem-to-numeric signed mean-correction readiness only; downstream bounded "
                "inverse, budget, spacetime, radial, damping and held-out NS gates remain required"
            ),
        }
