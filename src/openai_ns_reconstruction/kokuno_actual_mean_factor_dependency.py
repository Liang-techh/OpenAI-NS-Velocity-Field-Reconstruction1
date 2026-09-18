"""Pin the exact finite-head signed mean-defect dependency chain.

Agent 3 already has a fail-closed theorem-to-numeric gate in
:mod:`kokuno_actual_signed_mean_correction_gate`.  The pinned formal theorem
states

    meanBar(actualCross) - requestedStress
      = -missingWeight(N, physicalScale(n, x)) * requestedStress.

A tempting but incorrect next step is to invent a convenient smooth dyadic
partition in Python and call the resulting scalar the theorem's
``missingWeight``.  The source does not justify that promotion.  At the pinned
formal commit

* ``missingWeight`` is the finite sum
  ``sum_{m=-1}^{N-1} dyadicMask(m, q)^2``;
* ``physicalScale(n,x)`` is ``Q_n * coordinateQ(2h, x_slow)`` with
  ``Q_n = 2^{-n}``;
* ``requestedStress`` is built from the *actual cycle state* through the theta
  and axial residuals and the compact physical ``barSigma`` primitive;
* ``dyadicMask`` ultimately uses ``SquaredPartition.bump``.  In the pinned
  Mathlib revision, ``ContDiffBump.toFun`` evaluates through
  ``someContDiffBumpBase = Nonempty.some hb.out``.  No numeric bump/profile
  realization is exported by the formal repository into this Python project.

This module turns that provenance boundary into executable bookkeeping.  It
also exposes the parts that are genuinely numeric without choosing a surrogate:
the exact finite set of missing-weight indices and the algebraic physical-scale
map once ``coordinateQ`` is supplied.  It deliberately does *not* evaluate the
bump-dependent dyadic masks and therefore cannot materialize the finite-head
mean debt by itself.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral
from typing import Any

from .actual_signed_common_curl import PINNED_LEAN_COMMIT, PINNED_LEAN_REPOSITORY
from .actual_signed_mean_defect import (
    ActualSignedMeanDefectAdmission,
    PINNED_COVARIANCE_FILE,
    PINNED_REQUESTED_CROSS_DEFECT_THEOREM,
)
from .actual_signed_mean_tail import PINNED_MEAN_FILE


PINNED_PARTITION_FILE = "NavierStokes/SquaredPartition.lean"
PINNED_REQUEST_FILE = "NavierStokes/LocalSignedRequest.lean"
PINNED_INITIALIZATION_FILE = "NavierStokes/CorrectionInitialization.lean"
PINNED_MANIFEST_FILE = "lake-manifest.json"

PINNED_MATHLIB_REPOSITORY = "leanprover-community/mathlib4"
PINNED_MATHLIB_COMMIT = "85e3a25e006c35636f0e53b0e9296caca2685bc0"
PINNED_MATHLIB_BUMP_FILE = "Mathlib/Analysis/Calculus/BumpFunction/Basic.lean"

MISSING_WEIGHT_DEFINITION = "missingWeight(N,q)=sum_{m=-1}^{N-1} dyadicMask(m,q)^2"
PHYSICAL_SCALE_DEFINITION = "physicalScale(n,x)=2^(-n)*coordinateQ(2h,x_slow)"
REQUESTED_STRESS_DEFINITION = (
    "requestedStress=(physicalBarSigma(thetaResidual,exponent=2),"
    " physicalBarSigma(axialResidual,exponent=1))"
)
BUMP_NUMERIC_DEPENDENCY = (
    "dyadicMask->dyadicProfile->lineMask->translatedBump->SquaredPartition.bump"
    "->ContDiffBump.toFun->someContDiffBumpBase->Nonempty.some"
)


def _natural(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(value)


def missing_weight_indices(prepared_n: int) -> tuple[int, ...]:
    """Return the exact Lean ``Finset.Icc (-1) (N-1)`` indices.

    This is source-exact integer bookkeeping and does not evaluate any bump.
    For ``N=0`` the interval correctly contains only ``-1``.
    """

    n = _natural(prepared_n, "prepared_n")
    return tuple(range(-1, n))


def dyadic_q(band: int) -> float:
    """Floating realization of the exact algebraic scale ``Q_n = 2^{-n}``.

    The function fails rather than silently returning zero if IEEE-754
    underflow would erase the positive source scale.
    """

    n = _natural(band, "band")
    value = math.ldexp(1.0, -n)
    if not math.isfinite(value) or value <= 0.0:
        raise OverflowError("band is too large for a positive finite float dyadic scale")
    return value


def physical_scale_from_coordinate_q(band: int, coordinate_q: float) -> float:
    """Evaluate the source algebraic ``physicalScale`` once ``coordinateQ`` is known.

    ``coordinate_q`` must come from the same theorem point.  Supplying it does
    not materialize the bump-dependent ``missingWeight`` or the state-dependent
    requested stress.
    """

    q = float(coordinate_q)
    if not math.isfinite(q) or q <= 0.0:
        raise ValueError("coordinate_q must be strictly positive and finite")
    value = dyadic_q(band) * q
    if not math.isfinite(value) or value <= 0.0:
        raise OverflowError("physical scale is not representable as a positive finite float")
    return value


@dataclass(frozen=True)
class KokunoActualMeanFactorDependencyGate:
    """Expose what the pinned defect theorem does and does not numerically fix."""

    def evaluate(self, admission: ActualSignedMeanDefectAdmission) -> dict[str, Any]:
        if not isinstance(admission, ActualSignedMeanDefectAdmission):
            raise TypeError("admission must be an ActualSignedMeanDefectAdmission")
        if admission.status != "formal-structure":
            raise ValueError("mean-defect admission must retain formal-structure status")
        if admission.witness.lean_repository != PINNED_LEAN_REPOSITORY:
            raise ValueError("mean-defect repository drifted from the pinned formal source")
        if admission.witness.lean_commit != PINNED_LEAN_COMMIT:
            raise ValueError("mean-defect commit drifted from the pinned formal source")
        if admission.witness.theorem_symbol != PINNED_REQUESTED_CROSS_DEFECT_THEOREM:
            raise ValueError("mean-defect theorem symbol drifted from requested_cross_defect")
        if admission.witness.lean_file != PINNED_MEAN_FILE:
            raise ValueError("mean-defect theorem file drifted from the pinned source")

        prepared_n = admission.witness.prepared_N
        band = admission.witness.band
        component = admission.witness.component
        indices = missing_weight_indices(prepared_n)

        source_algebra_pinned = bool(
            admission.requested_cross_defect_identity_admitted
            and admission.finite_head_band_admitted
        )
        requested_stress_numeric_ready = bool(
            admission.actual_mean_cross_values_materialized
            and admission.finite_head_mean_debt_materialized
        )
        missing_weight_numeric_ready = bool(admission.missing_weight_values_materialized)
        finite_head_debt_numeric_ready = bool(
            source_algebra_pinned
            and requested_stress_numeric_ready
            and missing_weight_numeric_ready
        )

        blockers: list[str] = []
        if not missing_weight_numeric_ready:
            blockers.append(
                "no machine-linked numeric realization of the pinned dyadic bump/missingWeight is exported"
            )
        if not requested_stress_numeric_ready:
            blockers.append(
                "actual cycle-state theta/axial residual requestedStress values remain opaque"
            )
        if not finite_head_debt_numeric_ready:
            blockers.append("finite-head physical mean debt is not numerically materialized")

        return {
            "formal_repository": PINNED_LEAN_REPOSITORY,
            "formal_commit": PINNED_LEAN_COMMIT,
            "mean_file": PINNED_MEAN_FILE,
            "covariance_file": PINNED_COVARIANCE_FILE,
            "partition_file": PINNED_PARTITION_FILE,
            "request_file": PINNED_REQUEST_FILE,
            "initialization_file": PINNED_INITIALIZATION_FILE,
            "manifest_file": PINNED_MANIFEST_FILE,
            "mathlib_repository": PINNED_MATHLIB_REPOSITORY,
            "mathlib_commit": PINNED_MATHLIB_COMMIT,
            "mathlib_bump_file": PINNED_MATHLIB_BUMP_FILE,
            "requested_cross_defect_theorem": PINNED_REQUESTED_CROSS_DEFECT_THEOREM,
            "defect_formula": "DeltaC=-missingWeight(N,physicalScale(n,x))*requestedStress",
            "missing_weight_definition": MISSING_WEIGHT_DEFINITION,
            "physical_scale_definition": PHYSICAL_SCALE_DEFINITION,
            "requested_stress_definition": REQUESTED_STRESS_DEFINITION,
            "bump_numeric_dependency": BUMP_NUMERIC_DEPENDENCY,
            "prepared_N": prepared_n,
            "band": band,
            "component": "theta" if component == 0 else "axial",
            "missing_weight_indices": indices,
            "missing_weight_term_count": len(indices),
            "source_algebra_pinned": source_algebra_pinned,
            "physical_scale_algebra_executable_given_coordinate_q": True,
            "missing_weight_depends_on_pinned_mathlib_contdiff_bump": True,
            "mathlib_bump_base_selected_by_nonempty_some": True,
            "numeric_bump_profile_export_present": False,
            "missing_weight_values_materialized": missing_weight_numeric_ready,
            "requested_stress_actual_state_values_materialized": requested_stress_numeric_ready,
            "finite_head_mean_debt_materialized": finite_head_debt_numeric_ready,
            "arbitrary_python_bump_may_be_promoted_to_actual": False,
            "surrogate_defect_used": False,
            "real_candidate_defect_consumed": False,
            "signed_mean_inverse_input_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_rerun_allowed": False,
            "heldout_ns_residual_assessed": False,
            "residual_reduction_claimed": False,
            "pde_validated": False,
            "blockers": tuple(blockers),
            "scope": (
                "exact finite-head mean-factor provenance and partial algebraic evaluation only; "
                "no dyadic bump, actual cycle state, correction velocity, or NS residual is materialized"
            ),
        }
