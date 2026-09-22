"""Stable current-lineage angular-defect correction for Kokuno relative swirl.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
corrected 2026-09-09 reconstruction.

The public reconstruction defines

    H = sqrt(2X) E,
    I(X,eta) = int_0^X H dx,
    r_I = I/(XH),
    r_I' + (1+l) r_I = 1,        l = D_X log H.

On the eta-independent hold l=-lambda.  The two public relative-swirl bumps
then correct the entering angular defect so that I=XH/(1-lambda).

The repository's current leading lineage is not the hidden/source-exact field.
Its autonomous finite-N modulation plus frozen PA.17 repair leaves a small but
nonzero five-moment residual in row order (M,J,I,S,C_p).  Once the compact I1
perturbations end, their *difference* in I relative to the imported ideal outer
profile is a physical constant, because the current and imported profiles have
identical H thereafter.  This module turns that measured PA.17 I residual into
its contribution to the later relative-swirl angular target without forming

    r_I - 1/(1-lambda)

at the late hold.  That subtraction would erase the ~lambda^28-scale defect in
binary64.  If delta_I denotes the physical current-minus-imported angular
moment discrepancy, the first bump row is scaled at its center y1, so the exact
current correction to the row target is simply

    delta_target_current = - delta_I / (X H)|_{y1}.

PA.17 stores

    delta_I = r_I1 * X_1^(3/2) e_1,

where r_I1 is its normalized I-row residual.  The ratio above is evaluated in
logarithms, using the executable #1148 eta-independent hold, and therefore
never constructs enormous X or subtracts two O(1) r_I values.

Important provenance boundary
-----------------------------
PA.17's residual is a *difference from the imported ideal outer profile*, not
an absolute angular cumulative.  The corrected public reconstruction only
bounds the imported profile's pre-bump discrepancy by O(C_pre lambda^28) in the
passage used here; this repository has not yet materialized that absolute base
angular target.  Accordingly this module exposes only the exact current-lineage
*correction* to that target.  It deliberately has no zero default for the
missing base target and does not claim current absolute I/r_I is materialized.
It does not compose the relative-swirl bumps into Cartesian velocity and is not
Navier--Stokes validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_current_pulse_entry_moments import KokunoCurrentPulseEntryMoments
from .kokuno_pa16_current_cartesian_postpulse_eta_independent_hold_prefix import (
    KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix,
)
from .kokuno_public_relative_swirl_compensator import (
    SOURCE_BLOB,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_RELEASE,
    SOURCE_RELEASE_DATE,
    SOURCE_REPOSITORY,
    KokunoPublicRelativeSwirlCompensator,
)

SCHEMA = "kokuno-current-relative-swirl-angular-correction-v1"
PARENT_EXACT_HEAD = "d0e855e37c8f5d6d58f1a85b1dcb874651d20838"
MOMENT_ORDER = ("M", "J", "I", "S", "C_p")

_SOURCE_FORMULAS = {
    "angular_cumulative": "H=sqrt(2X)E; I=int_0^X H dx; r_I=I/(XH)",
    "angular_transport": "D_logX r_I + (1+l)r_I=1, l=D_logX log H",
    "hold_fixed_point": "on l=-lambda, r_*=1/(1-lambda)",
    "relative_swirl_target": "two bumps set I=XH/(1-lambda)",
    "PA17_I_scale": "delta I/(X_1^(3/2)e_1)=r_I1 in row order (M,J,I,S,C_p)",
    "stable_current_correction": "delta_target_current=-delta_I/(XH at first bump center)",
}

_NUMERICAL_REALIZATION = {
    "parent": "stack exactly on A1 #1154 public relative-swirl compensator algebra",
    "I_history": (
        "consume the actual frozen current PA.17 I-row residual already present in the A1 ancestry; "
        "never replace it by zero"
    ),
    "late_hold_stability": (
        "transport the physical delta_I directly into the row-scaled target; do not subtract two "
        "O(1) r_I values after ~30 log(1/lambda) decay"
    ),
    "denominator": (
        "evaluate log(XH) at the public first bump center from the executable #1148 "
        "eta-independent hold"
    ),
    "missing_absolute_anchor_policy": (
        "fail semantically closed: expose only current-minus-imported target correction until the "
        "imported/base absolute angular target is independently materialized"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_PA17_I_residual_consumed": True,
    "current_lineage_angular_target_correction_materialized": True,
    "late_hold_angular_cancellation_avoided": True,
    "imported_base_absolute_angular_target_materialized": False,
    "current_lineage_angular_entry_I_materialized": False,
    "current_lineage_absolute_r_I_materialized": False,
    "current_relative_swirl_total_target_materialized": False,
    "current_cartesian_relative_swirl_composed": False,
    "source_exact_pointwise_bump_shape_recovered": False,
    "source_hidden_parameters_recovered": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _eta_array(value: Any) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError("eta must contain only finite values")
    if np.any(np.abs(out) > 1.0):
        raise ValueError("eta must lie in [-1,1]")
    return out


@dataclass(frozen=True)
class KokunoCurrentRelativeSwirlAngularCorrection:
    """Materialize only the current PA.17 contribution to the late angular target."""

    moments: KokunoCurrentPulseEntryMoments = field(
        default_factory=KokunoCurrentPulseEntryMoments,
        repr=False,
        compare=False,
    )
    hold: KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix = field(
        default_factory=KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix,
        repr=False,
        compare=False,
    )
    compensator: KokunoPublicRelativeSwirlCompensator = field(
        default_factory=KokunoPublicRelativeSwirlCompensator,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.moments, KokunoCurrentPulseEntryMoments):
            raise TypeError("moments must be KokunoCurrentPulseEntryMoments")
        if not isinstance(
            self.hold, KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix
        ):
            raise TypeError(
                "hold must be KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix"
            )
        if not isinstance(self.compensator, KokunoPublicRelativeSwirlCompensator):
            raise TypeError("compensator must be KokunoPublicRelativeSwirlCompensator")
        lam = self.lambda_value
        if self.moments.lambda_value != lam or self.hold.lambda_value != lam:
            raise ValueError("A1 moment, hold and relative-swirl identities use different lambda")
        if not math.isclose(
            self.hold.hold_length,
            self.compensator.hold_length,
            rel_tol=0.0,
            abs_tol=2.0e-13,
        ):
            raise ValueError("#1148 hold and #1154 public compensator use different hold lengths")
        if not math.isclose(
            self.hold.first_relative_swirl_center_s,
            self.compensator.y1,
            rel_tol=0.0,
            abs_tol=2.0e-13,
        ):
            raise ValueError("#1148 and #1154 disagree on the first relative-swirl center")

    @property
    def lambda_value(self) -> float:
        return float(self.compensator.lambda_value)

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_TRUTH_BOUNDARY)

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    @property
    def log_delta_I_scale(self) -> float:
        """log(X_1^(3/2)e_1), the physical PA.17 I-row scale."""

        repair = self.moments.repair
        return float(1.5 * repair.log_X_1 + repair.log_e_1)

    def post_i1_I_residual_with_eta(self, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        """Return the normalized current-minus-imported PA.17 I residual and eta jet."""

        values = _eta_array(eta)
        residual, residual_eta = self.moments.post_i1_residual_normalized_with_eta(values)
        return (
            np.asarray(residual[..., 2], dtype=float),
            np.asarray(residual_eta[..., 2], dtype=float),
        )

    def log_XH_first_bump_center(self, eta: Any) -> np.ndarray:
        """Return log(XH) of the unedited #1148 hold at public first bump center."""

        values = _eta_array(eta)
        endpoint = self.hold._flatten_endpoint_state(values)
        E_flat = np.asarray(endpoint["E"], dtype=float)
        if np.any(E_flat <= 0.0) or np.any(~np.isfinite(E_flat)):
            raise RuntimeError("#1148 flattening endpoint must have finite positive E")
        y1 = self.compensator.y1
        return (
            1.5 * self.hold.log_X_flatten_end
            + 0.5 * math.log(2.0)
            + np.log(E_flat)
            + (1.0 - self.lambda_value) * y1
        )

    def correction_target_with_eta(self, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        """Return current PA.17 correction to the #1154 scaled angular target and eta jet.

        The denominator XH on the eta-independent hold has zero eta derivative at
        the flattening endpoint by construction, so the analytic eta jet is the
        same deterministic scale times the PA.17 residual eta jet.
        """

        values = _eta_array(eta)
        residual_I, residual_I_eta = self.post_i1_I_residual_with_eta(values)
        exponent = self.log_delta_I_scale - self.log_XH_first_bump_center(values)
        # This is a ratio of physical scales, not a free coefficient.  Refuse to
        # silently under/overflow a nonzero current residual.
        if np.any(exponent > math.log(np.finfo(float).max)):
            raise OverflowError("current angular correction scale exceeds float64 range")
        tiny_log = math.log(np.nextafter(0.0, 1.0))
        if np.any((exponent < tiny_log) & (residual_I != 0.0)):
            raise OverflowError(
                "current angular correction scale underflows float64; retain a log-domain target"
            )
        scale = np.exp(exponent)
        target = -scale * residual_I
        target_eta = -scale * residual_I_eta
        if np.any(~np.isfinite(target)) or np.any(~np.isfinite(target_eta)):
            raise RuntimeError("current angular correction target became non-finite")
        return np.asarray(target), np.asarray(target_eta)

    def apply_to_materialized_base_target(
        self,
        eta: Any,
        *,
        base_target: Any,
        base_target_eta: Any,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Add this measured current correction to an explicitly supplied base target.

        There is intentionally no default for either base quantity.  This method
        does not authenticate the external base target and therefore does not
        change this object's truth boundary.
        """

        values = _eta_array(eta)
        base, base_eta, values = np.broadcast_arrays(
            np.asarray(base_target, dtype=float),
            np.asarray(base_target_eta, dtype=float),
            values,
        )
        if np.any(~np.isfinite(base)) or np.any(~np.isfinite(base_eta)):
            raise ValueError("base_target and base_target_eta must be finite")
        correction, correction_eta = self.correction_target_with_eta(values)
        return np.asarray(base + correction), np.asarray(base_eta + correction_eta)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_path": SOURCE_PATH,
            "source_blob": SOURCE_BLOB,
            "source_release": SOURCE_RELEASE,
            "source_release_date": SOURCE_RELEASE_DATE,
            "moment_order": list(MOMENT_ORDER),
            "lambda_value": self.lambda_value,
            "source_formulas": self.source_formulas,
            "numerical_realization": self.numerical_realization,
            "truth_boundary": self.truth_boundary,
        }

    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def save_configuration(self, path: str | Path) -> None:
        Path(path).write_text(_canonical_json(self.configuration()) + "\n", encoding="utf-8")

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoCurrentRelativeSwirlAngularCorrection":
        expected = cls().configuration()
        actual = json.loads(_canonical_json(dict(payload)))
        if actual != expected:
            raise ValueError(
                "current angular-correction configuration/provenance drift detected; "
                "a changed realization requires a new semantic identity"
            )
        return cls()

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoCurrentRelativeSwirlAngularCorrection":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("configuration must decode to a JSON object")
        return cls.from_configuration(payload)
