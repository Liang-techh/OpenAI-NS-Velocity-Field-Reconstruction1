"""Absolute current angular target for Kokuno's late relative-swirl seam.

Pinned provenance is the KokunoYumeto corrected 2026-09-09 reconstruction at
commit 143f6773feb424ad9ed3a8d116653200f20346b7.  The public cumulative is

    H=sqrt(2X)E,  I(X)=int_0^X H dx,  r_I=I/(XH),
    D_logX r_I + (1+l)r_I = 1.

A1 #1159 deliberately exposed only the current-minus-imported PA.17 correction:
subtracting two O(1) late-hold r_I values would erase the O(lambda^28) defect
in binary64.  This child closes that missing base seam without inventing it.

First it computes the *absolute current* r_I at the end of eta flattening.  The
axis-side anchor is the actual A1 PA.10 cumulative I at X_i=110, already
materialized by KokunoPA10ActualXiPrefixMoments.  From log X_i to the flattening
endpoint it integrates the exact current A1 E profile in the stable normalized
identity

 r(b)=r(a) exp(log(XH)_a-log(XH)_b)
      + int_a^b exp(3(s-b)/2) E(s)/E(b) ds.

No enormous X or I is formed.  On the following unedited hold l=-lambda, so the
defect from r_*=1/(1-lambda) transports exactly.  At the first relative-swirl
row center y1, the scaled total target is therefore

    target_current = -(r_flat-r_*) exp(-(1-lambda)y1).

This is numerically stable even when the target is ~lambda^28.  Subtracting
#1159's independently transported current PA.17 contribution at the same
scaled row yields the previously missing imported/base target.  Both quantities
are O(lambda^28), so this subtraction does not suffer the forbidden O(1)-O(1)
cancellation.

The eta jet is a deterministic five-point finite-difference derivative of this
new cumulative adapter.  It is intentionally labelled numerical rather than
source/analytic.  This increment does not compose the relative-swirl bumps into
Cartesian velocity and is not Navier--Stokes validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_current_relative_swirl_angular_correction import (
    KokunoCurrentRelativeSwirlAngularCorrection,
)
from .kokuno_pa10_actual_xi_prefix_moments import KokunoPA10ActualXiPrefixMoments

SCHEMA = "kokuno-current-relative-swirl-absolute-target-v1"
PARENT_EXACT_HEAD = "6ddf35b9dc726c720033ee44161e29efc77d755e"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"

QUADRATURE_ORDER = 8
MAX_LOG_PANEL_WIDTH = 0.25
ETA_DERIVATIVE_STEP = 2.0e-5

_SOURCE_FORMULAS = {
    "angular_cumulative": "H=sqrt(2X)E; I=int_0^X H dx; r_I=I/(XH)",
    "normalized_transport": "D_logX r_I+(1+l)r_I=1",
    "stable_interval_identity": (
        "r(b)=r(a)exp(log(XH)_a-log(XH)_b)"
        "+int_a^b exp(3(s-b)/2)E(s)/E(b) ds"
    ),
    "hold_fixed_point": "l=-lambda => r_*=1/(1-lambda)",
    "scaled_total_target": "target=-(r_flat-r_*)exp(-(1-lambda)y1)",
    "relative_swirl_target": "two bumps set I=XH/(1-lambda)",
}

_NUMERICAL_REALIZATION = {
    "axis_side_anchor": (
        "actual current A1 PA.10 physical I at X_i=110 from "
        "KokunoPA10ActualXiPrefixMoments; no imported/base zero assumption"
    ),
    "outer_transport": (
        f"fixed composite Gauss-Legendre order {QUADRATURE_ORDER}, "
        f"maximum log-X panel width {MAX_LOG_PANEL_WIDTH}"
    ),
    "eta_jet": (
        f"deterministic five-point finite difference with h={ETA_DERIVATIVE_STEP}; "
        "numerical derivative, not source/analytic"
    ),
    "base_target": (
        "stable O(lambda^28)-scale total target minus #1159's independently "
        "transported current PA.17 correction"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "absolute_current_r_I_at_flattening_endpoint_materialized": True,
    "current_relative_swirl_total_target_materialized": True,
    "imported_base_absolute_angular_target_materialized": True,
    "current_PA17_target_correction_consumed": True,
    "late_hold_O1_subtraction_used": False,
    "eta_jet_materialized": True,
    "eta_jet_is_numerical_not_source_analytic": True,
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
class KokunoCurrentRelativeSwirlAbsoluteTarget:
    """Axis-anchored absolute angular target for the current A1 lineage."""

    correction: KokunoCurrentRelativeSwirlAngularCorrection = field(
        default_factory=KokunoCurrentRelativeSwirlAngularCorrection,
        repr=False,
        compare=False,
    )
    prefix: KokunoPA10ActualXiPrefixMoments = field(
        default_factory=KokunoPA10ActualXiPrefixMoments,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.correction, KokunoCurrentRelativeSwirlAngularCorrection):
            raise TypeError("correction must be KokunoCurrentRelativeSwirlAngularCorrection")
        if not isinstance(self.prefix, KokunoPA10ActualXiPrefixMoments):
            raise TypeError("prefix must be KokunoPA10ActualXiPrefixMoments")
        if not self.correction.truth_boundary[
            "current_lineage_angular_target_correction_materialized"
        ]:
            raise ValueError("exact #1159 current angular correction is required")
        lo, hi = self.prefix.eta_interval
        clo, chi = self.correction.hold.eta_interval
        if max(lo, clo) >= min(hi, chi):
            raise ValueError("PA.10 prefix and current outer lineage have disjoint eta intervals")
        if not self.prefix.X_i > 0.0:
            raise ValueError("absolute angular anchor requires X_i>0")
        if not self.log_X_i < self.log_X_flatten_end:
            raise ValueError("PA.10 angular anchor must precede eta-flattening endpoint")

    @property
    def hold(self):
        return self.correction.hold

    @property
    def compensator(self):
        return self.correction.compensator

    @property
    def lambda_value(self) -> float:
        return float(self.correction.lambda_value)

    @property
    def log_X_i(self) -> float:
        return math.log(float(self.prefix.X_i))

    @property
    def log_X_flatten_end(self) -> float:
        return float(self.hold.log_X_flatten_end)

    @property
    def r_star(self) -> float:
        return 1.0 / (1.0 - self.lambda_value)

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_TRUTH_BOUNDARY)

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    def _prefix_anchor(self, eta: float) -> tuple[float, float]:
        """Return absolute r_I(X_i) and E(X_i) from the actual current prefix."""
        moments = np.asarray(
            self.prefix.physical_prefix_moments_at_Xi(float(eta)), dtype=float
        )
        if moments.shape != (5,):
            raise RuntimeError("PA.10 physical prefix moments changed shape")
        # PA.10 public order is (M,I,J,S,C_p).
        I = float(moments[1])
        vals = self.prefix.profile_values(self.prefix.X_i, float(eta))
        E = float(np.asarray(vals["E_actual_prefix"]))
        H = math.sqrt(2.0 * self.prefix.X_i) * E
        denom = self.prefix.X_i * H
        if not (math.isfinite(I) and math.isfinite(E) and E > 0.0 and denom > 0.0):
            raise RuntimeError("invalid absolute PA.10 angular anchor")
        return I / denom, E

    def _outer_E(self, log_X: np.ndarray, eta: float) -> np.ndarray:
        """Evaluate exact current A1 E before/through eta flattening."""
        xx = np.asarray(log_X, dtype=float)
        ee = np.full(xx.shape, float(eta), dtype=float)
        vals = self.hold.parent.similarity_profile_values_logX(xx, ee)
        E = np.asarray(
            vals["E_current_leading_postpulse_eta_flattening"], dtype=float
        )
        if E.shape != xx.shape or np.any(~np.isfinite(E)) or np.any(E <= 0.0):
            raise RuntimeError("current A1 outer E evaluator returned invalid values")
        return E

    def _absolute_r_flat_scalar(self, eta: float) -> float:
        r0, E0 = self._prefix_anchor(float(eta))
        a = self.log_X_i
        b = self.log_X_flatten_end
        E_b = float(self._outer_E(np.asarray([b]), float(eta))[0])

        # Exact normalized anchor carry: XH = sqrt(2) X^(3/2) E.
        log_ratio = 1.5 * (a - b) + math.log(E0 / E_b)
        anchor = r0 * math.exp(log_ratio)

        span = b - a
        panels = max(1, int(math.ceil(span / MAX_LOG_PANEL_WIDTH)))
        edges = np.linspace(a, b, panels + 1, dtype=float)
        nodes, weights = np.polynomial.legendre.leggauss(QUADRATURE_ORDER)
        half = 0.5 * (edges[1:] - edges[:-1])
        mid = 0.5 * (edges[1:] + edges[:-1])
        s = (mid[:, None] + half[:, None] * nodes[None, :]).reshape(-1)
        E_s = self._outer_E(s, float(eta))
        kernel = np.exp(1.5 * (s - b)) * (E_s / E_b)
        integral = float(
            np.sum(
                half
                * np.sum(
                    kernel.reshape(panels, QUADRATURE_ORDER)
                    * weights[None, :],
                    axis=1,
                )
            )
        )
        result = anchor + integral
        if not math.isfinite(result) or result <= 0.0:
            raise RuntimeError("absolute normalized angular cumulative became invalid")
        return result

    def absolute_r_flat(self, eta: Any) -> np.ndarray:
        """Vectorized absolute current r_I at the eta-flattening endpoint."""
        values = _eta_array(eta)
        flat = values.reshape(-1)
        out = np.asarray(
            [self._absolute_r_flat_scalar(float(v)) for v in flat], dtype=float
        )
        return out.reshape(values.shape)

    def total_target(self, eta: Any) -> np.ndarray:
        """Stable current total scaled angular row target at the first bump."""
        values = _eta_array(eta)
        r_flat = self.absolute_r_flat(values)
        y1 = float(self.compensator.y1)
        scale = math.exp(-(1.0 - self.lambda_value) * y1)
        target = -(r_flat - self.r_star) * scale
        if np.any(~np.isfinite(target)):
            raise RuntimeError("current relative-swirl total target became non-finite")
        return np.asarray(target)

    def _derivative_scalar(self, eta: float) -> float:
        h = ETA_DERIVATIVE_STEP
        lo, hi = self.prefix.eta_interval
        x = float(eta)
        f = lambda z: float(self.total_target(float(z)))
        if x - 2.0 * h >= lo and x + 2.0 * h <= hi:
            return (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)) / (12*h)
        if x + 4.0 * h <= hi:
            return (-25*f(x) + 48*f(x+h) - 36*f(x+2*h) + 16*f(x+3*h) - 3*f(x+4*h)) / (12*h)
        if x - 4.0 * h >= lo:
            return (25*f(x) - 48*f(x-h) + 36*f(x-2*h) - 16*f(x-3*h) + 3*f(x-4*h)) / (12*h)
        raise ValueError("eta interval is too narrow for the frozen derivative stencil")

    def target_with_eta(self, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        """Return total target plus deterministic numerical eta derivative."""
        values = _eta_array(eta)
        target = self.total_target(values)
        derivative = np.asarray(
            [self._derivative_scalar(float(v)) for v in values.reshape(-1)],
            dtype=float,
        ).reshape(values.shape)
        if np.any(~np.isfinite(derivative)):
            raise RuntimeError("angular target eta derivative became non-finite")
        return target, derivative

    def decomposition_with_eta(self, eta: Any) -> dict[str, np.ndarray]:
        """Return total=current base + measured #1159 PA.17 correction."""
        values = _eta_array(eta)
        total, total_eta = self.target_with_eta(values)
        correction, correction_eta = self.correction.correction_target_with_eta(values)
        base = total - correction
        base_eta = total_eta - correction_eta
        arrays = (total, total_eta, correction, correction_eta, base, base_eta)
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("relative-swirl angular target decomposition is non-finite")
        return {
            "eta": values,
            "absolute_r_I_flatten_end": self.absolute_r_flat(values),
            "total_target_current": np.asarray(total),
            "total_target_current_eta": np.asarray(total_eta),
            "current_PA17_correction_target": np.asarray(correction),
            "current_PA17_correction_target_eta": np.asarray(correction_eta),
            "imported_base_target": np.asarray(base),
            "imported_base_target_eta": np.asarray(base_eta),
        }

    def solve_current_target(self, eta: Any):
        """Feed the materialized total target/jet into the existing #1154 solver."""
        total, total_eta = self.target_with_eta(eta)
        return self.compensator.solve_target_jet(total, total_eta)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_blob": SOURCE_BLOB,
            "source_path": SOURCE_PATH,
            "source_release": SOURCE_RELEASE,
            "source_release_date": SOURCE_RELEASE_DATE,
            "quadrature_order": QUADRATURE_ORDER,
            "max_log_panel_width": MAX_LOG_PANEL_WIDTH,
            "eta_derivative_step": ETA_DERIVATIVE_STEP,
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
    ) -> "KokunoCurrentRelativeSwirlAbsoluteTarget":
        expected = cls().configuration()
        actual = json.loads(_canonical_json(dict(payload)))
        if actual != expected:
            raise ValueError(
                "absolute angular-target configuration/provenance drift detected; "
                "a changed realization requires a new semantic identity"
            )
        return cls()

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoCurrentRelativeSwirlAbsoluteTarget":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("configuration must decode to a JSON object")
        return cls.from_configuration(payload)
