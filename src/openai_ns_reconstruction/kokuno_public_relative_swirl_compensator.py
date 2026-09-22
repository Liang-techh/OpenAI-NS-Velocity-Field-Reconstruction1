"""Executable low-dimensional relative-swirl compensator from the corrected Kokuno reconstruction.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

After eta flattening, the public construction retains the eta-independent
exponential for

    T_hold = 30 log(1/lambda).

Near the end it makes a *relative* edit

    E -> E (1 + c1 beta1 + c2 beta2),

using two fixed C-infinity bumps of width .3 centered three and one log-X
units before the end.  The coefficients meet two exact endpoint conditions:

1. angular: I = X H / (1-lambda), where I_X = H;
2. pressure neutrality: integral (E^2 - E_unedited^2) dy = 0.

Writing the unedited exponential relative to the first bump center y1 gives

    a_j = int exp((1-lambda)(y-y1)) beta_j(y) dy,

for the angular row, while exact pressure neutrality becomes

    2 b.c + c^T Q c = 0,

with

    b_j  = int exp((-1-2lambda)(y-y1)) beta_j(y) dy,
    Q_jk = int exp((-1-2lambda)(y-y1)) beta_j(y) beta_k(y) dy.

This module solves those two equations for an explicit angular target.  It
also provides a helper that converts an entering r_I=I/(XH) into that target
under the public unedited ODE r_I' + (1-lambda) r_I = 1.

Important provenance boundary
-----------------------------
The corrected public interface specifies the two smooth nonnegative bumps,
their common width and centers, but does not uniquely expose their pointwise
shape.  We therefore freeze the same repository-autonomous standard C-infinity
bump family already used elsewhere in this reconstruction.  It is NOT claimed
to be source-exact hidden data.

This module deliberately does not manufacture the current-lineage entering I.
It is executable algebra waiting for a separately materialized current r_I/I
state.  It is neither a Cartesian relative-swirl composition nor PDE
validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np


SCHEMA = "kokuno-public-relative-swirl-compensator-v1"
PARENT_EXACT_HEAD = "5bafa4199ef127bf50bce30a58cc69ded99e1aa9"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_RELEASE = "corrected 208-page reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"

CURRENT_LAMBDA = 0.05
PUBLIC_HOLD_MULTIPLIER = 30.0
PUBLIC_BUMP_WIDTH = 0.3
PUBLIC_FIRST_CENTER_FROM_END = 3.0
PUBLIC_SECOND_CENTER_FROM_END = 1.0
CURRENT_PARENT_PREFIX_GUARD_FROM_END = 3.5
QUADRATURE_ORDER = 96

_SOURCE_FORMULAS = {
    "hold_length": "T_hold=30*log(1/lambda)",
    "relative_edit": "E -> E*(1+c1*beta1+c2*beta2)",
    "bump_geometry": "width=.3; centers T_hold-3 and T_hold-1",
    "angular_target": "I=XH/(1-lambda) at the end",
    "angular_ode": "r_I'+(1+l)r_I=1; on the unedited hold l=-lambda",
    "angular_linear_weight": "exp((1-lambda)*y)",
    "pressure_target": "integral(E^2-E_unedited^2)dy=0",
    "pressure_linear_weight": "2*exp((-1-2lambda)*y)",
    "pressure_exact_remainder": "exp((-1-2lambda)*y)*(c1*beta1+c2*beta2)^2",
}

_NUMERICAL_REALIZATION = {
    "bump_shape": (
        "repository-autonomous standard C-infinity bump "
        "exp(1-1/(1-z^2)) for |z|<1, z=2(y-center)/0.3"
    ),
    "quadrature": "fixed 96-point Gauss-Legendre on exact compact supports",
    "row_scaling": "both public moment rows evaluated relative to first bump center",
    "quadratic_root": "real pressure root closest to the zero-target linearized branch",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_relative_swirl_two_bump_algebra_materialized": True,
    "public_hold_length_centers_width_materialized": True,
    "public_angular_endpoint_row_materialized": True,
    "public_pressure_neutrality_quadratic_row_materialized": True,
    "repository_autonomous_pointwise_bump_shape_materialized": True,
    "source_exact_pointwise_bump_shape_recovered": False,
    "current_lineage_angular_entry_I_materialized": False,
    "current_cartesian_relative_swirl_composed": False,
    "complete_terminal_hold_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _standard_bump(y: Any, center: float, width: float) -> np.ndarray:
    yy = _finite(y, "y")
    z = 2.0 * (yy - float(center)) / float(width)
    out = np.zeros_like(z, dtype=float)
    inside = np.abs(z) < 1.0
    if np.any(inside):
        q = z[inside]
        out[inside] = np.exp(1.0 - 1.0 / (1.0 - q * q))
    return out


def _standard_bump_prime(y: Any, center: float, width: float) -> np.ndarray:
    yy = _finite(y, "y")
    z = 2.0 * (yy - float(center)) / float(width)
    out = np.zeros_like(z, dtype=float)
    inside = np.abs(z) < 1.0
    if np.any(inside):
        q = z[inside]
        b = np.exp(1.0 - 1.0 / (1.0 - q * q))
        out[inside] = b * (-4.0 * q / (width * (1.0 - q * q) ** 2))
    return out


@dataclass(frozen=True)
class RelativeSwirlCompensatorSolution:
    """Vectorized exact-row solution plus analytic target derivative."""

    c1: np.ndarray
    c2: np.ndarray
    dc1_dtarget: np.ndarray
    dc2_dtarget: np.ndarray
    angular_residual: np.ndarray
    pressure_residual: np.ndarray
    discriminant: np.ndarray


class KokunoPublicRelativeSwirlCompensator:
    """Public two-bump relative-swirl algebra with explicit angular input."""

    def __init__(self, lambda_value: float = CURRENT_LAMBDA) -> None:
        lam = float(lambda_value)
        if not math.isfinite(lam) or lam <= 0.0 or lam >= 0.25:
            raise ValueError("lambda_value must be finite and satisfy 0<lambda<0.25")
        if lam != CURRENT_LAMBDA:
            raise ValueError(
                f"this current-lineage component freezes lambda={CURRENT_LAMBDA}; "
                "a different lambda requires a new semantic identity"
            )
        self.lambda_value = lam
        self.hold_length = PUBLIC_HOLD_MULTIPLIER * math.log(1.0 / lam)
        self.y1 = self.hold_length - PUBLIC_FIRST_CENTER_FROM_END
        self.y2 = self.hold_length - PUBLIC_SECOND_CENTER_FROM_END
        self.width = PUBLIC_BUMP_WIDTH
        self.angular_slope = 1.0 - lam
        self.pressure_slope = -1.0 - 2.0 * lam
        self._angular_row, self._pressure_linear_row, self._pressure_quadratic = (
            self._build_moments()
        )
        linear = np.vstack((self._angular_row, 2.0 * self._pressure_linear_row))
        if not np.all(np.isfinite(linear)) or abs(float(np.linalg.det(linear))) <= 1.0e-12:
            raise RuntimeError("relative-swirl linearized two-row system is singular")

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
    def angular_row_scaled(self) -> np.ndarray:
        return self._angular_row.copy()

    @property
    def pressure_linear_row_scaled(self) -> np.ndarray:
        return self._pressure_linear_row.copy()

    @property
    def pressure_quadratic_matrix_scaled(self) -> np.ndarray:
        return self._pressure_quadratic.copy()

    @property
    def linearized_condition_number(self) -> float:
        return float(
            np.linalg.cond(
                np.vstack((self._angular_row, 2.0 * self._pressure_linear_row))
            )
        )

    def beta1(self, y: Any) -> np.ndarray:
        return _standard_bump(y, self.y1, self.width)

    def beta2(self, y: Any) -> np.ndarray:
        return _standard_bump(y, self.y2, self.width)

    def beta1_prime(self, y: Any) -> np.ndarray:
        return _standard_bump_prime(y, self.y1, self.width)

    def beta2_prime(self, y: Any) -> np.ndarray:
        return _standard_bump_prime(y, self.y2, self.width)

    def _integral(self, func, a: float, b: float) -> float:
        if not (math.isfinite(a) and math.isfinite(b) and b >= a):
            raise ValueError("quadrature bounds must be finite and ordered")
        if a == b:
            return 0.0
        nodes, weights = np.polynomial.legendre.leggauss(QUADRATURE_ORDER)
        half = 0.5 * (b - a)
        mid = 0.5 * (b + a)
        x = half * nodes + mid
        return float(half * np.sum(weights * np.asarray(func(x), dtype=float)))

    def _build_moments(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        centers = (self.y1, self.y2)
        half_width = 0.5 * self.width
        angular = np.empty(2, dtype=float)
        pressure = np.empty(2, dtype=float)

        for j, center in enumerate(centers):
            lo, hi = center - half_width, center + half_width
            angular[j] = self._integral(
                lambda y, c=center: (
                    np.exp(self.angular_slope * (y - self.y1))
                    * _standard_bump(y, c, self.width)
                ),
                lo,
                hi,
            )
            pressure[j] = self._integral(
                lambda y, c=center: (
                    np.exp(self.pressure_slope * (y - self.y1))
                    * _standard_bump(y, c, self.width)
                ),
                lo,
                hi,
            )

        quadratic = np.zeros((2, 2), dtype=float)
        for i, ci in enumerate(centers):
            for j, cj in enumerate(centers):
                lo = max(ci - half_width, cj - half_width)
                hi = min(ci + half_width, cj + half_width)
                if hi <= lo:
                    continue
                quadratic[i, j] = self._integral(
                    lambda y, a=ci, b=cj: (
                        np.exp(self.pressure_slope * (y - self.y1))
                        * _standard_bump(y, a, self.width)
                        * _standard_bump(y, b, self.width)
                    ),
                    lo,
                    hi,
                )
        return angular, pressure, quadratic

    def angular_target_from_rI(
        self,
        r_I_entry: Any,
        *,
        y_entry: float | None = None,
    ) -> np.ndarray:
        """Convert an entering r_I=I/(XH) into the scaled angular bump target.

        The current A1 Cartesian prefix ends 3.5 log-X units before hold end.
        That value is used only as the default *location* if a future caller
        supplies the separately materialized current r_I there.  This helper
        never invents r_I itself.

        On the unedited eta-independent hold l=-lambda, hence

            r' + (1-lambda) r = 1.

        The bump row is scaled at y1, so a.c must equal

            exp((1-lambda)(T_hold-y1)) * (r_target-r_unedited_end).
        """
        r0 = _finite(r_I_entry, "r_I_entry")
        if y_entry is None:
            y0 = self.hold_length - CURRENT_PARENT_PREFIX_GUARD_FROM_END
        else:
            y0 = float(y_entry)
        first_support = self.y1 - 0.5 * self.width
        if not math.isfinite(y0) or y0 > first_support:
            raise ValueError("y_entry must be finite and no later than first bump support")
        r_star = 1.0 / self.angular_slope
        gap = self.hold_length - y0
        r_unedited_end = r_star + (r0 - r_star) * math.exp(-self.angular_slope * gap)
        return math.exp(
            self.angular_slope * (self.hold_length - self.y1)
        ) * (r_star - r_unedited_end)

    def _pressure_residual_for(self, coeff: np.ndarray) -> float:
        return float(
            2.0 * np.dot(self._pressure_linear_row, coeff)
            + coeff @ self._pressure_quadratic @ coeff
        )

    def _solve_scalar(self, target: float) -> tuple[np.ndarray, np.ndarray, float]:
        a = self._angular_row
        b = self._pressure_linear_row
        qmat = self._pressure_quadratic
        p = np.array([target / a[0], 0.0], dtype=float)
        q = np.array([-a[1] / a[0], 1.0], dtype=float)
        qa = float(q @ qmat @ q)
        qb = float(2.0 * np.dot(b, q) + 2.0 * p @ qmat @ q)
        qc = float(2.0 * np.dot(b, p) + p @ qmat @ p)
        disc = qb * qb - 4.0 * qa * qc

        scale = max(qb * qb, abs(4.0 * qa * qc), 1.0)
        if disc < -1.0e-14 * scale:
            raise ValueError(
                "angular target has no real pressure-neutral relative-swirl solution "
                "for this frozen autonomous bump realization"
            )
        disc = max(0.0, disc)

        linearized = np.linalg.solve(
            np.vstack((a, 2.0 * b)),
            np.array([target, 0.0], dtype=float),
        )
        candidates: list[np.ndarray] = []
        if abs(qa) <= 1.0e-14:
            if abs(qb) <= 1.0e-14:
                raise RuntimeError("degenerate pressure equation on angular line")
            candidates.append(p + q * (-qc / qb))
        else:
            root = math.sqrt(disc)
            candidates.append(p + q * ((-qb + root) / (2.0 * qa)))
            candidates.append(p + q * ((-qb - root) / (2.0 * qa)))

        coeff = min(candidates, key=lambda c: float(np.linalg.norm(c - linearized)))
        if float(np.sum(np.abs(coeff))) >= 1.0:
            raise ValueError(
                "solution leaves the source-small relative-edit regime; "
                "refusing to assert E>0 for this target"
            )

        pressure_gradient = 2.0 * b + 2.0 * (qmat @ coeff)
        jac = np.vstack((a, pressure_gradient))
        det = float(np.linalg.det(jac))
        if not math.isfinite(det) or abs(det) <= 1.0e-13:
            raise ValueError("relative-swirl solution branch has singular target derivative")
        dcoeff = np.linalg.solve(jac, np.array([1.0, 0.0], dtype=float))
        return coeff, dcoeff, disc

    def solve(self, angular_target: Any) -> RelativeSwirlCompensatorSolution:
        """Solve exact angular + pressure rows, vectorized over angular target."""
        target = _finite(angular_target, "angular_target")
        flat = target.reshape(-1)
        coeff = np.empty((2, flat.size), dtype=float)
        dcoeff = np.empty((2, flat.size), dtype=float)
        disc = np.empty(flat.size, dtype=float)
        angular_residual = np.empty(flat.size, dtype=float)
        pressure_residual = np.empty(flat.size, dtype=float)

        for k, value in enumerate(flat):
            c, dc, d = self._solve_scalar(float(value))
            coeff[:, k] = c
            dcoeff[:, k] = dc
            disc[k] = d
            angular_residual[k] = float(np.dot(self._angular_row, c) - value)
            pressure_residual[k] = self._pressure_residual_for(c)

        shape = target.shape
        return RelativeSwirlCompensatorSolution(
            c1=coeff[0].reshape(shape),
            c2=coeff[1].reshape(shape),
            dc1_dtarget=dcoeff[0].reshape(shape),
            dc2_dtarget=dcoeff[1].reshape(shape),
            angular_residual=angular_residual.reshape(shape),
            pressure_residual=pressure_residual.reshape(shape),
            discriminant=disc.reshape(shape),
        )

    def solve_target_jet(
        self,
        angular_target: Any,
        angular_target_eta: Any,
    ) -> tuple[RelativeSwirlCompensatorSolution, np.ndarray, np.ndarray]:
        """Return solution plus analytic eta jet via dc/dtarget * target_eta."""
        target, target_eta = np.broadcast_arrays(
            _finite(angular_target, "angular_target"),
            _finite(angular_target_eta, "angular_target_eta"),
        )
        solution = self.solve(target)
        return (
            solution,
            np.asarray(solution.dc1_dtarget * target_eta),
            np.asarray(solution.dc2_dtarget * target_eta),
        )

    def relative_edit(self, y: Any, c1: Any, c2: Any) -> np.ndarray:
        """Evaluate B=c1*beta1+c2*beta2 with NumPy broadcasting."""
        yy, a1, a2 = np.broadcast_arrays(
            _finite(y, "y"), _finite(c1, "c1"), _finite(c2, "c2")
        )
        return a1 * self.beta1(yy) + a2 * self.beta2(yy)

    def relative_edit_y(self, y: Any, c1: Any, c2: Any) -> np.ndarray:
        """Analytic y derivative of the relative edit."""
        yy, a1, a2 = np.broadcast_arrays(
            _finite(y, "y"), _finite(c1, "c1"), _finite(c2, "c2")
        )
        return a1 * self.beta1_prime(yy) + a2 * self.beta2_prime(yy)

    def edited_E(
        self,
        E_unedited: Any,
        y: Any,
        c1: Any,
        c2: Any,
    ) -> np.ndarray:
        """Apply the source-form relative edit to an explicit unedited E."""
        eu, yy, a1, a2 = np.broadcast_arrays(
            _finite(E_unedited, "E_unedited"),
            _finite(y, "y"),
            _finite(c1, "c1"),
            _finite(c2, "c2"),
        )
        factor = 1.0 + self.relative_edit(yy, a1, a2)
        if np.any(factor <= 0.0):
            raise ValueError("relative-swirl edit must preserve E>0")
        return eu * factor

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
            "lambda_value": self.lambda_value,
            "hold_multiplier": PUBLIC_HOLD_MULTIPLIER,
            "bump_width": self.width,
            "first_center_from_end": PUBLIC_FIRST_CENTER_FROM_END,
            "second_center_from_end": PUBLIC_SECOND_CENTER_FROM_END,
            "quadrature_order": QUADRATURE_ORDER,
            "source_formulas": self.source_formulas,
            "numerical_realization": self.numerical_realization,
            "truth_boundary": self.truth_boundary,
        }

    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def save_configuration(self, path: str | Path) -> None:
        Path(path).write_text(_canonical_json(self.configuration()) + "\n", encoding="utf-8")

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoPublicRelativeSwirlCompensator":
        expected = cls().configuration()
        actual = json.loads(_canonical_json(dict(payload)))
        if actual != expected:
            raise ValueError(
                "relative-swirl configuration/provenance drift detected; "
                "a changed realization requires a new semantic identity"
            )
        return cls(lambda_value=float(actual["lambda_value"]))

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPublicRelativeSwirlCompensator":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("configuration must decode to a JSON object")
        return cls.from_configuration(payload)
