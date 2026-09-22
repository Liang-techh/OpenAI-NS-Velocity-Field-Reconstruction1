"""Executable public M/J end-compensator algebra for Kokuno's axial pulse.

Pinned provenance
-----------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
corrected 2026-09-09 reconstruction.

For y=log(X/X_p), the public reconstruction writes

    R_b = Amp(eta) R_0(lambda y) + c_1(eta) beta_1(y) + c_2(eta) beta_2(y)

and chooses the two end coefficients by the exact two-row system

    [int exp(s1 y) beta_1, int exp(s1 y) beta_2] [c1] = -[m_p + Amp I1]
    [int exp(s2 y) beta_1, int exp(s2 y) beta_2] [c2]   [j_p + Amp I2]

with s1=1/2-lambda and s2=1/2-2 lambda.  Here ``m_p`` and ``j_p``
are *exactly the normalized entry quantities appearing in the public RHS*:

    m_p = M(X_p)/(X_p e_b f),
    j_p = J(X_p)/(X_p H(X_p) e_b f).

The source states that beta_1,beta_2 are non-negative C-infinity bumps of
identical width 0.3 centered at 13/lambda-3 and 13/lambda-1.  The corrected
reader does not uniquely prescribe their pointwise shape in this interface.
This module therefore uses one frozen repository-autonomous standard C-infinity
shape with full support width 0.3 and peak one.  That choice is part of the
semantic identity and is NOT claimed to recover a hidden source bump.

Numerically, the source explicitly recommends dividing each row by its
exponential value at the first bump center.  We implement that algebraically:
all matrix and main-pulse moments are evaluated directly in the row-scaled
coordinates, so no enormous intermediate exp(s*y) is formed.

This component deliberately accepts Amp, m_p and j_p as explicit inputs.  It
does not invent a current-lineage J value, does not identify the repository's
principal Amp proxy with the source-exact amplitude root, and does not compose
the compensators into the Cartesian field.  It is a low-dimensional executable
building block for that later composition, not PDE validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from .kokuno_source_main_axial_pulse_kernel import main_kernel_R0


SCHEMA = "kokuno-public-pulse-end-compensator-v1"
PARENT_EXACT_HEAD = "45da043dd2b4cd067f005a72c7e21fd0f2bcf309"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected 208-page reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"

# Current A1 outer schedule inherited from the #1107 lineage.  This is a
# repository-autonomous schedule value, not a newly fitted pulse parameter.
CURRENT_LAMBDA = 0.05
PUBLIC_BUMP_WIDTH = 0.3
MAIN_XI_END = 11.0
PULSE_XI_END = 13.0
QUADRATURE_ORDER = 96

_SOURCE_FORMULAS = {
    "pulse_ratio": "R_b=Amp(eta)R_0(lambda y)+c1(eta)beta1(y)+c2(eta)beta2(y)",
    "bump_centers": "y1=13/lambda-3, y2=13/lambda-1",
    "bump_width": "beta1,beta2 have identical width 0.3",
    "row_slopes": "s1=1/2-lambda, s2=1/2-2lambda",
    "row1": "sum_j c_j int exp(s1*y) beta_j dy = -m_p-Amp int exp(s1*y)R0(lambda*y)dy",
    "row2": "sum_j c_j int exp(s2*y) beta_j dy = -j_p-Amp int exp(s2*y)R0(lambda*y)dy",
    "m_entry": "m_p=M(X_p)/(X_p e_b f)",
    "j_entry": "j_p=J(X_p)/(X_p H(X_p) e_b f)",
    "source_stabilization": "divide each row by its exponential value at the first bump center",
}

_NUMERICAL_REALIZATION = {
    "bump_shape": (
        "repository-autonomous standard C-infinity bump: exp(1-1/(1-z^2)) "
        "for |z|<1, z=2(y-center)/0.3; full support width=0.3, peak=1"
    ),
    "quadrature": (
        "fixed 96-point Gauss-Legendre; bump moments integrated on their exact "
        "compact supports; R0 moments split at xi=0.02,10,11"
    ),
    "row_scaling": "publicly prescribed first-bump-center exponential scaling",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_end_compensation_linear_system_materialized": True,
    "public_bump_centers_and_width_materialized": True,
    "repository_autonomous_bump_shape_materialized": True,
    "source_exact_bump_shape_recovered": False,
    "source_exact_amplitude_root_materialized": False,
    "current_lineage_J_entry_materialized": False,
    "current_cartesian_end_compensation_composed": False,
    "source_hidden_parameters_recovered": False,
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


def _gauss_legendre_integral(func, a: float, b: float) -> float:
    if not (math.isfinite(a) and math.isfinite(b) and b >= a):
        raise ValueError("quadrature bounds must be finite and ordered")
    if b == a:
        return 0.0
    nodes, weights = np.polynomial.legendre.leggauss(QUADRATURE_ORDER)
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    x = half * nodes + mid
    return float(half * np.sum(weights * np.asarray(func(x), dtype=float)))


def _standard_bump(y: Any, center: float, width: float = PUBLIC_BUMP_WIDTH) -> np.ndarray:
    """Frozen autonomous C-infinity bump with *full* support width ``width``."""
    yy = _finite(y, "y")
    if not (math.isfinite(center) and math.isfinite(width) and width > 0.0):
        raise ValueError("center and width must be finite, with width>0")
    z = 2.0 * (yy - center) / width
    out = np.zeros_like(z, dtype=float)
    inside = np.abs(z) < 1.0
    if np.any(inside):
        q = z[inside]
        out[inside] = np.exp(1.0 - 1.0 / (1.0 - q * q))
    return out


def _standard_bump_prime(y: Any, center: float, width: float = PUBLIC_BUMP_WIDTH) -> np.ndarray:
    """Analytic y derivative of the frozen autonomous bump."""
    yy = _finite(y, "y")
    z = 2.0 * (yy - center) / width
    out = np.zeros_like(z, dtype=float)
    inside = np.abs(z) < 1.0
    if np.any(inside):
        q = z[inside]
        b = np.exp(1.0 - 1.0 / (1.0 - q * q))
        out[inside] = b * (-4.0 * q / (width * (1.0 - q * q) ** 2))
    return out


@dataclass(frozen=True)
class PulseEndCompensatorSolution:
    """Vectorized affine solution and deterministic algebra diagnostics."""

    c1: np.ndarray
    c2: np.ndarray
    residual_row1: np.ndarray
    residual_row2: np.ndarray


class KokunoPublicPulseEndCompensator:
    """Public two-row pulse-end closure algebra with frozen autonomous bumps."""

    def __init__(self, lambda_value: float = CURRENT_LAMBDA) -> None:
        lam = float(lambda_value)
        if not (math.isfinite(lam) and 0.0 < lam < 0.25):
            raise ValueError("lambda_value must be finite and satisfy 0<lambda<0.25")
        # The current A1 identity is deliberately frozen to the #1107 schedule.
        if lam != CURRENT_LAMBDA:
            raise ValueError(
                f"this current-lineage component freezes lambda={CURRENT_LAMBDA}; "
                "a different lambda requires a new semantic identity"
            )
        self.lambda_value = lam
        self.s1 = 0.5 - lam
        self.s2 = 0.5 - 2.0 * lam
        self.y1 = 13.0 / lam - 3.0
        self.y2 = 13.0 / lam - 1.0
        self.width = PUBLIC_BUMP_WIDTH
        self._matrix_scaled, self._main_scaled = self._build_scaled_system()
        det = float(np.linalg.det(self._matrix_scaled))
        if not (math.isfinite(det) and abs(det) > 1.0e-12):
            raise RuntimeError("row-scaled pulse end-compensator matrix is singular")

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
    def matrix_scaled(self) -> np.ndarray:
        return self._matrix_scaled.copy()

    @property
    def main_moments_scaled(self) -> np.ndarray:
        return self._main_scaled.copy()

    @property
    def row_condition_number(self) -> float:
        return float(np.linalg.cond(self._matrix_scaled))

    def beta1(self, y: Any) -> np.ndarray:
        return _standard_bump(y, self.y1, self.width)

    def beta2(self, y: Any) -> np.ndarray:
        return _standard_bump(y, self.y2, self.width)

    def beta1_prime(self, y: Any) -> np.ndarray:
        return _standard_bump_prime(y, self.y1, self.width)

    def beta2_prime(self, y: Any) -> np.ndarray:
        return _standard_bump_prime(y, self.y2, self.width)

    def _build_scaled_system(self) -> tuple[np.ndarray, np.ndarray]:
        slopes = (self.s1, self.s2)
        centers = (self.y1, self.y2)
        matrix = np.empty((2, 2), dtype=float)
        main = np.empty(2, dtype=float)
        half_width = 0.5 * self.width

        for i, slope in enumerate(slopes):
            # Source-prescribed row scaling: exp(-s_i*y1).
            for j, center in enumerate(centers):
                lo, hi = center - half_width, center + half_width
                matrix[i, j] = _gauss_legendre_integral(
                    lambda y, c=center, s=slope: (
                        np.exp(s * (y - self.y1)) * _standard_bump(y, c, self.width)
                    ),
                    lo,
                    hi,
                )

            # Substitute xi=lambda*y, then apply the same row scale directly.
            def integrand(xi: np.ndarray, s: float = slope) -> np.ndarray:
                return (
                    main_kernel_R0(xi)
                    * np.exp((s / self.lambda_value) * xi - s * self.y1)
                    / self.lambda_value
                )

            main[i] = sum(
                _gauss_legendre_integral(integrand, a, b)
                for a, b in ((0.0, 0.02), (0.02, 10.0), (10.0, MAIN_XI_END))
            )
        return matrix, main

    def _scaled_rhs(self, amplitude: Any, m_entry_ratio: Any, j_entry_ratio: Any) -> np.ndarray:
        amp, m0, j0 = np.broadcast_arrays(
            _finite(amplitude, "amplitude"),
            _finite(m_entry_ratio, "m_entry_ratio"),
            _finite(j_entry_ratio, "j_entry_ratio"),
        )
        rhs1 = -(
            m0 * math.exp(-self.s1 * self.y1)
            + amp * self._main_scaled[0]
        )
        rhs2 = -(
            j0 * math.exp(-self.s2 * self.y1)
            + amp * self._main_scaled[1]
        )
        return np.stack((rhs1, rhs2), axis=0)

    def solve(self, amplitude: Any, m_entry_ratio: Any, j_entry_ratio: Any) -> PulseEndCompensatorSolution:
        """Solve the public affine c1/c2 equations, vectorized over inputs."""
        rhs = self._scaled_rhs(amplitude, m_entry_ratio, j_entry_ratio)
        shape = rhs.shape[1:]
        coeff = np.linalg.solve(self._matrix_scaled, rhs.reshape(2, -1)).reshape((2,) + shape)
        residual = np.einsum("ij,j...->i...", self._matrix_scaled, coeff) - rhs
        return PulseEndCompensatorSolution(
            c1=np.asarray(coeff[0]),
            c2=np.asarray(coeff[1]),
            residual_row1=np.asarray(residual[0]),
            residual_row2=np.asarray(residual[1]),
        )

    def solve_eta_jet(
        self,
        amplitude: Any,
        m_entry_ratio: Any,
        j_entry_ratio: Any,
        amplitude_eta: Any,
        m_entry_ratio_eta: Any,
        j_entry_ratio_eta: Any,
    ) -> tuple[PulseEndCompensatorSolution, np.ndarray, np.ndarray]:
        """Return c plus analytic eta derivatives from the affine source system.

        The row matrix depends only on the frozen radial schedule, so eta
        differentiation acts only on Amp, m_p and j_p.
        """
        sol = self.solve(amplitude, m_entry_ratio, j_entry_ratio)
        deriv_rhs = self._scaled_rhs(amplitude_eta, m_entry_ratio_eta, j_entry_ratio_eta)
        shape = deriv_rhs.shape[1:]
        coeff_eta = np.linalg.solve(
            self._matrix_scaled, deriv_rhs.reshape(2, -1)
        ).reshape((2,) + shape)
        return sol, np.asarray(coeff_eta[0]), np.asarray(coeff_eta[1])

    def ratio_correction(self, y: Any, c1: Any, c2: Any) -> np.ndarray:
        """Evaluate c1*beta1+c2*beta2 with NumPy broadcasting."""
        yy = _finite(y, "y")
        a1, a2 = np.broadcast_arrays(_finite(c1, "c1"), _finite(c2, "c2"))
        return a1 * self.beta1(yy) + a2 * self.beta2(yy)

    def ratio_correction_y(self, y: Any, c1: Any, c2: Any) -> np.ndarray:
        """Analytic y derivative of the frozen end correction."""
        yy = _finite(y, "y")
        a1, a2 = np.broadcast_arrays(_finite(c1, "c1"), _finite(c2, "c2"))
        return a1 * self.beta1_prime(yy) + a2 * self.beta2_prime(yy)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": SOURCE_RELEASE,
                "release_date": SOURCE_RELEASE_DATE,
            },
            "lambda_value": self.lambda_value,
            "public_geometry": {
                "bump_width": self.width,
                "bump_centers": [self.y1, self.y2],
                "main_xi_end": MAIN_XI_END,
                "pulse_xi_end": PULSE_XI_END,
                "row_slopes": [self.s1, self.s2],
            },
            "numerical_realization": self.numerical_realization,
            "truth_boundary": self.truth_boundary,
        }

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    @classmethod
    def from_configuration(cls, config: Mapping[str, Any]) -> "KokunoPublicPulseEndCompensator":
        payload = dict(config)
        if payload.get("schema") != SCHEMA:
            raise ValueError("unexpected pulse end-compensator schema")
        if payload.get("parent_exact_head") != PARENT_EXACT_HEAD:
            raise ValueError("pulse end-compensator parent identity mismatch")
        expected_source = {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "release": SOURCE_RELEASE,
            "release_date": SOURCE_RELEASE_DATE,
        }
        if payload.get("source") != expected_source:
            raise ValueError("pulse end-compensator source provenance mismatch")
        if payload.get("lambda_value") != CURRENT_LAMBDA:
            raise ValueError("pulse end-compensator lambda identity mismatch")
        obj = cls(lambda_value=CURRENT_LAMBDA)
        if payload.get("public_geometry") != obj.configuration()["public_geometry"]:
            raise ValueError("pulse end-compensator public geometry mismatch")
        if payload.get("numerical_realization") != obj.numerical_realization:
            raise ValueError("pulse end-compensator numerical realization mismatch")
        if payload.get("truth_boundary") != obj.truth_boundary:
            raise ValueError("pulse end-compensator truth boundary mismatch")
        return obj
