"""Kokuno-native similarity coordinates for the public leading-field reconstruction.

The pinned public reconstruction uses

    tau = 1 - t = q (1 - eta**2),
    z = q**(1/2-h) eta,
    X = (x**2 + y**2) / (2 q),

hence q is the unique positive root above the source lower endpoint of

    q - z**2 q**(2h) = tau.

This module implements that convention directly.  It deliberately does not identify
it with the repository's separate user-Eq45 convention containing q**(-2h), and it
does not construct a complete three-dimensional Kokuno velocity by itself.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
CORRECTED_PDF = "released_ns_reader_corrected_20260909.pdf"
SCHEMA = "kokuno-native-similarity-coordinates-v1"

_SOURCE_FORMULAS = {
    "tau": "1-t",
    "A": "1/2+h",
    "D": "1/2-h",
    "z": "q^D*eta",
    "tau_similarity": "q*(1-eta^2)",
    "q_root": "q-z^2*q^(2h)=tau",
    "X": "(x^2+y^2)/(2q)",
    "d": "1-eta^2",
    "L": "1-2h*eta^2",
    "q_t": "-1/L",
    "eta_t": "D*eta/(q*L)",
    "X_t": "X/(q*L)",
    "q_z": "2*eta*q^(1-D)/L",
    "eta_z": "d/(q^D*L)",
    "X_z": "-2*eta*X/(q^D*L)",
}

_COORDINATE_CONTRACT = {
    "native_relation": "1-t=q(1-eta^2)=q-z^2*q^(2h)",
    "root_domain": "q>|z|^(1/D), with q_star=0 at z=0",
    "source_uniqueness_reason": "g'(q)=1-2h*z^2*q^(2h-1)>=1-2h>0 above q_star",
    "current_repo_eq45_user_relation": "q-z^2*q^(-2h)=1-t",
    "mapping_status": "native_kokuno_coordinates_executable; cross-convention_identification_not_claimed",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "native_q_eta_coordinates_executable": True,
    "cross_convention_equivalence_claimed": False,
    "full_leading_profile_reconstructed": False,
    "full_3d_velocity_candidate": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class KokunoNativeSimilarityCoordinates:
    """Vectorized solver and exact Jacobian for Kokuno's native coordinates.

    The source construction fixes ``0 < h < 1/100``.  The numerical root solve uses
    the fixed-point map ``q <- tau + z^2 q^(2h)`` starting above the source endpoint
    ``q_star``.  On the source domain its derivative is strictly below ``2h < .02``,
    so this is a strongly contractive, monotone-root-safe iteration rather than an
    unconstrained generic nonlinear solve.
    """

    h: float = 0.005
    rtol: float = 1.0e-13
    atol: float = 1.0e-14
    max_iterations: int = 64

    def __post_init__(self) -> None:
        h = float(self.h)
        rtol = float(self.rtol)
        atol = float(self.atol)
        max_iterations = int(self.max_iterations)
        if not np.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must be finite and satisfy 0 < h < 1e-2")
        if not np.isfinite(rtol) or not (0.0 < rtol <= 1.0e-8):
            raise ValueError("rtol must be finite and satisfy 0 < rtol <= 1e-8")
        if not np.isfinite(atol) or not (0.0 < atol <= 1.0e-10):
            raise ValueError("atol must be finite and satisfy 0 < atol <= 1e-10")
        if max_iterations < 8 or max_iterations > 256:
            raise ValueError("max_iterations must lie in [8, 256]")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "rtol", rtol)
        object.__setattr__(self, "atol", atol)
        object.__setattr__(self, "max_iterations", max_iterations)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def D(self) -> float:
        return 0.5 - self.h

    @staticmethod
    def _finite_array(value: Any, name: str) -> np.ndarray:
        array = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must contain only finite values")
        return array

    def q_star(self, z: Any) -> np.ndarray:
        z_array = self._finite_array(z, "z")
        result = np.power(np.abs(z_array), 1.0 / self.D)
        if not np.all(np.isfinite(result)):
            raise ValueError("z is too large for finite source-domain q_star")
        return result

    def solve_q(self, z: Any, t: Any) -> np.ndarray:
        """Return the unique Kokuno-native positive q for broadcastable ``z,t``."""

        z_array, t_array = np.broadcast_arrays(
            self._finite_array(z, "z"), self._finite_array(t, "t")
        )
        tau = 1.0 - t_array
        if not np.all(tau > 0.0):
            raise ValueError("Kokuno native coordinates require t < 1")

        q_star = np.power(np.abs(z_array), 1.0 / self.D)
        if not np.all(np.isfinite(q_star)):
            raise ValueError("z is too large for finite source-domain q_star")

        # q_star is the zero-tau endpoint. q_star + tau lies strictly inside the
        # monotone source domain, and the fixed-point map preserves that domain.
        q = q_star + tau
        z2 = z_array * z_array
        for _ in range(self.max_iterations):
            q_next = tau + z2 * np.power(q, 2.0 * self.h)
            tolerance = self.atol + self.rtol * np.abs(q_next)
            if np.all(np.abs(q_next - q) <= tolerance):
                q = q_next
                break
            q = q_next
        else:
            raise RuntimeError("Kokuno native q fixed-point solve did not converge")

        if not np.all(np.isfinite(q)) or not np.all(q > q_star) or not np.all(q > 0.0):
            raise RuntimeError("Kokuno native q solve left the unique positive source domain")

        residual = q - z2 * np.power(q, 2.0 * self.h) - tau
        residual_bound = 8.0 * (self.atol + self.rtol * np.maximum(np.abs(q), tau))
        if not np.all(np.abs(residual) <= residual_bound):
            raise RuntimeError("Kokuno native q solve failed its defining-equation residual check")
        return q

    def eta(self, z: Any, t: Any) -> np.ndarray:
        z_array, t_array = np.broadcast_arrays(
            self._finite_array(z, "z"), self._finite_array(t, "t")
        )
        q = self.solve_q(z_array, t_array)
        eta = z_array / np.power(q, self.D)
        if not np.all(np.abs(eta) < 1.0):
            raise RuntimeError("Kokuno native q solve did not produce |eta| < 1 for t < 1")
        return eta

    def evaluate(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Evaluate ``q,eta,X`` and the source exact fixed-physical-coordinate Jacobian."""

        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            self._finite_array(x, "x"),
            self._finite_array(y, "y"),
            self._finite_array(z, "z"),
            self._finite_array(t, "t"),
        )
        q = self.solve_q(z_array, t_array)
        q_D = np.power(q, self.D)
        eta = z_array / q_D
        X = (x_array * x_array + y_array * y_array) / (2.0 * q)
        d = 1.0 - eta * eta
        L = 1.0 - 2.0 * self.h * eta * eta

        if not np.all(d > 0.0) or not np.all(L > 0.0):
            raise RuntimeError("Kokuno native coordinate invariants left their source domain")

        return {
            "q": q,
            "eta": eta,
            "X": X,
            "d": d,
            "L": L,
            "q_t": -1.0 / L,
            "eta_t": self.D * eta / (q * L),
            "X_t": X / (q * L),
            "q_z": 2.0 * eta * np.power(q, 1.0 - self.D) / L,
            "eta_z": d / (q_D * L),
            "X_z": -2.0 * eta * X / (q_D * L),
        }

    def relation_residual(self, z: Any, t: Any) -> np.ndarray:
        z_array, t_array = np.broadcast_arrays(
            self._finite_array(z, "z"), self._finite_array(t, "t")
        )
        q = self.solve_q(z_array, t_array)
        return q - z_array * z_array * np.power(q, 2.0 * self.h) - (1.0 - t_array)

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "corrected_pdf": CORRECTED_PDF,
                "formula_evidence": "pinned_public_workbench_commit",
                "pdf_independently_parsed": False,
            },
            "parameters": {
                "h": self.h,
                "origin": "autonomous_demo_within_public_source_construction_bound",
                "source_bound": "0<h<1e-2",
            },
            "solver": {
                "method": "source-domain contractive fixed point q=tau+z^2*q^(2h)",
                "rtol": self.rtol,
                "atol": self.atol,
                "max_iterations": self.max_iterations,
                "derivative_bound": "0 <= 2h*z^2*q^(2h-1) < 2h < 0.02",
            },
            "formulas": dict(_SOURCE_FORMULAS),
            "coordinates": dict(_COORDINATE_CONTRACT),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoNativeSimilarityCoordinates":
        if not isinstance(payload, dict):
            raise ValueError("coordinate payload must be a JSON object")
        expected_keys = {"schema", "source", "parameters", "solver", "formulas", "coordinates", "truth_boundary"}
        if set(payload) - {"sha256"} != expected_keys:
            raise ValueError("coordinate payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported Kokuno native-coordinate schema")

        default_payload = cls().to_payload()
        if payload["source"] != default_payload["source"]:
            raise ValueError("Kokuno source/provenance metadata changed")
        if payload["formulas"] != _SOURCE_FORMULAS:
            raise ValueError("Kokuno coordinate formula metadata changed")
        if payload["coordinates"] != _COORDINATE_CONTRACT:
            raise ValueError("Kokuno coordinate-contract metadata changed")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("Kokuno truth-boundary metadata changed")

        parameters = payload["parameters"]
        if set(parameters) != {"h", "origin", "source_bound"}:
            raise ValueError("Kokuno coordinate parameter metadata changed")
        if parameters["origin"] != "autonomous_demo_within_public_source_construction_bound":
            raise ValueError("Kokuno coordinate parameter origin changed")
        if parameters["source_bound"] != "0<h<1e-2":
            raise ValueError("Kokuno coordinate source bound changed")

        solver = payload["solver"]
        expected_solver_keys = {"method", "rtol", "atol", "max_iterations", "derivative_bound"}
        if set(solver) != expected_solver_keys:
            raise ValueError("Kokuno coordinate solver metadata changed")
        if solver["method"] != default_payload["solver"]["method"]:
            raise ValueError("Kokuno coordinate solver method changed")
        if solver["derivative_bound"] != default_payload["solver"]["derivative_bound"]:
            raise ValueError("Kokuno coordinate contraction certificate changed")

        coordinates = cls(
            h=float(parameters["h"]),
            rtol=float(solver["rtol"]),
            atol=float(solver["atol"]),
            max_iterations=int(solver["max_iterations"]),
        )
        if "sha256" in payload and payload["sha256"] != coordinates.sha256:
            raise ValueError("Kokuno native-coordinate payload SHA mismatch")
        return coordinates

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoNativeSimilarityCoordinates":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_payload(payload)
