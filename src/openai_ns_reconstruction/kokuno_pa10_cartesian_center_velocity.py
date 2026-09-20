"""Executable Cartesian spacetime velocity for the public PA.10 contraction center.

Pinned public provenance:
KokunoYumeto/yang-mills-interacting-workbench@
143f6773feb424ad9ed3a8d116653200f20346b7,
navier-stokes/navier_stokes_workbench.tex, corrected 2026-09-09 reader.

The public reconstruction gives, for tau=1-t,

    z = q^D eta,       tau = q(1-eta^2),       X = (x^2+y^2)/(2q),
    u_1 = (v_0/(2q)) x - q^(-A-1/2) F y,
    u_2 = (v_0/(2q)) y + q^(-A-1/2) F x,
    u_3 = q^(-A) U.

This module applies those identities to the source-C-normalized PA.10
contraction center delivered by the preceding Agent-1 stack.  It is therefore
an executable *inner center* velocity, not the final corrected fixed point,
not the globally joined Kokuno field, and not an independently PDE-validated
candidate.

The physical-coordinate inverse uses the source uniqueness proof directly.
Eliminating eta gives

    g(q) = q - z^2 q^(2h) - tau = 0,

with q_* = |z|^(1/D), g(q_*)=-tau, and g'(q)>=1-2h for q>=q_*.
Hence [q_*, q_* + tau/(1-2h)] is a certified bracket.  A fixed bisection
realization is used only to evaluate this public coordinate map.

Field-semantic identity deliberately excludes the expensive source-C
certificate execution knobs (certificate_intervals/chunk_size).  Those knobs
affect evidence receipts, not the executable velocity values.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    KokunoPA10PhysicalCenterProfileContract,
)
from .kokuno_pa10_source_c_normalized_physical_center import (
    KokunoPA10SourceCNormalizedPhysicalCenter,
)


SCHEMA = "kokuno-pa10-cartesian-center-velocity-v1"
DEFAULT_TIME_INTERVAL = (0.25, 0.75)
DEFAULT_BISECTION_ITERATIONS = 80

_SOURCE_FORMULAS = {
    "coordinates": (
        "tau=1-t; z=q^D eta; tau=q(1-eta^2); "
        "X=(x^2+y^2)/(2q); A=1/2+h; D=1/2-h"
    ),
    "q_inverse": (
        "q solves q-z^2 q^(2h)-tau=0; "
        "q_*=|z|^(1/D); g'(q)>=1-2h for q>=q_*"
    ),
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "profile_input": "F=F_0(X,eta); U=U_0(X,eta); v0=v_0(X,eta)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_coordinate_map_executable": True,
    "source_coordinate_inverse_bisection_executable": True,
    "inner_cartesian_spacetime_center_velocity_executable": True,
    "inner_cartesian_spacetime_center_velocity_vectorized": True,
    "axis_regular_cartesian_formula_executable": True,
    "field_semantic_identity_separate_from_evidence_receipt_identity": True,
    "source_C_selected_from_residual": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_cartesian_spacetime_leading_velocity_materialized": False,
    "outer_join_localization_materialized": False,
    "matched_global_pressure_materialized": False,
    "complete_kokuno_composite_velocity": False,
    "unified_cartesian_velocity_export_ready": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_array(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoPA10CartesianCenterVelocity:
    """Source-native inner Cartesian velocity obtained from the PA.10 center."""

    source_center: KokunoPA10SourceCNormalizedPhysicalCenter = field(
        default_factory=KokunoPA10SourceCNormalizedPhysicalCenter
    )
    time_interval: tuple[float, float] = DEFAULT_TIME_INTERVAL
    bisection_iterations: int = DEFAULT_BISECTION_ITERATIONS

    def __post_init__(self) -> None:
        if not isinstance(self.source_center, KokunoPA10SourceCNormalizedPhysicalCenter):
            raise TypeError(
                "source_center must be KokunoPA10SourceCNormalizedPhysicalCenter"
            )
        if len(tuple(self.time_interval)) != 2:
            raise ValueError("time_interval must contain exactly two values")
        t0, t1 = (float(v) for v in self.time_interval)
        if not (math.isfinite(t0) and math.isfinite(t1) and 0.0 <= t0 < t1 < 1.0):
            raise ValueError("time_interval must satisfy 0<=t0<t1<1")
        iterations = int(self.bisection_iterations)
        if (
            isinstance(self.bisection_iterations, bool)
            or iterations != self.bisection_iterations
            or not 56 <= iterations <= 128
        ):
            raise ValueError("bisection_iterations must be an integer in [56,128]")
        object.__setattr__(self, "time_interval", (t0, t1))
        object.__setattr__(self, "bisection_iterations", iterations)

    @property
    def physical_profiles(self) -> KokunoPA10PhysicalCenterProfileContract:
        return self.source_center.physical_profiles

    @property
    def axis_domain(self):
        return self.source_center.axis_domain

    @property
    def h(self) -> float:
        return float(self.axis_domain.h)

    @property
    def A(self) -> float:
        return float(self.axis_domain.A)

    @property
    def D(self) -> float:
        return float(self.axis_domain.D)

    @property
    def source_X_interval(self) -> tuple[float, float]:
        return self.physical_profiles.source_X_interval

    def _broadcast_xyzt(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            _finite_array(t, "t"),
        )

    def _validate_time(self, t: np.ndarray) -> None:
        t0, t1 = self.time_interval
        if np.any((t < t0) | (t > t1)):
            raise ValueError(f"t must lie in the registered interval [{t0},{t1}]")

    def solve_q(self, z: Any, t: Any) -> np.ndarray:
        """Solve the source implicit q-coordinate on the registered time window."""
        z_arr, t_arr = np.broadcast_arrays(
            _finite_array(z, "z"), _finite_array(t, "t")
        )
        self._validate_time(t_arr)
        tau = 1.0 - t_arr
        derivative_floor = 1.0 - 2.0 * self.h
        if derivative_floor <= 0.0:
            raise RuntimeError("source q uniqueness derivative floor is nonpositive")

        q_star = np.power(np.abs(z_arr), 1.0 / self.D)
        if np.any(~np.isfinite(q_star)):
            raise OverflowError("q_* is outside binary64 range")
        low = q_star
        high = q_star + tau / derivative_floor
        high = np.nextafter(high, np.inf)

        def g(q: np.ndarray) -> np.ndarray:
            return q - z_arr * z_arr * np.power(q, 2.0 * self.h) - tau

        if np.any(g(low) > 64.0 * np.finfo(float).eps * (np.abs(low) + tau + 1.0)):
            raise RuntimeError("source q lower bracket failed")
        if np.any(g(high) < -256.0 * np.finfo(float).eps * (np.abs(high) + tau + 1.0)):
            raise RuntimeError("source q upper bracket failed")

        for _ in range(self.bisection_iterations):
            mid = low + 0.5 * (high - low)
            move_low = g(mid) <= 0.0
            low = np.where(move_low, mid, low)
            high = np.where(move_low, high, mid)
        q = low + 0.5 * (high - low)
        if np.any(~np.isfinite(q)) or np.any(q <= 0.0):
            raise RuntimeError("source q inverse returned a nonpositive/nonfinite value")
        return q

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        """Map Cartesian spacetime coordinates to the public (q,X,eta) variables."""
        x_arr, y_arr, z_arr, t_arr = self._broadcast_xyzt(x, y, z, t)
        self._validate_time(t_arr)
        tau = 1.0 - t_arr
        q = self.solve_q(z_arr, t_arr)
        eta = z_arr / np.power(q, self.D)
        X = (x_arr * x_arr + y_arr * y_arr) / (2.0 * q)
        residual = q - z_arr * z_arr * np.power(q, 2.0 * self.h) - tau

        if np.any(np.abs(eta) >= 1.0 + 2.0e-13):
            raise RuntimeError("source q inverse returned |eta|>=1")
        return {
            "x": x_arr,
            "y": y_arr,
            "z": z_arr,
            "t": t_arr,
            "tau": tau,
            "q": q,
            "X": X,
            "eta": eta,
            "q_equation_residual": residual,
        }

    def cartesian_from_similarity(
        self, X: Any, eta: Any, t: Any, theta: Any = 0.0
    ) -> dict[str, np.ndarray]:
        """Evaluate the public forward coordinate map on the inner source domain."""
        X_arr, eta_arr, t_arr, theta_arr = np.broadcast_arrays(
            _finite_array(X, "X"),
            _finite_array(eta, "eta"),
            _finite_array(t, "t"),
            _finite_array(theta, "theta"),
        )
        self._validate_time(t_arr)
        x0, x1 = self.source_X_interval
        if np.any((X_arr < x0) | (X_arr > x1)):
            raise ValueError(f"X must lie in [{x0},{x1}]")
        if np.any(np.abs(eta_arr) >= 1.0):
            raise ValueError("physical source eta must satisfy |eta|<1")
        tau = 1.0 - t_arr
        q = tau / (1.0 - eta_arr * eta_arr)
        z = np.power(q, self.D) * eta_arr
        r = np.sqrt(2.0 * q * X_arr)
        return {
            "x": r * np.cos(theta_arr),
            "y": r * np.sin(theta_arr),
            "z": z,
            "t": t_arr,
            "q": q,
            "X": X_arr,
            "eta": eta_arr,
            "theta": theta_arr,
        }

    def values(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Return coordinates, center profiles, and Cartesian velocity components."""
        coords = self.similarity_coordinates(x, y, z, t)
        X = coords["X"]
        eta = coords["eta"]
        x0, x1 = self.source_X_interval
        radial_tol = 64.0 * np.finfo(float).eps * max(1.0, abs(x1))
        if np.any(X < x0 - radial_tol) or np.any(X > x1 + radial_tol):
            raise ValueError(
                "point lies outside the source inner X interval; "
                "no outer/global continuation is materialized here"
            )
        X_eval = np.clip(X, x0, x1)
        profiles = self.source_center.values(X_eval, eta)
        q = coords["q"]
        q_radial = profiles["v_0"] / (2.0 * q)
        q_swirl = np.power(q, -self.A - 0.5) * profiles["F_0"]
        u1 = q_radial * coords["x"] - q_swirl * coords["y"]
        u2 = q_radial * coords["y"] + q_swirl * coords["x"]
        u3 = np.power(q, -self.A) * profiles["U_0"]
        return {
            **coords,
            "F_0": profiles["F_0"],
            "U_0": profiles["U_0"],
            "v_0": profiles["v_0"],
            "u": u1,
            "v": u2,
            "w": u3,
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        values = self.values(x, y, z, t)
        return np.stack((values["u"], values["v"], values["w"]), axis=-1)

    def field_configuration(self) -> dict[str, Any]:
        """Configuration whose hash changes iff executable field semantics change."""
        return {
            "schema": SCHEMA,
            "physical_profile_configuration": self.physical_profiles.configuration(),
            "time_interval": [float(v) for v in self.time_interval],
            "q_inverse": {
                "method": "source-monotone-certified-bracket-bisection",
                "iterations": int(self.bisection_iterations),
            },
        }

    @property
    def field_sha256(self) -> str:
        return hashlib.sha256(
            _canonical_json(self.field_configuration()).encode("utf-8")
        ).hexdigest()

    def evidence_configuration(self) -> dict[str, Any]:
        """Evidence settings are intentionally separate from field identity."""
        return {
            "field_sha256": self.field_sha256,
            "source_center_evidence_configuration": self.source_center.configuration(),
        }

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.field_configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10CartesianCenterVelocity":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        physical = KokunoPA10PhysicalCenterProfileContract.from_configuration(
            payload["physical_profile_configuration"]
        )
        q_inverse = dict(payload["q_inverse"])
        if q_inverse.get("method") != "source-monotone-certified-bracket-bisection":
            raise ValueError("unsupported q inverse method")
        return cls(
            source_center=KokunoPA10SourceCNormalizedPhysicalCenter(
                physical_profiles=physical
            ),
            time_interval=tuple(float(v) for v in payload["time_interval"]),
            bisection_iterations=int(q_inverse["iterations"]),
        )

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10CartesianCenterVelocity":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        _, x1 = self.source_X_interval
        X = np.asarray([0.0, 0.15 * x1, 0.35 * x1, 0.55 * x1], dtype=float)
        eta = np.asarray([0.0, 0.2, -0.35, 0.45], dtype=float)
        t = np.asarray([0.25, 0.4, 0.6, 0.75], dtype=float)
        theta = np.asarray([0.0, 0.4, 1.1, 2.0], dtype=float)
        physical = self.cartesian_from_similarity(X, eta, t, theta)
        values = self.values(
            physical["x"], physical["y"], physical["z"], physical["t"]
        )
        scale = (
            np.abs(values["q"])
            + np.abs(values["z"] * values["z"] * np.power(values["q"], 2.0 * self.h))
            + np.abs(values["tau"])
            + 1.0
        )
        q_relative_residual = np.max(
            np.abs(values["q_equation_residual"]) / scale
        )
        axis = self.cartesian_from_similarity(
            np.zeros(3),
            np.asarray([-0.3, 0.0, 0.3]),
            np.asarray([0.3, 0.5, 0.7]),
        )
        axis_velocity = self.velocity(
            axis["x"], axis["y"], axis["z"], axis["t"]
        )
        probe_velocity = self.velocity(
            physical["x"], physical["y"], physical["z"], physical["t"]
        )
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "field_configuration": self.field_configuration(),
            "field_sha256": self.field_sha256,
            "evidence_configuration": self.evidence_configuration(),
            "machine_checks": {
                "coordinate_roundtrip_max_abs_X": float(
                    np.max(np.abs(values["X"] - X))
                ),
                "coordinate_roundtrip_max_abs_eta": float(
                    np.max(np.abs(values["eta"] - eta))
                ),
                "q_equation_relative_residual_max": float(q_relative_residual),
                "axis_transverse_velocity_exact_zero": bool(
                    np.all(axis_velocity[..., :2] == 0.0)
                ),
                "velocity_nontrivial_on_probe": bool(
                    np.any(np.abs(probe_velocity) > 0.0)
                ),
                "all_probe_values_finite": bool(
                    all(np.all(np.isfinite(value)) for value in values.values())
                ),
            },
            "scientific_gates": {
                "momentum_max_l2": 1.0e-3,
                "divergence_max_l2": 1.0e-5,
                "free_residual_defined_forcing_forbidden": True,
            },
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    field = KokunoPA10CartesianCenterVelocity()
    payload = field.save_report(args.output)
    if args.config_output is not None:
        field.save_configuration(args.config_output)
    print("field_sha256=", payload["field_sha256"])
    print("receipt_sha256=", payload["receipt_sha256"])
    print("machine_checks=", payload["machine_checks"])


if __name__ == "__main__":
    _main()
