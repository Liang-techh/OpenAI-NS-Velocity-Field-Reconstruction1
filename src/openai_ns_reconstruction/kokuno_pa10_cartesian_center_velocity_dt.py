"""Analytic time derivative for the executable PA.10 Cartesian center velocity.

Pinned public provenance:
KokunoYumeto/yang-mills-interacting-workbench@
143f6773feb424ad9ed3a8d116653200f20346b7,
navier-stokes/navier_stokes_workbench.tex, corrected 2026-09-09 reader.

Agent-1 #803 materializes the public inner contraction-center velocity

    tau=1-t,
    q-z^2 q^(2h)-tau=0,
    X=(x^2+y^2)/(2q),       eta=z/q^D,

    u1=(v0/(2q))x-q^(-A-1/2)Fy,
    u2=(v0/(2q))y+q^(-A-1/2)Fx,
    u3=q^(-A)U.

For a fixed Cartesian point (x,y,z), implicit differentiation of the displayed
source coordinate equation gives, with L=1-2h eta^2,

    q_t=-1/L,
    X_t=X/(qL),
    eta_t=D eta/(qL).

This module uses those identities together with the already-executable analytic
physical-profile derivatives to expose ``velocity_dt``.  The only extra profile
quantity needed is ``(v0)_eta``.  At the PA.10 contraction center
``u0=Y s(eta)`` and therefore ``U0=U_*+X s``.  Differentiating the same public
incompressibility primitive used by #787 gives ``(v0)_eta`` analytically; no
finite-difference value is used in production.

Finite differences appear only in focused/report cross-checks.  This remains
an *inner contraction-center* derivative.  It does not materialize the final
fixed point, an outer/global join, pressure, restricted forcing, or an
independent Navier--Stokes validation.
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

from .kokuno_pa10_cartesian_center_velocity import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    KokunoPA10CartesianCenterVelocity,
)


SCHEMA = "kokuno-pa10-cartesian-center-velocity-dt-v1"

_SOURCE_FORMULAS = {
    "implicit_coordinate_time_derivatives": (
        "L=1-2h eta^2; q_t=-1/L; X_t=X/(qL); eta_t=D eta/(qL)"
    ),
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "center_radial_profile": (
        "u0=Y s(eta), s=-Z_*/(2L); U0=U_*+X s; "
        "v0 follows the displayed incompressibility primitive"
    ),
    "velocity_time_derivative": (
        "differentiate the displayed Cartesian velocity at fixed x,y,z using "
        "q_t,X_t,eta_t and analytic F_X,F_eta,U_X,U_eta,(v0)_X,(v0)_eta"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "inner_cartesian_center_velocity_inherited_from_agent1_803": True,
    "source_coordinate_time_derivatives_analytic": True,
    "center_v0_eta_analytic": True,
    "inner_cartesian_center_velocity_dt_executable": True,
    "inner_cartesian_center_velocity_dt_vectorized": True,
    "velocity_dt_uses_finite_difference_in_production": False,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_cartesian_spacetime_leading_velocity_materialized": False,
    "outer_join_localization_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "complete_candidate_api_ready": False,
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


def _relative_max(actual: np.ndarray, reference: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    reference = np.asarray(reference, dtype=float)
    scale = np.maximum(np.maximum(np.abs(actual), np.abs(reference)), 1.0)
    return float(np.max(np.abs(actual - reference) / scale))


@dataclass(frozen=True)
class KokunoPA10CartesianCenterVelocityTimeDerivative:
    """Add an analytic fixed-Cartesian-point time derivative to Agent-1 #803."""

    field: KokunoPA10CartesianCenterVelocity = field(
        default_factory=KokunoPA10CartesianCenterVelocity
    )

    def __post_init__(self) -> None:
        if not isinstance(self.field, KokunoPA10CartesianCenterVelocity):
            raise TypeError("field must be KokunoPA10CartesianCenterVelocity")

    @property
    def A(self) -> float:
        return float(self.field.A)

    @property
    def D(self) -> float:
        return float(self.field.D)

    @property
    def h(self) -> float:
        return float(self.field.h)

    @property
    def physical_profiles(self):
        return self.field.physical_profiles

    @property
    def axis_profiles(self):
        return self.physical_profiles.axis_profiles

    @property
    def axis_domain(self):
        return self.field.axis_domain

    @property
    def source_X_interval(self) -> tuple[float, float]:
        return self.field.source_X_interval

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Delegate the unchanged #803 velocity exactly."""
        return self.field.velocity(x, y, z, t)

    def coordinate_time_derivatives(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        """Return analytic q_t, X_t, eta_t at fixed Cartesian x,y,z."""
        coords = self.field.similarity_coordinates(x, y, z, t)
        eta = coords["eta"]
        q = coords["q"]
        L = 1.0 - 2.0 * self.h * eta * eta
        if np.any(L <= 0.0):
            raise RuntimeError("source coordinate L lost positivity")
        q_t = -1.0 / L
        X_t = coords["X"] / (q * L)
        eta_t = self.D * eta / (q * L)
        return {
            **coords,
            "L": L,
            "q_t": q_t,
            "X_t": X_t,
            "eta_t": eta_t,
        }

    def _slope_eta_eta(self, X: Any, eta: Any) -> np.ndarray:
        """Return s'' for u0(Y,eta)=Y*s(eta), using displayed source formulas."""
        X_arr, eta_arr = np.broadcast_arrays(
            np.asarray(X, dtype=float), np.asarray(eta, dtype=float)
        )
        Y = self.physical_profiles.Lambda * X_arr
        axis = self.axis_profiles.values(Y, eta_arr)
        axis_d = self.axis_profiles.derivatives(Y, eta_arr)

        A = float(self.axis_domain.A)
        D = float(self.axis_domain.D)
        h = float(self.axis_domain.h)
        j0 = float(self.axis_domain.j0)
        p2 = float(self.axis_domain.pressure_square)

        U = axis["U_star"]
        d = axis["d"]
        L = axis["L"]
        Z = axis["Z_star"]
        Z_eta = axis_d["Z_star_eta"]
        Pi_eta = axis["Pi_0_eta"]
        Pi_eta_eta = axis_d["Pi_0_eta_eta"]

        den = 1.0 + eta_arr * eta_arr
        Pi_eta_eta_eta = (
            24.0 * p2 * eta_arr * (5.0 * eta_arr * eta_arr - 3.0) / den**5
        )
        d_eta = -2.0 * eta_arr
        d_eta_eta = -2.0
        L_eta = -4.0 * h * eta_arr
        L_eta_eta = -4.0 * h
        H_eta_eta = -2.0 * j0 - 24.0 * eta_arr
        B_eta = -2.0 * U - 8.0 * eta_arr
        B_eta_eta = -16.0

        Z_eta_eta = (
            -A * (B_eta_eta * U + 8.0 * B_eta)
            - 4.0 * H_eta_eta
            - d_eta_eta * Pi_eta
            - 2.0 * d_eta * Pi_eta_eta
            - d * Pi_eta_eta_eta
            + 4.0 * A * (2.0 * Pi_eta + eta_arr * Pi_eta_eta)
        )

        quotient_eta_eta = (
            Z_eta_eta / L
            - Z * L_eta_eta / (L * L)
            - 2.0 * Z_eta * L_eta / (L * L)
            + 2.0 * Z * L_eta * L_eta / (L * L * L)
        )
        return -0.5 * quotient_eta_eta

    def v0_eta(self, X: Any, eta: Any) -> np.ndarray:
        """Analytic eta derivative of the PA.10 center radial profile v0."""
        values = self.physical_profiles.values(X, eta)
        X_arr = values["X"]
        eta_arr = values["eta"]
        Y = values["Y"]
        axis = self.axis_profiles.values(Y, eta_arr)
        axis_d = self.axis_profiles.derivatives(Y, eta_arr)

        L = axis["L"]
        L_eta = axis_d["L_eta"]
        d = axis["d"]
        d_eta = axis_d["d_eta"]
        slope = axis["u_0_Y"]
        slope_eta = axis_d["u_0_Y_eta"]
        slope_eta_eta = self._slope_eta_eta(X_arr, eta_arr)

        U = values["U_0"]
        U_eta = 4.0 + X_arr * slope_eta
        M_over_X = values["M_0_over_X"]
        M_over_X_eta = 4.0 + 0.5 * X_arr * slope_eta
        M_eta_over_X = 4.0 + 0.5 * X_arr * slope_eta
        M_eta_over_X_eta = 0.5 * X_arr * slope_eta_eta

        numerator = values["v_0"] * L
        numerator_eta = (
            2.0 * U
            + 2.0 * eta_arr * U_eta
            - 2.0 * self.D * (M_over_X + eta_arr * M_over_X_eta)
            - d_eta * M_eta_over_X
            - d * M_eta_over_X_eta
        )
        result = (numerator_eta * L - numerator * L_eta) / (L * L)
        if np.any(~np.isfinite(result)):
            raise RuntimeError("analytic v0_eta became nonfinite")
        return result

    def values(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Return #803 values plus analytic fixed-Cartesian-point time derivatives."""
        base = self.field.values(x, y, z, t)
        coord_dt = self.coordinate_time_derivatives(x, y, z, t)
        X = base["X"]
        eta = base["eta"]
        x0, x1 = self.source_X_interval
        radial_tol = 64.0 * np.finfo(float).eps * max(1.0, abs(x1))
        if np.any(X < x0 - radial_tol) or np.any(X > x1 + radial_tol):
            raise ValueError("point lies outside the source inner X interval")
        X_eval = np.clip(X, x0, x1)

        profile_d = self.physical_profiles.derivatives(X_eval, eta)
        v0_eta = self.v0_eta(X_eval, eta)
        q = base["q"]
        q_t = coord_dt["q_t"]
        X_t = coord_dt["X_t"]
        eta_t = coord_dt["eta_t"]

        F_t = profile_d["F_0_X"] * X_t + profile_d["F_0_eta"] * eta_t
        U_t = profile_d["U_0_X"] * X_t + profile_d["U_0_eta"] * eta_t
        v0_t = profile_d["v_0_X"] * X_t + v0_eta * eta_t

        radial = base["v_0"] / (2.0 * q)
        radial_t = v0_t / (2.0 * q) - base["v_0"] * q_t / (2.0 * q * q)

        swirl_scale = np.power(q, -self.A - 0.5)
        swirl = swirl_scale * base["F_0"]
        swirl_t = swirl_scale * (
            F_t + (-self.A - 0.5) * (q_t / q) * base["F_0"]
        )

        axial_scale = np.power(q, -self.A)
        axial_t = axial_scale * (U_t - self.A * (q_t / q) * base["U_0"])

        u_t = radial_t * base["x"] - swirl_t * base["y"]
        v_t = radial_t * base["y"] + swirl_t * base["x"]
        w_t = axial_t
        return {
            **base,
            "L_coordinate": coord_dt["L"],
            "q_t": q_t,
            "X_t": X_t,
            "eta_t": eta_t,
            "F_0_t": F_t,
            "U_0_t": U_t,
            "v_0_eta": v0_eta,
            "v_0_t": v0_t,
            "radial_coefficient": radial,
            "radial_coefficient_t": radial_t,
            "swirl_coefficient": swirl,
            "swirl_coefficient_t": swirl_t,
            "u_t": u_t,
            "v_t": v_t,
            "w_t": w_t,
        }

    def velocity_dt(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        values = self.values(x, y, z, t)
        return np.stack((values["u_t"], values["v_t"], values["w_t"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "field_configuration": self.field.field_configuration(),
            "derivative_realization": "analytic-source-chain-rule-v1",
        }

    @property
    def field_sha256(self) -> str:
        """The derivative wrapper does not change the underlying velocity identity."""
        return self.field.field_sha256

    @property
    def derivative_sha256(self) -> str:
        return hashlib.sha256(
            _canonical_json(self.configuration()).encode("utf-8")
        ).hexdigest()

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10CartesianCenterVelocityTimeDerivative":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        if payload.get("derivative_realization") != "analytic-source-chain-rule-v1":
            raise ValueError("unsupported derivative realization")
        return cls(
            field=KokunoPA10CartesianCenterVelocity.from_configuration(
                payload["field_configuration"]
            )
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10CartesianCenterVelocityTimeDerivative":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    @staticmethod
    def _fd4_time(field: KokunoPA10CartesianCenterVelocity, x, y, z, t, step: float):
        return (
            field.velocity(x, y, z, t - 2.0 * step)
            - 8.0 * field.velocity(x, y, z, t - step)
            + 8.0 * field.velocity(x, y, z, t + step)
            - field.velocity(x, y, z, t + 2.0 * step)
        ) / (12.0 * step)

    def report(self) -> dict[str, Any]:
        _, xmax = self.source_X_interval
        X = np.asarray([0.04, 0.13, 0.29, 0.47], dtype=float) * xmax
        eta = np.asarray([0.07, -0.19, 0.31, -0.43], dtype=float)
        t = np.asarray([0.31, 0.43, 0.57, 0.69], dtype=float)
        theta = np.asarray([0.2, 0.8, 1.4, 2.2], dtype=float)
        physical = self.field.cartesian_from_similarity(X, eta, t, theta)
        analytic = self.velocity_dt(
            physical["x"], physical["y"], physical["z"], physical["t"]
        )
        fd_coarse = self._fd4_time(
            self.field,
            physical["x"],
            physical["y"],
            physical["z"],
            physical["t"],
            2.0e-4,
        )
        fd_fine = self._fd4_time(
            self.field,
            physical["x"],
            physical["y"],
            physical["z"],
            physical["t"],
            1.0e-4,
        )

        profile_eta_step = 2.0e-5
        plus = self.physical_profiles.values(X, eta + profile_eta_step)["v_0"]
        minus = self.physical_profiles.values(X, eta - profile_eta_step)["v_0"]
        v0_eta_fd = (plus - minus) / (2.0 * profile_eta_step)
        v0_eta_analytic = self.v0_eta(X, eta)

        coord = self.coordinate_time_derivatives(
            physical["x"], physical["y"], physical["z"], physical["t"]
        )
        coord_step = 1.0e-5
        cplus = self.field.similarity_coordinates(
            physical["x"], physical["y"], physical["z"], physical["t"] + coord_step
        )
        cminus = self.field.similarity_coordinates(
            physical["x"], physical["y"], physical["z"], physical["t"] - coord_step
        )
        q_t_fd = (cplus["q"] - cminus["q"]) / (2.0 * coord_step)
        X_t_fd = (cplus["X"] - cminus["X"]) / (2.0 * coord_step)
        eta_t_fd = (cplus["eta"] - cminus["eta"]) / (2.0 * coord_step)

        axis = self.field.cartesian_from_similarity(
            np.zeros(3),
            np.asarray([-0.25, 0.0, 0.25]),
            np.asarray([0.32, 0.5, 0.68]),
        )
        axis_dt = self.velocity_dt(axis["x"], axis["y"], axis["z"], axis["t"])

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
            "configuration": self.configuration(),
            "field_sha256": self.field_sha256,
            "derivative_sha256": self.derivative_sha256,
            "machine_checks": {
                "velocity_dt_vs_fd4_fine_relative_max": _relative_max(analytic, fd_fine),
                "fd4_fine_vs_coarse_relative_max": _relative_max(fd_fine, fd_coarse),
                "v0_eta_vs_centered_fd_relative_max": _relative_max(
                    v0_eta_analytic, v0_eta_fd
                ),
                "q_t_vs_centered_fd_relative_max": _relative_max(coord["q_t"], q_t_fd),
                "X_t_vs_centered_fd_relative_max": _relative_max(coord["X_t"], X_t_fd),
                "eta_t_vs_centered_fd_relative_max": _relative_max(
                    coord["eta_t"], eta_t_fd
                ),
                "axis_transverse_velocity_dt_exact_zero": bool(
                    np.all(axis_dt[..., :2] == 0.0)
                ),
                "velocity_dt_nontrivial_on_probe": bool(np.any(np.abs(analytic) > 0.0)),
                "all_velocity_dt_probe_values_finite": bool(np.all(np.isfinite(analytic))),
            },
            "scientific_gates": {
                "momentum_max_l2": 1.0e-3,
                "divergence_max_l2": 1.0e-5,
                "free_residual_defined_forcing_forbidden": True,
            },
            "truth_boundary": self.truth_boundary,
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
    derivative = KokunoPA10CartesianCenterVelocityTimeDerivative()
    payload = derivative.save_report(args.output)
    if args.config_output is not None:
        derivative.save_configuration(args.config_output)
    print("field_sha256=", payload["field_sha256"])
    print("derivative_sha256=", payload["derivative_sha256"])
    print("receipt_sha256=", payload["receipt_sha256"])
    print("machine_checks=", payload["machine_checks"])


if __name__ == "__main__":
    _main()
