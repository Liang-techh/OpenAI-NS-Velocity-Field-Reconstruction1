"""Current-lineage Kokuno joined profile -> executable Cartesian leading velocity.

Pinned public provenance is the corrected 2026-09-09 KokunoYumeto
reconstruction at commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.

Agent 1 already has two governed ingredients on this ancestry:

* ``KokunoPA16CurrentJoinedProfile``: the current candidate-side ``F,U,E``
  profile on ``0 <= X <= X_h``;
* ``KokunoPA10CartesianCenterVelocity``: the public spacetime coordinate
  inverse and Cartesian assembly

      tau = 1-t,
      z = q^D eta,
      tau = q(1-eta^2),
      X = (x^2+y^2)/(2q),
      u1 = (v0/(2q)) x - q^(-A-1/2) F y,
      u2 = (v0/(2q)) y + q^(-A-1/2) F x,
      u3 = q^(-A) U.

This module binds those two existing components instead of inventing a second
coordinate map.  The only missing profile datum is the incompressibility
primitive needed for the public radial profile

    M(X,eta) = int_0^X U(s,eta) ds,
    v0 = [2 eta U - 2 D eta M/X - d M_eta/X] / L.

For ``X<=X_i`` the primitive is integrated piecewise across the already-fixed
PA.10 seams.  For ``X_i<X<=X_h`` it is continued in the scaled coordinate
``x=X/X_R``.  The long pre-restoration interval is integrated analytically,
while the public restoration / PA.16 intervals use deterministic
Gauss--Legendre quadrature.  ``partial_eta(M/X)`` is a frozen fourth-order
finite-difference realization because the current joined profile does not yet
expose an analytic eta jet.  This numerical derivative is explicitly not
source hidden data or an independent PDE validation.

The result is an executable *candidate-side leading velocity through X_h*.
It is not the outer/global completion beyond X_h, not matched pressure or
forcing, not the corrected fixed point, and not paper/OpenAI-field exactness.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_inner_join_exit import (
    LOG_JOIN_EXIT,
    LOG_REPAIR_START,
    LOG_RESTORE_END,
    LOG_RESTORE_START,
)
from .kokuno_pa10_cartesian_center_velocity import KokunoPA10CartesianCenterVelocity
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa16_current_joined_profile import KokunoPA16CurrentJoinedProfile


SCHEMA = "kokuno-pa16-current-cartesian-leading-velocity-v1"
DEFAULT_MOMENT_QUADRATURE_ORDER = 48
DEFAULT_ETA_FD_STEP = 2.0e-5

_SOURCE_FORMULAS = {
    "coordinates": (
        "tau=1-t; z=q^D eta; tau=q(1-eta^2); X=(x^2+y^2)/(2q); "
        "A=1/2+h; D=1/2-h"
    ),
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "radial_profile": (
        "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L"
    ),
    "joined_profile": "current Agent-1 F,U,E on 0<=X<=X_h",
}

_NUMERICAL_REALIZATION = {
    "coordinate_inverse": (
        "reuse KokunoPA10CartesianCenterVelocity certified monotone q-bisection; "
        "no duplicate inverse"
    ),
    "prefix_M": (
        "fixed piecewise Gauss-Legendre integration across current PA.10 seams"
    ),
    "joined_M": (
        "scaled-x continuation; pre-restoration/ideal plateaus analytic, "
        "restoration and PA16 repair fixed Gauss-Legendre"
    ),
    "M_eta_over_X": (
        "frozen fourth-order centered/one-sided eta finite difference of M/X"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "governed_source_coordinate_inverse_reused": True,
    "current_joined_F_U_E_consumed": True,
    "current_joined_incompressibility_primitive_materialized": True,
    "current_joined_radial_v0_materialized": True,
    "current_joined_cartesian_spacetime_leading_velocity_through_Xh_materialized": True,
    "axis_regular_cartesian_formula_executable": True,
    "velocity_interface_vectorized": True,
    "velocity_configuration_serializable": True,
    "eta_derivative_is_repository_numerical_realization": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_prepared_appendixA_pressure_stress_materialized": False,
    "source_admitted_kappa0_materialized": False,
    "source_global_inner_to_outer_join_admitted": False,
    "outer_global_leading_velocity_materialized": False,
    "velocity_beyond_Xh_materialized": False,
    "unified_global_cartesian_velocity_export_ready": False,
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


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@lru_cache(maxsize=8)
def _legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(order, bool) or not isinstance(order, int) or order < 16:
        raise ValueError("quadrature order must be an integer >=16")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes = np.asarray(nodes, dtype=float)
    weights = np.asarray(weights, dtype=float)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianLeadingVelocity:
    """Bind the current joined profile to the governed public Cartesian map."""

    joined: KokunoPA16CurrentJoinedProfile = field(
        default_factory=KokunoPA16CurrentJoinedProfile,
        repr=False,
        compare=False,
    )
    geometry: KokunoPA10CartesianCenterVelocity = field(
        default_factory=KokunoPA10CartesianCenterVelocity,
        repr=False,
        compare=False,
    )
    moment_quadrature_order: int = DEFAULT_MOMENT_QUADRATURE_ORDER
    eta_fd_step: float = DEFAULT_ETA_FD_STEP

    def __post_init__(self) -> None:
        if not isinstance(self.joined, KokunoPA16CurrentJoinedProfile):
            raise TypeError("joined must be KokunoPA16CurrentJoinedProfile")
        if not isinstance(self.geometry, KokunoPA10CartesianCenterVelocity):
            raise TypeError("geometry must be KokunoPA10CartesianCenterVelocity")
        q = int(self.moment_quadrature_order)
        if (
            isinstance(self.moment_quadrature_order, bool)
            or q != self.moment_quadrature_order
            or not 16 <= q <= 128
        ):
            raise ValueError("moment_quadrature_order must be an integer in [16,128]")
        step = float(self.eta_fd_step)
        if not math.isfinite(step) or not 1.0e-7 <= step <= 1.0e-3:
            raise ValueError("eta_fd_step must lie in [1e-7,1e-3]")
        if float(self.geometry.physical_profiles.C) != float(self.joined.moments.C):
            raise ValueError("Cartesian geometry and current joined profile must use the same source C")
        if not self.joined.route_ready:
            raise ValueError(
                "current numerical T_sh route does not fit the public source separation geometry"
            )
        object.__setattr__(self, "moment_quadrature_order", q)
        object.__setattr__(self, "eta_fd_step", step)

    @property
    def X_i(self) -> float:
        return float(self.joined.X_i)

    @property
    def X_h(self) -> float:
        return float(self.joined.X_h)

    @property
    def X_R(self) -> float:
        return float(self.joined.X_R)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.joined.eta_interval)

    @property
    def A(self) -> float:
        return float(self.geometry.A)

    @property
    def D(self) -> float:
        return float(self.geometry.D)

    def similarity_coordinates(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Reuse the governed public Cartesian -> ``(q,X,eta)`` inverse."""
        return self.geometry.similarity_coordinates(x, y, z, t)

    def _integrate_U_physical(self, a: float, b: float, eta: float) -> float:
        if not b > a:
            return 0.0
        nodes, weights = _legendre_rule(self.moment_quadrature_order)
        half = 0.5 * (b - a)
        mid = 0.5 * (a + b)
        X = mid + half * nodes
        vals = self.joined.profile_values(X, np.full_like(X, eta))
        U = np.asarray(vals["U_current_joined"], dtype=float)
        return float(half * np.dot(weights, U))

    def _prefix_M(self, X: float, eta: float) -> float:
        if X <= 0.0:
            return 0.0
        moments = self.joined.moments
        seams = (0.0, moments.X_0, moments.X_1, moments.X_b_ref, self.X_i)
        total = 0.0
        for left, right in zip(seams[:-1], seams[1:], strict=True):
            if X <= left:
                break
            stop = min(X, right)
            if stop > left:
                total += self._integrate_U_physical(left, stop, eta)
            if X <= right:
                break
        return float(total)

    def _selected_join(self):
        return self.joined.certificate.build_selected_inner_join(
            quadrature_points=self.joined.quadrature_points
        )

    def _integrate_scaled_U(
        self,
        a: float,
        b: float,
        *,
        eta: float,
        coefficients: np.ndarray | None = None,
    ) -> float:
        if not b > a:
            return 0.0
        nodes, weights = _legendre_rule(self.moment_quadrature_order)
        half = 0.5 * (b - a)
        mid = 0.5 * (a + b)
        x = mid + half * nodes
        data = self.joined.moments.pa16_input_at_eta(float(eta))
        vals = self._selected_join().profile_values_scaled(
            x,
            eta=float(eta),
            ell_i=float(data["ell_i"]),
            G_i=float(data["G_i"]),
            coefficients=coefficients,
        )
        U = np.asarray(vals["U"], dtype=float)
        return float(half * np.dot(weights, U))

    def _outer_scaled_U_integral(self, x: float, eta: float) -> float:
        """Return ``int_{x_i}^x U(s,eta) ds`` without integrating a huge X interval."""
        x_i = self.X_i / self.X_R
        if x < x_i - 32.0 * np.finfo(float).eps * max(x_i, 1.0e-300):
            raise ValueError("scaled x lies below x_i")
        x = max(x, x_i)
        x_restore0 = math.exp(LOG_RESTORE_START)
        x_restore1 = math.exp(LOG_RESTORE_END)
        x_repair0 = math.exp(LOG_REPAIR_START)
        x_exit = math.exp(LOG_JOIN_EXIT)
        if x > x_exit * (1.0 + 64.0 * np.finfo(float).eps):
            raise ValueError("scaled x lies beyond the current X_h exit")

        data = self.joined.moments.pa16_input_at_eta(float(eta))
        G_i = float(data["G_i"])
        total = 0.0

        first_stop = min(x, x_restore0)
        if first_stop > x_i:
            total += G_i * (first_stop - x_i)
        if x <= x_restore0:
            return float(total)

        restore_stop = min(x, x_restore1)
        total += self._integrate_scaled_U(x_restore0, restore_stop, eta=eta)
        if x <= x_restore1:
            return float(total)

        ideal_stop = min(x, x_repair0)
        if ideal_stop > x_restore1:
            total += 4.0 * eta * (ideal_stop - x_restore1)
        if x <= x_repair0:
            return float(total)

        repair_stop = min(x, x_exit)
        solved = self.joined.certificate.solve_selected_at_eta(
            float(eta),
            quadrature_points=self.joined.quadrature_points,
            absolute_tolerance=self.joined.absolute_tolerance,
        )
        if not solved.repair.success:
            raise RuntimeError("existing PA.16 solve did not return success")
        coefficients = np.asarray(solved.repair.coefficients, dtype=float)
        if coefficients.shape != (5,) or np.any(~np.isfinite(coefficients)):
            raise RuntimeError("existing PA.16 solve returned invalid coefficients")
        total += self._integrate_scaled_U(
            x_repair0,
            repair_stop,
            eta=eta,
            coefficients=coefficients,
        )
        return float(total)

    def _M_over_X_scalar(self, X: float, eta: float) -> float:
        X = float(X)
        eta = float(eta)
        if not math.isfinite(X) or not 0.0 <= X <= self.X_h:
            raise ValueError("X must lie in [0,X_h]")
        lo, hi = self.eta_interval
        if not math.isfinite(eta) or not lo <= eta <= hi:
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        if X == 0.0:
            vals = self.geometry.source_center.values(0.0, eta)
            return float(np.asarray(vals["M_0_over_X"]))
        if X <= self.X_i:
            return self._prefix_M(X, eta) / X

        x = X / self.X_R
        prefix_M_scaled = (
            float(self.joined.moments.physical_prefix_moments_at_Xi(eta)[0])
            / self.X_R
        )
        scaled_integral = self._outer_scaled_U_integral(x, eta)
        out = (prefix_M_scaled + scaled_integral) / x
        if not math.isfinite(out):
            raise RuntimeError("M/X continuation became non-finite")
        return float(out)

    def _M_ratio_and_eta_derivative(self, X: float, eta: float) -> tuple[float, float]:
        """Return ``M/X`` and ``partial_eta(M/X)`` using frozen FD4 stencils."""
        lo, hi = self.eta_interval
        h = self.eta_fd_step
        cache: dict[float, float] = {}

        def f(value: float) -> float:
            key = float(value)
            if key not in cache:
                cache[key] = self._M_over_X_scalar(float(X), key)
            return cache[key]

        center = f(float(eta))
        if eta - 2.0 * h >= lo and eta + 2.0 * h <= hi:
            derivative = (
                f(eta - 2.0 * h)
                - 8.0 * f(eta - h)
                + 8.0 * f(eta + h)
                - f(eta + 2.0 * h)
            ) / (12.0 * h)
        elif eta + 4.0 * h <= hi:
            derivative = (
                -25.0 * center
                + 48.0 * f(eta + h)
                - 36.0 * f(eta + 2.0 * h)
                + 16.0 * f(eta + 3.0 * h)
                - 3.0 * f(eta + 4.0 * h)
            ) / (12.0 * h)
        elif eta - 4.0 * h >= lo:
            derivative = (
                25.0 * center
                - 48.0 * f(eta - h)
                + 36.0 * f(eta - 2.0 * h)
                - 16.0 * f(eta - 3.0 * h)
                + 3.0 * f(eta - 4.0 * h)
            ) / (12.0 * h)
        else:
            raise ValueError("eta interval is too narrow for the frozen FD4 stencil")
        if not math.isfinite(derivative):
            raise RuntimeError("M_eta/X finite difference became non-finite")
        return center, float(derivative)

    def similarity_profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return current ``F,U,E,M/X,M_eta/X,v0`` on the joined source domain."""
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_h)):
            raise ValueError("X must lie in [0,X_h]")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")

        joined = self.joined.profile_values(X_arr, eta_arr)
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        m_ratio = np.empty_like(xf)
        m_eta_ratio = np.empty_like(xf)
        for index, (xv, ev) in enumerate(zip(xf, ef, strict=True)):
            m_ratio[index], m_eta_ratio[index] = self._M_ratio_and_eta_derivative(
                float(xv), float(ev)
            )

        # Reuse the already-governed source axis data for L,d rather than
        # re-encoding their eta polynomials here.  Y=0 is inside the center
        # profile domain and L,d are the source angular coefficients.
        axis = self.geometry.physical_profiles.axis_profiles.values(
            np.zeros_like(eta_arr), eta_arr
        )
        U = np.asarray(joined["U_current_joined"], dtype=float)
        numerator = (
            2.0 * eta_arr * U
            - 2.0 * self.D * eta_arr * m_ratio.reshape(shape)
            - np.asarray(axis["d"], dtype=float) * m_eta_ratio.reshape(shape)
        )
        v0 = numerator / np.asarray(axis["L"], dtype=float)
        if np.any(~np.isfinite(v0)):
            raise RuntimeError("current joined radial profile v0 became non-finite")
        return {
            "X": X_arr,
            "eta": eta_arr,
            "F_current_joined": np.asarray(joined["F_current_joined"], dtype=float),
            "U_current_joined": U,
            "E_current_joined": np.asarray(joined["E_current_joined"], dtype=float),
            "M_over_X_current_joined": m_ratio.reshape(shape),
            "M_eta_over_X_current_joined": m_eta_ratio.reshape(shape),
            "v0_current_joined": v0,
        }

    def values(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Return source coordinates, current profiles, and Cartesian leading velocity."""
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        radial_tol = 128.0 * np.finfo(float).eps * max(1.0, self.X_h)
        if np.any(X > self.X_h + radial_tol):
            raise ValueError(
                "point lies beyond current X_h; outer/global leading completion is not materialized"
            )
        X_eval = np.minimum(X, self.X_h)
        profile = self.similarity_profile_values(X_eval, coords["eta"])
        q = np.asarray(coords["q"], dtype=float)
        q_radial = profile["v0_current_joined"] / (2.0 * q)
        q_swirl = np.power(q, -self.A - 0.5) * profile["F_current_joined"]
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = np.power(q, -self.A) * profile["U_current_joined"]
        if np.any(~np.isfinite(u)) or np.any(~np.isfinite(v)) or np.any(~np.isfinite(w)):
            raise RuntimeError("current Cartesian leading velocity became non-finite")
        return {
            **coords,
            **profile,
            "u": u,
            "v": v,
            "w": w,
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Unified vectorized ``velocity(x,y,z,t)->[...,3]`` interface."""
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "joined": self.joined.configuration(),
            "geometry": self.geometry.field_configuration(),
            "moment_quadrature_order": self.moment_quadrature_order,
            "eta_fd_step": self.eta_fd_step,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianLeadingVelocity":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current Cartesian-leading schema")
        return cls(
            joined=KokunoPA16CurrentJoinedProfile.from_configuration(payload["joined"]),
            geometry=KokunoPA10CartesianCenterVelocity.from_configuration(payload["geometry"]),
            moment_quadrature_order=int(payload["moment_quadrature_order"]),
            eta_fd_step=float(payload["eta_fd_step"]),
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianLeadingVelocity":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        natural_X = 0.5 * float(self.joined.moments.X_0)
        eta = 0.25
        current = self.similarity_profile_values(natural_X, eta)
        old = self.geometry.source_center.values(natural_X, eta)
        axis_point = self.geometry.cartesian_from_similarity(0.0, eta, 0.5)
        axis_velocity = self.velocity(
            axis_point["x"], axis_point["y"], axis_point["z"], axis_point["t"]
        )
        return {
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "numerical_realization": copy.deepcopy(_NUMERICAL_REALIZATION),
            "configuration": self.configuration(),
            "parent_joined_semantic_sha256": self.joined.semantic_sha256,
            "governed_geometry_field_sha256": self.geometry.field_sha256,
            "machine_checks": {
                "natural_lane_v0_abs_replay_error": abs(
                    float(np.asarray(current["v0_current_joined"]))
                    - float(np.asarray(old["v_0"]))
                ),
                "axis_transverse_velocity_exact_zero": bool(
                    np.all(np.asarray(axis_velocity)[..., :2] == 0.0)
                ),
                "axis_velocity_finite": bool(np.all(np.isfinite(axis_velocity))),
                "velocity_interface_nontrivial": bool(np.any(np.abs(axis_velocity) > 0.0)),
            },
            "domain": {
                "X_i": self.X_i,
                "X_h": self.X_h,
                "eta_interval": list(self.eta_interval),
                "time_interval": list(self.geometry.time_interval),
            },
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    field = KokunoPA16CurrentCartesianLeadingVelocity()
    payload = field.report()
    text = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    if args.config_output is not None:
        field.save_configuration(args.config_output)


if __name__ == "__main__":
    _main()
