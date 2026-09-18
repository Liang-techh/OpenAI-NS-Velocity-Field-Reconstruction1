"""Executable first outer continuation of the Kokuno leading profile.

This module implements exactly one public reconstruction stage: the reference
continuation in ``continuation_exact.md`` / the corrected 2026-09-09 workbench.
It extends the finite nonlinear core profile through a short logarithmic-X
transition and then freezes the reference ``phi`` and ``U`` profiles.  Later
cone modulation, five-moment repair, pulse/release, heat compensation and the
final heat tail are deliberately *not* implemented here.

The resulting field is therefore a callable source-stage candidate useful for
closing the core->outer representation seam.  It is not the final Kokuno
leading field and is not an independently validated Navier--Stokes solution.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_leading_core_series import (
    CORE_Y_MAX,
    KokunoLeadingCoreSeriesCandidate,
)
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-reference-continuation-v1"

_SOURCE_MAP = {
    "source_stage": "continuation_exact/reference continuation before cone and moment repairs",
    "X0": "X0=4/Lambda",
    "log_coordinate": "y=log(X/X0)",
    "smooth_step": (
        "sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)); "
        "sigma=0 for s<=0, sigma=1 for s>=1"
    ),
    "transition_phi": (
        "D_y log(phi_r)=(1-sigma((y-t1)/t1))*D_X log(phi_nat), "
        "D_X=X*d/dX"
    ),
    "transition_U": (
        "D_y U_r=(1-sigma((y-t1)/t1))*D_X U_nat"
    ),
    "post_transition": "phi_r and U_r are constant in y for y>=2*t1",
    "swirl_map": "F=phi/C; E=sqrt(2X)*F",
    "pressure_identity": "Pi_X=F^2",
    "native_coordinates": "q-z^2*q^(2h)=1-t; X=(x^2+y^2)/(2q)",
    "cartesian_velocity": (
        "u1=x*v0/(2q)-y*q^(-1-h)*F; "
        "u2=y*v0/(2q)+x*q^(-1-h)*F; u3=q^(-1/2-h)*U"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "reference_continuation_executable": True,
    "velocity_api_compatible": True,
    "core_domain_only": False,
    "reference_profile_global_in_X": True,
    "cone_modulation_completed": False,
    "five_moment_repair_completed": False,
    "heat_compensation_completed": False,
    "core_to_heat_matching_completed": False,
    "global_leading_profile_reconstructed": False,
    "complete_kokuno_composite_velocity": False,
    "independent_full_momentum_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_array(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def source_smooth_step(s: Any) -> np.ndarray:
    """Source C-infinity step, evaluated without endpoint overflow."""

    values = _finite_array(s, "s")
    out = np.empty_like(values, dtype=float)
    out[values <= 0.0] = 0.0
    out[values >= 1.0] = 1.0
    mask = (values > 0.0) & (values < 1.0)
    if np.any(mask):
        x = values[mask]
        log_a = -1.0 / (x * x)
        log_b = -1.0 / ((1.0 - x) * (1.0 - x))
        delta = log_b - log_a
        local = np.empty_like(delta)
        high = delta > 50.0
        low = delta < -50.0
        mid = ~(high | low)
        local[high] = np.exp(-delta[high])
        local[low] = 1.0 - np.exp(delta[low])
        local[mid] = 1.0 / (1.0 + np.exp(delta[mid]))
        out[mask] = local
    return out


@dataclass(frozen=True)
class KokunoReferenceContinuationCandidate:
    """Source-stage reference continuation of the finite leading core.

    ``log_transition_width`` is the source's ``t1``.  It is an autonomous
    finite reconstruction choice constrained by the source requirement
    ``4*exp(2*t1)<4.1`` so every natural-profile sample needed by this stage
    remains inside the already executable core series.
    """

    h: float = 0.005
    j0: float = 0.02
    sigma: float = 0.5
    Lambda: float = 10.0
    C: float = 2.0
    pressure_scale: float = 1.0
    maxdegree: int = 14
    eta_nodes: int = 257
    quadrature_points: int = 24
    log_transition_width: float = 0.005

    _core: KokunoLeadingCoreSeriesCandidate = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(init=False, repr=False, compare=False)
    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        width = float(self.log_transition_width)
        max_width = 0.5 * math.log(CORE_Y_MAX / 4.0)
        if not math.isfinite(width) or not (0.0 < width < max_width):
            raise ValueError(
                "log_transition_width must satisfy 0<t1<0.5*log(4.1/4)"
            )
        if isinstance(self.quadrature_points, bool) or not isinstance(
            self.quadrature_points, (int, np.integer)
        ):
            raise TypeError("quadrature_points must be an integer")
        order = int(self.quadrature_points)
        if not 8 <= order <= 64:
            raise ValueError("quadrature_points must lie in [8,64]")

        core = KokunoLeadingCoreSeriesCandidate(
            h=float(self.h),
            j0=float(self.j0),
            sigma=float(self.sigma),
            Lambda=float(self.Lambda),
            C=float(self.C),
            pressure_scale=float(self.pressure_scale),
            maxdegree=int(self.maxdegree),
            eta_nodes=int(self.eta_nodes),
            quadrature_points=max(8, min(128, order)),
        )
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)

        object.__setattr__(self, "h", core.h)
        object.__setattr__(self, "j0", core.j0)
        object.__setattr__(self, "sigma", core.sigma)
        object.__setattr__(self, "Lambda", core.Lambda)
        object.__setattr__(self, "C", core.C)
        object.__setattr__(self, "pressure_scale", core.pressure_scale)
        object.__setattr__(self, "maxdegree", core.maxdegree)
        object.__setattr__(self, "eta_nodes", core.eta_nodes)
        object.__setattr__(self, "quadrature_points", order)
        object.__setattr__(self, "log_transition_width", width)
        object.__setattr__(self, "_core", core)
        object.__setattr__(self, "_coordinates", core.coordinates_model)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)

    @property
    def core(self) -> KokunoLeadingCoreSeriesCandidate:
        return self._core

    @property
    def X0(self) -> float:
        return 4.0 / self.Lambda

    @property
    def X1(self) -> float:
        return self.X0 * math.exp(self.log_transition_width)

    @property
    def X2(self) -> float:
        return self.X0 * math.exp(2.0 * self.log_transition_width)

    @property
    def max_source_transition_X(self) -> float:
        return CORE_Y_MAX / self.Lambda

    @staticmethod
    def _time_array(t: Any) -> np.ndarray:
        values = _finite_array(t, "t")
        if np.any((values < 0.0) | (values >= 1.0)):
            raise ValueError("time must satisfy 0 <= t < 1")
        return values

    def _integrate_y(self, lower: float, upper: float, values_fn) -> np.ndarray:
        if upper <= lower:
            return np.asarray(values_fn(np.empty(0, dtype=float)))[:0].sum(axis=0)
        center = 0.5 * (lower + upper)
        half = 0.5 * (upper - lower)
        y = center + half * self._nodes
        values = np.asarray(values_fn(y), dtype=float)
        if values.shape[0] != self._nodes.size:
            raise RuntimeError("continuation integrand returned an invalid quadrature shape")
        return half * np.tensordot(self._weights, values, axes=(0, 0))

    def _continued_at_X(self, X: float, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        """Return source-stage ``F,U`` at scalar X and arbitrary eta array."""

        X = float(X)
        if not math.isfinite(X) or X < 0.0:
            raise ValueError("X must be finite and nonnegative")
        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0 + 2.0e-13):
            raise ValueError("eta must satisfy |eta|<=1")
        eta_array = np.clip(eta_array, -1.0, 1.0)
        series = self._core.series
        if X <= self.X1:
            return (
                np.asarray(series.F(X, eta_array), dtype=float),
                np.asarray(series.U(X, eta_array), dtype=float),
            )

        y_upper = min(math.log(X / self.X0), 2.0 * self.log_transition_width)
        eta_flat = eta_array.reshape(-1)
        F1 = np.asarray(series.F(self.X1, eta_flat), dtype=float)
        U1 = np.asarray(series.U(self.X1, eta_flat), dtype=float)

        def logF_integrand(y: np.ndarray) -> np.ndarray:
            if y.size == 0:
                return np.empty((0, eta_flat.size), dtype=float)
            xq = self.X0 * np.exp(y)
            Xq = xq[:, None]
            Eq = eta_flat[None, :]
            Fnat = np.asarray(series.F(Xq, Eq), dtype=float)
            if np.any(Fnat <= 0.0):
                raise RuntimeError("source continuation requires positive natural phi/F")
            F_X = np.asarray(series.F_radial_derivative(Xq, Eq), dtype=float)
            s = (y - self.log_transition_width) / self.log_transition_width
            gate = 1.0 - source_smooth_step(s)
            return gate[:, None] * Xq * F_X / Fnat

        def U_integrand(y: np.ndarray) -> np.ndarray:
            if y.size == 0:
                return np.empty((0, eta_flat.size), dtype=float)
            xq = self.X0 * np.exp(y)
            Xq = xq[:, None]
            Eq = eta_flat[None, :]
            U_X = np.asarray(series.U_radial_derivative(Xq, Eq), dtype=float)
            s = (y - self.log_transition_width) / self.log_transition_width
            gate = 1.0 - source_smooth_step(s)
            return gate[:, None] * Xq * U_X

        log_increment = self._integrate_y(
            self.log_transition_width, y_upper, logF_integrand
        )
        U_increment = self._integrate_y(
            self.log_transition_width, y_upper, U_integrand
        )
        F = F1 * np.exp(log_increment)
        U = U1 + U_increment
        shape = eta_array.shape
        return F.reshape(shape), U.reshape(shape)

    def _broadcast_profile(self, X: Any, eta: Any, component: str) -> np.ndarray:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array < 0.0) or np.any(np.abs(eta_array) > 1.0 + 2.0e-13):
            raise ValueError("profile domain requires X>=0 and |eta|<=1")
        result = np.empty_like(X_array, dtype=float)
        flat = result.reshape(-1)
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            pair = self._continued_at_X(float(xx), float(ee))
            flat[index] = float(pair[0 if component == "F" else 1])
        return result

    def F(self, X: Any, eta: Any) -> np.ndarray:
        return self._broadcast_profile(X, eta, "F")

    def U(self, X: Any, eta: Any) -> np.ndarray:
        return self._broadcast_profile(X, eta, "U")

    def _eta_derivative(self, X: Any, eta: Any, component: str) -> np.ndarray:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array < 0.0) or np.any(np.abs(eta_array) > 1.0 + 2.0e-13):
            raise ValueError("profile domain requires X>=0 and |eta|<=1")
        grid = self._core.series.grid
        result = np.empty_like(X_array, dtype=float)
        flat = result.reshape(-1)
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            F_nodes, U_nodes = self._continued_at_X(float(xx), grid.eta)
            nodes = F_nodes if component == "F" else U_nodes
            flat[index] = float(grid.differentiate_and_interpolate(nodes, float(ee)))
        return result

    def F_eta(self, X: Any, eta: Any) -> np.ndarray:
        return self._eta_derivative(X, eta, "F")

    def U_eta(self, X: Any, eta: Any) -> np.ndarray:
        return self._eta_derivative(X, eta, "U")

    def F_X(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array < 0.0) or np.any(np.abs(eta_array) > 1.0 + 2.0e-13):
            raise ValueError("profile domain requires X>=0 and |eta|<=1")
        out = np.empty_like(X_array, dtype=float)
        series = self._core.series
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            x = float(xx)
            e = float(ee)
            if x <= self.X1:
                value = float(series.F_radial_derivative(x, e))
            elif x >= self.X2:
                value = 0.0
            else:
                F = float(self._continued_at_X(x, e)[0])
                Fnat = float(series.F(x, e))
                gate = 1.0 - float(
                    source_smooth_step(
                        (math.log(x / self.X0) - self.log_transition_width)
                        / self.log_transition_width
                    )
                )
                value = F * gate * float(series.F_radial_derivative(x, e)) / Fnat
            out.reshape(-1)[index] = value
        return out

    def U_X(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array < 0.0) or np.any(np.abs(eta_array) > 1.0 + 2.0e-13):
            raise ValueError("profile domain requires X>=0 and |eta|<=1")
        out = np.empty_like(X_array, dtype=float)
        series = self._core.series
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            x = float(xx)
            e = float(ee)
            if x <= self.X1:
                value = float(series.U_radial_derivative(x, e))
            elif x >= self.X2:
                value = 0.0
            else:
                gate = 1.0 - float(
                    source_smooth_step(
                        (math.log(x / self.X0) - self.log_transition_width)
                        / self.log_transition_width
                    )
                )
                value = gate * float(series.U_radial_derivative(x, e))
            out.reshape(-1)[index] = value
        return out

    def _average_U_and_eta(self, X: float, eta: float) -> tuple[float, float]:
        if X == 0.0:
            return float(self.U(0.0, eta)), float(self.U_eta(0.0, eta))
        # Average in X by Gauss-Legendre quadrature.  Differentiate the averaged
        # eta profile spectrally, matching the finite core's eta convention.
        radial_nodes = 0.5 * (self._nodes + 1.0) * X
        radial_weights = 0.5 * self._weights
        grid = self._core.series.grid
        sampled = np.stack(
            [self._continued_at_X(float(xx), grid.eta)[1] for xx in radial_nodes],
            axis=0,
        )
        average_nodes = np.tensordot(radial_weights, sampled, axes=(0, 0))
        avg = float(grid.interpolate(average_nodes, eta))
        avg_eta = float(grid.differentiate_and_interpolate(average_nodes, eta))
        return avg, avg_eta

    def v0(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        out = np.empty_like(X_array, dtype=float)
        D = 0.5 - self.h
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            x = float(xx)
            e = float(ee)
            U = float(self.U(x, e))
            average, average_eta = self._average_U_and_eta(x, e)
            d = 1.0 - e * e
            L = 1.0 - 2.0 * self.h * e * e
            out.reshape(-1)[index] = (
                2.0 * e * U - 2.0 * D * e * average - d * average_eta
            ) / L
        return out

    def _integrate_F2(self, lower: float, upper: float, eta: float) -> float:
        if upper <= lower:
            return 0.0
        center = 0.5 * (lower + upper)
        half = 0.5 * (upper - lower)
        points = center + half * self._nodes
        values = np.array([float(self.F(float(x), eta)) ** 2 for x in points])
        return float(half * (self._weights @ values))

    def Pi(self, X: Any, eta: Any) -> np.ndarray:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array < 0.0) or np.any(np.abs(eta_array) > 1.0 + 2.0e-13):
            raise ValueError("profile domain requires X>=0 and |eta|<=1")
        out = np.empty_like(X_array, dtype=float)
        series = self._core.series
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            x = float(xx)
            e = float(ee)
            if x <= self.X1:
                value = float(series.Pi(x, e))
            else:
                pi1 = float(series.Pi(self.X1, e))
                transition_end = min(x, self.X2)
                value = pi1 + self._integrate_F2(self.X1, transition_end, e)
                if x > self.X2:
                    F2 = float(self.F(self.X2, e)) ** 2
                    value += (x - self.X2) * F2
            out.reshape(-1)[index] = value
        return out

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        F = self.F(X, eta)
        return {
            "F": F,
            "F_X": self.F_X(X, eta),
            "F_eta": self.F_eta(X, eta),
            "U": self.U(X, eta),
            "U_X": self.U_X(X, eta),
            "U_eta": self.U_eta(X, eta),
            "v0": self.v0(X, eta),
            "Pi": self.Pi(X, eta),
            "Pi_X": F * F,
        }

    def coordinates(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        return self._coordinates.evaluate(x, y, z, t)

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            self._time_array(t),
        )
        coords = self.coordinates(x_array, y_array, z_array, t_array)
        q = coords["q"]
        X = coords["X"]
        eta = coords["eta"]
        F = self.F(X, eta)
        U = self.U(X, eta)
        v0 = self.v0(X, eta)
        A = 0.5 + self.h
        result = np.stack(
            (
                x_array * v0 / (2.0 * q) - y_array * q ** (-1.0 - self.h) * F,
                y_array * v0 / (2.0 * q) + x_array * q ** (-1.0 - self.h) * F,
                q ** (-A) * U,
            ),
            axis=-1,
        )
        if not np.all(np.isfinite(result)):
            raise RuntimeError("Kokuno reference-continuation velocity became non-finite")
        return result

    __call__ = velocity

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = _finite_array(points_xyz, "points_xyz")
        if points.ndim < 2 or points.shape[-1] != 3:
            raise ValueError("points_xyz must have shape (...,3)")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def pressure(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            self._time_array(t),
        )
        coords = self.coordinates(x_array, y_array, z_array, t_array)
        A = 0.5 + self.h
        result = coords["q"] ** (-2.0 * A) * self.Pi(coords["X"], coords["eta"])
        if not np.all(np.isfinite(result)):
            raise RuntimeError("Kokuno reference-continuation pressure became non-finite")
        return result

    def grid(self, x: Any, y: Any, z: Any, times: Any) -> np.ndarray:
        x_axis = _finite_array(x, "x")
        y_axis = _finite_array(y, "y")
        z_axis = _finite_array(z, "z")
        time_axis = self._time_array(times)
        for name, axis in (("x", x_axis), ("y", y_axis), ("z", z_axis), ("times", time_axis)):
            if axis.ndim != 1 or axis.size == 0:
                raise ValueError(f"{name} must be a nonempty one-dimensional axis")
        Xg, Yg, Zg = np.meshgrid(x_axis, y_axis, z_axis, indexing="ij")
        output = np.empty(
            (time_axis.size, x_axis.size, y_axis.size, z_axis.size, 3), dtype=float
        )
        for index, time in enumerate(time_axis):
            output[index] = self.velocity(Xg, Yg, Zg, float(time))
        return output

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_evidence": "pinned_public_workbench_commit",
            },
            "parameters": {
                "h": self.h,
                "j0": self.j0,
                "sigma": self.sigma,
                "Lambda": self.Lambda,
                "C": self.C,
                "pressure_scale": self.pressure_scale,
                "maxdegree": self.maxdegree,
                "eta_nodes": self.eta_nodes,
                "quadrature_points": self.quadrature_points,
                "log_transition_width": self.log_transition_width,
                "origin": "autonomous_finite_reconstruction_choices_not_hidden_parameters",
            },
            "transition": {
                "X0": self.X0,
                "X1": self.X1,
                "X2": self.X2,
                "source_core_limit": self.max_source_transition_X,
                "source_strict_bound_satisfied": self.X2 < self.max_source_transition_X,
            },
            "formula_map": dict(_SOURCE_MAP),
            "dependencies": {
                "leading_core_candidate_sha256": self._core.sha256,
                "native_coordinates_sha256": self._coordinates.sha256,
            },
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoReferenceContinuationCandidate":
        if not isinstance(payload, dict):
            raise ValueError("candidate payload must be a JSON object")
        expected = {
            "schema", "source", "parameters", "transition", "formula_map",
            "dependencies", "truth_boundary"
        }
        if set(payload) - {"sha256"} != expected:
            raise ValueError("reference-continuation payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported Kokuno reference-continuation schema")
        candidate = cls(**{
            key: value for key, value in payload["parameters"].items()
            if key != "origin"
        })
        canonical = candidate.to_payload()
        for key in expected - {"parameters"}:
            if payload[key] != canonical[key]:
                raise ValueError(f"reference-continuation {key} metadata changed")
        if payload["parameters"] != canonical["parameters"]:
            raise ValueError("reference-continuation parameter metadata changed")
        if "sha256" in payload and payload["sha256"] != candidate.sha256:
            raise ValueError("reference-continuation SHA mismatch")
        return candidate

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoReferenceContinuationCandidate":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
