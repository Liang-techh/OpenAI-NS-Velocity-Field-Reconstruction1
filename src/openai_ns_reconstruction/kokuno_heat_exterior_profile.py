"""Executable Kokuno heat-exterior leading-profile component.

This module implements the source heat factor from the pinned public reconstruction
and maps it through Kokuno's native ``q/X/eta`` coordinates.  It is deliberately a
component of the global leading field rather than a core-to-exterior matching rule:
the source outer construction and compensated splice remain separate dependencies.

For ``Z >= 0`` the pinned source defines

    H(Z) = Gamma(1+h)^(-1) int_0^inf exp(-v) v^h (1+Z v)^(-h) dv,

with ``A=1/2+h`` and

    E_heat(X,eta) = c_inf X^(-A) H(2(1-eta^2)/X),
    F_heat = E_heat / sqrt(2X),  U_heat = v0_heat = 0.

The resulting Cartesian velocity is the exact pure-swirl exterior heat field on
``X>0``.  The amplitude ``c_inf`` is explicit here because the source outer assembly,
which determines its final matched value, is not yet reconstructed in this lane.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import roots_genlaguerre

from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-heat-exterior-profile-v1"

_SOURCE_FORMULAS = {
    "heat_factor": (
        "H(Z)=Gamma(1+h)^(-1)*int_0^inf exp(-v)*v^h*(1+Z*v)^(-h) dv, Z>=0"
    ),
    "heat_derivative": (
        "H^(m)(Z)=(-1)^m*(h)_m/Gamma(1+h)*int_0^inf "
        "exp(-v)*v^(h+m)*(1+Z*v)^(-h-m) dv"
    ),
    "heat_ode": "Z^2*H_ZZ+(1+2*(1+h)*Z)*H_Z+h*(1+h)*H=0",
    "similarity_argument": "Z=2*(1-eta^2)/X=2*(1-t)/s",
    "E_heat": "c_inf*X^(-A)*H(Z), A=1/2+h",
    "F_heat": "E_heat/sqrt(2*X)",
    "exterior_components": "U_heat=0, v0_heat=0",
    "physical_swirl": "K(r,t)=c_inf*s^(-A)*H(2*(1-t)/s), s=r^2/2",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "heat_factor_formula_executable": True,
    "heat_exterior_velocity_component_executable": True,
    "source_outer_amplitude_recovered": False,
    "core_to_heat_matching_completed": False,
    "outer_compensation_completed": False,
    "axis_regular_complete_field": False,
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


@dataclass(frozen=True)
class KokunoHeatExteriorProfile:
    """Vectorized source heat factor and pure-swirl exterior velocity component.

    ``c_inf`` is an explicit autonomous amplitude placeholder pending the source outer
    assembly.  This class does not splice itself to the core series and therefore
    rejects the axis (``X=0``), where the exterior heat component alone is singular.
    """

    h: float = 0.005
    c_inf: float = 1.0
    quadrature_order: int = 96

    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )
    _laguerre_nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _laguerre_weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        h = float(self.h)
        c_inf = float(self.c_inf)
        if isinstance(self.quadrature_order, bool) or not isinstance(
            self.quadrature_order, (int, np.integer)
        ):
            raise TypeError("quadrature_order must be an integer")
        order = int(self.quadrature_order)
        if not np.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must be finite and satisfy 0 < h < 1e-2")
        if not np.isfinite(c_inf) or not (0.0 < c_inf <= 1.0e3):
            raise ValueError("c_inf must be finite and satisfy 0 < c_inf <= 1e3")
        if order < 16 or order > 256:
            raise ValueError("quadrature_order must lie in [16, 256]")

        nodes, weights = roots_genlaguerre(order, h)
        if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(weights)):
            raise ArithmeticError("generalized Gauss-Laguerre rule is non-finite")
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)

        object.__setattr__(self, "h", h)
        object.__setattr__(self, "c_inf", c_inf)
        object.__setattr__(self, "quadrature_order", order)
        object.__setattr__(self, "_coordinates", KokunoNativeSimilarityCoordinates(h=h))
        object.__setattr__(self, "_laguerre_nodes", nodes)
        object.__setattr__(self, "_laguerre_weights", weights)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @staticmethod
    def _rising(value: float, order: int) -> float:
        result = 1.0
        for k in range(order):
            result *= value + k
        return result

    def heat_factor(self, Z: Any, derivative: int = 0) -> np.ndarray:
        """Evaluate the source heat factor or one of its first four Z derivatives."""

        if isinstance(derivative, bool) or not isinstance(derivative, (int, np.integer)):
            raise TypeError("derivative must be an integer")
        derivative = int(derivative)
        if derivative < 0 or derivative > 4:
            raise ValueError("derivative must lie in [0, 4]")
        z = _finite_array(Z, "Z")
        if np.any(z < 0.0):
            raise ValueError("source heat factor requires Z >= 0")

        flat = z.reshape(-1)
        nodes = self._laguerre_nodes[:, None]
        kernel = np.power(nodes, derivative) * np.power(
            1.0 + nodes * flat[None, :], -self.h - derivative
        )
        factor = ((-1.0) ** derivative) * self._rising(self.h, derivative)
        values = factor * (self._laguerre_weights @ kernel) / math.gamma(1.0 + self.h)
        return np.asarray(values, dtype=float).reshape(z.shape)

    def heat_ode_residual(self, Z: Any) -> np.ndarray:
        z = _finite_array(Z, "Z")
        H = self.heat_factor(z, 0)
        H_Z = self.heat_factor(z, 1)
        H_ZZ = self.heat_factor(z, 2)
        return (
            z * z * H_ZZ
            + (1.0 + 2.0 * (1.0 + self.h) * z) * H_Z
            + self.h * (1.0 + self.h) * H
        )

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return the source exterior ``E,F,U,v0`` profile and first partials."""

        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array <= 0.0):
            raise ValueError("heat-exterior profile requires X > 0")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("heat-exterior profile requires |eta| <= 1")

        d = 1.0 - eta_array * eta_array
        Z = 2.0 * d / X_array
        H = self.heat_factor(Z, 0)
        H_Z = self.heat_factor(Z, 1)
        H_ZZ = self.heat_factor(Z, 2)
        Z_X = -Z / X_array
        Z_eta = -4.0 * eta_array / X_array

        E = self.c_inf * np.power(X_array, -self.A) * H
        E_X = self.c_inf * np.power(X_array, -self.A - 1.0) * (
            -self.A * H - Z * H_Z
        )
        E_eta = self.c_inf * np.power(X_array, -self.A) * H_Z * Z_eta

        B = 1.0 + self.h
        F = (self.c_inf / math.sqrt(2.0)) * np.power(X_array, -B) * H
        F_X = (self.c_inf / math.sqrt(2.0)) * np.power(X_array, -B - 1.0) * (
            -B * H - Z * H_Z
        )
        F_eta = (self.c_inf / math.sqrt(2.0)) * np.power(X_array, -B) * H_Z * Z_eta
        zeros = np.zeros_like(F)

        return {
            "Z": Z,
            "H": H,
            "H_Z": H_Z,
            "H_ZZ": H_ZZ,
            "E": E,
            "E_X": E_X,
            "E_eta": E_eta,
            "F": F,
            "F_X": F_X,
            "F_eta": F_eta,
            "U": zeros,
            "U_X": zeros,
            "U_eta": zeros,
            "v0": zeros,
            "Pi_X": F * F,
            "Z_X": Z_X,
            "Z_eta": Z_eta,
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate the pure-swirl heat-exterior component through native coordinates."""

        coordinates = self._coordinates.evaluate(x, y, z, t)
        X = coordinates["X"]
        if np.any(X <= 0.0):
            raise ValueError("heat-exterior velocity is not an axis field; require X > 0")
        profiles = self.profile_values(X, coordinates["eta"])
        factor = np.power(coordinates["q"], -1.0 - self.h) * profiles["F"]
        x_array, y_array = np.broadcast_arrays(_finite_array(x, "x"), _finite_array(y, "y"))
        x_array = np.broadcast_to(x_array, np.shape(factor))
        y_array = np.broadcast_to(y_array, np.shape(factor))
        zeros = np.zeros_like(factor)
        return np.stack((-y_array * factor, x_array * factor, zeros), axis=-1)

    def direct_heat_velocity(self, x: Any, y: Any, t: Any) -> np.ndarray:
        """Evaluate the source physical ``K(r,t)`` formula, independent of q/eta."""

        x_array, y_array, t_array = np.broadcast_arrays(
            _finite_array(x, "x"), _finite_array(y, "y"), _finite_array(t, "t")
        )
        tau = 1.0 - t_array
        if np.any(tau <= 0.0):
            raise ValueError("source heat field requires t < 1")
        s = 0.5 * (x_array * x_array + y_array * y_array)
        if np.any(s <= 0.0):
            raise ValueError("source heat exterior formula requires r > 0")
        Z = 2.0 * tau / s
        K = self.c_inf * np.power(s, -self.A) * self.heat_factor(Z)
        r = np.sqrt(2.0 * s)
        zeros = np.zeros_like(K)
        return np.stack((-y_array * K / r, x_array * K / r, zeros), axis=-1)

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = _finite_array(points_xyz, "points_xyz")
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points_xyz must have final dimension 3")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "heat factor and exact exterior pure-swirl component",
            },
            "parameters": {
                "h": self.h,
                "h_origin": "autonomous_choice_within_public_source_bound_0<h<1e-2",
                "c_inf": self.c_inf,
                "c_inf_origin": "autonomous_placeholder_pending_source_outer_assembly",
                "quadrature_order": self.quadrature_order,
                "quadrature": "generalized_gauss_laguerre_alpha_h",
            },
            "formulas": dict(_SOURCE_FORMULAS),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoHeatExteriorProfile":
        if not isinstance(payload, dict):
            raise ValueError("heat-exterior payload must be a JSON object")
        expected = {"schema", "source", "parameters", "formulas", "truth_boundary"}
        if set(payload) - {"sha256"} != expected:
            raise ValueError("heat-exterior payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported Kokuno heat-exterior schema")

        default = cls().to_payload()
        if payload["source"] != default["source"]:
            raise ValueError("Kokuno heat source/provenance metadata changed")
        if payload["formulas"] != _SOURCE_FORMULAS:
            raise ValueError("Kokuno heat formula metadata changed")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("Kokuno heat truth-boundary metadata changed")

        parameters = payload["parameters"]
        expected_parameter_keys = {
            "h",
            "h_origin",
            "c_inf",
            "c_inf_origin",
            "quadrature_order",
            "quadrature",
        }
        if set(parameters) != expected_parameter_keys:
            raise ValueError("Kokuno heat parameter metadata changed")
        if parameters["h_origin"] != default["parameters"]["h_origin"]:
            raise ValueError("Kokuno heat h-origin metadata changed")
        if parameters["c_inf_origin"] != default["parameters"]["c_inf_origin"]:
            raise ValueError("Kokuno heat amplitude provenance changed")
        if parameters["quadrature"] != default["parameters"]["quadrature"]:
            raise ValueError("Kokuno heat quadrature rule changed")

        result = cls(
            h=float(parameters["h"]),
            c_inf=float(parameters["c_inf"]),
            quadrature_order=int(parameters["quadrature_order"]),
        )
        if "sha256" in payload and payload["sha256"] != result.sha256:
            raise ValueError("Kokuno heat-exterior payload SHA mismatch")
        return result

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoHeatExteriorProfile":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
