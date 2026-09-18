"""Log-stable Appendix-B reference continuation at the selected pressure scale.

This module carries :class:`KokunoSourceRescaledCoreSeed` through the first
public continuation stage in KokunoYumeto's corrected 2026-09-09 reconstruction.
The source uses

    X0 = 4 / Lambda,          y = log(X / X0),

keeps the natural rescaled core through ``y=t1``, and on ``t1<y<2t1``
prescribes

    D_y log(phi_r) = (1-sigma((y-t1)/t1)) D_X log(phi_nat),
    D_y U_r        = (1-sigma((y-t1)/t1)) D_X U_nat,

before freezing ``phi_r`` and ``U_r``.  Since ``F=phi/C`` and ``C`` is
constant in X, the logarithmic equation is identical for ``F``.  The
implementation keeps ``log(F)`` as a primary quantity because the selected
source-compatible pressure datum makes ``Lambda`` enormous and ordinary
binary64 ``F`` legitimately underflows away from the phase maximum.

Only this reference stage is implemented.  The source fixed point/contraction
proof, stress activation, the Xi=110 Appendix-B trajectory, PA.10/PA.16 join,
global pressure matching, and the full leading field remain fail-closed.  The
selected pressure datum and the choice ``Lambda=P^2`` are autonomous
source-compatible numerical choices inherited from ``kokuno_rescaled_core_seed``;
they are not recovered hidden OpenAI/Kokuno parameters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_reference_continuation import source_smooth_step
from .kokuno_rescaled_core_seed import KokunoSourceRescaledCoreSeed
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-rescaled-reference-continuation-v1"

_SOURCE_FORMULAS = {
    "radial_origin": "X0=4/Lambda; y=log(X/X0)",
    "flat_step": (
        "sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), "
        "flat 0/1 extensions"
    ),
    "reference_angular": (
        "D_y log(phi_r)=(1-sigma((y-t1)/t1))*D_X log(phi_nat)"
    ),
    "reference_axial": (
        "D_y U_r=(1-sigma((y-t1)/t1))*D_X U_nat"
    ),
    "post_transition": "phi_r,U_r constant for y>=2*t1",
    "swirl_map": "F=phi/C; E=sqrt(2X)*F",
    "pressure_identity": "Pi_X=F^2",
    "incompressibility": (
        "V0=(2 eta X U-2D eta M-d M_eta)/L; M=int_0^X U dx; v0=V0/X"
    ),
    "cartesian_velocity": (
        "u1=x*v0/(2q)-y*q^(-1-h)*F; "
        "u2=y*v0/(2q)+x*q^(-1-h)*F; u3=q^(-1/2-h)*U"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "selected_source_pressure_datum_carried_into_reference_continuation": True,
    "source_reference_continuation_executable_at_selected_pressure_scale": True,
    "source_exact_radial_transition_equations_executable": True,
    "log_amplitude_representation_executable": True,
    "velocity_api_compatible": True,
    "selected_reference_profile_global_in_X": True,
    "continued_eta_derivatives_are_numerical": True,
    "source_hidden_pressure_scale_recovered": False,
    "source_hidden_rescaling_lambda_recovered": False,
    "source_fixed_point_solved": False,
    "source_contraction_threshold_verified": False,
    "source_complex_C_bound_verified": False,
    "source_appendix_B_activation_executed_at_selected_pressure_scale": False,
    "actual_source_appendix_B_trajectory_reconstructed": False,
    "source_pressure_datum_applied_to_selected_global_path": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_array(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoSourceRescaledReferenceContinuation:
    """Executable selected-pressure realization of the public reference stage."""

    seed: KokunoSourceRescaledCoreSeed = field(
        default_factory=KokunoSourceRescaledCoreSeed
    )
    log_transition_width: float = 0.005
    quadrature_points: int = 20
    eta_derivative_step: float = 2.0e-5

    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.seed, KokunoSourceRescaledCoreSeed):
            raise TypeError("seed must be a KokunoSourceRescaledCoreSeed")
        width = float(self.log_transition_width)
        max_width = 0.5 * math.log(4.1 / 4.0)
        if not math.isfinite(width) or not 0.0 < width < max_width:
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
        eta_step = float(self.eta_derivative_step)
        if not math.isfinite(eta_step) or not 1.0e-7 <= eta_step <= 1.0e-3:
            raise ValueError("eta_derivative_step must lie in [1e-7,1e-3]")
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "log_transition_width", width)
        object.__setattr__(self, "quadrature_points", order)
        object.__setattr__(self, "eta_derivative_step", eta_step)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)
        object.__setattr__(
            self, "_coordinates", KokunoNativeSimilarityCoordinates(h=self.h)
        )

    @property
    def h(self) -> float:
        return float(self.seed.h)

    @property
    def D(self) -> float:
        return 0.5 - self.h

    @property
    def pressure_scale(self) -> float:
        return float(self.seed.pressure_scale)

    @property
    def rescaling_lambda(self) -> float:
        return float(self.seed.rescaling_lambda)

    @property
    def X0(self) -> float:
        return 4.0 / self.rescaling_lambda

    @property
    def X1(self) -> float:
        return self.X0 * math.exp(self.log_transition_width)

    @property
    def X2(self) -> float:
        return self.X0 * math.exp(2.0 * self.log_transition_width)

    @property
    def max_source_transition_X(self) -> float:
        return 4.1 / self.rescaling_lambda

    def _integrate_y(self, lower: float, upper: float, fn: Callable[[np.ndarray], np.ndarray]) -> float:
        if upper <= lower:
            return 0.0
        center = 0.5 * (lower + upper)
        half = 0.5 * (upper - lower)
        y = center + half * self._nodes
        values = np.asarray(fn(y), dtype=float)
        if values.shape != y.shape or not np.all(np.isfinite(values)):
            raise RuntimeError("reference-continuation integrand returned invalid values")
        return float(half * np.dot(self._weights, values))

    def _natural_scalar(self, X: float, eta: float) -> dict[str, float]:
        values = self.seed.profile_values(np.asarray(X), np.asarray(eta))
        return {key: float(np.asarray(value)) for key, value in values.items()}

    def _natural_log_slope(self, X: Any, eta: float) -> np.ndarray:
        """Return D_X log(F_nat) without dividing by underflowed F."""
        X_array = np.asarray(X, dtype=float)
        Y = self.rescaling_lambda * X_array
        state = self.seed.axis_state(np.asarray(eta))
        chi = float(np.asarray(state["chi"]))
        z = Y * chi
        phi = self.seed.f0(z)
        phi_z = self.seed.f0_prime(z)
        if np.any(phi <= 0.0) or np.any(~np.isfinite(phi)):
            raise RuntimeError("natural f0 left the positive source-center range")
        return np.asarray((phi_z / phi) * Y * chi, dtype=float)

    def _natural_DU(self, X: Any, eta: float) -> np.ndarray:
        X_array = np.asarray(X, dtype=float)
        state = self.seed.axis_state(np.asarray(eta))
        B = -float(np.asarray(state["Z_star"])) / (2.0 * float(np.asarray(state["L"])))
        return X_array * B

    def _continued_scalar(self, X: float, eta: float) -> dict[str, float]:
        X = float(X)
        eta = float(eta)
        if not math.isfinite(X) or X < 0.0:
            raise ValueError("X must be finite and nonnegative")
        if not math.isfinite(eta) or abs(eta) > 1.0:
            raise ValueError("eta must lie in [-1,1]")
        if X <= self.X1:
            natural = self._natural_scalar(X, eta)
            F = natural["F"]
            return {
                "log_F": natural["log_F"],
                "F": F,
                "U": natural["U"],
                "D_X_log_F": float(self._natural_log_slope(np.asarray(X), eta)),
                "D_X_F": natural["D_X_F"],
                "D_X_U": natural["D_X_U"],
            }

        base = self._natural_scalar(self.X1, eta)
        y_upper = min(math.log(X / self.X0), 2.0 * self.log_transition_width)
        y_lower = self.log_transition_width

        def log_integrand(y: np.ndarray) -> np.ndarray:
            xq = self.X0 * np.exp(y)
            s = (y - self.log_transition_width) / self.log_transition_width
            gate = 1.0 - source_smooth_step(s)
            return gate * self._natural_log_slope(xq, eta)

        def u_integrand(y: np.ndarray) -> np.ndarray:
            xq = self.X0 * np.exp(y)
            s = (y - self.log_transition_width) / self.log_transition_width
            gate = 1.0 - source_smooth_step(s)
            return gate * self._natural_DU(xq, eta)

        log_F = base["log_F"] + self._integrate_y(y_lower, y_upper, log_integrand)
        U = base["U"] + self._integrate_y(y_lower, y_upper, u_integrand)
        with np.errstate(under="ignore"):
            F = float(math.exp(log_F)) if log_F > math.log(np.finfo(float).tiny) else 0.0

        if X < self.X2:
            y = math.log(X / self.X0)
            s = (y - self.log_transition_width) / self.log_transition_width
            gate = 1.0 - float(source_smooth_step(np.asarray(s)))
            D_X_log_F = gate * float(self._natural_log_slope(np.asarray(X), eta))
            D_X_U = gate * float(self._natural_DU(np.asarray(X), eta))
        else:
            D_X_log_F = 0.0
            D_X_U = 0.0
        return {
            "log_F": float(log_F),
            "F": F,
            "U": float(U),
            "D_X_log_F": float(D_X_log_F),
            "D_X_F": float(F * D_X_log_F),
            "D_X_U": float(D_X_U),
        }

    def _eta_fd(self, fn: Callable[[float], float], eta: float) -> float:
        h = self.eta_derivative_step
        if eta <= -1.0 + 2.0 * h:
            f0, f1, f2 = fn(eta), fn(eta + h), fn(eta + 2.0 * h)
            return (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h)
        if eta >= 1.0 - 2.0 * h:
            f0, f1, f2 = fn(eta), fn(eta - h), fn(eta - 2.0 * h)
            return (3.0 * f0 - 4.0 * f1 + f2) / (2.0 * h)
        return (fn(eta + h) - fn(eta - h)) / (2.0 * h)

    def _M_scalar(self, X: float, eta: float) -> float:
        X = float(X)
        eta = float(eta)
        if X <= 0.0:
            return 0.0
        state = self.seed.axis_state(np.asarray(eta))
        U_star = float(np.asarray(state["U_star"]))
        B = -float(np.asarray(state["Z_star"])) / (2.0 * float(np.asarray(state["L"])))
        if X <= self.X1:
            return X * U_star + 0.5 * X * X * B
        M1 = self.X1 * U_star + 0.5 * self.X1 * self.X1 * B
        y_upper = min(math.log(X / self.X0), 2.0 * self.log_transition_width)

        def integrand(y: np.ndarray) -> np.ndarray:
            out = np.empty_like(y)
            for i, yy in enumerate(y):
                xx = self.X0 * math.exp(float(yy))
                out[i] = self._continued_scalar(xx, eta)["U"] * xx
            return out

        value = M1 + self._integrate_y(self.log_transition_width, y_upper, integrand)
        if X > self.X2:
            value += (X - self.X2) * self._continued_scalar(self.X2, eta)["U"]
        return float(value)

    def _Pi_scalar(self, X: float, eta: float) -> float:
        X = float(X)
        eta = float(eta)
        if X <= self.X1:
            return self._natural_scalar(X, eta)["Pi"]
        Pi1 = self._natural_scalar(self.X1, eta)["Pi"]
        y_upper = min(math.log(X / self.X0), 2.0 * self.log_transition_width)

        def integrand(y: np.ndarray) -> np.ndarray:
            out = np.empty_like(y)
            for i, yy in enumerate(y):
                xx = self.X0 * math.exp(float(yy))
                F = self._continued_scalar(xx, eta)["F"]
                out[i] = F * F * xx
            return out

        value = Pi1 + self._integrate_y(self.log_transition_width, y_upper, integrand)
        if X > self.X2:
            F2 = self._continued_scalar(self.X2, eta)["F"]
            value += (X - self.X2) * F2 * F2
        if not math.isfinite(value):
            raise OverflowError("continued pressure became non-finite")
        return float(value)

    def _profile_scalar(self, X: float, eta: float) -> dict[str, float]:
        base = self._continued_scalar(X, eta)
        if X <= self.X1:
            natural = self._natural_scalar(X, eta)
            F_eta = natural["F_eta"]
            U_eta = natural["U_eta"]
            v0 = natural["v0"]
            M = self._M_scalar(X, eta)
        else:
            F_eta = self._eta_fd(lambda ee: self._continued_scalar(X, ee)["F"], eta)
            U_eta = self._eta_fd(lambda ee: self._continued_scalar(X, ee)["U"], eta)
            M = self._M_scalar(X, eta)
            M_eta = self._eta_fd(lambda ee: self._M_scalar(X, ee), eta)
            d = 1.0 - eta * eta
            L = 1.0 - 2.0 * self.h * eta * eta
            v0 = (
                2.0 * eta * base["U"]
                - 2.0 * self.D * eta * M / X
                - d * M_eta / X
            ) / L
        E = math.sqrt(2.0 * X) * base["F"] if X > 0.0 else 0.0
        D_X_E = 0.5 * E + math.sqrt(2.0 * X) * base["D_X_F"] if X > 0.0 else 0.0
        return {
            **base,
            "E": float(E),
            "D_X_E": float(D_X_E),
            "F_eta": float(F_eta),
            "U_eta": float(U_eta),
            "M": float(M),
            "v0": float(v0),
            "Pi": self._Pi_scalar(X, eta),
            "Pi_X": float(base["F"] * base["F"]),
        }

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array < 0.0) or np.any(np.abs(eta_array) > 1.0):
            raise ValueError("profile domain requires X>=0 and |eta|<=1")
        keys = (
            "log_F", "F", "E", "U", "v0", "Pi", "Pi_X", "M",
            "D_X_log_F", "D_X_F", "D_X_E", "D_X_U", "F_eta", "U_eta",
        )
        out = {key: np.empty_like(X_array, dtype=float) for key in keys}
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            values = self._profile_scalar(float(xx), float(ee))
            for key in keys:
                out[key].reshape(-1)[index] = values[key]
        if any(np.any(~np.isfinite(values)) for values in out.values()):
            raise OverflowError("rescaled reference profile produced non-finite values")
        return out

    def F(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["F"]

    def U(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["U"]

    def v0(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["v0"]

    def Pi(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["Pi"]

    def F_X(self, X: Any, eta: Any) -> np.ndarray:
        values = self.profile_values(X, eta)
        X_array = np.asarray(np.broadcast_arrays(_finite_array(X, "X"), _finite_array(eta, "eta"))[0])
        out = np.zeros_like(X_array, dtype=float)
        mask = X_array > 0.0
        out[mask] = values["D_X_F"][mask] / X_array[mask]
        return out

    def U_X(self, X: Any, eta: Any) -> np.ndarray:
        values = self.profile_values(X, eta)
        X_array = np.asarray(np.broadcast_arrays(_finite_array(X, "X"), _finite_array(eta, "eta"))[0])
        out = np.zeros_like(X_array, dtype=float)
        mask = X_array > 0.0
        out[mask] = values["D_X_U"][mask] / X_array[mask]
        if np.any(~mask):
            eta_array = np.asarray(np.broadcast_arrays(_finite_array(X, "X"), _finite_array(eta, "eta"))[1])
            for idx in np.argwhere(~mask):
                tup = tuple(idx)
                state = self.seed.axis_state(np.asarray(float(eta_array[tup])))
                out[tup] = -float(np.asarray(state["Z_star"])) / (2.0 * float(np.asarray(state["L"])))
        return out

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            _finite_array(x, "x"), _finite_array(y, "y"),
            _finite_array(z, "z"), _finite_array(t, "t")
        )
        coordinates = self._coordinates.evaluate(x_array, y_array, z_array, t_array)
        profiles = self.profile_values(coordinates["X"], coordinates["eta"])
        q = np.asarray(coordinates["q"], dtype=float)
        swirl_factor = np.power(q, -1.0 - self.h) * profiles["F"]
        radial_factor = profiles["v0"] / (2.0 * q)
        axial_factor = np.power(q, -0.5 - self.h) * profiles["U"]
        result = np.stack(
            (
                radial_factor * x_array - swirl_factor * y_array,
                radial_factor * y_array + swirl_factor * x_array,
                axial_factor,
            ),
            axis=-1,
        )
        if np.any(~np.isfinite(result)):
            raise OverflowError("rescaled reference velocity became non-finite")
        return result

    __call__ = velocity

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = np.asarray(points_xyz, dtype=float)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points_xyz must have final dimension 3")
        if not np.all(np.isfinite(points)):
            raise ValueError("points_xyz must contain only finite values")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def report(self) -> dict[str, Any]:
        eta = float(self.seed.phase_stationary_eta)
        mid = math.sqrt(self.X1 * self.X2)
        sample_X = np.asarray([self.X0, mid, 2.0 * self.X2])
        values = self.profile_values(sample_X, np.full(sample_X.shape, eta))
        return {
            "pressure_scale": self.pressure_scale,
            "rescaling_lambda": self.rescaling_lambda,
            "X0": self.X0,
            "X1": self.X1,
            "X2": self.X2,
            "max_source_transition_X": self.max_source_transition_X,
            "source_transition_fits_rescaled_core_domain": bool(
                self.X2 < self.max_source_transition_X
            ),
            "sample_X": sample_X.tolist(),
            "sample_log_F": values["log_F"].tolist(),
            "sample_U": values["U"].tolist(),
            "sample_v0": values["v0"].tolist(),
            "all_sample_fields_finite": bool(
                all(np.all(np.isfinite(v)) for v in values.values())
            ),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "seed": self.seed.to_payload(),
            "parameters": {
                "log_transition_width": self.log_transition_width,
                "quadrature_points": self.quadrature_points,
                "eta_derivative_step": self.eta_derivative_step,
            },
            "autonomous_numerics": {
                "transition_quadrature": "fixed Gauss-Legendre in source y=log(X/X0)",
                "continued_eta_derivatives": "second-order centered/one-sided finite differences",
                "underflow_policy": "materialized F may underflow to zero while log_F remains finite",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self._unsigned_payload()).encode("utf-8")).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceRescaledReferenceContinuation":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected rescaled-reference-continuation schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("rescaled reference source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("rescaled reference truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(_canonical_json(unsigned).encode("utf-8")).hexdigest()
        if claimed != expected:
            raise ValueError("rescaled reference payload SHA-256 mismatch")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing rescaled reference parameters")
        obj = cls(
            seed=KokunoSourceRescaledCoreSeed.from_payload(payload.get("seed")),
            log_transition_width=float(params["log_transition_width"]),
            quadrature_points=int(params["quadrature_points"]),
            eta_derivative_step=float(params["eta_derivative_step"]),
        )
        if obj.sha256 != claimed:
            raise ValueError("rescaled reference replay changed SHA-256")
        return obj

    def save_json(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return output

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceRescaledReferenceContinuation":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
