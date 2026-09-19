"""Executable selected-center pressure and radial-derivative primitives for PA.10.

Pinned public provenance is ``KokunoYumeto/yang-mills-interacting-workbench`` at
commit ``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

For the rescaled core the corrected reconstruction uses

    Phi_0(Y,eta) = f_0(Y chi(eta)),
    p(Y,eta) = integral_0^Y g(eta)^2 Phi_0(v,eta)^2 dv.

This module makes the pressure quantities that appear explicitly in PA.10's
R2 remainder executable on the repository-selected contraction center:

    p,
    p_eta,
    Y p_Y,

and, at the same time, exposes the selected-center factors ``Phi``,
``Phi_eta`` and ``Y Phi_Y`` needed by R1.  Differentiating the displayed
pressure primitive gives

    p_Y = g^2 Phi^2,
    p_eta = g^2 [2 Lambda zeta_* I + I_eta],
    I_eta = integral_0^Y 2 Phi(v) f_0'(v chi) v chi_eta dv.

The formula is evaluated directly rather than finite-differencing the huge
``Lambda`` scale.  Vectorized evaluation, deterministic engineering envelopes,
and SHA-bound save/load are provided.

Important truth boundary: this is the *selected contraction center*, not the
source fixed point and not a radius-one-ball coefficient-space certificate.
The sampled envelopes below therefore must not be inserted into the source
M/K proof as if they were source ball bounds.  They are executable candidate
inputs and diagnostics for the next rigorous bridge only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_rescaled_core_seed import KokunoSourceRescaledCoreSeed


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-selected-pressure-primitive-v1"

_SOURCE_FORMULAS = {
    "center_profile": "Phi_0(Y,eta)=f_0(Y chi(eta))",
    "pressure_primitive": "p(Y,eta)=integral_0^Y g(eta)^2 Phi_0(v,eta)^2 dv",
    "pressure_radial_derivative": "p_Y=g^2 Phi_0^2",
    "pressure_eta_derivative": (
        "p_eta=g^2[2 Lambda zeta_* I+I_eta], "
        "I_eta=integral_0^Y 2 Phi_0(v) f_0'(v chi) v chi_eta dv"
    ),
    "profile_radial_derivative": "Y Phi_Y=Y chi f_0'(Y chi)",
    "profile_eta_derivative": "Phi_eta=Y chi_eta f_0'(Y chi)",
    "center_domain": "0<=Y<=4.1, |eta|<=1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "selected_center_pressure_primitive_executable": True,
    "selected_center_pressure_eta_derivative_executable": True,
    "selected_center_Y_p_Y_executable": True,
    "selected_center_Y_Phi_Y_executable": True,
    "vectorized_evaluation_executable": True,
    "selected_center_engineering_envelope_executable": True,
    "selected_center_engineering_envelope_is_source_ball_bound": False,
    "source_fixed_point_pressure_bound_machine_bound": False,
    "source_pressure_radius_one_ball_norm_machine_bound": False,
    "source_pressure_radius_one_ball_lipschitz_machine_bound": False,
    "source_mixed_Y_Phi_Y_radius_one_ball_bound_machine_bound": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_B0_dependencies_machine_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "heldout_ns_residual_assessed": False,
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
class KokunoPA10SelectedPressurePrimitive:
    """Selected-center implementation of the PA.10 pressure primitive seam."""

    core: KokunoSourceRescaledCoreSeed = field(default_factory=KokunoSourceRescaledCoreSeed)
    quadrature_points: int = 48

    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.core, KokunoSourceRescaledCoreSeed):
            raise TypeError("core must be a KokunoSourceRescaledCoreSeed")
        if isinstance(self.quadrature_points, bool) or not isinstance(
            self.quadrature_points, (int, np.integer)
        ):
            raise TypeError("quadrature_points must be an integer")
        order = int(self.quadrature_points)
        if not 16 <= order <= 128:
            raise ValueError("quadrature_points must lie in [16,128]")
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "quadrature_points", order)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)

    def _integrals(
        self,
        Y: np.ndarray,
        chi: np.ndarray,
        chi_eta: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        flat_Y = np.asarray(Y, dtype=float).reshape(-1)
        flat_chi = np.asarray(chi, dtype=float).reshape(-1)
        flat_chi_eta = np.asarray(chi_eta, dtype=float).reshape(-1)
        I = np.empty_like(flat_Y)
        I_eta = np.empty_like(flat_Y)
        for index, (upper, local_chi, local_chi_eta) in enumerate(
            zip(flat_Y, flat_chi, flat_chi_eta, strict=True)
        ):
            if upper == 0.0:
                I[index] = 0.0
                I_eta[index] = 0.0
                continue
            v = 0.5 * upper * (self._nodes + 1.0)
            phi = self.core.f0(v * local_chi)
            phi_prime = self.core.f0_prime(v * local_chi)
            scale = 0.5 * upper
            I[index] = scale * float(np.dot(self._weights, phi * phi))
            eta_integrand = 2.0 * phi * phi_prime * v * local_chi_eta
            I_eta[index] = scale * float(np.dot(self._weights, eta_integrand))
        return I.reshape(np.shape(Y)), I_eta.reshape(np.shape(Y))

    def evaluate(self, Y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate selected-center profile/pressure primitives vectorially."""

        Y_array, eta_array = np.broadcast_arrays(
            _finite_array(Y, "Y"), _finite_array(eta, "eta")
        )
        if np.any((Y_array < 0.0) | (Y_array > 4.1)):
            raise ValueError("source contraction-center domain requires 0<=Y<=4.1")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")

        state = self.core.axis_state(eta_array)
        log_g = self.core.log_g(eta_array)
        with np.errstate(under="ignore"):
            g = np.exp(log_g)
        arg = Y_array * state["chi"]
        Phi = self.core.f0(arg)
        Phi_prime = self.core.f0_prime(arg)
        Y_Phi_Y = Y_array * state["chi"] * Phi_prime
        Phi_eta = Y_array * state["chi_eta"] * Phi_prime

        I, I_eta = self._integrals(Y_array, state["chi"], state["chi_eta"])
        g2 = g * g
        p = g2 * I
        p_Y = g2 * Phi * Phi
        Y_p_Y = Y_array * p_Y
        p_eta = g2 * (
            2.0 * self.core.rescaling_lambda * state["zeta_star"] * I + I_eta
        )

        result = {
            "Y": Y_array,
            "eta": eta_array,
            "Phi": Phi,
            "Phi_eta": Phi_eta,
            "Y_Phi_Y": Y_Phi_Y,
            "log_g": log_g,
            "g": g,
            "pressure_integral_I": I,
            "pressure_integral_I_eta": I_eta,
            "p": p,
            "p_Y": p_Y,
            "Y_p_Y": Y_p_Y,
            "p_eta": p_eta,
        }
        if not all(np.all(np.isfinite(value)) for value in result.values()):
            raise OverflowError("selected pressure primitive produced non-finite values")
        return result

    __call__ = evaluate

    @staticmethod
    def _chebyshev_lobatto(count: int) -> np.ndarray:
        if isinstance(count, bool) or not isinstance(count, int) or count < 3:
            raise ValueError("grid count must be an integer >=3")
        k = np.arange(count, dtype=float)
        return np.cos(np.pi * k / float(count - 1))[::-1]

    def engineering_envelope(
        self,
        *,
        y_count: int = 25,
        eta_count: int = 49,
        padding_fraction: float = 1.0e-9,
    ) -> dict[str, Any]:
        """Return a nested-grid selected-center engineering envelope.

        This is deliberately a numerical envelope, not a certified continuum
        supremum and not a source radius-one-ball coefficient-space bound.
        """

        if isinstance(y_count, bool) or not isinstance(y_count, int) or y_count < 5:
            raise ValueError("y_count must be an integer >=5")
        if isinstance(eta_count, bool) or not isinstance(eta_count, int) or eta_count < 5:
            raise ValueError("eta_count must be an integer >=5")
        pad = float(padding_fraction)
        if not math.isfinite(pad) or pad < 0.0 or pad > 0.1:
            raise ValueError("padding_fraction must lie in [0,0.1]")

        def sample(ny: int, ne: int) -> dict[str, float]:
            # Chebyshev-Lobatto includes endpoints and clusters near both ends.
            y_nodes = 2.05 * (self._chebyshev_lobatto(ny) + 1.0)
            eta_nodes = self._chebyshev_lobatto(ne)
            Y_grid, eta_grid = np.meshgrid(y_nodes, eta_nodes, indexing="ij")
            values = self.evaluate(Y_grid, eta_grid)
            names = ("Phi", "Phi_eta", "Y_Phi_Y", "p", "p_eta", "Y_p_Y")
            return {
                name: float(np.max(np.abs(values[name])))
                for name in names
            }

        coarse = sample(y_count, eta_count)
        fine = sample(2 * y_count - 1, 2 * eta_count - 1)
        envelope = {
            name: float(math.nextafter(max(coarse[name], fine[name]) * (1.0 + pad), math.inf))
            for name in fine
        }
        drift = {
            name: float(abs(fine[name] - coarse[name]) / max(fine[name], 1.0e-300))
            for name in fine
        }
        return {
            "coarse_grid": {"y_count": y_count, "eta_count": eta_count},
            "fine_grid": {"y_count": 2 * y_count - 1, "eta_count": 2 * eta_count - 1},
            "padding_fraction": pad,
            "coarse_sampled_abs_max": coarse,
            "fine_sampled_abs_max": fine,
            "relative_nested_grid_drift": drift,
            "engineering_abs_envelope": envelope,
            "continuum_supremum_certified": False,
            "source_radius_one_ball_bound": False,
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        envelope = self.engineering_envelope()
        eta0 = self.core.phase_stationary_eta
        sample = self.evaluate(np.asarray([0.0, 2.0, 4.1]), np.asarray([eta0, eta0, eta0]))
        return {
            "selected_lambda": self.core.rescaling_lambda,
            "phase_stationary_eta": eta0,
            "stationary_sample_p": sample["p"].tolist(),
            "stationary_sample_p_eta": sample["p_eta"].tolist(),
            "envelope": envelope,
            "truth_boundary": self.truth_boundary,
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
            "core": self.core.to_payload(),
            "parameters": {"quadrature_points": self.quadrature_points},
            "truth_boundary": self.truth_boundary,
        }

    @property
    def sha256(self) -> str:
        raw = _canonical_json(self._unsigned_payload()).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10SelectedPressurePrimitive":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected selected-pressure-primitive schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("selected-pressure source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("selected-pressure truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(_canonical_json(unsigned).encode("utf-8")).hexdigest()
        if claimed != expected:
            raise ValueError("selected-pressure payload SHA mismatch")
        obj = cls(
            core=KokunoSourceRescaledCoreSeed.from_payload(payload.get("core")),
            **dict(payload.get("parameters", {})),
        )
        if obj.sha256 != claimed:
            raise ValueError("selected-pressure reconstructed SHA mismatch")
        return obj

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n")

    @classmethod
    def load(cls, path: str | Path) -> "KokunoPA10SelectedPressurePrimitive":
        return cls.from_payload(json.loads(Path(path).read_text()))
