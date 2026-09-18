"""Log-stable executable center of Kokuno's Appendix-B rescaled core.

The corrected public reconstruction (2026-09-09, ``axis_exact.md`` in the
pinned 208-page reader) does *not* propagate the very large Appendix-A pressure
datum by a raw Taylor expansion in the physical radial coordinate.  It first
writes

    phi = phi_* Phi,
    U   = U_* + Lambda^-1 u,
    Pi  = Pi_0 + Lambda^-1 p,
    Y   = Lambda X,

and centers the contraction at

    Phi_0 = f_0(Y chi),        u_0 = -Y Z_*/(2 L),
    f_0(z) = sum_n (-z/2)^n / (n! (n+1)!).

This module makes that *center* executable at the source-compatible pressure
datum selected by :mod:`kokuno_pressure_datum_binding`.  It uses a log-amplitude
representation so the large ``Lambda`` scaling cannot overflow merely because
``phi_* = exp(Lambda integral zeta_*)`` is enormous.  The source fixed point,
its contraction threshold, the complex-domain C bound, and global matching are
not solved here.

The default ``Lambda = pressure_scale**2`` is an explicit autonomous numerical
conditioning choice.  It is not the source's hidden/existential threshold and
is deliberately recorded as such.  Likewise the real-axis amplitude is
normalized by the exact real-axis phase maximum; this is a finite executable
normalization, not a proof of the source complex-neighborhood C bound.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import hyp0f1

from .kokuno_pressure_datum_binding import KokunoPressureDatumBinding
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-rescaled-core-seed-v1"

_SOURCE_FORMULAS = {
    "axis_constants": "A=1/2+h; D=1/2-h; d=1-eta^2; L=1-2h*eta^2",
    "axis_U_H_W": (
        "U_*=4 eta+j_0; H_*=D eta+d U_*; "
        "W_*=1-d U_*'-2D eta U_*"
    ),
    "axis_Z": (
        "Z_*=-A(1-2 eta U_*)U_*-H_*U_*'-d Pi_0'+4A eta Pi_0"
    ),
    "chi_zeta": (
        "chi=H_*^2/(H_*^2+sigma_*^2); "
        "zeta_*=-L H_*/(H_*^2+sigma_*^2)"
    ),
    "source_rescaling": (
        "Y=Lambda X; phi=phi_* Phi; U=U_*+Lambda^-1 u; "
        "Pi=Pi_0+Lambda^-1 p"
    ),
    "angular_axis_factor": "phi_*=exp(Lambda*integral_0^eta zeta_*)",
    "pressure_primitive": "p=integral_0^Y g^2 Phi^2 dv; Pi_X=g^2 Phi^2",
    "contraction_center": (
        "Phi_0=f_0(Y chi); u_0=-Y Z_*/(2L); "
        "f_0(z)=sum_{n>=0}(-z/2)^n/(n!(n+1)!)"
    ),
    "center_domain": "0<=Y<=4.1, |eta|<=1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_exact_rescaling_executable": True,
    "source_contraction_center_executable": True,
    "source_scale_shifted_state_binary64_finite": True,
    "log_amplitude_representation_executable": True,
    "selected_real_axis_center_velocity_executable": True,
    "selected_rescaling_lambda_is_autonomous_conditioning_choice": True,
    "selected_real_axis_normalization_is_autonomous": True,
    "source_hidden_rescaling_lambda_recovered": False,
    "source_contraction_threshold_verified": False,
    "source_sigma_star_admissibility_verified": False,
    "source_complex_C_bound_verified": False,
    "source_fixed_point_solved": False,
    "source_pressure_datum_applied_to_selected_global_path": False,
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
class KokunoSourceRescaledCoreSeed:
    """Executable source contraction center in the scaled radial coordinate.

    ``binding`` supplies the selected public-pressure-bound-compatible axis
    datum ``Pi_0=-P^2(1+eta^2)^-2``.  ``rescaling_lambda_multiplier`` chooses
    ``Lambda = multiplier * P^2`` only to keep the shifted variables in a
    useful numerical range.  The corrected source merely requires Lambda to
    exceed a later existence threshold; this class does not certify that
    threshold.
    """

    binding: KokunoPressureDatumBinding = field(default_factory=KokunoPressureDatumBinding)
    rescaling_lambda_multiplier: float = 1.0
    pressure_quadrature_points: int = 32

    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.binding, KokunoPressureDatumBinding):
            raise TypeError("binding must be a KokunoPressureDatumBinding")
        multiplier = float(self.rescaling_lambda_multiplier)
        if not math.isfinite(multiplier) or not 1.0 <= multiplier <= 1.0e6:
            raise ValueError("rescaling_lambda_multiplier must lie in [1,1e6]")
        if isinstance(self.pressure_quadrature_points, bool) or not isinstance(
            self.pressure_quadrature_points, (int, np.integer)
        ):
            raise TypeError("pressure_quadrature_points must be an integer")
        order = int(self.pressure_quadrature_points)
        if not 16 <= order <= 128:
            raise ValueError("pressure_quadrature_points must lie in [16,128]")
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "rescaling_lambda_multiplier", multiplier)
        object.__setattr__(self, "pressure_quadrature_points", order)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)
        object.__setattr__(
            self,
            "_coordinates",
            KokunoNativeSimilarityCoordinates(h=self.h),
        )

    @property
    def h(self) -> float:
        return float(self.binding.outer_schedule.h)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def D(self) -> float:
        return 0.5 - self.h

    @property
    def j0(self) -> float:
        return float(self.binding.j0)

    @property
    def sigma_star(self) -> float:
        # Reuse the explicit selected B.13 sigma, but do not claim the source
        # compact-set admissibility inequality has been certified.
        return float(self.binding.sigma)

    @property
    def pressure_scale(self) -> float:
        return float(self.binding.selected_pressure_scale)

    @property
    def rescaling_lambda(self) -> float:
        value = self.rescaling_lambda_multiplier * self.pressure_scale**2
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("selected rescaling Lambda is outside float range")
        return value

    @staticmethod
    def f0(z: Any) -> np.ndarray:
        values = _finite_array(z, "z")
        return np.asarray(hyp0f1(2.0, -0.5 * values), dtype=float)

    @staticmethod
    def f0_prime(z: Any) -> np.ndarray:
        values = _finite_array(z, "z")
        return np.asarray(-0.25 * hyp0f1(3.0, -0.5 * values), dtype=float)

    def axis_state(self, eta: Any) -> dict[str, np.ndarray]:
        e = _finite_array(eta, "eta")
        if np.any(np.abs(e) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        d = 1.0 - e * e
        L = 1.0 - 2.0 * self.h * e * e
        L_eta = -4.0 * self.h * e
        U_star = 4.0 * e + self.j0
        U_star_eta = np.full_like(e, 4.0)
        H_star = self.D * e + d * U_star
        H_star_eta = self.D - 2.0 * e * U_star + 4.0 * d
        W_star = 1.0 - 4.0 * d - 2.0 * self.D * e * U_star

        p2 = self.pressure_scale**2
        denom = 1.0 + e * e
        Pi_0 = -p2 / np.power(denom, 2)
        Pi_0_eta = 4.0 * p2 * e / np.power(denom, 3)
        Pi_0_etaeta = 4.0 * p2 * (1.0 - 5.0 * e * e) / np.power(denom, 4)

        B = 1.0 - 2.0 * e * U_star
        B_eta = -2.0 * U_star - 8.0 * e
        Z_star = (
            -self.A * B * U_star
            - 4.0 * H_star
            - d * Pi_0_eta
            + 4.0 * self.A * e * Pi_0
        )
        Z_star_eta = (
            -self.A * (B_eta * U_star + 4.0 * B)
            - 4.0 * H_star_eta
            + 2.0 * e * Pi_0_eta
            - d * Pi_0_etaeta
            + 4.0 * self.A * (Pi_0 + e * Pi_0_eta)
        )
        sigma2 = self.sigma_star**2
        denominator = H_star * H_star + sigma2
        chi = H_star * H_star / denominator
        chi_eta = 2.0 * H_star * H_star_eta * sigma2 / (denominator * denominator)
        zeta = -L * H_star / denominator

        return {
            "d": d,
            "L": L,
            "L_eta": L_eta,
            "U_star": U_star,
            "U_star_eta": U_star_eta,
            "H_star": H_star,
            "H_star_eta": H_star_eta,
            "W_star": W_star,
            "Z_star": Z_star,
            "Z_star_eta": Z_star_eta,
            "chi": chi,
            "chi_eta": chi_eta,
            "zeta_star": zeta,
            "Pi_0": Pi_0,
            "Pi_0_eta": Pi_0_eta,
        }

    def _phase_scalar(self, eta: float) -> float:
        eta = float(eta)
        if not math.isfinite(eta) or abs(eta) > 1.0:
            raise ValueError("eta must lie in [-1,1]")
        if eta == 0.0:
            return 0.0

        def integrand(value: float) -> float:
            return float(self.axis_state(np.asarray(value))["zeta_star"])

        value, _ = quad(integrand, 0.0, eta, epsabs=2.0e-13, epsrel=2.0e-13, limit=100)
        return float(value)

    def phase(self, eta: Any) -> np.ndarray:
        values = _finite_array(eta, "eta")
        if np.any(np.abs(values) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        flat = values.reshape(-1)
        out = np.empty_like(flat)
        for index, value in enumerate(flat):
            out[index] = self._phase_scalar(float(value))
        return out.reshape(values.shape)

    @cached_property
    def phase_stationary_eta(self) -> float:
        def h_star(value: float) -> float:
            return float(self.axis_state(np.asarray(value))["H_star"])

        return float(brentq(h_star, -1.0, 1.0, xtol=1.0e-15, rtol=1.0e-15))

    @cached_property
    def real_axis_phase_max(self) -> float:
        candidates = (-1.0, 0.0, self.phase_stationary_eta, 1.0)
        return max(self._phase_scalar(value) for value in candidates)

    @property
    def log_C_real_axis(self) -> float:
        value = self.rescaling_lambda * self.real_axis_phase_max
        if not math.isfinite(value):
            raise OverflowError("real-axis log normalization is outside float range")
        return value

    def log_g(self, eta: Any) -> np.ndarray:
        phase = self.phase(eta)
        # Exact calculus puts the real-axis phase below its stationary maximum.
        # The min only suppresses positive roundoff after two independent
        # quadratures; it does not alter the mathematical profile.
        raw = self.rescaling_lambda * (phase - self.real_axis_phase_max)
        return np.minimum(raw, 0.0)

    def g(self, eta: Any) -> np.ndarray:
        with np.errstate(under="ignore"):
            return np.exp(self.log_g(eta))

    def _scaled_pressure_center(
        self, Y: np.ndarray, chi: np.ndarray, g: np.ndarray
    ) -> np.ndarray:
        flat_Y = np.asarray(Y, dtype=float).reshape(-1)
        flat_chi = np.asarray(chi, dtype=float).reshape(-1)
        flat_g = np.asarray(g, dtype=float).reshape(-1)
        out = np.empty_like(flat_Y)
        for index, (upper, local_chi, local_g) in enumerate(
            zip(flat_Y, flat_chi, flat_g, strict=True)
        ):
            if upper == 0.0 or local_g == 0.0:
                out[index] = 0.0
                continue
            points = 0.5 * upper * (self._nodes + 1.0)
            phi = self.f0(points * local_chi)
            out[index] = (
                local_g
                * local_g
                * 0.5
                * upper
                * float(np.dot(self._weights, phi * phi))
            )
        return out.reshape(np.shape(Y))

    def profile_values_scaled(self, Y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the source contraction center in ``(Y,eta)``.

        Returned ``F,E,U,v0,Pi`` are the selected real-axis-normalized center,
        not the fixed point.  ``log_g`` and ``log_F`` remain informative where
        binary64 necessarily underflows the materialized angular amplitude.
        """

        Y_array, eta_array = np.broadcast_arrays(
            _finite_array(Y, "Y"), _finite_array(eta, "eta")
        )
        if np.any((Y_array < 0.0) | (Y_array > 4.1)):
            raise ValueError("source contraction-center domain requires 0<=Y<=4.1")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")

        state = self.axis_state(eta_array)
        local_log_g = self.log_g(eta_array)
        with np.errstate(under="ignore"):
            local_g = np.exp(local_log_g)
        z = Y_array * state["chi"]
        Phi_0 = self.f0(z)
        Phi_z = self.f0_prime(z)
        if np.any(Phi_0 <= 0.0) or np.any(~np.isfinite(Phi_0)):
            raise RuntimeError("f_0 left its positive source center range")

        Lambda = self.rescaling_lambda
        X = Y_array / Lambda
        B = -state["Z_star"] / (2.0 * state["L"])
        B_eta = (
            -state["Z_star_eta"] / (2.0 * state["L"])
            + state["Z_star"] * state["L_eta"] / (2.0 * state["L"] ** 2)
        )
        u_0 = Y_array * B
        U = state["U_star"] + u_0 / Lambda
        U_eta = 4.0 + X * B_eta
        average_U = state["U_star"] + 0.5 * X * B
        average_U_eta = 4.0 + 0.5 * X * B_eta
        v0 = (
            2.0 * eta_array * U
            - 2.0 * self.D * eta_array * average_U
            - state["d"] * average_U_eta
        ) / state["L"]
        v0_X = (
            (2.0 - self.D) * eta_array * B
            - 0.5 * state["d"] * B_eta
        ) / state["L"]

        F = local_g * Phi_0
        log_F = local_log_g + np.log(Phi_0)
        D_X_F = local_g * Phi_z * Y_array * state["chi"]
        F_eta = local_g * (
            Lambda * state["zeta_star"] * Phi_0
            + Phi_z * Y_array * state["chi_eta"]
        )
        E = np.sqrt(2.0 * X) * F
        D_X_E = 0.5 * E + np.sqrt(2.0 * X) * D_X_F
        p_0 = self._scaled_pressure_center(Y_array, state["chi"], local_g)
        Pi = state["Pi_0"] + p_0 / Lambda

        return {
            "Y": Y_array,
            "X": X,
            "Phi_0": Phi_0,
            "u_0": u_0,
            "log_g": local_log_g,
            "g": local_g,
            "log_F": log_F,
            "F": F,
            "D_X_F": D_X_F,
            "F_eta": F_eta,
            "E": E,
            "D_X_E": D_X_E,
            "U": U,
            "D_X_U": X * B,
            "U_eta": U_eta,
            "average_U": average_U,
            "average_U_eta": average_U_eta,
            "v0": v0,
            "D_X_v0": X * v0_X,
            "Pi": Pi,
            "Pi_X": F * F,
            "scaled_pressure_correction": p_0,
            "physical_pressure_correction": p_0 / Lambda,
            **state,
        }

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array < 0.0):
            raise ValueError("X must be nonnegative")
        Y = self.rescaling_lambda * X_array
        return self.profile_values_scaled(Y, eta_array)

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate the selected center velocity in Kokuno native coordinates.

        Calls outside ``Lambda X<=4.1`` fail closed.  This velocity is useful as
        an executable large-pressure core seed only; it is not the source fixed
        point or a global candidate.
        """

        coordinates = self._coordinates.evaluate(x, y, z, t)
        profiles = self.profile_values(coordinates["X"], coordinates["eta"])
        x_array, y_array = np.broadcast_arrays(_finite_array(x, "x"), _finite_array(y, "y"))
        x_array = np.broadcast_to(x_array, np.shape(profiles["F"]))
        y_array = np.broadcast_to(y_array, np.shape(profiles["F"]))
        q = np.asarray(coordinates["q"], dtype=float)
        swirl_factor = np.power(q, -1.0 - self.h) * profiles["F"]
        radial_factor = profiles["v0"] / (2.0 * q)
        axial_factor = np.power(q, -0.5 - self.h) * profiles["U"]
        out = np.stack(
            (
                radial_factor * x_array - swirl_factor * y_array,
                radial_factor * y_array + swirl_factor * x_array,
                axial_factor,
            ),
            axis=-1,
        )
        if np.any(~np.isfinite(out)):
            raise OverflowError("selected rescaled-core center velocity became non-finite")
        return out

    __call__ = velocity

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = np.asarray(points_xyz, dtype=float)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points_xyz must have final dimension 3")
        if not np.all(np.isfinite(points)):
            raise ValueError("points_xyz must contain only finite values")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def report(self) -> dict[str, Any]:
        sample_eta = np.asarray([-1.0, self.phase_stationary_eta, 0.0, 1.0])
        center = self.profile_values_scaled(np.full(sample_eta.shape, 2.0), sample_eta)
        return {
            "pressure_scale": self.pressure_scale,
            "rescaling_lambda": self.rescaling_lambda,
            "lambda_to_pressure_square_ratio": (
                self.rescaling_lambda / (self.pressure_scale**2)
            ),
            "phase_stationary_eta": self.phase_stationary_eta,
            "real_axis_phase_max": self.real_axis_phase_max,
            "log_C_real_axis": self.log_C_real_axis,
            "sample_eta": sample_eta.tolist(),
            "sample_log_g": center["log_g"].tolist(),
            "sample_U": center["U"].tolist(),
            "sample_v0": center["v0"].tolist(),
            "all_sample_shifted_fields_finite": bool(
                np.all(np.isfinite(center["u_0"]))
                and np.all(np.isfinite(center["U"]))
                and np.all(np.isfinite(center["Pi"]))
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
            "binding": self.binding.to_payload(),
            "parameters": {
                "rescaling_lambda_multiplier": self.rescaling_lambda_multiplier,
                "pressure_quadrature_points": self.pressure_quadrature_points,
            },
            "autonomous_numerics": {
                "Lambda_choice": "rescaling_lambda_multiplier*selected_pressure_scale^2",
                "real_axis_C": "exp(Lambda*max_{eta in [-1,1]} integral_0^eta zeta_*) stored only as log_C",
                "phase_quadrature": "scipy adaptive Gauss-Kronrod on the real eta interval",
                "pressure_primitive": "fixed Gauss-Legendre in scaled Y",
                "positive_roundoff_log_g": "clipped to zero because exact real-axis phase is <= its maximum",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceRescaledCoreSeed":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno source-rescaled core-seed schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source-rescaled core formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("source-rescaled core truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(_canonical_json(unsigned).encode("utf-8")).hexdigest()
        if claimed != expected:
            raise ValueError("source-rescaled core payload SHA mismatch")
        obj = cls(
            binding=KokunoPressureDatumBinding.from_payload(payload.get("binding")),
            **dict(payload.get("parameters", {})),
        )
        if obj.to_payload() != payload:
            raise ValueError("source-rescaled core payload does not replay exactly")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceRescaledCoreSeed":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
