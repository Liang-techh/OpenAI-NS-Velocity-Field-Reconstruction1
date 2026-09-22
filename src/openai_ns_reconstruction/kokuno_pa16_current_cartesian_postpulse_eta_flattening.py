"""Carry the current A1 Cartesian leading field through Kokuno's eta-flattening stage.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
corrected 2026-09-09 reconstruction.

Immediately after the public pulse endpoint the reconstruction resets the local
coordinate y=log(X/X_end), sets U=0 permanently, and prescribes

    log E = log e_end -(1/2+lambda)y
            - theta_f(y) J_0 -(1-theta_f(y)) log 2,
    theta_f(y) = 1 - sigma(y/T_f),
    J_0(eta) = log(1+eta^2),

where sigma is the public flat C-infinity step.  The source only requires a
fixed sufficiently large T_f; it does not supply a unique numerical value.
This bounded repository realization therefore freezes T_f=100.0 as an
autonomous numerical choice.  It is not source-exact hidden data.

For an exact seam with the current #1133 candidate we use the equivalent form

    log E = log E_end -(1/2+lambda)y
            + sigma(y/T_f) (J_0-log 2),

with E_end read from the actual parent at xi=13.  Since that parent endpoint is
proportional to (1+eta^2)^(-1), this becomes eta-independent at y=T_f up to
floating-point roundoff, exactly as the public schedule requires.

The current M/M_eta primitive is not snapped to the idealized source value.
Because U=0 on this stage, physical M and M_eta are constant, hence

    M/X = exp(-y) (M/X)_end,
    M_eta/X = exp(-y) (M_eta/X)_end.

The radial profile is rebuilt from the full inherited primitive history.  No
pressure, forcing, terminal hold, relative-swirl correction, exterior heat,
complete held-out NS residual, or PDE validation is added here.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_pulse_end_compensated import (
    KokunoPA16CurrentCartesianPulseEndCompensated,
)

SCHEMA = "kokuno-pa16-current-cartesian-postpulse-eta-flattening-v1"
PARENT_EXACT_HEAD = "3b4af71c4ab547bc11f9f3e320afa4b62c25b930"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
T_F_AUTONOMOUS = 100.0
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "local_coordinate": "y=log(X/X_end), reset y=0 at the pulse endpoint",
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "J0": "J_0(eta)=log(1+eta^2)",
    "theta_f": "theta_f(y)=1-sigma(y/T_f)",
    "postpulse_E": (
        "log E=log e_end-(1/2+lambda)y-theta_f J_0-(1-theta_f)log 2"
    ),
    "public_l_window": "-lambda-0.1 <= l <= -lambda, l=X*d_X log(sqrt(2X)E)",
    "smooth_step": (
        "sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)) on 0<s<1; "
        "sigma=0 for s<=0 and sigma=1 for s>=1"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1133 Cartesian pulse endpoint identity",
    "T_f": "repository-autonomous frozen T_f=100.0; source supplies no unique numerical value",
    "seam": "use actual parent E/M/M_eta at xi=13 rather than idealized zeroing",
    "coordinate": "overflow-safe log-X/log-F continuation",
    "current_M_policy": "with public U=0, carry physical M/M_eta constantly so M/X scales exp(-y)",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "source_postpulse_eta_flattening_formula_materialized": True,
    "repository_autonomous_T_f_materialized": True,
    "source_exact_T_f_recovered": False,
    "current_cartesian_postpulse_eta_flattening_composed": True,
    "source_terminal_hold_after_eta_flattening_materialized": False,
    "source_relative_swirl_bumps_materialized": False,
    "source_exterior_heat_replacement_materialized": False,
    "outer_global_leading_velocity_materialized": False,
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


def _smooth_step_and_derivative(s: Any) -> tuple[np.ndarray, np.ndarray]:
    """Public flat step sigma and d sigma/ds, stably evaluated in float64."""

    values = np.asarray(s, dtype=float)
    sigma = np.zeros_like(values)
    deriv = np.zeros_like(values)
    sigma[values >= 1.0] = 1.0
    interior = (values > 0.0) & (values < 1.0)
    if not np.any(interior):
        return sigma, deriv

    u = values[interior]
    with np.errstate(over="ignore", divide="ignore", invalid="ignore", under="ignore"):
        logit = -1.0 / (u * u) + 1.0 / ((1.0 - u) * (1.0 - u))
        sig = np.empty_like(u)
        positive = logit >= 0.0
        sig[positive] = 1.0 / (1.0 + np.exp(-logit[positive]))
        exp_logit = np.exp(logit[~positive])
        sig[~positive] = exp_logit / (1.0 + exp_logit)
        local_deriv = np.zeros_like(u)
        live = (sig > 0.0) & (sig < 1.0)
        if np.any(live):
            ul = u[live]
            logit_prime = 2.0 / (ul**3) + 2.0 / ((1.0 - ul) ** 3)
            local_deriv[live] = sig[live] * (1.0 - sig[live]) * logit_prime

    sigma[interior] = sig
    deriv[interior] = local_deriv
    return sigma, deriv


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianPostPulseEtaFlattening:
    """Current Cartesian leading candidate through the immediate eta-flattening stage."""

    parent: KokunoPA16CurrentCartesianPulseEndCompensated = field(
        default_factory=KokunoPA16CurrentCartesianPulseEndCompensated,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianPulseEndCompensated):
            raise TypeError("parent must be KokunoPA16CurrentCartesianPulseEndCompensated")
        truth = self.parent.truth_boundary
        if not truth["current_pulse_endpoint_xi13_materialized"]:
            raise ValueError("post-pulse flattening requires the current xi=13 pulse endpoint")
        if truth["source_terminal_tail_schedule_bound_into_current_velocity"]:
            raise ValueError("parent unexpectedly already contains a post-pulse terminal schedule")
        if not math.isfinite(T_F_AUTONOMOUS) or T_F_AUTONOMOUS <= 0.0:
            raise ValueError("autonomous T_f must be finite and positive")
        if self.log_radius_q1_flatten_end >= math.log(np.finfo(float).max):
            raise ValueError("eta-flattening endpoint has non-representable physical radius")

    @property
    def leading(self):
        return self.parent

    @property
    def geometry(self):
        return self.parent.geometry

    @property
    def A(self) -> float:
        return float(self.parent.A)

    @property
    def D(self) -> float:
        return float(self.parent.D)

    @property
    def lambda_value(self) -> float:
        return float(self.parent.lambda_value)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.parent.eta_interval)

    @property
    def T_f(self) -> float:
        return T_F_AUTONOMOUS

    @property
    def log_X_pulse_end(self) -> float:
        return float(self.parent.log_X_pulse_end)

    @property
    def log_X_flatten_end(self) -> float:
        return self.log_X_pulse_end + self.T_f

    @property
    def log_radius_q1_flatten_end(self) -> float:
        return 0.5 * (_LOG2 + self.log_X_flatten_end)

    def similarity_coordinates_logX(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.parent.similarity_coordinates_logX(x, y, z, t)

    def _broadcast_log_similarity(
        self, log_X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_flatten_end))
        if np.any(log_arr > self.log_X_flatten_end + tol):
            raise ValueError("log_X lies beyond the bounded eta-flattening endpoint")
        log_arr = np.minimum(log_arr, self.log_X_flatten_end)
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return log_arr, eta_arr

    def _endpoint_state(self, eta: np.ndarray) -> dict[str, np.ndarray]:
        log_end = np.full(np.asarray(eta).shape, self.log_X_pulse_end, dtype=float)
        p = self.parent.similarity_profile_values_logX(log_end, eta)
        E = np.asarray(p["E_current_leading_with_pulse_end"], dtype=float)
        M_ratio = np.asarray(p["M_over_X_current_leading_with_pulse_end"], dtype=float)
        M_eta_ratio = np.asarray(
            p["M_eta_over_X_current_leading_with_pulse_end"], dtype=float
        )
        if np.any(E <= 0.0) or np.any(~np.isfinite(E)):
            raise RuntimeError("parent pulse endpoint must have finite positive E")
        return {"E": E, "M_over_X": M_ratio, "M_eta_over_X": M_eta_ratio}

    def _postpulse_profile_logX(
        self, log_X: np.ndarray, eta: np.ndarray
    ) -> dict[str, np.ndarray]:
        y = np.asarray(log_X, dtype=float) - self.log_X_pulse_end
        if np.any((y < -2.0e-12) | (y > self.T_f + 2.0e-12)):
            raise ValueError("post-pulse eta-flattening profile requires 0<=y<=T_f")
        y = np.clip(y, 0.0, self.T_f)
        state = self._endpoint_state(eta)

        s = y / self.T_f
        sigma, sigma_prime = _smooth_step_and_derivative(s)
        sigma_y = sigma_prime / self.T_f
        J0 = np.log1p(eta * eta)
        J0_eta = 2.0 * eta / (1.0 + eta * eta)
        log_E = (
            np.log(state["E"])
            - (0.5 + self.lambda_value) * y
            + sigma * (J0 - _LOG2)
        )
        E = np.exp(log_E)
        E_eta = -(1.0 - sigma) * J0_eta * E
        log_F = log_E - 0.5 * (_LOG2 + log_X)
        F = np.exp(log_F)
        U = np.zeros_like(E)
        U_eta = np.zeros_like(E)

        carry = np.exp(-y)
        M_ratio = carry * state["M_over_X"]
        M_eta_ratio = carry * state["M_eta_over_X"]

        axis = self.geometry.physical_profiles.axis_profiles.values(
            np.zeros_like(eta), eta
        )
        L = np.asarray(axis["L"], dtype=float)
        d = np.asarray(axis["d"], dtype=float)
        v0 = (-2.0 * self.D * eta * M_ratio - d * M_eta_ratio) / L

        dlogE = (
            -(0.5 + self.lambda_value)
            + sigma_y * (J0 - _LOG2)
        )
        ell = 0.5 + dlogE
        arrays = (E, E_eta, F, U, M_ratio, M_eta_ratio, v0, dlogE, ell)
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("post-pulse eta-flattening profile produced non-finite values")
        if np.any(E <= 0.0) or np.any(F < 0.0):
            raise RuntimeError("post-pulse eta-flattening lost positive E/nonnegative F")

        return {
            "log_X": np.asarray(log_X, dtype=float),
            "eta": np.asarray(eta, dtype=float),
            "y": y,
            "sigma": sigma,
            "sigma_prime": sigma_prime,
            "J0": J0,
            "E": E,
            "E_eta": E_eta,
            "F": F,
            "log_F": log_F,
            "U": U,
            "U_eta": U_eta,
            "M_over_X": M_ratio,
            "M_eta_over_X": M_eta_ratio,
            "v0": v0,
            "E_DlogX": E * dlogE,
            "F_DlogX": F * (dlogE - 0.5),
            "log_F_DlogX": dlogE - 0.5,
            "U_DlogX": np.zeros_like(U),
            "ell": ell,
        }

    def similarity_profile_values_logX(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)

        F = np.empty_like(lf)
        log_F = np.empty_like(lf)
        U = np.empty_like(lf)
        E = np.empty_like(lf)
        M_ratio = np.empty_like(lf)
        M_eta_ratio = np.empty_like(lf)
        v0 = np.empty_like(lf)
        stage_y = np.zeros_like(lf)
        sigma = np.zeros_like(lf)
        region = np.empty(lf.shape, dtype=object)

        inherited = lf <= self.log_X_pulse_end
        if np.any(inherited):
            p = self.parent.similarity_profile_values_logX(lf[inherited], ef[inherited])
            F[inherited] = np.asarray(p["F_current_leading_with_pulse_end"], dtype=float)
            log_F[inherited] = np.asarray(
                p["log_F_current_leading_with_pulse_end"], dtype=float
            )
            U[inherited] = np.asarray(p["U_current_leading_with_pulse_end"], dtype=float)
            E[inherited] = np.asarray(p["E_current_leading_with_pulse_end"], dtype=float)
            M_ratio[inherited] = np.asarray(
                p["M_over_X_current_leading_with_pulse_end"], dtype=float
            )
            M_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_leading_with_pulse_end"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_leading_with_pulse_end"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        active = ~inherited
        if np.any(active):
            p = self._postpulse_profile_logX(lf[active], ef[active])
            F[active] = p["F"]
            log_F[active] = p["log_F"]
            U[active] = p["U"]
            E[active] = p["E"]
            M_ratio[active] = p["M_over_X"]
            M_eta_ratio[active] = p["M_eta_over_X"]
            v0[active] = p["v0"]
            stage_y[active] = p["y"]
            sigma[active] = p["sigma"]
            region[active] = "current_cartesian_postpulse_eta_flattening_logX"

        arrays = (F, log_F, U, E, M_ratio, M_eta_ratio, v0)
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("post-pulse Cartesian candidate produced non-finite values")
        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "postpulse_y": stage_y.reshape(shape),
            "sigma_postpulse_eta_flattening": sigma.reshape(shape),
            "F_current_leading_postpulse_eta_flattening": F.reshape(shape),
            "log_F_current_leading_postpulse_eta_flattening": log_F.reshape(shape),
            "U_current_leading_postpulse_eta_flattening": U.reshape(shape),
            "E_current_leading_postpulse_eta_flattening": E.reshape(shape),
            "M_over_X_current_leading_postpulse_eta_flattening": M_ratio.reshape(shape),
            "M_eta_over_X_current_leading_postpulse_eta_flattening": M_eta_ratio.reshape(shape),
            "v0_current_leading_postpulse_eta_flattening": v0.reshape(shape),
        }

    def similarity_log_radial_derivatives(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)
        F_D = np.empty_like(lf)
        log_F_D = np.empty_like(lf)
        U_D = np.empty_like(lf)
        E_D = np.empty_like(lf)
        ell = np.empty_like(lf)

        inherited = lf <= self.log_X_pulse_end
        if np.any(inherited):
            d = self.parent.similarity_log_radial_derivatives(lf[inherited], ef[inherited])
            F_D[inherited] = np.asarray(
                d["F_DlogX_current_leading_with_pulse_end"], dtype=float
            )
            log_F_D[inherited] = np.asarray(
                d["log_F_DlogX_current_leading_with_pulse_end"], dtype=float
            )
            U_D[inherited] = np.asarray(
                d["U_DlogX_current_leading_with_pulse_end"], dtype=float
            )
            E_D[inherited] = np.asarray(
                d["E_DlogX_current_leading_with_pulse_end"], dtype=float
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                ell[inherited] = 0.5 + E_D[inherited] / np.asarray(
                    self.parent.similarity_profile_values_logX(
                        lf[inherited], ef[inherited]
                    )["E_current_leading_with_pulse_end"],
                    dtype=float,
                )

        active = ~inherited
        if np.any(active):
            p = self._postpulse_profile_logX(lf[active], ef[active])
            F_D[active] = p["F_DlogX"]
            log_F_D[active] = p["log_F_DlogX"]
            U_D[active] = p["U_DlogX"]
            E_D[active] = p["E_DlogX"]
            ell[active] = p["ell"]

        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_DlogX_current_leading_postpulse_eta_flattening": F_D.reshape(shape),
            "log_F_DlogX_current_leading_postpulse_eta_flattening": log_F_D.reshape(shape),
            "U_DlogX_current_leading_postpulse_eta_flattening": U_D.reshape(shape),
            "E_DlogX_current_leading_postpulse_eta_flattening": E_D.reshape(shape),
            "ell_current_leading_postpulse_eta_flattening": ell.reshape(shape),
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        coords = self.similarity_coordinates_logX(x, y, z, t)
        xb = np.asarray(coords["x"], dtype=float)
        yb = np.asarray(coords["y"], dtype=float)
        zb = np.asarray(coords["z"], dtype=float)
        tb = np.asarray(coords["t"], dtype=float)
        q = np.asarray(coords["q"], dtype=float)
        eta = np.asarray(coords["eta"], dtype=float)
        radius = np.asarray(coords["radius"], dtype=float)
        log_X = np.asarray(coords["log_X"], dtype=float)

        out = np.empty(radius.shape + (3,), dtype=float)
        axis = radius == 0.0
        if np.any(axis):
            out[axis] = self.parent.velocity(xb[axis], yb[axis], zb[axis], tb[axis])

        nonaxis = ~axis
        if np.any(nonaxis):
            p = self.similarity_profile_values_logX(log_X[nonaxis], eta[nonaxis])
            log_F = np.asarray(
                p["log_F_current_leading_postpulse_eta_flattening"], dtype=float
            )
            U = np.asarray(p["U_current_leading_postpulse_eta_flattening"], dtype=float)
            v0 = np.asarray(p["v0_current_leading_postpulse_eta_flattening"], dtype=float)
            unit_x = xb[nonaxis] / radius[nonaxis]
            unit_y = yb[nonaxis] / radius[nonaxis]
            radial_radius = v0 * radius[nonaxis] / (2.0 * q[nonaxis])
            log_swirl_radius = (
                (-self.A - 0.5) * np.log(q[nonaxis])
                + log_F
                + np.log(radius[nonaxis])
            )
            swirl_radius = np.exp(log_swirl_radius)
            u1 = radial_radius * unit_x - swirl_radius * unit_y
            u2 = radial_radius * unit_y + swirl_radius * unit_x
            axial = q[nonaxis] ** (-self.A) * U
            out[nonaxis] = np.stack((u1, u2, axial), axis=-1)

        if np.any(~np.isfinite(out)):
            raise RuntimeError("post-pulse eta-flattening Cartesian velocity became non-finite")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        truth = dict(self.parent.truth_boundary)
        truth.update(_TRUTH_UPDATES)
        return truth

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_current_cartesian_pulse_end": self.parent.configuration(),
            "bound_scope": {
                "parent_exact_head": PARENT_EXACT_HEAD,
                "source": {
                    "repository": SOURCE_REPOSITORY,
                    "commit": SOURCE_COMMIT,
                    "blob": SOURCE_BLOB,
                    "path": SOURCE_PATH,
                    "release": SOURCE_RELEASE,
                    "release_date": SOURCE_RELEASE_DATE,
                },
                "postpulse_stage": "eta_flattening_only",
                "T_f": self.T_f,
                "T_f_role": "repository_autonomous_not_source_exact",
                "current_M_policy": "carry_actual_parent_M_and_M_eta_with_U_zero",
                "terminal_hold_materialized": False,
                "relative_swirl_materialized": False,
                "exterior_heat_materialized": False,
                "mutable": False,
            },
            "truth_boundary": self.truth_boundary,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianPostPulseEtaFlattening":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected post-pulse eta-flattening schema")
        parent_payload = payload.get("parent_current_cartesian_pulse_end")
        if not isinstance(parent_payload, Mapping):
            raise ValueError("missing current pulse-end parent configuration")
        candidate = cls(
            parent=KokunoPA16CurrentCartesianPulseEndCompensated.from_configuration(
                parent_payload
            )
        )
        if _canonical_json(dict(payload)) != _canonical_json(candidate.configuration()):
            raise ValueError("serialized post-pulse source/truth binding changed")
        return candidate

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianPostPulseEtaFlattening":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "source_formulas": self.source_formulas,
            "numerical_realization": self.numerical_realization,
            "domain": {
                "lambda": self.lambda_value,
                "T_f": self.T_f,
                "log_X_pulse_end": self.log_X_pulse_end,
                "log_X_flatten_end": self.log_X_flatten_end,
            },
            "truth_boundary": self.truth_boundary,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    def flattening_report(
        self,
        eta: Any = (-1.0, -0.75, -0.25, 0.0, 0.4, 0.8, 1.0),
        *,
        radial_samples: int = 2001,
    ) -> dict[str, Any]:
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        if radial_samples < 101:
            raise ValueError("radial_samples must be at least 101")

        y = np.linspace(0.0, self.T_f, radial_samples, dtype=float)
        yy, ee = np.meshgrid(y, eta_arr, indexing="ij")
        p = self._postpulse_profile_logX(self.log_X_pulse_end + yy, ee)
        ell = np.asarray(p["ell"], dtype=float)
        end = self._postpulse_profile_logX(
            np.full(eta_arr.shape, self.log_X_flatten_end), eta_arr
        )
        end_E = np.asarray(end["E"], dtype=float)
        endpoint_relative_spread = (
            float((np.max(end_E) - np.min(end_E)) / np.max(end_E))
            if np.max(end_E) > 0.0
            else 0.0
        )
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "semantic_sha256": self.semantic_sha256,
            "eta": eta_arr.tolist(),
            "lambda_value": self.lambda_value,
            "T_f": self.T_f,
            "log_X_pulse_end": self.log_X_pulse_end,
            "log_X_flatten_end": self.log_X_flatten_end,
            "min_ell_on_deterministic_grid": float(np.min(ell)),
            "max_ell_on_deterministic_grid": float(np.max(ell)),
            "public_ell_lower_bound": -self.lambda_value - 0.1,
            "public_ell_upper_bound": -self.lambda_value,
            "endpoint_E_relative_eta_spread": endpoint_relative_spread,
            "max_abs_M_over_X_flatten_endpoint": float(
                np.max(np.abs(end["M_over_X"]))
            ),
            "max_abs_M_eta_over_X_flatten_endpoint": float(
                np.max(np.abs(end["M_eta_over_X"]))
            ),
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "T_f=100 is a frozen repository-autonomous choice; the public reconstruction supplies no unique numerical T_f",
                "the ell-window check in this receipt is deterministic numerical evidence, not a rigorous global proof",
                "U=0 is the public post-pulse schedule, not residual minimization by amplitude collapse",
                "terminal hold, relative-swirl bumps, exterior heat, matched pressure, forcing, and held-out complete NS residual remain absent",
            ],
        }


__all__ = ["KokunoPA16CurrentCartesianPostPulseEtaFlattening"]
