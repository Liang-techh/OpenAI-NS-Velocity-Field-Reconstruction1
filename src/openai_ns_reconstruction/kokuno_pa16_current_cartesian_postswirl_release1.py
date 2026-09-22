"""First public post-relative-swirl exterior transition for the current Kokuno A1 lineage.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

After the post-pulse eta flattening, the 30 log(1/lambda) hold, and the two
relative-swirl bumps, the corrected public schedule says:

    Starting at X_rel, vary l from -lambda to -1 in a unit interval,
    hold l=-1 for 4 log(1/h), and vary from -1 to -h in a unit interval,
    all using sigma.

This increment implements ONLY the first unit interval.  With
s = log(X/X_rel) in [0,1],

    l(s) = -lambda - (1-lambda) sigma(s),
    U = 0,

and d_s log E = l - 1/2.  The exact current endpoint from A1 #1179 is
retained; no hidden source coefficient is inferred and no later hold/release,
terminal multiplier, exterior heat replacement, pressure, forcing or NS
validation is included.

The parent #1179 uses a fixed 96-significant-digit Decimal total only because
the preceding relative-swirl edit is below binary64 relative resolution.  The
two compact bumps are already zero at X_rel, so this stage has no hidden
sub-epsilon additive correction.  Nevertheless this class keeps the same
Decimal-valued public velocity representation at the seam and throughout the
new interval, avoiding an unannounced representation downgrade.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_relative_swirl_decimal_composed import (
    DECIMAL_DIGITS,
    KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed,
)

SCHEMA = "kokuno-pa16-current-cartesian-postswirl-release1-v1"
PARENT_EXACT_HEAD = "8c5c5b6d55a68de8285ed1dd0ebb55128f3078a4"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
RELEASE1_LENGTH = 1.0
SIGMA_QUADRATURE_ORDER = 96

_SOURCE_FORMULAS = {
    "stage": "starting at X_rel, vary l from -lambda to -1 in a unit interval using sigma",
    "coordinate": "s=log(X/X_rel), 0<=s<=1",
    "slope": "l(s)=-lambda-(1-lambda)*sigma(s)",
    "E_ode": "D_logX log(E)=l-1/2",
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "swirl_relation": "F=E/sqrt(2X)",
    "primitive_transport": "U=0 => M,M_eta constant => D_logX(M/X)=-(M/X)",
    "cartesian_velocity": "u_r=(v0/(2q))r; u_theta=q^(-A-1/2) r F; u3=q^(-A)U",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1179 precision-qualified relative-swirl endpoint",
    "sigma_integral": "fixed 96-point Gauss-Legendre on [0,s]",
    "public_velocity_representation": (
        "fixed 96-significant-decimal-digit object array, ROUND_HALF_EVEN; "
        "precision inherited from #1179 and not caller-tunable"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "post_relative_swirl_release1_materialized": True,
    "post_relative_swirl_release1_public_velocity_materialized": True,
    "post_relative_swirl_release1_decimal_representation_materialized": True,
    "source_release1_exact_schedule_materialized": True,
    "source_l_minus1_hold_materialized": False,
    "source_l_minus1_to_minus_h_transition_materialized": False,
    "source_terminal_multiplier_materialized": False,
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


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _configure_context(ctx) -> None:
    ctx.prec = DECIMAL_DIGITS
    ctx.rounding = ROUND_HALF_EVEN
    ctx.Emax = 999999999
    ctx.Emin = -999999999


def _embed_decimal(value: float) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("only finite binary64 values may enter Decimal stage representation")
    with localcontext() as ctx:
        _configure_context(ctx)
        return +Decimal.from_float(value)


def _sigma(value: Any) -> np.ndarray:
    """Public flat C-infinity step, evaluated stably in binary64."""
    y = _finite(value, "sigma argument")
    out = np.zeros_like(y, dtype=float)
    out[y >= 1.0] = 1.0
    inside = (y > 0.0) & (y < 1.0)
    if np.any(inside):
        z = y[inside]
        logit = -1.0 / (z * z) + 1.0 / ((1.0 - z) ** 2)
        vals = np.empty_like(logit)
        pos = logit >= 0.0
        if np.any(pos):
            e = np.exp(-logit[pos])
            vals[pos] = 1.0 / (1.0 + e)
        if np.any(~pos):
            e = np.exp(logit[~pos])
            vals[~pos] = e / (1.0 + e)
        out[inside] = vals
    return out


def _sigma_prime(value: Any) -> np.ndarray:
    y = _finite(value, "sigma argument")
    out = np.zeros_like(y, dtype=float)
    inside = (y > 0.0) & (y < 1.0)
    if np.any(inside):
        z = y[inside]
        sig = _sigma(z)
        logit_prime = 2.0 / (z ** 3) + 2.0 / ((1.0 - z) ** 3)
        out[inside] = sig * (1.0 - sig) * logit_prime
    return out


_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(SIGMA_QUADRATURE_ORDER)


def _sigma_integral(value: Any) -> np.ndarray:
    """Integral int_0^s sigma(v) dv for s in [0,1]."""
    s = _finite(value, "release coordinate")
    tol = 64.0 * np.finfo(float).eps
    if np.any((s < -tol) | (s > 1.0 + tol)):
        raise ValueError("release coordinate must lie in [0,1]")
    sc = np.clip(s, 0.0, 1.0)
    flat = sc.reshape(-1)
    out = np.empty_like(flat)
    for i, upper in enumerate(flat):
        if upper == 0.0:
            out[i] = 0.0
            continue
        x = 0.5 * upper * (_GL_NODES + 1.0)
        out[i] = 0.5 * upper * float(np.dot(_GL_WEIGHTS, _sigma(x)))
    return out.reshape(sc.shape)


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianPostSwirlRelease1:
    """Current Cartesian leading on the first public post-swirl release interval."""

    parent: KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed = field(
        default_factory=KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(
            self.parent, KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed
        ):
            raise TypeError(
                "parent must be KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed"
            )
        truth = self.parent.truth_boundary
        if not truth["current_cartesian_relative_swirl_composed"]:
            raise ValueError("exact #1179 relative-swirl total is required")
        if not truth["precision_qualified_public_velocity_materialized"]:
            raise ValueError("exact #1179 precision-qualified public velocity is required")
        seam = self.parent.split.split_profile_logX(
            self.log_X_rel, np.asarray([-0.5, 0.0, 0.5], dtype=float)
        )
        if np.any(np.asarray(seam["relative_edit"], dtype=float) != 0.0):
            raise ValueError("relative-swirl edit is not exactly zero at X_rel")

    @property
    def split(self):
        return self.parent.split

    @property
    def geometry(self):
        return self.parent.geometry

    @property
    def A(self) -> float:
        return float(self.split.A)

    @property
    def lambda_value(self) -> float:
        return float(self.split.lambda_value)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.split.eta_interval)

    @property
    def log_X_rel(self) -> float:
        return float(self.split.log_X_hold_end)

    @property
    def log_X_release1_end(self) -> float:
        return self.log_X_rel + RELEASE1_LENGTH

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

    def similarity_coordinates_logX(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.split.similarity_coordinates_logX(x, y, z, t)

    def _broadcast_release(
        self, log_X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        s = log_arr - self.log_X_rel
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_release1_end))
        if np.any((s < -tol) | (s > RELEASE1_LENGTH + tol)):
            raise ValueError("post-swirl release1 requires 0<=log(X/X_rel)<=1")
        s = np.clip(s, 0.0, RELEASE1_LENGTH)
        return self.log_X_rel + s, eta_arr, s

    def profile_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Vectorized low-dimensional release profile and analytic log-X derivatives."""
        log_arr, eta_arr, s = self._broadcast_release(log_X, eta)
        seam = self.split.split_profile_logX(
            np.full_like(log_arr, self.log_X_rel, dtype=float), eta_arr
        )
        E_rel = np.asarray(seam["E_base"], dtype=float)
        F_rel = np.asarray(seam["F_base"], dtype=float)
        log_F_rel = np.asarray(seam["log_F_base"], dtype=float)
        M_rel = np.asarray(seam["M_over_X_base"], dtype=float)
        M_eta_rel = np.asarray(seam["M_eta_over_X_base"], dtype=float)
        v0_rel = np.asarray(seam["v0_base"], dtype=float)

        sig = _sigma(s)
        sig_i = _sigma_integral(s)
        lam = self.lambda_value
        ell = -lam - (1.0 - lam) * sig
        log_E_ratio = -(0.5 + lam) * s - (1.0 - lam) * sig_i
        E = E_rel * np.exp(log_E_ratio)
        log_F = log_F_rel + log_E_ratio - 0.5 * s
        F = F_rel * np.exp(log_E_ratio - 0.5 * s)
        decay = np.exp(-s)
        M_ratio = M_rel * decay
        M_eta_ratio = M_eta_rel * decay
        v0 = v0_rel * decay
        U = np.zeros_like(E)

        DlogX_E = (ell - 0.5) * E
        DlogX_F = (ell - 1.0) * F
        DlogX_U = np.zeros_like(E)
        DlogX_M_ratio = -M_ratio
        DlogX_M_eta_ratio = -M_eta_ratio
        DlogX_v0 = -v0

        arrays = (
            E, F, log_F, M_ratio, M_eta_ratio, v0, ell,
            DlogX_E, DlogX_F, DlogX_M_ratio, DlogX_M_eta_ratio, DlogX_v0,
        )
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("post-swirl release1 profile became non-finite")

        return {
            "log_X": log_arr,
            "eta": eta_arr,
            "release_s": s,
            "sigma": sig,
            "sigma_prime": _sigma_prime(s),
            "sigma_integral": sig_i,
            "ell": ell,
            "E": E,
            "F": F,
            "log_F": log_F,
            "U": U,
            "M_over_X": M_ratio,
            "M_eta_over_X": M_eta_ratio,
            "v0": v0,
            "DlogX_E": DlogX_E,
            "DlogX_F": DlogX_F,
            "DlogX_U": DlogX_U,
            "DlogX_M_over_X": DlogX_M_ratio,
            "DlogX_M_eta_over_X": DlogX_M_eta_ratio,
            "DlogX_v0": DlogX_v0,
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return the release-stage Cartesian velocity as fixed-precision Decimals."""
        coords = self.similarity_coordinates_logX(x, y, z, t)
        xb = np.asarray(coords["x"], dtype=float)
        yb = np.asarray(coords["y"], dtype=float)
        zb = np.asarray(coords["z"], dtype=float)
        tb = np.asarray(coords["t"], dtype=float)
        q = np.asarray(coords["q"], dtype=float)
        eta = np.asarray(coords["eta"], dtype=float)
        radius = np.asarray(coords["radius"], dtype=float)
        log_X = np.asarray(coords["log_X"], dtype=float)
        shape = radius.shape

        xb_f, yb_f = xb.reshape(-1), yb.reshape(-1)
        zb_f, tb_f = zb.reshape(-1), tb.reshape(-1)
        q_f, eta_f = q.reshape(-1), eta.reshape(-1)
        r_f, lx_f = radius.reshape(-1), log_X.reshape(-1)
        out_float = np.zeros((r_f.size, 3), dtype=float)

        axis = r_f == 0.0
        if np.any(axis):
            inherited = self.parent.split.parent.velocity(
                xb_f[axis], yb_f[axis], zb_f[axis], tb_f[axis]
            )
            out_float[axis] = np.asarray(inherited, dtype=float).reshape((-1, 3))

        nonaxis = ~axis
        if np.any(nonaxis):
            p = self.profile_logX(lx_f[nonaxis], eta_f[nonaxis])
            r = r_f[nonaxis]
            xx = xb_f[nonaxis]
            yy = yb_f[nonaxis]
            qq = q_f[nonaxis]
            if np.any(qq <= 0.0):
                raise ValueError("similarity q must stay positive")
            v0 = np.asarray(p["v0"], dtype=float)
            radial_scale = v0 / (2.0 * qq)
            utheta = np.exp(
                -(self.A + 0.5) * np.log(qq)
                + np.asarray(p["log_F"], dtype=float)
                + np.log(r)
            )
            tx = -yy / r
            ty = xx / r
            out_float[nonaxis, 0] = radial_scale * xx + tx * utheta
            out_float[nonaxis, 1] = radial_scale * yy + ty * utheta
            out_float[nonaxis, 2] = 0.0

        if np.any(~np.isfinite(out_float)):
            raise RuntimeError("post-swirl release1 Cartesian velocity became non-finite")

        out = np.empty(out_float.shape, dtype=object)
        for idx in np.ndindex(out_float.shape):
            out[idx] = _embed_decimal(out_float[idx])
        return out.reshape(shape + (3,))

    def release_report(self) -> dict[str, Any]:
        eta = np.asarray([-0.5, 0.0, 0.5], dtype=float)
        p0 = self.profile_logX(self.log_X_rel, eta)
        p1 = self.profile_logX(self.log_X_release1_end, eta)
        e_ratio = np.asarray(p1["E"], dtype=float) / np.asarray(p0["E"], dtype=float)
        f_ratio = np.asarray(p1["F"], dtype=float) / np.asarray(p0["F"], dtype=float)
        seam_radius = math.exp(0.5 * (math.log(2.0) + self.log_X_rel))
        parent_v = self.parent.velocity(seam_radius, 0.0, 0.0, 0.0)
        child_v = self.velocity(seam_radius, 0.0, 0.0, 0.0)
        seam_exact = bool(
            all(a == b for a, b in zip(parent_v.reshape(-1), child_v.reshape(-1)))
        )
        return {
            "schema": "kokuno-agent1-postswirl-release1-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "lambda": self.lambda_value,
            "release_length": RELEASE1_LENGTH,
            "sigma_integral_end": float(_sigma_integral(np.asarray(1.0))),
            "ell_start": float(np.asarray(p0["ell"])[0]),
            "ell_end": float(np.asarray(p1["ell"])[0]),
            "E_ratio_end_min": float(np.min(e_ratio)),
            "E_ratio_end_max": float(np.max(e_ratio)),
            "F_ratio_end_min": float(np.min(f_ratio)),
            "F_ratio_end_max": float(np.max(f_ratio)),
            "parent_child_decimal_seam_exact": seam_exact,
            "truth_boundary": self.truth_boundary,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_blob": SOURCE_BLOB,
            "source_path": SOURCE_PATH,
            "source_release": SOURCE_RELEASE,
            "source_release_date": SOURCE_RELEASE_DATE,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "decimal_digits": DECIMAL_DIGITS,
            "release1_length": RELEASE1_LENGTH,
            "sigma_quadrature_order": SIGMA_QUADRATURE_ORDER,
            "source_formulas": self.source_formulas,
            "numerical_realization": self.numerical_realization,
            "truth_boundary": self.truth_boundary,
        }

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]):
        obj = cls()
        if dict(payload) != obj.configuration():
            raise ValueError("configuration/provenance drift")
        return obj

    @classmethod
    def load_configuration(cls, path: str | Path):
        return cls.from_configuration(json.loads(Path(path).read_text()))
