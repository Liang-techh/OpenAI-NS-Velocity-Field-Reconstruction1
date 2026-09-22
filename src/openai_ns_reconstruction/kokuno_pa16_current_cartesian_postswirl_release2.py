"""Current Kokuno A1 leading through the public post-swirl -1 -> -h transition.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

The corrected public outer schedule states that, after the l=-1 hold of length
4 log(1/h), l varies from -1 to -h in one logarithmic-X unit, using the same
flat C-infinity step sigma.

This increment materializes ONLY that one-unit transition.  It consumes the
exact A1 #1196 hold endpoint and the current lineage's already-bound h value;
it does not recover a hidden/source-exact numerical h.  With

    s = log(X/X_release2_start),  0 <= s <= 1,

the public schedule is realized as

    l(s) = -1 + (1-h_current) sigma(s),
    U = 0,
    D_logX log(E) = l - 1/2.

Hence

    log(E/E_start)
      = -3 s/2 + (1-h_current) int_0^s sigma(v) dv,

and F=E/sqrt(2X).  Physical M and M_eta stay constant because U=0, so
M/X, M_eta/X and v0 decay as exp(-s).

Scalar binary64 F may already be subnormal or zero at these enormous X values,
so log_F remains a first-class underflow-safe profile coordinate and Cartesian
rF is assembled from log_F.  Public velocity keeps the inherited fixed
96-significant-digit Decimal output encoding, while new profile/Cartesian
arithmetic remains binary64 before Decimal embedding.  This module therefore
does not claim end-to-end 96-digit arithmetic.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_postswirl_lminus1_hold import (
    KokunoPA16CurrentCartesianPostSwirlLMinus1Hold,
)
from .kokuno_pa16_current_cartesian_postswirl_release1 import (
    DECIMAL_DIGITS,
    _embed_decimal,
    _sigma,
    _sigma_integral,
    _sigma_prime,
)

SCHEMA = "kokuno-pa16-current-cartesian-postswirl-release2-v1"
PARENT_EXACT_HEAD = "e0cf691884fd371c477de5ae55a18cdabbd72e81"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
SOURCE_H_UPPER_BOUND = 1.0 / 100.0
RELEASE2_LENGTH = 1.0
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "coordinate_exponents": "A=1/2+h; D=1/2-h; completed construction fixes 0<h<1/100",
    "stage": "after the 4*log(1/h) l=-1 hold, vary l from -1 to -h in one unit interval using sigma",
    "coordinate": "s=log(X/X_release2_start), 0<=s<=1",
    "slope": "l(s)=-1+(1-h)*sigma(s)",
    "E_ode": "D_logX log(E)=l-1/2",
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "swirl_relation": "F=E/sqrt(2X); D_logX log(F)=l-1",
    "primitive_transport": "U=0 => M,M_eta constant => D_logX(M/X)=-(M/X)",
    "cartesian_velocity": "u_r=(v0/(2q))r; u_theta=q^(-A-1/2) r F; u3=q^(-A)U",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1196 l=-1 hold endpoint",
    "h_binding": "inherit h_current=A_current-1/2=1/2-D_current from #1196; no caller h knob",
    "sigma_integral": "reuse fixed 96-point Gauss-Legendre public-step integral from A1 #1188",
    "profile_arithmetic": "binary64 NumPy/Python arithmetic",
    "swirl_scale_representation": (
        "keep log_F as the underflow-safe scale coordinate and assemble Cartesian rF from log_F; "
        "scalar binary64 F may legitimately underflow to zero"
    ),
    "public_velocity_encoding": (
        "fixed 96-significant-decimal-digit object-array encoding inherited from the current A1 lineage; "
        "Decimal encoding is not an end-to-end 96-digit arithmetic claim"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "source_l_minus1_hold_materialized": True,
    "source_l_minus1_to_minus_h_transition_materialized": True,
    "current_cartesian_l_minus1_to_minus_h_composed": True,
    "current_h_bound_from_similarity_exponent": True,
    "source_exact_h_recovered": False,
    "overflow_safe_log_F_materialized": True,
    "decimal_output_encoding_materialized": True,
    "end_to_end_96digit_arithmetic_materialized": False,
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


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianPostSwirlRelease2:
    """Executable current leading on the public one-unit l=-1 -> -h transition."""

    parent: KokunoPA16CurrentCartesianPostSwirlLMinus1Hold = field(
        default_factory=KokunoPA16CurrentCartesianPostSwirlLMinus1Hold,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianPostSwirlLMinus1Hold):
            raise TypeError("parent must be KokunoPA16CurrentCartesianPostSwirlLMinus1Hold")
        truth = self.parent.truth_boundary
        if not truth["source_l_minus1_hold_materialized"]:
            raise ValueError("exact #1196 l=-1 hold stage is required")
        if not truth["current_cartesian_l_minus1_hold_composed"]:
            raise ValueError("exact #1196 Cartesian hold composition is required")
        h_a = self.A - 0.5
        h_d = 0.5 - self.D
        if not math.isclose(h_a, h_d, rel_tol=0.0, abs_tol=64.0 * np.finfo(float).eps):
            raise ValueError("current A,D do not encode one consistent h")
        if not (0.0 < h_a < SOURCE_H_UPPER_BOUND):
            raise ValueError("current h must satisfy the public construction range 0<h<1/100")
        if not math.isclose(h_a, self.parent.h_value, rel_tol=0.0, abs_tol=0.0):
            raise ValueError("release2 h identity drifted from exact #1196 parent")
        if self.log_radius_q1_release2_end >= math.log(np.finfo(float).max):
            raise ValueError("release2 endpoint has non-representable physical radius")

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
    def D(self) -> float:
        return float(self.split.D)

    @property
    def lambda_value(self) -> float:
        return float(self.split.lambda_value)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.split.eta_interval)

    @property
    def h_value(self) -> float:
        return self.A - 0.5

    @property
    def log_X_release2_start(self) -> float:
        return float(self.parent.log_X_hold_end)

    @property
    def log_X_release2_end(self) -> float:
        return self.log_X_release2_start + RELEASE2_LENGTH

    @property
    def log_radius_q1_release2_end(self) -> float:
        return 0.5 * (_LOG2 + self.log_X_release2_end)

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
        return self.parent.similarity_coordinates_logX(x, y, z, t)

    def _broadcast_release2(
        self, log_X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        s = log_arr - self.log_X_release2_start
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_release2_end))
        if np.any((s < -tol) | (s > RELEASE2_LENGTH + tol)):
            raise ValueError("post-swirl release2 requires 0<=log(X/X_release2_start)<=1")
        s = np.clip(s, 0.0, RELEASE2_LENGTH)
        return self.log_X_release2_start + s, eta_arr, s

    def profile_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Vectorized -1 -> -h profile with analytic log-X derivatives."""
        log_arr, eta_arr, s = self._broadcast_release2(log_X, eta)
        seam = self.parent.profile_logX(
            np.full_like(log_arr, self.log_X_release2_start, dtype=float), eta_arr
        )
        E0 = np.asarray(seam["E"], dtype=float)
        F0 = np.asarray(seam["F"], dtype=float)
        log_F0 = np.asarray(seam["log_F"], dtype=float)
        M0 = np.asarray(seam["M_over_X"], dtype=float)
        M_eta0 = np.asarray(seam["M_eta_over_X"], dtype=float)
        v00 = np.asarray(seam["v0"], dtype=float)

        sig = _sigma(s)
        sig_i = _sigma_integral(s)
        h = self.h_value
        ell = -1.0 + (1.0 - h) * sig
        log_E_ratio = -1.5 * s + (1.0 - h) * sig_i
        log_F_ratio = log_E_ratio - 0.5 * s

        E = E0 * np.exp(log_E_ratio)
        F = F0 * np.exp(log_F_ratio)
        log_F = log_F0 + log_F_ratio
        primitive_decay = np.exp(-s)
        M_ratio = M0 * primitive_decay
        M_eta_ratio = M_eta0 * primitive_decay
        v0 = v00 * primitive_decay
        U = np.zeros_like(E)

        DlogX_E = (ell - 0.5) * E
        DlogX_F = (ell - 1.0) * F
        DlogX_log_F = ell - 1.0
        DlogX_U = np.zeros_like(E)
        DlogX_M_ratio = -M_ratio
        DlogX_M_eta_ratio = -M_eta_ratio
        DlogX_v0 = -v0

        arrays = (
            E,
            F,
            log_F,
            M_ratio,
            M_eta_ratio,
            v0,
            ell,
            DlogX_E,
            DlogX_F,
            DlogX_log_F,
            DlogX_M_ratio,
            DlogX_M_eta_ratio,
            DlogX_v0,
        )
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("post-swirl release2 profile became non-finite")
        if np.any(E <= 0.0) or np.any(F < 0.0):
            raise RuntimeError("post-swirl release2 lost positive E/nonnegative F")

        return {
            "log_X": log_arr,
            "eta": eta_arr,
            "release2_s": s,
            "h_current": np.full_like(E, h),
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
            "DlogX_log_F": DlogX_log_F,
            "DlogX_U": DlogX_U,
            "DlogX_M_over_X": DlogX_M_ratio,
            "DlogX_M_eta_over_X": DlogX_M_eta_ratio,
            "DlogX_v0": DlogX_v0,
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return release2 Cartesian velocity with inherited Decimal encoding."""
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
            inherited = self.parent.velocity(
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
            raise RuntimeError("release2 Cartesian velocity became non-finite")

        out = np.empty(out_float.shape, dtype=object)
        for idx in np.ndindex(out_float.shape):
            out[idx] = _embed_decimal(out_float[idx])
        return out.reshape(shape + (3,))

    def release2_report(self) -> dict[str, Any]:
        eta = np.asarray([-0.5, 0.0, 0.5], dtype=float)
        p0 = self.profile_logX(self.log_X_release2_start, eta)
        p1 = self.profile_logX(self.log_X_release2_end, eta)
        e_ratio = np.asarray(p1["E"], dtype=float) / np.asarray(p0["E"], dtype=float)
        f_ratio_log_chart = np.exp(
            np.asarray(p1["log_F"], dtype=float) - np.asarray(p0["log_F"], dtype=float)
        )
        m0 = np.asarray(p0["M_over_X"], dtype=float)
        m1 = np.asarray(p1["M_over_X"], dtype=float)
        nz = np.abs(m0) > 0.0
        m_ratio = m1[nz] / m0[nz] if np.any(nz) else np.asarray([], dtype=float)

        seam_radius = math.exp(0.5 * (_LOG2 + self.log_X_release2_start))
        parent_v = self.parent.velocity(seam_radius, 0.0, 0.0, 0.0)
        child_v = self.velocity(seam_radius, 0.0, 0.0, 0.0)
        seam_exact = bool(
            all(a == b for a, b in zip(parent_v.reshape(-1), child_v.reshape(-1)))
        )

        h = self.h_value
        expected_e = math.exp(-1.0 - 0.5 * h)
        expected_f = math.exp(-1.5 - 0.5 * h)
        expected_primitive = math.exp(-1.0)

        return {
            "schema": "kokuno-agent1-postswirl-release2-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "A": self.A,
            "D": self.D,
            "h_current": h,
            "h_source_upper_bound": SOURCE_H_UPPER_BOUND,
            "release2_length": RELEASE2_LENGTH,
            "ell_start": float(np.asarray(p0["ell"])[0]),
            "ell_end": float(np.asarray(p1["ell"])[0]),
            "sigma_integral_end": float(np.asarray(p1["sigma_integral"])[0]),
            "E_ratio_end_min": float(np.min(e_ratio)),
            "E_ratio_end_max": float(np.max(e_ratio)),
            "F_ratio_end_min": float(np.min(f_ratio_log_chart)),
            "F_ratio_end_max": float(np.max(f_ratio_log_chart)),
            "M_over_X_ratio_end_min": float(np.min(m_ratio)) if m_ratio.size else None,
            "M_over_X_ratio_end_max": float(np.max(m_ratio)) if m_ratio.size else None,
            "expected_E_ratio": expected_e,
            "expected_F_ratio": expected_f,
            "expected_primitive_ratio": expected_primitive,
            "parent_child_decimal_seam_exact": seam_exact,
            "decimal_digits": DECIMAL_DIGITS,
            "end_to_end_96digit_arithmetic_materialized": False,
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
            "A": self.A,
            "D": self.D,
            "h_current_from_A_minus_half": self.h_value,
            "source_h_upper_bound": SOURCE_H_UPPER_BOUND,
            "release2_length": RELEASE2_LENGTH,
            "decimal_digits": DECIMAL_DIGITS,
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
