"""Current Kokuno A1 leading through the public post-swirl l=-1 hold.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

The public coordinate map fixes

    A = 1/2 + h,  D = 1/2 - h,

with the completed construction choosing 0 < h < 1/100.  After the first
post-relative-swirl unit transition (A1 #1188), the public outer schedule says
to hold

    l = D_logX log(H) = -1

for 4 log(1/h), before a later one-unit transition from -1 to -h.

This increment materializes ONLY that constant-l hold.  It does not invent a
new h: the current candidate value is bound from its already frozen similarity
exponent by h_current = A - 1/2 and cross-checked against D.  This is a
repository candidate binding of a public parameter relation, not recovery of a
hidden/source-exact numerical h.

With s=log(X/X_hold_start), 0<=s<=4 log(1/h_current), public post-pulse U=0
and l=-1 give

    E = E_start exp(-3 s/2),
    F = F_start exp(-2 s),
    U = 0.

Physical M and M_eta remain constant, hence M/X, M_eta/X and v0 decay as
exp(-s).  The exact #1188 current primitive history is therefore carried rather
than reset.  The public velocity keeps #1188's Decimal output encoding, while
all new profile/Cartesian arithmetic remains binary64 before Decimal embedding;
this module does NOT claim end-to-end 96-digit evaluation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_postswirl_release1 import (
    DECIMAL_DIGITS,
    KokunoPA16CurrentCartesianPostSwirlRelease1,
    _embed_decimal,
)

SCHEMA = "kokuno-pa16-current-cartesian-postswirl-lminus1-hold-v1"
PARENT_EXACT_HEAD = "f5b5b8b532e0802e181bfead5e748cd4008906bf"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
SOURCE_H_UPPER_BOUND = 1.0 / 100.0
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "coordinate_exponents": "A=1/2+h; D=1/2-h; completed construction fixes 0<h<1/100",
    "stage": "after release1 hold l=-1 for 4*log(1/h)",
    "coordinate": "s=log(X/X_hold_start), 0<=s<=4*log(1/h)",
    "slope": "l=D_logX log(H)=-1",
    "E_ode": "D_logX log(E)=l-1/2=-3/2",
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "swirl_relation": "F=E/sqrt(2X)",
    "primitive_transport": "U=0 => M,M_eta constant => D_logX(M/X)=-(M/X)",
    "cartesian_velocity": "u_r=(v0/(2q))r; u_theta=q^(-A-1/2) r F; u3=q^(-A)U",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1188 first-release endpoint",
    "h_binding": "h_current=A_current-1/2, cross-checked with 1/2-D_current; no caller h knob",
    "profile_arithmetic": "binary64 NumPy/Python arithmetic",
    "public_velocity_encoding": (
        "fixed 96-significant-decimal-digit object-array encoding inherited from #1188; "
        "Decimal encoding is not an end-to-end 96-digit arithmetic claim"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "post_relative_swirl_release1_materialized": True,
    "source_l_minus1_hold_materialized": True,
    "current_cartesian_l_minus1_hold_composed": True,
    "current_h_bound_from_similarity_exponent": True,
    "source_exact_h_recovered": False,
    "decimal_output_encoding_materialized": True,
    "end_to_end_96digit_arithmetic_materialized": False,
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


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianPostSwirlLMinus1Hold:
    """Executable current leading on the public constant-l=-1 hold."""

    parent: KokunoPA16CurrentCartesianPostSwirlRelease1 = field(
        default_factory=KokunoPA16CurrentCartesianPostSwirlRelease1,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianPostSwirlRelease1):
            raise TypeError("parent must be KokunoPA16CurrentCartesianPostSwirlRelease1")
        truth = self.parent.truth_boundary
        if not truth["post_relative_swirl_release1_materialized"]:
            raise ValueError("exact #1188 first-release stage is required")
        h_a = self.A - 0.5
        h_d = 0.5 - self.D
        if not math.isclose(h_a, h_d, rel_tol=0.0, abs_tol=64.0 * np.finfo(float).eps):
            raise ValueError("current A,D do not encode one consistent h")
        if not (0.0 < h_a < SOURCE_H_UPPER_BOUND):
            raise ValueError("current h must satisfy the public construction range 0<h<1/100")
        if not math.isfinite(self.hold_length) or self.hold_length <= 0.0:
            raise ValueError("public l=-1 hold length must be finite and positive")
        if self.log_radius_q1_hold_end >= math.log(np.finfo(float).max):
            raise ValueError("l=-1 hold endpoint has non-representable physical radius")

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
    def log_X_hold_start(self) -> float:
        return float(self.parent.log_X_release1_end)

    @property
    def hold_length(self) -> float:
        return 4.0 * math.log(1.0 / self.h_value)

    @property
    def log_X_hold_end(self) -> float:
        return self.log_X_hold_start + self.hold_length

    @property
    def log_radius_q1_hold_end(self) -> float:
        return 0.5 * (_LOG2 + self.log_X_hold_end)

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

    def _broadcast_hold(
        self, log_X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        s = log_arr - self.log_X_hold_start
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_hold_end))
        if np.any((s < -tol) | (s > self.hold_length + tol)):
            raise ValueError("l=-1 hold requires 0<=log(X/X_hold_start)<=4*log(1/h)")
        s = np.clip(s, 0.0, self.hold_length)
        return self.log_X_hold_start + s, eta_arr, s

    def profile_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Vectorized l=-1 profile with analytic log-X derivatives."""
        log_arr, eta_arr, s = self._broadcast_hold(log_X, eta)
        seam = self.parent.profile_logX(
            np.full_like(log_arr, self.log_X_hold_start, dtype=float), eta_arr
        )
        E0 = np.asarray(seam["E"], dtype=float)
        F0 = np.asarray(seam["F"], dtype=float)
        log_F0 = np.asarray(seam["log_F"], dtype=float)
        M0 = np.asarray(seam["M_over_X"], dtype=float)
        M_eta0 = np.asarray(seam["M_eta_over_X"], dtype=float)
        v00 = np.asarray(seam["v0"], dtype=float)

        E_decay = np.exp(-1.5 * s)
        F_decay = np.exp(-2.0 * s)
        primitive_decay = np.exp(-s)
        E = E0 * E_decay
        F = F0 * F_decay
        log_F = log_F0 - 2.0 * s
        M_ratio = M0 * primitive_decay
        M_eta_ratio = M_eta0 * primitive_decay
        v0 = v00 * primitive_decay
        U = np.zeros_like(E)
        ell = np.full_like(E, -1.0)

        DlogX_E = -1.5 * E
        DlogX_F = -2.0 * F
        DlogX_U = np.zeros_like(E)
        DlogX_M_ratio = -M_ratio
        DlogX_M_eta_ratio = -M_eta_ratio
        DlogX_v0 = -v0

        arrays = (
            E, F, log_F, M_ratio, M_eta_ratio, v0,
            DlogX_E, DlogX_F, DlogX_M_ratio, DlogX_M_eta_ratio, DlogX_v0,
        )
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("l=-1 hold profile became non-finite")
        if np.any(E <= 0.0) or np.any(F < 0.0):
            raise RuntimeError("l=-1 hold lost positive E/nonnegative F")

        return {
            "log_X": log_arr,
            "eta": eta_arr,
            "hold_s": s,
            "h_current": np.full_like(E, self.h_value),
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
        """Return the hold-stage Cartesian velocity with inherited Decimal encoding."""
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
            raise RuntimeError("l=-1 hold Cartesian velocity became non-finite")

        out = np.empty(out_float.shape, dtype=object)
        for idx in np.ndindex(out_float.shape):
            out[idx] = _embed_decimal(out_float[idx])
        return out.reshape(shape + (3,))

    def hold_report(self) -> dict[str, Any]:
        eta = np.asarray([-0.5, 0.0, 0.5], dtype=float)
        p0 = self.profile_logX(self.log_X_hold_start, eta)
        p1 = self.profile_logX(self.log_X_hold_end, eta)
        e_ratio = np.asarray(p1["E"], dtype=float) / np.asarray(p0["E"], dtype=float)
        f_ratio = np.asarray(p1["F"], dtype=float) / np.asarray(p0["F"], dtype=float)
        m0 = np.asarray(p0["M_over_X"], dtype=float)
        m1 = np.asarray(p1["M_over_X"], dtype=float)
        nz = np.abs(m0) > 0.0
        m_ratio = m1[nz] / m0[nz] if np.any(nz) else np.asarray([], dtype=float)
        seam_radius = math.exp(0.5 * (_LOG2 + self.log_X_hold_start))
        parent_v = self.parent.velocity(seam_radius, 0.0, 0.0, 0.0)
        child_v = self.velocity(seam_radius, 0.0, 0.0, 0.0)
        seam_exact = bool(
            all(a == b for a, b in zip(parent_v.reshape(-1), child_v.reshape(-1)))
        )
        return {
            "schema": "kokuno-agent1-postswirl-lminus1-hold-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "A": self.A,
            "D": self.D,
            "h_current": self.h_value,
            "h_source_upper_bound": SOURCE_H_UPPER_BOUND,
            "hold_length": self.hold_length,
            "ell": -1.0,
            "E_ratio_end_min": float(np.min(e_ratio)),
            "E_ratio_end_max": float(np.max(e_ratio)),
            "F_ratio_end_min": float(np.min(f_ratio)),
            "F_ratio_end_max": float(np.max(f_ratio)),
            "M_over_X_ratio_end_min": float(np.min(m_ratio)) if m_ratio.size else None,
            "M_over_X_ratio_end_max": float(np.max(m_ratio)) if m_ratio.size else None,
            "expected_E_ratio_h6": self.h_value ** 6,
            "expected_F_ratio_h8": self.h_value ** 8,
            "expected_primitive_ratio_h4": self.h_value ** 4,
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
            "hold_length": self.hold_length,
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
