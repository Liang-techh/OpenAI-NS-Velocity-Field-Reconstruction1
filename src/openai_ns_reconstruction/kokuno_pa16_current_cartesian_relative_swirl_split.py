"""Compensated relative-swirl correction channel for the current Kokuno A1 lineage.

Pinned public provenance is KokunoYumeto's corrected 2026-09-09 reconstruction
at commit 143f6773feb424ad9ed3a8d116653200f20346b7.  After eta flattening and
an eta-independent hold, the public construction makes the relative edit

    E = E_unedited * (1 + c1 beta1 + c2 beta2),

with two smooth bumps.  The coefficients are fixed by the public angular
cumulative endpoint condition and exact pressure neutrality.  A1 #1168
materializes the absolute current angular target and feeds it into the #1154
solver.

For the current lambda=0.05 lineage the row-scaled target is O(lambda^28).
Consequently c1,c2 and the pointwise relative edit are far below binary64
relative precision.  Forming ``1 + edit`` and claiming that ordinary float64
Cartesian velocity now contains the correction would silently erase the new
field.  This increment therefore makes the mathematically nonzero correction
executable as a *split* channel:

    base profile/velocity + delta profile/velocity.

The delta is itself representable in float64 and carries analytic log-X and
eta jets.  The ordinary current Cartesian candidate is intentionally not
promoted yet; a later integration layer must use a precision/representation
that can retain the split sum.  This is not PDE validation and does not claim
source-exact hidden bump data.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_current_relative_swirl_absolute_target import (
    KokunoCurrentRelativeSwirlAbsoluteTarget,
)

SCHEMA = "kokuno-pa16-current-cartesian-relative-swirl-split-v1"
PARENT_EXACT_HEAD = "bb246140c6c2451b22b59aad9ba90604bc2acb21"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "relative_edit": "E=E_unedited*(1+c1*beta1+c2*beta2)",
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "angular_target": "I=XH/(1-lambda) after the two relative-swirl bumps",
    "pressure_neutrality": "integral(E^2-E_unedited^2)dy=0",
    "cartesian_velocity": (
        "u_r=(v0/(2q))r; u_theta=q^(-A-1/2) r F; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1168 absolute current target and #1154 solver",
    "pointwise_bump_shape": (
        "repository-autonomous standard C-infinity bump already frozen by #1154"
    ),
    "representation": (
        "preserve base and relative-swirl delta as separate float64 channels; "
        "do not round the O(lambda^28) edit away in 1+edit"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "public_relative_swirl_two_bump_algebra_materialized": True,
    "absolute_current_r_I_at_flattening_endpoint_materialized": True,
    "current_relative_swirl_total_target_materialized": True,
    "relative_swirl_profile_delta_channel_materialized": True,
    "relative_swirl_cartesian_delta_channel_materialized": True,
    "binary64_total_relative_swirl_sum_is_resolved": False,
    "current_cartesian_relative_swirl_composed": False,
    "source_terminal_hold_after_eta_flattening_materialized": False,
    "source_exact_pointwise_bump_shape_recovered": False,
    "source_hidden_parameters_recovered": False,
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


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianRelativeSwirlSplit:
    """Executable split base+delta representation of the late relative-swirl edit."""

    target: KokunoCurrentRelativeSwirlAbsoluteTarget = field(
        default_factory=KokunoCurrentRelativeSwirlAbsoluteTarget,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.target, KokunoCurrentRelativeSwirlAbsoluteTarget):
            raise TypeError("target must be KokunoCurrentRelativeSwirlAbsoluteTarget")
        truth = self.target.truth_boundary
        if not truth["current_relative_swirl_total_target_materialized"]:
            raise ValueError("exact #1168 current angular target is required")
        if truth["current_cartesian_relative_swirl_composed"]:
            raise ValueError("parent unexpectedly already composes relative swirl")
        if abs(self.parent.hold_length - self.compensator.hold_length) > 2e-13:
            raise ValueError("hold/compensator geometry mismatch")
        if self.compensator.width <= 0.0:
            raise ValueError("relative-swirl bump width must be positive")

    @property
    def parent(self):
        return self.target.hold

    @property
    def compensator(self):
        return self.target.compensator

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
    def log_X_flatten_end(self) -> float:
        return float(self.parent.log_X_flatten_end)

    @property
    def log_X_hold_end(self) -> float:
        return float(self.parent.log_X_hold_end)

    @property
    def hold_length(self) -> float:
        return float(self.parent.hold_length)

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

    def _broadcast_hold(self, log_X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        s = log_arr - self.log_X_flatten_end
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_hold_end))
        if np.any((s < -tol) | (s > self.hold_length + tol)):
            raise ValueError("relative-swirl split requires the public eta-independent hold")
        return np.clip(log_arr, self.log_X_flatten_end, self.log_X_hold_end), eta_arr

    def _coefficients(self, eta: np.ndarray):
        target = self.target.total_target(eta)
        return self.compensator.solve(target)

    def _coefficient_eta_jet(self, eta: np.ndarray):
        target, target_eta = self.target.target_with_eta(eta)
        return self.compensator.solve_target_jet(target, target_eta)

    def split_profile_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return unedited base and representable relative-swirl delta separately."""
        log_arr, eta_arr = self._broadcast_hold(log_X, eta)
        base = self.parent.unedited_hold_profile_logX(log_arr, eta_arr)
        s = np.asarray(base["hold_s"], dtype=float)
        sol = self._coefficients(eta_arr)
        c1 = np.asarray(sol.c1, dtype=float)
        c2 = np.asarray(sol.c2, dtype=float)
        edit = np.asarray(self.compensator.relative_edit(s, c1, c2), dtype=float)
        edit_y = np.asarray(self.compensator.relative_edit_y(s, c1, c2), dtype=float)

        E0 = np.asarray(base["E_unedited"], dtype=float)
        F0 = np.asarray(base["F_unedited"], dtype=float)
        log_F0 = np.asarray(base["log_F_unedited"], dtype=float)
        M_ratio = np.asarray(base["M_over_X_unedited"], dtype=float)
        M_eta_ratio = np.asarray(base["M_eta_over_X_unedited"], dtype=float)
        v0 = np.asarray(base["v0_unedited"], dtype=float)

        delta_E = E0 * edit
        delta_F = F0 * edit
        delta_log_F_first_order = edit
        delta_E_DlogX = E0 * (
            -(0.5 + self.lambda_value) * edit + edit_y
        )
        delta_F_DlogX = F0 * (-(1.0 + self.lambda_value) * edit + edit_y)

        arrays = (
            c1, c2, edit, edit_y, delta_E, delta_F,
            delta_E_DlogX, delta_F_DlogX, M_ratio, M_eta_ratio, v0,
        )
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("relative-swirl split profile became non-finite")

        return {
            "log_X": log_arr,
            "eta": eta_arr,
            "hold_s": s,
            "c1": c1,
            "c2": c2,
            "relative_edit": edit,
            "relative_edit_DlogX": edit_y,
            "E_base": E0,
            "F_base": F0,
            "log_F_base": log_F0,
            "U_base": np.asarray(base["U_unedited"], dtype=float),
            "M_over_X_base": M_ratio,
            "M_eta_over_X_base": M_eta_ratio,
            "v0_base": v0,
            "delta_E": delta_E,
            "delta_F": delta_F,
            "delta_log_F_first_order": delta_log_F_first_order,
            "delta_U": np.zeros_like(delta_E),
            "delta_M_over_X": np.zeros_like(delta_E),
            "delta_M_eta_over_X": np.zeros_like(delta_E),
            "delta_v0": np.zeros_like(delta_E),
            "delta_E_DlogX": delta_E_DlogX,
            "delta_F_DlogX": delta_F_DlogX,
            "delta_U_DlogX": np.zeros_like(delta_E),
        }

    def split_profile_eta_jet_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_hold(log_X, eta)
        base = self.parent.unedited_hold_profile_logX(log_arr, eta_arr)
        s = np.asarray(base["hold_s"], dtype=float)
        sol, c1_eta, c2_eta = self._coefficient_eta_jet(eta_arr)
        c1 = np.asarray(sol.c1, dtype=float)
        c2 = np.asarray(sol.c2, dtype=float)
        c1_eta = np.asarray(c1_eta, dtype=float)
        c2_eta = np.asarray(c2_eta, dtype=float)
        beta1 = np.asarray(self.compensator.beta1(s), dtype=float)
        beta2 = np.asarray(self.compensator.beta2(s), dtype=float)
        edit = np.asarray(self.compensator.relative_edit(s, c1, c2), dtype=float)
        edit_eta = c1_eta * beta1 + c2_eta * beta2
        E0 = np.asarray(base["E_unedited"], dtype=float)
        F0 = np.asarray(base["F_unedited"], dtype=float)
        return {
            "relative_edit": edit,
            "relative_edit_eta": edit_eta,
            "c1_eta": c1_eta,
            "c2_eta": c2_eta,
            "delta_E_eta": E0 * edit_eta,
            "delta_F_eta": F0 * edit_eta,
            "delta_U_eta": np.zeros_like(edit_eta),
        }

    def velocity_split(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return (unedited_base_velocity, relative_swirl_delta_velocity).

        The two arrays have identical ``(...,3)`` shape.  They are deliberately
        not added in float64 because the relative correction is below ordinary
        relative precision for the current lineage.
        """
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
        base_f = np.zeros((r_f.size, 3), dtype=float)
        delta_f = np.zeros((r_f.size, 3), dtype=float)

        axis = r_f == 0.0
        if np.any(axis):
            inherited = self.parent.velocity(
                xb_f[axis], yb_f[axis], zb_f[axis], tb_f[axis]
            )
            base_f[axis] = np.asarray(inherited, dtype=float).reshape((-1, 3))

        nonaxis = ~axis
        if np.any(nonaxis):
            p = self.split_profile_logX(lx_f[nonaxis], eta_f[nonaxis])
            r = r_f[nonaxis]
            xx = xb_f[nonaxis]
            yy = yb_f[nonaxis]
            qq = q_f[nonaxis]
            v0 = np.asarray(p["v0_base"], dtype=float)
            radial_scale = v0 / (2.0 * qq)
            log_r = np.log(r)
            utheta_base = np.exp(
                -(self.A + 0.5) * np.log(qq) + p["log_F_base"] + log_r
            )
            # Compute the correction multiplicatively from the representable
            # relative edit without ever forming 1+edit.
            utheta_delta = utheta_base * np.asarray(p["relative_edit"], dtype=float)
            tx = -yy / r
            ty = xx / r
            base_f[nonaxis, 0] = radial_scale * xx + tx * utheta_base
            base_f[nonaxis, 1] = radial_scale * yy + ty * utheta_base
            base_f[nonaxis, 2] = 0.0
            delta_f[nonaxis, 0] = tx * utheta_delta
            delta_f[nonaxis, 1] = ty * utheta_delta
            delta_f[nonaxis, 2] = 0.0

        if np.any(~np.isfinite(base_f)) or np.any(~np.isfinite(delta_f)):
            raise RuntimeError("relative-swirl split Cartesian velocity became non-finite")
        return base_f.reshape(shape + (3,)), delta_f.reshape(shape + (3,))

    def representation_report(self) -> dict[str, Any]:
        eta = np.asarray([-0.2, 0.0, 0.2], dtype=float)
        s = np.asarray(
            [self.compensator.y1, self.compensator.y2, self.hold_length], dtype=float
        )
        # Broadcast eta across representative source locations.
        ee, ss = np.meshgrid(eta, s, indexing="ij")
        p = self.split_profile_logX(self.log_X_flatten_end + ss, ee)
        edit = np.asarray(p["relative_edit"], dtype=float)
        nonzero = np.abs(edit[edit != 0.0])
        eps = np.finfo(float).eps
        rounded = (1.0 + edit) - 1.0
        sol = self._coefficients(eta)
        return {
            "schema": "kokuno-agent1-current-relative-swirl-split-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "eta_probe": eta.tolist(),
            "c1": np.asarray(sol.c1).tolist(),
            "c2": np.asarray(sol.c2).tolist(),
            "max_abs_relative_edit": float(np.max(np.abs(edit))),
            "min_nonzero_abs_relative_edit": (
                float(np.min(nonzero)) if nonzero.size else 0.0
            ),
            "float64_epsilon": float(eps),
            "max_abs_rounded_1_plus_edit_minus_1": float(np.max(np.abs(rounded))),
            "edit_is_below_float64_relative_epsilon": bool(
                np.max(np.abs(edit)) < eps
            ),
            "max_abs_angular_residual": float(np.max(np.abs(sol.angular_residual))),
            "max_abs_pressure_residual": float(np.max(np.abs(sol.pressure_residual))),
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
            "target_semantic_sha256": _semantic_value(self.target),
            "hold_semantic_sha256": _semantic_value(self.parent),
            "compensator_semantic_sha256": _semantic_value(self.compensator),
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
