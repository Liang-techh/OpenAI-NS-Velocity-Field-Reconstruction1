"""Materialize Kokuno's unedited post-pulse hold and a safe Cartesian prefix.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
corrected 2026-09-09 reconstruction.

After the public eta-flattening stage, the reconstruction says to retain the
now eta-independent exponential for

    T_hold = 30 log(1/lambda).

It then places two relative-swirl bumps of width 0.3 at three and one log-X
units before the end of that hold.  Therefore the *unedited* exponential is a
valid public baseline over the full hold, but it is not the complete source
field inside the later bump neighborhoods.

This module makes that distinction executable.  It exposes the full unedited
baseline as a low-dimensional helper, while the unified current Cartesian
velocity is extended only to a conservative pre-bump prefix ending 3.5 units
before the public hold endpoint.  That guard is a repository-autonomous domain
safety margin, not a source parameter and not a residual-tuning knob.  The
relative-swirl bumps themselves remain unmaterialized.

The post-pulse source schedule has U=0 permanently.  With
s=log(X/X_flatten_end), the unedited baseline is

    E = E_flat exp[-(1/2+lambda)s],
    F = F_flat exp[-(1+lambda)s],
    U = 0,

and physical M,M_eta stay constant, so M/X and M_eta/X decay as exp(-s).
The exact current #1140 endpoint is preserved; its tiny binary64 eta-spread is
not snapped to invented hidden data.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_postpulse_eta_flattening import (
    KokunoPA16CurrentCartesianPostPulseEtaFlattening,
)

SCHEMA = "kokuno-pa16-current-cartesian-postpulse-eta-independent-hold-prefix-v1"
PARENT_EXACT_HEAD = "c5442c11165d7f0976889b826bf58dce38a3d3dd"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
RELATIVE_SWIRL_WIDTH_PUBLIC = 0.3
FIRST_RELATIVE_SWIRL_CENTER_FROM_END = 3.0
SECOND_RELATIVE_SWIRL_CENTER_FROM_END = 1.0
PRE_BUMP_GUARD_FROM_HOLD_END_AUTONOMOUS = 3.5
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "hold_length": "T_hold=30*log(1/lambda)",
    "hold_unedited_E": "E_unedited=E_flat*exp[-(1/2+lambda)*s]",
    "hold_unedited_F": "F_unedited=F_flat*exp[-(1+lambda)*s]",
    "hold_unedited_ell": "ell_unedited=-lambda",
    "relative_swirl_geometry": (
        "two relative-swirl bumps of width 0.3 centered three and one "
        "log-X units before the hold endpoint"
    ),
    "relative_swirl_constraint": (
        "the two bumps set I=XH/(1-lambda) there and have total "
        "integral int(E^2-E_unedited^2)dy=0"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1140 post-pulse eta-flattening Cartesian identity",
    "hold_length": "evaluate public 30*log(1/lambda) in float64",
    "full_hold": "materialize only the public unedited exponential baseline",
    "candidate_domain": (
        "extend unified Cartesian candidate only to T_hold-3.5, a conservative "
        "repository-autonomous guard before the first relative-swirl bump"
    ),
    "seam": "use actual #1140 E/F/M/M_eta at the flattening endpoint",
    "coordinate": "continue inherited overflow-safe log-X/log-F chart",
    "current_M_policy": (
        "with public U=0, physical M/M_eta remain constant so normalized "
        "primitives scale exp(-s)"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "source_postpulse_eta_flattening_formula_materialized": True,
    "repository_autonomous_T_f_materialized": True,
    "source_exact_T_f_recovered": False,
    "current_cartesian_postpulse_eta_flattening_composed": True,
    "source_unedited_eta_independent_hold_baseline_materialized": True,
    "current_cartesian_pre_relative_swirl_hold_prefix_composed": True,
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


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix:
    """Current Cartesian leading candidate up to the pre-relative-swirl hold prefix."""

    parent: KokunoPA16CurrentCartesianPostPulseEtaFlattening = field(
        default_factory=KokunoPA16CurrentCartesianPostPulseEtaFlattening,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianPostPulseEtaFlattening):
            raise TypeError("parent must be KokunoPA16CurrentCartesianPostPulseEtaFlattening")
        truth = self.parent.truth_boundary
        if not truth["current_cartesian_postpulse_eta_flattening_composed"]:
            raise ValueError("hold prefix requires exact #1140 eta-flattening parent")
        if truth["source_terminal_hold_after_eta_flattening_materialized"]:
            raise ValueError("parent unexpectedly already contains the terminal hold")
        if truth["source_relative_swirl_bumps_materialized"]:
            raise ValueError("parent unexpectedly already contains relative-swirl bumps")
        if not (0.0 < self.lambda_value < 1.0):
            raise ValueError("public hold length requires 0<lambda<1")
        if self.prefix_length <= 0.0:
            raise ValueError("pre-relative-swirl prefix must have positive length")
        if PRE_BUMP_GUARD_FROM_HOLD_END_AUTONOMOUS <= (
            FIRST_RELATIVE_SWIRL_CENTER_FROM_END + RELATIVE_SWIRL_WIDTH_PUBLIC
        ):
            raise ValueError("autonomous pre-bump guard must stay strictly before bump width")
        if self.log_radius_q1_prefix_end >= math.log(np.finfo(float).max):
            raise ValueError("pre-relative-swirl prefix endpoint has non-representable physical radius")

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
        return float(self.parent.T_f)

    @property
    def log_X_flatten_end(self) -> float:
        return float(self.parent.log_X_flatten_end)

    @property
    def hold_length(self) -> float:
        return 30.0 * math.log(1.0 / self.lambda_value)

    @property
    def first_relative_swirl_center_s(self) -> float:
        return self.hold_length - FIRST_RELATIVE_SWIRL_CENTER_FROM_END

    @property
    def second_relative_swirl_center_s(self) -> float:
        return self.hold_length - SECOND_RELATIVE_SWIRL_CENTER_FROM_END

    @property
    def prefix_length(self) -> float:
        return self.hold_length - PRE_BUMP_GUARD_FROM_HOLD_END_AUTONOMOUS

    @property
    def log_X_hold_end(self) -> float:
        return self.log_X_flatten_end + self.hold_length

    @property
    def log_X_prefix_end(self) -> float:
        return self.log_X_flatten_end + self.prefix_length

    @property
    def log_radius_q1_prefix_end(self) -> float:
        return 0.5 * (_LOG2 + self.log_X_prefix_end)

    def similarity_coordinates_logX(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.parent.similarity_coordinates_logX(x, y, z, t)

    def _broadcast_candidate_log_similarity(
        self, log_X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_prefix_end))
        if np.any(log_arr > self.log_X_prefix_end + tol):
            raise ValueError(
                "log_X lies beyond the bounded pre-relative-swirl Cartesian prefix"
            )
        log_arr = np.minimum(log_arr, self.log_X_prefix_end)
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return log_arr, eta_arr

    def _flatten_endpoint_state(self, eta: np.ndarray) -> dict[str, np.ndarray]:
        log_end = np.full(np.asarray(eta).shape, self.log_X_flatten_end, dtype=float)
        p = self.parent.similarity_profile_values_logX(log_end, eta)
        E = np.asarray(
            p["E_current_leading_postpulse_eta_flattening"], dtype=float
        )
        log_F = np.asarray(
            p["log_F_current_leading_postpulse_eta_flattening"], dtype=float
        )
        M_ratio = np.asarray(
            p["M_over_X_current_leading_postpulse_eta_flattening"], dtype=float
        )
        M_eta_ratio = np.asarray(
            p["M_eta_over_X_current_leading_postpulse_eta_flattening"], dtype=float
        )
        if np.any(E <= 0.0) or np.any(~np.isfinite(E)) or np.any(~np.isfinite(log_F)):
            raise RuntimeError("parent flattening endpoint must have finite positive E/log_F")
        return {
            "E": E,
            "log_F": log_F,
            "M_over_X": M_ratio,
            "M_eta_over_X": M_eta_ratio,
        }

    def unedited_hold_profile_logX(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        """Evaluate the public unedited exponential baseline over the full hold."""

        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        s = log_arr - self.log_X_flatten_end
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_hold_end))
        if np.any((s < -tol) | (s > self.hold_length + tol)):
            raise ValueError("unedited hold baseline requires 0<=s<=T_hold")
        s = np.clip(s, 0.0, self.hold_length)
        state = self._flatten_endpoint_state(eta_arr)

        E = state["E"] * np.exp(-(0.5 + self.lambda_value) * s)
        log_F = state["log_F"] - (1.0 + self.lambda_value) * s
        F = np.exp(log_F)
        U = np.zeros_like(E)
        carry = np.exp(-s)
        M_ratio = carry * state["M_over_X"]
        M_eta_ratio = carry * state["M_eta_over_X"]

        axis = self.geometry.physical_profiles.axis_profiles.values(
            np.zeros_like(eta_arr), eta_arr
        )
        L = np.asarray(axis["L"], dtype=float)
        d = np.asarray(axis["d"], dtype=float)
        v0 = (-2.0 * self.D * eta_arr * M_ratio - d * M_eta_ratio) / L

        arrays = (E, F, U, M_ratio, M_eta_ratio, v0)
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("unedited eta-independent hold baseline became non-finite")
        if np.any(E <= 0.0) or np.any(F < 0.0):
            raise RuntimeError("unedited eta-independent hold lost positive E/nonnegative F")

        return {
            "log_X": log_arr,
            "eta": eta_arr,
            "hold_s": s,
            "E_unedited": E,
            "E_eta_unedited": np.zeros_like(E),
            "F_unedited": F,
            "log_F_unedited": log_F,
            "U_unedited": U,
            "U_eta_unedited": np.zeros_like(U),
            "M_over_X_unedited": M_ratio,
            "M_eta_over_X_unedited": M_eta_ratio,
            "v0_unedited": v0,
            "E_DlogX_unedited": -(0.5 + self.lambda_value) * E,
            "F_DlogX_unedited": -(1.0 + self.lambda_value) * F,
            "log_F_DlogX_unedited": np.full_like(E, -(1.0 + self.lambda_value)),
            "U_DlogX_unedited": np.zeros_like(U),
            "ell_unedited": np.full_like(E, -self.lambda_value),
        }

    def similarity_profile_values_logX(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_candidate_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)

        F = np.empty_like(lf)
        log_F = np.empty_like(lf)
        U = np.empty_like(lf)
        E = np.empty_like(lf)
        M_ratio = np.empty_like(lf)
        M_eta_ratio = np.empty_like(lf)
        v0 = np.empty_like(lf)
        hold_s = np.zeros_like(lf)
        region = np.empty(lf.shape, dtype=object)

        inherited = lf <= self.log_X_flatten_end
        if np.any(inherited):
            p = self.parent.similarity_profile_values_logX(lf[inherited], ef[inherited])
            F[inherited] = np.asarray(
                p["F_current_leading_postpulse_eta_flattening"], dtype=float
            )
            log_F[inherited] = np.asarray(
                p["log_F_current_leading_postpulse_eta_flattening"], dtype=float
            )
            U[inherited] = np.asarray(
                p["U_current_leading_postpulse_eta_flattening"], dtype=float
            )
            E[inherited] = np.asarray(
                p["E_current_leading_postpulse_eta_flattening"], dtype=float
            )
            M_ratio[inherited] = np.asarray(
                p["M_over_X_current_leading_postpulse_eta_flattening"], dtype=float
            )
            M_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_leading_postpulse_eta_flattening"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_leading_postpulse_eta_flattening"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        active = ~inherited
        if np.any(active):
            p = self.unedited_hold_profile_logX(lf[active], ef[active])
            F[active] = p["F_unedited"]
            log_F[active] = p["log_F_unedited"]
            U[active] = p["U_unedited"]
            E[active] = p["E_unedited"]
            M_ratio[active] = p["M_over_X_unedited"]
            M_eta_ratio[active] = p["M_eta_over_X_unedited"]
            v0[active] = p["v0_unedited"]
            hold_s[active] = p["hold_s"]
            region[active] = "current_cartesian_pre_relative_swirl_hold_prefix_logX"

        arrays = (F, log_F, U, E, M_ratio, M_eta_ratio, v0)
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("pre-relative-swirl Cartesian prefix produced non-finite values")
        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "postpulse_hold_s": hold_s.reshape(shape),
            "F_current_leading_pre_relative_swirl": F.reshape(shape),
            "log_F_current_leading_pre_relative_swirl": log_F.reshape(shape),
            "U_current_leading_pre_relative_swirl": U.reshape(shape),
            "E_current_leading_pre_relative_swirl": E.reshape(shape),
            "M_over_X_current_leading_pre_relative_swirl": M_ratio.reshape(shape),
            "M_eta_over_X_current_leading_pre_relative_swirl": M_eta_ratio.reshape(shape),
            "v0_current_leading_pre_relative_swirl": v0.reshape(shape),
        }

    def similarity_log_radial_derivatives(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_candidate_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)
        F_D = np.empty_like(lf)
        log_F_D = np.empty_like(lf)
        U_D = np.empty_like(lf)
        E_D = np.empty_like(lf)
        ell = np.empty_like(lf)

        inherited = lf <= self.log_X_flatten_end
        if np.any(inherited):
            d = self.parent.similarity_log_radial_derivatives(
                lf[inherited], ef[inherited]
            )
            F_D[inherited] = np.asarray(
                d["F_DlogX_current_leading_postpulse_eta_flattening"], dtype=float
            )
            log_F_D[inherited] = np.asarray(
                d["log_F_DlogX_current_leading_postpulse_eta_flattening"], dtype=float
            )
            U_D[inherited] = np.asarray(
                d["U_DlogX_current_leading_postpulse_eta_flattening"], dtype=float
            )
            E_D[inherited] = np.asarray(
                d["E_DlogX_current_leading_postpulse_eta_flattening"], dtype=float
            )
            ell[inherited] = np.asarray(
                d["ell_current_leading_postpulse_eta_flattening"], dtype=float
            )

        active = ~inherited
        if np.any(active):
            p = self.unedited_hold_profile_logX(lf[active], ef[active])
            F_D[active] = p["F_DlogX_unedited"]
            log_F_D[active] = p["log_F_DlogX_unedited"]
            U_D[active] = p["U_DlogX_unedited"]
            E_D[active] = p["E_DlogX_unedited"]
            ell[active] = p["ell_unedited"]

        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_DlogX_current_leading_pre_relative_swirl": F_D.reshape(shape),
            "log_F_DlogX_current_leading_pre_relative_swirl": log_F_D.reshape(shape),
            "U_DlogX_current_leading_pre_relative_swirl": U_D.reshape(shape),
            "E_DlogX_current_leading_pre_relative_swirl": E_D.reshape(shape),
            "ell_current_leading_pre_relative_swirl": ell.reshape(shape),
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
                p["log_F_current_leading_pre_relative_swirl"], dtype=float
            )
            U = np.asarray(p["U_current_leading_pre_relative_swirl"], dtype=float)
            v0 = np.asarray(p["v0_current_leading_pre_relative_swirl"], dtype=float)
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
            raise RuntimeError("pre-relative-swirl Cartesian velocity became non-finite")
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
            "parent_current_cartesian_postpulse_eta_flattening": self.parent.configuration(),
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
                "postpulse_stage": "unedited_eta_independent_hold_baseline_plus_pre_bump_prefix",
                "T_f": self.T_f,
                "T_f_role": "inherited_repository_autonomous_not_source_exact",
                "hold_length": self.hold_length,
                "hold_length_role": "public_formula_30_log_inverse_lambda",
                "relative_swirl_width_public": RELATIVE_SWIRL_WIDTH_PUBLIC,
                "first_relative_swirl_center_s": self.first_relative_swirl_center_s,
                "second_relative_swirl_center_s": self.second_relative_swirl_center_s,
                "pre_bump_guard_from_hold_end": PRE_BUMP_GUARD_FROM_HOLD_END_AUTONOMOUS,
                "pre_bump_guard_role": "repository_autonomous_domain_safety_not_source_parameter",
                "candidate_prefix_length": self.prefix_length,
                "relative_swirl_materialized": False,
                "exterior_heat_materialized": False,
                "mutable": False,
            },
            "truth_boundary": self.truth_boundary,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected post-pulse hold-prefix schema")
        parent_payload = payload.get("parent_current_cartesian_postpulse_eta_flattening")
        if not isinstance(parent_payload, Mapping):
            raise ValueError("missing current post-pulse eta-flattening parent configuration")
        candidate = cls(
            parent=KokunoPA16CurrentCartesianPostPulseEtaFlattening.from_configuration(
                parent_payload
            )
        )
        if _canonical_json(dict(payload)) != _canonical_json(candidate.configuration()):
            raise ValueError("serialized hold-prefix source/truth binding changed")
        return candidate

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix":
        return cls.from_configuration(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

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
                "hold_length": self.hold_length,
                "prefix_length": self.prefix_length,
                "log_X_flatten_end": self.log_X_flatten_end,
                "log_X_prefix_end": self.log_X_prefix_end,
            },
            "truth_boundary": self.truth_boundary,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    def hold_prefix_report(
        self,
        eta: Any = (-1.0, -0.75, -0.25, 0.0, 0.4, 0.8, 1.0),
    ) -> dict[str, Any]:
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")

        start = self.unedited_hold_profile_logX(
            np.full(eta_arr.shape, self.log_X_flatten_end), eta_arr
        )
        prefix_end = self.unedited_hold_profile_logX(
            np.full(eta_arr.shape, self.log_X_prefix_end), eta_arr
        )
        unedited_full_end = self.unedited_hold_profile_logX(
            np.full(eta_arr.shape, self.log_X_hold_end), eta_arr
        )

        def relative_spread(values: np.ndarray) -> float:
            vmax = float(np.max(values))
            return (
                float((np.max(values) - np.min(values)) / vmax)
                if vmax > 0.0
                else 0.0
            )

        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "semantic_sha256": self.semantic_sha256,
            "eta": eta_arr.tolist(),
            "lambda_value": self.lambda_value,
            "T_f_inherited_autonomous": self.T_f,
            "hold_length_public_formula": self.hold_length,
            "relative_swirl_width_public": RELATIVE_SWIRL_WIDTH_PUBLIC,
            "first_relative_swirl_center_s": self.first_relative_swirl_center_s,
            "second_relative_swirl_center_s": self.second_relative_swirl_center_s,
            "pre_bump_guard_from_hold_end_autonomous": (
                PRE_BUMP_GUARD_FROM_HOLD_END_AUTONOMOUS
            ),
            "candidate_prefix_length": self.prefix_length,
            "log_X_flatten_end": self.log_X_flatten_end,
            "log_X_prefix_end": self.log_X_prefix_end,
            "log_X_unedited_hold_end": self.log_X_hold_end,
            "start_E_relative_eta_spread": relative_spread(
                np.asarray(start["E_unedited"], dtype=float)
            ),
            "prefix_end_E_relative_eta_spread": relative_spread(
                np.asarray(prefix_end["E_unedited"], dtype=float)
            ),
            "unedited_full_end_E_relative_eta_spread": relative_spread(
                np.asarray(unedited_full_end["E_unedited"], dtype=float)
            ),
            "ell_unedited": -self.lambda_value,
            "max_abs_U_prefix_endpoint": float(
                np.max(np.abs(prefix_end["U_unedited"]))
            ),
            "max_abs_M_over_X_prefix_endpoint": float(
                np.max(np.abs(prefix_end["M_over_X_unedited"]))
            ),
            "max_abs_M_eta_over_X_prefix_endpoint": float(
                np.max(np.abs(prefix_end["M_eta_over_X_unedited"]))
            ),
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "the full 30log(1/lambda) exponential is an unedited baseline, not the complete source field inside the later relative-swirl bump neighborhoods",
                "unified Cartesian velocity stops at a conservative autonomous guard 3.5 units before the hold endpoint",
                "relative-swirl bumps, their I condition and zero-integral compensation are not materialized",
                "T_f=100 is inherited from exact #1140 and remains repository-autonomous rather than source-exact",
                "this receipt is source-schedule/executable evidence, not independent PDE validation",
                "exterior/global leading, pressure, forcing, and held-out complete NS residual remain absent",
            ],
        }


__all__ = ["KokunoPA16CurrentCartesianPostPulseEtaIndependentHoldPrefix"]
