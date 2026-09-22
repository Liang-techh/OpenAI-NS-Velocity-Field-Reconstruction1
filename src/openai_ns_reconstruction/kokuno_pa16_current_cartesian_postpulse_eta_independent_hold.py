"""Carry the current A1 Cartesian leading field through Kokuno's eta-independent hold.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
corrected 2026-09-09 reconstruction.

Immediately after the public post-pulse eta-flattening stage, the corrected
reconstruction says to retain the now eta-independent exponential for

    T_hold = 30 log(1/lambda),

before the later pair of relative-swirl bumps.  The public schedule has already
set U=0 permanently after the axial pulse.  This module implements only that
bounded hold.  It does not materialize the relative-swirl bumps, their
zero-integral adjustment, the later l:-lambda -> -1 transition, exterior heat,
pressure, forcing, or a complete Navier--Stokes residual.

The exact current #1140 endpoint is preserved.  On the hold, with
s=log(X/X_flatten_end),

    E(s,eta) = E_flat(eta) exp[-(1/2+lambda)s],
    F(s,eta) = F_flat(eta) exp[-(1+lambda)s],
    U = 0,

and physical M,M_eta remain constant, hence their normalized versions decay as
exp(-s).  The parent endpoint is eta-independent up to its documented
float64 roundoff.  We do not snap it to a new hidden/source-exact constant.
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

SCHEMA = "kokuno-pa16-current-cartesian-postpulse-eta-independent-hold-v1"
PARENT_EXACT_HEAD = "c5442c11165d7f0976889b826bf58dce38a3d3dd"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "eta_flattening_handoff": "after eta flattening, retain the now eta-independent exponential",
    "hold_length": "T_hold=30*log(1/lambda)",
    "hold_E": "E=E_flat*exp[-(1/2+lambda)*s], s=log(X/X_flatten_end)",
    "hold_F": "F=F_flat*exp[-(1+lambda)*s]",
    "hold_ell": "ell=X*d_X log(sqrt(2X)E)=-lambda",
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1140 post-pulse eta-flattening Cartesian identity",
    "hold_length": "evaluate the public exact expression 30*log(1/lambda) in float64",
    "seam": "use actual #1140 E/F/M/M_eta at the flattening endpoint; do not snap tiny eta-roundoff",
    "coordinate": "continue the inherited overflow-safe log-X/log-F chart",
    "current_M_policy": "with public U=0, carry physical M/M_eta constantly so normalized primitives scale exp(-s)",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "source_postpulse_eta_flattening_formula_materialized": True,
    "repository_autonomous_T_f_materialized": True,
    "source_exact_T_f_recovered": False,
    "current_cartesian_postpulse_eta_flattening_composed": True,
    "source_terminal_hold_after_eta_flattening_materialized": True,
    "current_cartesian_postpulse_eta_independent_hold_composed": True,
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
class KokunoPA16CurrentCartesianPostPulseEtaIndependentHold:
    """Current Cartesian leading candidate through the public eta-independent hold."""

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
            raise ValueError("eta-independent hold requires exact #1140 eta-flattening parent")
        if truth["source_terminal_hold_after_eta_flattening_materialized"]:
            raise ValueError("parent unexpectedly already contains the post-flattening hold")
        if not (0.0 < self.lambda_value < 1.0):
            raise ValueError("public hold length requires 0<lambda<1")
        if not math.isfinite(self.hold_length) or self.hold_length <= 0.0:
            raise ValueError("eta-independent hold length must be finite and positive")
        if self.log_radius_q1_hold_end >= math.log(np.finfo(float).max):
            raise ValueError("eta-independent hold endpoint has non-representable physical radius")

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
    def log_X_hold_end(self) -> float:
        return self.log_X_flatten_end + self.hold_length

    @property
    def log_radius_q1_hold_end(self) -> float:
        return 0.5 * (_LOG2 + self.log_X_hold_end)

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
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_hold_end))
        if np.any(log_arr > self.log_X_hold_end + tol):
            raise ValueError("log_X lies beyond the bounded eta-independent hold")
        log_arr = np.minimum(log_arr, self.log_X_hold_end)
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

    def _hold_profile_logX(
        self, log_X: np.ndarray, eta: np.ndarray
    ) -> dict[str, np.ndarray]:
        s = np.asarray(log_X, dtype=float) - self.log_X_flatten_end
        if np.any((s < -2.0e-12) | (s > self.hold_length + 2.0e-12)):
            raise ValueError("eta-independent hold profile requires 0<=s<=T_hold")
        s = np.clip(s, 0.0, self.hold_length)
        state = self._flatten_endpoint_state(eta)

        E = state["E"] * np.exp(-(0.5 + self.lambda_value) * s)
        log_F = state["log_F"] - (1.0 + self.lambda_value) * s
        F = np.exp(log_F)
        U = np.zeros_like(E)
        U_eta = np.zeros_like(E)

        carry = np.exp(-s)
        M_ratio = carry * state["M_over_X"]
        M_eta_ratio = carry * state["M_eta_over_X"]

        axis = self.geometry.physical_profiles.axis_profiles.values(
            np.zeros_like(eta), eta
        )
        L = np.asarray(axis["L"], dtype=float)
        d = np.asarray(axis["d"], dtype=float)
        v0 = (-2.0 * self.D * eta * M_ratio - d * M_eta_ratio) / L

        arrays = (E, F, U, M_ratio, M_eta_ratio, v0)
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("eta-independent hold profile produced non-finite values")
        if np.any(E <= 0.0) or np.any(F < 0.0):
            raise RuntimeError("eta-independent hold lost positive E/nonnegative F")

        return {
            "log_X": np.asarray(log_X, dtype=float),
            "eta": np.asarray(eta, dtype=float),
            "hold_s": s,
            "E": E,
            "E_eta": np.zeros_like(E),
            "F": F,
            "log_F": log_F,
            "U": U,
            "U_eta": U_eta,
            "M_over_X": M_ratio,
            "M_eta_over_X": M_eta_ratio,
            "v0": v0,
            "E_DlogX": -(0.5 + self.lambda_value) * E,
            "F_DlogX": -(1.0 + self.lambda_value) * F,
            "log_F_DlogX": np.full_like(E, -(1.0 + self.lambda_value)),
            "U_DlogX": np.zeros_like(U),
            "ell": np.full_like(E, -self.lambda_value),
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
            p = self._hold_profile_logX(lf[active], ef[active])
            F[active] = p["F"]
            log_F[active] = p["log_F"]
            U[active] = p["U"]
            E[active] = p["E"]
            M_ratio[active] = p["M_over_X"]
            M_eta_ratio[active] = p["M_eta_over_X"]
            v0[active] = p["v0"]
            hold_s[active] = p["hold_s"]
            region[active] = "current_cartesian_postpulse_eta_independent_hold_logX"

        arrays = (F, log_F, U, E, M_ratio, M_eta_ratio, v0)
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("post-pulse hold Cartesian candidate produced non-finite values")
        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "postpulse_hold_s": hold_s.reshape(shape),
            "F_current_leading_postpulse_hold": F.reshape(shape),
            "log_F_current_leading_postpulse_hold": log_F.reshape(shape),
            "U_current_leading_postpulse_hold": U.reshape(shape),
            "E_current_leading_postpulse_hold": E.reshape(shape),
            "M_over_X_current_leading_postpulse_hold": M_ratio.reshape(shape),
            "M_eta_over_X_current_leading_postpulse_hold": M_eta_ratio.reshape(shape),
            "v0_current_leading_postpulse_hold": v0.reshape(shape),
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
            p = self._hold_profile_logX(lf[active], ef[active])
            F_D[active] = p["F_DlogX"]
            log_F_D[active] = p["log_F_DlogX"]
            U_D[active] = p["U_DlogX"]
            E_D[active] = p["E_DlogX"]
            ell[active] = p["ell"]

        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_DlogX_current_leading_postpulse_hold": F_D.reshape(shape),
            "log_F_DlogX_current_leading_postpulse_hold": log_F_D.reshape(shape),
            "U_DlogX_current_leading_postpulse_hold": U_D.reshape(shape),
            "E_DlogX_current_leading_postpulse_hold": E_D.reshape(shape),
            "ell_current_leading_postpulse_hold": ell.reshape(shape),
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
            log_F = np.asarray(p["log_F_current_leading_postpulse_hold"], dtype=float)
            U = np.asarray(p["U_current_leading_postpulse_hold"], dtype=float)
            v0 = np.asarray(p["v0_current_leading_postpulse_hold"], dtype=float)
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
            raise RuntimeError("post-pulse eta-independent hold Cartesian velocity became non-finite")
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
                "postpulse_stage": "eta_independent_exponential_hold_only",
                "T_f": self.T_f,
                "T_f_role": "inherited_repository_autonomous_not_source_exact",
                "hold_length": self.hold_length,
                "hold_length_role": "public_formula_30_log_inverse_lambda",
                "relative_swirl_materialized": False,
                "exterior_heat_materialized": False,
                "mutable": False,
            },
            "truth_boundary": self.truth_boundary,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianPostPulseEtaIndependentHold":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected post-pulse eta-independent hold schema")
        parent_payload = payload.get("parent_current_cartesian_postpulse_eta_flattening")
        if not isinstance(parent_payload, Mapping):
            raise ValueError("missing current post-pulse eta-flattening parent configuration")
        candidate = cls(
            parent=KokunoPA16CurrentCartesianPostPulseEtaFlattening.from_configuration(
                parent_payload
            )
        )
        if _canonical_json(dict(payload)) != _canonical_json(candidate.configuration()):
            raise ValueError("serialized post-pulse hold source/truth binding changed")
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
    ) -> "KokunoPA16CurrentCartesianPostPulseEtaIndependentHold":
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
                "log_X_flatten_end": self.log_X_flatten_end,
                "log_X_hold_end": self.log_X_hold_end,
            },
            "truth_boundary": self.truth_boundary,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    def hold_report(
        self,
        eta: Any = (-1.0, -0.75, -0.25, 0.0, 0.4, 0.8, 1.0),
    ) -> dict[str, Any]:
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")

        start = self._hold_profile_logX(
            np.full(eta_arr.shape, self.log_X_flatten_end), eta_arr
        )
        end = self._hold_profile_logX(
            np.full(eta_arr.shape, self.log_X_hold_end), eta_arr
        )
        start_E = np.asarray(start["E"], dtype=float)
        end_E = np.asarray(end["E"], dtype=float)
        start_spread = (
            float((np.max(start_E) - np.min(start_E)) / np.max(start_E))
            if np.max(start_E) > 0.0
            else 0.0
        )
        end_spread = (
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
            "T_f_inherited_autonomous": self.T_f,
            "hold_length_public_formula": self.hold_length,
            "expected_hold_length": 30.0 * math.log(1.0 / self.lambda_value),
            "log_X_flatten_end": self.log_X_flatten_end,
            "log_X_hold_end": self.log_X_hold_end,
            "start_E_relative_eta_spread": start_spread,
            "end_E_relative_eta_spread": end_spread,
            "ell_on_hold": -self.lambda_value,
            "max_abs_U_hold_endpoint": float(np.max(np.abs(end["U"]))),
            "max_abs_M_over_X_hold_endpoint": float(
                np.max(np.abs(end["M_over_X"]))
            ),
            "max_abs_M_eta_over_X_hold_endpoint": float(
                np.max(np.abs(end["M_eta_over_X"]))
            ),
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "this implements only the public eta-independent exponential hold; the relative-swirl bumps are not materialized",
                "T_f=100 is inherited from exact #1140 and remains repository-autonomous rather than source-exact",
                "tiny eta spread inherited at the flattening seam is preserved rather than snapped to hidden/source-exact data",
                "the receipt is source-schedule/executable evidence, not independent PDE validation",
                "relative-swirl, exterior heat/global leading, pressure, forcing, and held-out complete NS residual remain absent",
            ],
        }


__all__ = ["KokunoPA16CurrentCartesianPostPulseEtaIndependentHold"]
