"""Executable public ideal-profile target for the Kokuno final radial join.

Pinned public provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.

The public reconstruction introduces

    f(eta) = (1+eta^2)^-1,
    X_R = 110 (C P_*)^10,
    x = X/X_R,
    U_id = 4 eta,
    E_id = P_* f x^(1/10),
    X_h = X_R exp(-5),

and states that the repaired inner/exterior profiles agree with this same
ideal pair in a neighborhood of X_h.  It also displays five normalized moment
targets used by the nonlinear five-bump repair.

This module materializes only that public ideal target and the displayed
moments.  The corrected source does not identify a numerical P_* in the
material used by Agent 1, so ``P_star=1`` is an explicit repository-autonomous
positive realization.  It is not inferred from figures, residual fitting, or
hidden OpenAI data.  Likewise, this module does not fabricate the two-U/three-E
bump coefficients or the unspecified width of the source join neighborhood.
It therefore is not an outer/global velocity and is not PDE validation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_source_c_normalized_physical_center import (
    SELECTED_SOURCE_C,
    KokunoPA10SourceCNormalizedPhysicalCenter,
)


SCHEMA = "kokuno-pa10-ideal-join-profile-contract-v1"
DEFAULT_P_STAR = 1.0

_SOURCE_FORMULAS = {
    "f": "f=(1+eta^2)^-1",
    "X_R": "X_R=110(C P_*)^10",
    "x": "x=X/X_R",
    "ideal_U": "U_id=4 eta",
    "ideal_E": "E_id=P_* f x^(1/10)",
    "physical_F": "E=sqrt(2X)F, so F_id=E_id/sqrt(2X)",
    "join_radius": "X_h=X_R exp(-5)",
    "hat_M": "hat M=4 eta x",
    "hat_I": "hat I=(5 sqrt(2)/8) P_* f x^(8/5)",
    "hat_J": "hat J=4 eta hat I",
    "hat_S": "hat S=16 eta^2 x-(5/12)P_*^2 f^2 x^(6/5)",
    "C_p": "C_p=(5/2)P_*^2 f^2 x^(1/5)",
    "repair": "two U-bumps plus three E-bumps solve the five moment equations on -6<log x<-5",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_ideal_join_target_executable": True,
    "public_five_normalized_moment_targets_executable": True,
    "ideal_join_first_derivatives_executable": True,
    "source_C_inherited_from_agent1_normalization": True,
    "P_star_is_repository_autonomous_positive_realization": True,
    "source_P_star_numerically_identified": False,
    "P_star_inferred_from_figure": False,
    "five_bump_nonlinear_repair_materialized": False,
    "source_join_neighborhood_width_known": False,
    "source_join_neighborhood_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "source_center_is_final_corrected_fixed_point": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
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
class KokunoPA10IdealJoinProfileContract:
    """Public ideal pair and five moment targets at the source final-join scale."""

    P_star: float = DEFAULT_P_STAR
    source_normalized: KokunoPA10SourceCNormalizedPhysicalCenter = field(
        default_factory=KokunoPA10SourceCNormalizedPhysicalCenter,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        p = float(self.P_star)
        if not math.isfinite(p) or p <= 0.0:
            raise ValueError("P_star must be finite and strictly positive")
        if not isinstance(self.source_normalized, KokunoPA10SourceCNormalizedPhysicalCenter):
            raise TypeError("source_normalized must be KokunoPA10SourceCNormalizedPhysicalCenter")
        if float(self.source_normalized.C) != SELECTED_SOURCE_C:
            raise ValueError(f"source-normalized C must remain {SELECTED_SOURCE_C:g}")
        object.__setattr__(self, "P_star", p)

    @property
    def C(self) -> float:
        return float(self.source_normalized.C)

    @property
    def X_R(self) -> float:
        return 110.0 * (self.C * self.P_star) ** 10

    @property
    def X_h(self) -> float:
        return self.X_R * math.exp(-5.0)

    @property
    def inner_X_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.source_normalized.physical_profiles.source_X_interval)

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any(X_arr <= 0.0):
            raise ValueError("X must be strictly positive for the ideal outer profile")
        return X_arr, eta_arr

    @staticmethod
    def f(eta: Any) -> np.ndarray:
        eta_arr = _finite(eta, "eta")
        return 1.0 / (1.0 + eta_arr * eta_arr)

    def values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast(X, eta)
        f = self.f(eta_arr)
        x = X_arr / self.X_R
        U = 4.0 * eta_arr
        E = self.P_star * f * np.power(x, 0.1)
        F = E / np.sqrt(2.0 * X_arr)
        return {
            "X": X_arr,
            "eta": eta_arr,
            "x": x,
            "f": f,
            "U_ideal": U,
            "E_ideal": E,
            "F_ideal": F,
        }

    def derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        values = self.values(X, eta)
        X_arr = values["X"]
        eta_arr = values["eta"]
        E = values["E_ideal"]
        F = values["F_ideal"]
        eta_factor = -2.0 * eta_arr / (1.0 + eta_arr * eta_arr)
        return {
            "U_ideal_X": np.zeros_like(X_arr),
            "U_ideal_eta": np.full_like(X_arr, 4.0),
            "E_ideal_X": 0.1 * E / X_arr,
            "E_ideal_eta": eta_factor * E,
            "F_ideal_X": -0.4 * F / X_arr,
            "F_ideal_eta": eta_factor * F,
        }

    def normalized_moments(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast(X, eta)
        x = X_arr / self.X_R
        f = self.f(eta_arr)
        hat_I = (5.0 * math.sqrt(2.0) / 8.0) * self.P_star * f * np.power(x, 8.0 / 5.0)
        return {
            "hat_M": 4.0 * eta_arr * x,
            "hat_I": hat_I,
            "hat_J": 4.0 * eta_arr * hat_I,
            "hat_S": 16.0 * eta_arr * eta_arr * x
            - (5.0 / 12.0) * self.P_star**2 * f * f * np.power(x, 6.0 / 5.0),
            "C_p": (5.0 / 2.0) * self.P_star**2 * f * f * np.power(x, 1.0 / 5.0),
        }

    def join_gap_diagnostic(self) -> dict[str, float | bool]:
        inner_max = float(self.inner_X_interval[1])
        if not math.isfinite(inner_max) or inner_max <= 0.0:
            raise RuntimeError("current PA.10 source-X upper endpoint is not positive finite")
        ratio = self.X_h / inner_max
        return {
            "inner_X_max": inner_max,
            "X_R": self.X_R,
            "X_h": self.X_h,
            "X_h_over_inner_X_max": ratio,
            "log_X_h_over_inner_X_max": math.log(ratio),
            "source_join_is_outside_current_inner_domain": self.X_h > inner_max,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "P_star": self.P_star,
            "source_C": self.C,
            "P_star_semantics": "repository-autonomous-positive-realization-not-source-identified",
        }

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoPA10IdealJoinProfileContract":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        if float(payload.get("source_C")) != SELECTED_SOURCE_C:
            raise ValueError(f"serialized source_C must remain {SELECTED_SOURCE_C:g}")
        if payload.get("P_star_semantics") != "repository-autonomous-positive-realization-not-source-identified":
            raise ValueError("P_star semantics mismatch")
        return cls(P_star=float(payload["P_star"]))

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA10IdealJoinProfileContract":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_C": self.C,
            "P_star": self.P_star,
            "source_formulas": _SOURCE_FORMULAS,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        join_values = self.values(self.X_h, 0.5)
        moments = self.normalized_moments(self.X_h, 0.5)
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
            "join_gap": self.join_gap_diagnostic(),
            "probe_at_X_h_eta_0p5": {
                "U_ideal": float(join_values["U_ideal"]),
                "E_ideal": float(join_values["E_ideal"]),
                "F_ideal": float(join_values["F_ideal"]),
                **{name: float(value) for name, value in moments.items()},
            },
            "truth_boundary": self.truth_boundary,
        }
