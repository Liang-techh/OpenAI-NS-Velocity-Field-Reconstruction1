"""Source-compatible Kokuno oscillatory phase-label selection.

The corrected 2026-09-09 reconstruction selects the angular/axial phase label
from source data by

    B_s = sqrt(lambda0 / (epsilon*k^2*(1+u_*^2)^(3/2))),
    x0 = sigma*B_s*u_*/2,
    q = B_s*(K - sigma*u_* g0/(L_s*|g0|^2)),
    tilde_p = R0*q_r,
    j = a specified nearest nonzero integer to k*tilde_p,
    p = j/k,
    p_z = q_z,
    k = ceil(epsilon^(-1/2)).

This module executes exactly those algebraic formulas once ``R0``, ``g0`` and
``K`` are supplied.  The source says "specified nearest nonzero integer" but
does not fix the half-tie convention in the retained statement.  We therefore
use a deterministic *autonomous numerical convention*: half ties away from
zero, and when zero would be nearest choose the sign of the target (with +1 at
an exactly-zero target).  That convention is provenance-labelled and is not
claimed to recover any hidden Kokuno/OpenAI label.

The actual positive-order background, auxiliary-torus data, support partition,
and direct physical velocity remain outside this increment.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from math import ceil, floor, sqrt
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-phase-label-selection-v1"

_SOURCE_FORMULAS = {
    "k": "ceil(epsilon^(-1/2))",
    "B_s": "sqrt(lambda0/(epsilon*k^2*(1+u_*^2)^(3/2)))",
    "x0": "sigma*B_s*u_*/2",
    "q": "B_s*(K-sigma*u_*g0/(L_s*|g0|^2))",
    "tilde_p": "R0*q_r",
    "j": "specified nearest nonzero integer to k*tilde_p",
    "p": "j/k",
    "p_z": "q_z",
}

_TRUTH_BOUNDARY = {
    "source_phase_label_selection_formulas_executable": True,
    "source_requires_p_equals_j_over_k": True,
    "source_requires_k_p_integer_nonzero": True,
    "nearest_nonzero_integer_tie_rule_source_specified": False,
    "autonomous_deterministic_tie_rule_used": True,
    "source_actual_discrete_label_recovered": False,
    "source_actual_background_path_instantiated": False,
    "positive_order_background_corrections_included": False,
    "source_actual_auxiliary_torus_evaluation_instantiated": False,
    "source_actual_support_partition_instantiated": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _bounded_finite(value: Any, name: str, lower: float, upper: float) -> float:
    out = float(value)
    if not np.isfinite(out) or not (lower <= out <= upper):
        raise ValueError(f"{name} must be finite in [{lower},{upper}]")
    return out


def _finite_pair(value: Any, name: str) -> tuple[float, float]:
    arr = np.asarray(value, dtype=float)
    if arr.shape != (2,) or not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be a finite 2-vector (r,z)")
    if np.max(np.abs(arr)) > 1.0e6:
        raise ValueError(f"{name} entries must satisfy |component|<=1e6")
    return float(arr[0]), float(arr[1])


@dataclass(frozen=True)
class KokunoSourcePhaseLabelSelection:
    """Evaluate one source-compatible phase label from supplied source data.

    Defaults are bounded autonomous replay values only.  They are not a claim
    about the source's hidden/discrete choice.
    """

    epsilon: float = 0.25
    lambda0: float = 1.0
    u_star: float = 2.0
    L_s: float = 1.0
    sigma: int = 1
    R0: float = 1.0
    g0: tuple[float, float] = (1.0, 0.5)
    K: tuple[float, float] = (0.5, 1.0)

    def __post_init__(self) -> None:
        epsilon = _bounded_finite(self.epsilon, "epsilon", 1.0e-6, 1.0)
        lambda0 = _bounded_finite(self.lambda0, "lambda0", 1.0e-8, 100.0)
        u_star = _bounded_finite(self.u_star, "u_star", 0.05, 20.0)
        L_s = _bounded_finite(self.L_s, "L_s", 1.0e-8, 1.0e4)
        R0 = _bounded_finite(self.R0, "R0", 1.0e-8, 1.0e4)
        sigma = int(self.sigma)
        if sigma != self.sigma or sigma not in (-1, 1):
            raise ValueError("sigma must be exactly -1 or +1")
        g0 = _finite_pair(self.g0, "g0")
        K = _finite_pair(self.K, "K")
        g0_norm_sq = g0[0] * g0[0] + g0[1] * g0[1]
        if g0_norm_sq <= 1.0e-24:
            raise ValueError("g0 must be nonzero with |g0|^2>1e-24")
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "lambda0", lambda0)
        object.__setattr__(self, "u_star", u_star)
        object.__setattr__(self, "L_s", L_s)
        object.__setattr__(self, "R0", R0)
        object.__setattr__(self, "sigma", sigma)
        object.__setattr__(self, "g0", g0)
        object.__setattr__(self, "K", K)

    @staticmethod
    def nearest_nonzero_integer(target: float) -> int:
        """Autonomous deterministic realization of source's unspecified choice.

        Half ties are sent away from zero.  If ordinary nearest-integer rounding
        would give zero, select the target's sign; exact zero maps to +1.
        """
        x = float(target)
        if not np.isfinite(x):
            raise ValueError("nearest-nonzero target must be finite")
        sign = -1 if x < 0.0 else 1
        magnitude = int(floor(abs(x) + 0.5))
        return sign * max(1, magnitude)

    @property
    def k(self) -> int:
        return int(ceil(self.epsilon ** -0.5))

    @property
    def B_s(self) -> float:
        denominator = self.epsilon * self.k * self.k * (1.0 + self.u_star * self.u_star) ** 1.5
        return sqrt(self.lambda0 / denominator)

    @property
    def x0(self) -> float:
        return self.sigma * self.B_s * self.u_star / 2.0

    @property
    def q(self) -> tuple[float, float]:
        g_r, g_z = self.g0
        K_r, K_z = self.K
        g_norm_sq = g_r * g_r + g_z * g_z
        factor = self.sigma * self.u_star / (self.L_s * g_norm_sq)
        return self.B_s * (K_r - factor * g_r), self.B_s * (K_z - factor * g_z)

    @property
    def tilde_p(self) -> float:
        return self.R0 * self.q[0]

    @property
    def j(self) -> int:
        return self.nearest_nonzero_integer(self.k * self.tilde_p)

    @property
    def p(self) -> float:
        return self.j / self.k

    @property
    def p_z(self) -> float:
        return self.q[1]

    def phase_contract(self, *, m: int = 1) -> KokunoOscillatoryPhaseContract:
        """Build the fixed-label phase contract without inventing base-field data."""
        return KokunoOscillatoryPhaseContract(
            p=self.p,
            p_z=self.p_z,
            x0=self.x0,
            epsilon=self.epsilon,
            m=m,
        )

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "source oscillatory discrete phase-label selection",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "inputs": {
                "epsilon": self.epsilon,
                "lambda0": self.lambda0,
                "u_star": self.u_star,
                "L_s": self.L_s,
                "sigma": self.sigma,
                "R0": self.R0,
                "g0": list(self.g0),
                "K": list(self.K),
                "origin": "bounded replay inputs; not hidden-label recovery",
            },
            "selection": {
                "k": self.k,
                "B_s": self.B_s,
                "x0": self.x0,
                "q": list(self.q),
                "tilde_p": self.tilde_p,
                "j": self.j,
                "p": self.p,
                "p_z": self.p_z,
                "k_times_p": self.k * self.p,
                "tie_rule": "autonomous: half-away-from-zero; zero-nearest uses target sign; exact zero -> +1",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourcePhaseLabelSelection":
        if not isinstance(payload, dict):
            raise ValueError("phase-label payload must be a JSON object")
        if set(payload) - {"sha256"} != {"schema", "source", "inputs", "selection", "truth_boundary"}:
            raise ValueError("phase-label payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported source phase-label schema")
        inputs = payload["inputs"]
        if set(inputs) != {"epsilon", "lambda0", "u_star", "L_s", "sigma", "R0", "g0", "K", "origin"}:
            raise ValueError("phase-label input metadata changed")
        obj = cls(
            epsilon=float(inputs["epsilon"]),
            lambda0=float(inputs["lambda0"]),
            u_star=float(inputs["u_star"]),
            L_s=float(inputs["L_s"]),
            sigma=int(inputs["sigma"]),
            R0=float(inputs["R0"]),
            g0=tuple(inputs["g0"]),
            K=tuple(inputs["K"]),
        )
        expected = obj.to_payload()
        for key in ("source", "inputs", "selection", "truth_boundary"):
            if payload[key] != expected[key]:
                raise ValueError(f"phase-label {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("phase-label payload SHA mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourcePhaseLabelSelection":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
