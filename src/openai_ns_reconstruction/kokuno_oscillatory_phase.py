"""Executable source phase/covector contract for Kokuno oscillatory labels.

The corrected 2026-09-09 reconstruction retains the full normalized cylindrical
phase

    Phi = p*theta + p_z*Z/epsilon + x0*R - v*H_Phi,
    H_Phi = p*F + p_z*G,    F = V/R,

and the corresponding (non-unit-normalized) covector

    n_Phi = (x0 - v*d_R H_Phi,
             p/R,
             p_z - epsilon*v*d_Z H_Phi).

The source discrete label is not generally integer ``p``: it is ``p=j/k`` with
``j`` a nonzero integer and ``k=ceil(epsilon^(-1/2))``.  Periodicity is thus
encoded by ``k*p=j in Z\{0}``.  This module enforces that source-compatible
rational grid but does not itself choose ``j``; see ``kokuno_source_phase_label``
for the source label-selection algebra.

This module implements only that source algebra. The base-field values V, G and
their R/Z derivatives are supplied by the caller. It deliberately does not map
Agent-1 profile symbols onto source V/G, does not build an oscillatory velocity,
and does not identify any numerical label as a hidden source/OpenAI parameter.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from math import ceil
from pathlib import Path
from typing import Any

import numpy as np

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-oscillatory-phase-contract-v2"

_SOURCE_FORMULAS = {
    "phase": "Phi=p*theta+p_z*Z/epsilon+x0*R-v*H_Phi",
    "H_Phi": "H_Phi=p*F+p_z*G",
    "F": "F=V/R",
    "n_Phi": "(x0-v*d_R(H_Phi), p/R, p_z-epsilon*v*d_Z(H_Phi))",
    "k": "ceil(epsilon^(-1/2))",
    "angular_label_grid": "p=j/k with j integer nonzero (equivalently k*p in Z\\{0})",
    "harmonic": "k_m=k*m, m integer nonzero",
}

_TRUTH_BOUNDARY = {
    "source_phase_formula_executable": True,
    "source_covector_formula_executable": True,
    "source_rational_angular_label_grid_enforced": True,
    "caller_supplies_source_base_field_symbols": True,
    "agent1_profile_to_source_V_G_mapping_completed": False,
    "complete_curl_velocity_changed": False,
    "phase_label_selected": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True)
class KokunoOscillatoryPhaseContract:
    """Bounded executable adapter for one fixed source oscillatory label.

    ``R, theta, Z`` are source normalized cylindrical coordinates. The evaluator
    requires ``R>0`` because the source covector contains ``p/R``; the source
    wave boxes themselves are supported away from the axis.

    ``p`` must lie on the source grid ``p=j/k`` with nonzero integer ``j``.
    Numeric defaults remain autonomous diagnostic label values, not recovered
    source parameters.
    """

    p: float = 1.0
    p_z: float = 1.0
    x0: float = 1.0
    epsilon: float = 0.25
    m: int = 1

    def __post_init__(self) -> None:
        p = float(self.p)
        m = int(self.m)
        p_z = float(self.p_z)
        x0 = float(self.x0)
        epsilon = float(self.epsilon)
        if not np.isfinite(p):
            raise ValueError("p must be finite")
        if m != self.m or m == 0 or abs(m) > 32:
            raise ValueError("m must be a nonzero integer with |m|<=32")
        if not np.isfinite(p_z) or abs(p_z) > 64.0:
            raise ValueError("p_z must be finite with |p_z|<=64")
        if not np.isfinite(x0) or abs(x0) > 64.0:
            raise ValueError("x0 must be finite with |x0|<=64")
        if not np.isfinite(epsilon) or not (1.0e-6 <= epsilon <= 1.0):
            raise ValueError("epsilon must lie in [1e-6,1]")
        k = int(ceil(epsilon ** -0.5))
        j_float = k * p
        j = int(round(j_float))
        if j == 0 or abs(j) > 256 or not np.isclose(j_float, j, rtol=0.0, atol=2.0e-12):
            raise ValueError("p must satisfy p=j/k for a nonzero integer j with |j|<=256")
        object.__setattr__(self, "p", p)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "p_z", p_z)
        object.__setattr__(self, "x0", x0)
        object.__setattr__(self, "epsilon", epsilon)

    @property
    def k(self) -> int:
        return int(ceil(self.epsilon ** -0.5))

    @property
    def j(self) -> int:
        return int(round(self.k * self.p))

    @property
    def k_m(self) -> int:
        return self.k * self.m

    def _broadcast(
        self,
        R: Any,
        theta: Any,
        Z: Any,
        pulse_v: Any,
        V: Any,
        V_R: Any,
        V_Z: Any,
        G: Any,
        G_R: Any,
        G_Z: Any,
    ) -> tuple[np.ndarray, ...]:
        arrays = np.broadcast_arrays(
            _finite(R, "R"),
            _finite(theta, "theta"),
            _finite(Z, "Z"),
            _finite(pulse_v, "pulse_v"),
            _finite(V, "V"),
            _finite(V_R, "V_R"),
            _finite(V_Z, "V_Z"),
            _finite(G, "G"),
            _finite(G_R, "G_R"),
            _finite(G_Z, "G_Z"),
        )
        if not np.all(arrays[0] > 0.0):
            raise ValueError("source oscillatory phase requires R>0 (wave shell is away from axis)")
        return tuple(np.asarray(a, dtype=float) for a in arrays)

    def evaluate_from_V(
        self,
        R: Any,
        theta: Any,
        Z: Any,
        pulse_v: Any,
        V: Any,
        V_R: Any,
        V_Z: Any,
        G: Any,
        G_R: Any,
        G_Z: Any,
    ) -> dict[str, np.ndarray | int]:
        """Evaluate ``Phi`` and ``n_Phi`` from caller-supplied ``V,G`` data.

        The quotient derivatives are exact:
        ``F=V/R``, ``F_R=V_R/R-V/R^2``, ``F_Z=V_Z/R``.
        """
        R, theta, Z, pulse_v, V, V_R, V_Z, G, G_R, G_Z = self._broadcast(
            R, theta, Z, pulse_v, V, V_R, V_Z, G, G_R, G_Z
        )
        F = V / R
        F_R = V_R / R - V / (R * R)
        F_Z = V_Z / R
        H = self.p * F + self.p_z * G
        H_R = self.p * F_R + self.p_z * G_R
        H_Z = self.p * F_Z + self.p_z * G_Z

        Phi = self.p * theta + self.p_z * Z / self.epsilon + self.x0 * R - pulse_v * H
        n_R = self.x0 - pulse_v * H_R
        n_theta = self.p / R
        n_Z = self.p_z - self.epsilon * pulse_v * H_Z
        n = np.stack((n_R, n_theta, n_Z), axis=-1)
        n_norm_sq = np.sum(n * n, axis=-1)
        if not np.all(np.isfinite(Phi)) or not np.all(np.isfinite(n)):
            raise RuntimeError("source phase/covector evaluation produced a nonfinite result")
        if not np.all(n_norm_sq > 0.0):
            raise RuntimeError("source phase covector unexpectedly vanished")
        return {
            "Phi": Phi,
            "H_Phi": H,
            "H_Phi_R": H_R,
            "H_Phi_Z": H_Z,
            "F": F,
            "F_R": F_R,
            "F_Z": F_Z,
            "n_Phi": n,
            "n_Phi_norm_sq": n_norm_sq,
            "k": self.k,
            "j": self.j,
            "k_m": self.k_m,
        }

    def phase_from_V(
        self,
        R: Any,
        theta: Any,
        Z: Any,
        pulse_v: Any,
        V: Any,
        G: Any,
    ) -> np.ndarray:
        """Evaluate only the source phase when derivative data are unnecessary."""
        R, theta, Z, pulse_v, V, G = np.broadcast_arrays(
            _finite(R, "R"),
            _finite(theta, "theta"),
            _finite(Z, "Z"),
            _finite(pulse_v, "pulse_v"),
            _finite(V, "V"),
            _finite(G, "G"),
        )
        if not np.all(R > 0.0):
            raise ValueError("source oscillatory phase requires R>0 (wave shell is away from axis)")
        F = V / R
        H = self.p * F + self.p_z * G
        return self.p * theta + self.p_z * Z / self.epsilon + self.x0 * R - pulse_v * H

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "stage-9 fixed-label oscillatory phase/covector",
                "source_formulas": dict(_SOURCE_FORMULAS),
                "axis_contract": "source wave shell is away from R=0; evaluator requires R>0",
            },
            "parameters": {
                "p": self.p,
                "j": self.j,
                "p_z": self.p_z,
                "x0": self.x0,
                "epsilon": self.epsilon,
                "m": self.m,
                "k": self.k,
                "k_m": self.k_m,
                "origin": "bounded source-grid-compatible label values; no hidden-parameter recovery",
            },
            "caller_contract": {
                "required_base_data": ["V", "V_R", "V_Z", "G", "G_R", "G_Z"],
                "derived": ["F=V/R", "F_R=V_R/R-V/R^2", "F_Z=V_Z/R"],
                "angular_label_grid": "p=j/k, j integer nonzero",
                "agent1_symbol_mapping": "not implemented in this increment",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoOscillatoryPhaseContract":
        if not isinstance(payload, dict):
            raise ValueError("phase payload must be a JSON object")
        if set(payload) - {"sha256"} != {"schema", "source", "parameters", "caller_contract", "truth_boundary"}:
            raise ValueError("phase payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported oscillatory phase schema")
        parameters = payload["parameters"]
        if set(parameters) != {"p", "j", "p_z", "x0", "epsilon", "m", "k", "k_m", "origin"}:
            raise ValueError("phase parameter metadata changed")
        obj = cls(
            p=float(parameters["p"]),
            p_z=float(parameters["p_z"]),
            x0=float(parameters["x0"]),
            epsilon=float(parameters["epsilon"]),
            m=int(parameters["m"]),
        )
        expected = obj.to_payload()
        if payload["source"] != expected["source"]:
            raise ValueError("phase source/provenance metadata changed")
        if payload["caller_contract"] != expected["caller_contract"]:
            raise ValueError("phase caller contract changed")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("phase truth boundary changed")
        for key in ("j", "k", "k_m", "origin"):
            if parameters[key] != expected["parameters"][key]:
                raise ValueError("phase derived parameter metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("phase payload SHA mismatch")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoOscillatoryPhaseContract":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
