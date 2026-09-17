"""Bridge Agent-1 leading profiles into Kokuno's oscillatory chart variables.

The corrected 2026-09-09 reader uses the dyadic chart

    R = r / Q**(1/2),  Z = z / Q**D,  T = (1-t)/Q,

with physical velocity ``Q**(-A) u_*`` and writes the chart tangential
components as ``(V,G)`` with ``F_chart = V/R``.  Agent-1's leading profile is
expressed in the native similarity variables ``s=q/Q``, ``X`` and ``eta``:

    u_theta = q**(-A) * sqrt(2X) * F_profile,
    u_z     = q**(-A) * U.

Consequently the *leading-only* chart map is

    F_chart = s**(-(A+1/2)) F_profile,
    V       = R F_chart,
    G       = s**(-A) U.

This module implements that source-implied scaling and its R/Z derivatives,
then feeds it into :class:`KokunoOscillatoryPhaseContract`.  It does NOT claim
to reconstruct Proposition 5.5's actual background: the source states that the
actual chart V,G differ from their leading profiles by positive-order terms.
No complete-curl velocity is changed here.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from math import exp2
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-agent1-leading-chart-phase-bridge-v1"

_SOURCE_FORMULAS = {
    "chart": "R=r/Q^(1/2), Z=z/Q^D, T=(1-t)/Q, s=q/Q",
    "physical_scaling": "physical velocity = Q^(-A) u_*",
    "agent1_leading": "u_theta=q^(-A)*sqrt(2X)*F_profile; u_z=q^(-A)*U",
    "leading_chart": "F_chart=V/R=s^(-(A+1/2))*F_profile; G=s^(-A)*U",
    "epsilon": "epsilon=Q^h=2^(-ell*h)",
    "source_background_boundary": (
        "actual chart b=O(epsilon); actual V,G differ from leading profiles by O(epsilon^2)"
    ),
}

_TRUTH_BOUNDARY = {
    "source_chart_scaling_executable": True,
    "agent1_leading_profile_to_source_leading_V_G_mapping_completed": True,
    "source_phase_composed_with_agent1_leading_profile": True,
    "actual_background_V_G_mapping_completed": False,
    "positive_order_background_corrections_included": False,
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
class KokunoAgent1LeadingChartPhaseBridge:
    """Map an Agent-1 leading core profile into source oscillatory chart data.

    ``ell`` is a bounded dyadic-band label used only to impose the source
    relation ``Q=2^-ell`` and ``epsilon=Q^h`` on the phase contract.  Its
    numerical default is an autonomous replay value, not a recovered label.
    """

    ell: int = 64

    def __post_init__(self) -> None:
        if isinstance(self.ell, bool) or not isinstance(self.ell, (int, np.integer)):
            raise TypeError("ell must be an integer")
        ell = int(self.ell)
        if ell < 1 or ell > 1022:
            raise ValueError("ell must lie in [1,1022]")
        object.__setattr__(self, "ell", ell)

    def epsilon(self, candidate: KokunoLeadingCoreSeriesCandidate) -> float:
        value = float(exp2(-self.ell * float(candidate.h)))
        if not np.isfinite(value) or not (1.0e-6 <= value <= 1.0):
            raise ValueError("Q^h lies outside the phase-contract epsilon bounds")
        return value

    def phase_contract(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        *,
        p: int = 1,
        p_z: float = 1.0,
        x0: float = 1.0,
        m: int = 1,
    ) -> KokunoOscillatoryPhaseContract:
        return KokunoOscillatoryPhaseContract(
            p=p, p_z=p_z, x0=x0, epsilon=self.epsilon(candidate), m=m
        )

    @staticmethod
    def _broadcast_chart(R: Any, Z: Any, T: Any) -> tuple[np.ndarray, ...]:
        R, Z, T = np.broadcast_arrays(
            _finite(R, "R"), _finite(Z, "Z"), _finite(T, "T")
        )
        if not np.all(R > 0.0):
            raise ValueError("oscillatory chart bridge requires R>0")
        if not np.all(T > 0.0):
            raise ValueError("oscillatory chart bridge requires T>0")
        return R, Z, T

    def leading_chart_fields(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        R: Any,
        Z: Any,
        T: Any,
    ) -> dict[str, np.ndarray]:
        """Return leading-only ``V,G`` and exact chart R/Z derivatives.

        The normalized source relation for ``s=q/Q`` is
        ``T=s-Z^2*s^(2h)``.  Reusing Agent-1's native coordinate solver at
        normalized inputs ``(R,0,Z,1-T)`` evaluates this same relation without
        introducing a second root solver.
        """
        R, Z, T = self._broadcast_chart(R, Z, T)
        coords = candidate.coordinates_model.evaluate(R, 0.0, Z, 1.0 - T)
        s = coords["q"]
        eta = coords["eta"]
        X = coords["X"]
        if np.any(candidate.Lambda * X > 4.1 + 2.0e-12):
            raise ValueError("leading chart bridge requires Agent-1 core Lambda*X <= 4.1")

        profile = candidate.profile_values(X, eta)
        A = 0.5 + candidate.h
        alpha = A + 0.5
        s_inv = 1.0 / s
        s_minus_alpha = np.power(s, -alpha)
        s_minus_A = np.power(s, -A)

        Fp = profile["F"]
        Fp_X = profile["F_X"]
        Fp_eta = profile["F_eta"]
        U = profile["U"]
        U_X = profile["U_X"]
        U_eta = profile["U_eta"]

        # Agent-1 normalized coordinate Jacobian at (R,0,Z,1-T) is exactly
        # the source chart Jacobian of (s,eta,X) with respect to Z at fixed T.
        s_Z = coords["q_z"]
        eta_Z = coords["eta_z"]
        X_Z = coords["X_z"]
        X_R = R * s_inv

        F_chart = s_minus_alpha * Fp
        F_chart_R = s_minus_alpha * Fp_X * X_R
        F_chart_Z = s_minus_alpha * (
            Fp_X * X_Z + Fp_eta * eta_Z - alpha * Fp * s_Z * s_inv
        )
        V = R * F_chart
        V_R = F_chart + R * F_chart_R
        V_Z = R * F_chart_Z

        G = s_minus_A * U
        G_R = s_minus_A * U_X * X_R
        G_Z = s_minus_A * (
            U_X * X_Z + U_eta * eta_Z - A * U * s_Z * s_inv
        )

        arrays = (F_chart, F_chart_R, F_chart_Z, V, V_R, V_Z, G, G_R, G_Z)
        if not all(np.all(np.isfinite(a)) for a in arrays):
            raise RuntimeError("leading chart bridge produced a nonfinite value")
        return {
            "s": s,
            "eta": eta,
            "X": X,
            "F_profile": Fp,
            "U_profile": U,
            "F_chart": F_chart,
            "F_chart_R": F_chart_R,
            "F_chart_Z": F_chart_Z,
            "V": V,
            "V_R": V_R,
            "V_Z": V_Z,
            "G": G,
            "G_R": G_R,
            "G_Z": G_Z,
        }

    def evaluate_phase_on_leading(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        contract: KokunoOscillatoryPhaseContract,
        R: Any,
        theta: Any,
        Z: Any,
        T: Any,
        pulse_v: Any,
    ) -> dict[str, np.ndarray | int]:
        """Evaluate #262 source phase/covector on the leading-only chart map."""
        expected_epsilon = self.epsilon(candidate)
        if not np.isclose(contract.epsilon, expected_epsilon, rtol=2.0e-15, atol=0.0):
            raise ValueError("phase epsilon must equal Q^h for this dyadic bridge")
        fields = self.leading_chart_fields(candidate, R, Z, T)
        return contract.evaluate_from_V(
            R,
            theta,
            Z,
            pulse_v,
            fields["V"],
            fields["V_R"],
            fields["V_Z"],
            fields["G"],
            fields["G_R"],
            fields["G_Z"],
        )

    def to_payload(self, candidate: KokunoLeadingCoreSeriesCandidate) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "dyadic-chart leading tangential fields feeding stage-9 phase",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "ell": self.ell,
                "Q": float(exp2(-self.ell)),
                "epsilon": self.epsilon(candidate),
                "h": float(candidate.h),
                "origin": "ell is a bounded autonomous replay label; no hidden/source label recovered",
            },
            "dependency": {
                "agent1_candidate_schema": candidate.to_payload()["schema"],
                "agent1_candidate_sha256": candidate.sha256,
                "agent2_phase_schema": "kokuno-oscillatory-phase-contract-v1",
                "mapping_scope": "leading chart V,G only; not Proposition-5.5 corrected background",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    def sha256(self, candidate: KokunoLeadingCoreSeriesCandidate) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload(candidate)).encode("utf-8")).hexdigest()

    def save_json(
        self, path: str | Path, candidate: KokunoLeadingCoreSeriesCandidate
    ) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload(candidate)
        payload["sha256"] = self.sha256(candidate)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target
