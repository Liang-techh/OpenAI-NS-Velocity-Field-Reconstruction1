"""Bind Kokuno's source-normalized axial derivative to the pulse/curl path.

This module closes one narrow Agent-2 seam without claiming the missing actual
background or pulse source.  The corrected 2026-09-09 reader supplies

    D_z = epsilon * partial_Z,
    Phi = p*theta + p_z*Z/epsilon + x0*R - v*H_Phi,
    n_Phi = (x0-v*(H_Phi)_R, p/R, p_z-epsilon*v*(H_Phi)_Z),
    n_Phi' = (-(H_Phi)_R, 0, -epsilon*(H_Phi)_Z),

and, on the tangential background used by the projected pulse equation,

    K = [[0, -2F, 0],
         [2F + R F_R, 0, 0],
         [G_R, 0, 0]].

The Agent-1 bridge already evaluates the *leading-only* chart fields F,G and
the source phase/covector.  Here those sourced formulas are assembled along an
increasing pulse coordinate v.  Source-normalized D_z derivatives of n_Phi,
n_Phi' and K are evaluated with a fixed fourth-order centered Z stencil; that
stencil is an autonomous numerical approximation, not a paper formula.

A caller must still supply f_m and D_z f_m.  With those inputs the existing
projected-pulse sensitivity solver can identify xi with the source operator D_z
and return D_z C_m for the already tested complete-curl coefficient.  Actual
positive-order background corrections, actual source pulse forcing, D_r C_m,
real conjugate-pair/Q-scaled velocity, and public [u,v,w] materialization all
remain unresolved.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_oscillatory_leading_bridge import KokunoAgent1LeadingChartPhaseBridge
from .kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract
from .kokuno_source_pulse_sensitivity import KokunoProjectedPulseSensitivityContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-leading-source-dz-pulse-bridge-v1"

_SOURCE_FORMULAS = {
    "source_D_z": "D_z=epsilon*partial_Z",
    "phase": "Phi=p*theta+p_z*Z/epsilon+x0*R-v*H_Phi; H_Phi=p*F+p_z*G",
    "covector": "n_Phi=(x0-v*(H_Phi)_R,p/R,p_z-epsilon*v*(H_Phi)_Z)",
    "pulse_derivative": "n_Phi'=(-(H_Phi)_R,0,-epsilon*(H_Phi)_Z)",
    "background_matrix": "K=[[0,-2F,0],[2F+R*F_R,0,0],[G_R,0,0]]",
    "pulse_coordinate": "D_r v=D_z v=0",
}

_TRUTH_BOUNDARY = {
    "source_D_z_operator_bound_to_leading_background": True,
    "source_n_Phi_prime_formula_executable": True,
    "source_background_K_formula_executable_on_leading_fields": True,
    "leading_only_D_z_C_m_executable_with_caller_supplied_forcing_derivative": True,
    "source_actual_pulse_forcing_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "source_actual_D_z_pulse_path_instantiated": False,
    "source_actual_D_r_path_instantiated": False,
    "actual_background_V_G_mapping_completed": False,
    "positive_order_background_corrections_included": False,
    "real_conjugate_pair_materialized": False,
    "physical_Q_scaled_velocity_materialized": False,
    "public_velocity_correction_materialized": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_scalar(value: Any, name: str) -> float:
    out = float(value)
    if not np.isfinite(out):
        raise ValueError(f"{name} must be finite")
    return out


def _finite_path(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim != 1 or out.size < 2:
        raise ValueError(f"{name} must be a one-dimensional path with at least two nodes")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _complex_path3(value: Any, name: str, nodes: int) -> np.ndarray:
    out = np.asarray(value, dtype=np.complex128)
    if out.shape != (nodes, 3):
        raise ValueError(f"{name} must have shape ({nodes},3)")
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoLeadingSourceDzPulseBridge:
    """Instantiate source ``D_z`` geometry on the Agent-1 leading-only chart.

    ``ell`` is the same bounded autonomous dyadic replay label used by the
    leading-chart bridge. ``dz_step`` controls only the repository's centered
    finite-difference approximation to source-normalized coefficient
    derivatives.  Neither value is a recovered Kokuno/OpenAI hidden parameter.
    """

    ell: int = 64
    dz_step: float = 2.0e-4

    def __post_init__(self) -> None:
        leading = KokunoAgent1LeadingChartPhaseBridge(ell=self.ell)
        object.__setattr__(self, "ell", leading.ell)
        step = float(self.dz_step)
        if not np.isfinite(step) or not (1.0e-6 <= step <= 5.0e-3):
            raise ValueError("dz_step must lie in [1e-6,5e-3]")
        object.__setattr__(self, "dz_step", step)

    @property
    def leading_bridge(self) -> KokunoAgent1LeadingChartPhaseBridge:
        return KokunoAgent1LeadingChartPhaseBridge(ell=self.ell)

    def _validate_binding(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
    ) -> None:
        expected = self.leading_bridge.epsilon(candidate)
        if not np.isclose(phase_contract.epsilon, expected, rtol=2.0e-15, atol=0.0):
            raise ValueError("phase epsilon must equal Q^h for the bound leading chart")

    def _raw_geometry(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
        R: float,
        Z: float,
        T: float,
        pulse_v: np.ndarray,
    ) -> dict[str, np.ndarray]:
        self._validate_binding(candidate, phase_contract)
        R = _finite_scalar(R, "R")
        Z = _finite_scalar(Z, "Z")
        T = _finite_scalar(T, "T")
        if R <= 0.0:
            raise ValueError("source oscillatory chart requires R>0")
        if T <= 0.0:
            raise ValueError("source oscillatory chart requires T>0")
        v = _finite_path(pulse_v, "pulse_v")
        if np.any(np.diff(v) <= 0.0):
            raise ValueError("pulse_v nodes must be strictly increasing")

        fields = self.leading_bridge.leading_chart_fields(candidate, R, Z, T)
        phase = phase_contract.evaluate_from_V(
            R,
            0.0,
            Z,
            v,
            fields["V"],
            fields["V_R"],
            fields["V_Z"],
            fields["G"],
            fields["G_R"],
            fields["G_Z"],
        )
        N = v.size
        H_R = np.broadcast_to(np.asarray(phase["H_Phi_R"], dtype=float), (N,))
        H_Z = np.broadcast_to(np.asarray(phase["H_Phi_Z"], dtype=float), (N,))
        n_prime = np.stack(
            (-H_R, np.zeros(N, dtype=float), -phase_contract.epsilon * H_Z),
            axis=-1,
        )

        F = np.broadcast_to(np.asarray(phase["F"], dtype=float), (N,))
        F_R = np.broadcast_to(np.asarray(phase["F_R"], dtype=float), (N,))
        G_R = np.broadcast_to(np.asarray(fields["G_R"], dtype=float), (N,))
        K = np.zeros((N, 3, 3), dtype=float)
        K[:, 0, 1] = -2.0 * F
        K[:, 1, 0] = 2.0 * F + R * F_R
        K[:, 2, 0] = G_R

        return {
            "pulse_v": v.copy(),
            "Phi_theta0": np.asarray(phase["Phi"], dtype=float),
            "n_Phi": np.asarray(phase["n_Phi"], dtype=float),
            "n_Phi_prime": n_prime,
            "K": K,
            "F": F,
            "F_R": F_R,
            "G_R": G_R,
        }

    def source_geometry(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
        R: float,
        Z: float,
        T: float,
        pulse_v: Any,
    ) -> dict[str, np.ndarray | float | str]:
        """Return sourced pulse geometry and source-normalized ``D_z`` derivatives.

        The source formulas for ``n_Phi'`` and ``K`` are exact on the supplied
        leading-only fields.  ``D_z n_Phi``, ``D_z n_Phi'`` and ``D_z K`` use
        the declared fourth-order chart-Z stencil and are therefore numerical
        approximations to the sourced operator ``epsilon*partial_Z``.
        """
        v = _finite_path(pulse_v, "pulse_v")
        center = self._raw_geometry(candidate, phase_contract, R, Z, T, v)
        h = self.dz_step
        fm2 = self._raw_geometry(candidate, phase_contract, R, Z - 2.0 * h, T, v)
        fm1 = self._raw_geometry(candidate, phase_contract, R, Z - h, T, v)
        fp1 = self._raw_geometry(candidate, phase_contract, R, Z + h, T, v)
        fp2 = self._raw_geometry(candidate, phase_contract, R, Z + 2.0 * h, T, v)
        factor = phase_contract.epsilon / (12.0 * h)

        def dz(name: str) -> np.ndarray:
            return factor * (fm2[name] - 8.0 * fm1[name] + 8.0 * fp1[name] - fp2[name])

        return {
            **center,
            "D_z_n_Phi": dz("n_Phi"),
            "D_z_n_Phi_prime": dz("n_Phi_prime"),
            "D_z_K": dz("K"),
            "dz_step": h,
            "D_z_definition": "epsilon*partial_Z",
            "D_z_discretization": "fourth-order centered chart-Z stencil",
        }

    def solve_Dz_coefficient(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
        R: float,
        Z: float,
        T: float,
        pulse_v: Any,
        forcing: Any,
        D_z_forcing: Any,
    ) -> dict[str, Any]:
        """Solve the leading-only pulse path and return executable ``D_z C_m``.

        ``forcing`` and ``D_z_forcing`` are caller-supplied source-normalized
        pulse-source data.  Supplying them does not mark the actual source
        forcing as reconstructed.  No public velocity is assembled here.
        """
        geometry = self.source_geometry(candidate, phase_contract, R, Z, T, pulse_v)
        v = np.asarray(geometry["pulse_v"], dtype=float)
        f = _complex_path3(forcing, "forcing", v.size)
        f_Dz = _complex_path3(D_z_forcing, "D_z_forcing", v.size)
        sensitivity = KokunoProjectedPulseSensitivityContract(
            epsilon=phase_contract.epsilon,
            m=phase_contract.m,
        )
        solved = sensitivity.solve_path(
            v,
            geometry["n_Phi"],
            geometry["n_Phi_prime"],
            geometry["K"],
            f,
            geometry["D_z_n_Phi"],
            geometry["D_z_n_Phi_prime"],
            geometry["D_z_K"],
            f_Dz,
        )
        return {
            **solved,
            "D_z_C_m": solved["C_m_xi"],
            "D_z_t_m": solved["t_m_xi"],
            "geometry": geometry,
            "derivative_identification": "xi is source-normalized D_z=epsilon*partial_Z",
            "forcing_origin": "caller supplied; not asserted to be the actual Kokuno pulse forcing",
        }

    def to_payload(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
    ) -> dict[str, Any]:
        self._validate_binding(candidate, phase_contract)
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "formula_scope": "stage-9 pulse geometry, n_Phi prime, K, and source-normalized D_z",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "ell": self.ell,
                "dz_step": self.dz_step,
                "epsilon": phase_contract.epsilon,
                "m": phase_contract.m,
                "origin": (
                    "ell and dz_step are bounded autonomous replay/numerical choices; "
                    "they are not recovered hidden/source values"
                ),
            },
            "dependency": {
                "agent1_candidate_sha256": candidate.sha256,
                "leading_bridge_sha256": self.leading_bridge.sha256(candidate),
                "phase_contract_sha256": phase_contract.sha256,
                "pulse_sensitivity_schema": "kokuno-source-projected-pulse-sensitivity-v1",
            },
            "numerical_method": {
                "D_z": "source operator epsilon*partial_Z",
                "slow_derivative": "fourth-order centered finite difference in chart Z",
                "pulse_solve": "delegated unchanged to KokunoProjectedPulseSensitivityContract",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    def sha256(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
    ) -> str:
        return hashlib.sha256(
            _canonical_json(self.to_payload(candidate, phase_contract)).encode("utf-8")
        ).hexdigest()

    def save_json(
        self,
        path: str | Path,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
    ) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload(candidate, phase_contract)
        payload["sha256"] = self.sha256(candidate, phase_contract)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(
        cls,
        payload: dict[str, Any],
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
    ) -> "KokunoLeadingSourceDzPulseBridge":
        if not isinstance(payload, dict):
            raise ValueError("source D_z pulse bridge payload must be an object")
        required = {"schema", "source", "parameters", "dependency", "numerical_method", "truth_boundary"}
        if set(payload) - {"sha256"} != required:
            raise ValueError("source D_z pulse bridge payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported source D_z pulse bridge schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {"ell", "dz_step", "epsilon", "m", "origin"}:
            raise ValueError("source D_z pulse bridge parameters changed")
        obj = cls(ell=params["ell"], dz_step=params["dz_step"])
        expected = obj.to_payload(candidate, phase_contract)
        for key in ("source", "parameters", "dependency", "numerical_method", "truth_boundary"):
            if payload[key] != expected[key]:
                raise ValueError(f"source D_z pulse bridge {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256(candidate, phase_contract):
            raise ValueError("source D_z pulse bridge SHA256 mismatch")
        return obj

    @classmethod
    def load_json(
        cls,
        path: str | Path,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
    ) -> "KokunoLeadingSourceDzPulseBridge":
        return cls.from_payload(
            json.loads(Path(path).read_text(encoding="utf-8")), candidate, phase_contract
        )
