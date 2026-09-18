"""Bind Kokuno's source-normalized radial derivative to the pulse/curl path.

This module closes one narrow Agent-2 seam after the source-normalized axial
bridge.  The corrected 2026-09-09 reader supplies, on an extended torus chart,

    D_r = partial_R + M_i*d_r*R**(d_r-1)*L_i,
    D_z = epsilon*partial_Z,
    D_r v = D_z v = 0,

with ``L_i=v_r dot partial_{Y_i}``.  The Agent-1 leading-only chart fields used
by the existing phase bridge have no auxiliary ``Y_i`` dependence.  Therefore
``L_i`` annihilates those *geometry coefficients* and source ``D_r`` reduces
exactly to ordinary ``partial_R`` on ``n_Phi``, ``n_Phi'`` and ``K`` in this
restricted leading-only sector.

The radial derivatives below are evaluated with a fixed fourth-order centered
R stencil.  That stencil is an autonomous numerical approximation.  A caller
must still supply ``f_m`` and its full source-normalized ``D_r f_m``.  With
those inputs the existing projected-pulse sensitivity solver identifies its
parameter derivative with ``D_r`` and returns ``D_r C_m``.  This does not
instantiate the actual positive-order background, the actual source pulse
forcing, the auxiliary-dependent source path, a real conjugate pair, or a
public physical ``[u,v,w]`` correction.
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
from .kokuno_source_dz_pulse_bridge import KokunoLeadingSourceDzPulseBridge
from .kokuno_source_pulse_sensitivity import KokunoProjectedPulseSensitivityContract

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-leading-source-dr-pulse-bridge-v1"

_SOURCE_FORMULAS = {
    "source_D_r": "D_r=partial_R+M_i*d_r*R^(d_r-1)*L_i",
    "auxiliary_direction": "L_i=v_r dot partial_{Y_i}",
    "pulse_coordinate": "D_r v=D_z v=0",
    "phase": "Phi=p*theta+p_z*Z/epsilon+x0*R-v*H_Phi; H_Phi=p*F+p_z*G",
    "covector": "n_Phi=(x0-v*(H_Phi)_R,p/R,p_z-epsilon*v*(H_Phi)_Z)",
    "pulse_derivative": "n_Phi'=(-(H_Phi)_R,0,-epsilon*(H_Phi)_Z)",
    "background_matrix": "K=[[0,-2F,0],[2F+R*F_R,0,0],[G_R,0,0]]",
}

_TRUTH_BOUNDARY = {
    "source_D_r_operator_decomposed": True,
    "leading_geometry_auxiliary_independent": True,
    "source_D_r_geometry_reduces_to_partial_R_on_leading_only_bridge": True,
    "leading_only_D_r_C_m_executable_with_caller_supplied_forcing_derivative": True,
    "source_actual_pulse_forcing_instantiated": False,
    "source_actual_background_path_instantiated": False,
    "source_actual_D_r_pulse_path_instantiated": False,
    "source_actual_auxiliary_dependent_D_r_path_instantiated": False,
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
class KokunoLeadingSourceDrPulseBridge:
    """Instantiate source ``D_r`` on auxiliary-independent leading geometry.

    ``ell`` is the bounded autonomous dyadic replay label already used by the
    leading-chart bridge. ``dr_step`` controls only this repository's centered
    finite-difference approximation to ``partial_R``.  Neither is a recovered
    Kokuno/OpenAI hidden parameter.
    """

    ell: int = 64
    dr_step: float = 2.0e-4

    def __post_init__(self) -> None:
        leading = KokunoAgent1LeadingChartPhaseBridge(ell=self.ell)
        object.__setattr__(self, "ell", leading.ell)
        step = float(self.dr_step)
        if not np.isfinite(step) or not (1.0e-6 <= step <= 5.0e-3):
            raise ValueError("dr_step must lie in [1e-6,5e-3]")
        object.__setattr__(self, "dr_step", step)

    @property
    def leading_bridge(self) -> KokunoAgent1LeadingChartPhaseBridge:
        return KokunoAgent1LeadingChartPhaseBridge(ell=self.ell)

    @property
    def _geometry_bridge(self) -> KokunoLeadingSourceDzPulseBridge:
        # Reuse only the already tested source n_Phi/n_Phi'/K assembly.  Its
        # D_z finite-difference method is not used by this module.
        return KokunoLeadingSourceDzPulseBridge(ell=self.ell)

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
        return self._geometry_bridge._raw_geometry(
            candidate,
            phase_contract,
            R,
            Z,
            T,
            pulse_v,
        )

    def source_geometry(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
        R: float,
        Z: float,
        T: float,
        pulse_v: Any,
    ) -> dict[str, np.ndarray | float | str | bool]:
        """Return source pulse geometry and restricted source ``D_r`` derivatives.

        The source operator contains an auxiliary derivative.  The Agent-1
        leading-only geometry has no auxiliary coordinate, so that term is
        exactly zero for the geometry coefficients differentiated here.  The
        remaining ``partial_R`` derivative is approximated by the declared
        fourth-order centered stencil.
        """
        R = _finite_scalar(R, "R")
        Z = _finite_scalar(Z, "Z")
        T = _finite_scalar(T, "T")
        v = _finite_path(pulse_v, "pulse_v")
        h = self.dr_step
        if R <= 2.0 * h:
            raise ValueError("R must exceed 2*dr_step for the away-from-axis D_r stencil")

        center = self._raw_geometry(candidate, phase_contract, R, Z, T, v)
        rm2 = self._raw_geometry(candidate, phase_contract, R - 2.0 * h, Z, T, v)
        rm1 = self._raw_geometry(candidate, phase_contract, R - h, Z, T, v)
        rp1 = self._raw_geometry(candidate, phase_contract, R + h, Z, T, v)
        rp2 = self._raw_geometry(candidate, phase_contract, R + 2.0 * h, Z, T, v)
        factor = 1.0 / (12.0 * h)

        def dr(name: str) -> np.ndarray:
            return factor * (rm2[name] - 8.0 * rm1[name] + 8.0 * rp1[name] - rp2[name])

        return {
            **center,
            "D_r_n_Phi": dr("n_Phi"),
            "D_r_n_Phi_prime": dr("n_Phi_prime"),
            "D_r_K": dr("K"),
            "dr_step": h,
            "D_r_definition": "partial_R+M_i*d_r*R^(d_r-1)*L_i",
            "D_r_geometry_reduction": "L_i=0 on auxiliary-independent leading-only geometry; D_r=partial_R",
            "leading_geometry_auxiliary_independent": True,
            "D_r_discretization": "fourth-order centered chart-R stencil",
        }

    def solve_Dr_coefficient(
        self,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
        R: float,
        Z: float,
        T: float,
        pulse_v: Any,
        forcing: Any,
        D_r_forcing: Any,
    ) -> dict[str, Any]:
        """Solve the leading-only pulse path and return executable ``D_r C_m``.

        ``D_r_forcing`` is the caller's *full source-normalized* radial
        derivative, including any auxiliary contribution if the supplied source
        has one.  This module only reduces ``D_r`` to ``partial_R`` on the
        auxiliary-independent leading geometry coefficients.
        """
        geometry = self.source_geometry(candidate, phase_contract, R, Z, T, pulse_v)
        v = np.asarray(geometry["pulse_v"], dtype=float)
        f = _complex_path3(forcing, "forcing", v.size)
        f_Dr = _complex_path3(D_r_forcing, "D_r_forcing", v.size)
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
            geometry["D_r_n_Phi"],
            geometry["D_r_n_Phi_prime"],
            geometry["D_r_K"],
            f_Dr,
        )
        return {
            **solved,
            "D_r_C_m": solved["C_m_xi"],
            "D_r_t_m": solved["t_m_xi"],
            "geometry": geometry,
            "derivative_identification": (
                "xi is source D_r; geometry sector has L_i=0, while caller supplies full D_r forcing"
            ),
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
                "formula_scope": "Section-6 extended radial operator plus stage-9 pulse geometry",
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "ell": self.ell,
                "dr_step": self.dr_step,
                "epsilon": phase_contract.epsilon,
                "m": phase_contract.m,
                "origin": (
                    "ell and dr_step are bounded autonomous replay/numerical choices; "
                    "they are not recovered hidden/source values"
                ),
            },
            "dependency": {
                "agent1_candidate_sha256": candidate.sha256,
                "leading_bridge_sha256": self.leading_bridge.sha256(candidate),
                "phase_contract_sha256": phase_contract.sha256,
                "source_geometry_schema": "kokuno-leading-source-dz-pulse-bridge-v1",
                "pulse_sensitivity_schema": "kokuno-source-projected-pulse-sensitivity-v1",
            },
            "numerical_method": {
                "D_r": "source operator partial_R+M_i*d_r*R^(d_r-1)*L_i",
                "geometry_sector": "auxiliary-independent leading-only coefficients, hence L_i=0",
                "slow_derivative": "fourth-order centered finite difference in chart R",
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
    ) -> "KokunoLeadingSourceDrPulseBridge":
        if not isinstance(payload, dict):
            raise ValueError("source D_r pulse bridge payload must be an object")
        required = {"schema", "source", "parameters", "dependency", "numerical_method", "truth_boundary"}
        if set(payload) - {"sha256"} != required:
            raise ValueError("source D_r pulse bridge payload keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported source D_r pulse bridge schema")
        params = payload["parameters"]
        if not isinstance(params, dict) or set(params) != {"ell", "dr_step", "epsilon", "m", "origin"}:
            raise ValueError("source D_r pulse bridge parameters changed")
        obj = cls(ell=params["ell"], dr_step=params["dr_step"])
        expected = obj.to_payload(candidate, phase_contract)
        for key in ("source", "parameters", "dependency", "numerical_method", "truth_boundary"):
            if payload[key] != expected[key]:
                raise ValueError(f"source D_r pulse bridge {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256(candidate, phase_contract):
            raise ValueError("source D_r pulse bridge SHA256 mismatch")
        return obj

    @classmethod
    def load_json(
        cls,
        path: str | Path,
        candidate: KokunoLeadingCoreSeriesCandidate,
        phase_contract: KokunoOscillatoryPhaseContract,
    ) -> "KokunoLeadingSourceDrPulseBridge":
        return cls.from_payload(
            json.loads(Path(path).read_text(encoding="utf-8")), candidate, phase_contract
        )
