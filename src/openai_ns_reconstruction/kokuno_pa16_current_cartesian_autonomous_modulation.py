"""Bind the current A1 RF40 candidate to the existing autonomous Appendix-C modulation.

Pinned provenance: KokunoYumeto corrected 2026-09-09 reconstruction,
commit 143f6773feb424ad9ed3a8d116653200f20346b7, Appendix C.

The public reconstruction inserts
    E_N = E exp(A/N),  U_N = U + B/N,  phi = N log X (mod 1),
before the reserved I1 repair.  The repository already has an explicitly
AUTONOMOUS source-form realization in ``KokunoRadialModulationDiscrepancy``;
the source proves existence of a suitable loop but does not publish its hidden
numerical realization.  This adapter does not rename that autonomous loop as
source-exact.  It only binds the existing realization to the *current* A1
candidate from ``KokunoPA16CurrentCartesianRF40PowerLaw``.

The important executable seam is incompressibility memory.  The current A1
candidate already carries its own prefix M=int_0^X U ds.  The modulation adds
Delta U only on a compact pre-I1 interval, so this class adds
    Delta(M/X) = int Delta U(exp s,eta) exp(s-log X) ds
and the corresponding eta derivative to the current primitive before rebuilding
v0.  Evaluation fails closed at the start of I1; no repair is silently skipped.

This is a candidate-side autonomous overlay, not paper-exact reconstruction,
not independent PDE validation, and not a residual-improvement claim.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa16_current_cartesian_rf40_power_law import (
    KokunoPA16CurrentCartesianRF40PowerLaw,
)
from .kokuno_radial_modulation_discrepancy import (
    KokunoRadialModulationDiscrepancy,
)

SCHEMA = "kokuno-pa16-current-cartesian-autonomous-modulation-v1"

_SOURCE_FORMULAS = {
    "finite_N_insertion": "E_N=E exp(A/N); U_N=U+B/N; phi=N log X mod 1",
    "prefix_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
}
_NUMERICAL_REALIZATION = {
    "loop": (
        "reuse repository-autonomous compact periodic A/B realization from "
        "KokunoRadialModulationDiscrepancy; it is not a recovered source loop"
    ),
    "memory_quadrature": (
        "deterministic Gauss-Legendre panels in log X, aligned at "
        "1/(frequency*panels_per_period)"
    ),
    "current_base": "consume exact A1 #1005 current primitive and Cartesian map",
    "domain_stop": "fail closed at the left endpoint of reserved I1",
    "new_residual_tuning": "none",
}
_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_RF40_base_backbone_consumed": True,
    "source_finite_frequency_insertion_formula_reused": True,
    "repository_autonomous_modulation_bound_to_current_lineage": True,
    "modulation_prefix_M_memory_carried": True,
    "cartesian_velocity_executable_through_pre_I1_domain": True,
    "velocity_interface_vectorized": True,
    "configuration_serializable": True,
    "source_hidden_loop_parameters_recovered": False,
    "source_admissible_loop_reconstructed": False,
    "actual_source_admissible_loop_discrepancy_supplied": False,
    "cone_modulation_completed": False,
    "I1_repair_applied_to_current_modulation": False,
    "I2_overlay_applied_to_current_lineage": False,
    "I3_overlay_applied_to_current_lineage": False,
    "I4_overlay_applied_to_current_lineage": False,
    "source_terminal_tail_schedule_bound": False,
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
class KokunoPA16CurrentCartesianAutonomousModulation:
    """Current A1 Cartesian leading candidate with the selected autonomous pre-I1 loop."""

    current: KokunoPA16CurrentCartesianRF40PowerLaw = field(
        default_factory=KokunoPA16CurrentCartesianRF40PowerLaw,
        repr=False,
        compare=False,
    )
    frequency: int = 16
    alpha: float = 5.0e-4
    beta: float = 5.0e-4
    support_left_offset: float = -35.3
    support_right_offset: float = -29.1
    quadrature_order: int = 8
    panels_per_period: int = 3

    _modulation: KokunoRadialModulationDiscrepancy = field(
        init=False, repr=False, compare=False
    )
    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.current, KokunoPA16CurrentCartesianRF40PowerLaw):
            raise TypeError("current must be KokunoPA16CurrentCartesianRF40PowerLaw")
        modulation = KokunoRadialModulationDiscrepancy(
            outer_schedule=self.current.outer_base.outer_schedule,
            frequency=self.frequency,
            alpha=self.alpha,
            beta=self.beta,
            support_left_offset=self.support_left_offset,
            support_right_offset=self.support_right_offset,
            quadrature_order=self.quadrature_order,
            panels_per_period=self.panels_per_period,
        )
        nodes, weights = leggauss(int(modulation.quadrature_order))
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)

        left, right = modulation.support_log_interval
        i1_left, _ = modulation.outer_schedule.reserved_log_intervals()["I1"]
        if not (
            self.current.outer_base.log_X_w
            < left
            < right
            < i1_left
            < self.current.log_X_4
        ):
            raise ValueError("autonomous support/I1 ordering is inconsistent with current RF40")
        if not math.isclose(
            self.current.log_X_4,
            modulation._outer_base.log_X_end,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("current RF40 endpoint and modulation outer schedule disagree")

        object.__setattr__(self, "frequency", int(modulation.frequency))
        object.__setattr__(self, "alpha", float(modulation.alpha))
        object.__setattr__(self, "beta", float(modulation.beta))
        object.__setattr__(self, "support_left_offset", float(modulation.support_left_offset))
        object.__setattr__(self, "support_right_offset", float(modulation.support_right_offset))
        object.__setattr__(self, "quadrature_order", int(modulation.quadrature_order))
        object.__setattr__(self, "panels_per_period", int(modulation.panels_per_period))
        object.__setattr__(self, "_modulation", modulation)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)

    @property
    def modulation(self) -> KokunoRadialModulationDiscrepancy:
        return self._modulation

    @property
    def outer_base(self):
        return self.current.outer_base

    @property
    def geometry(self):
        return self.current.geometry

    @property
    def A(self) -> float:
        return float(self.current.A)

    @property
    def D(self) -> float:
        return float(self.current.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.current.eta_interval)

    @property
    def support_log_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.modulation.support_log_interval)

    @property
    def log_X_I1_start(self) -> float:
        return float(self.modulation.outer_schedule.reserved_log_intervals()["I1"][0])

    @property
    def X_I1_start(self) -> float:
        return float(math.exp(self.log_X_I1_start))

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.current.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_I1_start)):
            raise ValueError(
                "X must lie in the current autonomous-modulation domain "
                f"[0,{self.X_I1_start:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _panel_integral_pair(self, y: float, eta: float) -> tuple[float, float]:
        """Return Delta(M/X), Delta(M_eta/X) at one positive log-radius."""
        left, right = self.support_log_interval
        if y <= left:
            return 0.0, 0.0
        upper = min(float(y), right)
        if upper <= left:
            return 0.0, 0.0

        panel = 1.0 / (float(self.frequency) * float(self.panels_per_period))
        count = max(1, int(math.ceil((upper - left) / panel)))
        edges = np.linspace(left, upper, count + 1)
        total = 0.0
        total_eta = 0.0
        for a, b in zip(edges[:-1], edges[1:]):
            mid = 0.5 * (a + b)
            half = 0.5 * (b - a)
            samples = mid + half * self._nodes
            etas = np.full(samples.shape, float(eta), dtype=float)
            modified = self.modulation.profile_values_logX(samples, etas)
            base = self.outer_base.profile_values_logX(samples, etas)
            weight = np.exp(samples - float(y))
            delta_u = np.asarray(modified["U"], dtype=float) - np.asarray(
                base["U"], dtype=float
            )
            delta_u_eta = np.asarray(modified["U_eta"], dtype=float) - np.asarray(
                base["U_eta"], dtype=float
            )
            total += half * float(np.dot(self._weights, delta_u * weight))
            total_eta += half * float(
                np.dot(self._weights, delta_u_eta * weight)
            )
        return total, total_eta

    def _delta_memory_ratio(
        self, X: np.ndarray, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        out = np.zeros_like(X, dtype=float)
        out_eta = np.zeros_like(X, dtype=float)
        flat_X = X.reshape(-1)
        flat_eta = eta.reshape(-1)
        flat_out = out.reshape(-1)
        flat_out_eta = out_eta.reshape(-1)
        left, _ = self.support_log_interval
        for i, (xv, ev) in enumerate(zip(flat_X, flat_eta)):
            if xv <= 0.0:
                continue
            y = math.log(float(xv))
            if y <= left:
                continue
            dm, dm_eta = self._panel_integral_pair(y, float(ev))
            flat_out[i] = dm
            flat_out_eta[i] = dm_eta
        return out, out_eta

    def similarity_profile_values(
        self, X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)

        base = self.current.similarity_profile_values(xf, ef)
        F = np.asarray(base["F_current_rf40_power_law"], dtype=float).copy()
        U = np.asarray(base["U_current_rf40_power_law"], dtype=float).copy()
        E = np.asarray(base["E_current_rf40_power_law"], dtype=float).copy()
        m_ratio = np.asarray(
            base["M_over_X_current_rf40_power_law"], dtype=float
        ).copy()
        m_eta_ratio = np.asarray(
            base["M_eta_over_X_current_rf40_power_law"], dtype=float
        ).copy()
        region = np.full(xf.shape, "current_RF40_before_autonomous_modulation", dtype=object)

        positive = xf > 0.0
        log_x = np.full_like(xf, -math.inf)
        log_x[positive] = np.log(xf[positive])
        left, right = self.support_log_interval
        selected = positive & (log_x >= left)
        active = positive & (log_x >= left) & (log_x <= right)

        if np.any(selected):
            modified = self.modulation.profile_values_logX(
                log_x[selected], ef[selected]
            )
            F[selected] = np.asarray(modified["F"], dtype=float)
            U[selected] = np.asarray(modified["U"], dtype=float)
            E[selected] = np.asarray(modified["E"], dtype=float)
            region[selected] = "post_autonomous_modulation_memory"
        if np.any(active):
            region[active] = "selected_autonomous_modulation"

        dm, dm_eta = self._delta_memory_ratio(xf, ef)
        m_ratio += dm
        m_eta_ratio += dm_eta

        changed = selected
        v0 = np.asarray(base["v0_current_rf40_power_law"], dtype=float).copy()
        if np.any(changed):
            axis = self.geometry.physical_profiles.axis_profiles.values(
                np.zeros(np.count_nonzero(changed), dtype=float), ef[changed]
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            v0[changed] = (
                2.0 * ef[changed] * U[changed]
                - 2.0 * self.D * ef[changed] * m_ratio[changed]
                - d * m_eta_ratio[changed]
            ) / L

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current autonomous-modulation profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current autonomous-modulation profile lost positive F/E")

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_autonomous_modulation": F.reshape(shape),
            "U_current_autonomous_modulation": U.reshape(shape),
            "E_current_autonomous_modulation": E.reshape(shape),
            "M_over_X_current_autonomous_modulation": m_ratio.reshape(shape),
            "M_eta_over_X_current_autonomous_modulation": m_eta_ratio.reshape(shape),
            "delta_M_over_X_autonomous_modulation": dm.reshape(shape),
            "delta_M_eta_over_X_autonomous_modulation": dm_eta.reshape(shape),
            "v0_current_autonomous_modulation": v0.reshape(shape),
        }

    def similarity_radial_derivatives(
        self, X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        if np.any(X_arr <= 0.0):
            raise ValueError("similarity_radial_derivatives requires X>0")
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        base = self.current.similarity_radial_derivatives(xf, ef)
        F_X = np.asarray(base["F_current_rf40_power_law_X"], dtype=float).copy()
        U_X = np.asarray(base["U_current_rf40_power_law_X"], dtype=float).copy()
        E_X = np.asarray(base["E_current_rf40_power_law_X"], dtype=float).copy()

        log_x = np.log(xf)
        left, _ = self.support_log_interval
        selected = log_x >= left
        if np.any(selected):
            modified = self.modulation.profile_values_logX(
                log_x[selected], ef[selected]
            )
            inv_x = 1.0 / xf[selected]
            F_X[selected] = np.asarray(modified["D_X_F"], dtype=float) * inv_x
            U_X[selected] = np.asarray(modified["D_X_U"], dtype=float) * inv_x
            E_X[selected] = np.asarray(modified["D_X_E"], dtype=float) * inv_x

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_current_autonomous_modulation_X": F_X.reshape(shape),
            "U_current_autonomous_modulation_X": U_X.reshape(shape),
            "E_current_autonomous_modulation_X": E_X.reshape(shape),
        }

    def values(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        tol = 128.0 * np.finfo(float).eps * self.X_I1_start
        if np.any(X > self.X_I1_start + tol):
            raise ValueError(
                "Cartesian point lies at/after the reserved I1 repair; "
                "current autonomous-modulation candidate fails closed there"
            )
        X_eval = np.minimum(X, self.X_I1_start)
        profile = self.similarity_profile_values(X_eval, coords["eta"])
        q = np.asarray(coords["q"], dtype=float)
        q_radial = profile["v0_current_autonomous_modulation"] / (2.0 * q)
        q_swirl = (
            np.power(q, -self.A - 0.5)
            * profile["F_current_autonomous_modulation"]
        )
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = np.power(q, -self.A) * profile["U_current_autonomous_modulation"]
        return {**coords, **profile, "u": u, "v": v, "w": w}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "current": self.current.configuration(),
            "autonomous_modulation": {
                "frequency": self.frequency,
                "alpha": self.alpha,
                "beta": self.beta,
                "support_left_offset": self.support_left_offset,
                "support_right_offset": self.support_right_offset,
                "quadrature_order": self.quadrature_order,
                "panels_per_period": self.panels_per_period,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianAutonomousModulation":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current autonomous-modulation schema")
        mod = payload.get("autonomous_modulation")
        if not isinstance(mod, Mapping):
            raise ValueError("missing autonomous_modulation configuration")
        return cls(
            current=KokunoPA16CurrentCartesianRF40PowerLaw.from_configuration(
                payload["current"]
            ),
            frequency=int(mod["frequency"]),
            alpha=float(mod["alpha"]),
            beta=float(mod["beta"]),
            support_left_offset=float(mod["support_left_offset"]),
            support_right_offset=float(mod["support_right_offset"]),
            quadrature_order=int(mod["quadrature_order"]),
            panels_per_period=int(mod["panels_per_period"]),
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianAutonomousModulation":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "configuration": self.configuration(),
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "truth_boundary": _TRUTH_BOUNDARY,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        left, right = self.support_log_interval
        eta = 0.2
        mid = 0.5 * (left + right)
        profile = self.similarity_profile_values(math.exp(mid), eta)
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "parent_exact_head": "2c76ebdc41d6c566f43a1305034ba2ff9dce410b",
            "current_RF40_log_X_end": self.current.log_X_4,
            "autonomous_modulation_support_log_X": [left, right],
            "reserved_I1_start_log_X": self.log_X_I1_start,
            "sample_eta": eta,
            "sample_log_X": mid,
            "sample_delta_M_over_X": float(
                np.asarray(profile["delta_M_over_X_autonomous_modulation"])
            ),
            "sample_velocity_magnitude_profile": {
                "F": float(np.asarray(profile["F_current_autonomous_modulation"])),
                "U": float(np.asarray(profile["U_current_autonomous_modulation"])),
                "v0": float(np.asarray(profile["v0_current_autonomous_modulation"])),
            },
            "semantic_sha256": self.semantic_sha256,
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "the selected A/B loop is repository-autonomous, not a recovered source loop",
                "I1 and all later reserved-patch repairs are not applied",
                "evaluation fails closed at the start of I1",
                "no pressure/forcing or held-out NS residual is supplied",
            ],
        }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()

    candidate = KokunoPA16CurrentCartesianAutonomousModulation()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
