"""Bind the current A1 pre-I1 candidate to the existing PA.17 I1 repair family.

Pinned provenance is the public KokunoYumeto corrected 2026-09-09 reconstruction
(`yang-mills-interacting-workbench@143f6773...`).  The public reconstruction
inserts a finite-frequency radial modulation and restores its five induced
prefix-moment discrepancies on reserved interval I1 with the later PA.17 map.

This module does not recover the unpublished admissible source loop.  It reuses
the repository-autonomous modulation already bound to the current A1 lineage by
``KokunoPA16CurrentCartesianAutonomousModulation`` and the already-reviewed
eta-smooth PA.17 coefficient family from ``KokunoEtaSmoothShearRepairFamily``.

The executable seam added here is current-lineage composition.  Before I1 the
candidate is replayed exactly.  On I1 the current unmodulated RF40 base is
combined with the PA.17 profile correction while the *actual* modulation
prefix-M memory arriving from the current pre-I1 candidate is propagated as a
physical constant (therefore as 1/X in M/X).  The resulting total M and M_eta
are used to rebuild v0 and hence the Cartesian velocity.

This is a provenance-labelled candidate realization.  It is not paper-exact,
not an independent PDE validation, and not a residual-improvement claim.
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

from .kokuno_eta_smooth_shear_repair_family import (
    KokunoEtaSmoothShearRepairFamily,
)
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa16_current_cartesian_autonomous_modulation import (
    KokunoPA16CurrentCartesianAutonomousModulation,
)

SCHEMA = "kokuno-pa16-current-cartesian-i1-repair-v1"
PARENT_EXACT_HEAD = "ddf71bab2a2038bab2fba5b854ad713c6b6c7e95"

_SOURCE_FORMULAS = {
    "finite_N_insertion": "E_N=E exp(A/N); U_N=U+B/N; phi=N log X mod 1",
    "repair_location": "five modulation-induced prefix-moment discrepancies are restored on reserved I1",
    "PA17_background": "U=0; E=K(eta) X^(-1/2-lambda)",
    "PA17_rows": "row order (M,J,I,S,C_p)",
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
}
_NUMERICAL_REALIZATION = {
    "incoming_modulation": "consume exact current #1035 autonomous pre-I1 candidate",
    "repair_family": (
        "reuse repository-autonomous #399 eta-smooth Chebyshev PA.17 family "
        "on the exact current modulation dependency"
    ),
    "memory_handoff": (
        "measure actual Delta(M/X),Delta(M_eta/X) at I1 start from #1035; "
        "carry physical Delta M,Delta M_eta as constants through I1"
    ),
    "domain_stop": "fail closed at the upper endpoint of reserved I1",
    "new_residual_tuning": "none",
}
_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_pre_I1_modulation_consumed": True,
    "eta_smooth_PA17_family_consumed": True,
    "I1_repair_applied_to_current_modulation": True,
    "current_modulation_prefix_M_memory_carried_through_I1": True,
    "cartesian_velocity_executable_through_I1": True,
    "velocity_interface_vectorized": True,
    "configuration_serializable": True,
    "child_semantic_identity_binds_parent_and_modulation_dependency": True,
    "parent_pre_I1_semantic_dependency_binding_repaired_by_this_increment": False,
    "source_hidden_loop_parameters_recovered": False,
    "source_admissible_loop_reconstructed": False,
    "actual_source_admissible_loop_discrepancy_supplied": False,
    "I1_cone_repair_applied_to_actual_source_modulation": False,
    "cone_modulation_completed": False,
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
class KokunoPA16CurrentCartesianI1Repair:
    """Current A1 Cartesian leading candidate with the selected PA.17 I1 repair."""

    pre_i1: KokunoPA16CurrentCartesianAutonomousModulation = field(
        default_factory=KokunoPA16CurrentCartesianAutonomousModulation,
        repr=False,
        compare=False,
    )
    interpolation_nodes: int = 11
    closure_tolerance: float = 5.0e-7

    _i1_family: KokunoEtaSmoothShearRepairFamily = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(
            self.pre_i1, KokunoPA16CurrentCartesianAutonomousModulation
        ):
            raise TypeError(
                "pre_i1 must be KokunoPA16CurrentCartesianAutonomousModulation"
            )
        family = KokunoEtaSmoothShearRepairFamily(
            modulation=self.pre_i1.modulation,
            interpolation_nodes=self.interpolation_nodes,
            closure_tolerance=self.closure_tolerance,
        )
        i1_low, i1_high = self.pre_i1.modulation.outer_schedule.reserved_log_intervals()[
            "I1"
        ]
        if not math.isclose(
            float(i1_low),
            self.pre_i1.log_X_I1_start,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("pre-I1 endpoint and reserved I1 start disagree")
        repair = family.repair
        support_low = repair.log_X_1 + math.log(0.80)
        support_high = repair.log_X_1 + math.log(1.70)
        if not (float(i1_low) < support_low < support_high < float(i1_high)):
            raise ValueError("PA.17 repair support is not strictly inside reserved I1")
        if family.modulation.sha256 != self.pre_i1.modulation.sha256:
            raise ValueError("I1 family detached from the exact current modulation identity")
        object.__setattr__(self, "interpolation_nodes", int(family.interpolation_nodes))
        object.__setattr__(self, "closure_tolerance", float(family.closure_tolerance))
        object.__setattr__(self, "_i1_family", family)

    @property
    def i1_family(self) -> KokunoEtaSmoothShearRepairFamily:
        return self._i1_family

    @property
    def current(self):
        return self.pre_i1.current

    @property
    def geometry(self):
        return self.pre_i1.geometry

    @property
    def A(self) -> float:
        return float(self.pre_i1.A)

    @property
    def D(self) -> float:
        return float(self.pre_i1.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.pre_i1.eta_interval)

    @property
    def log_X_I1_start(self) -> float:
        return float(self.pre_i1.log_X_I1_start)

    @property
    def log_X_I1_end(self) -> float:
        return float(
            self.pre_i1.modulation.outer_schedule.reserved_log_intervals()["I1"][1]
        )

    @property
    def X_I1_start(self) -> float:
        return float(math.exp(self.log_X_I1_start))

    @property
    def X_I1_end(self) -> float:
        return float(math.exp(self.log_X_I1_end))

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.pre_i1.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_I1_end)):
            raise ValueError(
                "X must lie in the current I1-repair domain "
                f"[0,{self.X_I1_end:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _incoming_modulation_memory_ratio(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        eta_arr = np.asarray(eta, dtype=float)
        X0 = np.full(eta_arr.shape, self.X_I1_start, dtype=float)
        incoming = self.pre_i1.similarity_profile_values(X0, eta_arr)
        return (
            np.asarray(
                incoming["delta_M_over_X_autonomous_modulation"], dtype=float
            ),
            np.asarray(
                incoming["delta_M_eta_over_X_autonomous_modulation"], dtype=float
            ),
        )

    def _carried_modulation_memory_ratio(
        self, X: np.ndarray, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        incoming, incoming_eta = self._incoming_modulation_memory_ratio(eta)
        factor = self.X_I1_start / np.asarray(X, dtype=float)
        return factor * incoming, factor * incoming_eta

    def similarity_profile_values(
        self, X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)

        F = np.empty_like(xf)
        U = np.empty_like(xf)
        E = np.empty_like(xf)
        m_ratio = np.empty_like(xf)
        m_eta_ratio = np.empty_like(xf)
        mod_m_ratio = np.zeros_like(xf)
        mod_m_eta_ratio = np.zeros_like(xf)
        repair_m_ratio = np.zeros_like(xf)
        repair_m_eta_ratio = np.zeros_like(xf)
        v0 = np.empty_like(xf)
        region = np.empty(xf.shape, dtype=object)

        inherited = xf <= self.X_I1_start
        if np.any(inherited):
            p = self.pre_i1.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(
                p["F_current_autonomous_modulation"], dtype=float
            )
            U[inherited] = np.asarray(
                p["U_current_autonomous_modulation"], dtype=float
            )
            E[inherited] = np.asarray(
                p["E_current_autonomous_modulation"], dtype=float
            )
            m_ratio[inherited] = np.asarray(
                p["M_over_X_current_autonomous_modulation"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_autonomous_modulation"], dtype=float
            )
            mod_m_ratio[inherited] = np.asarray(
                p["delta_M_over_X_autonomous_modulation"], dtype=float
            )
            mod_m_eta_ratio[inherited] = np.asarray(
                p["delta_M_eta_over_X_autonomous_modulation"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_autonomous_modulation"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        on_i1 = ~inherited
        if np.any(on_i1):
            Xs = xf[on_i1]
            etas = ef[on_i1]
            base = self.current.similarity_profile_values(Xs, etas)
            repair = self.i1_family.profile_correction_logX(np.log(Xs), etas)
            dm, dm_eta = self._carried_modulation_memory_ratio(Xs, etas)

            F[on_i1] = np.asarray(
                base["F_current_rf40_power_law"], dtype=float
            ) + np.asarray(repair["delta_F"], dtype=float)
            U[on_i1] = np.asarray(
                base["U_current_rf40_power_law"], dtype=float
            ) + np.asarray(repair["delta_U"], dtype=float)
            E[on_i1] = np.asarray(
                base["E_current_rf40_power_law"], dtype=float
            ) + np.asarray(repair["delta_E"], dtype=float)

            rep_m = np.asarray(repair["delta_M"], dtype=float) / Xs
            rep_m_eta = np.asarray(repair["delta_M_eta"], dtype=float) / Xs
            m_ratio[on_i1] = (
                np.asarray(base["M_over_X_current_rf40_power_law"], dtype=float)
                + dm
                + rep_m
            )
            m_eta_ratio[on_i1] = (
                np.asarray(
                    base["M_eta_over_X_current_rf40_power_law"], dtype=float
                )
                + dm_eta
                + rep_m_eta
            )
            mod_m_ratio[on_i1] = dm
            mod_m_eta_ratio[on_i1] = dm_eta
            repair_m_ratio[on_i1] = rep_m
            repair_m_eta_ratio[on_i1] = rep_m_eta

            axis = self.geometry.physical_profiles.axis_profiles.values(
                np.zeros_like(etas), etas
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            v0[on_i1] = (
                2.0 * etas * U[on_i1]
                - 2.0 * self.D * etas * m_ratio[on_i1]
                - d * m_eta_ratio[on_i1]
            ) / L
            region[on_i1] = "current_eta_smooth_PA17_I1_repair"

        arrays = (
            F,
            U,
            E,
            m_ratio,
            m_eta_ratio,
            mod_m_ratio,
            mod_m_eta_ratio,
            repair_m_ratio,
            repair_m_eta_ratio,
            v0,
        )
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current I1-repair profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current I1-repair profile lost positive F/E")

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_i1_repair": F.reshape(shape),
            "U_current_i1_repair": U.reshape(shape),
            "E_current_i1_repair": E.reshape(shape),
            "M_over_X_current_i1_repair": m_ratio.reshape(shape),
            "M_eta_over_X_current_i1_repair": m_eta_ratio.reshape(shape),
            "delta_M_over_X_autonomous_modulation": mod_m_ratio.reshape(shape),
            "delta_M_eta_over_X_autonomous_modulation": mod_m_eta_ratio.reshape(shape),
            "delta_M_over_X_I1_repair": repair_m_ratio.reshape(shape),
            "delta_M_eta_over_X_I1_repair": repair_m_eta_ratio.reshape(shape),
            "v0_current_i1_repair": v0.reshape(shape),
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
        F_X = np.empty_like(xf)
        U_X = np.empty_like(xf)
        E_X = np.empty_like(xf)

        inherited = xf <= self.X_I1_start
        if np.any(inherited):
            d0 = self.pre_i1.similarity_radial_derivatives(
                xf[inherited], ef[inherited]
            )
            F_X[inherited] = np.asarray(
                d0["F_current_autonomous_modulation_X"], dtype=float
            )
            U_X[inherited] = np.asarray(
                d0["U_current_autonomous_modulation_X"], dtype=float
            )
            E_X[inherited] = np.asarray(
                d0["E_current_autonomous_modulation_X"], dtype=float
            )

        on_i1 = ~inherited
        if np.any(on_i1):
            Xs = xf[on_i1]
            etas = ef[on_i1]
            base = self.current.similarity_radial_derivatives(Xs, etas)
            repair = self.i1_family.profile_correction_logX(np.log(Xs), etas)
            F_X[on_i1] = np.asarray(
                base["F_current_rf40_power_law_X"], dtype=float
            ) + np.asarray(repair["delta_F_X"], dtype=float)
            U_X[on_i1] = np.asarray(
                base["U_current_rf40_power_law_X"], dtype=float
            ) + np.asarray(repair["delta_U_X"], dtype=float)
            E_X[on_i1] = np.asarray(
                base["E_current_rf40_power_law_X"], dtype=float
            ) + np.asarray(repair["delta_E_X"], dtype=float)

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_current_i1_repair_X": F_X.reshape(shape),
            "U_current_i1_repair_X": U_X.reshape(shape),
            "E_current_i1_repair_X": E_X.reshape(shape),
        }

    def values(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        tol = 128.0 * np.finfo(float).eps * self.X_I1_end
        if np.any(X > self.X_I1_end + tol):
            raise ValueError(
                "Cartesian point lies after reserved I1; "
                "current I1-repair candidate fails closed there"
            )
        X_eval = np.minimum(X, self.X_I1_end)
        profile = self.similarity_profile_values(X_eval, coords["eta"])
        q = np.asarray(coords["q"], dtype=float)
        q_radial = profile["v0_current_i1_repair"] / (2.0 * q)
        q_swirl = np.power(q, -self.A - 0.5) * profile["F_current_i1_repair"]
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = np.power(q, -self.A) * profile["U_current_i1_repair"]
        return {**coords, **profile, "u": u, "v": v, "w": w}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def i1_exit_report(self, eta: Any) -> dict[str, Any]:
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        X = np.full(eta_arr.shape, self.X_I1_end, dtype=float)
        total = self.similarity_profile_values(X, eta_arr)
        base = self.current.similarity_profile_values(X, eta_arr)
        outgoing_m = np.asarray(
            total["M_over_X_current_i1_repair"], dtype=float
        ) - np.asarray(base["M_over_X_current_rf40_power_law"], dtype=float)
        outgoing_m_eta = np.asarray(
            total["M_eta_over_X_current_i1_repair"], dtype=float
        ) - np.asarray(
            base["M_eta_over_X_current_rf40_power_law"], dtype=float
        )
        incoming_m, incoming_m_eta = self._incoming_modulation_memory_ratio(eta_arr)
        return {
            "eta": np.asarray(eta_arr, dtype=float).tolist(),
            "family_closure": self.i1_family.closure_report(eta_arr),
            "incoming_autonomous_M_over_X_at_I1_start": np.asarray(
                incoming_m, dtype=float
            ).tolist(),
            "incoming_autonomous_M_eta_over_X_at_I1_start": np.asarray(
                incoming_m_eta, dtype=float
            ).tolist(),
            "outgoing_overlay_M_over_X_at_I1_exit": outgoing_m.tolist(),
            "outgoing_overlay_M_eta_over_X_at_I1_exit": outgoing_m_eta.tolist(),
            "max_abs_coefficient": self.i1_family.family_report()[
                "max_abs_coefficient"
            ],
            "coefficient_limit": self.i1_family.family_report()[
                "coefficient_limit"
            ],
            "independent_PDE_validation": False,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "pre_i1": self.pre_i1.configuration(),
            "i1_repair": {
                "interpolation_nodes": self.interpolation_nodes,
                "closure_tolerance": self.closure_tolerance,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianI1Repair":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current I1-repair schema")
        repair = payload.get("i1_repair")
        if not isinstance(repair, Mapping):
            raise ValueError("missing i1_repair configuration")
        return cls(
            pre_i1=KokunoPA16CurrentCartesianAutonomousModulation.from_configuration(
                payload["pre_i1"]
            ),
            interpolation_nodes=int(repair["interpolation_nodes"]),
            closure_tolerance=float(repair["closure_tolerance"]),
        )

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
    ) -> "KokunoPA16CurrentCartesianI1Repair":
        return cls.from_configuration(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "configuration": self.configuration(),
            "dependencies": {
                "pre_i1_parent_exact_head": PARENT_EXACT_HEAD,
                "pre_i1_semantic_sha256": self.pre_i1.semantic_sha256,
                "modulation_dependency_sha256": self.pre_i1.modulation.sha256,
                "eta_smooth_i1_family_sha256": self.i1_family.sha256,
            },
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
        return hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta = np.asarray([0.0, 0.25, 0.6, 1.0], dtype=float)
        exit_report = self.i1_exit_report(eta)
        active_log_X = self.i1_family.repair.log_X_1
        active = self.similarity_profile_values(math.exp(active_log_X), 0.25)
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "parent_exact_head": PARENT_EXACT_HEAD,
            "reserved_I1_log_X": [self.log_X_I1_start, self.log_X_I1_end],
            "repair_center_log_X": active_log_X,
            "family_report": self.i1_family.family_report(),
            "exit_report": exit_report,
            "sample_active_profile": {
                "F": float(np.asarray(active["F_current_i1_repair"])),
                "U": float(np.asarray(active["U_current_i1_repair"])),
                "v0": float(np.asarray(active["v0_current_i1_repair"])),
            },
            "semantic_dependencies": {
                "pre_i1_semantic_sha256": self.pre_i1.semantic_sha256,
                "modulation_dependency_sha256": self.pre_i1.modulation.sha256,
                "eta_smooth_i1_family_sha256": self.i1_family.sha256,
            },
            "semantic_sha256": self.semantic_sha256,
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "the selected finite-N A/B loop remains repository-autonomous",
                "the PA.17 interpolation family is repository-autonomous numerics",
                "the standalone #1035 parent semantic-dependency finding is not relabelled fixed",
                "I2/I3/I4 and terminal/global outer completion are not applied",
                "no matched pressure, restricted forcing, or held-out complete NS residual is supplied",
            ],
        }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()
    candidate = KokunoPA16CurrentCartesianI1Repair()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
