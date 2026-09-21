"""Current Kokuno Cartesian leading field through the first post-X_R RF40 turn.

Pinned provenance is the corrected 2026-09-09 KokunoYumeto public
reconstruction at commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.

Agent 1 already materializes the current candidate-side Cartesian leading field
through ``X_R``.  The corrected outer schedule then resets the logarithmic
coordinate ``y=log(X/X_R)`` and, for ``0<=y<=1``, prescribes

    U = 4 eta,
    ell = d_y log(sqrt(2X) E) = (3/5) (1-sigma(y)),

hence

    log(E/f) = log(P_*) + y/10 - (3/5) int_0^y sigma(s) ds,
    F = E/sqrt(2X),
    f(eta) = (1+eta^2)^(-1).

The public ``KokunoOuterBaseSchedule`` already implements this RF40 first-turn
formula and its radial derivatives.  This adapter does not duplicate that
schedule.  It binds the *current* A1 ``X_R`` handoff to it and preserves the
actual current incompressibility primitive instead of silently resetting the
floating-point PA.16 history to the ideal value.

Because ``U=4 eta`` throughout this first turn, the carried primitive is exact:

    (M/X)(X) = 4 eta + (X_R/X) [(M/X)_R - 4 eta],
    (M_eta/X)(X) = 4 + (X_R/X) [(M_eta/X)_R - 4].

The Cartesian radial profile remains

    v0 = (2 eta U - 2 D eta M/X - d M_eta/X) / L.

This is one deliberately small executable increment.  It stops at
``X_1=e X_R`` and fails closed beyond that point.  The later RF40 axial
shutdown, lambda turn, power-law interval, cone modulation, I1--I4 overlays,
globally matched pressure/forcing, complete held-out NS residual, paper exact
field and OpenAI field identity remain unmaterialized.
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

from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa16_current_cartesian_exterior_to_xr import (
    KokunoPA16CurrentCartesianExteriorToXR,
)


SCHEMA = "kokuno-pa16-current-cartesian-rf40-first-turn-v1"

_SOURCE_FORMULAS = {
    "coordinates": "y=log(X/X_R), reset to y=0 at X=X_R",
    "first_turn_domain": "0<=y<=1, equivalently X_R<=X<=e X_R",
    "first_turn_U": "U=4 eta",
    "first_turn_ell": "ell=d_y log(sqrt(2X)E)=(3/5)(1-sigma(y))",
    "first_turn_E": (
        "log(E/f)=log(P_*)+y/10-(3/5) int_0^y sigma(s) ds; "
        "f=(1+eta^2)^(-1)"
    ),
    "F": "F=E/sqrt(2X)",
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "XR_seam": "consume the exact current #980 candidate-side X_R handoff",
    "RF40_profile": "reuse KokunoOuterBaseSchedule; do not duplicate the public schedule",
    "primitive_memory": (
        "carry current floating-point M/X and d_eta(M/X) analytically under U=4 eta; "
        "do not reset the current candidate to ideal moments at X_R"
    ),
    "coordinates": "reuse the governed current A1 q/X/eta inverse and Cartesian map",
    "new_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_cartesian_leading_through_XR_consumed": True,
    "public_RF40_first_turn_formula_reused": True,
    "current_lineage_RF40_first_turn_materialized": True,
    "current_cartesian_spacetime_leading_velocity_through_RF40_first_turn_materialized": True,
    "current_incompressibility_memory_carried_through_RF40_first_turn": True,
    "axis_regular_cartesian_formula_reused": True,
    "velocity_interface_vectorized": True,
    "velocity_configuration_serializable": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_exact_current_lineage_certified": False,
    "full_post_XR_RF40_current_lineage_materialized": False,
    "RF40_axial_shutdown_current_lineage_materialized": False,
    "RF40_lambda_turn_current_lineage_materialized": False,
    "RF40_power_law_current_lineage_materialized": False,
    "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_prepared_appendixA_pressure_stress_materialized": False,
    "source_admitted_kappa0_materialized": False,
    "source_global_inner_to_outer_join_admitted": False,
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
class KokunoPA16CurrentCartesianRF40FirstTurn:
    """Extend the current candidate-side Cartesian leading field one RF40 stage."""

    current: KokunoPA16CurrentCartesianExteriorToXR = field(
        default_factory=KokunoPA16CurrentCartesianExteriorToXR,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.current, KokunoPA16CurrentCartesianExteriorToXR):
            raise TypeError("current must be KokunoPA16CurrentCartesianExteriorToXR")
        if not math.isclose(
            math.log(self.X_R), self.outer_base.log_X_R, rel_tol=0.0, abs_tol=2.0e-12
        ):
            raise ValueError("current X_R and RF40 schedule X_R disagree")
        if not math.isclose(
            self.log_X_1, self.outer_base.log_stage1_end, rel_tol=0.0, abs_tol=2.0e-12
        ):
            raise ValueError("RF40 first-turn endpoint disagrees with public outer schedule")
        if not (0.0 < self.X_R < self.X_1 < np.finfo(float).max):
            raise ValueError("materialized current RF40 first-turn interval is invalid")

    @property
    def outer_base(self):
        return self.current.outer_base

    @property
    def X_R(self) -> float:
        return float(self.current.X_R)

    @property
    def log_X_1(self) -> float:
        return float(self.outer_base.log_stage1_end)

    @property
    def X_1(self) -> float:
        return float(math.exp(self.log_X_1))

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.current.eta_interval)

    @property
    def A(self) -> float:
        return float(self.current.A)

    @property
    def D(self) -> float:
        return float(self.current.D)

    @property
    def h(self) -> float:
        return self.A - 0.5

    def similarity_coordinates(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        return self.current.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_1)):
            raise ValueError(
                f"X must lie in the current RF40-first-turn domain [0,{self.X_1:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _XR_memory(self, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return current ``(M/X)_R`` and ``d_eta(M/X)_R`` at the exact X_R seam."""
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        seam = self.current.similarity_profile_values(
            np.full(unique.shape, self.X_R, dtype=float), unique
        )
        m = np.asarray(seam["M_over_X_current_exterior"], dtype=float)
        m_eta = np.asarray(seam["M_eta_over_X_current_exterior"], dtype=float)
        return m[inverse].reshape(values.shape), m_eta[inverse].reshape(values.shape)

    def similarity_profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return unified ``F,U,E,M/X,M_eta/X,v0`` through the first RF40 turn."""
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)

        F = np.empty_like(xf)
        U = np.empty_like(xf)
        E = np.empty_like(xf)
        m_ratio = np.empty_like(xf)
        m_eta_ratio = np.empty_like(xf)
        v0 = np.empty_like(xf)
        region = np.empty(xf.shape, dtype=object)

        inherited = xf <= self.X_R
        if np.any(inherited):
            vals = self.current.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(vals["F_current_exterior"], dtype=float)
            U[inherited] = np.asarray(vals["U_current_exterior"], dtype=float)
            E[inherited] = np.asarray(vals["E_current_exterior"], dtype=float)
            m_ratio[inherited] = np.asarray(
                vals["M_over_X_current_exterior"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                vals["M_eta_over_X_current_exterior"], dtype=float
            )
            v0[inherited] = np.asarray(vals["v0_current_exterior"], dtype=float)
            region[inherited] = "current_cartesian_leading_through_XR"

        first_turn = ~inherited
        if np.any(first_turn):
            Xo = xf[first_turn]
            etao = ef[first_turn]
            base = self.outer_base.profile_values(Xo, etao)
            stages = np.asarray(base["stage"], dtype=object)
            if np.any(stages != "first_l_turn"):
                raise RuntimeError("X_R<X<=e X_R must remain on the public RF40 first_l_turn")

            seam_m, seam_m_eta = self._XR_memory(etao)
            ratio = self.X_R / Xo
            mo = 4.0 * etao + ratio * (seam_m - 4.0 * etao)
            meo = 4.0 + ratio * (seam_m_eta - 4.0)

            axis = self.current.current.geometry.physical_profiles.axis_profiles.values(
                np.zeros_like(etao), etao
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            Uo = np.asarray(base["U"], dtype=float)
            v0o = (2.0 * etao * Uo - 2.0 * self.D * etao * mo - d * meo) / L

            F[first_turn] = np.asarray(base["F"], dtype=float)
            U[first_turn] = Uo
            E[first_turn] = np.asarray(base["E"], dtype=float)
            m_ratio[first_turn] = mo
            m_eta_ratio[first_turn] = meo
            v0[first_turn] = v0o
            region[first_turn] = "public_RF40_first_l_turn"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current RF40 first-turn profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current RF40 first-turn profile lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_rf40_first_turn": F.reshape(shape),
            "U_current_rf40_first_turn": U.reshape(shape),
            "E_current_rf40_first_turn": E.reshape(shape),
            "M_over_X_current_rf40_first_turn": m_ratio.reshape(shape),
            "M_eta_over_X_current_rf40_first_turn": m_eta_ratio.reshape(shape),
            "v0_current_rf40_first_turn": v0.reshape(shape),
        }

    def similarity_radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return first physical-X derivatives of ``F,U,E`` through ``e X_R``."""
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        if np.any(X_arr <= 0.0):
            raise ValueError("similarity_radial_derivatives requires X>0")
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        F_X = np.empty_like(xf)
        U_X = np.empty_like(xf)
        E_X = np.empty_like(xf)

        inherited = xf <= self.X_R
        if np.any(inherited):
            vals = self.current.similarity_radial_derivatives(xf[inherited], ef[inherited])
            F_X[inherited] = np.asarray(vals["F_current_exterior_X"], dtype=float)
            U_X[inherited] = np.asarray(vals["U_current_exterior_X"], dtype=float)
            E_X[inherited] = np.asarray(vals["E_current_exterior_X"], dtype=float)

        first_turn = ~inherited
        if np.any(first_turn):
            vals = self.outer_base.profile_values(xf[first_turn], ef[first_turn])
            stages = np.asarray(vals["stage"], dtype=object)
            if np.any(stages != "first_l_turn"):
                raise RuntimeError("RF40 derivative evaluation escaped the first_l_turn")
            F_X[first_turn] = np.asarray(vals["F_X"], dtype=float)
            U_X[first_turn] = np.asarray(vals["U_X"], dtype=float)
            E_X[first_turn] = np.asarray(vals["E_X"], dtype=float)

        if any(np.any(~np.isfinite(arr)) for arr in (F_X, U_X, E_X)):
            raise RuntimeError("current RF40 first-turn radial derivatives are non-finite")
        return {
            "F_current_rf40_first_turn_X": F_X.reshape(shape),
            "U_current_rf40_first_turn_X": U_X.reshape(shape),
            "E_current_rf40_first_turn_X": E_X.reshape(shape),
        }

    def values(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Return source coordinates, current profiles and Cartesian velocity through ``e X_R``."""
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        radial_tol = 128.0 * np.finfo(float).eps * self.X_1
        if np.any(X > self.X_1 + radial_tol):
            raise ValueError(
                "point lies beyond current RF40 first-turn endpoint; later RF40 is not materialized"
            )
        X_eval = np.minimum(X, self.X_1)
        profile = self.similarity_profile_values(X_eval, coords["eta"])
        q = np.asarray(coords["q"], dtype=float)
        q_radial = profile["v0_current_rf40_first_turn"] / (2.0 * q)
        q_swirl = np.power(q, -self.A - 0.5) * profile["F_current_rf40_first_turn"]
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = np.power(q, -self.A) * profile["U_current_rf40_first_turn"]
        return {**coords, **profile, "u": u, "v": v, "w": w}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Unified vectorized ``velocity(x,y,z,t)->[...,3]`` through RF40 first turn."""
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {"schema": SCHEMA, "current": self.current.configuration()}

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianRF40FirstTurn":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current Cartesian RF40 first-turn schema")
        return cls(
            current=KokunoPA16CurrentCartesianExteriorToXR.from_configuration(
                payload["current"]
            )
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA16CurrentCartesianRF40FirstTurn":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "truth_boundary": _TRUTH_BOUNDARY,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta = 0.25
        at_R = self.current.similarity_profile_values(self.X_R, eta)
        right_R = self.outer_base.profile_values(self.X_R, eta)
        at_R_jet = self.current.similarity_radial_derivatives(self.X_R, eta)
        right_R_jet = self.outer_base.profile_values(self.X_R, eta)
        at_1 = self.similarity_profile_values(self.X_1, eta)
        at_1_jet = self.similarity_radial_derivatives(self.X_1, eta)

        def scalar(value: Any) -> float:
            return float(np.asarray(value))

        value_jump = {
            "F": scalar(right_R["F"]) - scalar(at_R["F_current_exterior"]),
            "U": scalar(right_R["U"]) - scalar(at_R["U_current_exterior"]),
            "E": scalar(right_R["E"]) - scalar(at_R["E_current_exterior"]),
        }
        jet_jump = {
            "F_X": scalar(right_R_jet["F_X"]) - scalar(at_R_jet["F_current_exterior_X"]),
            "U_X": scalar(right_R_jet["U_X"]) - scalar(at_R_jet["U_current_exterior_X"]),
            "E_X": scalar(right_R_jet["E_X"]) - scalar(at_R_jet["E_current_exterior_X"]),
        }
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "numerical_realization": copy.deepcopy(_NUMERICAL_REALIZATION),
            "domain": {
                "X_R": self.X_R,
                "X_1": self.X_1,
                "log_X_R": math.log(self.X_R),
                "log_X_1": self.log_X_1,
                "eta_interval": list(self.eta_interval),
                "time_interval": list(self.current.current.geometry.time_interval),
            },
            "XR_seam_probe_eta": eta,
            "XR_value_jump_public_right_minus_current_left": value_jump,
            "XR_first_radial_jet_jump_public_right_minus_current_left": jet_jump,
            "X1_probe": {
                "F": scalar(at_1["F_current_rf40_first_turn"]),
                "U": scalar(at_1["U_current_rf40_first_turn"]),
                "E": scalar(at_1["E_current_rf40_first_turn"]),
                "M_over_X_minus_ideal": scalar(at_1["M_over_X_current_rf40_first_turn"]) - 4.0 * eta,
                "M_eta_over_X_minus_ideal": scalar(at_1["M_eta_over_X_current_rf40_first_turn"]) - 4.0,
                "v0": scalar(at_1["v0_current_rf40_first_turn"]),
                "X_F_X_over_F": self.X_1
                * scalar(at_1_jet["F_current_rf40_first_turn_X"])
                / scalar(at_1["F_current_rf40_first_turn"]),
                "X_E_X_over_E": self.X_1
                * scalar(at_1_jet["E_current_rf40_first_turn_X"])
                / scalar(at_1["E_current_rf40_first_turn"]),
                "X_U_X": self.X_1 * scalar(at_1_jet["U_current_rf40_first_turn_X"]),
            },
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="write deterministic scoped receipt JSON")
    parser.add_argument("--config-output", required=True, help="write deterministic configuration JSON")
    args = parser.parse_args()

    candidate = KokunoPA16CurrentCartesianRF40FirstTurn()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    _main()
