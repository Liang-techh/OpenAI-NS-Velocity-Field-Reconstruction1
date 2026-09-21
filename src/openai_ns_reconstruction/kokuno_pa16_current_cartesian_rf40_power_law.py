"""Current Kokuno Cartesian leading field through the RF40c power-law stage.

Pinned provenance is the corrected 2026-09-09 KokunoYumeto public
reconstruction at commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.

The parent Agent-1 candidate reaches the RF40 lambda-turn endpoint

    X_3 = X_w = X_R exp(T_d+2).

The corrected RF40c base schedule then keeps, for logarithmic length ``T_w``,

    U = 0,
    ell = d_logX log(sqrt(2X) E) = -lambda,

on

    X_w < X <= X_4 = X_w exp(T_w).

The repository already implements this public branch in
``KokunoOuterBaseSchedule``.  This adapter reuses it and carries the actual
current-lineage incompressibility primitive.  Since U=0, physical M and M_eta
remain constant, so

    M/X = (X_w/X) (M/X)_w,
    M_eta/X = (X_w/X) (M_eta/X)_w.

No hidden source parameter is inferred and no residual knob is introduced.
This minimal increment completes only the current RF40 *base backbone* through
its power-law endpoint.  Cone modulation, I1--I4 overlays, terminal-tail/global
assembly, matched pressure/forcing and held-out PDE validation remain separate.
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
from .kokuno_pa16_current_cartesian_rf40_lambda_turn import (
    KokunoPA16CurrentCartesianRF40LambdaTurn,
)


SCHEMA = "kokuno-pa16-current-cartesian-rf40-power-law-v1"

_SOURCE_FORMULAS = {
    "coordinates": "RF40c uses log(X/X_w) after X_w=X_R exp(T_d+2)",
    "power_law_domain": "X_w<X<=X_w exp(T_w); T_w=60 log(1/lambda)",
    "power_law_U": "U=0",
    "power_law_ell": "ell=d_logX log(sqrt(2X)E)=-lambda",
    "power_law_E": "E=e_w f(eta) exp[(-1/2-lambda) log(X/X_w)]",
    "F": "F=E/sqrt(2X)",
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "current_memory_propagation": (
        "because U=0, M and M_eta are constant: "
        "M/X=(X_w/X)(M/X)_w and M_eta/X=(X_w/X)(M_eta/X)_w"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "Xw_seam": "consume the exact current Agent-1 lambda-turn X_w handoff",
    "RF40_profile": (
        "reuse KokunoOuterBaseSchedule power_law values/derivatives; "
        "do not duplicate the public RF40c schedule"
    ),
    "primitive_memory": (
        "carry current floating-point M/X and d_eta(M/X) exactly under U=0; "
        "do not reset to the public base primitive"
    ),
    "coordinates": "reuse the governed current A1 q/X/eta inverse and Cartesian map",
    "new_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_cartesian_leading_through_RF40_lambda_turn_consumed": True,
    "public_RF40_power_law_formula_reused": True,
    "current_lineage_RF40_power_law_materialized": True,
    "current_RF40_base_backbone_through_power_law_endpoint_materialized": True,
    "current_cartesian_spacetime_leading_velocity_through_RF40_power_law_materialized": True,
    "current_incompressibility_memory_carried_through_RF40_power_law": True,
    "axis_regular_cartesian_formula_reused": True,
    "velocity_interface_vectorized": True,
    "velocity_configuration_serializable": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_exact_current_lineage_certified": False,
    "full_post_XR_RF40_current_lineage_materialized": False,
    "current_lineage_cone_I1_I2_I3_I4_outer_overlays_completed": False,
    "source_terminal_tail_schedule_bound": False,
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
class KokunoPA16CurrentCartesianRF40PowerLaw:
    """Extend the current candidate through the RF40c constant-lambda base stage."""

    current: KokunoPA16CurrentCartesianRF40LambdaTurn = field(
        default_factory=KokunoPA16CurrentCartesianRF40LambdaTurn,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.current, KokunoPA16CurrentCartesianRF40LambdaTurn):
            raise TypeError("current must be KokunoPA16CurrentCartesianRF40LambdaTurn")
        if not math.isclose(
            self.log_X_3,
            self.outer_base.log_X_w,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("current X_3 and public RF40 X_w endpoint disagree")
        if not math.isclose(
            self.log_X_4,
            self.outer_base.log_X_end,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("RF40 power-law endpoint disagrees with public schedule")
        if not math.isclose(
            self.log_X_4 - self.log_X_3,
            self.T_w,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("RF40 power-law logarithmic length must equal T_w")
        if not (0.0 < self.X_3 < self.X_4 < np.finfo(float).max):
            raise ValueError("materialized current RF40 power-law interval is invalid")

    @property
    def outer_base(self):
        return self.current.outer_base

    @property
    def X_R(self) -> float:
        return float(self.current.X_R)

    @property
    def X_3(self) -> float:
        return float(self.current.X_3)

    @property
    def log_X_3(self) -> float:
        return float(self.current.log_X_3)

    @property
    def log_X_4(self) -> float:
        return float(self.outer_base.log_X_end)

    @property
    def X_4(self) -> float:
        return float(math.exp(self.log_X_4))

    @property
    def T_w(self) -> float:
        return float(self.outer_base.outer_schedule.T_w)

    @property
    def lambda_outer(self) -> float:
        return float(self.outer_base.outer_schedule.lambda_outer)

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

    @property
    def geometry(self):
        return self.current.geometry

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.current.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_4)):
            raise ValueError(
                "X must lie in the current RF40-power-law domain "
                f"[0,{self.X_4:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _X3_memory(self, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        seam = self.current.similarity_profile_values(
            np.full(unique.shape, self.X_3, dtype=float), unique
        )
        m = np.asarray(seam["M_over_X_current_rf40_lambda_turn"], dtype=float)
        m_eta = np.asarray(
            seam["M_eta_over_X_current_rf40_lambda_turn"], dtype=float
        )
        return (
            m[inverse].reshape(values.shape),
            m_eta[inverse].reshape(values.shape),
        )

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
        v0 = np.empty_like(xf)
        region = np.empty(xf.shape, dtype=object)

        inherited = xf <= self.X_3
        if np.any(inherited):
            vals = self.current.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(vals["F_current_rf40_lambda_turn"], dtype=float)
            U[inherited] = np.asarray(vals["U_current_rf40_lambda_turn"], dtype=float)
            E[inherited] = np.asarray(vals["E_current_rf40_lambda_turn"], dtype=float)
            m_ratio[inherited] = np.asarray(
                vals["M_over_X_current_rf40_lambda_turn"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                vals["M_eta_over_X_current_rf40_lambda_turn"], dtype=float
            )
            v0[inherited] = np.asarray(vals["v0_current_rf40_lambda_turn"], dtype=float)
            region[inherited] = "current_cartesian_leading_through_RF40_lambda_turn"

        power = ~inherited
        if np.any(power):
            Xo = xf[power]
            etao = ef[power]
            base = self.outer_base.profile_values(Xo, etao)
            stages = np.asarray(base["stage"], dtype=object)
            if np.any(stages != "power_law"):
                raise RuntimeError(
                    "X_w<X<=X_4 must remain on the public RF40 power_law stage"
                )

            seam_m, seam_m_eta = self._X3_memory(etao)
            ratio = self.X_3 / Xo
            mo = ratio * seam_m
            meo = ratio * seam_m_eta

            axis = self.geometry.physical_profiles.axis_profiles.values(
                np.zeros_like(etao), etao
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            Uo = np.asarray(base["U"], dtype=float)
            v0o = (
                2.0 * etao * Uo
                - 2.0 * self.D * etao * mo
                - d * meo
            ) / L

            F[power] = np.asarray(base["F"], dtype=float)
            U[power] = Uo
            E[power] = np.asarray(base["E"], dtype=float)
            m_ratio[power] = mo
            m_eta_ratio[power] = meo
            v0[power] = v0o
            region[power] = "public_RF40_power_law"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current RF40 power-law profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current RF40 power-law profile lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_rf40_power_law": F.reshape(shape),
            "U_current_rf40_power_law": U.reshape(shape),
            "E_current_rf40_power_law": E.reshape(shape),
            "M_over_X_current_rf40_power_law": m_ratio.reshape(shape),
            "M_eta_over_X_current_rf40_power_law": m_eta_ratio.reshape(shape),
            "v0_current_rf40_power_law": v0.reshape(shape),
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

        inherited = xf <= self.X_3
        if np.any(inherited):
            vals = self.current.similarity_radial_derivatives(xf[inherited], ef[inherited])
            F_X[inherited] = np.asarray(vals["F_current_rf40_lambda_turn_X"], dtype=float)
            U_X[inherited] = np.asarray(vals["U_current_rf40_lambda_turn_X"], dtype=float)
            E_X[inherited] = np.asarray(vals["E_current_rf40_lambda_turn_X"], dtype=float)

        power = ~inherited
        if np.any(power):
            vals = self.outer_base.profile_values(xf[power], ef[power])
            stages = np.asarray(vals["stage"], dtype=object)
            if np.any(stages != "power_law"):
                raise RuntimeError(
                    "RF40 derivative evaluation escaped the power_law stage"
                )
            F_X[power] = np.asarray(vals["F_X"], dtype=float)
            U_X[power] = np.asarray(vals["U_X"], dtype=float)
            E_X[power] = np.asarray(vals["E_X"], dtype=float)

        if any(np.any(~np.isfinite(arr)) for arr in (F_X, U_X, E_X)):
            raise RuntimeError("current RF40 power-law radial derivatives are non-finite")
        return {
            "F_current_rf40_power_law_X": F_X.reshape(shape),
            "U_current_rf40_power_law_X": U_X.reshape(shape),
            "E_current_rf40_power_law_X": E_X.reshape(shape),
        }

    def values(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        radial_tol = 128.0 * np.finfo(float).eps * self.X_4
        if np.any(X > self.X_4 + radial_tol):
            raise ValueError(
                "point lies beyond current RF40 power-law endpoint; "
                "cone/I1-I4 overlays and later outer assembly are not materialized"
            )
        X_eval = np.minimum(X, self.X_4)
        profile = self.similarity_profile_values(X_eval, coords["eta"])
        q = np.asarray(coords["q"], dtype=float)
        q_radial = profile["v0_current_rf40_power_law"] / (2.0 * q)
        q_swirl = np.power(q, -self.A - 0.5) * profile["F_current_rf40_power_law"]
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = np.power(q, -self.A) * profile["U_current_rf40_power_law"]
        return {**coords, **profile, "u": u, "v": v, "w": w}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {"schema": SCHEMA, "current": self.current.configuration()}

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianRF40PowerLaw":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current Cartesian RF40 power-law schema")
        return cls(
            current=KokunoPA16CurrentCartesianRF40LambdaTurn.from_configuration(
                payload["current"]
            )
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianRF40PowerLaw":
        return cls.from_configuration(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

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
        left = self.current.similarity_profile_values(self.X_3, eta)
        public_3 = self.outer_base.profile_values(self.X_3, eta)
        left_jet = self.current.similarity_radial_derivatives(self.X_3, eta)
        at_4 = self.similarity_profile_values(self.X_4, eta)
        at_4_jet = self.similarity_radial_derivatives(self.X_4, eta)
        public_4 = self.outer_base.profile_values(self.X_4, eta)

        def scalar(value: Any) -> float:
            return float(np.asarray(value))

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
                "X_3": self.X_3,
                "X_4": self.X_4,
                "log_X_3": self.log_X_3,
                "log_X_4": self.log_X_4,
                "T_w": self.T_w,
                "lambda_outer": self.lambda_outer,
                "eta_interval": list(self.eta_interval),
                "time_interval": list(self.geometry.time_interval),
            },
            "X3_seam_probe_eta": eta,
            "X3_value_jump_public_minus_current": {
                "F": scalar(public_3["F"]) - scalar(left["F_current_rf40_lambda_turn"]),
                "U": scalar(public_3["U"]) - scalar(left["U_current_rf40_lambda_turn"]),
                "E": scalar(public_3["E"]) - scalar(left["E_current_rf40_lambda_turn"]),
            },
            "X3_first_radial_jet_jump_public_minus_current": {
                "F_X": scalar(public_3["F_X"]) - scalar(left_jet["F_current_rf40_lambda_turn_X"]),
                "U_X": scalar(public_3["U_X"]) - scalar(left_jet["U_current_rf40_lambda_turn_X"]),
                "E_X": scalar(public_3["E_X"]) - scalar(left_jet["E_current_rf40_lambda_turn_X"]),
            },
            "X4_probe": {
                "F": scalar(at_4["F_current_rf40_power_law"]),
                "U": scalar(at_4["U_current_rf40_power_law"]),
                "E": scalar(at_4["E_current_rf40_power_law"]),
                "public_stage": str(np.asarray(public_4["stage"]).item()),
                "public_ell": scalar(public_4["ell"]),
                "M_over_X": scalar(at_4["M_over_X_current_rf40_power_law"]),
                "M_eta_over_X": scalar(at_4["M_eta_over_X_current_rf40_power_law"]),
                "v0": scalar(at_4["v0_current_rf40_power_law"]),
                "X_F_X_over_F": self.X_4
                * scalar(at_4_jet["F_current_rf40_power_law_X"])
                / scalar(at_4["F_current_rf40_power_law"]),
                "X_E_X_over_E": self.X_4
                * scalar(at_4_jet["E_current_rf40_power_law_X"])
                / scalar(at_4["E_current_rf40_power_law"]),
                "X_U_X": self.X_4 * scalar(at_4_jet["U_current_rf40_power_law_X"]),
            },
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()

    candidate = KokunoPA16CurrentCartesianRF40PowerLaw()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
