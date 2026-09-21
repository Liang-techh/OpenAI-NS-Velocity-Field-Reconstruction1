"""Current Kokuno Cartesian leading field through the RF40 lambda-turn stage.

Pinned provenance is the corrected 2026-09-09 KokunoYumeto public
reconstruction at commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.

The parent Agent-1 candidate reaches the RF40 axial-shutdown endpoint

    X_2 = X_R exp(1+T_d).

The corrected RF40b schedule then resets the local logarithmic coordinate

    y = log(X/X_2)

and, for ``0 <= y <= 1``, prescribes

    U = 0,
    ell = d_y log(sqrt(2X) E) = -lambda sigma(y).

The repository already implements this public schedule in
``KokunoOuterBaseSchedule``.  This adapter reuses that implementation and
carries the actual current-lineage incompressibility primitive.  Since U=0 on
this interval, physical M and M_eta are constant, hence

    M/X = (X_2/X) (M/X)_2,
    M_eta/X = (X_2/X) (M_eta/X)_2.

No hidden source parameter is inferred and no residual knob is introduced.
This minimal increment stops at ``X_3=X_w=X_R exp(T_d+2)``.  The following
constant-lambda power-law stage and cone/I1--I4 overlays remain unmaterialized.
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
from .kokuno_pa16_current_cartesian_rf40_axial_shutdown import (
    KokunoPA16CurrentCartesianRF40AxialShutdown,
)


SCHEMA = "kokuno-pa16-current-cartesian-rf40-lambda-turn-v1"

_SOURCE_FORMULAS = {
    "coordinates": "RF40b resets y=log(X/X_2) at the lambda-turn left endpoint",
    "lambda_turn_domain": "0<=y<=1; X_2<=X<=X_3 with X_3=X_w=X_R exp(T_d+2)",
    "lambda_turn_U": "U=0",
    "lambda_turn_ell": "ell=d_y log(sqrt(2X)E)=-lambda*sigma(y)",
    "log_E_evolution": "d_y log E=ell-1/2",
    "F": "F=E/sqrt(2X)",
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "current_memory_propagation": (
        "because U=0, M and M_eta are constant: "
        "M/X=(X_2/X)(M/X)_2 and M_eta/X=(X_2/X)(M_eta/X)_2"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "X2_seam": "consume the exact current Agent-1 axial-shutdown X_2 handoff",
    "RF40_profile": (
        "reuse KokunoOuterBaseSchedule lambda_turn values/derivatives and its "
        "fixed source-step quadrature; do not duplicate the public schedule"
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
    "current_cartesian_leading_through_RF40_axial_shutdown_consumed": True,
    "public_RF40_lambda_turn_formula_reused": True,
    "current_lineage_RF40_lambda_turn_materialized": True,
    "current_cartesian_spacetime_leading_velocity_through_RF40_lambda_turn_materialized": True,
    "current_incompressibility_memory_carried_through_RF40_lambda_turn": True,
    "axis_regular_cartesian_formula_reused": True,
    "velocity_interface_vectorized": True,
    "velocity_configuration_serializable": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_exact_current_lineage_certified": False,
    "full_post_XR_RF40_current_lineage_materialized": False,
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
class KokunoPA16CurrentCartesianRF40LambdaTurn:
    """Extend the current candidate-side Cartesian leading field one RF40 stage."""

    current: KokunoPA16CurrentCartesianRF40AxialShutdown = field(
        default_factory=KokunoPA16CurrentCartesianRF40AxialShutdown,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.current, KokunoPA16CurrentCartesianRF40AxialShutdown):
            raise TypeError(
                "current must be KokunoPA16CurrentCartesianRF40AxialShutdown"
            )
        if not math.isclose(
            self.log_X_2,
            self.outer_base.log_stage2_end,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("current X_2 and public RF40 stage-2 endpoint disagree")
        if not math.isclose(
            self.log_X_3,
            self.outer_base.log_X_w,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("RF40 lambda-turn endpoint disagrees with public schedule")
        if not math.isclose(
            self.log_X_3 - self.log_X_2,
            1.0,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("RF40 lambda-turn must have unit logarithmic length")
        if not (0.0 < self.X_2 < self.X_3 < np.finfo(float).max):
            raise ValueError("materialized current RF40 lambda-turn interval is invalid")

    @property
    def outer_base(self):
        return self.current.outer_base

    @property
    def X_R(self) -> float:
        return float(self.current.X_R)

    @property
    def X_2(self) -> float:
        return float(self.current.X_2)

    @property
    def log_X_2(self) -> float:
        return float(self.current.log_X_2)

    @property
    def log_X_3(self) -> float:
        return float(self.outer_base.log_X_w)

    @property
    def X_3(self) -> float:
        return float(math.exp(self.log_X_3))

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
        return self.current.current.current.current.geometry

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.current.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_3)):
            raise ValueError(
                "X must lie in the current RF40-lambda-turn domain "
                f"[0,{self.X_3:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _X2_memory(self, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        seam = self.current.similarity_profile_values(
            np.full(unique.shape, self.X_2, dtype=float), unique
        )
        m = np.asarray(
            seam["M_over_X_current_rf40_axial_shutdown"], dtype=float
        )
        m_eta = np.asarray(
            seam["M_eta_over_X_current_rf40_axial_shutdown"], dtype=float
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

        inherited = xf <= self.X_2
        if np.any(inherited):
            vals = self.current.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(
                vals["F_current_rf40_axial_shutdown"], dtype=float
            )
            U[inherited] = np.asarray(
                vals["U_current_rf40_axial_shutdown"], dtype=float
            )
            E[inherited] = np.asarray(
                vals["E_current_rf40_axial_shutdown"], dtype=float
            )
            m_ratio[inherited] = np.asarray(
                vals["M_over_X_current_rf40_axial_shutdown"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                vals["M_eta_over_X_current_rf40_axial_shutdown"], dtype=float
            )
            v0[inherited] = np.asarray(
                vals["v0_current_rf40_axial_shutdown"], dtype=float
            )
            region[inherited] = "current_cartesian_leading_through_RF40_axial_shutdown"

        turn = ~inherited
        if np.any(turn):
            Xo = xf[turn]
            etao = ef[turn]
            base = self.outer_base.profile_values(Xo, etao)
            stages = np.asarray(base["stage"], dtype=object)
            if np.any(stages != "lambda_turn"):
                raise RuntimeError(
                    "X_2<X<=X_3 must remain on the public RF40 lambda_turn stage"
                )

            seam_m, seam_m_eta = self._X2_memory(etao)
            ratio = self.X_2 / Xo
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

            F[turn] = np.asarray(base["F"], dtype=float)
            U[turn] = Uo
            E[turn] = np.asarray(base["E"], dtype=float)
            m_ratio[turn] = mo
            m_eta_ratio[turn] = meo
            v0[turn] = v0o
            region[turn] = "public_RF40_lambda_turn"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current RF40 lambda-turn profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current RF40 lambda-turn profile lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_rf40_lambda_turn": F.reshape(shape),
            "U_current_rf40_lambda_turn": U.reshape(shape),
            "E_current_rf40_lambda_turn": E.reshape(shape),
            "M_over_X_current_rf40_lambda_turn": m_ratio.reshape(shape),
            "M_eta_over_X_current_rf40_lambda_turn": m_eta_ratio.reshape(shape),
            "v0_current_rf40_lambda_turn": v0.reshape(shape),
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

        inherited = xf <= self.X_2
        if np.any(inherited):
            vals = self.current.similarity_radial_derivatives(
                xf[inherited], ef[inherited]
            )
            F_X[inherited] = np.asarray(
                vals["F_current_rf40_axial_shutdown_X"], dtype=float
            )
            U_X[inherited] = np.asarray(
                vals["U_current_rf40_axial_shutdown_X"], dtype=float
            )
            E_X[inherited] = np.asarray(
                vals["E_current_rf40_axial_shutdown_X"], dtype=float
            )

        turn = ~inherited
        if np.any(turn):
            vals = self.outer_base.profile_values(xf[turn], ef[turn])
            stages = np.asarray(vals["stage"], dtype=object)
            if np.any(stages != "lambda_turn"):
                raise RuntimeError(
                    "RF40 derivative evaluation escaped the lambda_turn stage"
                )
            F_X[turn] = np.asarray(vals["F_X"], dtype=float)
            U_X[turn] = np.asarray(vals["U_X"], dtype=float)
            E_X[turn] = np.asarray(vals["E_X"], dtype=float)

        if any(np.any(~np.isfinite(arr)) for arr in (F_X, U_X, E_X)):
            raise RuntimeError("current RF40 lambda-turn radial derivatives are non-finite")
        return {
            "F_current_rf40_lambda_turn_X": F_X.reshape(shape),
            "U_current_rf40_lambda_turn_X": U_X.reshape(shape),
            "E_current_rf40_lambda_turn_X": E_X.reshape(shape),
        }

    def values(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        radial_tol = 128.0 * np.finfo(float).eps * self.X_3
        if np.any(X > self.X_3 + radial_tol):
            raise ValueError(
                "point lies beyond current RF40 lambda-turn endpoint; "
                "the constant-lambda power-law stage is not materialized"
            )
        X_eval = np.minimum(X, self.X_3)
        profile = self.similarity_profile_values(X_eval, coords["eta"])
        q = np.asarray(coords["q"], dtype=float)
        q_radial = profile["v0_current_rf40_lambda_turn"] / (2.0 * q)
        q_swirl = (
            np.power(q, -self.A - 0.5)
            * profile["F_current_rf40_lambda_turn"]
        )
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = np.power(q, -self.A) * profile["U_current_rf40_lambda_turn"]
        return {**coords, **profile, "u": u, "v": v, "w": w}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {"schema": SCHEMA, "current": self.current.configuration()}

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianRF40LambdaTurn":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current Cartesian RF40 lambda-turn schema")
        return cls(
            current=KokunoPA16CurrentCartesianRF40AxialShutdown.from_configuration(
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
    ) -> "KokunoPA16CurrentCartesianRF40LambdaTurn":
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
        return hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta = 0.25
        left = self.current.similarity_profile_values(self.X_2, eta)
        public_2 = self.outer_base.profile_values(self.X_2, eta)
        left_jet = self.current.similarity_radial_derivatives(self.X_2, eta)
        at_3 = self.similarity_profile_values(self.X_3, eta)
        at_3_jet = self.similarity_radial_derivatives(self.X_3, eta)
        public_3 = self.outer_base.profile_values(self.X_3, eta)

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
                "X_2": self.X_2,
                "X_3": self.X_3,
                "log_X_2": self.log_X_2,
                "log_X_3": self.log_X_3,
                "lambda_outer": self.lambda_outer,
                "eta_interval": list(self.eta_interval),
                "time_interval": list(self.geometry.time_interval),
            },
            "X2_seam_probe_eta": eta,
            "X2_value_jump_public_minus_current": {
                "F": scalar(public_2["F"])
                - scalar(left["F_current_rf40_axial_shutdown"]),
                "U": scalar(public_2["U"])
                - scalar(left["U_current_rf40_axial_shutdown"]),
                "E": scalar(public_2["E"])
                - scalar(left["E_current_rf40_axial_shutdown"]),
            },
            "X2_first_radial_jet_jump_public_minus_current": {
                "F_X": scalar(public_2["F_X"])
                - scalar(left_jet["F_current_rf40_axial_shutdown_X"]),
                "U_X": scalar(public_2["U_X"])
                - scalar(left_jet["U_current_rf40_axial_shutdown_X"]),
                "E_X": scalar(public_2["E_X"])
                - scalar(left_jet["E_current_rf40_axial_shutdown_X"]),
            },
            "X3_probe": {
                "F": scalar(at_3["F_current_rf40_lambda_turn"]),
                "U": scalar(at_3["U_current_rf40_lambda_turn"]),
                "E": scalar(at_3["E_current_rf40_lambda_turn"]),
                "public_ell": scalar(public_3["ell"]),
                "M_over_X": scalar(at_3["M_over_X_current_rf40_lambda_turn"]),
                "M_eta_over_X": scalar(
                    at_3["M_eta_over_X_current_rf40_lambda_turn"]
                ),
                "v0": scalar(at_3["v0_current_rf40_lambda_turn"]),
                "X_F_X_over_F": self.X_3
                * scalar(at_3_jet["F_current_rf40_lambda_turn_X"])
                / scalar(at_3["F_current_rf40_lambda_turn"]),
                "X_E_X_over_E": self.X_3
                * scalar(at_3_jet["E_current_rf40_lambda_turn_X"])
                / scalar(at_3["E_current_rf40_lambda_turn"]),
                "X_U_X": self.X_3
                * scalar(at_3_jet["U_current_rf40_lambda_turn_X"]),
            },
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()

    candidate = KokunoPA16CurrentCartesianRF40LambdaTurn()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
