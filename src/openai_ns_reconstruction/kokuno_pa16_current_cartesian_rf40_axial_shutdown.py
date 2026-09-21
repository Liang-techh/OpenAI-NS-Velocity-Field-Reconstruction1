"""Current Kokuno Cartesian leading field through the RF40 axial-shutdown stage.

Pinned provenance is the corrected 2026-09-09 KokunoYumeto public
reconstruction at commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.

Agent 1 already materializes the current candidate-side Cartesian leading field
through the first post-``X_R`` RF40 turn, ending at

    X_1 = e X_R.

The corrected RF40b schedule then resets the local logarithmic coordinate

    y = log(X/X_1)

and on ``0 <= y <= T_d`` prescribes

    U = k(y) eta,
    ell = d_y log(sqrt(2X) E) = 0,
    k(y) = 4 [1-sigma(log(1+y)/M_d)].

The repository already contains this public schedule in
``KokunoOuterBaseSchedule``.  This adapter reuses it rather than creating a
second RF40 implementation.

The only current-lineage bookkeeping required here is the incompressibility
primitive.  Let ``b(y)`` be the public base coefficient returned by
``KokunoOuterBaseSchedule`` so that the public ideal-start primitive is

    (M/X)_base = eta b(y),        (M_eta/X)_base = b(y).

Because the current and public fields use exactly the same ``U`` throughout
this stage, ``M_current-M_base`` is constant in physical X.  Therefore the
actual #986 seam memory propagates exactly as

    M/X = eta b(y) + (X_1/X) [(M/X)_1 - 4 eta],
    M_eta/X = b(y) + (X_1/X) [(M_eta/X)_1 - 4].

The Cartesian radial profile remains

    v0 = (2 eta U - 2 D eta M/X - d M_eta/X) / L.

This is one deliberately small executable increment.  It stops at

    X_2 = X_R exp(1+T_d)

and fails closed beyond that point.  The later RF40 lambda turn, power-law
stage, cone modulation, I1--I4 overlays, globally matched pressure/forcing,
complete held-out NS residual, paper-exact field and OpenAI-field identity
remain unmaterialized.
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
from .kokuno_pa16_current_cartesian_rf40_first_turn import (
    KokunoPA16CurrentCartesianRF40FirstTurn,
)


SCHEMA = "kokuno-pa16-current-cartesian-rf40-axial-shutdown-v1"

_SOURCE_FORMULAS = {
    "coordinates": (
        "RF40b resets y=log(X/X_1) at X_1=e X_R for the axial-shutdown stage"
    ),
    "axial_shutdown_domain": (
        "0<=y<=T_d; X_1<=X<=X_2 with X_2=X_R exp(1+T_d)"
    ),
    "axial_shutdown_U": (
        "U=k(y) eta; k(y)=4[1-sigma(log(1+y)/M_d)]"
    ),
    "axial_shutdown_ell": "ell=d_y log(sqrt(2X)E)=0",
    "log_E_evolution": "d_y log E=ell-1/2=-1/2",
    "F": "F=E/sqrt(2X)",
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "current_memory_propagation": (
        "if b=(M_base/eta)/X from the public schedule, then "
        "M_current/X=eta*b+(X_1/X)[(M/X)_1-4eta] and "
        "M_eta_current/X=b+(X_1/X)[(M_eta/X)_1-4]"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "X1_seam": "consume the exact current #986 candidate-side X_1 handoff",
    "RF40_profile": (
        "reuse KokunoOuterBaseSchedule axial_shutdown values/derivatives and "
        "its fixed source-step quadrature; do not duplicate the public schedule"
    ),
    "primitive_memory": (
        "carry the #986 floating-point M/X and d_eta(M/X) seam by an exact "
        "difference-from-public-base identity; do not reset to ideal moments"
    ),
    "coordinates": "reuse the governed current A1 q/X/eta inverse and Cartesian map",
    "new_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_cartesian_leading_through_RF40_first_turn_consumed": True,
    "public_RF40_axial_shutdown_formula_reused": True,
    "current_lineage_RF40_axial_shutdown_materialized": True,
    "current_cartesian_spacetime_leading_velocity_through_RF40_axial_shutdown_materialized": True,
    "current_incompressibility_memory_carried_through_RF40_axial_shutdown": True,
    "axis_regular_cartesian_formula_reused": True,
    "velocity_interface_vectorized": True,
    "velocity_configuration_serializable": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_exact_current_lineage_certified": False,
    "full_post_XR_RF40_current_lineage_materialized": False,
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
class KokunoPA16CurrentCartesianRF40AxialShutdown:
    """Extend the current candidate-side Cartesian leading field one RF40 stage."""

    current: KokunoPA16CurrentCartesianRF40FirstTurn = field(
        default_factory=KokunoPA16CurrentCartesianRF40FirstTurn,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.current, KokunoPA16CurrentCartesianRF40FirstTurn):
            raise TypeError(
                "current must be KokunoPA16CurrentCartesianRF40FirstTurn"
            )
        if not math.isclose(
            self.log_X_1, self.outer_base.log_stage1_end,
            rel_tol=0.0, abs_tol=2.0e-12,
        ):
            raise ValueError("current X_1 and public RF40 stage-1 endpoint disagree")
        if not math.isclose(
            self.log_X_2, self.outer_base.log_stage2_end,
            rel_tol=0.0, abs_tol=2.0e-12,
        ):
            raise ValueError("RF40 axial-shutdown endpoint disagrees with public schedule")
        if not (0.0 < self.X_1 < self.X_2 < np.finfo(float).max):
            raise ValueError("materialized current RF40 axial-shutdown interval is invalid")

    @property
    def outer_base(self):
        return self.current.outer_base

    @property
    def X_R(self) -> float:
        return float(self.current.X_R)

    @property
    def X_1(self) -> float:
        return float(self.current.X_1)

    @property
    def log_X_1(self) -> float:
        return float(self.current.log_X_1)

    @property
    def log_X_2(self) -> float:
        return float(self.outer_base.log_stage2_end)

    @property
    def X_2(self) -> float:
        return float(math.exp(self.log_X_2))

    @property
    def T_d(self) -> float:
        return float(self.outer_base.outer_schedule.T_d)

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

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.current.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(
            _finite(X, "X"), _finite(eta, "eta")
        )
        if np.any((X_arr < 0.0) | (X_arr > self.X_2)):
            raise ValueError(
                "X must lie in the current RF40-axial-shutdown domain "
                f"[0,{self.X_2:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _X1_memory(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return current ``(M/X)_1`` and ``d_eta(M/X)_1`` at the exact X_1 seam."""
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        seam = self.current.similarity_profile_values(
            np.full(unique.shape, self.X_1, dtype=float), unique
        )
        m = np.asarray(
            seam["M_over_X_current_rf40_first_turn"], dtype=float
        )
        m_eta = np.asarray(
            seam["M_eta_over_X_current_rf40_first_turn"], dtype=float
        )
        return (
            m[inverse].reshape(values.shape),
            m_eta[inverse].reshape(values.shape),
        )

    def similarity_profile_values(
        self, X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        """Return unified ``F,U,E,M/X,M_eta/X,v0`` through RF40 axial shutdown."""
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

        inherited = xf <= self.X_1
        if np.any(inherited):
            vals = self.current.similarity_profile_values(
                xf[inherited], ef[inherited]
            )
            F[inherited] = np.asarray(
                vals["F_current_rf40_first_turn"], dtype=float
            )
            U[inherited] = np.asarray(
                vals["U_current_rf40_first_turn"], dtype=float
            )
            E[inherited] = np.asarray(
                vals["E_current_rf40_first_turn"], dtype=float
            )
            m_ratio[inherited] = np.asarray(
                vals["M_over_X_current_rf40_first_turn"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                vals["M_eta_over_X_current_rf40_first_turn"], dtype=float
            )
            v0[inherited] = np.asarray(
                vals["v0_current_rf40_first_turn"], dtype=float
            )
            region[inherited] = "current_cartesian_leading_through_RF40_first_turn"

        axial = ~inherited
        if np.any(axial):
            Xo = xf[axial]
            etao = ef[axial]
            base = self.outer_base.profile_values(Xo, etao)
            stages = np.asarray(base["stage"], dtype=object)
            if np.any(stages != "axial_shutdown"):
                raise RuntimeError(
                    "X_1<X<=X_2 must remain on the public RF40 axial_shutdown stage"
                )

            seam_m, seam_m_eta = self._X1_memory(etao)
            b = np.asarray(base["m_ratio"], dtype=float)
            ratio = self.X_1 / Xo
            mo = etao * b + ratio * (seam_m - 4.0 * etao)
            meo = b + ratio * (seam_m_eta - 4.0)

            axis = (
                self.current.current.current.geometry.physical_profiles
                .axis_profiles.values(np.zeros_like(etao), etao)
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            Uo = np.asarray(base["U"], dtype=float)
            v0o = (
                2.0 * etao * Uo
                - 2.0 * self.D * etao * mo
                - d * meo
            ) / L

            F[axial] = np.asarray(base["F"], dtype=float)
            U[axial] = Uo
            E[axial] = np.asarray(base["E"], dtype=float)
            m_ratio[axial] = mo
            m_eta_ratio[axial] = meo
            v0[axial] = v0o
            region[axial] = "public_RF40_axial_shutdown"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError(
                "current RF40 axial-shutdown profile produced non-finite values"
            )
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError(
                "current RF40 axial-shutdown profile lost positive F/E"
            )
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_rf40_axial_shutdown": F.reshape(shape),
            "U_current_rf40_axial_shutdown": U.reshape(shape),
            "E_current_rf40_axial_shutdown": E.reshape(shape),
            "M_over_X_current_rf40_axial_shutdown": m_ratio.reshape(shape),
            "M_eta_over_X_current_rf40_axial_shutdown": m_eta_ratio.reshape(shape),
            "v0_current_rf40_axial_shutdown": v0.reshape(shape),
        }

    def similarity_radial_derivatives(
        self, X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        """Return first physical-X derivatives of ``F,U,E`` through ``X_2``."""
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        if np.any(X_arr <= 0.0):
            raise ValueError("similarity_radial_derivatives requires X>0")
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        F_X = np.empty_like(xf)
        U_X = np.empty_like(xf)
        E_X = np.empty_like(xf)

        inherited = xf <= self.X_1
        if np.any(inherited):
            vals = self.current.similarity_radial_derivatives(
                xf[inherited], ef[inherited]
            )
            F_X[inherited] = np.asarray(
                vals["F_current_rf40_first_turn_X"], dtype=float
            )
            U_X[inherited] = np.asarray(
                vals["U_current_rf40_first_turn_X"], dtype=float
            )
            E_X[inherited] = np.asarray(
                vals["E_current_rf40_first_turn_X"], dtype=float
            )

        axial = ~inherited
        if np.any(axial):
            vals = self.outer_base.profile_values(
                xf[axial], ef[axial]
            )
            stages = np.asarray(vals["stage"], dtype=object)
            if np.any(stages != "axial_shutdown"):
                raise RuntimeError(
                    "RF40 derivative evaluation escaped the axial_shutdown stage"
                )
            F_X[axial] = np.asarray(vals["F_X"], dtype=float)
            U_X[axial] = np.asarray(vals["U_X"], dtype=float)
            E_X[axial] = np.asarray(vals["E_X"], dtype=float)

        if any(
            np.any(~np.isfinite(arr)) for arr in (F_X, U_X, E_X)
        ):
            raise RuntimeError(
                "current RF40 axial-shutdown radial derivatives are non-finite"
            )
        return {
            "F_current_rf40_axial_shutdown_X": F_X.reshape(shape),
            "U_current_rf40_axial_shutdown_X": U_X.reshape(shape),
            "E_current_rf40_axial_shutdown_X": E_X.reshape(shape),
        }

    def values(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        """Return source coordinates, profiles and Cartesian velocity through ``X_2``."""
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        radial_tol = 128.0 * np.finfo(float).eps * self.X_2
        if np.any(X > self.X_2 + radial_tol):
            raise ValueError(
                "point lies beyond current RF40 axial-shutdown endpoint; "
                "later RF40 is not materialized"
            )
        X_eval = np.minimum(X, self.X_2)
        profile = self.similarity_profile_values(
            X_eval, coords["eta"]
        )
        q = np.asarray(coords["q"], dtype=float)
        q_radial = (
            profile["v0_current_rf40_axial_shutdown"] / (2.0 * q)
        )
        q_swirl = (
            np.power(q, -self.A - 0.5)
            * profile["F_current_rf40_axial_shutdown"]
        )
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = (
            np.power(q, -self.A)
            * profile["U_current_rf40_axial_shutdown"]
        )
        return {**coords, **profile, "u": u, "v": v, "w": w}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Unified vectorized ``velocity(x,y,z,t)->[...,3]`` through axial shutdown."""
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {"schema": SCHEMA, "current": self.current.configuration()}

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianRF40AxialShutdown":
        if (
            not isinstance(payload, Mapping)
            or payload.get("schema") != SCHEMA
        ):
            raise ValueError(
                "unexpected current Cartesian RF40 axial-shutdown schema"
            )
        return cls(
            current=KokunoPA16CurrentCartesianRF40FirstTurn.from_configuration(
                payload["current"]
            )
        )

    def save_configuration(
        self, path: str | Path
    ) -> dict[str, Any]:
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
    ) -> "KokunoPA16CurrentCartesianRF40AxialShutdown":
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
        at_1 = self.current.similarity_profile_values(self.X_1, eta)
        right_1 = self.outer_base.profile_values(self.X_1, eta)
        at_1_jet = self.current.similarity_radial_derivatives(
            self.X_1, eta
        )
        right_1_jet = self.outer_base.profile_values(self.X_1, eta)
        at_2 = self.similarity_profile_values(self.X_2, eta)
        at_2_jet = self.similarity_radial_derivatives(self.X_2, eta)
        public_2 = self.outer_base.profile_values(self.X_2, eta)

        def scalar(value: Any) -> float:
            return float(np.asarray(value))

        value_jump = {
            "F": scalar(right_1["F"])
            - scalar(at_1["F_current_rf40_first_turn"]),
            "U": scalar(right_1["U"])
            - scalar(at_1["U_current_rf40_first_turn"]),
            "E": scalar(right_1["E"])
            - scalar(at_1["E_current_rf40_first_turn"]),
        }
        jet_jump = {
            "F_X": scalar(right_1_jet["F_X"])
            - scalar(at_1_jet["F_current_rf40_first_turn_X"]),
            "U_X": scalar(right_1_jet["U_X"])
            - scalar(at_1_jet["U_current_rf40_first_turn_X"]),
            "E_X": scalar(right_1_jet["E_X"])
            - scalar(at_1_jet["E_current_rf40_first_turn_X"]),
        }

        b2 = scalar(public_2["m_ratio"])
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
            "numerical_realization": copy.deepcopy(
                _NUMERICAL_REALIZATION
            ),
            "domain": {
                "X_1": self.X_1,
                "X_2": self.X_2,
                "log_X_1": self.log_X_1,
                "log_X_2": self.log_X_2,
                "T_d": self.T_d,
                "eta_interval": list(self.eta_interval),
                "time_interval": list(
                    self.current.current.current.geometry.time_interval
                ),
            },
            "X1_seam_probe_eta": eta,
            "X1_value_jump_public_right_minus_current_left": value_jump,
            "X1_first_radial_jet_jump_public_right_minus_current_left": jet_jump,
            "X2_probe": {
                "F": scalar(
                    at_2["F_current_rf40_axial_shutdown"]
                ),
                "U": scalar(
                    at_2["U_current_rf40_axial_shutdown"]
                ),
                "E": scalar(
                    at_2["E_current_rf40_axial_shutdown"]
                ),
                "public_k": scalar(public_2["k"]),
                "public_m_ratio_coefficient": b2,
                "M_over_X_minus_public_base": (
                    scalar(
                        at_2["M_over_X_current_rf40_axial_shutdown"]
                    )
                    - eta * b2
                ),
                "M_eta_over_X_minus_public_base": (
                    scalar(
                        at_2[
                            "M_eta_over_X_current_rf40_axial_shutdown"
                        ]
                    )
                    - b2
                ),
                "v0": scalar(
                    at_2["v0_current_rf40_axial_shutdown"]
                ),
                "X_F_X_over_F": self.X_2
                * scalar(
                    at_2_jet[
                        "F_current_rf40_axial_shutdown_X"
                    ]
                )
                / scalar(
                    at_2["F_current_rf40_axial_shutdown"]
                ),
                "X_E_X_over_E": self.X_2
                * scalar(
                    at_2_jet[
                        "E_current_rf40_axial_shutdown_X"
                    ]
                )
                / scalar(
                    at_2["E_current_rf40_axial_shutdown"]
                ),
                "X_U_X": self.X_2
                * scalar(
                    at_2_jet[
                        "U_current_rf40_axial_shutdown_X"
                    ]
                ),
            },
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", required=True,
        help="write deterministic scoped receipt JSON",
    )
    parser.add_argument(
        "--config-output", required=True,
        help="write deterministic configuration JSON",
    )
    args = parser.parse_args()

    candidate = KokunoPA16CurrentCartesianRF40AxialShutdown()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
