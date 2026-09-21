"""Current Kokuno Cartesian leading field through the public exterior-preservation interval.

Pinned provenance is the corrected 2026-09-09 KokunoYumeto reconstruction at
commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.

Agent 1 already materializes the current candidate-side joined profile and its
Cartesian spacetime velocity through

    X_h = X_R exp(-5).

The corrected reconstruction arranges the PA.16 exit so that, on a
neighborhood of ``X_h``, the field agrees with the ideal outer pair and the
five prefix moments agree.  The public temporary-reference branch continuing
toward ``X_R`` is

    U = 4 eta,
    E = P_* (1+eta^2)^(-1) (X/X_R)^(1/10),
    F = E / sqrt(2X).

This module makes only that ``X_h <= X <= X_R`` exterior-preservation bridge
executable on the *current* A1 ancestry.  It reuses the existing
``KokunoOuterBaseSchedule`` temporary-reference evaluator and the governed
Cartesian ``q/X/eta`` map from the current #965 field.

A small but important numerical detail is kept explicit.  The current PA.16
solve is a floating-point candidate realization, so this adapter does not
silently reset the incompressibility primitive to its ideal value at ``X_h``.
Instead it carries the actual seam memory exactly under ``U=4 eta``:

    (M/X)(X) = 4 eta + (X_h/X) [(M/X)_h - 4 eta],
    (M_eta/X)(X) = 4 + (X_h/X) [(M_eta/X)_h - 4].

Thus the Cartesian radial component is value-compatible at ``X_h`` even if a
tiny numerical PA.16 closure remainder is present, while that memory decays
under the public exterior branch.  This is repository numerical bookkeeping,
not hidden source data and not independent PDE validation.

The result is still not the post-``X_R`` RF40 schedule with cone/I1--I4
overlays, not a globally matched pressure/forcing field, and not a paper- or
OpenAI-exact Navier--Stokes solution.
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

from .kokuno_outer_base_schedule import KokunoOuterBaseSchedule
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa16_current_cartesian_leading_velocity import (
    KokunoPA16CurrentCartesianLeadingVelocity,
)


SCHEMA = "kokuno-pa16-current-cartesian-exterior-to-xr-v1"

_SOURCE_FORMULAS = {
    "exterior_preservation": (
        "at X_h=X_R exp(-5), equality of the field on a neighborhood plus "
        "equality of (M,I,J,S,C_p) enables exterior preservation"
    ),
    "ideal_outer_pair": (
        "U=4 eta; E=P_* (1+eta^2)^(-1) (X/X_R)^(1/10); "
        "F=E/sqrt(2X), for the temporary-reference branch X<=X_R"
    ),
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "Xh_seam": "consume the exact current #965 candidate-side X_h seam",
    "primitive_memory": (
        "carry current floating-point M/X and d_eta(M/X) analytically under U=4 eta; "
        "do not reset them to ideal values"
    ),
    "outer_profile": (
        "reuse KokunoOuterBaseSchedule temporary-reference branch with the exact same "
        "current A1 outer_schedule/C/P_*/X_R identity"
    ),
    "coordinates": "reuse the governed current #965 q/X/eta inverse and Cartesian map",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_cartesian_leading_through_Xh_consumed": True,
    "public_ideal_temporary_reference_branch_reused": True,
    "candidate_side_exterior_preservation_Xh_to_XR_materialized": True,
    "current_cartesian_spacetime_leading_velocity_through_XR_materialized": True,
    "velocity_beyond_Xh_materialized": True,
    "current_incompressibility_memory_carried": True,
    "axis_regular_cartesian_formula_reused": True,
    "velocity_interface_vectorized": True,
    "velocity_configuration_serializable": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_exact_exterior_preservation_certified": False,
    "post_XR_RF40_current_lineage_materialized": False,
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
class KokunoPA16CurrentCartesianExteriorToXR:
    """Extend the current candidate-side Cartesian leading field through ``X_R``."""

    current: KokunoPA16CurrentCartesianLeadingVelocity = field(
        default_factory=KokunoPA16CurrentCartesianLeadingVelocity,
        repr=False,
        compare=False,
    )
    _outer: KokunoOuterBaseSchedule = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.current, KokunoPA16CurrentCartesianLeadingVelocity):
            raise TypeError("current must be KokunoPA16CurrentCartesianLeadingVelocity")
        if not self.current.joined.route_ready:
            raise ValueError("current joined route must be numerically geometry-feasible")

        outer = KokunoOuterBaseSchedule(
            outer_schedule=self.current.joined.moments.outer_schedule
        )
        object.__setattr__(self, "_outer", outer)

        if not math.isclose(self.h, outer.h, rel_tol=0.0, abs_tol=2.0e-15):
            raise ValueError("current Cartesian scaling h and outer schedule h disagree")
        if not math.isclose(
            math.log(self.X_R), outer.log_X_R, rel_tol=0.0, abs_tol=2.0e-12
        ):
            raise ValueError("current X_R and outer temporary-reference X_R disagree")
        if not math.isclose(
            self.X_h / self.X_R, math.exp(-5.0), rel_tol=2.0e-14, abs_tol=0.0
        ):
            raise ValueError("current X_h is not the public X_R exp(-5) join radius")
        if not (0.0 < self.X_h < self.X_R < np.finfo(float).max):
            raise ValueError("materialized current X_h/X_R interval is invalid")

    @property
    def outer_base(self) -> KokunoOuterBaseSchedule:
        return self._outer

    @property
    def X_h(self) -> float:
        return float(self.current.X_h)

    @property
    def X_R(self) -> float:
        return float(self.current.X_R)

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
        if np.any((X_arr < 0.0) | (X_arr > self.X_R)):
            raise ValueError(f"X must lie in the current exterior domain [0,{self.X_R:.17e}]")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _seam_memory(self, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return current ``(M/X)_h`` and ``d_eta(M/X)_h`` for requested etas."""
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        seam = self.current.similarity_profile_values(
            np.full(unique.shape, self.X_h, dtype=float), unique
        )
        m = np.asarray(seam["M_over_X_current_joined"], dtype=float)
        m_eta = np.asarray(seam["M_eta_over_X_current_joined"], dtype=float)
        return m[inverse].reshape(values.shape), m_eta[inverse].reshape(values.shape)

    def similarity_profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return unified ``F,U,E,M/X,M_eta/X,v0`` on ``0<=X<=X_R``."""
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

        inner = xf <= self.X_h
        if np.any(inner):
            vals = self.current.similarity_profile_values(xf[inner], ef[inner])
            F[inner] = np.asarray(vals["F_current_joined"], dtype=float)
            U[inner] = np.asarray(vals["U_current_joined"], dtype=float)
            E[inner] = np.asarray(vals["E_current_joined"], dtype=float)
            m_ratio[inner] = np.asarray(vals["M_over_X_current_joined"], dtype=float)
            m_eta_ratio[inner] = np.asarray(vals["M_eta_over_X_current_joined"], dtype=float)
            v0[inner] = np.asarray(vals["v0_current_joined"], dtype=float)
            region[inner] = "current_joined_through_Xh"

        exterior = ~inner
        if np.any(exterior):
            Xo = xf[exterior]
            etao = ef[exterior]
            base = self.outer_base.profile_values(Xo, etao)
            stages = np.asarray(base["stage"], dtype=object)
            if np.any(stages != "temporary_reference"):
                raise RuntimeError("X_h<X<=X_R must remain on the public temporary-reference branch")

            seam_m, seam_m_eta = self._seam_memory(etao)
            ratio = self.X_h / Xo
            mo = 4.0 * etao + ratio * (seam_m - 4.0 * etao)
            meo = 4.0 + ratio * (seam_m_eta - 4.0)

            axis = self.current.geometry.physical_profiles.axis_profiles.values(
                np.zeros_like(etao), etao
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            Uo = np.asarray(base["U"], dtype=float)
            v0o = (2.0 * etao * Uo - 2.0 * self.D * etao * mo - d * meo) / L

            F[exterior] = np.asarray(base["F"], dtype=float)
            U[exterior] = Uo
            E[exterior] = np.asarray(base["E"], dtype=float)
            m_ratio[exterior] = mo
            m_eta_ratio[exterior] = meo
            v0[exterior] = v0o
            region[exterior] = "public_exterior_preservation_to_XR"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current exterior profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current exterior profile lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_exterior": F.reshape(shape),
            "U_current_exterior": U.reshape(shape),
            "E_current_exterior": E.reshape(shape),
            "M_over_X_current_exterior": m_ratio.reshape(shape),
            "M_eta_over_X_current_exterior": m_eta_ratio.reshape(shape),
            "v0_current_exterior": v0.reshape(shape),
        }

    def similarity_radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return first physical-X derivatives of ``F,U,E`` through ``X_R``."""
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        if np.any(X_arr <= 0.0):
            raise ValueError("similarity_radial_derivatives requires X>0")
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        F_X = np.empty_like(xf)
        U_X = np.empty_like(xf)
        E_X = np.empty_like(xf)

        inner = xf <= self.X_h
        if np.any(inner):
            vals = self.current.joined.radial_derivatives(xf[inner], ef[inner])
            F_X[inner] = np.asarray(vals["F_current_joined_X"], dtype=float)
            U_X[inner] = np.asarray(vals["U_current_joined_X"], dtype=float)
            E_X[inner] = np.asarray(vals["E_current_joined_X"], dtype=float)

        exterior = ~inner
        if np.any(exterior):
            vals = self.outer_base.profile_values(xf[exterior], ef[exterior])
            F_X[exterior] = np.asarray(vals["F_X"], dtype=float)
            U_X[exterior] = np.asarray(vals["U_X"], dtype=float)
            E_X[exterior] = np.asarray(vals["E_X"], dtype=float)

        if np.any(~np.isfinite(F_X)) or np.any(~np.isfinite(U_X)) or np.any(~np.isfinite(E_X)):
            raise RuntimeError("current exterior radial derivatives are non-finite")
        return {
            "F_current_exterior_X": F_X.reshape(shape),
            "U_current_exterior_X": U_X.reshape(shape),
            "E_current_exterior_X": E_X.reshape(shape),
        }

    def values(self, x: Any, y: Any, z: Any, t: Any) -> dict[str, np.ndarray]:
        """Return source coordinates, unified profiles, and Cartesian velocity through ``X_R``."""
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        radial_tol = 128.0 * np.finfo(float).eps * self.X_R
        if np.any(X > self.X_R + radial_tol):
            raise ValueError("point lies beyond current X_R; post-X_R global leading is not materialized")
        X_eval = np.minimum(X, self.X_R)
        profile = self.similarity_profile_values(X_eval, coords["eta"])
        q = np.asarray(coords["q"], dtype=float)
        q_radial = profile["v0_current_exterior"] / (2.0 * q)
        q_swirl = np.power(q, -self.A - 0.5) * profile["F_current_exterior"]
        u = q_radial * coords["x"] - q_swirl * coords["y"]
        v = q_radial * coords["y"] + q_swirl * coords["x"]
        w = np.power(q, -self.A) * profile["U_current_exterior"]
        return {**coords, **profile, "u": u, "v": v, "w": w}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Unified vectorized ``velocity(x,y,z,t)->[...,3]`` through ``X_R``."""
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def configuration(self) -> dict[str, Any]:
        return {"schema": SCHEMA, "current": self.current.configuration()}

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianExteriorToXR":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current Cartesian exterior schema")
        return cls(
            current=KokunoPA16CurrentCartesianLeadingVelocity.from_configuration(
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
    def load_configuration(cls, path: str | Path) -> "KokunoPA16CurrentCartesianExteriorToXR":
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
        seam_current = self.current.similarity_profile_values(self.X_h, eta)
        seam_ideal = self.outer_base.profile_values(self.X_h, eta)
        seam_jet_current = self.current.joined.radial_derivatives(self.X_h, eta)
        seam_jet_ideal = self.outer_base.profile_values(self.X_h, eta)
        at_XR = self.similarity_profile_values(self.X_R, eta)

        def scalar(value: Any) -> float:
            return float(np.asarray(value))

        value_jump = {
            "F": scalar(seam_ideal["F"]) - scalar(seam_current["F_current_joined"]),
            "U": scalar(seam_ideal["U"]) - scalar(seam_current["U_current_joined"]),
            "E": scalar(seam_ideal["E"]) - scalar(seam_current["E_current_joined"]),
        }
        jet_jump = {
            "F_X": scalar(seam_jet_ideal["F_X"]) - scalar(seam_jet_current["F_current_joined_X"]),
            "U_X": scalar(seam_jet_ideal["U_X"]) - scalar(seam_jet_current["U_current_joined_X"]),
            "E_X": scalar(seam_jet_ideal["E_X"]) - scalar(seam_jet_current["E_current_joined_X"]),
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
                "X_h": self.X_h,
                "X_R": self.X_R,
                "log_X_h": math.log(self.X_h),
                "log_X_R": math.log(self.X_R),
                "eta_interval": list(self.eta_interval),
                "time_interval": list(self.current.geometry.time_interval),
            },
            "seam_probe_eta": eta,
            "seam_value_jump_ideal_minus_current": value_jump,
            "seam_first_radial_jet_jump_ideal_minus_current": jet_jump,
            "XR_probe": {
                "F": scalar(at_XR["F_current_exterior"]),
                "U": scalar(at_XR["U_current_exterior"]),
                "E": scalar(at_XR["E_current_exterior"]),
                "M_over_X_minus_ideal": scalar(at_XR["M_over_X_current_exterior"]) - 4.0 * eta,
                "M_eta_over_X_minus_ideal": scalar(at_XR["M_eta_over_X_current_exterior"]) - 4.0,
                "v0": scalar(at_XR["v0_current_exterior"]),
            },
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="write deterministic scoped receipt JSON")
    parser.add_argument("--config-output", required=True, help="write deterministic configuration JSON")
    args = parser.parse_args()

    candidate = KokunoPA16CurrentCartesianExteriorToXR()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    _main()
