"""Preserve the current A1 *leading* velocity through reserved interval I4.

This is one narrow Kokuno Agent-1 increment stacked exactly on A1 PR #1072
head ``5fb7b583062b4a991db86e39fdcdb231022b70c9``.

The pinned corrected 2026-09-09 public reconstruction assigns

    I_pos = I3,   I_mean = I4,

and, crucially for Agent-1 scope, states that throughout the RF40 power-law
stage including I4 the leading profiles remain

    U_0 = 0,
    E_0 = c_patch (1+eta^2)^(-1) X^(-1/2-lambda).

It further states that the subsequent modifications preserve this exact leading
statement on I4: cone repair is in I1, heat compensation is in I2, and later
pulse/exterior-heat edits start after I4.  Therefore I4's role as the reserved
*mean-correction* interval must not be misread as a leading E_0/U_0 patch.

Accordingly this adapter does not invent an I4 leading overlay.  It consumes
#1072 unchanged through I3, preserves the current RF40 leading F/U/E through
the gap and I4, and carries the actual current-lineage incompressibility-memory
remainder.  Since leading U_0=0 after I3, the physical Delta M and Delta M_eta
remain constant, so their normalized ratios decay exactly as 1/X.  v0 is then
rebuilt from the complete current primitive history.

The public reconstruction is provenance/formula source only.  This module does
not materialize the source mean-correction construction itself, does not bind
the later terminal/exterior schedule, and does not claim paper/OpenAI exactness
or PDE validation.
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
from .kokuno_pa16_current_cartesian_i3_leading_preservation import (
    KokunoPA16CurrentCartesianI3LeadingPreservation,
)

SCHEMA = "kokuno-pa16-current-cartesian-i4-leading-preservation-v1"
PARENT_EXACT_HEAD = "5fb7b583062b4a991db86e39fdcdb231022b70c9"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"

_SOURCE_FORMULAS = {
    "reserved_intervals": (
        "I1=X_w exp((T_w-25,T_w-20)); I2=X_w exp((T_w-20,T_w-15)); "
        "I3=X_w exp((T_w-14,T_w-9)); I4=X_w exp((T_w-8,T_w-3))"
    ),
    "assignment": "I_pos=I3; I_mean=I4",
    "I4_scope": "I4 is the reserved mean-correction interval",
    "leading_on_I4": (
        "throughout the RF40 power-law stage including I4, "
        "U_0=0 and E_0=c_patch (1+eta^2)^(-1) X^(-1/2-lambda)"
    ),
    "leading_preservation_statement": (
        "subsequent modifications preserve that exact leading statement on I4; "
        "none changes E_0,U_0 on I4"
    ),
    "F": "F=E/sqrt(2X)",
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "post_I3_leading_memory": (
        "leading U_0=0, so the current physical overlay Delta M and Delta M_eta "
        "are constant and Delta(M/X),Delta(M_eta/X) decay as X_I3_end/X"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact #1072 current leading-through-I3 candidate",
    "I4_leading_action": (
        "no leading overlay is applied: reuse the current RF40 power-law F/U/E "
        "through I4 because the public I4 role is mean correction while E_0,U_0 "
        "remain unchanged there"
    ),
    "primitive_handoff": (
        "measure total-minus-RF40 M/X and M_eta/X at exact current I3 exit and "
        "carry the associated physical remainder with exact 1/X scaling"
    ),
    "derivatives": "reuse current RF40 analytic first radial F/U/E derivatives",
    "coordinates": "reuse governed current q/X/eta inverse and Cartesian map",
    "new_tuning_parameters": "none",
    "domain_stop": "fail closed at the I4 upper endpoint before later terminal edits",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "exact_current_I3_parent_consumed": True,
    "source_I4_is_reserved_mean_correction_interval": True,
    "source_leading_E0_U0_unchanged_on_I4": True,
    "source_I4_is_not_a_leading_E0_U0_patch": True,
    "current_leading_preserved_through_I4": True,
    "current_incompressibility_memory_carried_through_I4": True,
    "cartesian_velocity_executable_through_I4": True,
    "velocity_interface_vectorized": True,
    "configuration_serializable": True,
    "source_positive_order_I3_profiles_materialized": False,
    "I3_positive_order_correction_materialized": False,
    "source_I4_mean_correction_materialized": False,
    "current_I4_mean_correction_materialized": False,
    "leading_I4_overlay_invented": False,
    "source_terminal_tail_schedule_bound_into_current_velocity": False,
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
class KokunoPA16CurrentCartesianI4LeadingPreservation:
    """Current leading Cartesian candidate through the reserved I4 interval."""

    parent: KokunoPA16CurrentCartesianI3LeadingPreservation = field(
        default_factory=KokunoPA16CurrentCartesianI3LeadingPreservation,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianI3LeadingPreservation):
            raise TypeError(
                "parent must be KokunoPA16CurrentCartesianI3LeadingPreservation"
            )
        intervals = self.outer_schedule.reserved_log_intervals()
        i3_low, i3_high = intervals["I3"]
        i4_low, i4_high = intervals["I4"]
        if not math.isclose(
            float(i3_high), self.parent.log_X_I3_end, rel_tol=0.0, abs_tol=2.0e-12
        ):
            raise ValueError("parent I3 endpoint and reserved schedule disagree")
        if not float(i3_low) < float(i3_high) < float(i4_low) < float(i4_high):
            raise ValueError("reserved I3/I4 ordering is inconsistent")
        if float(i4_high) >= math.log(np.finfo(float).max):
            raise ValueError("current I4 endpoint is not materializable in float64 X")
        if self.X_I4_end > self.current.X_4:
            raise ValueError("reserved I4 must remain inside the current RF40 power-law domain")

    @property
    def current(self):
        return self.parent.current

    @property
    def outer_schedule(self):
        return self.parent.outer_schedule

    @property
    def geometry(self):
        return self.parent.geometry

    @property
    def A(self) -> float:
        return float(self.parent.A)

    @property
    def D(self) -> float:
        return float(self.parent.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.parent.eta_interval)

    @property
    def log_X_I3_end(self) -> float:
        return float(self.parent.log_X_I3_end)

    @property
    def X_I3_end(self) -> float:
        return float(self.parent.X_I3_end)

    @property
    def log_X_I4_start(self) -> float:
        return float(self.outer_schedule.reserved_log_intervals()["I4"][0])

    @property
    def log_X_I4_end(self) -> float:
        return float(self.outer_schedule.reserved_log_intervals()["I4"][1])

    @property
    def X_I4_start(self) -> float:
        return float(math.exp(self.log_X_I4_start))

    @property
    def X_I4_end(self) -> float:
        return float(math.exp(self.log_X_I4_end))

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.parent.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_I4_end)):
            raise ValueError(
                "X must lie in the current leading-through-I4 domain "
                f"[0,{self.X_I4_end:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _i3_exit_overlay_memory(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        X = np.full(unique.shape, self.X_I3_end, dtype=float)
        total = self.parent.similarity_profile_values(X, unique)
        base = self.current.similarity_profile_values(X, unique)
        delta = np.asarray(
            total["M_over_X_current_leading_through_I3"], dtype=float
        ) - np.asarray(base["M_over_X_current_rf40_power_law"], dtype=float)
        delta_eta = np.asarray(
            total["M_eta_over_X_current_leading_through_I3"], dtype=float
        ) - np.asarray(base["M_eta_over_X_current_rf40_power_law"], dtype=float)
        return (
            delta[inverse].reshape(values.shape),
            delta_eta[inverse].reshape(values.shape),
        )

    def _carried_post_i3_memory(
        self, X: np.ndarray, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        dm, dm_eta = self._i3_exit_overlay_memory(eta)
        factor = self.X_I3_end / np.asarray(X, dtype=float)
        return factor * dm, factor * dm_eta

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
        carried_m = np.zeros_like(xf)
        carried_m_eta = np.zeros_like(xf)
        region = np.empty(xf.shape, dtype=object)

        inherited = xf <= self.X_I3_end
        if np.any(inherited):
            p = self.parent.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(p["F_current_leading_through_I3"], dtype=float)
            U[inherited] = np.asarray(p["U_current_leading_through_I3"], dtype=float)
            E[inherited] = np.asarray(p["E_current_leading_through_I3"], dtype=float)
            m_ratio[inherited] = np.asarray(
                p["M_over_X_current_leading_through_I3"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_leading_through_I3"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_leading_through_I3"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        post_i3 = ~inherited
        if np.any(post_i3):
            Xo = xf[post_i3]
            etao = ef[post_i3]
            base = self.current.similarity_profile_values(Xo, etao)
            dm, dm_eta = self._carried_post_i3_memory(Xo, etao)
            mo = np.asarray(base["M_over_X_current_rf40_power_law"], dtype=float) + dm
            meo = np.asarray(base["M_eta_over_X_current_rf40_power_law"], dtype=float) + dm_eta
            Uo = np.asarray(base["U_current_rf40_power_law"], dtype=float)
            if np.any(Uo != 0.0):
                raise RuntimeError("leading U_0 must remain zero from I3 exit through I4")

            axis = self.geometry.physical_profiles.axis_profiles.values(
                np.zeros_like(etao), etao
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            v0o = (
                2.0 * etao * Uo
                - 2.0 * self.D * etao * mo
                - d * meo
            ) / L

            F[post_i3] = np.asarray(base["F_current_rf40_power_law"], dtype=float)
            U[post_i3] = Uo
            E[post_i3] = np.asarray(base["E_current_rf40_power_law"], dtype=float)
            m_ratio[post_i3] = mo
            m_eta_ratio[post_i3] = meo
            v0[post_i3] = v0o
            carried_m[post_i3] = dm
            carried_m_eta[post_i3] = dm_eta
            in_i4 = np.log(Xo) >= self.log_X_I4_start
            labels = np.where(
                in_i4,
                "reserved_I4_leading_preserved",
                "post_I3_pre_I4_gap",
            )
            region[post_i3] = labels

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current leading-through-I4 profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current leading-through-I4 profile lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_leading_through_I4": F.reshape(shape),
            "U_current_leading_through_I4": U.reshape(shape),
            "E_current_leading_through_I4": E.reshape(shape),
            "M_over_X_current_leading_through_I4": m_ratio.reshape(shape),
            "M_eta_over_X_current_leading_through_I4": m_eta_ratio.reshape(shape),
            "v0_current_leading_through_I4": v0.reshape(shape),
            "delta_M_over_X_carried_from_I3_exit": carried_m.reshape(shape),
            "delta_M_eta_over_X_carried_from_I3_exit": carried_m_eta.reshape(shape),
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

        inherited = xf <= self.X_I3_end
        if np.any(inherited):
            p = self.parent.similarity_radial_derivatives(xf[inherited], ef[inherited])
            F_X[inherited] = np.asarray(
                p["F_current_leading_through_I3_X"], dtype=float
            )
            U_X[inherited] = np.asarray(
                p["U_current_leading_through_I3_X"], dtype=float
            )
            E_X[inherited] = np.asarray(
                p["E_current_leading_through_I3_X"], dtype=float
            )

        post_i3 = ~inherited
        if np.any(post_i3):
            p = self.current.similarity_radial_derivatives(xf[post_i3], ef[post_i3])
            F_X[post_i3] = np.asarray(
                p["F_current_rf40_power_law_X"], dtype=float
            )
            U_X[post_i3] = np.asarray(
                p["U_current_rf40_power_law_X"], dtype=float
            )
            E_X[post_i3] = np.asarray(
                p["E_current_rf40_power_law_X"], dtype=float
            )

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_current_leading_through_I4_X": F_X.reshape(shape),
            "U_current_leading_through_I4_X": U_X.reshape(shape),
            "E_current_leading_through_I4_X": E_X.reshape(shape),
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        tol = 128.0 * np.finfo(float).eps * self.X_I4_end
        if np.any(X > self.X_I4_end + tol):
            raise ValueError("Cartesian point lies beyond the current leading-through-I4 domain")
        X = np.minimum(X, self.X_I4_end)
        eta = np.asarray(coords["eta"], dtype=float)
        p = self.similarity_profile_values(X, eta)

        xb, yb, qb = np.broadcast_arrays(
            _finite(x, "x"), _finite(y, "y"), np.asarray(coords["q"], dtype=float)
        )
        F = np.asarray(p["F_current_leading_through_I4"], dtype=float)
        U = np.asarray(p["U_current_leading_through_I4"], dtype=float)
        v0 = np.asarray(p["v0_current_leading_through_I4"], dtype=float)
        radial = v0 / (2.0 * qb)
        swirl = qb ** (-self.A - 0.5) * F
        axial = qb ** (-self.A) * U
        return np.stack(
            (radial * xb - swirl * yb, radial * yb + swirl * xb, axial), axis=-1
        )

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_i3": self.parent.configuration(),
            "bound_source_scope": {
                "source_commit": SOURCE_COMMIT,
                "source_blob": SOURCE_BLOB,
                "I4_role": "mean_correction_not_leading_E0_U0",
                "leading_action": "preserve_RF40_power_law",
                "mutable": False,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianI4LeadingPreservation":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current leading-through-I4 schema")
        parent_payload = payload.get("parent_i3")
        bound = payload.get("bound_source_scope")
        if not isinstance(parent_payload, Mapping) or not isinstance(bound, Mapping):
            raise ValueError("missing current I3 parent or I4 source-scope binding")
        expected = {
            "source_commit": SOURCE_COMMIT,
            "source_blob": SOURCE_BLOB,
            "I4_role": "mean_correction_not_leading_E0_U0",
            "leading_action": "preserve_RF40_power_law",
            "mutable": False,
        }
        if dict(bound) != expected:
            raise ValueError("serialized I4 source-scope binding differs from the frozen contract")
        parent = KokunoPA16CurrentCartesianI3LeadingPreservation.from_configuration(
            parent_payload
        )
        return cls(parent=parent)

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
    ) -> "KokunoPA16CurrentCartesianI4LeadingPreservation":
        return cls.from_configuration(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "blob": SOURCE_BLOB,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "domain": {
                "log_X_I3_end": self.log_X_I3_end,
                "log_X_I4_start": self.log_X_I4_start,
                "log_X_I4_end": self.log_X_I4_end,
            },
            "truth_boundary": _TRUTH_BOUNDARY,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta_probe = np.asarray(0.37)
        seam = self.similarity_profile_values(
            np.asarray(self.X_I3_end), eta_probe
        )
        end = self.similarity_profile_values(
            np.asarray(self.X_I4_end), eta_probe
        )
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "semantic_sha256": self.semantic_sha256,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "blob": SOURCE_BLOB,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formula_mapping": copy.deepcopy(_SOURCE_FORMULAS),
            "numerical_realization": copy.deepcopy(_NUMERICAL_REALIZATION),
            "domain": {
                "log_X_I3_end": self.log_X_I3_end,
                "log_X_I4_start": self.log_X_I4_start,
                "log_X_I4_end": self.log_X_I4_end,
            },
            "eta_0p37_receipt": {
                "M_over_X_at_I3_exit": float(
                    np.asarray(seam["M_over_X_current_leading_through_I4"])
                ),
                "M_over_X_at_I4_exit": float(
                    np.asarray(end["M_over_X_current_leading_through_I4"])
                ),
                "U_at_I4_exit": float(
                    np.asarray(end["U_current_leading_through_I4"])
                ),
                "region_at_I4_exit": str(np.asarray(end["region"]).item()),
            },
            "truth_boundary": self.truth_boundary,
            "note": (
                "I4 preservation is a leading-profile scope statement. The source mean-correction "
                "construction itself is not materialized here and this report is not PDE validation."
            ),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()
    candidate = KokunoPA16CurrentCartesianI4LeadingPreservation()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
