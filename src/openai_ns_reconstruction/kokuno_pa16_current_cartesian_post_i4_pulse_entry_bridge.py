"""Bridge the current Kokuno leading velocity from I4 exit to pulse entry.

This is one bounded Kokuno Agent-1 increment stacked exactly on A1 PR #1079
head ``b06742ca6e189499192ede3cce40f62cdc1e35ca``.

Pinned public provenance is KokunoYumeto's corrected 2026-09-09
reconstruction.  Its complete outer schedule holds

    ell = d_logX log(sqrt(2X) E) = -lambda

for ``T_w = 60 log(1/lambda)``.  The reserved I4 interval occupies
``(T_w-8,T_w-3)`` inside that hold, while the leading axial pulse starts only
at the hold endpoint ``X_p``.  The existing current RF40 power-law endpoint is

    X_4 = X_w exp(T_w),

so this module identifies the source pulse-entry seam with ``X_p = X_4`` and
closes only the remaining three-log-unit interval from I4 exit to that seam.
On this bridge the leading profile is still the same RF40 power law,

    U_0 = 0,
    ell = -lambda,
    F_0 = E_0 / sqrt(2X).

The exact #1079 current candidate may carry a nonzero repository-side
incompressibility-memory remainder from its earlier autonomous modulation and
repairs.  We measure that total-minus-RF40 remainder at the I4 exit and preserve
physical ``Delta M`` and ``Delta M_eta`` exactly.  Because ``U_0=0`` on the
bridge, their normalized values decay as ``X_I4_end / X``.  ``v0`` is rebuilt
from the full current primitive history.

The source pulse ``U=E R_b`` and its ``Amp(eta), c1, c2`` solve are deliberately
not materialized here.  The candidate fails closed immediately beyond pulse
entry.  Kokuno remains a public reconstruction/provenance source only; this is
not a paper-exact/OpenAI-exact field and is not independent PDE validation.
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
from .kokuno_pa16_current_cartesian_i4_leading_preservation import (
    KokunoPA16CurrentCartesianI4LeadingPreservation,
)


SCHEMA = "kokuno-pa16-current-cartesian-post-i4-pulse-entry-bridge-v1"
PARENT_EXACT_HEAD = "b06742ca6e189499192ede3cce40f62cdc1e35ca"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"

_SOURCE_FORMULAS = {
    "RF40_constant_lambda_hold": (
        "after the unit lambda turn, hold ell=-lambda for "
        "T_w=60 log(1/lambda)"
    ),
    "reserved_I4": "I4=(T_w-8,T_w-3) in log(X/X_w)",
    "pulse_entry": (
        "at X=X_p start the pulse after the T_w hold; therefore the current "
        "RF40 endpoint X_4=X_w exp(T_w) is the pulse-entry seam X_p"
    ),
    "post_I4_gap": (
        "I4 ends at T_w-3, so the pre-pulse bridge has exact logarithmic "
        "length 3"
    ),
    "pre_pulse_leading": "U_0=0 and ell=-lambda on the bridge",
    "power_law_E": (
        "E_0 continues the constant-lambda power law, proportional to "
        "f(eta) X^(-1/2-lambda)"
    ),
    "F": "F=E/sqrt(2X)",
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "bridge_memory": (
        "because U_0=0, physical Delta M and Delta M_eta stay constant; "
        "normalized remainders decay as X_I4_end/X"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
    "next_unmaterialized_source_stage": (
        "at pulse entry the source next uses U=E R_b with R_b depending on "
        "Amp(eta), c1, c2; that pulse is outside this increment"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact #1079 current leading-through-I4 candidate",
    "bridge_profile": (
        "reuse the exact current RF40 power-law F/U/E values and analytic "
        "first radial derivatives from I4 exit through X_4"
    ),
    "primitive_handoff": (
        "measure exact #1079 total-minus-current-RF40 M/X and M_eta/X at I4 "
        "exit and carry the associated physical remainder with exact 1/X scaling"
    ),
    "pulse_entry_identity": "X_p := current RF40 endpoint X_4",
    "new_tuning_parameters": "none",
    "domain_stop": "fail closed for X>X_p before evaluating the axial pulse",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "exact_current_I4_parent_consumed": True,
    "source_post_I4_pre_pulse_constant_lambda_bridge_identified": True,
    "source_post_I4_pre_pulse_bridge_log_length_is_three": True,
    "current_leading_preserved_from_I4_exit_to_pulse_entry": True,
    "current_incompressibility_memory_carried_to_pulse_entry": True,
    "current_cartesian_velocity_executable_to_pulse_entry": True,
    "velocity_interface_vectorized": True,
    "configuration_serializable": True,
    "source_positive_order_I3_profiles_materialized": False,
    "I3_positive_order_correction_materialized": False,
    "source_I4_mean_correction_materialized": False,
    "current_I4_mean_correction_materialized": False,
    "source_axial_pulse_materialized": False,
    "source_pulse_amplitude_root_materialized": False,
    "source_pulse_M_J_end_correction_materialized": False,
    "source_terminal_tail_schedule_bound_into_current_velocity": False,
    "source_exterior_heat_replacement_materialized": False,
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
class KokunoPA16CurrentCartesianPostI4PulseEntryBridge:
    """Current leading Cartesian candidate through the pre-pulse entry seam."""

    parent: KokunoPA16CurrentCartesianI4LeadingPreservation = field(
        default_factory=KokunoPA16CurrentCartesianI4LeadingPreservation,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianI4LeadingPreservation):
            raise TypeError(
                "parent must be KokunoPA16CurrentCartesianI4LeadingPreservation"
            )
        if not math.isclose(
            self.log_X_p - self.parent.log_X_I4_end,
            3.0,
            rel_tol=0.0,
            abs_tol=4.0e-12,
        ):
            raise ValueError(
                "corrected schedule requires a three-log-unit I4-exit to pulse-entry gap"
            )
        if not math.isclose(
            self.log_X_p,
            self.current.log_X_4,
            rel_tol=0.0,
            abs_tol=2.0e-12,
        ):
            raise ValueError("pulse-entry seam must equal current RF40 endpoint X_4")
        if not (0.0 < self.parent.X_I4_end < self.X_p < np.finfo(float).max):
            raise ValueError("post-I4 pulse-entry bridge is not materializable")

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
    def log_X_I4_end(self) -> float:
        return float(self.parent.log_X_I4_end)

    @property
    def X_I4_end(self) -> float:
        return float(self.parent.X_I4_end)

    @property
    def log_X_p(self) -> float:
        # The corrected schedule starts the pulse at the end of the RF40
        # constant-lambda hold.  In the current executable base that endpoint
        # is named X_4.
        return float(self.current.log_X_4)

    @property
    def X_p(self) -> float:
        return float(self.current.X_4)

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.parent.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_p)):
            raise ValueError(
                "X must lie in the current post-I4-to-pulse-entry domain "
                f"[0,{self.X_p:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _i4_exit_overlay_memory(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        X = np.full(unique.shape, self.X_I4_end, dtype=float)
        total = self.parent.similarity_profile_values(X, unique)
        base = self.current.similarity_profile_values(X, unique)
        delta = np.asarray(
            total["M_over_X_current_leading_through_I4"], dtype=float
        ) - np.asarray(base["M_over_X_current_rf40_power_law"], dtype=float)
        delta_eta = np.asarray(
            total["M_eta_over_X_current_leading_through_I4"], dtype=float
        ) - np.asarray(base["M_eta_over_X_current_rf40_power_law"], dtype=float)
        return (
            delta[inverse].reshape(values.shape),
            delta_eta[inverse].reshape(values.shape),
        )

    def _carried_post_i4_memory(
        self, X: np.ndarray, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        dm, dm_eta = self._i4_exit_overlay_memory(eta)
        factor = self.X_I4_end / np.asarray(X, dtype=float)
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

        inherited = xf <= self.X_I4_end
        if np.any(inherited):
            p = self.parent.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(p["F_current_leading_through_I4"], dtype=float)
            U[inherited] = np.asarray(p["U_current_leading_through_I4"], dtype=float)
            E[inherited] = np.asarray(p["E_current_leading_through_I4"], dtype=float)
            m_ratio[inherited] = np.asarray(
                p["M_over_X_current_leading_through_I4"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_leading_through_I4"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_leading_through_I4"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        bridge = ~inherited
        if np.any(bridge):
            Xo = xf[bridge]
            etao = ef[bridge]
            base = self.current.similarity_profile_values(Xo, etao)
            dm, dm_eta = self._carried_post_i4_memory(Xo, etao)
            mo = np.asarray(base["M_over_X_current_rf40_power_law"], dtype=float) + dm
            meo = np.asarray(base["M_eta_over_X_current_rf40_power_law"], dtype=float) + dm_eta
            Uo = np.asarray(base["U_current_rf40_power_law"], dtype=float)
            if np.any(Uo != 0.0):
                raise RuntimeError("leading U_0 must remain zero before pulse entry")

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

            F[bridge] = np.asarray(base["F_current_rf40_power_law"], dtype=float)
            U[bridge] = Uo
            E[bridge] = np.asarray(base["E_current_rf40_power_law"], dtype=float)
            m_ratio[bridge] = mo
            m_eta_ratio[bridge] = meo
            v0[bridge] = v0o
            carried_m[bridge] = dm
            carried_m_eta[bridge] = dm_eta
            region[bridge] = "post_I4_pre_pulse_constant_lambda_bridge"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("post-I4 pulse-entry bridge produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("post-I4 pulse-entry bridge lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_leading_to_pulse_entry": F.reshape(shape),
            "U_current_leading_to_pulse_entry": U.reshape(shape),
            "E_current_leading_to_pulse_entry": E.reshape(shape),
            "M_over_X_current_leading_to_pulse_entry": m_ratio.reshape(shape),
            "M_eta_over_X_current_leading_to_pulse_entry": m_eta_ratio.reshape(shape),
            "v0_current_leading_to_pulse_entry": v0.reshape(shape),
            "delta_M_over_X_carried_from_I4_exit": carried_m.reshape(shape),
            "delta_M_eta_over_X_carried_from_I4_exit": carried_m_eta.reshape(shape),
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

        inherited = xf <= self.X_I4_end
        if np.any(inherited):
            p = self.parent.similarity_radial_derivatives(xf[inherited], ef[inherited])
            F_X[inherited] = np.asarray(
                p["F_current_leading_through_I4_X"], dtype=float
            )
            U_X[inherited] = np.asarray(
                p["U_current_leading_through_I4_X"], dtype=float
            )
            E_X[inherited] = np.asarray(
                p["E_current_leading_through_I4_X"], dtype=float
            )

        bridge = ~inherited
        if np.any(bridge):
            p = self.current.similarity_radial_derivatives(xf[bridge], ef[bridge])
            F_X[bridge] = np.asarray(
                p["F_current_rf40_power_law_X"], dtype=float
            )
            U_X[bridge] = np.asarray(
                p["U_current_rf40_power_law_X"], dtype=float
            )
            E_X[bridge] = np.asarray(
                p["E_current_rf40_power_law_X"], dtype=float
            )

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_current_leading_to_pulse_entry_X": F_X.reshape(shape),
            "U_current_leading_to_pulse_entry_X": U_X.reshape(shape),
            "E_current_leading_to_pulse_entry_X": E_X.reshape(shape),
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        tol = 128.0 * np.finfo(float).eps * self.X_p
        if np.any(X > self.X_p + tol):
            raise ValueError("Cartesian point lies beyond the current pulse-entry domain")
        X = np.minimum(X, self.X_p)
        eta = np.asarray(coords["eta"], dtype=float)
        p = self.similarity_profile_values(X, eta)

        xb, yb, qb = np.broadcast_arrays(
            _finite(x, "x"), _finite(y, "y"), np.asarray(coords["q"], dtype=float)
        )
        F = np.asarray(p["F_current_leading_to_pulse_entry"], dtype=float)
        U = np.asarray(p["U_current_leading_to_pulse_entry"], dtype=float)
        v0 = np.asarray(p["v0_current_leading_to_pulse_entry"], dtype=float)
        radial = v0 / (2.0 * qb)
        swirl = qb ** (-self.A - 0.5) * F
        axial = qb ** (-self.A) * U
        return np.stack(
            (radial * xb - swirl * yb, radial * yb + swirl * xb, axial), axis=-1
        )

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_i4": self.parent.configuration(),
            "bound_source_scope": {
                "source_commit": SOURCE_COMMIT,
                "source_blob": SOURCE_BLOB,
                "bridge_role": "post_I4_constant_lambda_to_pulse_entry",
                "pulse_entry_identity": "X_p_equals_current_X_4",
                "source_axial_pulse_materialized": False,
                "mutable": False,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianPostI4PulseEntryBridge":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current post-I4 pulse-entry bridge schema")
        parent_payload = payload.get("parent_i4")
        bound = payload.get("bound_source_scope")
        if not isinstance(parent_payload, Mapping) or not isinstance(bound, Mapping):
            raise ValueError("missing current I4 parent or post-I4 source-scope binding")
        expected = {
            "source_commit": SOURCE_COMMIT,
            "source_blob": SOURCE_BLOB,
            "bridge_role": "post_I4_constant_lambda_to_pulse_entry",
            "pulse_entry_identity": "X_p_equals_current_X_4",
            "source_axial_pulse_materialized": False,
            "mutable": False,
        }
        if dict(bound) != expected:
            raise ValueError(
                "serialized post-I4 source-scope binding differs from the frozen contract"
            )
        parent = KokunoPA16CurrentCartesianI4LeadingPreservation.from_configuration(
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
    ) -> "KokunoPA16CurrentCartesianPostI4PulseEntryBridge":
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
                "log_X_I4_end": self.log_X_I4_end,
                "log_X_p": self.log_X_p,
                "post_I4_gap_log_length": self.log_X_p - self.log_X_I4_end,
            },
            "truth_boundary": _TRUTH_BOUNDARY,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta_probe = np.asarray(0.37)
        seam = self.similarity_profile_values(np.asarray(self.X_I4_end), eta_probe)
        entry = self.similarity_profile_values(np.asarray(self.X_p), eta_probe)
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
                "log_X_I4_end": self.log_X_I4_end,
                "log_X_p": self.log_X_p,
                "post_I4_gap_log_length": self.log_X_p - self.log_X_I4_end,
            },
            "eta_0p37_receipt": {
                "M_over_X_at_I4_exit": float(
                    np.asarray(seam["M_over_X_current_leading_to_pulse_entry"])
                ),
                "M_over_X_at_pulse_entry": float(
                    np.asarray(entry["M_over_X_current_leading_to_pulse_entry"])
                ),
                "U_at_pulse_entry": float(
                    np.asarray(entry["U_current_leading_to_pulse_entry"])
                ),
                "region_at_pulse_entry": str(np.asarray(entry["region"]).item()),
            },
            "truth_boundary": self.truth_boundary,
            "note": (
                "This receipt closes only the source-specified constant-lambda bridge from "
                "I4 exit to pulse entry. The axial pulse itself is not materialized and "
                "this report is not independent PDE validation."
            ),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()
    candidate = KokunoPA16CurrentCartesianPostI4PulseEntryBridge()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
