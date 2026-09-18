"""Source-scheduled multi-band Kokuno localized real-pair velocity family.

This module closes one specific gap left by the earlier supplied-data family:
Kokuno's corrected 2026-09-09 reconstruction assigns every dyadic band its own

    Q_ell = 2^(-ell),      epsilon_ell = Q_ell^h.

A single shared ``epsilon`` is therefore not a source-compatible realization of
an actually multi-band family.  Here the whole slow squared partition is first
validated, then every beta=(ell,a) label is curled with its *own* source-derived
``epsilon_ell`` while retaining the complete support-gradient curl terms.  The
real m=+/-1 pair is reconstructed by conjugation, each label receives the
source physical factor ``Q_ell^(-A)``, A=1/2+h, and the labels are summed only
after those per-band operations.

The caller still supplies Phi, n_Phi, t_plus, D_r C_plus, D_z C_plus and the
partition values/derivatives.  Actual positive-order/background data and the
actual auxiliary-torus mode labels are not recovered here.  Thus the output is
a source-scheduled supplied-mode correction family, not a paper-exact field or
a public velocity(x,y,z,t) map.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates
from .kokuno_source_band_covering import KokunoSourceBandCovering
from .kokuno_source_partition_family import KokunoSourcePartitionFamilyContract
from .kokuno_source_support_localized_curl import KokunoSourceSupportLocalizedCurl

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-multiband-real-pair-family-v1"

_SOURCE_FORMULAS = {
    "band_scales": "Q_ell=2^(-ell); epsilon_ell=Q_ell^h",
    "slow_partition": "eta_beta=chi_ell*chi_{ell,a}; sum_beta eta_beta^2=1",
    "localized_potential": "A_{beta,m}=eta_beta*C_{beta,m}*exp(i*k_{ell,m}*Phi_beta)",
    "band_wavenumber": "k_{ell,m}=m/epsilon_ell",
    "real_harmonics": "u_{beta,-}=conjugate(u_{beta,+}); u_beta=2 Re(u_{beta,+})",
    "physical_scaling": "u_phys,beta=Q_ell^(-A) u_beta, A=1/2+h",
}

_TRUTH_BOUNDARY = {
    "source_per_band_Q_and_epsilon_schedule_executable": True,
    "whole_slow_partition_guard_reused": True,
    "per_band_complete_curl_executable": True,
    "per_band_real_conjugate_pair_executable": True,
    "per_band_Q_physical_scaling_executable": True,
    "distinct_dyadic_band_columns_materialized_from_supplied_modes": True,
    "bounded_two_coordinate_handoff_available": True,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_actual_partition_labels_instantiated": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "genuinely_independent_second_covariance_column_ready": False,
    "agent3_rank_screen_still_required": True,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_real(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def _real_vector3(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim < 2 or out.shape[-1] != 3 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must have finite shape sample_shape+(beta,3)")
    return out


def _complex_vector3(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=np.complex128)
    if out.ndim < 2 or out.shape[-1] != 3:
        raise ValueError(f"{name} must have shape sample_shape+(beta,3)")
    if not np.all(np.isfinite(out.real)) or not np.all(np.isfinite(out.imag)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _rotate_cylindrical_by_beta(v: np.ndarray, theta: np.ndarray) -> np.ndarray:
    c = np.cos(theta)[..., None]
    s = np.sin(theta)[..., None]
    return np.stack(
        (
            v[..., 0] * c - v[..., 1] * s,
            v[..., 0] * s + v[..., 1] * c,
            v[..., 2],
        ),
        axis=-1,
    )


@dataclass(frozen=True)
class KokunoSourceMultiBandRealPairFamily:
    """Evaluate >=2 interacting dyadic bands with source-derived fast scales."""

    h: float = 0.005
    partition_atol: float = 2.0e-12
    derivative_rtol: float = 2.0e-10
    reality_atol: float = 2.0e-12
    coefficient_max_l1_update: float = 0.125

    def __post_init__(self) -> None:
        h = float(self.h)
        if not np.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
        for name, value, upper in (
            ("partition_atol", self.partition_atol, 1.0e-8),
            ("derivative_rtol", self.derivative_rtol, 1.0e-6),
            ("reality_atol", self.reality_atol, 1.0e-8),
        ):
            val = float(value)
            if not np.isfinite(val) or not (0.0 < val <= upper):
                raise ValueError(f"{name} must lie in (0,{upper:g}]")
            object.__setattr__(self, name, val)
        object.__setattr__(self, "h", h)
        # Reuse #400's bounded-unit guard at construction time.
        KokunoBoundedFamilyCoefficientCoordinates(
            max_l1_update=self.coefficient_max_l1_update
        )

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def partition_contract(self) -> KokunoSourcePartitionFamilyContract:
        return KokunoSourcePartitionFamilyContract(
            partition_atol=self.partition_atol,
            derivative_rtol=self.derivative_rtol,
        )

    def _validated_inputs(
        self,
        phase: Any,
        n_phi: Any,
        t_plus: Any,
        D_r_C_plus: Any,
        D_z_C_plus: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        beta_labels: Sequence[Any],
    ) -> dict[str, Any]:
        family = self.partition_contract.validate_family(
            eta, D_r_eta, D_z_eta, beta_labels
        )
        eta_array = family["eta"]
        phase_array = _finite_real(phase, "phase")
        n_array = _real_vector3(n_phi, "n_phi")
        t_array = _complex_vector3(t_plus, "t_plus")
        Dr_array = _complex_vector3(D_r_C_plus, "D_r_C_plus")
        Dz_array = _complex_vector3(D_z_C_plus, "D_z_C_plus")
        if phase_array.shape != eta_array.shape:
            raise ValueError("phase must have exactly the eta shape, including beta")
        expected_vector = eta_array.shape + (3,)
        for name, value in (
            ("n_phi", n_array),
            ("t_plus", t_array),
            ("D_r_C_plus", Dr_array),
            ("D_z_C_plus", Dz_array),
        ):
            if value.shape != expected_vector:
                raise ValueError(f"{name} must have shape eta.shape+(3,)")

        labels = family["beta_labels"]
        ell = tuple(int(label[0]) for label in labels)
        active = tuple(sorted(set(ell)))
        if len(active) < 2:
            raise ValueError(
                "multi-band handoff requires at least two distinct ell bands; a same-band phase copy is not a second band"
            )
        if active[-1] - active[0] > 4:
            raise ValueError("active interacting band window must satisfy max(ell)-min(ell)<=4")
        schedules = tuple(KokunoSourceBandCovering(e, self.h) for e in ell)
        return {
            **family,
            "phase": phase_array,
            "n_phi": n_array,
            "t_plus": t_array,
            "D_r_C_plus": Dr_array,
            "D_z_C_plus": Dz_array,
            "ell_by_beta": ell,
            "active_ell_bands": active,
            "band_schedules": schedules,
        }

    def physical_family(
        self,
        R: Any,
        theta: Any,
        phase: Any,
        n_phi: Any,
        t_plus: Any,
        D_r_C_plus: Any,
        D_z_C_plus: Any,
        eta: Any,
        D_r_eta: Any,
        D_z_eta: Any,
        beta_labels: Sequence[Any],
    ) -> dict[str, Any]:
        """Return a real Q-scaled supplied-mode family with per-band epsilon."""
        data = self._validated_inputs(
            phase,
            n_phi,
            t_plus,
            D_r_C_plus,
            D_z_C_plus,
            eta,
            D_r_eta,
            D_z_eta,
            beta_labels,
        )
        family_shape = data["eta"].shape
        sample_shape = family_shape[:-1]
        n_beta = family_shape[-1]

        R_array = _finite_real(R, "R")
        try:
            R_family = np.broadcast_to(R_array, family_shape)
        except ValueError:
            if R_array.shape == sample_shape:
                R_family = np.broadcast_to(R_array[..., None], family_shape)
            else:
                raise ValueError("R must broadcast to the eta family shape") from None
        if np.any(R_family <= 0.0):
            raise ValueError("R must stay strictly away from the cylindrical axis")

        theta_array = _finite_real(theta, "theta")
        try:
            theta_sample = np.broadcast_to(theta_array, sample_shape)
        except ValueError:
            raise ValueError("theta must broadcast to sample axes only") from None

        by_beta: list[np.ndarray] = []
        q_label: list[float] = []
        eps_label: list[float] = []
        covering_level: list[int] = []
        plus_by_beta: list[np.ndarray] = []
        minus_by_beta: list[np.ndarray] = []

        for j, schedule in enumerate(data["band_schedules"]):
            plus = KokunoSourceSupportLocalizedCurl(
                epsilon=schedule.epsilon, m=1, partition_atol=self.partition_atol
            ).localized_mode(
                R_family[..., j],
                data["phase"][..., j],
                data["n_phi"][..., j, :],
                data["t_plus"][..., j, :],
                data["D_r_C_plus"][..., j, :],
                data["D_z_C_plus"][..., j, :],
                data["eta"][..., j],
                data["D_r_eta"][..., j],
                data["D_z_eta"][..., j],
            )
            minus = KokunoSourceSupportLocalizedCurl(
                epsilon=schedule.epsilon, m=-1, partition_atol=self.partition_atol
            ).localized_mode(
                R_family[..., j],
                data["phase"][..., j],
                data["n_phi"][..., j, :],
                np.conjugate(data["t_plus"][..., j, :]),
                np.conjugate(data["D_r_C_plus"][..., j, :]),
                np.conjugate(data["D_z_C_plus"][..., j, :]),
                data["eta"][..., j],
                data["D_r_eta"][..., j],
                data["D_z_eta"][..., j],
            )
            pair_complex = plus["velocity"] + minus["velocity"]
            scale = max(1.0, float(np.max(np.abs(pair_complex.real), initial=0.0)))
            if float(np.max(np.abs(pair_complex.imag), initial=0.0)) > self.reality_atol * scale:
                raise RuntimeError("per-band conjugate pair failed the declared reality guard")
            pair = pair_complex.real
            if not np.allclose(
                pair,
                2.0 * plus["velocity"].real,
                rtol=0.0,
                atol=self.reality_atol * scale,
            ):
                raise RuntimeError("per-band m=+/-1 pair disagrees with 2*Re(m=+1)")
            physical = (schedule.Q ** (-self.A)) * pair
            by_beta.append(physical)
            plus_by_beta.append(plus["velocity"])
            minus_by_beta.append(minus["velocity"])
            q_label.append(schedule.Q)
            eps_label.append(schedule.epsilon)
            covering_level.append(schedule.covering_level)

        physical_by_beta = np.stack(by_beta, axis=-2)
        plus_velocity = np.stack(plus_by_beta, axis=-2)
        minus_velocity = np.stack(minus_by_beta, axis=-2)
        cyl_total = np.sum(physical_by_beta, axis=-2)
        cart_by_beta = _rotate_cylindrical_by_beta(physical_by_beta, theta_sample)
        cart_total = np.sum(cart_by_beta, axis=-2)

        band_columns = []
        band_column_norms = []
        for ell in data["active_ell_bands"]:
            indices = [j for j, value in enumerate(data["ell_by_beta"]) if value == ell]
            column = np.sum(np.take(cart_by_beta, indices, axis=-2), axis=-2)
            if not np.any(column != 0.0):
                raise ValueError(f"active ell={ell} collapses to an exactly zero physical band column")
            band_columns.append(column)
            band_column_norms.append(float(np.linalg.norm(column.reshape(-1, 3))))
        cart_by_band = np.stack(band_columns, axis=-2)
        if not np.allclose(np.sum(cart_by_band, axis=-2), cart_total, rtol=0.0, atol=2.0e-12):
            raise RuntimeError("band aggregation does not reproduce the total Cartesian correction")

        q_source = np.asarray(q_label, dtype=float)
        eps_source = np.asarray(eps_label, dtype=float)
        q_by_beta = np.broadcast_to(q_source, family_shape)
        eps_by_beta = np.broadcast_to(eps_source, family_shape)

        return {
            "beta_labels": data["beta_labels"],
            "ell_by_beta": data["ell_by_beta"],
            "active_ell_bands": data["active_ell_bands"],
            "Q_source_by_label": q_source,
            "epsilon_source_by_label": eps_source,
            "covering_level_by_label": np.asarray(covering_level, dtype=int),
            "Q_by_beta": q_by_beta,
            "epsilon_by_beta": eps_by_beta,
            "A": self.A,
            "partition_error": data["partition_error"],
            "partition_radial_closure": data["radial_closure"],
            "partition_axial_closure": data["axial_closure"],
            "plus_velocity_by_beta": plus_velocity,
            "minus_velocity_by_beta": minus_velocity,
            "velocity_physical_cylindrical_by_beta": physical_by_beta,
            "velocity_physical_cylindrical_total": cyl_total,
            "velocity_physical_cartesian_by_beta": cart_by_beta,
            "velocity_physical_cartesian_by_band": cart_by_band,
            "velocity_physical_cartesian_total": cart_total,
            "band_column_norms": np.asarray(band_column_norms, dtype=float),
            "distinct_dyadic_band_columns_materialized": True,
            "genuinely_independent_second_covariance_column_ready": False,
            "agent3_rank_screen_still_required": True,
        }

    def agent3_handoff(
        self,
        physical_family_result: dict[str, Any],
        *,
        theta: Any,
        delta_common: float = 0.0,
        delta_band: float = 0.0,
    ) -> dict[str, Any]:
        """Attach #400 bounded units without claiming covariance/rank success."""
        required = {
            "beta_labels",
            "active_ell_bands",
            "velocity_physical_cylindrical_by_beta",
            "velocity_physical_cylindrical_total",
            "velocity_physical_cartesian_by_band",
            "Q_source_by_label",
            "epsilon_source_by_label",
        }
        missing = sorted(required - set(physical_family_result))
        if missing:
            raise ValueError(f"multi-band physical result is missing required keys: {missing}")
        if len(tuple(physical_family_result["active_ell_bands"])) < 2:
            raise ValueError("Agent-3 handoff requires at least two distinct dyadic bands")
        coordinates = KokunoBoundedFamilyCoefficientCoordinates(
            max_l1_update=self.coefficient_max_l1_update
        ).consume_physical_family(
            physical_family_result,
            delta_common=delta_common,
            delta_band=delta_band,
            theta=theta,
        )
        return {
            "active_ell_bands": tuple(physical_family_result["active_ell_bands"]),
            "Q_source_by_label": np.asarray(physical_family_result["Q_source_by_label"], dtype=float),
            "epsilon_source_by_label": np.asarray(
                physical_family_result["epsilon_source_by_label"], dtype=float
            ),
            "rank_screen_columns_cartesian": np.asarray(
                physical_family_result["velocity_physical_cartesian_by_band"], dtype=float
            ),
            "bounded_coordinates": coordinates,
            "parameter_count": coordinates["parameter_count"],
            "parameter_units": coordinates["parameter_units"],
            "genuinely_independent_second_covariance_column_ready": False,
            "agent3_rank_screen_still_required": True,
            "handoff_scope": "distinct source-scheduled dyadic-band columns from caller-supplied mode data",
        }

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "source_formulas": dict(_SOURCE_FORMULAS),
            },
            "parameters": {
                "h": self.h,
                "A": self.A,
                "partition_atol": self.partition_atol,
                "derivative_rtol": self.derivative_rtol,
                "reality_atol": self.reality_atol,
                "coefficient_max_l1_update": self.coefficient_max_l1_update,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSourceMultiBandRealPairFamily":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected multi-band real-pair family schema")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("multi-band parameters are missing")
        obj = cls(
            h=params.get("h"),
            partition_atol=params.get("partition_atol"),
            derivative_rtol=params.get("derivative_rtol"),
            reality_atol=params.get("reality_atol"),
            coefficient_max_l1_update=params.get("coefficient_max_l1_update"),
        )
        expected = obj.to_payload()
        for key in ("schema", "source", "parameters", "truth_boundary"):
            if payload.get(key) != expected[key]:
                raise ValueError(f"multi-band {key} metadata changed")
        if "sha256" in payload and payload["sha256"] != obj.sha256:
            raise ValueError("multi-band sha256 mismatch")
        extra = set(payload) - {"schema", "source", "parameters", "truth_boundary", "sha256"}
        if extra:
            raise ValueError("multi-band payload contains unexpected keys")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSourceMultiBandRealPairFamily":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
