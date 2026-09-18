"""Bind a certified autonomous rectangle witness to Kokuno slow labels.

The corrected 2026-09-09 reader fixes the *structure* needed by the oscillatory
lane: beta=(ell,a) slow labels, two auxiliary-rectangle signs sigma=+/- inside
each beta, finite rational torus centers separated under the relevant J_g
iterates, and the per-band pulse length L_s=2*r0/c_i.  It does not publish a
unique numerical center tuple, color assignment, grid origin, or r0.

This module therefore joins two already-audited Agent-2 components without
upgrading repository choices to source data:

* ``KokunoRationalRectangleSeparation`` supplies an exact, deterministic,
  repository-autonomous source-compatible center/r0 witness.
* ``KokunoSourceCompatiblePartitionRealization`` supplies autonomous beta
  labels and squared-partition weights obeying the displayed source contract.

For each beta this module deterministically assigns two distinct certified
centers to sigma=+/- and emits Q, epsilon and L_s aligned with that beta.  The
assignment is a repository convention, stable under permutation of the beta
axis.  It is useful input geometry for the source-phase-bound complete-curl
adapter, but it is NOT the unreleased Kokuno rectangle/mode realization and
cannot by itself materialize a paper-exact velocity field.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from typing import Any

import numpy as np

from .kokuno_rational_rectangle_separation import KokunoRationalRectangleSeparation
from .kokuno_source_band_covering import KokunoSourceBandCovering

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-autonomous-signed-rectangle-geometry-v1"

_SOURCE_CONTRACT = {
    "slow_label": "beta=(ell,a)",
    "signed_label": "gamma=(ell,a,sigma), sigma in {+,-}, assigned after the beta partition",
    "signed_rectangles": "the two sigma prototypes use disjoint auxiliary-torus rectangles",
    "center_separation": "c_mu != J_g^Delta c_nu mod Z^2 on the finite interaction range",
    "radius_guards": "4*r0*C_v<1 and 2*r0*C_v*(C_J+1)<d_c",
    "band_scales": "Q=2^(-ell), epsilon=Q^h",
    "pulse_length": "L_s=2*r0/c_i",
}

_AUTONOMOUS_CONVENTIONS = {
    "center_pool": "deterministic rational witness from KokunoRationalRectangleSeparation",
    "beta_color": "SHA256 canonical beta label modulo the certified center count",
    "signed_pair": "sigma_plus uses beta_color; sigma_minus uses (beta_color+1) mod color_count",
    "purpose": "one deterministic source-compatible geometry witness for executable diagnostics",
}

_TRUTH_BOUNDARY = {
    "source_signed_rectangle_structure_preserved": True,
    "source_rectangle_separation_guards_certified": True,
    "source_band_Q_epsilon_Ls_schedule_executable": True,
    "autonomous_source_compatible_signed_rectangle_witness_instantiated": True,
    "autonomous_source_compatible_partition_labels_consumable": True,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_actual_rectangle_labels_instantiated": False,
    "source_actual_partition_labels_instantiated": False,
    "actual_source_h_sigma_pulse_integrals_bound": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "actual_source_physical_covariance_rank_two_assessed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _canonical_beta(label: Any) -> tuple[int, tuple[int, int, int]]:
    try:
        ell_raw, a_raw = label
        a_tuple = tuple(a_raw)
    except (TypeError, ValueError) as exc:
        raise ValueError("each beta label must be (ell,(a1,a2,a3))") from exc
    if isinstance(ell_raw, bool) or not isinstance(ell_raw, (int, np.integer)):
        raise ValueError("beta ell must be an integer")
    ell = int(ell_raw)
    if ell <= 4:
        raise ValueError("beta ell must exceed 4 on the source separated-band regime")
    if len(a_tuple) != 3:
        raise ValueError("beta grid label a must contain exactly three integers")
    a: list[int] = []
    for value in a_tuple:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError("beta grid coordinates must be integers")
        a.append(int(value))
    return ell, (a[0], a[1], a[2])


def _center_payload(center: tuple[Fraction, Fraction]) -> list[dict[str, int]]:
    return [
        {"numerator": center[0].numerator, "denominator": center[0].denominator},
        {"numerator": center[1].numerator, "denominator": center[1].denominator},
    ]


@dataclass(frozen=True)
class KokunoAutonomousSignedRectangleGeometry:
    """Deterministic, provenance-labelled signed rectangle geometry witness."""

    h: float = 0.005
    ell0: int = 5
    color_count: int = 4
    denominator: int = 17
    safety: float = 0.5

    def __post_init__(self) -> None:
        h = float(self.h)
        if not math.isfinite(h) or not (0.0 < h < 0.01):
            raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
        if isinstance(self.ell0, bool) or not isinstance(self.ell0, (int, np.integer)):
            raise ValueError("ell0 must be an integer")
        ell0 = int(self.ell0)
        if ell0 <= 4:
            raise ValueError("ell0 must exceed 4")
        if isinstance(self.color_count, bool) or not isinstance(self.color_count, (int, np.integer)):
            raise ValueError("color_count must be an integer")
        color_count = int(self.color_count)
        if color_count < 2:
            raise ValueError("color_count must be at least two so sigma signs are disjoint")
        if isinstance(self.denominator, bool) or not isinstance(self.denominator, (int, np.integer)):
            raise ValueError("denominator must be an integer")
        denominator = int(self.denominator)
        if denominator < 2:
            raise ValueError("denominator must be at least two")
        safety = float(self.safety)
        if not math.isfinite(safety) or not (0.0 < safety < 1.0):
            raise ValueError("safety must lie strictly between zero and one")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "ell0", ell0)
        object.__setattr__(self, "color_count", color_count)
        object.__setattr__(self, "denominator", denominator)
        object.__setattr__(self, "safety", safety)
        # Construct once here so invalid autonomous choices fail closed at creation.
        _ = self.witness

    @property
    def witness(self) -> KokunoRationalRectangleSeparation:
        return KokunoRationalRectangleSeparation.autonomous_for_band_window(
            h=self.h,
            ell0=self.ell0,
            color_count=self.color_count,
            denominator=self.denominator,
            safety=self.safety,
        )

    def _center_pair_indices(self, beta: tuple[int, tuple[int, int, int]]) -> tuple[int, int]:
        encoded = _canonical_json([beta[0], list(beta[1])]).encode("utf-8")
        plus = int.from_bytes(hashlib.sha256(encoded).digest()[:8], "big") % self.color_count
        minus = (plus + 1) % self.color_count
        return plus, minus

    def instantiate(self, partition: dict[str, Any]) -> dict[str, Any]:
        """Attach certified signed auxiliary geometry to one partition result.

        ``partition`` is expected to be the output of the existing source-
        compatible partition contract/realization.  We recheck the key shape and
        squared-partition identities here so forged label/weight bookkeeping
        cannot silently acquire rectangle provenance.
        """
        if not isinstance(partition, dict):
            raise TypeError("partition must be a partition-result dictionary")
        required = ("eta", "D_r_eta", "D_z_eta", "beta_labels")
        if any(key not in partition for key in required):
            raise ValueError("partition result is missing eta/derivative/label fields")

        eta = np.asarray(partition["eta"], dtype=float)
        D_r_eta = np.asarray(partition["D_r_eta"], dtype=float)
        D_z_eta = np.asarray(partition["D_z_eta"], dtype=float)
        if eta.ndim < 1 or D_r_eta.shape != eta.shape or D_z_eta.shape != eta.shape:
            raise ValueError("partition eta and derivative arrays must share a nonempty beta axis")
        if not np.all(np.isfinite(eta)) or not np.all(np.isfinite(D_r_eta)) or not np.all(np.isfinite(D_z_eta)):
            raise ValueError("partition arrays must contain only finite real values")

        labels = tuple(_canonical_beta(label) for label in partition["beta_labels"])
        if len(labels) != eta.shape[-1]:
            raise ValueError("beta label count must equal the partition last-axis length")
        if len(set(labels)) != len(labels):
            raise ValueError("beta labels must be unique")
        if min(label[0] for label in labels) < self.ell0:
            raise ValueError("partition reaches a band below this geometry witness ell0")

        squared_sum = np.sum(eta * eta, axis=-1)
        if not np.allclose(squared_sum, 1.0, rtol=0.0, atol=2.0e-12):
            raise ValueError("partition fails sum_beta eta_beta^2=1")
        radial = np.sum(eta * D_r_eta, axis=-1)
        axial = np.sum(eta * D_z_eta, axis=-1)
        rscale = np.maximum(1.0, np.linalg.norm(eta, axis=-1) * np.linalg.norm(D_r_eta, axis=-1))
        zscale = np.maximum(1.0, np.linalg.norm(eta, axis=-1) * np.linalg.norm(D_z_eta, axis=-1))
        if np.any(np.abs(radial) > 2.0e-10 * rscale):
            raise ValueError("partition fails differentiated radial squared-partition closure")
        if np.any(np.abs(axial) > 2.0e-10 * zscale):
            raise ValueError("partition fails differentiated axial squared-partition closure")

        witness = self.witness
        q_values: list[float] = []
        epsilon_values: list[float] = []
        L_s_values: list[float] = []
        covering_levels: list[int] = []
        center_indices = np.empty((len(labels), 2), dtype=np.int64)
        centers_float = np.empty((len(labels), 2, 2), dtype=float)
        gamma_labels: list[tuple[int, tuple[int, int, int], str]] = []
        exact_centers: list[list[list[dict[str, int]]]] = []

        for j, beta in enumerate(labels):
            band = KokunoSourceBandCovering(ell=beta[0], h=self.h)
            pulse = band.pulse_length_from_r0(witness.r0)
            plus_index, minus_index = self._center_pair_indices(beta)
            if plus_index == minus_index:
                raise AssertionError("signed autonomous center assignment must be distinct")
            center_indices[j] = (plus_index, minus_index)
            pair_exact = (witness.centers[plus_index], witness.centers[minus_index])
            centers_float[j, 0] = tuple(float(x) for x in pair_exact[0])
            centers_float[j, 1] = tuple(float(x) for x in pair_exact[1])
            exact_centers.append([_center_payload(pair_exact[0]), _center_payload(pair_exact[1])])
            q_values.append(band.Q)
            epsilon_values.append(band.epsilon)
            L_s_values.append(float(pulse["L_s"]))
            covering_levels.append(band.covering_level)
            gamma_labels.extend(
                (
                    (beta[0], beta[1], "sigma_plus"),
                    (beta[0], beta[1], "sigma_minus"),
                )
            )

        if np.any(center_indices[:, 0] == center_indices[:, 1]):
            raise AssertionError("sigma signs must never share one auxiliary center")

        return {
            "beta_labels": labels,
            "gamma_labels": tuple(gamma_labels),
            "eta": eta,
            "D_r_eta": D_r_eta,
            "D_z_eta": D_z_eta,
            "Q_by_beta": np.asarray(q_values, dtype=float),
            "epsilon_by_beta": np.asarray(epsilon_values, dtype=float),
            "L_s_by_beta": np.asarray(L_s_values, dtype=float),
            "covering_level_by_beta": np.asarray(covering_levels, dtype=np.int64),
            "center_index_by_beta_sign": center_indices,
            "center_float_by_beta_sign": centers_float,
            "center_exact_by_beta_sign": exact_centers,
            "sign_labels": ("sigma_plus", "sigma_minus"),
            "r0": witness.r0,
            "r0_binding": witness.r0_binding,
            "center_binding": witness.center_binding,
            "rectangle_geometry_is_autonomous": True,
            "source_actual_rectangle_labels_instantiated": False,
            "source_actual_partition_labels_instantiated": False,
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    def receipt(self) -> dict[str, Any]:
        witness = self.witness
        payload = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_contract": dict(_SOURCE_CONTRACT),
            "autonomous_conventions": dict(_AUTONOMOUS_CONVENTIONS),
            "parameters": {
                "h": self.h,
                "ell0": self.ell0,
                "color_count": self.color_count,
                "denominator": self.denominator,
                "safety": self.safety,
            },
            "witness": {
                "delta_max": witness.delta_max,
                "centers": [_center_payload(center) for center in witness.centers],
                "d_c_squared": {
                    "numerator": witness.d_c_squared_exact.numerator,
                    "denominator": witness.d_c_squared_exact.denominator,
                },
                "r0": witness.r0,
                "injectivity_guard_passed": witness.injectivity_guard_passed,
                "separation_guard_passed": witness.separation_guard_passed,
                "center_binding": witness.center_binding,
                "r0_binding": witness.r0_binding,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
        return payload
