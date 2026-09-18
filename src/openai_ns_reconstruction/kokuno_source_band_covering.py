"""Deterministic source band/covering scales for Kokuno oscillatory fields.

The corrected 2026-09-09 reconstruction displays the per-dyadic-band relations

    Q = 2^(-ell),  epsilon = Q^h,  S_* = ell^2,
    i = floor(log_{T_g}(Q^(-1-h) / S_*)),
    c_i = T_g^i Q^(1+h),
    M_i = Lambda_g^i Q^(d_r/2),

with the exact bounds

    T_g^(-1) S_*^(-1) < c_i <= S_*^(-1),
    T_g^(-rho_g) epsilon^(-kappa_s) S_*^(-rho_g)
        < M_i <= epsilon^(-kappa_s) S_*^(-rho_g).

It also uses a product-grid mesh S_*^(-3) and proves that interacting bands
separated by at most four levels have a uniformly bounded covering-level
difference.  These displayed formulas are deterministic once ``ell`` and ``h``
are fixed, so this module makes that part executable.

The reader does *not* display a unique rational rectangle-center tuple, grid
origin, or rectangle radius r0.  This module therefore does not invent them and
does not claim to materialize the actual auxiliary-torus mode family or public
oscillatory velocity.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from typing import Any

from .kokuno_physical_evaluation import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_KAPPA_S,
    SOURCE_LAMBDA_G,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    SOURCE_T_G,
    KokunoPhysicalEvaluationMap,
)

SCHEMA = "kokuno-source-band-covering-v1"

_SOURCE_IDENTITIES = {
    "band_scales": "Q=2^(-ell); epsilon=Q^h; S_*=ell^2",
    "covering_level": "i=floor(log_{T_g}(Q^(-1-h)/S_*))",
    "covering_coefficients": "c_i=T_g^i Q^(1+h); M_i=Lambda_g^i Q^(d_r/2)",
    "c_i_bound": "T_g^(-1) S_*^(-1) < c_i <= S_*^(-1)",
    "M_i_bound": "T_g^(-rho_g) epsilon^(-kappa_s) S_*^(-rho_g) < M_i <= epsilon^(-kappa_s) S_*^(-rho_g)",
    "grid_mesh": "S_*^(-3)=ell^(-6)",
    "interacting_bands": "|ell-ell'|<=4 for the source enlarged-box interaction bound",
    "covering_level_difference": "|i(ell)-i(ell')| <= 1 + (4(1+h)log2 + 8/ell0)/log(T_g)",
}

_TRUTH_BOUNDARY = {
    "source_band_covering_schedule_executable": True,
    "source_per_band_Q_epsilon_Sstar_i_ci_Mi_bound": True,
    "source_interacting_band_level_bound_executable": True,
    "source_product_grid_mesh_bound": True,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_partition_grid_origin_recovered": False,
    "actual_positive_order_background_bound": False,
    "actual_auxiliary_torus_mode_family_bound": False,
    "public_xyz_t_velocity_correction_materialized": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _validate_ell(ell: int) -> int:
    if isinstance(ell, bool) or not isinstance(ell, int):
        raise ValueError("ell must be an integer dyadic-band index")
    if ell <= 4:
        raise ValueError("ell must exceed 4 on the source separated-band regime")
    return ell


def _validate_h(h: float) -> float:
    out = float(h)
    if not math.isfinite(out) or not (0.0 < out < 0.01):
        raise ValueError("h must satisfy the corrected-reader range 0<h<1/100")
    return out


@dataclass(frozen=True)
class KokunoSourceBandCovering:
    """Evaluate the corrected reader's deterministic per-band covering scales."""

    ell: int
    h: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "ell", _validate_ell(self.ell))
        object.__setattr__(self, "h", _validate_h(self.h))
        if self.covering_level < 0:
            raise ValueError("source covering level i must be nonnegative")
        diagnostics = self.bound_diagnostics()
        if not diagnostics["c_i_bound_passed"] or not diagnostics["M_i_bound_passed"]:
            raise ValueError("computed band scales violate the corrected-reader floor bounds")

    @property
    def log_Q(self) -> float:
        return -self.ell * math.log(2.0)

    @property
    def Q(self) -> float:
        return math.exp(self.log_Q)

    @property
    def epsilon(self) -> float:
        return math.exp(self.h * self.log_Q)

    @property
    def S_star(self) -> float:
        return float(self.ell * self.ell)

    @property
    def product_grid_mesh(self) -> float:
        return self.S_star ** -3.0

    @property
    def rho_g(self) -> float:
        return math.log(SOURCE_LAMBDA_G) / math.log(SOURCE_T_G)

    @property
    def d_r(self) -> float:
        return KokunoPhysicalEvaluationMap.source_radial_exponent(
            SOURCE_LAMBDA_G,
            SOURCE_T_G,
            self.h,
            kappa_s=SOURCE_KAPPA_S,
        )["d_r"]

    @property
    def covering_level_real(self) -> float:
        log_target = -(1.0 + self.h) * self.log_Q - math.log(self.S_star)
        return log_target / math.log(SOURCE_T_G)

    @property
    def covering_level(self) -> int:
        return math.floor(self.covering_level_real)

    @property
    def c_i(self) -> float:
        log_value = self.covering_level * math.log(SOURCE_T_G) + (1.0 + self.h) * self.log_Q
        return math.exp(log_value)

    @property
    def M_i(self) -> float:
        log_value = self.covering_level * math.log(SOURCE_LAMBDA_G) + 0.5 * self.d_r * self.log_Q
        return math.exp(log_value)

    def bound_diagnostics(self) -> dict[str, float | bool]:
        """Return scale-free checks of the two exact floor inequalities."""
        c_scaled = self.c_i * self.S_star
        m_scaled = self.M_i * (self.epsilon ** SOURCE_KAPPA_S) * (self.S_star ** self.rho_g)
        c_lower = SOURCE_T_G ** -1.0
        m_lower = SOURCE_T_G ** (-self.rho_g)
        slack = 64.0 * math.ulp(1.0)
        return {
            "c_i_scaled": c_scaled,
            "c_i_lower_scaled": c_lower,
            "M_i_scaled": m_scaled,
            "M_i_lower_scaled": m_lower,
            "c_i_bound_passed": (c_scaled > c_lower - slack) and (c_scaled <= 1.0 + slack),
            "M_i_bound_passed": (m_scaled > m_lower - slack) and (m_scaled <= 1.0 + slack),
        }

    @staticmethod
    def interacting_level_difference_bound(*, h: float, ell0: int) -> float:
        """Return the source four-band upper bound on |i(ell)-i(ell')|."""
        hh = _validate_h(h)
        e0 = _validate_ell(ell0)
        return 1.0 + (4.0 * (1.0 + hh) * math.log(2.0) + 8.0 / e0) / math.log(SOURCE_T_G)

    def interaction_receipt(self, other: "KokunoSourceBandCovering", *, ell0: int | None = None) -> dict[str, Any]:
        """Check the displayed covering-level bound for one interacting band pair."""
        if not isinstance(other, KokunoSourceBandCovering):
            raise TypeError("other must be KokunoSourceBandCovering")
        if self.h != other.h:
            raise ValueError("interacting source bands must use the same h")
        band_difference = abs(self.ell - other.ell)
        if band_difference > 4:
            raise ValueError("source interaction receipt is restricted to |ell-ell'|<=4")
        lower_band = min(self.ell, other.ell)
        e0 = lower_band if ell0 is None else _validate_ell(ell0)
        if e0 > lower_band:
            raise ValueError("ell0 cannot exceed the lower interacting band")
        bound = self.interacting_level_difference_bound(h=self.h, ell0=e0)
        actual = abs(self.covering_level - other.covering_level)
        return {
            "ell": self.ell,
            "ell_prime": other.ell,
            "band_difference": band_difference,
            "i": self.covering_level,
            "i_prime": other.covering_level,
            "level_difference": actual,
            "source_upper_bound": bound,
            "bound_passed": actual <= bound + 64.0 * math.ulp(bound),
        }

    def pulse_length_from_r0(self, r0: float) -> dict[str, Any]:
        """Apply L_s=2 r0/c_i while keeping r0 explicitly caller supplied.

        The corrected reader proves existence of a sufficiently small rectangle
        radius but does not publish a unique numerical value.  This helper never
        upgrades that choice to source-recovered data.
        """
        radius = float(r0)
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError("r0 must be finite and positive")
        return {
            "r0": radius,
            "r0_binding": "caller_supplied_not_source_recovered",
            "L_s": 2.0 * radius / self.c_i,
        }

    def receipt(self) -> dict[str, Any]:
        payload = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_identities": dict(_SOURCE_IDENTITIES),
            "band": {
                "ell": self.ell,
                "h": self.h,
                "Q": self.Q,
                "epsilon": self.epsilon,
                "S_star": self.S_star,
                "product_grid_mesh": self.product_grid_mesh,
                "covering_level_i": self.covering_level,
                "covering_level_real": self.covering_level_real,
                "c_i": self.c_i,
                "M_i": self.M_i,
                "rho_g": self.rho_g,
                "d_r": self.d_r,
            },
            "bounds": self.bound_diagnostics(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
        return payload
