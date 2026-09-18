"""Source-compatible rational rectangle separation for Kokuno oscillatory labels.

The corrected 2026-09-09 reconstruction proves that the finite color labels can
be assigned rational torus centers c_mu with

    c_mu != J_g**Delta c_nu  (mod Z^2),  0 <= Delta <= Delta_max,

apart from the tautology Delta=0, mu=nu.  For such a finite family it defines a
strictly positive minimum torus separation d_c and chooses r0 small enough that

    4 r0 C_v < 1,
    2 r0 C_v (C_J + 1) < d_c,

where C_J=max ||J_g**Delta|| and C_v=|v_r|+|v_t|.  Those inequalities make the
enlarged rectangles injective and keep joined-label auxiliary preimages
disjoint.

The source proves existence; it does not publish a unique rational center tuple,
color count, denominator, or r0.  This module therefore provides an exact
certificate for caller-supplied rational centers and a deterministic,
repository-autonomous witness generator.  Autonomous witnesses are
source-compatible realizations, never recovered source data or paper-exact
parameters.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from typing import Any, Iterable

from .kokuno_physical_evaluation import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_B_G,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
    SOURCE_T_G,
    SOURCE_V_R,
    SOURCE_V_T,
)
from .kokuno_source_band_covering import KokunoSourceBandCovering

SCHEMA = "kokuno-rational-rectangle-separation-v1"

_SOURCE_IDENTITIES = {
    "center_constraints": "c_mu != J_g^Delta c_nu mod Z^2 for 0<=Delta<=Delta_max, excluding Delta=0 and mu=nu",
    "positive_separation": "d_c=min dist_T2(c_mu-J_g^Delta c_nu, Z^2)>0 over the finite non-tautological constraints",
    "C_J": "C_J=max_{0<=Delta<=Delta_max} ||J_g^Delta||",
    "C_v": "C_v=|v_r|+|v_t|",
    "injectivity_guard": "4*r0*C_v < 1",
    "separation_guard": "2*r0*C_v*(C_J+1) < d_c",
    "pulse_length": "L_s=2*r0/c_i",
}

_TRUTH_BOUNDARY = {
    "source_finite_rectangle_separation_structure_executable": True,
    "source_small_r0_inequalities_executable": True,
    "exact_rational_torus_collision_certificate": True,
    "source_rectangle_centers_recovered": False,
    "source_rectangle_radius_r0_recovered": False,
    "source_color_count_recovered": False,
    "source_center_denominator_recovered": False,
    "autonomous_rational_center_witness_available": True,
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


def _mod1(x: Fraction) -> Fraction:
    return x - math.floor(x)


def _coerce_fraction(value: Any) -> Fraction:
    if isinstance(value, bool):
        raise ValueError("torus coordinates must be rational numbers, not booleans")
    if isinstance(value, Fraction):
        return _mod1(value)
    if isinstance(value, int):
        return _mod1(Fraction(value, 1))
    raise ValueError("torus coordinates must be exact fractions or integers")


def _normalize_center(value: Iterable[Any]) -> tuple[Fraction, Fraction]:
    pair = tuple(value)
    if len(pair) != 2:
        raise ValueError("each torus center must contain exactly two rational coordinates")
    return (_coerce_fraction(pair[0]), _coerce_fraction(pair[1]))


def _apply_j_power(
    center: tuple[Fraction, Fraction], power: int
) -> tuple[Fraction, Fraction]:
    if isinstance(power, bool) or not isinstance(power, int) or power < 0:
        raise ValueError("covering power Delta must be a nonnegative integer")
    x, y = center
    for _ in range(power):
        x, y = _mod1(3 * x + y), _mod1(x + 5 * y)
    return x, y


def _wrapped_abs(delta: Fraction) -> Fraction:
    d = _mod1(delta)
    return min(d, 1 - d)


def _torus_distance_squared(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> Fraction:
    dx = _wrapped_abs(left[0] - right[0])
    dy = _wrapped_abs(left[1] - right[1])
    return dx * dx + dy * dy


def _constraint_minimum(
    centers: tuple[tuple[Fraction, Fraction], ...],
    delta_max: int,
) -> tuple[Fraction, tuple[int, int, int]]:
    minimum: Fraction | None = None
    argmin: tuple[int, int, int] | None = None
    for mu, c_mu in enumerate(centers):
        for nu, c_nu in enumerate(centers):
            for delta in range(delta_max + 1):
                if delta == 0 and mu == nu:
                    continue
                transformed = _apply_j_power(c_nu, delta)
                distance_sq = _torus_distance_squared(c_mu, transformed)
                if distance_sq == 0:
                    raise ValueError(
                        f"rational centers violate source separation at mu={mu}, nu={nu}, Delta={delta}"
                    )
                if minimum is None or distance_sq < minimum:
                    minimum = distance_sq
                    argmin = (mu, nu, delta)
    if minimum is None or argmin is None:
        raise ValueError("at least two rational centers are required")
    return minimum, argmin


def _source_c_v() -> float:
    return math.hypot(*SOURCE_V_R) + math.hypot(*SOURCE_V_T)


def _source_c_j(delta_max: int) -> float:
    # J_g is real symmetric positive definite with largest eigenvalue T_g.
    # Hence ||J_g^Delta||_2 = T_g^Delta and the finite maximum occurs at Delta_max.
    return SOURCE_T_G**delta_max


def _validate_delta_max(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("delta_max must be a nonnegative integer")
    return value


@dataclass(frozen=True)
class KokunoRationalRectangleSeparation:
    """Certificate the corrected reader's finite rectangle-separation guards."""

    centers: tuple[tuple[Fraction, Fraction], ...]
    delta_max: int
    r0: float
    center_binding: str = "caller_supplied_rational_not_source_recovered"
    r0_binding: str = "caller_supplied_not_source_recovered"
    witness_denominator: int | None = None
    witness_safety: float | None = None

    def __post_init__(self) -> None:
        centers = tuple(_normalize_center(center) for center in self.centers)
        if len(centers) < 2:
            raise ValueError("at least two rational centers are required")
        object.__setattr__(self, "centers", centers)
        object.__setattr__(self, "delta_max", _validate_delta_max(self.delta_max))

        radius = float(self.r0)
        if not math.isfinite(radius) or radius <= 0.0:
            raise ValueError("r0 must be finite and positive")
        object.__setattr__(self, "r0", radius)

        minimum_sq, _ = _constraint_minimum(centers, self.delta_max)
        if minimum_sq <= 0:
            raise ValueError("finite rational center family must have positive torus separation")
        if not self.injectivity_guard_passed:
            raise ValueError("r0 violates source injectivity guard 4*r0*C_v < 1")
        if not self.separation_guard_passed:
            raise ValueError("r0 violates source separation guard 2*r0*C_v*(C_J+1) < d_c")

        if self.witness_denominator is not None:
            if (
                isinstance(self.witness_denominator, bool)
                or not isinstance(self.witness_denominator, int)
                or self.witness_denominator < 2
            ):
                raise ValueError("witness_denominator must be an integer >=2")
        if self.witness_safety is not None:
            safety = float(self.witness_safety)
            if not math.isfinite(safety) or not (0.0 < safety < 1.0):
                raise ValueError("witness_safety must lie strictly between 0 and 1")

    @property
    def d_c_squared_exact(self) -> Fraction:
        return _constraint_minimum(self.centers, self.delta_max)[0]

    @property
    def closest_constraint(self) -> tuple[int, int, int]:
        return _constraint_minimum(self.centers, self.delta_max)[1]

    @property
    def d_c(self) -> float:
        return math.sqrt(float(self.d_c_squared_exact))

    @property
    def C_v(self) -> float:
        return _source_c_v()

    @property
    def C_J(self) -> float:
        return _source_c_j(self.delta_max)

    @property
    def injectivity_radius_upper(self) -> float:
        return 1.0 / (4.0 * self.C_v)

    @property
    def separation_radius_upper(self) -> float:
        return self.d_c / (2.0 * self.C_v * (self.C_J + 1.0))

    @property
    def strict_radius_upper(self) -> float:
        return min(self.injectivity_radius_upper, self.separation_radius_upper)

    @property
    def injectivity_guard_passed(self) -> bool:
        return 4.0 * self.r0 * self.C_v < 1.0

    @property
    def separation_guard_passed(self) -> bool:
        return 2.0 * self.r0 * self.C_v * (self.C_J + 1.0) < self.d_c

    @staticmethod
    def source_delta_max(*, h: float, ell0: int) -> int:
        """Choose an integer upper bound from #430's displayed four-band estimate."""
        bound = KokunoSourceBandCovering.interacting_level_difference_bound(h=h, ell0=ell0)
        return int(math.ceil(bound))

    @classmethod
    def autonomous_rational_witness(
        cls,
        *,
        color_count: int,
        denominator: int,
        delta_max: int,
        safety: float = 0.5,
    ) -> "KokunoRationalRectangleSeparation":
        """Build one deterministic source-compatible rational witness.

        The concrete grid, color count, greedy ordering and safety fraction are
        repository choices.  They are deliberately tagged autonomous and are not
        source-recovered hidden parameters.
        """
        if isinstance(color_count, bool) or not isinstance(color_count, int) or color_count < 2:
            raise ValueError("color_count must be an integer >=2")
        if isinstance(denominator, bool) or not isinstance(denominator, int) or denominator < 2:
            raise ValueError("denominator must be an integer >=2")
        dm = _validate_delta_max(delta_max)
        sf = float(safety)
        if not math.isfinite(sf) or not (0.0 < sf < 1.0):
            raise ValueError("safety must lie strictly between 0 and 1")

        selected: list[tuple[Fraction, Fraction]] = []
        for a in range(denominator):
            for b in range(denominator):
                candidate = (Fraction(a, denominator), Fraction(b, denominator))
                trial = tuple(selected + [candidate])
                if len(trial) == 1:
                    collision = False
                    for delta in range(1, dm + 1):
                        if _torus_distance_squared(candidate, _apply_j_power(candidate, delta)) == 0:
                            collision = True
                            break
                    if collision:
                        continue
                else:
                    try:
                        _constraint_minimum(trial, dm)
                    except ValueError:
                        continue
                selected.append(candidate)
                if len(selected) == color_count:
                    break
            if len(selected) == color_count:
                break
        if len(selected) != color_count:
            raise ValueError(
                "requested autonomous rational witness was not found on the chosen denominator grid"
            )

        centers = tuple(selected)
        minimum_sq, _ = _constraint_minimum(centers, dm)
        d_c = math.sqrt(float(minimum_sq))
        c_v = _source_c_v()
        c_j = _source_c_j(dm)
        upper = min(1.0 / (4.0 * c_v), d_c / (2.0 * c_v * (c_j + 1.0)))
        r0 = sf * upper
        return cls(
            centers=centers,
            delta_max=dm,
            r0=r0,
            center_binding="repository_autonomous_source_compatible_rational_witness",
            r0_binding="repository_autonomous_source_compatible_safety_fraction",
            witness_denominator=denominator,
            witness_safety=sf,
        )

    @classmethod
    def autonomous_for_band_window(
        cls,
        *,
        h: float,
        ell0: int,
        color_count: int,
        denominator: int,
        safety: float = 0.5,
    ) -> "KokunoRationalRectangleSeparation":
        return cls.autonomous_rational_witness(
            color_count=color_count,
            denominator=denominator,
            delta_max=cls.source_delta_max(h=h, ell0=ell0),
            safety=safety,
        )

    def constraint_rows(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for mu, c_mu in enumerate(self.centers):
            for nu, c_nu in enumerate(self.centers):
                for delta in range(self.delta_max + 1):
                    if delta == 0 and mu == nu:
                        continue
                    transformed = _apply_j_power(c_nu, delta)
                    distance_sq = _torus_distance_squared(c_mu, transformed)
                    rows.append(
                        {
                            "mu": mu,
                            "nu": nu,
                            "Delta": delta,
                            "distance_squared_numerator": distance_sq.numerator,
                            "distance_squared_denominator": distance_sq.denominator,
                            "separated": distance_sq > 0,
                        }
                    )
        return rows

    def pulse_lengths(
        self, bands: Iterable[KokunoSourceBandCovering]
    ) -> list[dict[str, Any]]:
        output: list[dict[str, Any]] = []
        for band in bands:
            if not isinstance(band, KokunoSourceBandCovering):
                raise TypeError("bands must contain KokunoSourceBandCovering instances")
            pulse = band.pulse_length_from_r0(self.r0)
            output.append(
                {
                    "ell": band.ell,
                    "covering_level_i": band.covering_level,
                    "c_i": band.c_i,
                    "r0": self.r0,
                    "L_s": pulse["L_s"],
                    "r0_binding": self.r0_binding,
                }
            )
        return output

    def receipt(self) -> dict[str, Any]:
        minimum = self.d_c_squared_exact
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
            "source_constants_reused": {
                "b_g": SOURCE_B_G,
                "v_r": list(SOURCE_V_R),
                "v_t": list(SOURCE_V_T),
                "T_g": SOURCE_T_G,
            },
            "realization": {
                "center_binding": self.center_binding,
                "r0_binding": self.r0_binding,
                "color_count": len(self.centers),
                "delta_max": self.delta_max,
                "centers": [
                    [
                        {"numerator": x.numerator, "denominator": x.denominator},
                        {"numerator": y.numerator, "denominator": y.denominator},
                    ]
                    for x, y in self.centers
                ],
                "witness_denominator": self.witness_denominator,
                "witness_safety": self.witness_safety,
            },
            "separation": {
                "d_c_squared_exact": {
                    "numerator": minimum.numerator,
                    "denominator": minimum.denominator,
                },
                "d_c": self.d_c,
                "closest_constraint": list(self.closest_constraint),
                "C_v": self.C_v,
                "C_J": self.C_J,
                "r0": self.r0,
                "injectivity_radius_upper": self.injectivity_radius_upper,
                "separation_radius_upper": self.separation_radius_upper,
                "strict_radius_upper": self.strict_radius_upper,
                "injectivity_guard_passed": self.injectivity_guard_passed,
                "separation_guard_passed": self.separation_guard_passed,
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()
        return payload
