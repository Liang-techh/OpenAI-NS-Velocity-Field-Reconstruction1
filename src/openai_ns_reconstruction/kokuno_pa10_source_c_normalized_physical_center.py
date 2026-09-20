"""Bind one source-complex normalization C to the executable PA.10 center.

Pinned public provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader
(Zenodo 22678406).

The corrected reconstruction requires

    phi_*(eta) = exp(Lambda * int_0^eta zeta_*(w) dw),
    C >= sup_{eta in Omega} |phi_*(eta)|,
    g = phi_*/C,
    F = g Phi.

Agent-1 #787 made ``F=(phi_*/C)Phi`` executable but deliberately retained the
older autonomous outer-schedule ``C=2`` and left the complex-domain source
normalization fail-closed.  This module closes exactly that numerical barrier
for one explicit repository existence-choice, ``C=1000``, and propagates the
same value through the already-existing outer schedule and physical-center
profile contract.  The choice is made only from the public source condition;
no NS residual, CR001 energy normalization, or held-out data is consulted.

The certificate uses the complex tube already machine-certified by the PA.10
source-axis-domain object.  On the real enlarged interval it partitions the
public rational

    zeta_* = -L H_*/(H_*^2+sigma_*^2)

and encloses each cell directly through interval bounds for ``H_*`` and ``L``.
This preserves the sign of zeta_* instead of replacing the real primitive by
an unusably large integral of |zeta_*|.  For the complex part, every point in
the tube is reached from its real projection by a vertical segment (or from an
endpoint by a cap segment); the same polynomial derivative bounds give an
absolute zeta_* envelope on that segment.  Therefore

    Re int_0^z zeta_* <= real_primitive_upper + tube_radius*zeta_abs_upper.

The resulting bound is an Agent-1 machine self-certificate for the configured
binary64 datum.  Independent Agent-4 admission is still required.  This is the
PA.10 contraction center, not the final corrected fixed point, and it is not
yet a Cartesian spacetime/global Kokuno velocity or PDE validation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from functools import cached_property
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
    KokunoPA10PhysicalCenterProfileContract,
)


SCHEMA = "kokuno-pa10-source-c-normalized-physical-center-v1"
SELECTED_SOURCE_C = 1000.0
DEFAULT_CERTIFICATE_INTERVALS = 2**21
DEFAULT_CHUNK_SIZE = 2**15
_PRIMITIVE_ROUNDOFF_PADDING = 1.0e-8

_SOURCE_FORMULAS = {
    "phi_star": "phi_*=exp(Lambda int_0^eta zeta_*(w)dw)",
    "zeta_star": "zeta_*=-L H_*/(H_*^2+sigma_*^2)",
    "source_C_normalization": "C>=sup_{eta in Omega}|phi_*(eta)|",
    "normalized_g": "g=phi_*/C; sup_Omega|g|<=1",
    "physical_F": "F=(phi_*/C)Phi=g Phi",
    "same_outer_C": "the same source-normalization C enters the outer schedule",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_complex_C_normalization_agent1_machine_certified": True,
    "source_complex_C_normalization_independently_admitted": False,
    "selected_C_is_repository_autonomous_source_choice": True,
    "selected_C_uses_same_outer_schedule_C": True,
    "selected_C_chosen_from_ns_residual": False,
    "selected_C_is_CR001_energy_normalizer": False,
    "source_normalized_physical_center_profiles_executable": True,
    "source_center_is_final_corrected_fixed_point": False,
    "fixed_point_correction_materialized": False,
    "global_leading_profile_reconstructed": False,
    "cartesian_spacetime_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "complete_kokuno_composite_velocity": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _round_guard(value: Any, multiplier: float = 128.0) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    return multiplier * np.finfo(float).eps * (np.abs(array) + 1.0)


def _next_up(value: float) -> float:
    return math.nextafter(float(value), math.inf)


def _normalized_default_profile() -> KokunoPA10PhysicalCenterProfileContract:
    """Rebuild #787 with exactly one changed autonomous scalar: outer/source C."""
    base = KokunoPA10PhysicalCenterProfileContract()
    payload = copy.deepcopy(base.configuration())
    axis = payload["center_profile_configuration"]["axis_profile_configuration"]
    axis["outer_schedule"]["C"] = SELECTED_SOURCE_C
    return KokunoPA10PhysicalCenterProfileContract.from_configuration(payload)


@dataclass(frozen=True)
class KokunoPA10SourceCNormalizedPhysicalCenter:
    """Executable PA.10 center with one source-complex-normalized C choice."""

    physical_profiles: KokunoPA10PhysicalCenterProfileContract = field(
        default_factory=_normalized_default_profile
    )
    certificate_intervals: int = DEFAULT_CERTIFICATE_INTERVALS
    chunk_size: int = DEFAULT_CHUNK_SIZE

    def __post_init__(self) -> None:
        if not isinstance(self.physical_profiles, KokunoPA10PhysicalCenterProfileContract):
            raise TypeError("physical_profiles must be KokunoPA10PhysicalCenterProfileContract")
        count = int(self.certificate_intervals)
        chunk = int(self.chunk_size)
        if isinstance(self.certificate_intervals, bool) or count != self.certificate_intervals:
            raise TypeError("certificate_intervals must be an integer")
        if isinstance(self.chunk_size, bool) or chunk != self.chunk_size:
            raise TypeError("chunk_size must be an integer")
        if count < 2**18 or count % 2 != 0:
            raise ValueError("certificate_intervals must be an even integer >=2^18")
        if not 1024 <= chunk <= 2**18:
            raise ValueError("chunk_size must lie in [1024,2^18]")
        if float(self.physical_profiles.C) != SELECTED_SOURCE_C:
            raise ValueError(
                f"this source-normalized contract requires outer/profile C={SELECTED_SOURCE_C:g}"
            )
        object.__setattr__(self, "certificate_intervals", count)
        object.__setattr__(self, "chunk_size", chunk)

    @property
    def axis_domain(self):
        return self.physical_profiles.axis_profiles.domain

    @property
    def Lambda(self) -> float:
        return float(self.physical_profiles.Lambda)

    @property
    def C(self) -> float:
        return float(self.physical_profiles.C)

    def _global_polynomial_bounds(self) -> dict[str, float]:
        domain = self.axis_domain
        R = 1.0 + float(domain.enlarged_real_margin)
        tube = float(domain.complex_tube_radius)
        h = float(domain.h)
        D = float(domain.D)
        j0 = float(domain.j0)
        z_abs = R + tube
        H_prime_real = _next_up(abs(D + 4.0) + 2.0 * abs(j0) * R + 12.0 * R * R)
        L_prime_real = _next_up(4.0 * abs(h) * R)
        H_prime_complex = _next_up(
            abs(D + 4.0) + 2.0 * abs(j0) * z_abs + 12.0 * z_abs * z_abs
        )
        L_prime_complex = _next_up(4.0 * abs(h) * z_abs)
        return {
            "R": R,
            "tube": tube,
            "h": h,
            "D": D,
            "j0": j0,
            "sigma": float(domain.sigma_star),
            "H_prime_real_upper": H_prime_real,
            "L_prime_real_upper": L_prime_real,
            "H_prime_complex_upper": H_prime_complex,
            "L_prime_complex_upper": L_prime_complex,
        }

    def _cell_bounds(
        self, indices: np.ndarray, *, constants: Mapping[str, float], step: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Enclose real zeta_* and complex vertical/cap |zeta_*| cellwise."""
        R = float(constants["R"])
        tube = float(constants["tube"])
        h = float(constants["h"])
        D = float(constants["D"])
        j0 = float(constants["j0"])
        sigma = abs(float(constants["sigma"]))
        half = 0.5 * step
        eps = np.finfo(float).eps

        x = -R + (indices.astype(float) + 0.5) * step
        x2 = x * x
        x3 = x2 * x
        H_mid = j0 + (D + 4.0) * x - j0 * x2 - 4.0 * x3
        L_mid = 1.0 - 2.0 * h * x2

        # Standard gamma_n-style evaluation padding for the configured binary64
        # coefficients, then analytic cell variation from the polynomial
        # derivative majorants.  The padding is intentionally much larger than
        # the actual Horner operation count.
        gamma_H = (32.0 * eps) / (1.0 - 32.0 * eps)
        gamma_L = (16.0 * eps) / (1.0 - 16.0 * eps)
        H_majorant = (
            abs(j0)
            + abs(D + 4.0) * np.abs(x)
            + abs(j0) * x2
            + 4.0 * np.abs(x3)
        )
        L_majorant = 1.0 + 2.0 * abs(h) * x2
        H_eval_error = gamma_H * H_majorant + _round_guard(H_mid, 64.0)
        L_eval_error = gamma_L * L_majorant + _round_guard(L_mid, 64.0)
        H_radius = float(constants["H_prime_real_upper"]) * half + H_eval_error
        L_radius = float(constants["L_prime_real_upper"]) * half + L_eval_error

        H_lower = np.nextafter(H_mid - H_radius, -np.inf)
        H_upper = np.nextafter(H_mid + H_radius, np.inf)
        L_lower = np.nextafter(L_mid - L_radius, -np.inf)
        L_upper = np.nextafter(L_mid + L_radius, np.inf)
        if np.any(L_lower <= 0.0):
            raise RuntimeError("real-cell L interval lost positivity")

        # f(H)=H/(H^2+sigma^2) has only the two critical points H=+-sigma.
        # Enclosing f on H_lower..H_upper therefore needs only endpoints plus
        # those critical values when they lie inside the interval.
        def ratio(value: np.ndarray) -> np.ndarray:
            return value / (value * value + sigma * sigma)

        f_lo_raw = ratio(H_lower)
        f_hi_raw = ratio(H_upper)
        f_lo = np.minimum(
            f_lo_raw - _round_guard(f_lo_raw),
            f_hi_raw - _round_guard(f_hi_raw),
        )
        f_hi = np.maximum(
            f_lo_raw + _round_guard(f_lo_raw),
            f_hi_raw + _round_guard(f_hi_raw),
        )
        critical = 1.0 / (2.0 * sigma)
        critical_guard = float(_round_guard(critical))
        contains_minus_sigma = (H_lower <= -sigma) & (H_upper >= -sigma)
        contains_plus_sigma = (H_lower <= sigma) & (H_upper >= sigma)
        f_lo = np.where(
            contains_minus_sigma,
            np.minimum(f_lo, -critical - critical_guard),
            f_lo,
        )
        f_hi = np.where(
            contains_plus_sigma,
            np.maximum(f_hi, critical + critical_guard),
            f_hi,
        )

        products = np.stack(
            (
                -L_lower * f_lo,
                -L_lower * f_hi,
                -L_upper * f_lo,
                -L_upper * f_hi,
            ),
            axis=0,
        )
        product_guard = _round_guard(products)
        zeta_lower = np.min(products - product_guard, axis=0)
        zeta_upper = np.max(products + product_guard, axis=0)

        # For z=x+iy with |y|<=tube, compare H(z),L(z) with their real
        # projections.  The denominator factorization
        # H^2+sigma^2=(H-i sigma)(H+i sigma) gives a positive lower bound for
        # both factors after subtracting the complex vertical variation.
        min_abs_H_real = np.where(
            (H_lower <= 0.0) & (H_upper >= 0.0),
            0.0,
            np.minimum(np.abs(H_lower), np.abs(H_upper)),
        )
        max_abs_H_real = np.maximum(np.abs(H_lower), np.abs(H_upper))
        H_vertical = _next_up(float(constants["H_prime_complex_upper"]) * tube)
        L_vertical = _next_up(float(constants["L_prime_complex_upper"]) * tube)
        factor_lower = (
            np.sqrt(min_abs_H_real * min_abs_H_real + sigma * sigma) - H_vertical
        )
        factor_lower = factor_lower - _round_guard(factor_lower)
        if np.any(factor_lower <= 0.0):
            raise RuntimeError("complex zeta denominator lower bound failed")
        numerator_upper = (L_upper + L_vertical) * (max_abs_H_real + H_vertical)
        zeta_abs_upper = numerator_upper / (factor_lower * factor_lower)
        zeta_abs_upper = np.nextafter(
            zeta_abs_upper + _round_guard(zeta_abs_upper, 256.0), np.inf
        )
        if np.any(~np.isfinite(zeta_abs_upper)):
            raise RuntimeError("complex zeta upper bound became non-finite")
        return zeta_lower, zeta_upper, zeta_abs_upper

    @cached_property
    def source_C_certificate(self) -> dict[str, Any]:
        """Return a deterministic upper bound for sup_Omega |phi_*|."""
        domain_certificate = self.axis_domain.complex_domain_certificate()
        if not domain_certificate["source_complex_neighborhood_machine_certified"]:
            raise RuntimeError("upstream complex source tube is not certified")

        constants = self._global_polynomial_bounds()
        R = float(constants["R"])
        tube = float(constants["tube"])
        count = int(self.certificate_intervals)
        step = (2.0 * R) / count
        mid = count // 2

        max_real_primitive = 0.0
        max_tube_primitive = 0.0

        def scan_positive() -> float:
            nonlocal max_real_primitive, max_tube_primitive
            total = np.longdouble(0.0)
            for start in range(mid, count, self.chunk_size):
                stop = min(start + self.chunk_size, count)
                indices = np.arange(start, stop, dtype=np.int64)
                _, zeta_upper, zeta_abs = self._cell_bounds(
                    indices, constants=constants, step=step
                )
                increments = np.asarray(zeta_upper * step, dtype=np.longdouble)
                if increments.size:
                    prefix = np.empty(increments.size, dtype=np.longdouble)
                    prefix[0] = total
                    if increments.size > 1:
                        prefix[1:] = total + np.cumsum(increments[:-1], dtype=np.longdouble)
                    cell_real = prefix + np.maximum(
                        np.asarray(zeta_upper, dtype=np.longdouble), 0.0
                    ) * np.longdouble(step)
                    cell_tube = cell_real + np.longdouble(tube) * np.asarray(
                        zeta_abs, dtype=np.longdouble
                    )
                    max_real_primitive = max(max_real_primitive, float(np.max(cell_real)))
                    max_tube_primitive = max(max_tube_primitive, float(np.max(cell_tube)))
                    total += np.sum(increments, dtype=np.longdouble)
            return float(total)

        def scan_negative() -> float:
            nonlocal max_real_primitive, max_tube_primitive
            total = np.longdouble(0.0)
            for high in range(mid, 0, -self.chunk_size):
                low = max(0, high - self.chunk_size)
                indices = np.arange(high - 1, low - 1, -1, dtype=np.int64)
                zeta_lower, _, zeta_abs = self._cell_bounds(
                    indices, constants=constants, step=step
                )
                oriented = -zeta_lower
                increments = np.asarray(oriented * step, dtype=np.longdouble)
                if increments.size:
                    prefix = np.empty(increments.size, dtype=np.longdouble)
                    prefix[0] = total
                    if increments.size > 1:
                        prefix[1:] = total + np.cumsum(increments[:-1], dtype=np.longdouble)
                    cell_real = prefix + np.maximum(
                        np.asarray(oriented, dtype=np.longdouble), 0.0
                    ) * np.longdouble(step)
                    cell_tube = cell_real + np.longdouble(tube) * np.asarray(
                        zeta_abs, dtype=np.longdouble
                    )
                    max_real_primitive = max(max_real_primitive, float(np.max(cell_real)))
                    max_tube_primitive = max(max_tube_primitive, float(np.max(cell_tube)))
                    total += np.sum(increments, dtype=np.longdouble)
            return float(total)

        right_endpoint_upper = scan_positive()
        left_endpoint_upper = scan_negative()

        # The Euclidean tube has two semicircular end caps.  Every cap point is
        # connected to the corresponding real endpoint by a straight segment
        # of length <=tube.  Bound zeta on that endpoint disk independently.
        h = float(constants["h"])
        D = float(constants["D"])
        j0 = float(constants["j0"])
        sigma = abs(float(constants["sigma"]))
        H_vertical = _next_up(float(constants["H_prime_complex_upper"]) * tube)
        L_vertical = _next_up(float(constants["L_prime_complex_upper"]) * tube)
        cap_upper = -math.inf
        for x, base_upper in ((-R, left_endpoint_upper), (R, right_endpoint_upper)):
            H_x = j0 + (D + 4.0) * x - j0 * x * x - 4.0 * x**3
            L_x = 1.0 - 2.0 * h * x * x
            factor = math.hypot(H_x, sigma) - H_vertical
            factor -= float(_round_guard(factor))
            if factor <= 0.0:
                raise RuntimeError("complex end-cap denominator lower bound failed")
            zeta_cap = (
                (abs(L_x) + L_vertical) * (abs(H_x) + H_vertical) / (factor * factor)
            )
            zeta_cap = _next_up(zeta_cap + float(_round_guard(zeta_cap, 256.0)))
            cap_upper = max(cap_upper, base_upper + tube * zeta_cap)

        primitive_upper = max(0.0, max_tube_primitive, cap_upper)
        primitive_upper = _next_up(primitive_upper + _PRIMITIVE_ROUNDOFF_PADDING)
        log_phi_upper = _next_up(self.Lambda * primitive_upper)
        if log_phi_upper >= math.log(np.finfo(float).max):
            raise OverflowError("source phi_* complex upper exceeds binary64")
        phi_upper = _next_up(math.exp(log_phi_upper))
        g_upper = _next_up(phi_upper / self.C)
        certified = bool(phi_upper < self.C and g_upper < 1.0)
        return {
            "source_complex_domain": "certified Euclidean tube around enlarged real I",
            "enlarged_real_half_width": R,
            "complex_tube_radius": tube,
            "certificate_intervals": count,
            "chunk_size": int(self.chunk_size),
            "cell_width": step,
            "real_primitive_max_upper": _next_up(max_real_primitive),
            "complex_tube_primitive_max_upper": primitive_upper,
            "log_phi_star_complex_sup_upper": log_phi_upper,
            "phi_star_complex_sup_upper": phi_upper,
            "configured_C": self.C,
            "C_over_phi_star_upper": self.C / phi_upper,
            "normalized_g_complex_sup_upper": g_upper,
            "upstream_complex_neighborhood_machine_certified": True,
            "source_complex_C_normalization_agent1_machine_certified": certified,
            "source_complex_C_normalization_independently_admitted": False,
            "certificate_interpretation": (
                "Agent-1 source-compatible machine self-certificate; independent A4 audit required"
            ),
        }

    def values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        return self.physical_profiles.values(X, eta)

    def derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        return self.physical_profiles.derivatives(X, eta)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "physical_profile_configuration": self.physical_profiles.configuration(),
            "certificate_intervals": int(self.certificate_intervals),
            "chunk_size": int(self.chunk_size),
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10SourceCNormalizedPhysicalCenter":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        if payload.get("schema") != SCHEMA:
            raise ValueError("configuration schema mismatch")
        return cls(
            physical_profiles=KokunoPA10PhysicalCenterProfileContract.from_configuration(
                payload["physical_profile_configuration"]
            ),
            certificate_intervals=int(payload["certificate_intervals"]),
            chunk_size=int(payload["chunk_size"]),
        )

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10SourceCNormalizedPhysicalCenter":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        certificate = copy.deepcopy(self.source_C_certificate)
        if not certificate["source_complex_C_normalization_agent1_machine_certified"]:
            raise RuntimeError("configured source C failed its own certificate")

        x0, x1 = self.physical_profiles.source_X_interval
        X = np.linspace(x0, x1, 17)
        eta = np.linspace(-1.0, 1.0, 19)
        XX, EE = np.meshgrid(X, eta, indexing="ij")
        values = self.values(XX, EE)
        schedule = self.axis_domain.binding.outer_schedule
        payload: dict[str, Any] = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "source_C_certificate": certificate,
            "normalization_contract": {
                "selected_C": self.C,
                "same_C_in_physical_F": bool(self.physical_profiles.C == self.C),
                "same_C_in_outer_schedule": bool(float(schedule.C) == self.C),
                "selection_basis": (
                    "public source complex-normalization condition only; no NS residual or CR001 energy data"
                ),
                "CR001_energy_normalization_used_as_source_C_proof": False,
                "source_C_used_as_CR001_energy_normalizer": False,
            },
            "configuration": self.configuration(),
            "machine_checks": {
                "source_C_condition_closed_at_agent1_self_certificate_level": bool(
                    certificate["source_complex_C_normalization_agent1_machine_certified"]
                ),
                "normalized_g_complex_sup_upper_lt_one": bool(
                    certificate["normalized_g_complex_sup_upper"] < 1.0
                ),
                "same_C_propagated_into_outer_schedule": bool(float(schedule.C) == self.C),
                "physical_center_F_nontrivial_on_probe": bool(
                    np.any(np.abs(values["F_0"]) > 0.0)
                ),
                "physical_center_U_nontrivial_on_probe": bool(
                    np.any(np.abs(values["U_0"]) > 0.0)
                ),
                "all_physical_center_values_finite": bool(
                    all(np.all(np.isfinite(value)) for value in values.values())
                ),
            },
            "scientific_gates": {
                "momentum_max_l2": 1.0e-3,
                "divergence_max_l2": 1.0e-5,
                "free_residual_defined_forcing_forbidden": True,
            },
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }
        payload["receipt_sha256"] = hashlib.sha256(
            _canonical_json(payload).encode("utf-8")
        ).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.report()["receipt_sha256"])

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    profile = KokunoPA10SourceCNormalizedPhysicalCenter()
    payload = profile.save_report(args.output)
    if args.config_output is not None:
        profile.save_configuration(args.config_output)
    print("receipt_sha256=", payload["receipt_sha256"])
    print("source_C_certificate=", payload["source_C_certificate"])
    print("machine_checks=", payload["machine_checks"])


if __name__ == "__main__":
    _main()
