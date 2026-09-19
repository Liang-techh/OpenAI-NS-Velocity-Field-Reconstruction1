"""Executable universal operator primitives behind Kokuno's PA.10 M/K bounds.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected reader dated
2026-09-09 (Zenodo 22678406).

The corrected reconstruction introduces the coefficient norm

    a_{alpha,beta} = 20^{-alpha} rho^{-beta} beta!
                     binom(alpha+beta,beta)
                     / ((alpha+1)^2 (beta+1)^2)

and proves a finite collection of universal operator estimates before defining
the radius-one-ball constants ``M`` and ``K``.  In particular, with

    C_sq = 8 sum_{n>=1} n^{-2},

one has ``||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho`` and, for a sequence
vanishing below radial degree ``b``,

    ||J_nu F||_rho <= 80 / ((b+1)(b+nu)) ||F||_rho,   nu in {1,2}.

For ``T = J_2 chi / 2`` and multiplier norm ``M_chi`` the source obtains

    ||T^k|| <= (40 M_chi)^k / (k! (k+1)!),

so the absolute Neumann series gives a computable upper bound for
``||(1+T)^-1||`` once an actual ``M_chi`` bound is supplied.

This module makes those *universal* primitives executable and can propagate
caller-supplied radius-one-ball bounds for R1/R2 into conditional M/K envelopes.
It deliberately does not invent the source-dependent values rho, M_chi, the
R1/R2 ball norms, or their Lipschitz constants.  Consequently it does not
promote M, K, the fixed-point radius, B_0/T_sh, PA.16, or PDE validation.
"""

from __future__ import annotations

from dataclasses import dataclass
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_pa10_fixed_point_radius import KokunoPA10FixedPointRadiusGate


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-operator-primitives-v1"

_SOURCE_FORMULAS = {
    "coefficient_weight": (
        "a_{alpha,beta}=20^{-alpha} rho^{-beta} beta! "
        "binom(alpha+beta,beta)/((alpha+1)^2(beta+1)^2)"
    ),
    "coefficient_norm": (
        "||F||_rho=sup_{alpha,beta,eta in I} "
        "|partial_eta^beta F_alpha(eta)|/a_{alpha,beta}"
    ),
    "C_sq": "C_sq=8 sum_{n>=1} n^-2",
    "product_bound": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
    "J_nu_coefficients": (
        "(J_nu F)_{alpha+1}=F_alpha/((alpha+1)(alpha+nu)); "
        "(J_nu F)_0=0"
    ),
    "J_nu_vanishing_bound": (
        "||J_nu F||_rho <= 80/((b+1)(b+nu)) ||F||_rho, nu in {1,2}"
    ),
    "mixed_derivative_template": (
        "J_nu[(partial_eta F)(Y partial_Y G)] <= "
        "(80/rho) C_sq^2 ||F||_rho ||G||_rho before extra factors"
    ),
    "T_definition": "T=J_2 chi/2, multiplication by chi before integration",
    "T_power_bound": "||T^k|| <= (40 M_chi)^k/(k!(k+1)!)",
    "inverse_bound": (
        "||(1+T)^-1|| <= sum_{k>=0}(40 M_chi)^k/(k!(k+1)!)"
    ),
    "M_pair": "M bounds ((1+T)^-1 J_2 R_1/2, J_1 R_2/2)",
    "K_pair": "K is the Lipschitz constant of the same pair on the radius-one ball",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_coefficient_space_formula_recorded": True,
    "source_universal_operator_primitives_executable": True,
    "source_product_constant_executable": True,
    "source_J_nu_bound_executable": True,
    "source_T_power_bound_executable": True,
    "source_inverse_absolute_series_bound_executable": True,
    "conditional_MK_envelope_executable": True,
    "source_rho_machine_bound": False,
    "source_chi_multiplier_norm_machine_bound": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
    "source_R1_radius_one_ball_lipschitz_machine_bound": False,
    "source_R2_radius_one_ball_lipschitz_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_contraction_factor_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
    "source_fixed_point_solved": False,
    "source_B0_dependencies_machine_bound": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _nonnegative_finite(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and >=0")
    return out


def _positive_finite(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be finite and >0")
    return out


@dataclass(frozen=True)
class KokunoPA10OperatorPrimitiveBounds:
    """Universal source-displayed coefficient-space bounds.

    No instance field is a reconstructed source constant.  Source-dependent
    quantities are explicit method arguments so a caller cannot accidentally
    serialize a diagnostic number as if it were part of the pinned source.
    """

    @property
    def C_sq(self) -> float:
        # 8 * zeta(2) = 8*pi^2/6 = 4*pi^2/3.
        return float(4.0 * math.pi * math.pi / 3.0)

    @property
    def product_constant(self) -> float:
        return float(self.C_sq * self.C_sq)

    def j_nu_bound(self, nu: int, vanishing_order: int = 0) -> float:
        """Return the coefficient-space norm factor for ``J_nu``.

        The source prints the ``nu=2`` instance explicitly.  The ``nu=1``
        instance follows from the displayed coefficient formula for ``J_nu``
        and the same exact weight-ratio estimate; using 80 keeps the shared
        source bound conservative.
        """
        if isinstance(nu, bool) or nu not in (1, 2):
            raise ValueError("nu must be 1 or 2")
        if isinstance(vanishing_order, bool) or not isinstance(vanishing_order, int):
            raise TypeError("vanishing_order must be an integer")
        if vanishing_order < 0:
            raise ValueError("vanishing_order must be >=0")
        b = int(vanishing_order)
        return float(80.0 / ((b + 1.0) * (b + float(nu))))

    def mixed_derivative_template_constant(self, rho: float) -> float:
        """Universal prefactor before extra undifferentiated factors.

        This evaluates the source's ``(80/rho)`` derivative-loss factor and
        the two convolution sums, each bounded by ``C_sq``.  A numerical rho
        supplied here is caller evidence, not a source-certified radius.
        """
        r = _positive_finite(rho, "rho")
        return float((80.0 / r) * self.product_constant)

    def t_power_bound(self, k: int, multiplier_norm_chi: float) -> float:
        """Evaluate ``(40 M_chi)^k/(k!(k+1)!)`` by stable recurrence."""
        if isinstance(k, bool) or not isinstance(k, int):
            raise TypeError("k must be an integer")
        if k < 0:
            raise ValueError("k must be >=0")
        m = _nonnegative_finite(multiplier_norm_chi, "multiplier_norm_chi")
        term = 1.0
        for j in range(1, k + 1):
            term *= (40.0 * m) / (float(j) * float(j + 1))
            if not math.isfinite(term):
                raise OverflowError("T-power bound is outside float range")
        return float(term)

    def inverse_one_plus_t_bound(
        self,
        multiplier_norm_chi: float,
        *,
        relative_tolerance: float = 2.0e-15,
        max_terms: int = 100000,
    ) -> float:
        """Absolute-series upper bound for ``||(1+T)^-1||``.

        The summand recurrence is
        ``a_k/a_{k-1}=40 M_chi/(k(k+1))``.  Summation stops only after the
        positive tail term is below a relative floating-point tolerance once
        the terms are decreasing.  This is numerical evaluation of the
        displayed convergent majorant, not an estimate of ``M_chi`` itself.
        """
        m = _nonnegative_finite(multiplier_norm_chi, "multiplier_norm_chi")
        tol = _positive_finite(relative_tolerance, "relative_tolerance")
        if isinstance(max_terms, bool) or not isinstance(max_terms, int):
            raise TypeError("max_terms must be an integer")
        if max_terms < 2:
            raise ValueError("max_terms must be >=2")
        total = 1.0
        term = 1.0
        previous = math.inf
        for k in range(1, max_terms + 1):
            term *= (40.0 * m) / (float(k) * float(k + 1))
            if not math.isfinite(term) or not math.isfinite(total):
                raise OverflowError("inverse absolute-series bound is outside float range")
            total += term
            decreasing = term <= previous
            if decreasing and term <= tol * total:
                return float(total)
            previous = term
        raise RuntimeError("inverse absolute-series bound did not converge within max_terms")

    def conditional_pair_envelopes(
        self,
        *,
        multiplier_norm_chi: float,
        R1_ball_norm: float,
        R2_ball_norm: float,
        R1_lipschitz_norm: float,
        R2_lipschitz_norm: float,
        R1_vanishing_order: int = 0,
        R2_vanishing_order: int = 0,
    ) -> dict[str, float]:
        """Propagate supplied remainder bounds into conditional M/K envelopes.

        The pair norm is conservatively taken as the maximum of its two
        component bounds.  All supplied remainder and multiplier bounds remain
        caller-owned evidence; this method never changes the source truth flags.
        """
        mchi = _nonnegative_finite(multiplier_norm_chi, "multiplier_norm_chi")
        r1 = _nonnegative_finite(R1_ball_norm, "R1_ball_norm")
        r2 = _nonnegative_finite(R2_ball_norm, "R2_ball_norm")
        l1 = _nonnegative_finite(R1_lipschitz_norm, "R1_lipschitz_norm")
        l2 = _nonnegative_finite(R2_lipschitz_norm, "R2_lipschitz_norm")
        inverse = self.inverse_one_plus_t_bound(mchi)
        j2 = self.j_nu_bound(2, R1_vanishing_order)
        j1 = self.j_nu_bound(1, R2_vanishing_order)

        M_phi = 0.5 * inverse * j2 * r1
        M_u = 0.5 * j1 * r2
        K_phi = 0.5 * inverse * j2 * l1
        K_u = 0.5 * j1 * l2
        values = {
            "inverse_one_plus_T_absolute_bound": float(inverse),
            "J2_R1_operator_factor": float(j2),
            "J1_R2_operator_factor": float(j1),
            "conditional_M_phi_component": float(M_phi),
            "conditional_M_u_component": float(M_u),
            "conditional_M_pair_max_envelope": float(max(M_phi, M_u)),
            "conditional_K_phi_component": float(K_phi),
            "conditional_K_u_component": float(K_u),
            "conditional_K_pair_max_envelope": float(max(K_phi, K_u)),
        }
        if not all(math.isfinite(value) for value in values.values()):
            raise OverflowError("conditional M/K envelope is outside float range")
        return values

    def diagnostic_fixed_point_gate_from_remainder_bounds(
        self,
        *,
        multiplier_norm_chi: float,
        R1_ball_norm: float,
        R2_ball_norm: float,
        R1_lipschitz_norm: float,
        R2_lipschitz_norm: float,
        R1_vanishing_order: int = 0,
        R2_vanishing_order: int = 0,
    ) -> KokunoPA10FixedPointRadiusGate:
        """Feed conditional envelopes into #560's formula-only theorem gate.

        The returned gate still reports ``source_operator_constant_*`` false;
        this helper is for arithmetic planning until all source-dependent input
        bounds above are independently machine-bound.
        """
        bounds = self.conditional_pair_envelopes(
            multiplier_norm_chi=multiplier_norm_chi,
            R1_ball_norm=R1_ball_norm,
            R2_ball_norm=R2_ball_norm,
            R1_lipschitz_norm=R1_lipschitz_norm,
            R2_lipschitz_norm=R2_lipschitz_norm,
            R1_vanishing_order=R1_vanishing_order,
            R2_vanishing_order=R2_vanishing_order,
        )
        return KokunoPA10FixedPointRadiusGate(
            operator_constant_M=bounds["conditional_M_pair_max_envelope"],
            operator_constant_K=bounds["conditional_K_pair_max_envelope"],
        )

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
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
            "universal_execution": {
                "C_sq": self.C_sq,
                "product_constant_C_sq_squared": self.product_constant,
                "J1_bound_b0": self.j_nu_bound(1, 0),
                "J2_bound_b0": self.j_nu_bound(2, 0),
                "interpretation": (
                    "source-displayed universal coefficient-space arithmetic only; "
                    "rho, M_chi, R1/R2 ball norms and Lipschitz norms are not supplied"
                ),
            },
            "truth_boundary": self.truth_boundary,
        }

    def to_payload(self) -> dict[str, Any]:
        payload = self.report()
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode()).hexdigest()
        return payload

    @property
    def sha256(self) -> str:
        return str(self.to_payload()["sha256"])

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n")

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10OperatorPrimitiveBounds":
        raw = copy.deepcopy(payload)
        supplied_sha = raw.pop("sha256", None)
        expected_sha = hashlib.sha256(_canonical_json(raw).encode()).hexdigest()
        if supplied_sha != expected_sha:
            raise ValueError("payload sha256 mismatch")
        if raw.get("schema") != SCHEMA:
            raise ValueError("schema mismatch")
        baseline = cls().report()
        if raw.get("source") != baseline["source"]:
            raise ValueError("source provenance mismatch")
        if raw.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source formula metadata mismatch")
        if raw.get("universal_execution") != baseline["universal_execution"]:
            raise ValueError("universal execution mismatch")
        if raw.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth boundary mismatch")
        return cls()

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoPA10OperatorPrimitiveBounds":
        return cls.from_payload(json.loads(Path(path).read_text()))
