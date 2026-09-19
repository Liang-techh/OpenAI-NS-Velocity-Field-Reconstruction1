"""Conditional coefficient-space envelopes for Kokuno's exact PA.10 remainders.

Pinned public provenance is ``KokunoYumeto/yang-mills-interacting-workbench`` at
commit ``143f6773feb424ad9ed3a8d116653200f20346b7``, file
``navier-stokes/navier_stokes_workbench.tex`` (corrected reader dated
2026-09-09, Zenodo 22678406).

The corrected reconstruction writes, with ``d=1-eta^2``,

    U   = U_* + Lambda^-1 u,
    B   = -2 D eta A(u) - d partial_eta A(u),
    W   = W_* + Lambda^-1 B,
    H_c = H_* + Lambda^-1 d u,

and the exact remainders

    L R1 = [W + h(1-2 eta U) + d u zeta_*] Phi
           + W Y Phi_Y + H_c Phi_eta,

    L R2 = [A(1-4 eta U_*) + d U_*'] u
           - 2 A eta Lambda^-1 u^2
           + W Y u_Y + H_* u_eta + d Lambda^-1 u u_eta
           - 4 A eta p + d p_eta - 2 eta Y p_Y.

This module turns *that displayed monomial bookkeeping* into an executable,
fail-closed bound calculator.  It deliberately does not invent the numerical
coefficient-space bounds still missing from the public reconstruction.  A
caller must supply nonnegative norm/Lipschitz envelopes for the primitive
factors on the radius-one ball.  In particular ``Y Phi_Y`` and ``Y u_Y`` are
not silently assumed bounded merely because ``Phi`` and ``u`` are bounded;
the source uses radial inverses / mixed derivative estimates for those terms,
and a future source certificate must bind them accordingly.

The product rule is the source coefficient-space algebra estimate

    ||F G||_rho <= C_sq^2 ||F||_rho ||G||_rho,

with ``C_sq=4*pi^2/3``.  For a product of several factors the pair estimate is
iterated.  Lipschitz envelopes use the source-described one-factor-at-a-time
replacement argument.  The resulting R1/R2 values may be fed into Agent-1's
existing conditional ``M/K`` arithmetic only as *diagnostic planning inputs*
until all primitive bounds are independently source-certified.

No velocity, forcing, pressure fit, residual threshold, or ST006 datum is
changed here.  This is not a PDE validation and not a paper-exact recovery.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

from .kokuno_pa10_operator_primitives import KokunoPA10OperatorPrimitiveBounds


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-remainder-ball-bounds-v1"

_SOURCE_FORMULAS = {
    "rescaled_U": "U=U_*+Lambda^-1 u",
    "B": "B=-2D eta A(u)-d partial_eta A(u)",
    "W": "W=W_*+Lambda^-1 B",
    "H_c": "H_c=H_*+Lambda^-1 d u",
    "R1": (
        "L R1=[W+h(1-2 eta U)+d u zeta_*] Phi"
        "+W Y Phi_Y+H_c Phi_eta"
    ),
    "R2": (
        "L R2=[A(1-4 eta U_*)+d U_*']u"
        "-2A eta Lambda^-1 u^2+W Y u_Y+H_*u_eta"
        "+d Lambda^-1 u u_eta-4A eta p+d p_eta-2eta Y p_Y"
    ),
    "product": "||FG||_rho <= C_sq^2 ||F||_rho ||G||_rho",
    "lipschitz": "replace factors one at a time in each remainder monomial",
    "pair_M": "M bounds ((1+T)^-1 J_2 R1/2, J_1 R2/2)",
    "pair_K": "K is the Lipschitz constant of the same pair on the radius-one ball",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_R1_R2_formulas_recorded": True,
    "exact_displayed_R1_R2_monomial_bookkeeping_executable": True,
    "source_product_algebra_bound_used": True,
    "one_factor_at_a_time_lipschitz_bookkeeping_executable": True,
    "conditional_remainder_ball_envelope_executable": True,
    "conditional_MK_bridge_executable": True,
    "diagnostic_primitive_bounds_are_source_bounds": False,
    "source_rho_machine_bound": False,
    "source_chi_multiplier_norm_machine_bound": False,
    "source_R1_radius_one_ball_norm_machine_bound": False,
    "source_R2_radius_one_ball_norm_machine_bound": False,
    "source_R1_radius_one_ball_lipschitz_machine_bound": False,
    "source_R2_radius_one_ball_lipschitz_machine_bound": False,
    "source_operator_constant_M_machine_bound": False,
    "source_operator_constant_K_machine_bound": False,
    "source_contraction_invariant_ball_machine_verified": False,
    "source_fixed_point_distance_machine_bound": False,
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
class BallFactorBound:
    """Norm and local Lipschitz envelopes for one coefficient-space factor."""

    norm: float
    lipschitz: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "norm", _nonnegative_finite(self.norm, "norm"))
        object.__setattr__(
            self,
            "lipschitz",
            _nonnegative_finite(self.lipschitz, "lipschitz"),
        )


@dataclass(frozen=True)
class KokunoPA10RemainderPrimitiveBounds:
    """Caller-owned primitive envelopes used by the exact R1/R2 algebra.

    Fixed axis multipliers should normally have Lipschitz value zero.  Moving
    quantities may have nonzero Lipschitz values with respect to the radius-one
    ball variables.  These values are *not* source-certified merely by being
    placed in this dataclass.

    The two logarithmic-radial derivative factors ``Y_Phi_Y`` and ``Y_u_Y``
    are explicit rather than inferred from ``Phi``/``u``.  This prevents the
    implementation from hiding the source's derivative-loss seam.
    """

    eta: BallFactorBound
    d: BallFactorBound
    L_inverse: BallFactorBound
    U_star: BallFactorBound
    U_star_eta: BallFactorBound
    W_star: BallFactorBound
    H_star: BallFactorBound
    zeta_star: BallFactorBound
    Phi: BallFactorBound
    Y_Phi_Y: BallFactorBound
    Phi_eta: BallFactorBound
    u: BallFactorBound
    Y_u_Y: BallFactorBound
    u_eta: BallFactorBound
    A_u: BallFactorBound
    A_u_eta: BallFactorBound
    p: BallFactorBound
    p_eta: BallFactorBound
    Y_p_Y: BallFactorBound

    def __post_init__(self) -> None:
        for item in fields(self):
            value = getattr(self, item.name)
            if not isinstance(value, BallFactorBound):
                raise TypeError(f"{item.name} must be a BallFactorBound")

    def to_dict(self) -> dict[str, dict[str, float]]:
        return {name: asdict(getattr(self, name)) for name in self.__dataclass_fields__}

    @classmethod
    def diagnostic_unit_fixture(cls) -> "KokunoPA10RemainderPrimitiveBounds":
        """Return a deterministic arithmetic-only fixture.

        Every fixed multiplier gets ``(1,0)`` and every moving factor gets
        ``(1,1)``.  This fixture has no source-estimate status; it exists only
        so CI can execute every algebraic path without smuggling in hidden
        numerical reconstruction data.
        """

        fixed = BallFactorBound(1.0, 0.0)
        moving = BallFactorBound(1.0, 1.0)
        return cls(
            eta=fixed,
            d=fixed,
            L_inverse=fixed,
            U_star=fixed,
            U_star_eta=fixed,
            W_star=fixed,
            H_star=fixed,
            zeta_star=fixed,
            Phi=moving,
            Y_Phi_Y=moving,
            Phi_eta=moving,
            u=moving,
            Y_u_Y=moving,
            u_eta=moving,
            A_u=moving,
            A_u_eta=moving,
            p=moving,
            p_eta=moving,
            Y_p_Y=moving,
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "KokunoPA10RemainderPrimitiveBounds":
        expected = set(cls.__dataclass_fields__)
        if set(payload) != expected:
            missing = sorted(expected - set(payload))
            extra = sorted(set(payload) - expected)
            raise ValueError(f"primitive-bound fields mismatch; missing={missing}, extra={extra}")
        return cls(
            **{
                name: BallFactorBound(
                    norm=float(payload[name]["norm"]),
                    lipschitz=float(payload[name]["lipschitz"]),
                )
                for name in expected
            }
        )


@dataclass(frozen=True)
class KokunoPA10RemainderBallBounds:
    """Exact displayed R1/R2 monomial envelope calculator."""

    h: float
    A: float
    D: float
    rescaling_lambda: float

    def __post_init__(self) -> None:
        h = _nonnegative_finite(self.h, "h")
        A = _positive_finite(self.A, "A")
        D = _positive_finite(self.D, "D")
        lam = _positive_finite(self.rescaling_lambda, "rescaling_lambda")
        if h >= 0.5:
            raise ValueError("h must be <0.5")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "A", A)
        object.__setattr__(self, "D", D)
        object.__setattr__(self, "rescaling_lambda", lam)

    @property
    def primitive_operator_bounds(self) -> KokunoPA10OperatorPrimitiveBounds:
        return KokunoPA10OperatorPrimitiveBounds()

    @property
    def product_constant(self) -> float:
        return self.primitive_operator_bounds.product_constant

    @staticmethod
    def scale(factor: BallFactorBound, scalar: float) -> BallFactorBound:
        s = abs(float(scalar))
        if not math.isfinite(s):
            raise ValueError("scalar must be finite")
        return BallFactorBound(s * factor.norm, s * factor.lipschitz)

    @staticmethod
    def add(*factors: BallFactorBound) -> BallFactorBound:
        return BallFactorBound(
            sum(item.norm for item in factors),
            sum(item.lipschitz for item in factors),
        )

    def product(self, *factors: BallFactorBound) -> BallFactorBound:
        """Iterate the source product algebra and one-factor Lipschitz rule."""
        if not factors:
            return BallFactorBound(1.0, 0.0)
        if not all(isinstance(item, BallFactorBound) for item in factors):
            raise TypeError("all product factors must be BallFactorBound values")
        n = len(factors)
        algebra = self.product_constant ** max(n - 1, 0)
        norm_product = math.prod(item.norm for item in factors)
        norm = algebra * norm_product
        lip_sum = 0.0
        for index, item in enumerate(factors):
            if item.lipschitz == 0.0:
                continue
            other_product = math.prod(
                other.norm for j, other in enumerate(factors) if j != index
            )
            lip_sum += item.lipschitz * other_product
        lipschitz = algebra * lip_sum
        if not math.isfinite(norm) or not math.isfinite(lipschitz):
            raise OverflowError("product envelope is outside float range")
        return BallFactorBound(norm, lipschitz)

    def _derived_factors(
        self, primitives: KokunoPA10RemainderPrimitiveBounds
    ) -> dict[str, BallFactorBound]:
        one = BallFactorBound(1.0, 0.0)
        inv_lambda = 1.0 / self.rescaling_lambda

        U = self.add(primitives.U_star, self.scale(primitives.u, inv_lambda))
        B = self.add(
            self.scale(self.product(primitives.eta, primitives.A_u), 2.0 * self.D),
            self.product(primitives.d, primitives.A_u_eta),
        )
        W = self.add(primitives.W_star, self.scale(B, inv_lambda))
        H_c = self.add(
            primitives.H_star,
            self.scale(self.product(primitives.d, primitives.u), inv_lambda),
        )
        one_minus_2_eta_U = self.add(
            one,
            self.scale(self.product(primitives.eta, U), 2.0),
        )
        R1_coefficient = self.add(
            W,
            self.scale(one_minus_2_eta_U, self.h),
            self.product(primitives.d, primitives.u, primitives.zeta_star),
        )
        R2_linear_coefficient = self.add(
            self.scale(
                self.add(
                    one,
                    self.scale(
                        self.product(primitives.eta, primitives.U_star),
                        4.0,
                    ),
                ),
                self.A,
            ),
            self.product(primitives.d, primitives.U_star_eta),
        )
        return {
            "U": U,
            "B": B,
            "W": W,
            "H_c": H_c,
            "one_minus_2_eta_U": one_minus_2_eta_U,
            "R1_coefficient": R1_coefficient,
            "R2_linear_coefficient": R2_linear_coefficient,
        }

    def remainder_envelopes(
        self, primitives: KokunoPA10RemainderPrimitiveBounds
    ) -> dict[str, Any]:
        """Return termwise and total conditional norm/Lipschitz envelopes."""
        if not isinstance(primitives, KokunoPA10RemainderPrimitiveBounds):
            raise TypeError("primitives must be KokunoPA10RemainderPrimitiveBounds")
        d = self._derived_factors(primitives)
        inv_lambda = 1.0 / self.rescaling_lambda

        r1_terms = {
            "coefficient_times_Phi": self.product(d["R1_coefficient"], primitives.Phi),
            "W_times_Y_Phi_Y": self.product(d["W"], primitives.Y_Phi_Y),
            "H_c_times_Phi_eta": self.product(d["H_c"], primitives.Phi_eta),
        }
        r1_numerator = self.add(*r1_terms.values())
        r1 = self.product(primitives.L_inverse, r1_numerator)

        r2_terms = {
            "linear_coefficient_times_u": self.product(
                d["R2_linear_coefficient"], primitives.u
            ),
            "2A_eta_Lambda_inverse_u_squared": self.scale(
                self.product(primitives.eta, primitives.u, primitives.u),
                2.0 * self.A * inv_lambda,
            ),
            "W_times_Y_u_Y": self.product(d["W"], primitives.Y_u_Y),
            "H_star_times_u_eta": self.product(primitives.H_star, primitives.u_eta),
            "d_Lambda_inverse_u_u_eta": self.scale(
                self.product(primitives.d, primitives.u, primitives.u_eta),
                inv_lambda,
            ),
            "4A_eta_p": self.scale(
                self.product(primitives.eta, primitives.p), 4.0 * self.A
            ),
            "d_times_p_eta": self.product(primitives.d, primitives.p_eta),
            "2_eta_Y_p_Y": self.scale(
                self.product(primitives.eta, primitives.Y_p_Y), 2.0
            ),
        }
        r2_numerator = self.add(*r2_terms.values())
        r2 = self.product(primitives.L_inverse, r2_numerator)

        result: dict[str, Any] = {
            "derived": {name: asdict(value) for name, value in d.items()},
            "R1_terms": {name: asdict(value) for name, value in r1_terms.items()},
            "R2_terms": {name: asdict(value) for name, value in r2_terms.items()},
            "R1_numerator": asdict(r1_numerator),
            "R2_numerator": asdict(r2_numerator),
            "R1": asdict(r1),
            "R2": asdict(r2),
        }
        flat = []
        for section in (result["derived"], result["R1_terms"], result["R2_terms"]):
            for value in section.values():
                flat.extend(value.values())
        flat.extend(result["R1"].values())
        flat.extend(result["R2"].values())
        if not all(math.isfinite(float(value)) for value in flat):
            raise OverflowError("remainder envelope contains a non-finite value")
        return result

    def diagnostic_conditional_MK(
        self,
        primitives: KokunoPA10RemainderPrimitiveBounds,
        *,
        multiplier_norm_chi: float,
        R1_vanishing_order: int = 0,
        R2_vanishing_order: int = 0,
    ) -> dict[str, float]:
        """Bridge conditional R1/R2 envelopes into #569's M/K arithmetic.

        This method is intentionally named ``diagnostic``.  It is not a route
        for upgrading caller-provided primitive bounds into source-certified
        M/K values.
        """
        remainder = self.remainder_envelopes(primitives)
        return self.primitive_operator_bounds.conditional_pair_envelopes(
            multiplier_norm_chi=multiplier_norm_chi,
            R1_ball_norm=remainder["R1"]["norm"],
            R2_ball_norm=remainder["R2"]["norm"],
            R1_lipschitz_norm=remainder["R1"]["lipschitz"],
            R2_lipschitz_norm=remainder["R2"]["lipschitz"],
            R1_vanishing_order=R1_vanishing_order,
            R2_vanishing_order=R2_vanishing_order,
        )

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def to_payload(
        self,
        primitives: KokunoPA10RemainderPrimitiveBounds,
        *,
        multiplier_norm_chi: float | None = None,
    ) -> dict[str, Any]:
        remainder = self.remainder_envelopes(primitives)
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
            "parameters": {
                "h": self.h,
                "A": self.A,
                "D": self.D,
                "rescaling_lambda": self.rescaling_lambda,
                "product_constant": self.product_constant,
            },
            "primitive_bounds": primitives.to_dict(),
            "remainder_envelopes": remainder,
            "truth_boundary": self.truth_boundary,
        }
        if multiplier_norm_chi is not None:
            mchi = _nonnegative_finite(multiplier_norm_chi, "multiplier_norm_chi")
            payload["diagnostic_conditional_MK"] = self.diagnostic_conditional_MK(
                primitives,
                multiplier_norm_chi=mchi,
            )
            payload["diagnostic_multiplier_norm_chi"] = mchi
        payload["sha256"] = hashlib.sha256(_canonical_json(payload).encode()).hexdigest()
        return payload

    @classmethod
    def from_payload(
        cls, payload: dict[str, Any]
    ) -> tuple["KokunoPA10RemainderBallBounds", KokunoPA10RemainderPrimitiveBounds]:
        data = copy.deepcopy(payload)
        expected_sha = data.pop("sha256", None)
        if not isinstance(expected_sha, str):
            raise ValueError("payload is missing sha256")
        actual_sha = hashlib.sha256(_canonical_json(data).encode()).hexdigest()
        if actual_sha != expected_sha:
            raise ValueError("payload sha256 mismatch")
        if data.get("schema") != SCHEMA:
            raise ValueError("unexpected schema")
        if data.get("source") != {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "corrected_release": CORRECTED_RELEASE,
            "corrected_release_date": CORRECTED_RELEASE_DATE,
        }:
            raise ValueError("source provenance mismatch")
        if data.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source formula map mismatch")
        if data.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth boundary mismatch")
        parameters = data["parameters"]
        obj = cls(
            h=float(parameters["h"]),
            A=float(parameters["A"]),
            D=float(parameters["D"]),
            rescaling_lambda=float(parameters["rescaling_lambda"]),
        )
        if not math.isclose(
            float(parameters["product_constant"]),
            obj.product_constant,
            rel_tol=0.0,
            abs_tol=0.0,
        ):
            raise ValueError("product constant mismatch")
        primitives = KokunoPA10RemainderPrimitiveBounds.from_dict(data["primitive_bounds"])
        replay = obj.to_payload(
            primitives,
            multiplier_norm_chi=data.get("diagnostic_multiplier_norm_chi"),
        )
        if replay != payload:
            raise ValueError("payload replay mismatch")
        return obj, primitives

    def save_json(
        self,
        path: str | Path,
        primitives: KokunoPA10RemainderPrimitiveBounds,
        *,
        multiplier_norm_chi: float | None = None,
    ) -> Path:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps(
                self.to_payload(primitives, multiplier_norm_chi=multiplier_norm_chi),
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        return destination

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> tuple["KokunoPA10RemainderBallBounds", KokunoPA10RemainderPrimitiveBounds]:
        payload = json.loads(Path(path).read_text())
        return cls.from_payload(payload)


__all__ = [
    "BallFactorBound",
    "KokunoPA10RemainderPrimitiveBounds",
    "KokunoPA10RemainderBallBounds",
]
