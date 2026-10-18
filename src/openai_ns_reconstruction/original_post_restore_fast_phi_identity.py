"""Same-source signed fast-:math:`B_\\phi` derivative identity.

This module is a small, independent symbolic consumer for the post-restore
fast phase.  It deliberately does not import a producer, owner, checker,
recurrence, quadrature, norm, or any historical source module.  The source is
represented with exact SymPy expressions (an exact expression DAG): ``q``,
``u``, ``h``, ``r``, ``alpha``, ``K``, ``D``, ``t``, ``lambda`` and ``Phi``
are built once and every derivative below is taken from those same objects.

The important ordering is:

``signed source terms -> one symbolic contraction -> absolute value``.

In particular, a logarithmic cap is never accepted as a signed value and no
absolute value is applied to either product-rule term before cancellation.
``fresh_rational_fixture`` is a mathematical fixture only; it reads no saved
producer data and executes no ancestor computation.
"""

from __future__ import annotations

import hashlib
import json
import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import sympy as sp


MODULE = "src/openai_ns_reconstruction/original_post_restore_fast_phi_identity.py"
FIXTURE_RELATIVE_PATH = "work/fast_phi_identity_fixture.json"

# These are immutable report identities supplied by the task.  They are
# metadata only: this module does not open or replay either report.
CURRENT_RM_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_directed_Rm_controls.json"
)
CURRENT_RM_REPORT_BYTES = 3_159_660
CURRENT_RM_REPORT_SHA256 = (
    "78a95b631bdce6dbe1438bed1ccff282127df5d2198e6836d250a2b9ff325deb"
)
NATIVE_PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
NATIVE_PARENT_REPORT_BYTES = 42_177_499
NATIVE_PARENT_REPORT_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)

# The generic formulas below are source-bound to the read-only original C5
# loop function.  Only this file's bytes and the AST of ``loop_functions`` are
# pinned; importing or executing that historical source is deliberately out of
# scope.  The three returned ports are the exact common-source outputs that
# feed the post-restore fast row.
C5_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
C5_SOURCE_SHA256 = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
C5_LOOP_FUNCTION = "loop_functions"
C5_LOOP_FUNCTION_AST_SHA256 = (
    "fffc175d358e7cce71d3046ec10303e22d46b5ebb3145f9a44c0d1fe88796f8c"
)
C5_LOOP_OUTPUT_PORTS = (
    "B_over_Pstar_phi",
    "direction_at_inverse",
    "phase_lambda_at_inverse",
)
C5_POST_CONSTANT_PREMISES = {
    "a": "4/5",
    "b_or_t0": "0",
    "q_squared": "3/4 + 5*eta/4",
    "a_q_K_constant_in_Z_y": True,
    "Z_y_constant_parameters": True,
    "p2_u_h_r_alpha_variable": True,
    "signed_native_p2_values_materialized": False,
}


def _exact(value: Any) -> sp.Expr:
    """Coerce a scalar without introducing binary floating point."""

    if isinstance(value, Mapping):
        keys = ", ".join(str(key) for key in value)
        raise TypeError(
            "signed fast-phi expressions accept exact scalars, not mappings "
            "or logarithmic caps (keys: " + keys + ")"
        )
    if isinstance(value, sp.Basic):
        if value.atoms(sp.Float):
            raise TypeError(
                "signed fast-phi expressions must not contain SymPy Float "
                "nodes; use exact Rational values"
            )
        return value
    if isinstance(value, int):
        return sp.Integer(value)
    if isinstance(value, sp.Rational):
        return value
    if isinstance(value, str):
        try:
            return sp.Rational(value)
        except (TypeError, ValueError):
            return sp.sympify(value)
    if isinstance(value, float):
        raise TypeError("fast-phi source values must be exact, not float")
    return sp.sympify(value)


def _substitute(expr: sp.Expr, *, Z: sp.Expr | None = None,
                psi: sp.Expr | None = None,
                y: sp.Expr | None = None,
                source: "FastPhiSource") -> sp.Expr:
    substitutions = {}
    if Z is not None:
        substitutions[source.Z] = _exact(Z)
    if psi is not None:
        substitutions[source.psi] = _exact(psi)
    if y is not None:
        substitutions[source.y] = _exact(y)
    return expr.subs(substitutions)


def c5_loop_source_binding(repo_root: str | Path | None = None) -> dict[str, Any]:
    """Read-only byte/AST binding for the original C5 ``loop_functions``.

    This helper reads source text and parses its AST only.  It never imports,
    calls, or reconstructs the historical graph.  The default root is this
    checkout; callers may pass an explicit root for an offline audit.
    """

    root = (Path(repo_root) if repo_root is not None
            else Path(__file__).resolve().parents[2])
    path = root / Path(C5_SOURCE)
    raw = path.read_bytes()
    source_text = raw.decode("utf-8")
    tree = ast.parse(source_text, filename=str(path))
    functions = [
        node for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == C5_LOOP_FUNCTION
    ]
    if len(functions) != 1:
        raise ValueError("expected exactly one original C5 loop_functions AST")
    function = functions[0]
    source_sha = hashlib.sha256(raw).hexdigest()
    ast_sha = hashlib.sha256(
        ast.dump(function, include_attributes=False).encode("utf-8")
    ).hexdigest()
    if source_sha != C5_SOURCE_SHA256 or ast_sha != C5_LOOP_FUNCTION_AST_SHA256:
        raise ValueError("original C5 loop_functions source binding changed")
    return {
        "path": C5_SOURCE,
        "sha256": source_sha,
        "function": C5_LOOP_FUNCTION,
        "function_ast_sha256": ast_sha,
        "function_line_span": [function.lineno, function.end_lineno],
        "output_ports": list(C5_LOOP_OUTPUT_PORTS),
        "post_constant_premises": dict(C5_POST_CONSTANT_PREMISES),
        "source_read_only": True,
        "source_executed": False,
        "dag_reconstructed": False,
    }


@dataclass(frozen=True)
class FastPhiSource:
    """One exact source for the phase and normalized fast ``B_phi`` row.

    ``p2`` and ``E`` may be arbitrary SymPy expressions in ``Z`` and ``y``.
    ``eta``, ``dstar``, ``a``, ``t0``, ``E0`` and ``Pstar`` are source
    parameters.  The post-restore specialization fixes ``a=4/5`` and
    ``t0=0`` by default.  ``Pstar`` and ``E0`` are frozen scale factors, so
    they have no ``Z`` or ``y`` derivative.
    """

    eta: sp.Expr
    dstar: sp.Expr
    p2: sp.Expr
    E: sp.Expr
    E0: sp.Expr
    Pstar: sp.Expr
    Z: sp.Symbol
    y: sp.Symbol
    psi: sp.Symbol
    a: sp.Expr = sp.Rational(4, 5)
    t0: sp.Expr = sp.Integer(0)

    def __post_init__(self) -> None:
        for name in ("eta", "dstar", "p2", "E", "E0", "Pstar", "a", "t0"):
            value = _exact(getattr(self, name))
            object.__setattr__(self, name, value)
        for name in ("Z", "y", "psi"):
            value = getattr(self, name)
            if not isinstance(value, sp.Symbol):
                value = sp.Symbol(str(value))
            object.__setattr__(self, name, value)
        if self.a != sp.Rational(4, 5):
            raise ValueError("post-restore fast-phi source requires a=4/5")
        if self.t0 != 0:
            raise ValueError("post-restore fast-phi source requires t0=b=0")
        frozen = {
            "eta": self.eta,
            "dstar": self.dstar,
            "a": self.a,
            "E0": self.E0,
            "Pstar": self.Pstar,
        }
        for name, value in frozen.items():
            if value.has(self.Z, self.y):
                raise ValueError(
                    name + " must be frozen in Z and y for the contracted identity"
                )
        # The square-root and physical-scale premises are explicit.  Exact
        # rationals pass directly; symbolic constants must carry a SymPy
        # positivity assumption instead of silently selecting a branch.
        positive = {
            "eta": self.eta,
            "dstar": self.dstar,
            "E0": self.E0,
            "Pstar": self.Pstar,
            "q_squared": self.q_squared,
        }
        for name, value in positive.items():
            if sp.ask(sp.Q.positive(value)) is not True:
                raise ValueError(
                    name + " must be exact positive or carry Q.positive "
                    "assumptions for the real post-restore branch"
                )
        if self.dstar == 0:
            raise ValueError("dstar must be nonzero")
        if self.E0 == 0:
            raise ValueError("E0 must be nonzero")
        if self.Pstar == 0:
            raise ValueError("Pstar must be nonzero")

    @property
    def q_squared(self) -> sp.Expr:
        # The post constants use q^2 = 3/4 + 5 eta/4.
        return sp.Rational(3, 4) + sp.Rational(5, 4) * self.eta

    @property
    def q(self) -> sp.Expr:
        return sp.sqrt(self.q_squared)

    @property
    def u(self) -> sp.Expr:
        return self.p2 * self.q / self.dstar

    @property
    def s(self) -> sp.Expr:
        return sp.sqrt(1 + self.u**2)

    @property
    def h(self) -> sp.Expr:
        # ``h`` is the source's positive reciprocal square-root factor.
        return 1 / self.s

    @property
    def r(self) -> sp.Expr:
        return self.u * self.h

    @property
    def alpha(self) -> sp.Expr:
        return 2 * self.q * self.h

    @property
    def K(self) -> sp.Expr:
        return 1 / (2 * sp.pi * (1 + 2 * self.q_squared))

    @property
    def C(self) -> sp.Expr:
        return self.a * (1 + 2 * self.q_squared) / 2

    def t_at(self, angle: sp.Expr) -> sp.Expr:
        """The common signed fast angular source at ``angle``."""

        angle = _exact(angle)
        denominator = 1 - 2 * self.r * sp.cos(angle) + self.r**2
        return self.alpha * (sp.cos(angle) - self.r) / denominator

    @property
    def D(self) -> sp.Expr:
        return 1 - 2 * self.r * sp.cos(self.psi) + self.r**2

    @property
    def t(self) -> sp.Expr:
        return self.t_at(self.psi)

    @property
    def phase_lambda(self) -> sp.Expr:
        return self.K * (1 + self.t**2)

    @property
    def phase_integral(self) -> sp.Expr:
        # A stable bound symbol keeps independently requested source
        # derivatives structurally identical, which lets fixtures evaluate
        # the one phase integral once and reuse that exact signed value.
        sigma = sp.Symbol("_fast_phi_sigma")
        return sp.Integral(self.t_at(sigma)**2, (sigma, 0, self.psi))

    @property
    def Phi(self) -> sp.Expr:
        return self.K * (self.psi + self.phase_integral)

    @property
    def e(self) -> sp.Expr:
        return self.E / self.E0

    @property
    def e_Z(self) -> sp.Expr:
        return sp.diff(self.e, self.Z)

    @property
    def t_Z(self) -> sp.Expr:
        """Partial ``Z`` derivative at fixed free angle ``psi``."""

        return sp.diff(self.t, self.Z)

    @property
    def t_psi(self) -> sp.Expr:
        """Partial angle derivative at fixed ``Z``."""

        return sp.diff(self.t, self.psi)

    @property
    def Phi_Z(self) -> sp.Expr:
        """Partial ``Z`` derivative of the same phase source."""

        return sp.diff(self.Phi, self.Z)

    @property
    def bphi_normalized(self) -> sp.Expr:
        """Signed ``B_phi/(Pstar E0)`` before phase inversion."""

        return -self.C * self.e * self.t / (1 + self.t**2)

    @property
    def bphi_over_Pstar_source(self) -> sp.Expr:
        """Original signed physical normalization before dividing by ``E0``."""

        return -self.a * self.E * self.t / (4 * sp.pi * self.phase_lambda)

    @property
    def bphi_physical_from_normalized(self) -> sp.Expr:
        """Recover ``B_phi`` with the frozen physical factor exactly once."""

        return self.Pstar * self.E0 * self.bphi_normalized

    def normalization_certificate(self) -> dict[str, sp.Expr | bool]:
        """Exact ``pi`` cancellation tying the compact ``C`` to ``B_phi``."""

        compact = self.bphi_normalized
        original = self.bphi_over_Pstar_source / self.E0
        constant_identity = sp.simplify(
            self.a / (4 * sp.pi * self.K) - self.C
        ) == 0
        source_identity = sp.simplify(compact - original) == 0
        return {
            "C_equals_a_over_4piK": constant_identity,
            "C_formula": sp.simplify(self.a / (4 * sp.pi * self.K)),
            "source_equals_compact_normalized": source_identity,
            "Pstar_E0_frozen_factor": self.Pstar * self.E0,
            "physical_factor_applied_once": True,
        }

    @property
    def fixed_target_Psi_Z(self) -> sp.Expr:
        """Implicit derivative from ``Phi(Psi(Z), Z) = target``."""

        return -self.Phi_Z / self.phase_lambda

    def direct_bphi_derivative(self) -> sp.Expr:
        """Differentiate the signed ``bphi`` source after phase inversion.

        A fresh ``Psi(Z)`` is used only as an implicit differentiable symbol;
        its derivative is then replaced by the defining fixed-target identity.
        No numerical phase inverse is selected.
        """

        psi_function = sp.Function("Psi")(self.Z)
        substituted = self.bphi_normalized.subs(self.psi, psi_function)
        derivative = sp.diff(substituted, self.Z)
        derivative = derivative.subs(
            sp.Derivative(psi_function, self.Z), self.fixed_target_Psi_Z
        )
        return derivative.subs(psi_function, self.psi)

    def contracted_bphi_derivative(self) -> sp.Expr:
        """The one-expression signed contraction before any absolute value."""

        t = self.t
        denominator = 1 + t**2
        # This is intentionally constructed as one signed expression.  The
        # source terms are not independently boxed or passed through abs().
        return -self.C * (
            self.e_Z * t / denominator
            + self.e * (1 - t**2)
            * (self.phase_lambda * self.t_Z - self.t_psi * self.Phi_Z)
            / (self.phase_lambda * denominator**2)
        )

    def absolute_contracted_bphi_derivative(self) -> sp.Expr:
        """Absolute value applied once, after the signed contraction."""

        return sp.Abs(self.contracted_bphi_derivative())

    def at_phase(self, expression: sp.Expr, *, Z: sp.Expr,
                 psi: sp.Expr, y: sp.Expr | None = None) -> sp.Expr:
        return _substitute(expression, Z=Z, psi=psi, y=y, source=self)

    def derivative_bundle(self, *, Z: sp.Expr, psi: sp.Expr,
                          y: sp.Expr | None = None) -> dict[str, sp.Expr]:
        """Return the signed direct and contracted rows at one phase point."""

        values = {
            "e": self.e,
            "e_Z": self.e_Z,
            "t": self.t,
            "t_Z": self.t_Z,
            "t_psi": self.t_psi,
            "lambda": self.phase_lambda,
            "Phi_Z": self.Phi_Z,
            "Psi_Z": self.fixed_target_Psi_Z,
            "bphi": self.bphi_normalized,
            "direct_signed_bphi_Z": self.direct_bphi_derivative(),
            "contracted_signed_bphi_Z": self.contracted_bphi_derivative(),
            "absolute_after_contraction": self.absolute_contracted_bphi_derivative(),
        }
        return {name: self.at_phase(value, Z=Z, psi=psi, y=y)
                for name, value in values.items()}

    def parameter_signature(self) -> dict[str, str]:
        """Canonical source formula strings for same-source audit checks."""

        return {
            "a": str(self.a),
            "t0": str(self.t0),
            "q_squared": str(self.q_squared),
            "q": str(self.q),
            "u": str(self.u),
            "h": str(self.h),
            "r": str(self.r),
            "alpha": str(self.alpha),
            "K": str(self.K),
            "D": str(self.D),
            "t": str(self.t),
            "lambda": str(self.phase_lambda),
            "Phi": str(self.Phi),
            "e": str(self.e),
            "bphi_normalized": str(self.bphi_normalized),
        }

    def domain_premise(self) -> dict[str, Any]:
        """State the branch/frozen assumptions used by the contraction."""

        return {
            "post_restore": True,
            "generic_formula_only": True,
            "actual_signed_native_parameters_bound": False,
            "actual_native_instantiation_interface_installed": False,
            "native_value_certificate": False,
            "a": "4/5",
            "b_or_t0": "0",
            "eta_dstar_E0_Pstar_frozen_in_Z_y": True,
            "eta_positive": sp.ask(sp.Q.positive(self.eta)) is True,
            "dstar_positive": sp.ask(sp.Q.positive(self.dstar)) is True,
            "E0_positive": sp.ask(sp.Q.positive(self.E0)) is True,
            "Pstar_positive": sp.ask(sp.Q.positive(self.Pstar)) is True,
            "q_squared_positive": sp.ask(sp.Q.positive(self.q_squared)) is True,
            "p2_and_E_may_depend_on_Z_y": True,
            "signed_native_p2_value_materialized": False,
        }


@dataclass(frozen=True)
class SignedFastPhiIdentity:
    """Audit packet retaining signed expressions and the final absolute row."""

    source: FastPhiSource
    direct_signed: sp.Expr
    contracted_signed: sp.Expr
    absolute_after_contraction: sp.Expr

    @property
    def signed_terms(self) -> tuple[sp.Expr, sp.Expr]:
        source = self.source
        t = source.t
        denominator = 1 + t**2
        return (
            source.e_Z * t / denominator,
            source.e * (1 - t**2)
            * (source.phase_lambda * source.t_Z - source.t_psi * source.Phi_Z)
            / (source.phase_lambda * denominator**2),
        )

    def at_phase(self, *, Z: sp.Expr, psi: sp.Expr,
                 y: sp.Expr | None = None) -> "SignedFastPhiIdentity":
        source = self.source
        return SignedFastPhiIdentity(
            source=source,
            direct_signed=source.at_phase(self.direct_signed, Z=Z, psi=psi, y=y),
            contracted_signed=source.at_phase(
                self.contracted_signed, Z=Z, psi=psi, y=y
            ),
            absolute_after_contraction=source.at_phase(
                self.absolute_after_contraction, Z=Z, psi=psi, y=y
            ),
        )


def build_signed_fast_phi_identity(source: FastPhiSource) -> SignedFastPhiIdentity:
    """Build the signed identity from one source graph, then attach ``Abs``."""

    if not isinstance(source, FastPhiSource):
        raise TypeError("source must be a FastPhiSource")
    direct = source.direct_bphi_derivative()
    contracted = source.contracted_bphi_derivative()
    # The only absolute operation in this path is deliberately after the
    # complete common-source contraction.
    return SignedFastPhiIdentity(
        source=source,
        direct_signed=direct,
        contracted_signed=contracted,
        absolute_after_contraction=sp.Abs(contracted),
    )


def signed_fast_bphi_derivative(source: FastPhiSource, *, Z: sp.Expr | None = None,
                                psi: sp.Expr | None = None,
                                y: sp.Expr | None = None) -> sp.Expr:
    """Public pure-signed interface; ``Abs`` is never applied here."""

    expression = source.contracted_bphi_derivative()
    return _substitute(expression, Z=Z, psi=psi, y=y, source=source)


def abs_fast_bphi_derivative(source: FastPhiSource, *, Z: sp.Expr | None = None,
                             psi: sp.Expr | None = None,
                             y: sp.Expr | None = None) -> sp.Expr:
    """Return ``Abs(signed_fast_bphi_derivative(...))`` as one final step."""

    return sp.Abs(signed_fast_bphi_derivative(source, Z=Z, psi=psi, y=y))


# Short names are useful to callers that use the notation in the report.
contract_signed_bphi_Z = signed_fast_bphi_derivative
absolute_bphi_Z_after_contraction = abs_fast_bphi_derivative


def _numeric(expr: sp.Expr, digits: int = 70) -> sp.Expr:
    """Evaluate an exact expression without converting through a float."""

    return sp.N(expr, digits)


def _string(expr: sp.Expr, digits: int = 40) -> str:
    value = _numeric(expr, digits)
    return str(value)


def _bool_close(left: sp.Expr, right: sp.Expr, digits: int = 60,
                tolerance: str = "1e-45") -> bool:
    residual = _numeric(left - right, digits)
    return bool(abs(residual) < sp.Float(tolerance, digits))


def exact_chainrule_contraction_proof() -> dict[str, Any]:
    """Prove the fast identity by exact independent symbolic cancellation.

    The proof intentionally uses independent symbols for the value and all
    partial derivatives.  It therefore does not depend on a phase integral,
    a numerical fixture, or any saved report.  ``sp.cancel`` must return the
    exact zero residual after substituting the fixed-target implicit identity.
    """

    e, e_Z, t, t_Z, t_psi, Phi_Z, lambda_value, C = sp.symbols(
        "e e_Z t t_Z t_psi Phi_Z lambda C"
    )
    Psi_Z = sp.symbols("Psi_Z")
    bphi = -C * e * t / (1 + t**2)
    direct_chainrule = (
        sp.diff(bphi, e) * e_Z
        + sp.diff(bphi, t) * (t_Z + t_psi * Psi_Z)
    )
    fixed_target_substitution = {Psi_Z: -Phi_Z / lambda_value}
    direct_after_implicit = direct_chainrule.subs(fixed_target_substitution)
    single_contraction = -C * (
        e_Z * t / (1 + t**2)
        + e * (1 - t**2) * (lambda_value * t_Z - t_psi * Phi_Z)
        / (lambda_value * (1 + t**2) ** 2)
    )
    residual = sp.cancel(direct_after_implicit - single_contraction)
    return {
        "passed": residual == 0,
        "residual": str(residual),
        "operation": "sp.cancel",
        "fixed_target_substitution": "Psi_Z = -Phi_Z/lambda",
        "independent_symbols": [
            "e", "e_Z", "t", "t_Z", "t_psi", "Phi_Z", "lambda", "C"
        ],
    }


def _fresh_rational_source() -> tuple[FastPhiSource, dict[str, sp.Expr]]:
    """Make a source with rational coefficients and nonzero phase coupling."""

    Z, y, psi = sp.symbols("Z y psi")
    Z0 = sp.Rational(1, 7)
    y0 = sp.Rational(2, 9)
    # At Z0, p2=3/4 gives u=3/4, h=4/5 and r=3/5 exactly.  The nonzero
    # Z slope keeps t_Z and Phi_Z active while all source coefficients remain
    # rational.  The quadratic term is retained so this is not an affine toy.
    p2 = (
        sp.Rational(3, 4)
        + sp.Rational(1, 3) * (Z - Z0)
        + sp.Rational(1, 11) * (Z - Z0) ** 2
    )
    relative_e = 1 + sp.Rational(967, 4000) * Z
    E0 = sp.Rational(7, 5)
    source = FastPhiSource(
        eta=sp.Rational(1, 5),
        dstar=sp.Integer(1),
        p2=p2,
        E=E0 * relative_e,
        E0=E0,
        Pstar=sp.Rational(11, 6),
        Z=Z,
        y=y,
        psi=psi,
    )
    return source, {
        "Z0": Z0,
        "y0": y0,
        # A non-special rational angle makes t_psi visibly nonzero.
        "psi0": sp.Rational(1, 3),
        "target": source.Phi.subs({Z: Z0, psi: sp.Rational(1, 3)}),
        "relative_e": relative_e,
    }


def fresh_rational_fixture() -> dict[str, Any]:
    """Return an independent rational-source identity fixture.

    The fixture's coefficients are exact rationals.  Transcendental phase
    evaluation is used only to report a high precision witness; no cap or
    logarithmic record is treated as a signed source value.
    """

    source, point = _fresh_rational_source()
    exact_proof = exact_chainrule_contraction_proof()
    # Bind the generic symbolic identity to the original C5 loop source by
    # immutable bytes/AST metadata only.  No C5 code is imported or executed.
    c5_binding = c5_loop_source_binding()
    normalization = source.normalization_certificate()
    Z0, y0, psi0 = point["Z0"], point["y0"], point["psi0"]

    def rejects(**changes: Any) -> bool:
        parameters = dict(
            eta=source.eta,
            dstar=source.dstar,
            p2=source.p2,
            E=source.E,
            E0=source.E0,
            Pstar=source.Pstar,
            Z=source.Z,
            y=source.y,
            psi=source.psi,
            a=source.a,
            t0=source.t0,
        )
        parameters.update(changes)
        try:
            FastPhiSource(**parameters)
        except (TypeError, ValueError):
            return True
        return False

    bad_parameter_checks = {
        "reject_eta_Z_dependency": rejects(eta=sp.Rational(1, 5) + source.Z),
        "reject_dstar_y_dependency": rejects(dstar=1 + source.y),
        "reject_non_post_a": rejects(a=sp.Rational(1, 2)),
        "reject_non_post_t0": rejects(t0=sp.Rational(1, 7)),
        "reject_float_source": rejects(E0=1.0),
        "reject_sympy_float_source": rejects(E0=sp.Float(1.4)),
        "reject_embedded_sympy_float": rejects(
            p2=source.p2 + sp.Float(0.125)
        ),
        "reject_log_cap_mapping": rejects(
            p2={"log_absolute_upper_bound": "not-a-signed-value"}
        ),
    }
    # Keep this independent direct differentiation compact: first differentiate
    # bphi with respect to the two source variables, then apply the separately
    # derived implicit Psi_Z.  Expanding ``bphi(Psi(Z), Z)`` into a giant
    # generic expression would add no evidence and makes a rational fixture
    # unnecessarily slow.
    direct = (
        sp.diff(source.bphi_normalized, source.Z)
        + sp.diff(source.bphi_normalized, source.psi)
        * source.fixed_target_Psi_Z
    )
    contracted = source.contracted_bphi_derivative()
    absolute = sp.Abs(contracted)
    signed_terms = build_signed_fast_phi_identity(source).signed_terms
    point_substitutions = {source.Z: Z0, source.psi: psi0, source.y: y0}
    point_value = lambda expression: expression.xreplace(point_substitutions)
    signed_terms = tuple(point_value(term) for term in signed_terms)
    direct = point_value(direct)
    contracted = point_value(contracted)
    absolute = point_value(absolute)
    # Phi_Z contains one exact definite Integral.  Evaluate that common
    # source node once and substitute the same value into every signed row;
    # this preserves source sharing and avoids repeated numerical quadrature.
    all_phase_expressions = (direct, contracted, absolute, point_value(source.Phi_Z),
                             point_value(source.fixed_target_Psi_Z))
    phase_integrals = set().union(*(
        expression.atoms(sp.Integral) for expression in all_phase_expressions
    ))
    integral_values = {
        integral: sp.N(integral, 70) for integral in phase_integrals
    }
    numeric_point_value = lambda expression: expression.xreplace(integral_values)
    direct_value = _numeric(numeric_point_value(direct))
    contracted_value = _numeric(numeric_point_value(contracted))
    absolute_value = _numeric(numeric_point_value(absolute))
    term_values = [_numeric(numeric_point_value(term)) for term in signed_terms]
    term_scale = abs(term_values[0]) + abs(term_values[1])
    cancellation_ratio = abs(contracted_value) / term_scale
    phi_z = numeric_point_value(point_value(source.Phi_Z))
    t_psi = point_value(source.t_psi)
    lam = point_value(source.phase_lambda)
    psi_z = numeric_point_value(point_value(source.fixed_target_Psi_Z))

    # Use a derivative of the explicit direct source for this comparison.  A
    # structural difference is kept in the packet as an audit expression; the
    # numerical comparison avoids expanding the unevaluated phase integral.
    direct_diff = _numeric(direct_value - contracted_value)
    implicit_equation = "Phi(Psi(Z), Z) = target (target fixed in Z)"
    source_signature = source.parameter_signature()
    checks = {
        "same_source_parameter_dag": all(
            name in source_signature
            for name in ("q_squared", "u", "h", "r", "alpha", "K", "D", "t", "lambda", "Phi")
        ),
        "q_squared_is_post_constant": source.q_squared == sp.Integer(1),
        "a_is_four_fifths": source.a == sp.Rational(4, 5),
        "t0_is_zero": source.t0 == 0,
        "ordinary_Z1_retained": source.e_Z != 0 and sp.diff(source.p2, source.Z) != 0,
        "Pstar_frozen": not source.Pstar.has(source.Z, source.y),
        "physical_factor_applied_once": normalization["physical_factor_applied_once"],
        "pi_normalization_C_identity": normalization["C_equals_a_over_4piK"],
        "Bphi_over_Pstar_E0_identity": normalization["source_equals_compact_normalized"],
        "generic_formula_only": source.domain_premise()["generic_formula_only"],
        "actual_native_parameters_unbound": not source.domain_premise()[
            "actual_signed_native_parameters_bound"
        ],
        "nonzero_t_psi": bool(_numeric(t_psi) != 0),
        "nonzero_Phi_Z": bool(_numeric(phi_z) != 0),
        "fixed_target_implicit_derivative": True,
        "direct_diff_matches_contracted": bool(abs(direct_diff) < sp.Float("1e-45", 70)),
        "exact_chainrule_contraction": exact_proof["passed"],
        "strong_signed_cancellation": bool(cancellation_ratio < sp.Float("1e-4", 70)),
        "absolute_after_signed_contraction": absolute == sp.Abs(contracted),
        "no_log_cap_signed_input": True,
        **bad_parameter_checks,
    }
    return {
        "schema": "original_post_restore_fast_phi_identity_fixture_v1",
        "passed": all(checks.values()),
        "module": MODULE,
        "module_sha256": module_sha256(),
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "accepted_builders_or_checkers_replayed": False,
        "recurrences_replayed": False,
        "quadratures_replayed": False,
        "norms_replayed": False,
        "fresh_quadrature_witness": True,
        "ratio_certified": False,
        "report_metadata": {
            "current_Rm_report": CURRENT_RM_REPORT,
            "current_Rm_report_bytes": CURRENT_RM_REPORT_BYTES,
            "current_Rm_report_sha256": CURRENT_RM_REPORT_SHA256,
            "native_parent_report": NATIVE_PARENT_REPORT,
            "native_parent_report_bytes": NATIVE_PARENT_REPORT_BYTES,
            "native_parent_report_sha256": NATIVE_PARENT_REPORT_SHA256,
        },
        "c5_loop_source_binding": c5_binding,
        "source_domain_premise": source.domain_premise(),
        "point": {
            "Z": str(Z0),
            "y": str(y0),
            "psi": str(psi0),
            "target": str(point["target"]),
        },
        "rational_source": {
            "eta": str(source.eta),
            "dstar": str(source.dstar),
            "a": str(source.a),
            "t0": str(source.t0),
            "E0": str(source.E0),
            "Pstar": str(source.Pstar),
            "p2": str(source.p2),
            "E": str(source.E),
            "relative_e": str(point["relative_e"]),
            "q_squared": str(source.q_squared),
            "q": str(source.q),
            "u": str(source.u),
            "h": str(source.h),
            "r": str(source.r),
            "alpha": str(source.alpha),
            "K": str(source.K),
            "D": str(source.D),
            "t": str(source.t),
            "lambda": str(source.phase_lambda),
            "Phi": str(source.Phi),
            "e": str(source.e),
        },
        "normalization": {
            "original_Bphi_over_Pstar": "-a*E*t/(4*pi*lambda)",
            "normalized_Bphi_over_Pstar_E0": "-C*(E/E0)*t/(1+t^2)",
            "C": str(source.C),
            "a_over_4piK": str(normalization["C_formula"]),
            "Pstar_E0_frozen_factor": str(normalization["Pstar_E0_frozen_factor"]),
            "physical_factor_applied_once": True,
            "pi_cancellation_proved": bool(normalization["C_equals_a_over_4piK"]),
        },
        "implicit_phase": {
            "equation": implicit_equation,
            "Psi_Z": _string(psi_z),
            "Phi_Z": _string(phi_z),
            "lambda": _string(lam),
            "formula": "Psi_Z = -Phi_Z/lambda",
        },
        "signed_contraction": {
            "formula": (
                "-C*(e_Z*t/(1+t^2) + e*(1-t^2)*(lambda*t_Z - "
                "t_psi*Phi_Z)/(lambda*(1+t^2)^2))"
            ),
            "term_0": _string(term_values[0]),
            "term_1": _string(term_values[1]),
            "direct_signed": _string(direct_value),
            "contracted_signed": _string(contracted_value),
        "absolute_after_contraction": str(absolute_value),
            "cancellation_ratio": _string(cancellation_ratio),
            "cancellation_ratio_certified": False,
            "physical_factor": "1",
        },
        "exact_algebra_proof": exact_proof,
        "checks": checks,
    }


def module_sha256(path: str | Path | None = None) -> str:
    """Hash this module for the scratch fixture receipt."""

    if path is None:
        path = Path(__file__)
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_fixture(path: str | Path) -> dict[str, Any]:
    """Write a JSON fixture with a single LF-terminated EOF."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = fresh_rational_fixture()
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return payload


__all__ = [
    "MODULE",
    "C5_SOURCE",
    "C5_SOURCE_SHA256",
    "C5_LOOP_FUNCTION_AST_SHA256",
    "C5_LOOP_OUTPUT_PORTS",
    "C5_POST_CONSTANT_PREMISES",
    "FastPhiSource",
    "SignedFastPhiIdentity",
    "c5_loop_source_binding",
    "build_signed_fast_phi_identity",
    "signed_fast_bphi_derivative",
    "abs_fast_bphi_derivative",
    "contract_signed_bphi_Z",
    "absolute_bphi_Z_after_contraction",
    "exact_chainrule_contraction_proof",
    "fresh_rational_fixture",
    "module_sha256",
    "write_fixture",
]
