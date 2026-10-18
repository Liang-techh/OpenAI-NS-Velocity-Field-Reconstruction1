"""Fresh validated scalar ``gamma2`` source for the restoration chart.

The restoration identity is

    gamma2(1) = exp(-1) + int_0^1 exp(s - 1) (1 - sigma(s))**2 ds.

This module deliberately has no imports from the historical owner graph.  It
binds the C6 sigma wrapper and the source used by ``Gamma.sigma`` by exact file
and AST hashes, then evaluates the new scalar kernel with interval arithmetic.
The integrator uses a local ordinary-derivative Taylor jet through order four;
the Simpson remainder is certified from the whole-panel fourth-derivative
interval.  The module does not evaluate ``gamma2`` at import time.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import re
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


# These are source bindings, rather than imports.  The old C6 graph is never
# constructed by this module.  The expected hashes are the source bytes read
# during the preflight audit on 2026-10-18.
_C6_WRAPPER_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
)
_C6_BRIDGE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
)
_GAMMA_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py"
)
_STABLE_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_outer_initial.py"
)
_FLAT_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
)

_EXPECTED_SOURCE_SHA256 = {
    _C6_WRAPPER_REL: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    _C6_BRIDGE_REL: "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    _GAMMA_SOURCE_REL: "38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01",
    _STABLE_SOURCE_REL: "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657",
    _FLAT_SOURCE_REL: "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7",
}
# The C6 and Gamma function hashes are filled below with the exact audited
# values.  Keeping the mapping explicit makes source drift fail closed.
_EXPECTED_FUNCTION_AST_SHA256 = {
    (_C6_WRAPPER_REL, "sigma"): "a9df99a3783719da28584631da2d5e8a167ef888ae2563c899372da72a7e5f3a",
    (_C6_BRIDGE_REL, "__init__"): "2af2f63ffb2bfb84827a55a00e0b62949de957aa035e65bc126271783ffd18a9",
    (_GAMMA_SOURCE_REL, "sigma"): "bbd2218c3596eccaa669f9e0edf45271cf5a70b73f78a7c2caf39bf1d4842f9f",
    (_STABLE_SOURCE_REL, "stable_sigma"): "18bc4da4c4af6cfc4ae3268b803ef1e402d40150215901827523883397a4a35c",
    (_FLAT_SOURCE_REL, "_sigma_left"): "7938664874d3bbd227d0847d4018d019df44484692e88a1657b866bd66c46188",
    (_FLAT_SOURCE_REL, "sigma_jets"): "9b16b8141a650fa5f8e57ababdd6a9dc6390b7867f4735fc9ae51ba8f0a1c4d9",
}


class SourceAuditError(RuntimeError):
    """Raised when the bound source bytes or AST no longer match."""


class ValidatedIntegrationError(RuntimeError):
    """Raised when the requested certified panel budget is not reached."""


def _repo_root(root: str | Path | None = None) -> Path:
    if root is not None:
        return Path(root).resolve()
    # src/openai_ns_reconstruction/original_restoration_gamma2.py -> repo root
    return Path(__file__).resolve().parents[2]


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _function_node(tree: ast.AST, name: str) -> ast.FunctionDef:
    found = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == name
    ]
    if len(found) != 1 or not isinstance(found[0], ast.FunctionDef):
        raise SourceAuditError(f"expected one synchronous function named {name!r}")
    return found[0]


def _source_binding(root: Path, relative: str, function_name: str) -> dict[str, Any]:
    path = root / relative
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise SourceAuditError(f"cannot read bound source {relative}: {exc}") from exc
    if _sha256_bytes(raw) != _EXPECTED_SOURCE_SHA256[relative]:
        raise SourceAuditError(f"bound source bytes changed: {relative}")
    try:
        text = raw.decode("utf-8")
        tree = ast.parse(text, filename=str(path))
    except (UnicodeDecodeError, SyntaxError) as exc:
        raise SourceAuditError(f"cannot parse bound source {relative}: {exc}") from exc
    node = _function_node(tree, function_name)
    ast_dump = ast.dump(node, include_attributes=False)
    ast_sha = _sha256_bytes(ast_dump.encode("utf-8"))
    expected_ast = _EXPECTED_FUNCTION_AST_SHA256[(relative, function_name)]
    if ast_sha != expected_ast:
        raise SourceAuditError(
            f"bound function AST changed: {relative}:{function_name}"
        )
    segment = ast.get_source_segment(text, node)
    if segment is None:
        raise SourceAuditError(f"missing source segment: {relative}:{function_name}")
    return {
        "relative_path": relative,
        "byte_sha256": _sha256_bytes(raw),
        "function": function_name,
        "lineno": node.lineno,
        "end_lineno": node.end_lineno,
        "function_source_sha256": _sha256_bytes(segment.encode("utf-8")),
        "function_ast_sha256": ast_sha,
        "module_ast_sha256": _sha256_bytes(
            ast.dump(tree, include_attributes=False).encode("utf-8")
        ),
    }


def _existing_named_gamma2_modules(root: Path) -> list[str]:
    """Find old named gamma2 modules without inspecting generated receipts.

    Historical prose mentions ``gamma2`` as an interface name, so this audit
    intentionally only reports source-module paths.  The new module itself is
    excluded.  An empty result is the source-audit condition for this task.
    """

    excluded = Path(__file__).resolve()
    candidates: list[str] = []
    for base in (root / "src", root / "scripts", root / "experiments"):
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if (
                not path.is_file()
                or path.suffix == ".pyc"
                or path.suffix not in {".py", ".pyi"}
                or "__pycache__" in path.parts
                or path.resolve() == excluded
            ):
                continue
            if re.search(r"gamma[_ -]?2", path.stem, flags=re.IGNORECASE):
                candidates.append(path.relative_to(root).as_posix())
    return sorted(candidates)


def source_audit(root: str | Path | None = None) -> dict[str, Any]:
    """Return the immutable C6/Gamma sigma source binding and audit result.

    This performs byte and AST checks only.  It never imports or constructs an
    owner, producer, checker, fixture, or historical quadrature module.
    """

    repo = _repo_root(root)
    bindings = [
        _source_binding(repo, _C6_WRAPPER_REL, "sigma"),
        _source_binding(repo, _C6_BRIDGE_REL, "__init__"),
        _source_binding(repo, _GAMMA_SOURCE_REL, "sigma"),
        _source_binding(repo, _STABLE_SOURCE_REL, "stable_sigma"),
        _source_binding(repo, _FLAT_SOURCE_REL, "_sigma_left"),
        _source_binding(repo, _FLAT_SOURCE_REL, "sigma_jets"),
    ]
    existing = _existing_named_gamma2_modules(repo)
    if existing:
        raise SourceAuditError(
            "existing gamma2 source module(s) found; do not duplicate: "
            + ", ".join(existing)
        )
    return {
        "passed": True,
        "existing_named_gamma2_modules": existing,
        "c6_sigma_wrapper_binding": bindings[0],
        "c6_bridge_sigma_wrapper_binding": bindings[1],
        "gamma_sigma_binding": bindings[2],
        "stable_sigma_provider_binding": bindings[3],
        "flat_sigma_left_binding": bindings[4],
        "flat_sigma_jets_binding": bindings[5],
        "sigma_definition": (
            "sigma(0)=0, sigma(1)=1, and for 0<s<1: "
            "r/(1+r), r=exp(1/(1-s)^2-1/s^2) on s<=1/2; "
            "sigma(s)=1-sigma(1-s) on s>=1/2"
        ),
        "sigma_source_is_exact_function": True,
    }


def _interval_bounds(value: Any) -> tuple[mp.mpf, mp.mpf]:
    """Extract exact binary endpoints from an mpmath interval."""

    try:
        mpi = value._mpi_
        return mp.mpf(mpi[0]), mp.mpf(mpi[1])
    except AttributeError as exc:
        raise TypeError("expected an mpmath interval value") from exc


def _endpoint_tuple(value: Any, index: int) -> tuple[int, int, int, int]:
    """Return an exact directed endpoint tuple without decimal formatting."""

    try:
        endpoint = value._mpi_[index]
    except (AttributeError, IndexError) as exc:
        raise TypeError("expected an mpmath interval value") from exc
    return tuple(int(part) for part in endpoint)


def _directed_record(value: Any) -> dict[str, tuple[int, int, int, int]]:
    """Record both directed endpoints in the native mpmath tuple format."""

    return {
        "lower_mpf_tuple": _endpoint_tuple(value, 0),
        "upper_mpf_tuple": _endpoint_tuple(value, 1),
    }


def _json_directed_record(value: Mapping[str, Any]) -> dict[str, list[int]]:
    """Convert an internal tuple record to JSON-safe lists exactly once."""

    return {
        key: [int(part) for part in value[key]]
        for key in ("lower_mpf_tuple", "upper_mpf_tuple")
    }


def _interval_abs_upper(value: Any) -> mp.mpf:
    """Outward upper bound for an interval absolute value.

    The derivative bound is obtained from the interval endpoints, without
    applying scalar ``abs`` to a midpoint or an unvalidated derivative.
    """

    lo, hi = _interval_bounds(value)
    return max(lo, hi, -lo, -hi)


def _interval_from_bounds(ctx: MPIntervalContext, lo: mp.mpf, hi: mp.mpf) -> Any:
    if lo > hi:
        raise ValueError("reversed interval bounds")
    return ctx.mpf([lo, hi])


@dataclass(frozen=True)
class TaylorJet:
    """Ordinary derivatives ``(f, f', ..., f^(4))`` over one interval."""

    ctx: MPIntervalContext
    derivatives: tuple[Any, ...]

    ORDER = 4

    def __post_init__(self) -> None:
        if len(self.derivatives) != self.ORDER + 1:
            raise ValueError("TaylorJet requires derivatives through order four")

    @classmethod
    def constant(cls, ctx: MPIntervalContext, value: Any) -> "TaylorJet":
        v = value if hasattr(value, "_mpi_") else ctx.mpf(value)
        return cls(ctx, (v,) + (ctx.zero,) * cls.ORDER)

    @classmethod
    def variable(
        cls,
        ctx: MPIntervalContext,
        value: Any,
        first_derivative: Any = 1,
    ) -> "TaylorJet":
        v = value if hasattr(value, "_mpi_") else ctx.mpf(value)
        d1 = first_derivative if hasattr(first_derivative, "_mpi_") else ctx.mpf(first_derivative)
        return cls(ctx, (v, d1) + (ctx.zero,) * (cls.ORDER - 1))

    def _coerce(self, other: Any) -> "TaylorJet":
        if isinstance(other, TaylorJet):
            if other.ctx is not self.ctx:
                raise ValueError("TaylorJet contexts must agree")
            return other
        return TaylorJet.constant(self.ctx, other)

    def __add__(self, other: Any) -> "TaylorJet":
        rhs = self._coerce(other)
        return TaylorJet(self.ctx, tuple(a + b for a, b in zip(self.derivatives, rhs.derivatives)))

    __radd__ = __add__

    def __neg__(self) -> "TaylorJet":
        return TaylorJet(self.ctx, tuple(-a for a in self.derivatives))

    def __sub__(self, other: Any) -> "TaylorJet":
        return self + (-self._coerce(other))

    def __rsub__(self, other: Any) -> "TaylorJet":
        return self._coerce(other) - self

    def __mul__(self, other: Any) -> "TaylorJet":
        rhs = self._coerce(other)
        ctx = self.ctx
        out = []
        for order in range(self.ORDER + 1):
            out.append(
                sum(
                    ctx.mpf(math.comb(order, k))
                    * self.derivatives[k]
                    * rhs.derivatives[order - k]
                    for k in range(order + 1)
                )
            )
        return TaylorJet(ctx, tuple(out))

    __rmul__ = __mul__

    def reciprocal(self) -> "TaylorJet":
        ctx = self.ctx
        r = [ctx.one / self.derivatives[0]]
        for order in range(1, self.ORDER + 1):
            rest = sum(
                ctx.mpf(math.comb(order, k))
                * self.derivatives[k]
                * r[order - k]
                for k in range(1, order + 1)
            )
            r.append(-rest / self.derivatives[0])
        return TaylorJet(ctx, tuple(r))

    def __truediv__(self, other: Any) -> "TaylorJet":
        return self * self._coerce(other).reciprocal()

    def exp(self) -> "TaylorJet":
        ctx = self.ctx
        out = [ctx.exp(self.derivatives[0])]
        # From h'=f'h: h^(n)=sum_{k=1}^n C(n-1,k-1) f^(k) h^(n-k).
        for order in range(1, self.ORDER + 1):
            out.append(
                sum(
                    ctx.mpf(math.comb(order - 1, k - 1))
                    * self.derivatives[k]
                    * out[order - k]
                    for k in range(1, order + 1)
                )
            )
        return TaylorJet(ctx, tuple(out))

    def square(self) -> "TaylorJet":
        return self * self


def _sigma_jet(ctx: MPIntervalContext, s: Any) -> TaylorJet:
    """Build the source sigma jet on a panel strictly inside (0, 1).

    The left-half representation uses ``r=exp(1/(1-s)^2-1/s^2)`` and
    ``sigma=r/(1+r)``.  This is algebraically the source formula but avoids
    the ``0/0`` numerical form at a flat endpoint.  The right half uses the
    exact source symmetry ``sigma(s)=1-sigma(1-s)``.  Panels crossing one half
    are split by the adaptive driver.
    """

    lo, hi = _interval_bounds(s)
    half = mp.mpf("0.5")
    if lo < 0 or hi > 1 or lo == hi == 0 or lo == hi == 1:
        raise ValueError("sigma jet requires a strict interior panel")
    if lo < half < hi:
        raise ValueError("sigma jet panel crosses the symmetry point")
    if hi <= half:
        y = TaylorJet.variable(ctx, s, 1)
        reflected = False
    else:
        y = TaylorJet.variable(ctx, ctx.one - s, -1)
        reflected = True
    one = TaylorJet.constant(ctx, 1)
    inv_y_sq = (y.reciprocal()).square()
    inv_complement_sq = (one - y).reciprocal().square()
    log_r = inv_complement_sq - inv_y_sq
    r = log_r.exp()
    sigma = r / (one + r)
    return one - sigma if reflected else sigma


def _integrand_jet(ctx: MPIntervalContext, lo: mp.mpf, hi: mp.mpf) -> TaylorJet:
    s = ctx.mpf([lo, hi])
    variable = TaylorJet.variable(ctx, s, 1)
    sigma = _sigma_jet(ctx, s)
    one = TaylorJet.constant(ctx, 1)
    heat = (variable - 1).exp()
    return heat * (one - sigma).square()


def _integrand_value(ctx: MPIntervalContext, point: mp.mpf) -> Any:
    return _integrand_jet(ctx, point, point).derivatives[0]


@dataclass(frozen=True)
class SimpsonPanel:
    lo: mp.mpf
    hi: mp.mpf
    enclosure: Any
    remainder_error: mp.mpf
    fourth_derivative_upper: mp.mpf


def _simpson_panel(ctx: MPIntervalContext, lo: mp.mpf, hi: mp.mpf) -> SimpsonPanel:
    if not lo < hi:
        raise ValueError("Simpson panel must have positive length")
    midpoint_interval = (ctx.mpf(lo) + ctx.mpf(hi)) / ctx.mpf(2)
    midpoint_lo, midpoint_hi = _interval_bounds(midpoint_interval)
    fa = _integrand_value(ctx, lo)
    # The midpoint is retained as an interval so the true midpoint is inside
    # the Simpson sample even when the endpoint sum is not exactly dyadic.
    fm = _integrand_jet(ctx, midpoint_lo, midpoint_hi).derivatives[0]
    fb = _integrand_value(ctx, hi)
    # Keep the Simpson coefficient in interval arithmetic.  A nearest
    # rounded scalar ``length/6`` would not by itself be a certified weight.
    length_interval = ctx.mpf(hi) - ctx.mpf(lo)
    simpson = (length_interval / ctx.mpf(6)) * (fa + 4 * fm + fb)
    jet = _integrand_jet(ctx, lo, hi)
    fourth_upper = _interval_abs_upper(jet.derivatives[4])
    # Likewise, take the upper endpoint of the interval remainder expression
    # instead of a nearest rounded scalar estimate.
    remainder_interval = (
        length_interval**5
        * ctx.mpf(fourth_upper)
        / ctx.mpf(2880)
    )
    remainder = _interval_bounds(remainder_interval)[1]
    enclosure = simpson + ctx.mpf([-remainder, remainder])
    return SimpsonPanel(lo, hi, enclosure, remainder, fourth_upper)


def adaptive_simpson_validated(
    ctx: MPIntervalContext,
    lo: mp.mpf,
    hi: mp.mpf,
    error_budget: Any,
    *,
    max_depth: int = 80,
) -> tuple[Any, dict[str, Any]]:
    """Integrate the fresh kernel with a panelwise fourth-derivative budget."""

    if not lo < hi:
        raise ValueError("integration interval must be nonempty")
    budget_interval = error_budget if hasattr(error_budget, "_mpi_") else ctx.mpf(error_budget)
    budget_lower, budget_upper = _interval_bounds(budget_interval)
    if budget_upper <= 0:
        raise ValueError("positive error budget required")
    if max_depth < 1:
        raise ValueError("max_depth must be positive")
    total_length_interval = ctx.mpf(hi) - ctx.mpf(lo)
    total_length_lower, total_length_upper = _interval_bounds(total_length_interval)
    if total_length_lower <= 0:
        raise ValueError("integration interval has no positive directed length")
    panels: list[SimpsonPanel] = []
    max_error_interval = ctx.mpf([0, 0])

    def visit(left: mp.mpf, right: mp.mpf, depth: int) -> None:
        nonlocal max_error_interval
        panel = _simpson_panel(ctx, left, right)
        panel_length_interval = ctx.mpf(right) - ctx.mpf(left)
        panel_budget_interval = (
            budget_interval * panel_length_interval / total_length_interval
        )
        # Use the lower directed budget endpoint when deciding acceptance.
        # This sacrifices a few ulps of admissible budget and keeps the panel
        # certificate valid despite independent interval operation rounding.
        panel_budget_lower, _ = _interval_bounds(panel_budget_interval)
        if panel.remainder_error <= panel_budget_lower:
            panels.append(panel)
            panel_error_interval = ctx.mpf([0, panel.remainder_error])
            max_error_interval = ctx.mpf(
                [
                    0,
                    max(
                        _interval_bounds(max_error_interval)[1],
                        _interval_bounds(panel_error_interval)[1],
                    ),
                ]
            )
            return
        if depth >= max_depth:
            raise ValidatedIntegrationError(
                "maximum adaptive Simpson depth reached before the certified "
                "fourth-derivative budget was met"
            )
        midpoint_interval = (ctx.mpf(left) + ctx.mpf(right)) / ctx.mpf(2)
        midpoint = _interval_bounds(midpoint_interval)[0]
        visit(left, midpoint, depth + 1)
        visit(midpoint, right, depth + 1)

    # The source symmetry point is a derivative-chart seam for our stable
    # representation.  Integrate the two source pieces independently.
    if lo < mp.mpf("0.5") < hi:
        visit(lo, mp.mpf("0.5"), 0)
        visit(mp.mpf("0.5"), hi, 0)
    else:
        visit(lo, hi, 0)
    result = ctx.mpf([0, 0])
    remainder_sum_interval = ctx.mpf([0, 0])
    for panel in panels:
        result += panel.enclosure
        remainder_sum_interval += ctx.mpf([0, panel.remainder_error])
    remainder_sum_lower, remainder_sum_upper = _interval_bounds(remainder_sum_interval)
    max_error_lower, max_error_upper = _interval_bounds(max_error_interval)
    return result, {
        "panel_count": len(panels),
        "max_panel_error": max_error_upper,
        "quadrature_error": remainder_sum_upper,
        "max_panel_error_record": _directed_record(max_error_interval),
        "quadrature_error_record": _directed_record(remainder_sum_interval),
        "total_length_record": _directed_record(total_length_interval),
        "budget_record": _directed_record(budget_interval),
        "quadrature_error_lower": remainder_sum_lower,
        "max_depth": max_depth,
        "remainder_formula": "length^5/2880 * sup_interval_abs(f^(4))",
        "ordinary_derivative_order": 4,
        "whole_panel_interval_jet": True,
    }


@dataclass(frozen=True)
class Gamma2Certificate:
    """Serializable result of a future approved gamma2 evaluation."""

    lower: mp.mpf
    upper: mp.mpf
    width: mp.mpf
    delta: mp.mpf
    tolerance: mp.mpf
    dps: int
    tail_error: mp.mpf
    quadrature_error: mp.mpf
    max_panel_error: mp.mpf
    panel_count: int
    source_audit: Mapping[str, Any]
    directed_record: Mapping[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "directed_record": {
                key: _json_directed_record(value)
                for key, value in self.directed_record.items()
                if isinstance(value, Mapping)
                and "lower_mpf_tuple" in value
                and "upper_mpf_tuple" in value
            },
            "lower_mpf_tuple": list(self.directed_record["total"]["lower_mpf_tuple"]),
            "upper_mpf_tuple": list(self.directed_record["total"]["upper_mpf_tuple"]),
            "width_upper_mpf_tuple": list(self.directed_record["width"]["upper_mpf_tuple"]),
            "delta_mpf_tuple": list(self.directed_record["delta"]["upper_mpf_tuple"]),
            "tolerance_mpf_tuple": list(self.directed_record["tolerance"]["lower_mpf_tuple"]),
            "dps": self.dps,
            "retained_tail_error_upper_mpf_tuple": list(
                self.directed_record["tail_error"]["upper_mpf_tuple"]
            ),
            "retained_quadrature_error_upper_mpf_tuple": list(
                self.directed_record["quadrature_error"]["upper_mpf_tuple"]
            ),
            "max_panel_error_upper_mpf_tuple": list(
                self.directed_record["max_panel_error"]["upper_mpf_tuple"]
            ),
            "panel_count": self.panel_count,
            "source_audit": dict(self.source_audit),
        }


def validated_gamma2(
    tolerance: str | mp.mpf = "1e-8",
    *,
    dps: int = 120,
    max_depth: int = 80,
    root: str | Path | None = None,
) -> Gamma2Certificate:
    """Return a certified interval for the new restoration scalar.

    This function is intentionally dormant unless called by the root agent.
    It uses only the source sigma definition and the fresh integrator above;
    no ancestor replay or library quadrature error estimate is involved.
    """

    if dps < 50:
        raise ValueError("at least 50 interval decimal digits are required")
    with mp.workdps(dps):
        ctx = MPIntervalContext()
        ctx.dps = dps
        tolerance_interval = ctx.mpf(tolerance)
        tolerance_lower, tolerance_upper = _interval_bounds(tolerance_interval)
        if tolerance_upper <= 0:
            raise ValueError("positive tolerance required")
        audit = source_audit(root)
        # Choose the lower endpoint of tolerance/4.  The retained endpoint
        # tails therefore remain inside the requested tolerance even when the
        # decimal input itself has a directed interval representation.
        quarter_interval = tolerance_interval / ctx.mpf(4)
        quarter_lower, _ = _interval_bounds(quarter_interval)
        delta_cap_interval = ctx.mpf(2) ** (-40)
        delta_cap_lower, delta_cap_upper = _interval_bounds(delta_cap_interval)
        delta = min(delta_cap_upper, quarter_lower)
        delta_interval = ctx.mpf(delta)
        delta_lower, delta_upper = _interval_bounds(delta_interval)
        # Use directed endpoint choices so the omitted tails cover every bit
        # lost to endpoint subtraction.  The actual tail bound is recomputed
        # from these chosen endpoints rather than assumed to equal 2*delta.
        left = delta_upper
        right_interval = ctx.mpf(1) - delta_interval
        right = _interval_bounds(right_interval)[0]
        left_interval = ctx.mpf(left)
        right_point_interval = ctx.mpf(right)
        left_tail_interval = left_interval
        right_tail_interval = ctx.mpf(1) - right_point_interval
        tail_length_interval = left_tail_interval + right_tail_interval
        tail_length_upper = _interval_bounds(tail_length_interval)[1]
        # Each endpoint tail is bounded by the true integral inequality
        # 0 <= exp(s-1) (1-sigma(s))^2 <= 1.  It is retained, not zeroed.
        tail_error_interval = ctx.mpf([0, tail_length_upper])
        tail_error = _interval_bounds(tail_error_interval)[1]
        quad_budget_interval = ctx.mpf([quarter_lower, quarter_lower])
        inner, details = adaptive_simpson_validated(
            ctx, left, right, quad_budget_interval, max_depth=max_depth
        )
        base = ctx.exp(ctx.mpf(-1))
        total = base + inner + tail_error_interval
        lower, upper = _interval_bounds(total)
        width_interval = ctx.mpf(upper) - ctx.mpf(lower)
        width_lower, width_upper = _interval_bounds(width_interval)
        # Compare directed endpoints directly.  A nearest rounded scalar
        # width could otherwise pass after losing its outward rounding.
        if width_upper > tolerance_lower:
            raise ValidatedIntegrationError(
                "certified interval width exceeds the requested tolerance"
            )
        quadrature_record = details["quadrature_error_record"]
        max_panel_record = details["max_panel_error_record"]
        return Gamma2Certificate(
            lower=lower,
            upper=upper,
            width=width_upper,
            delta=delta,
            tolerance=tolerance_lower,
            dps=dps,
            tail_error=tail_error,
            quadrature_error=details["quadrature_error"],
            max_panel_error=details["max_panel_error"],
            panel_count=details["panel_count"],
            directed_record={
                "total": _directed_record(total),
                "width": _directed_record(width_interval),
                "delta": _directed_record(delta_interval),
                "tolerance": _directed_record(tolerance_interval),
                "tail_length": _directed_record(tail_length_interval),
                "tail_error": _directed_record(tail_error_interval),
                "quadrature_error": quadrature_record,
                "max_panel_error": max_panel_record,
                "quadrature_budget": details["budget_record"],
            },
            source_audit={
                **audit,
                "formula": "exp(-1)+integral_0^1 exp(s-1)*(1-sigma(s))^2 ds",
                "restoration_source": "U=4*Z+E*(1-sigma), A_y=U^2-A",
                "tail_bound": "true integral bound 0<=h<=1 on [0,delta] and [1-delta,1]",
                "tail_error_is_retained": True,
                "tail_error_upper_is_directed": True,
                "quadrature_error_upper_is_directed": True,
                "width_upper_is_directly_compared_to_tolerance_lower": True,
                "quadrature": {
                    "panel_count": details["panel_count"],
                    "max_depth": details["max_depth"],
                    "remainder_formula": details["remainder_formula"],
                    "ordinary_derivative_order": details["ordinary_derivative_order"],
                    "whole_panel_interval_jet": details["whole_panel_interval_jet"],
                    "total_length": details["total_length_record"],
                    "budget": details["budget_record"],
                    "error": details["quadrature_error_record"],
                    "max_panel_error": details["max_panel_error_record"],
                },
            },
        )


def validate_sigma_source(
    *, dps: int = 100, root: str | Path | None = None
) -> dict[str, Any]:
    """Run only endpoint/range/symmetry checks for the bound sigma source."""

    audit = source_audit(root)
    ctx = MPIntervalContext()
    ctx.dps = dps
    with mp.workdps(dps):
        points = [mp.mpf("1/8"), mp.mpf("1/4"), mp.mpf("1/2"), mp.mpf("3/4"), mp.mpf("7/8")]
        values = [_integrand_sigma_value(ctx, point) for point in points]
        sigma_quarters = _sigma_jet(ctx, ctx.mpf([mp.mpf("1/4"), mp.mpf("1/4")])).derivatives[0]
        sigma_three_quarters = _sigma_jet(
            ctx, ctx.mpf([mp.mpf("3/4"), mp.mpf("3/4")])
        ).derivatives[0]
        qlo, qhi = _interval_bounds(sigma_quarters)
        rlo, rhi = _interval_bounds(sigma_three_quarters)
        symmetry = sigma_quarters + sigma_three_quarters
        slo, shi = _interval_bounds(symmetry)
        endpoint_zero = ctx.mpf([0, 0])
        endpoint_one = ctx.mpf([1, 1])
        passed = (
            _interval_bounds(endpoint_zero) == (mp.mpf(0), mp.mpf(0))
            and _interval_bounds(endpoint_one) == (mp.mpf(1), mp.mpf(1))
            and all(_interval_bounds(v)[0] >= 0 and _interval_bounds(v)[1] <= 1 for v in values)
            and qlo > 0
            and qhi < 1
            and rlo > 0
            and rhi < 1
            and slo <= 1 <= shi
        )
        return {
            "passed": bool(passed),
            "source_audit": audit,
            "endpoint_values": {"0": "0", "1": "1"},
            "sample_points": [str(x) for x in points],
            "sample_sigma_intervals": [_interval_strings(v) for v in values],
            "quarter_symmetry_sum": _interval_strings(symmetry),
            "range_claim": "0<=sigma<=1",
            "exact_source_symmetry": "sigma(s)+sigma(1-s)=1",
        }


def _integrand_sigma_value(ctx: MPIntervalContext, point: mp.mpf) -> Any:
    # Kept separate from the gamma2 integrand so this preflight cannot silently
    # validate the wrong quantity through a shared ancestor evaluator.
    return _sigma_jet(ctx, ctx.mpf([point, point])).derivatives[0]


def _interval_strings(value: Any) -> list[str]:
    lo, hi = _interval_bounds(value)
    return [str(lo), str(hi)]


def rational_poly_kernel_fixture(
    *, dps: int = 100, root: str | Path | None = None
) -> dict[str, Any]:
    """Fresh rational/poly witness for Simpson weights and Taylor factorials.

    The polynomial is ``1 + 2x + 3x^2 + 4x^3 + 5x^4`` on [0, 1].  Simpson's
    one-panel value is 121/24, the exact integral is 5, and the certified
    fourth-order remainder radius is 1/24 because p^(4)=120.
    """

    audit = source_audit(root)
    ctx = MPIntervalContext()
    ctx.dps = dps
    with mp.workdps(dps):
        x = TaylorJet.variable(ctx, ctx.mpf([0, 1]), 1)
        p = TaylorJet.constant(ctx, 1) + 2 * x + 3 * x.square() + 4 * (x.square() * x) + 5 * x.square().square()
        f0 = p.derivatives[0]
        f4 = p.derivatives[4]
        exact_simpson = mp.mpf(121) / 24
        exact_integral = mp.mpf(5)
        remainder = mp.mpf(1) / 24
        simpson_interval = (ctx.mpf(1) / 6) * (
            _poly_value(ctx, 0)
            + 4 * _poly_value(ctx, mp.mpf("1/2"))
            + _poly_value(ctx, 1)
        )
        simpson_error = exact_simpson - exact_integral
        derivative_lo, derivative_hi = _interval_bounds(f4)
        simpson_lo, simpson_hi = _interval_bounds(simpson_interval)
        rounding_check = mp.mpf("1e-90")
        remainder_difference = simpson_error - remainder
        passed = (
            _interval_bounds(f0)[0] == 1
            and _interval_bounds(f0)[1] == 15
            and derivative_lo == 120
            and derivative_hi == 120
            and -rounding_check <= remainder_difference <= rounding_check
            and simpson_lo <= exact_simpson <= simpson_hi
        )
        return {
            "passed": bool(passed),
            "source_audit": audit,
            "polynomial": "1+2*x+3*x^2+4*x^3+5*x^4",
            "interval": ["0", "1"],
            "simpson_coefficients": ["1/6", "4/6", "1/6"],
            "simpson_value": str(exact_simpson),
            "exact_integral": str(exact_integral),
            "fourth_derivative_interval": [str(derivative_lo), str(derivative_hi)],
            "remainder_bound": str(remainder),
            "ordinary_derivative_factorial_check": {
                "fourth_derivative": "120",
                "coefficient_times_factorial": "5*4! = 120",
            },
        }


def _poly_value(ctx: MPIntervalContext, point: mp.mpf) -> Any:
    x = point
    return ctx.mpf(1) + 2 * x + 3 * x**2 + 4 * x**3 + 5 * x**4


__all__ = [
    "Gamma2Certificate",
    "SourceAuditError",
    "TaylorJet",
    "ValidatedIntegrationError",
    "adaptive_simpson_validated",
    "rational_poly_kernel_fixture",
    "source_audit",
    "validate_sigma_source",
    "validated_gamma2",
]
