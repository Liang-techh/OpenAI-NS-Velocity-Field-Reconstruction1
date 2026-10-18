"""Source-bound normalized axial pressure jets on the current real Z tube.

This module consumes the accepted fourteen-stage pressure archive only.  It
does not instantiate ``CompliantPressureDatum`` or execute any source
producer.  The returned list contains ordinary (unfactorialized) derivative
intervals for ``p0=P0/Pstar**2`` over the whole interval
``Z in [37/100, 93/250]``.
"""
from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext

from .original_finite_N_signed_source import require
from .original_ordinary_jet_bounds import jet_mul, jet_reciprocal


BASE = "experiments/root_st073/"
PRESSURE_JSON = BASE + "lei_ren_part1_paper_compliant_pressure_source.json"
PRESSURE_SOURCE = "lei_ren_part1_paper_compliant_pressure_source.py"
DATUM_SOURCE = "lei_ren_part1_paper_logarithmic_pressure_datum.py"
Q_SOURCE = "lei_ren_part1_paper_candidate_pressure_function.py"
CORE_NUMERIC = "src/openai_ns_reconstruction/original_core_numeric.py"

Z_LO = "37/100"
Z_HI = "93/250"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _interval(ctx: MPIntervalContext, value: Any) -> Any:
    """Decode an archive interval without going through decimal endpoints."""
    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, Mapping):
        exact = value.get("source_enclosure_exact_mpf_tuples")
        if isinstance(exact, (list, tuple)) and len(exact) == 2:
            lo, hi = exact
        else:
            lo = value.get("lower_exact_mpf_tuple")
            hi = value.get("upper_exact_mpf_tuple")
        if lo is not None and hi is not None:
            return ctx.mpf([
                mp.make_mpf(tuple(int(v) for v in lo)),
                mp.make_mpf(tuple(int(v) for v in hi)),
            ])
        if "lower" in value and "upper" in value:
            return ctx.mpf([str(value["lower"]), str(value["upper"])])
        if "upper" in value:
            return ctx.mpf(str(value["upper"]))
    return ctx.mpf(str(value))


def _record(value: Any) -> dict[str, Any]:
    lo = mp.make_mpf(value._mpi_[0])
    hi = mp.make_mpf(value._mpi_[1])
    return {
        "lower_exact_mpf_tuple": list(map(int, lo._mpf_)),
        "upper_exact_mpf_tuple": list(map(int, hi._mpf_)),
        "lower_decimal_display": mp.nstr(lo, 28),
        "upper_decimal_display": mp.nstr(hi, 28),
        "directed_interval": True,
    }


def _ast_guard(path: Path, function: str, fragments: tuple[str, ...],
               *, class_name: str | None = None) -> dict[str, Any]:
    """Check source syntax without importing or executing the source module."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    scope: ast.AST = tree
    if class_name is not None:
        classes = [node for node in tree.body
                   if isinstance(node, ast.ClassDef) and node.name == class_name]
        require(len(classes) == 1, "Unique source class required: " + class_name)
        scope = classes[0]
    functions = [node for node in ast.walk(scope)
                 if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                 and node.name == function]
    require(len(functions) == 1, "Unique source function required: " + function)
    node = functions[0]
    nodes = {ast.dump(part, include_attributes=False) for part in ast.walk(node)}
    checks: dict[str, bool] = {}
    for fragment in fragments:
        wanted = ast.parse(fragment).body[0]
        if isinstance(wanted, ast.Expr):
            wanted = wanted.value
        checks[fragment] = ast.dump(wanted, include_attributes=False) in nodes
    return {
        "path": str(path),
        "class": class_name,
        "function": function,
        "required_AST": checks,
        "function_AST_sha256": hashlib.sha256(
            ast.dump(node, include_attributes=False).encode()
        ).hexdigest(),
        "source_executed": False,
        "passed": all(checks.values()),
    }


def _pin(root: Path, parent: Mapping[str, Any], relative: str,
         *, parent_section: str = "input_hashes") -> dict[str, Any]:
    """Bind a source file to the accepted parent report hash."""
    section = parent.get(parent_section, {})
    key = relative[len(BASE):] if relative.startswith(BASE) else relative
    expected = section.get(key)
    require(isinstance(expected, str), "Missing accepted source hash: " + key)
    actual = _sha(root / relative)
    require(actual == expected, "Changed accepted source: " + relative)
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "bound_to_parent": True,
    }


def _symmetric_bound(ctx: MPIntervalContext, bound: Any) -> Any:
    """Return the symmetric interval [-upper(bound), upper(bound)]."""
    upper = mp.make_mpf(bound._mpi_[1])
    return ctx.mpf([-upper, upper])


def _q_jet(ctx: MPIntervalContext, z: Any, order: int) -> list[Any]:
    """Compute ordinary derivatives of q=(1+Z**2)^(-2) by jet arithmetic."""
    zero = ctx.mpf(0)
    one = ctx.mpf(1)
    zjet = [z] + ([one] if order >= 1 else []) + [zero] * max(order - 1, 0)
    z2 = jet_mul(ctx, zjet, zjet)
    denominator = [one + z2[0]] + z2[1:]
    inverse = jet_reciprocal(ctx, denominator)
    return jet_mul(ctx, inverse, inverse)


def _flat_jet(ctx: MPIntervalContext, source: Mapping[str, Any], order: int) -> tuple[list[Any], Any, Any]:
    """Build the saved F_flat remainder using ordinary Cauchy derivatives."""
    stages = source.get("stages")
    require(isinstance(stages, Mapping) and isinstance(stages.get("z_flatten"), Mapping),
            "Saved z_flatten stage is required")
    rho = _interval(ctx, source["complex_strip_half_width"])
    flat_upper = _interval(ctx, source["flatten_complex_mass_upper"])
    flat_zero = _interval(ctx, stages["z_flatten"]["mass"])
    flat_zero_lower = mp.make_mpf(flat_zero._mpi_[0])
    flat_upper_lower = mp.make_mpf(flat_upper._mpi_[0])
    rho_lower = mp.make_mpf(rho._mpi_[0])
    require(flat_zero_lower >= 0 and flat_upper_lower > 0 and rho_lower > 0,
            "Flatten mass bounds must be nonnegative")
    rows = [flat_zero]
    for k in range(1, order + 1):
        rows.append(_symmetric_bound(ctx, math.factorial(k) * flat_upper / rho**k))
    return rows, rho, flat_upper


def build(c: MPIntervalContext, root: Any, parent: Mapping[str, Any], order: int = 5):
    """Return ``(p0jets, proof, pins)`` for normalized axial pressure.

    ``p0jets[k]`` is a directed interval enclosing the ordinary derivative
    ``D_Z**k(P0/Pstar**2)`` at every real ``Z`` in the requested tube.  The
    flattening source is retained as a signed remainder at order zero and as
    a symmetric Cauchy derivative remainder at higher orders.
    """
    require(isinstance(c, MPIntervalContext), "An MPIntervalContext is required")
    require(isinstance(parent, Mapping), "Accepted parent mapping is required")
    require(type(order) is int and not isinstance(order, bool) and order >= 0,
            "Nonnegative ordinary jet order required")
    root = Path(root)

    family = parent.get("source_family")
    require(isinstance(family, Mapping), "Accepted parent source family is required")
    require(parent.get("checks", {}).get("passed") is True,
            "Passed accepted parent report required")

    pins = {"pressure_source_json": _pin(root, parent, PRESSURE_JSON)}
    for label, name in (("pressure_source", PRESSURE_SOURCE),
                        ("datum_source", DATUM_SOURCE),
                        ("q_source", Q_SOURCE)):
        pins[label] = _pin(root, parent, BASE + name)

    # The scale AST is a source-only guard.  The normalized output below never
    # materializes Pstar or exp(5); it records the exact logPstar formula.
    core_expected = parent.get("implementation_hashes", {}).get(CORE_NUMERIC)
    require(isinstance(core_expected, str), "Accepted core scale source hash required")
    core_path = root / CORE_NUMERIC
    require(_sha(core_path) == core_expected, "Changed accepted core scale source")
    pins["core_scale_source"] = {
        "path": CORE_NUMERIC,
        "sha256": core_expected,
        "expected_sha256": core_expected,
        "bound_to_parent": True,
    }

    pressure_path = root / PRESSURE_JSON
    report = _read(pressure_path)
    source = report.get("compliant_source")
    require(isinstance(source, Mapping), "compliant_source archive is required")
    for key in ("implicit_source_sha256", "datum_enclosure_sha256"):
        require(source.get(key) == family.get(key), "Pressure source family mismatch: " + key)
    require(source.get("stage_count_total") == 14 and
            source.get("global_real_axis_jets_supported") is True and
            source.get("complete_preheat_pressure_enclosed") is True,
            "Complete global pressure source receipt required")

    nested_hashes = source.get("input_hashes")
    require(isinstance(nested_hashes, Mapping), "Pressure nested input hashes are required")
    for name in (PRESSURE_SOURCE, DATUM_SOURCE, Q_SOURCE):
        require(nested_hashes.get(name) == pins[
            {PRESSURE_SOURCE: "pressure_source", DATUM_SOURCE: "datum_source",
             Q_SOURCE: "q_source"}[name]
        ]["sha256"], "Nested pressure source hash mismatch: " + name)

    guards = {
        "compliant_constructor": _ast_guard(
            root / (BASE + PRESSURE_SOURCE), "__init__",
            ("self.rho = c.mpf('.25')",
             "self.stages['z_flatten'] = dict(beta='variable [0,2]', mass=box)",
             "self.flatten_complex_upper = self.tail_upper/(1-self.rho**2)**2"),
            class_name="CompliantPressureDatum"),
        "normalized_jets_cauchy": _ast_guard(
            root / (BASE + DATUM_SOURCE), "normalized_jets",
            ("remainder=self.stages['z_flatten']['mass']",
             "bound=self.flatten_complex_upper/self.rho**n",
             "remainder=c.mpf([endpoints(-bound)[0],endpoints(bound)[1]])",
             "rows.append(-(self.m2*a+(self.m0 if n==0 else 0)+remainder))"),
            class_name="LogarithmicPressureDatum"),
        "q_source": _ast_guard(
            root / (BASE + Q_SOURCE), "q_jets",
            ("a = 1 + center**2", "b = 2*center",
             "coefficients = [1/a**2]",
             "coefficients.append(-(b*(n+2)*current+(n+3)*previous)/(a*(n+1)))")),
        "core_logPstar": _ast_guard(
            core_path, "__init__", ("self.logP=c.exp(40)+11",),
            class_name="OriginalCoreNumeric"),
        "core_normalized_pressure_formula": _ast_guard(
            core_path, "_seed",
            ("pressure.append(-(M2*q+(M0 if k == 0 else 0)+remainder))",)),
    }
    require(all(item["passed"] for item in guards.values()),
            "Pinned pressure source AST guard failed")

    z = c.mpf([c.mpf(37) / 100, c.mpf(93) / 250])
    q = _q_jet(c, z, order)
    flat, rho, flat_upper = _flat_jet(c, source, order)
    m2 = _interval(c, source["fixed_beta2_mass_normalized"])
    m0 = _interval(c, source["fixed_beta0_mass_normalized"])
    require(mp.make_mpf(m2._mpi_[0]) >= 0 and mp.make_mpf(m0._mpi_[0]) >= 0,
            "Saved normalized pressure masses must be nonnegative")
    zero = c.mpf(0)
    p0jets = [-(m2 * q[k] + flat[k] + (m0 if k == 0 else zero))
              for k in range(order + 1)]
    controls = _ordinary_jet_controls(c)
    require(controls["passed"], "Exact ordinary jet controls failed")

    proof = {
        "passed": True,
        "source_family": dict(family),
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "domain": {"Z_lower": Z_LO, "Z_upper": Z_HI,
                   "exact_real_interval": True},
        "ordinary_Z_orders": list(range(order + 1)),
        "ordinary_derivative_convention": "D_Z^k; no factorial divisor",
        "normalized_units": "P0/Pstar^2",
        "defining_function": "p0(Z)=-(M2*(1+Z^2)^(-2)+M0+F_flat(Z))",
        "M2": _record(m2),
        "M0": _record(m0),
        "F_flat_order_zero_remainder": _record(flat[0]),
        "F_flat_higher_order_remainders": [_record(v) for v in flat[1:]],
        "flatten_complex_mass_upper": _record(flat_upper),
        "complex_strip_half_width": _record(rho),
        "q_jets": [_record(v) for v in q],
        "p0_jets": [_record(v) for v in p0jets],
        "ordinary_jet_controls": controls,
        "cauchy_derivative_bound": "k!*flatten_complex_mass_upper/rho^k for k>=1",
        "source_family_hashes": {
            "implicit_source_sha256": source["implicit_source_sha256"],
            "datum_enclosure_sha256": source["datum_enclosure_sha256"],
        },
        "physical_scale": {
            "logPstar_formula": "exp(40)+11",
            "Pstar_squared_not_applied": True,
            "physical_P0_formula": "P0=Pstar^2*p0",
        },
        "source_AST": guards,
        "numeric_cap_not_promoted_to_value": True,
        "flatten_remainder_retained": True,
    }
    return p0jets, proof, pins


def _ordinary_jet_controls(ctx: MPIntervalContext) -> dict[str, Any]:
    """Exact rational controls for ordinary, unfactorialized jets through Z5."""
    x = sp.symbols("x")
    f = (sp.Rational(7, 5) + 2*x - sp.Rational(5, 3)*x**2
         + sp.Rational(7, 4)*x**3 - sp.Rational(2, 5)*x**4
         + sp.Rational(3, 8)*x**5)
    g = (sp.Rational(5, 3) - x + sp.Rational(4, 7)*x**2
         + sp.Rational(3, 2)*x**3 - sp.Rational(1, 3)*x**4
         + sp.Rational(2, 9)*x**5)
    order = 5
    derivative = lambda h: [sp.diff(h, x, k).subs(x, 0) for k in range(order + 1)]
    F, G = derivative(f), derivative(g)
    product = [sum((sp.binomial(k, j) * F[j] * G[k-j]
                    for j in range(k + 1)), sp.Integer(0))
               for k in range(order + 1)]
    reciprocal = [1 / F[0]]
    for k in range(1, order + 1):
        reciprocal.append(-sum((sp.binomial(k, j) * F[j] * reciprocal[k-j]
                                for j in range(1, k + 1)), sp.Integer(0)) / F[0])
    expected_product = derivative(f * g)
    expected_reciprocal = derivative(1 / f)
    # Exercise the actual q implementation as an algorithm-level guard.  A
    # generic product/reciprocal control does not detect composing the
    # reciprocal operation the wrong number of times.  At Z=0 the ordinary
    # derivatives of (1+Z**2)**(-2) are [1, 0, -4, 0, 72, 0].
    q_at_zero = _q_jet(ctx, ctx.mpf(0), order)
    q_expected = [1, 0, -4, 0, 72, 0]
    q_algorithm = all(
        mp.make_mpf(value._mpi_[0]) == mp.mpf(expected) and
        mp.make_mpf(value._mpi_[1]) == mp.mpf(expected)
        for value, expected in zip(q_at_zero, q_expected)
    )
    checks = {
        "mul_binomial": product == expected_product,
        "reciprocal": all(sp.simplify(reciprocal[k] - expected_reciprocal[k]) == 0
                           for k in range(order + 1)),
        "ordinary_factorial": derivative(x**2)[2] == 2,
        "q_jet_algorithm": q_algorithm,
    }
    return {
        "passed": all(checks.values()),
        "order": order,
        "factorial_normalized": False,
        "checks": checks,
        "q_at_zero_expected": q_expected,
        "q_at_zero_exact": q_algorithm,
        "no_source_claim": True,
    }


__all__ = ["build"]
