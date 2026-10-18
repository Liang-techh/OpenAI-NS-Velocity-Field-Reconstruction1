"""Analytic global bounds for the original flat transition ``sigma``.

The actual source is zero on ``x <= 0``, one on ``x >= 1``, and

    exp(-1/x**2) / (exp(-1/x**2) + exp(-1/(1-x)**2))

on ``0 < x < 1``.  This module does not evaluate that function or sample a
grid.  It pins the accepted source and returns integer majorants for the
absolute ordinary derivatives through order four.  The proof uses the
reflection ``sigma(1-x)=1-sigma(x)`` and the exponential flat-tail bound on
the half interval ``0 < x <= 1/2``.
"""
from __future__ import annotations

import ast
import hashlib
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping

import sympy as sp

from .original_finite_N_signed_source import require


BASE = "experiments/root_st073/"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_outer_initial.py"
SIGMA_WRAPPER = BASE + "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py"

# These are deliberately integer caps.  They are obtained from the displayed
# term table below, rather than from interval samples or an optimization call.
DERIVATIVE_UPPER_BOUNDS = (1, 1296, 38880, 19290960, 2251925712)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _normal(text: str) -> str:
    return "".join(text.split())


def _ast_guard(path: Path, function: str,
               fragments: tuple[str, ...]) -> dict[str, Any]:
    """Pin source syntax without importing or executing the source module."""
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    functions = [node for node in ast.walk(tree)
                 if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                 and node.name == function]
    require(len(functions) == 1, "Unique sigma source function required: " + function)
    node = functions[0]
    body = _normal(ast.get_source_segment(source, node) or "")
    checks = {_normal(fragment): _normal(fragment) in body
              for fragment in fragments}
    return {
        "passed": all(checks.values()),
        "path": str(path),
        "function": function,
        "required_source_fragments": checks,
        "function_AST_sha256": hashlib.sha256(
            ast.dump(node, include_attributes=False).encode()).hexdigest(),
        "source_executed": False,
    }


def _pin(root: Path, parent: Mapping[str, Any], relative: str) -> dict[str, Any]:
    inputs = parent.get("input_hashes")
    require(isinstance(inputs, Mapping), "Parent input_hashes are required")
    key = Path(relative).name
    expected = inputs.get(relative) or inputs.get(key)
    require(isinstance(expected, str), "Missing parent source pin: " + key)
    actual = _sha(root / relative)
    require(actual == expected, "Changed pinned sigma source: " + relative)
    return {
        "path": relative,
        "key": key,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": True,
        "source_executed": False,
    }


def _symbolic_controls() -> dict[str, Any]:
    """Verify the logistic derivative identities by formal differentiation."""
    s = sp.symbols("s")
    l1, l2, l3, l4 = sp.symbols("L1 L2 L3 L4")
    q = s * (1 - s)
    a = 1 - 2 * s

    def total_derivative(expr: sp.Expr) -> sp.Expr:
        variables = (s, l1, l2, l3, l4)
        rates = (q * l1, l2, l3, l4, 0)
        return sp.expand(sum(sp.diff(expr, var) * rate
                             for var, rate in zip(variables, rates)))

    formulas = [
        s,
        q * l1,
        q * (l2 + a * l1**2),
        q * (l3 + 3 * a * l1 * l2 + (1 - 6 * q) * l1**3),
        q * (l4 + a * (4 * l1 * l3 + 3 * l2**2)
             + 6 * (1 - 6 * q) * l1**2 * l2
             + a * (1 - 12 * q) * l1**4),
    ]
    generated = [s]
    for _ in range(4):
        generated.append(total_derivative(generated[-1]))
    checks = [sp.simplify(generated[k] - formulas[k]) == 0
              for k in range(5)]
    return {
        "passed": all(checks),
        "checks": {f"order_{k}": checks[k] for k in range(5)},
        "ordinary_derivatives": True,
        "formulas": [str(expr) for expr in formulas],
        "q_definition": "s*(1-s)",
        "a_definition": "1-2*s",
    }


def _power_cap(power: Fraction) -> int:
    """Integer cap for exp(4-t)t**power, t>=4.

    Since e < 3, exp(4) < 81.  Since the maximum of t**p exp(-t) on
    t>0 is p**p, and p <= ceil(p), the returned integer is valid.
    """
    ceiling = (power.numerator + power.denominator - 1) // power.denominator
    return 81 * ceiling ** ceiling


def _term_table() -> tuple[dict[str, Any], ...]:
    # On m=min(x,1-x)<=1/2, |L^(k)| <= c_k*m^(-(k+2)) with
    # c=(4,12,48,240).  Each row is one term in the absolute logistic
    # derivative formula, after replacing q by exp(4-m^-2).
    rows = (
        (1, ((4, Fraction(3, 2)),)),
        (2, ((12, Fraction(2)), (16, Fraction(3)))),
        (3, ((48, Fraction(5, 2)), (144, Fraction(7, 2)),
             (64, Fraction(9, 2)))),
        (4, ((240, Fraction(3)), (768, Fraction(4)), (432, Fraction(4)),
             (1152, Fraction(5)), (512, Fraction(6)))),
    )
    output = []
    for order, terms in rows:
        term_rows = []
        total = 0
        for coefficient, power in terms:
            cap = _power_cap(power)
            contribution = coefficient * cap
            total += contribution
            term_rows.append({
                "coefficient": coefficient,
                "m_power": str(power),
                "exp_power_cap": cap,
                "contribution_cap": contribution,
            })
        output.append({"order": order, "terms": term_rows,
                       "sum_cap": total,
                       "selected_integer_cap": DERIVATIVE_UPPER_BOUNDS[order],
                       "cap_dominates_sum": DERIVATIVE_UPPER_BOUNDS[order] >= total})
    return tuple(output)


def _analytic_proof() -> dict[str, Any]:
    table = _term_table()
    require(all(row["cap_dominates_sum"] for row in table),
            "Analytic sigma cap table does not close")
    return {
        "passed": True,
        "derivative_upper_bounds": list(DERIVATIVE_UPPER_BOUNDS),
        "domain_split": {
            "left": "0<x<=1/2",
            "right": "1/2<=x<1 via sigma(1-x)=1-sigma(x)",
            "outside": "sigma is constant on x<=0 and x>=1",
        },
        "logistic_log_ratio": (
            "L=1/(1-x)^2-1/x^2; sigma=exp(L)/(1+exp(L)) on the left"),
        "left_tail_bound": "q=sigma(1-sigma)<=exp(4-m^-2), m=x<=1/2",
        "log_derivative_bounds": {
            "L1": "4*m^-3",
            "L2": "12*m^-4",
            "L3": "48*m^-5",
            "L4": "240*m^-6",
        },
        "power_majorant": (
            "for t=m^-2>=4 and p>0: exp(4-t)t^p <= 81*ceil(p)^ceil(p)"),
        "term_table": list(table),
        "reflection_identity": (
            "sigma(1-x)=1-sigma(x); |sigma^(k)(1-x)|=|sigma^(k)(x)|"),
        "seam_flatness": (
            "m^-n*exp(-m^-2) tends to zero for every n<=12; hence the "
            "interior derivatives through order four match the constant "
            "outside derivatives at x=0 and x=1"),
        "no_sampling_or_supremum": True,
        "no_numerical_source_execution": True,
    }


def build(root: Any, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Return a source-pinned proof of global sigma derivative caps."""
    root = Path(root)
    require(isinstance(parent, Mapping), "Parent mapping is required")
    pins = {
        "stable_sigma_source": _pin(root, parent, SIGMA_SOURCE),
        "Gamma_sigma_wrapper": _pin(root, parent, SIGMA_WRAPPER),
    }
    guards = {
        "stable_sigma": _ast_guard(
            root / SIGMA_SOURCE, "stable_sigma",
            ("if v<=0:return c.mpf(0),c.mpf(0)",
             "if v>=1:return c.mpf(1),c.mpf(0)",
             "a=c.exp(-1/x**2);b=c.exp(-1/(1-x)**2)",
             "value=a/(a+b)",
             "derivative=value*(b/(a+b))*(2/x**3+2/(1-x)**3)")),
        "Gamma_sigma_binding": _ast_guard(
            root / SIGMA_WRAPPER, "sigma",
            ("return g.node('exact_original_source_flat_sigma', argument=argument.node",
             "definition='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside'",
             "source_binding=parent.active.current.ast_binding(flat_source.stable_sigma)",
             "exact_function_not_its_interval_cap=True")),
    }
    require(all(guard["passed"] for guard in guards.values()),
            "Pinned sigma source AST guard failed")
    symbolic = _symbolic_controls()
    require(symbolic["passed"], "Symbolic logistic derivative controls failed")
    analytic = _analytic_proof()
    return {
        "passed": True,
        "source_only": True,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "source_pins": pins,
        "source_pins_ok": True,
        "AST": guards,
        "symbolic_controls": symbolic,
        "analytic_proof": analytic,
        "derivative_upper_bounds": list(DERIVATIVE_UPPER_BOUNDS),
        "ordinary_derivative_convention": "absolute D_x^k sigma, no factorial divisor",
        "caps_are_global_scalar_majorants": True,
        "numeric_sampling": False,
        "numerical_supremum": False,
    }


__all__ = ["build"]
