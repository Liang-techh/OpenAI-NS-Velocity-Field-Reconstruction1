"""Fresh synthetic rejection fixture for the v2 Green scale source leaf.

This fixture deliberately does not call the v1 fixture and does not consume
any caller-provided repository or source text.  It builds three new AST-only
inputs: a changed first row tuple, a wrong Volterra ``FIELDS`` binding, and a
Z-dependent first-chart boundary.  Each case must be rejected by the supplied
v2 text auditor; no canonical or v1 mutation case is re-entered here.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping
from pathlib import Path
from typing import Any


SCHEMA = "n1_green_scale_commutator_fixture_v2"
GOOD_ROWS = "('5','6','6','4','5','6')"
BAD_FIRST_ROW_ROWS = "('6','6','6','4','5','6')"
GOOD_FIELDS_BINDING = "FIELDS=regular.FIELDS"
WRONG_FIELDS_BINDING = "FIELDS=regular.WRONG"
GOOD_FIRST_BOUNDARY = "g.zero"
WRONG_FIRST_BOUNDARY = "g.symbol('Z')"


def _parse(text: str, label: str) -> None:
    if not isinstance(text, str):
        raise TypeError(f"{label} must be source text")
    try:
        ast.parse(text, filename=label, mode="exec")
    except SyntaxError as exc:
        raise ValueError(f"v2 synthetic source failed AST parse: {label}: {exc}") from exc


def _synthetic_texts(
    rows: str,
    fields_binding: str,
    first_boundary: str,
    first_rows: str | None = None,
) -> tuple[str, str, str]:
    """Build one independent v2 source triple without reading any file."""

    first_rows = rows if first_rows is None else first_rows
    volterra = '''
import synthetic_regular as regular
FIELDS_BINDING_SENTINEL
def green_apply(g,vectors,chart,coordinates,parameters):
    c=g.constant;a=c3.C3Algebra(g);hb=parameters['hb'];r=coordinates['core'];s=coordinates['first']
    t=g.symbol('v2_green_t');v=g.symbol('v2_green_phase_v')
    integral=lambda q,var,lo,hi: c3.C3Function(*(controls.integral(g,row,g.nodes[var.node]['name'],lo,hi,
        actual_raw_regular_n1_Green_body=True,ordinary_slow_Z_derivative_order=k,
        endpoint_independent_of_Z=True,zero_endpoint_integrand_continuation=0,
        complete_history_not_reset_at_core_seam=True,signed_body_not_magnitude_cap=True)
        for k,row in enumerate(c3.rows(q))))
    if chart=='core':
        X=g.unary('sqrt',r);sample=g.mul(r,t,t);measure=X
        K=g.mul(c('1/2'),X,t,g.sub(g.one,g.mul(t,t)))
        log=g.mul(X,t,g.unary('log',g.quotient(g.one,t,'positive t inside unit integral')))
        kernels=(K,log,g.neg(K),g.one,g.mul(t,t,t),t)
    elif chart=='first':
        X=g.mul(c(2),g.unary('exp',g.mul(c('1/2'),hb,s)));sample=g.mul(c(4),t,t);measure=c(2)
        distance=g.mul(hb,s);tail=g.mul(distance,g.unary('exprel',g.neg(distance)))
        K=g.mul(t,g.add(g.sub(g.one,g.mul(t,t)),g.mul(t,t,tail)))
        log=g.mul(c(2),t,g.add(g.mul(c('1/2'),hb,s),g.unary('log',g.quotient(g.one,t,'positive t inside unit integral'))))
        ratio=g.mul(t,g.unary('exp',g.neg(g.mul(c('1/2'),hb,s))))
        kernels=(K,log,g.neg(K),g.one,g.mul(ratio,ratio,ratio),ratio)
    else:raise ValueError('Regular n1 owns only core/first charts')
    out=[]
    for row,kernel in zip(CORE_ROWS_SENTINEL,kernels):
        body=collar.at(g,vectors['core'][row],r,sample)
        out.append(integral(a.scale(body,g.mul(kernel,measure)),t,g.zero,g.one))
    if chart=='first':
        x=g.mul(c(2),g.unary('exp',g.mul(c('1/2'),hb,v)))
        distance=g.mul(hb,g.sub(s,v));ratio=g.unary('exp',g.neg(g.mul(c('1/2'),distance)))
        K=g.mul(c('1/2'),x,distance,g.unary('exprel',g.neg(distance)))
        log=g.mul(c('1/2'),x,distance);measure=g.mul(c('1/2'),hb,x)
        kernels=(K,log,g.neg(K),g.one,g.mul(ratio,ratio,ratio),ratio)
        for i,(row,kernel) in enumerate(zip(FIRST_ROWS_SENTINEL,kernels)):
            body=collar.at(g,vectors['first'][row],s,v)
            out[i]=a.add(out[i],integral(a.scale(body,g.mul(kernel,measure)),v,FIRST_BOUNDARY_SENTINEL,s))
    return dict(zip(FIELDS,out))
'''
    volterra = volterra.replace("FIELDS_BINDING_SENTINEL", fields_binding)
    volterra = volterra.replace("CORE_ROWS_SENTINEL", rows)
    volterra = volterra.replace("FIRST_ROWS_SENTINEL", first_rows)
    volterra = volterra.replace("FIRST_BOUNDARY_SENTINEL", first_boundary)
    regular = """
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
"""
    first_raw = """
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
CANONICAL_STATE = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
SCALE_DIAGONAL = (
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    "epsilon",
)
"""
    return volterra, regular, first_raw


def _auditor(auditor: Any) -> Any:
    candidate = getattr(auditor, "audit_source_texts", None)
    if not callable(candidate):
        candidate = auditor
    if callable(candidate) and getattr(candidate, "__name__", "") == "build_source_bound_proof":
        candidate = getattr(candidate, "__globals__", {}).get("audit_source_texts")
    if not callable(candidate):
        raise TypeError("auditor must be callable or expose audit_source_texts")
    return candidate


def _invoke_case(auditor: Any, texts: tuple[str, str, str]) -> Mapping[str, Any] | None:
    try:
        for label, text in zip(("volterra_v2", "regular_v2", "first_raw_v2"), texts):
            _parse(text, label)
        candidate = _auditor(auditor)
        try:
            result = candidate(*texts, parent=None)
        except TypeError as first_error:
            try:
                result = candidate(*texts)
            except TypeError:
                raise first_error
        return result if isinstance(result, Mapping) else None
    except (TypeError, ValueError, KeyError, OSError):
        return None


def _checks(result: Mapping[str, Any] | None) -> Mapping[str, Any]:
    if not isinstance(result, Mapping):
        return {}
    checks = result.get("checks")
    return checks if isinstance(checks, Mapping) else {}


def _rejected(result: Mapping[str, Any] | None, required_false: str) -> bool:
    checks = _checks(result)
    return (
        isinstance(result, Mapping)
        and result.get("passed") is False
        and checks.get(required_false) is False
    )


def _unexpected_false_keys(result: Mapping[str, Any] | None) -> set[str]:
    """Return structural false keys, excluding provenance-only gates."""

    checks = _checks(result)
    return {
        str(key)
        for key, value in checks.items()
        if value is False and not str(key).startswith(("parent_report_", "source_pin_"))
    }


def run_fixture(
    auditor: Any,
    repo: str | Path | None = None,
    parent: Mapping[str, Any] | str | Path | None = None,
    source_texts: Mapping[str, str] | tuple[str, str, str] | list[str] | None = None,
    **_: Any,
) -> dict[str, Any]:
    """Run only three new v2 synthetic rejection cases."""

    # These compatibility arguments are intentionally not consumed.  The
    # actual source-only proof is a later, independent runner invocation.
    del repo, parent, source_texts

    row_case = _synthetic_texts(
        BAD_FIRST_ROW_ROWS,
        GOOD_FIELDS_BINDING,
        GOOD_FIRST_BOUNDARY,
        first_rows=GOOD_ROWS,
    )
    fields_case = _synthetic_texts(
        GOOD_ROWS,
        WRONG_FIELDS_BINDING,
        GOOD_FIRST_BOUNDARY,
    )
    boundary_case = _synthetic_texts(
        GOOD_ROWS,
        GOOD_FIELDS_BINDING,
        WRONG_FIRST_BOUNDARY,
    )

    row_result = _invoke_case(auditor, row_case)
    fields_result = _invoke_case(auditor, fields_case)
    boundary_result = _invoke_case(auditor, boundary_case)
    checks = {
        "new_row_tuple_case_rejected": _rejected(
            row_result,
            "green_exact_rows_present_in_both_for_loops",
        ),
        "new_fields_binding_case_rejected": _rejected(
            fields_result,
            "green_green_FIELDS_module_bound_to_regular_FIELDS",
        ),
        "new_first_boundary_case_rejected": _rejected(
            boundary_result,
            "green_first_loop_out_index_add_accumulation",
        ),
    }
    false_by_case = {
        "row_tuple": _unexpected_false_keys(row_result),
        "fields_binding": _unexpected_false_keys(fields_result),
        "first_boundary": _unexpected_false_keys(boundary_result),
    }
    shared_unexpected_false_checks = sorted(
        false_by_case["row_tuple"]
        & false_by_case["fields_binding"]
        & false_by_case["first_boundary"]
    )
    checks["no_shared_unexpected_false_checks"] = not shared_unexpected_false_checks
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "shared_unexpected_false_checks": shared_unexpected_false_checks,
        "unexpected_false_checks_by_case": {
            key: sorted(value) for key, value in false_by_case.items()
        },
        "synthetic_only": True,
        "source_texts_only": True,
        "source_texts_consumed": False,
        "actual_source_called_here": False,
        "v1_fixture_called": False,
        "v1_canonical_case_called": False,
        "v1_old_mutation_cases_called": False,
        "implementation_imported": False,
        "old_modules_called": False,
        "provider_executed": False,
        "numerical_residual_evaluated": False,
        "cases": {
            "first_row_tuple": BAD_FIRST_ROW_ROWS,
            "fields_binding": WRONG_FIELDS_BINDING,
            "first_chart_boundary": WRONG_FIRST_BOUNDARY,
        },
        "scope": {
            "only_rejection_cases": True,
            "actual_source_proof_deferred": True,
            "old_v1_samples_reentered": False,
        },
    }


def run(auditor: Any, repo: str | Path | None = None, parent: Mapping[str, Any] | str | Path | None = None, **kwargs: Any) -> dict[str, Any]:
    return run_fixture(auditor, repo=repo, parent=parent, **kwargs)


__all__ = ["SCHEMA", "run_fixture", "run"]
