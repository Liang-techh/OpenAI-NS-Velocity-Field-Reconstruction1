"""Fresh v3 synthetic guards for the corrected Green first-loop target check.

Only two new negative source triples are sent to the supplied v3 text
auditor.  One changes the first-chart accumulation target to ``out[i+1]``;
the other keeps that target correct and changes scale column 6.  No v1/v2
fixture input, canonical case, or historical mutation is re-entered.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping
from pathlib import Path
from typing import Any


SCHEMA = "n1_green_scale_commutator_fixture_v3"
GOOD_OUT_TARGET = "out[i]"
WRONG_OUT_TARGET = "out[i+1]"
GOOD_SCALE_COLUMN_6 = '"epsilon",'
WRONG_SCALE_COLUMN_6 = '"epsilon*column6_wrong",'
EXPECTED_TARGET_FALSE = "green_first_loop_out_index_add_accumulation"
EXPECTED_SCALE_FALSE = "first_raw_SCALE_DIAGONAL_module_literal"


def _parse(text: str, label: str) -> None:
    if not isinstance(text, str):
        raise TypeError(f"{label} must be source text")
    try:
        ast.parse(text, filename=label, mode="exec")
    except SyntaxError as exc:
        raise ValueError(f"v3 synthetic source failed AST parse: {label}: {exc}") from exc


def _synthetic_texts(out_target: str, scale_column_6: str) -> tuple[str, str, str]:
    """Build a new source triple without reading any repository text."""

    volterra = '''
import synthetic_regular_v3 as regular
FIELDS=regular.FIELDS
def green_apply(g,vectors,chart,coordinates,parameters):
    c=g.constant;a=c3.C3Algebra(g);hb=parameters['hb'];r=coordinates['core'];s=coordinates['first']
    t=g.symbol('v3_green_t');v=g.symbol('v3_green_phase_v')
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
    for row,kernel in zip(('5','6','6','4','5','6'),kernels):
        body=collar.at(g,vectors['core'][row],r,sample)
        out.append(integral(a.scale(body,g.mul(kernel,measure)),t,g.zero,g.one))
    if chart=='first':
        x=g.mul(c(2),g.unary('exp',g.mul(c('1/2'),hb,v)))
        distance=g.mul(hb,g.sub(s,v));ratio=g.unary('exp',g.neg(g.mul(c('1/2'),distance)))
        K=g.mul(c('1/2'),x,distance,g.unary('exprel',g.neg(distance)))
        log=g.mul(c('1/2'),x,distance);measure=g.mul(c('1/2'),hb,x)
        kernels=(K,log,g.neg(K),g.one,g.mul(ratio,ratio,ratio),ratio)
        for i,(row,kernel) in enumerate(zip(('5','6','6','4','5','6'),kernels)):
            body=collar.at(g,vectors['first'][row],s,v)
            OUT_TARGET_SENTINEL=a.add(out[i],integral(a.scale(body,g.mul(kernel,measure)),v,g.zero,s))
    return dict(zip(FIELDS,out))
'''
    volterra = volterra.replace("OUT_TARGET_SENTINEL", out_target)
    first_raw = '''
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
CANONICAL_STATE = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
SCALE_DIAGONAL = (
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    SCALE_COLUMN_6_SENTINEL
)
'''.replace("SCALE_COLUMN_6_SENTINEL", scale_column_6)
    regular = """
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
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
        for label, text in zip(("volterra_v3", "regular_v3", "first_raw_v3"), texts):
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


def _rejected(result: Mapping[str, Any] | None, *required_false: str) -> bool:
    checks = _checks(result)
    return (
        isinstance(result, Mapping)
        and result.get("passed") is False
        and any(checks.get(key) is False for key in required_false)
    )


def _unexpected_false_keys(result: Mapping[str, Any] | None) -> set[str]:
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
    """Run only the two new v3 negative source triples."""

    del repo, parent, source_texts
    target_case = _synthetic_texts(WRONG_OUT_TARGET, GOOD_SCALE_COLUMN_6)
    scale_case = _synthetic_texts(GOOD_OUT_TARGET, WRONG_SCALE_COLUMN_6)
    target_result = _invoke_case(auditor, target_case)
    scale_result = _invoke_case(auditor, scale_case)

    target_checks = _checks(target_result)
    scale_checks = _checks(scale_result)
    false_by_case = {
        "first_target_out_i_plus_1": _unexpected_false_keys(target_result),
        "first_raw_scale_column_6": _unexpected_false_keys(scale_result),
    }
    shared_unexpected_false_checks = sorted(
        false_by_case["first_target_out_i_plus_1"]
        & false_by_case["first_raw_scale_column_6"]
    )
    checks = {
        "first_target_case_rejected": _rejected(
            target_result,
            EXPECTED_TARGET_FALSE,
        ),
        "scale_column_6_case_rejected": _rejected(
            scale_result,
            "first_raw_SCALE_DIAGONAL_module_literal",
            "scale_mapping_by_literal_input_row",
        ),
        "target_case_scale_and_mapping_otherwise_pass": (
            target_checks.get("first_raw_SCALE_DIAGONAL_module_literal") is True
            and target_checks.get("scale_mapping_by_literal_input_row") is True
        ),
        "scale_case_first_accumulation_target_passes": (
            scale_checks.get(EXPECTED_TARGET_FALSE) is True
        ),
        "case_specific_failure_keys_differ": (
            target_checks.get(EXPECTED_TARGET_FALSE) is False
            and scale_checks.get("first_raw_SCALE_DIAGONAL_module_literal") is False
            and EXPECTED_TARGET_FALSE not in false_by_case["first_raw_scale_column_6"]
            and "first_raw_SCALE_DIAGONAL_module_literal"
            not in false_by_case["first_target_out_i_plus_1"]
        ),
        "no_shared_unexpected_false_checks": not shared_unexpected_false_checks,
    }
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
        "v2_fixture_called": False,
        "v1_or_v2_synthetic_cases_reentered": False,
        "implementation_imported": False,
        "old_modules_called": False,
        "provider_executed": False,
        "numerical_residual_evaluated": False,
        "synthetic_records": {
            "first_target_out_i_plus_1": "v3_new_first_target_out_i_plus_1",
            "first_raw_scale_column_6": "v3_new_first_raw_scale_column_6",
        },
        "scope": {
            "only_new_negative_cases": True,
            "actual_source_proof_deferred": True,
            "unknown_common_structure_stops_gate": True,
        },
    }


def run(auditor: Any, repo: str | Path | None = None, parent: Mapping[str, Any] | str | Path | None = None, **kwargs: Any) -> dict[str, Any]:
    return run_fixture(auditor, repo=repo, parent=parent, **kwargs)


__all__ = ["SCHEMA", "run_fixture", "run"]
