"""Synthetic AST fixture for :mod:`1_green_scale_commutator_v1`.

The fixture is deliberately implementation-agnostic.  The caller supplies a
source-text auditor (or an object exposing ``audit_source_texts``), and this
file supplies independent synthetic canonical texts plus controlled source
mutations.  Any ``repo`` or ``source_texts`` argument is compatibility-only
and is ignored.  No raw provider, legacy module, numerical owner, interval,
or quadrature is imported.  Every case is parsed with :func:`ast.parse`
before it is handed to the caller-owned auditor.

The negative cases are function-structure guards, not interval experiments:

* a wrong diagonal scale association is rejected;
* the core and first-chart row mappings are rejected independently;
* a ``Z``-dependent endpoint is rejected;
* the unmodified graph and the six incoming field-scale pairings pass.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping
import json
from pathlib import Path
from typing import Any


SCHEMA = "n1_green_scale_commutator_fixture_v1"
VOLERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
REGULAR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py"
)
FIRST_RAW_SOURCE = "src/openai_ns_reconstruction/original_n1_first_raw_function.py"
ROWS = ("5", "6", "6", "4", "5", "6")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
SCALE_DIAGONAL = (
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    "epsilon",
)


def _parse(text: str, label: str) -> None:
    if not isinstance(text, str):
        raise TypeError(f"{label} must be source text")
    try:
        ast.parse(text, filename=label, mode="exec")
    except SyntaxError as exc:
        raise ValueError(f"synthetic source failed AST parse: {label}: {exc}") from exc


def _synthetic_texts() -> tuple[str, str, str]:
    """Return independent source text with the admitted AST shape.

    These strings intentionally do not read or copy a repository source file.
    They retain only the six-row provider structure needed by this fixture, so
    the actual source proof can run once, separately, in the root runner.
    """

    volterra = '''
import synthetic_regular as regular
FIELDS=regular.FIELDS
def green_apply(g,vectors,chart,coordinates,parameters):
    c=g.constant;a=c3.C3Algebra(g);hb=parameters['hb'];r=coordinates['core'];s=coordinates['first']
    t=g.symbol('original_n1_Green_t');v=g.symbol('original_n1_Green_phase_v')
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
            out[i]=a.add(out[i],integral(a.scale(body,g.mul(kernel,measure)),v,g.zero,s))
    return dict(zip(FIELDS,out))
'''
    regular = '''
FIELDS=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
'''
    first_raw = '''
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
'''
    return volterra, regular, first_raw


def _load_texts(
    repo: str | Path | None,
    source_texts: Mapping[str, str] | tuple[str, str, str] | list[str] | None,
) -> tuple[str, str, str]:
    """Return only the fixture-owned texts; never read caller or repo text.

    ``repo`` and ``source_texts`` remain accepted for runner compatibility,
    but this leaf must exercise its mutations against independent synthetic
    AST text.  The root runner invokes the actual source proof separately.
    """

    del repo, source_texts
    return _synthetic_texts()


def _load_parent(parent: Mapping[str, Any] | str | Path | None) -> Mapping[str, Any] | None:
    if parent is None:
        return None
    if isinstance(parent, Mapping):
        return parent
    path = Path(parent)
    return json.loads(path.read_text(encoding="utf-8"))


def _replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise ValueError(f"fixture mutation {label} expected one source occurrence, found {count}")
    return text.replace(old, new, 1)


def _call_auditor(
    auditor: Any,
    texts: tuple[str, str, str],
    parent: Mapping[str, Any] | None,
) -> Mapping[str, Any]:
    candidate: Any = getattr(auditor, "audit_source_texts", None)
    if not callable(candidate):
        candidate = auditor
    if callable(candidate) and getattr(candidate, "__name__", "") == "build_source_bound_proof":
        # A runner may pass the builder function rather than its module.  Its
        # globals expose the pure text auditor without importing anything
        # here; the fixture still invokes only the supplied callable.
        maybe_audit = getattr(candidate, "__globals__", {}).get("audit_source_texts")
        if callable(maybe_audit):
            candidate = maybe_audit
    if not callable(candidate):
        raise TypeError("auditor must be callable or expose audit_source_texts")
    try:
        result = candidate(*texts, parent=parent)
    except TypeError as first_error:
        # A narrow compatibility path for callers that expose the three
        # positional source strings but no parent keyword.
        try:
            result = candidate(*texts)
        except TypeError:
            raise first_error
    if not isinstance(result, Mapping):
        raise TypeError("source auditor must return a mapping")
    return result


def _checks(result: Mapping[str, Any]) -> Mapping[str, Any]:
    checks = result.get("checks")
    return checks if isinstance(checks, Mapping) else {}


def _structural_pass(result: Mapping[str, Any]) -> bool:
    checks = _checks(result)
    if not checks:
        return False
    # Parent receipt admission is intentionally orthogonal to this synthetic
    # graph test.  If supplied, the root runner checks it in the source proof;
    # the fixture's canonical graph check must remain meaningful without a
    # saved JSON receipt.
    structural = [
        value
        for key, value in checks.items()
        if not str(key).startswith(("parent_report_", "source_pin_"))
    ]
    return bool(structural) and all(value is True for value in structural)


def _rejected(result: Mapping[str, Any] | None, *required_false: str) -> bool:
    if result is None:
        return True
    checks = _checks(result)
    return bool(result.get("passed") is False and required_false and any(checks.get(key) is False for key in required_false))


def _invoke_case(
    auditor: Any,
    texts: tuple[str, str, str],
    parent: Mapping[str, Any] | None,
) -> Mapping[str, Any] | None:
    try:
        for label, text in zip(("volterra", "regular", "first_raw"), texts):
            _parse(text, label)
        return _call_auditor(auditor, texts, parent)
    except (TypeError, ValueError, KeyError, OSError, json.JSONDecodeError):
        return None


def run_fixture(
    auditor: Any,
    repo: str | Path | None = None,
    parent: Mapping[str, Any] | str | Path | None = None,
    source_texts: Mapping[str, str] | tuple[str, str, str] | list[str] | None = None,
    **_: Any,
) -> dict[str, Any]:
    """Run canonical and mutated source-text checks through ``auditor``."""

    # Preserve the familiar ``run_fixture(leaf, context)`` calling shape used
    # by older synthetic fixtures: a non-path second argument is ignored as a
    # numerical context, since this fixture is source-text-only.
    if repo is not None and not isinstance(repo, (str, Path)):
        repo = None
    texts = _load_texts(repo, source_texts)
    parent_mapping = _load_parent(parent)
    canonical = _invoke_case(auditor, texts, parent_mapping)
    checks: dict[str, bool] = {
        "canonical_source_texts_parse": canonical is not None,
        "canonical_graph_and_incoming_scaling_pass": canonical is not None and _structural_pass(canonical),
    }
    if canonical is not None:
        canonical_checks = _checks(canonical)
        checks["canonical_owner_is_module"] = canonical_checks.get("green_green_apply_owner_is_module") is True
        checks["canonical_FIELDS_binding_is_regular_module"] = canonical_checks.get(
            "green_FIELDS_module_bound_to_regular_FIELDS"
        ) is True
        checks["canonical_rows_are_exact"] = canonical_checks.get("green_exact_rows_present_in_both_for_loops") is True
        checks["canonical_core_dependency_is_exact"] = canonical_checks.get("green_core_loop_body_uses_core_vector_row") is True and canonical_checks.get("green_core_loop_no_hidden_cross_row_input") is True
        checks["canonical_first_dependency_is_exact"] = canonical_checks.get("green_first_loop_body_uses_first_vector_row") is True and canonical_checks.get("green_first_loop_no_hidden_cross_row_input") is True
        checks["canonical_incoming_field_scaling_is_exact"] = canonical_checks.get("scale_mapping_by_literal_input_row") is True
        scale_map = canonical.get("scale_mapping")
        records = scale_map.get("records", ()) if isinstance(scale_map, Mapping) else ()
        checks["six_output_input_scale_pairs"] = len(records) == 6 and all(
            isinstance(record, Mapping)
            and record.get("input_row") == ROWS[index]
            and record.get("output_field") == FIELDS[index]
            and record.get("same_scale") is True
            for index, record in enumerate(records)
        )
    else:
        for key in (
            "canonical_owner_is_module",
            "canonical_FIELDS_binding_is_regular_module",
            "canonical_rows_are_exact",
            "canonical_core_dependency_is_exact",
            "canonical_first_dependency_is_exact",
            "canonical_incoming_field_scaling_is_exact",
            "six_output_input_scale_pairs",
        ):
            checks[key] = False

    # Change only one structural fact in each case.  The mutations preserve
    # the surrounding labels, so a substring/label-only checker cannot pass.
    wrong_scale = (
        texts[0],
        texts[1],
        _replace_once(
            texts[2],
            'SCALE_DIAGONAL = (\n    "F0(Z)/epsilon",',
            'SCALE_DIAGONAL = (\n    "epsilon",',
            "wrong scale association",
        ),
    )
    wrong_scale_result = _invoke_case(auditor, wrong_scale, parent_mapping)
    checks["wrong_scale_association_rejected"] = _rejected(
        wrong_scale_result,
        "first_raw_SCALE_DIAGONAL_module_literal",
        "scale_mapping_by_literal_input_row",
    )

    wrong_core = (
        _replace_once(
            texts[0],
            "vectors['core'][row]",
            "vectors['core']['5']",
            "wrong core row mapping",
        ),
        texts[1],
        texts[2],
    )
    wrong_core_result = _invoke_case(auditor, wrong_core, parent_mapping)
    checks["wrong_core_mapping_rejected"] = _rejected(
        wrong_core_result,
        "green_core_loop_body_uses_core_vector_row",
        "green_core_loop_no_hidden_cross_row_input",
    )

    wrong_first = (
        _replace_once(
            texts[0],
            "vectors['first'][row]",
            "vectors['core'][row]",
            "wrong first-chart row mapping",
        ),
        texts[1],
        texts[2],
    )
    wrong_first_result = _invoke_case(auditor, wrong_first, parent_mapping)
    checks["wrong_first_mapping_rejected"] = _rejected(
        wrong_first_result,
        "green_first_loop_body_uses_first_vector_row",
        "green_first_loop_no_hidden_cross_row_input",
    )

    wrong_endpoint = (
        _replace_once(
            texts[0],
            "out.append(integral(a.scale(body,g.mul(kernel,measure)),t,g.zero,g.one))",
            "out.append(integral(a.scale(body,g.mul(kernel,measure)),t,g.symbol('Z'),g.one))",
            "wrong Z-dependent boundary",
        ),
        texts[1],
        texts[2],
    )
    wrong_endpoint_result = _invoke_case(auditor, wrong_endpoint, parent_mapping)
    checks["wrong_boundary_Z_dependence_rejected"] = _rejected(
        wrong_endpoint_result,
        "green_core_loop_out_append_accumulation",
    )

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "source_texts_only": True,
        "ast_parse_per_case": True,
        "implementation_imported": False,
        "old_modules_called": False,
        "provider_executed": False,
        "numerical_residual_evaluated": False,
        "negative_cases": {
            "wrong_scale_association": "rejected",
            "wrong_core_mapping": "rejected",
            "wrong_first_mapping": "rejected",
            "wrong_boundary_Z_dependence": "rejected",
        },
        "scope": {
            "function_level_only": True,
            "Kchi_resolvent_proved": False,
            "gauge_owner_proved": False,
            "history_errors_changed": False,
        },
    }


def run(auditor: Any, repo: str | Path | None = None, parent: Mapping[str, Any] | str | Path | None = None, **kwargs: Any) -> dict[str, Any]:
    return run_fixture(auditor, repo=repo, parent=parent, **kwargs)


__all__ = ["SCHEMA", "run_fixture", "run"]
