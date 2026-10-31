"""Static provenance admission for one original Q4 provider callable.

This module checks a previously authenticated descriptor against caller-owned
Python source bytes.  It deliberately has no provider import, source loading,
callable identity check, source execution, graph access, or numerical/error
claim.  The only source operation is UTF-8 decoding followed by ``ast.parse``.

The target is intentionally narrow.  A descriptor for another method is not
accepted merely because it happens to have a valid hash: the admitted
qualified name is the one original component method named below.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping
import hashlib
from typing import Any


SCHEMA = "original_q4_provider_provenance_gate_v1"
EXPECTED_QUALNAME = "CurrentOriginalN1CoreFixedPointFunctions.component_enclosure"
EXPECTED_CLASS = "CurrentOriginalN1CoreFixedPointFunctions"
EXPECTED_METHOD = "component_enclosure"

_REQUIRED_DESCRIPTOR_FIELDS = (
    "function",
    "module",
    "module_sha256",
    "AST_sha256",
)


def _check(actual: Any, claimed: Any, match: bool) -> dict[str, Any]:
    """Create the uniform actual/claimed/match shape used by every check."""

    return {"actual": actual, "claimed": claimed, "match": bool(match)}


def _basename_candidate(value: Any) -> str | None:
    """Return the lexical final component without accepting it as valid.

    Both slash conventions are handled explicitly so that a Windows path is
    rejected even when the gate is called from a POSIX process, and vice versa.
    """

    if not isinstance(value, str):
        return None
    return value.replace("\\", "/").rsplit("/", 1)[-1]


def _is_strict_basename(value: Any) -> bool:
    """Accept only a non-empty literal basename, never a path alias."""

    if not isinstance(value, str) or not value:
        return False
    if value in {".", ".."}:
        return False
    if "\x00" in value or "/" in value or "\\" in value or ":" in value:
        return False
    return _basename_candidate(value) == value


def _descriptor_value(descriptor: Mapping[str, Any], key: str) -> tuple[Any, bool]:
    """Read one claimed field while preserving the missing-field distinction."""

    if key not in descriptor:
        return None, False
    return descriptor[key], True


def _target_methods(tree: ast.Module) -> tuple[int, list[ast.FunctionDef | ast.AsyncFunctionDef]]:
    """Find direct methods of the one permitted top-level class.

    Nested classes, functions, aliases, and arbitrary qualified names are not
    promoted to this target.  Requiring one class definition and one direct
    method definition makes duplicate or unknown target ownership rejectable.
    """

    classes = [
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == EXPECTED_CLASS
    ]
    methods: list[ast.FunctionDef | ast.AsyncFunctionDef] = []
    for class_node in classes:
        methods.extend(
            node
            for node in class_node.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == EXPECTED_METHOD
        )
    return len(classes), methods


def _source_precision_decorator(node: ast.AST) -> bool:
    """Recognize exactly the implementation decorator used by the wrapper."""

    return (
        isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "common"
        and node.attr == "source_precision"
        and isinstance(node.ctx, ast.Load)
    )


def _implementation_target(
    tree: ast.Module,
) -> tuple[int, list[ast.FunctionDef | ast.AsyncFunctionDef], bool]:
    """Find the target and report whether its decorator AST is exact."""

    class_count, methods = _target_methods(tree)
    decorator_ok = (
        class_count == 1
        and len(methods) == 1
        and len(methods[0].decorator_list) == 1
        and _source_precision_decorator(methods[0].decorator_list[0])
    )
    return class_count, methods, decorator_ok


def _forwarding_call(node: ast.Return) -> bool:
    """Check the static ``return fn(*args, **kwargs)`` forwarding shape."""

    value = node.value
    if not isinstance(value, ast.Call):
        return False
    if not isinstance(value.func, ast.Name) or value.func.id != "fn":
        return False
    if len(value.args) != 1 or not isinstance(value.args[0], ast.Starred):
        return False
    starred = value.args[0].value
    if not isinstance(starred, ast.Name) or starred.id != "args":
        return False
    if len(value.keywords) != 1:
        return False
    keyword = value.keywords[0]
    return (
        keyword.arg is None
        and isinstance(keyword.value, ast.Name)
        and keyword.value.id == "kwargs"
    )


def _wrapper_contract(tree: ast.Module) -> tuple[bool, dict[str, Any]]:
    """Audit the source-level ``source_precision``/``call`` wrapper recipe."""

    source_precision = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "source_precision"
    ]
    details: dict[str, Any] = {
        "source_precision_count": len(source_precision),
        "call_count": 0,
        "wraps_fn": False,
        "returns_call": False,
        "forwards_fn": False,
    }
    if len(source_precision) != 1:
        return False, details

    producer = source_precision[0]
    positional = [*producer.args.posonlyargs, *producer.args.args]
    signature_ok = (
        len(positional) == 1
        and positional[0].arg == "fn"
        and not producer.args.vararg
        and not producer.args.kwarg
        and not producer.args.kwonlyargs
    )
    calls = [
        node
        for node in producer.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "call"
    ]
    details["call_count"] = len(calls)
    if len(calls) != 1:
        return False, details
    call = calls[0]
    wraps_decorators = [
        decorator
        for decorator in call.decorator_list
        if isinstance(decorator, ast.Call)
        and isinstance(decorator.func, ast.Name)
        and decorator.func.id == "wraps"
        and len(decorator.args) == 1
        and isinstance(decorator.args[0], ast.Name)
        and decorator.args[0].id == "fn"
        and not decorator.keywords
    ]
    details["wraps_fn"] = len(call.decorator_list) == 1 and len(wraps_decorators) == 1

    producer_returns = [
        statement
        for statement in producer.body
        if isinstance(statement, ast.Return)
        and isinstance(statement.value, ast.Name)
        and statement.value.id == "call"
    ]
    details["returns_call"] = len(producer_returns) == 1
    forwarding_returns = [
        node
        for node in ast.walk(call)
        if isinstance(node, ast.Return)
    ]
    details["forwards_fn"] = (
        len(forwarding_returns) == 1 and _forwarding_call(forwarding_returns[0])
    )
    return (
        signature_ok
        and details["wraps_fn"]
        and details["returns_call"]
        and details["forwards_fn"],
        details,
    )


def _rejected_report(
    checks: dict[str, dict[str, Any]],
    *,
    reasons: list[str] | None = None,
) -> dict[str, Any]:
    """Return the common fail-closed report without adding semantic claims."""

    return {
        "schema": SCHEMA,
        "status": "rejected",
        "passed": False,
        "checks": checks,
        "errors": list(reasons or ()),
        "source_eval": False,
        "cross_owner_identity": False,
        "certified_errors": False,
        "scope": {
            "source_eval": False,
            "source_values_read": False,
            "cross_owner_identity": False,
            "cross_owner_identity_proved": False,
            "certified_errors": False,
            "static_descriptor_only": True,
        },
    }


def verify(
    descriptor: Mapping[str, Any],
    source_bytes: bytes,
    source_filename: str,
    *,
    expected_qualname: str,
    expected_family: Any = None,
) -> dict[str, Any]:
    """Verify a provider descriptor against caller-supplied raw source bytes.

    The return value is always a report.  Malformed inputs, invalid UTF-8,
    syntax errors, missing fields, unknown targets, and duplicate targets all
    become ``status == "rejected"``; no malformed input is repaired or
    inferred.  ``expected_family`` is an optional metadata comparison only;
    this gate does not derive a family from source or assert owner identity.
    """

    checks: dict[str, dict[str, Any]] = {}
    reasons: list[str] = []

    # Keep the malformed-input path reportable rather than allowing a caller
    # supplied object to turn a failed admission into an exception.
    if not isinstance(descriptor, Mapping):
        return _rejected_report(
            {
                "descriptor": _check(None, "mapping", False),
            },
            reasons=["descriptor_not_mapping"],
        )

    source_is_bytes = isinstance(source_bytes, bytes)
    filename_is_basename = _is_strict_basename(source_filename)
    checks["basename"] = _check(
        _basename_candidate(source_filename),
        source_filename,
        filename_is_basename,
    )
    checks["source_filename"] = _check(
        _basename_candidate(source_filename),
        source_filename,
        filename_is_basename,
    )
    if not filename_is_basename:
        reasons.append("source_filename_not_strict_basename")

    expected_name_ok = expected_qualname == EXPECTED_QUALNAME
    checks["expected_qualname"] = _check(
        EXPECTED_QUALNAME,
        expected_qualname,
        expected_name_ok,
    )
    if not expected_name_ok:
        reasons.append("unexpected_expected_qualname")

    claimed: dict[str, Any] = {}
    present: dict[str, bool] = {}
    try:
        for key in _REQUIRED_DESCRIPTOR_FIELDS:
            claimed[key], present[key] = _descriptor_value(descriptor, key)
            if not present[key]:
                reasons.append(f"missing_{key}")
    except Exception:
        return _rejected_report(
            {"descriptor": _check(None, "readable mapping", False)},
            reasons=["descriptor_not_readable"],
        )

    # The module digest is always over the exact bytes supplied by the caller;
    # a claimed path or claimed digest cannot substitute for this computation.
    actual_module_sha256: str | None = None
    if source_is_bytes:
        actual_module_sha256 = hashlib.sha256(source_bytes).hexdigest()
    checks["module_sha256"] = _check(
        actual_module_sha256,
        claimed["module_sha256"],
        source_is_bytes
        and type(claimed["module_sha256"]) is str
        and present["module_sha256"]
        and claimed["module_sha256"] == actual_module_sha256,
    )
    if not checks["module_sha256"]["match"]:
        reasons.append("module_sha256_mismatch")

    checks["module"] = _check(
        source_filename if isinstance(source_filename, str) else None,
        claimed["module"],
        filename_is_basename
        and type(claimed["module"]) is str
        and present["module"]
        and claimed["module"] == source_filename,
    )
    if not checks["module"]["match"]:
        reasons.append("module_mismatch")

    actual_function: str | None = None
    actual_ast_sha256: str | None = None
    class_count = 0
    method_count = 0
    syntax_ok = False

    if not source_is_bytes:
        reasons.append("source_bytes_not_bytes")
    else:
        try:
            source_text = source_bytes.decode("utf-8")
            tree = ast.parse(source_text)
            syntax_ok = True
            class_count, methods = _target_methods(tree)
            method_count = len(methods)
            if class_count == 1 and method_count == 1:
                actual_function = EXPECTED_QUALNAME
                method_node = methods[0]
                method_module = ast.Module(body=[method_node], type_ignores=[])
                serialized = ast.dump(method_module).encode()
                actual_ast_sha256 = hashlib.sha256(serialized).hexdigest()
            elif class_count == 0:
                reasons.append("target_class_not_found")
            elif class_count != 1:
                reasons.append("duplicate_target_class")
            elif method_count == 0:
                reasons.append("target_method_not_found")
            else:
                reasons.append("duplicate_target_method")
        except Exception:
            reasons.append("source_not_valid_utf8_python")

    checks["target"] = _check(
        {
            "class": EXPECTED_CLASS if class_count == 1 else None,
            "method": EXPECTED_METHOD if method_count == 1 else None,
            "class_count": class_count,
            "method_count": method_count,
            "syntax": syntax_ok,
        },
        EXPECTED_QUALNAME,
        syntax_ok
        and class_count == 1
        and method_count == 1
        and expected_name_ok,
    )
    checks["function"] = _check(
        actual_function,
        claimed["function"],
        expected_name_ok
        and actual_function == EXPECTED_QUALNAME
        and type(claimed["function"]) is str
        and present["function"]
        and claimed["function"] == EXPECTED_QUALNAME,
    )
    if not checks["function"]["match"]:
        reasons.append("function_mismatch")

    checks["AST_sha256"] = _check(
        actual_ast_sha256,
        claimed["AST_sha256"],
        actual_ast_sha256 is not None
        and type(claimed["AST_sha256"]) is str
        and present["AST_sha256"]
        and claimed["AST_sha256"] == actual_ast_sha256,
    )
    if not checks["AST_sha256"]["match"]:
        reasons.append("AST_sha256_mismatch")

    if expected_family is not None:
        # Family is metadata only.  ``family`` is canonical; ``source_family``
        # is accepted as the established receipt spelling.  If both spellings
        # occur, they must agree rather than allowing an alias to override it.
        try:
            family_claims = [
                descriptor[key]
                for key in ("family", "source_family")
                if key in descriptor
            ]
        except Exception:
            family_claims = []
        family_present = bool(family_claims)
        claimed_family = family_claims[0] if family_claims else None
        family_aliases_agree = len(family_claims) <= 1 or all(
            value == claimed_family for value in family_claims[1:]
        )
        checks["family"] = _check(
            claimed_family if family_present else None,
            expected_family,
            family_present
            and family_aliases_agree
            and claimed_family == expected_family,
        )
        if not checks["family"]["match"]:
            reasons.append("family_mismatch")

    matched = all(item["match"] for item in checks.values())
    if not matched:
        # Preserve first occurrence order while avoiding repeated diagnostics
        # from independent checks that fail for the same malformed input.
        reasons = list(dict.fromkeys(reasons))
        return _rejected_report(checks, reasons=reasons)

    return {
        "schema": SCHEMA,
        "status": "verified",
        "passed": True,
        "checks": checks,
        "errors": [],
        "source_eval": False,
        "cross_owner_identity": False,
        "certified_errors": False,
        "scope": {
            "source_eval": False,
            "source_values_read": False,
            "cross_owner_identity": False,
            "cross_owner_identity_proved": False,
            "certified_errors": False,
            "static_descriptor_only": True,
        },
    }


def _wrapper_report(
    status: str,
    checks: dict[str, dict[str, Any]],
    *,
    reasons: list[str] | None = None,
) -> dict[str, Any]:
    """Return a wrapper-evidence-only report."""

    matched = status == "wrapper_evidence_matched"
    return {
        "schema": SCHEMA,
        "status": status,
        "passed": matched,
        "checks": checks,
        "errors": list(reasons or ()),
        "source_eval": False,
        "live_callable_identity": False,
        "source_identity": False,
        "cross_owner_identity": False,
        "certified_errors": False,
        "scope": {
            "wrapper_evidence_only": True,
            "source_eval": False,
            "source_values_read": False,
            "live_callable_identity": False,
            "source_identity": False,
            "cross_owner_identity": False,
            "cross_owner_identity_proved": False,
            "certified_errors": False,
            "static_descriptor_only": True,
        },
    }


def verify_wrapped(
    descriptor: Mapping[str, Any],
    implementation_bytes: bytes,
    implementation_filename: str,
    wrapper_bytes: bytes,
    wrapper_filename: str,
    *,
    expected_qualname: str,
) -> dict[str, Any]:
    """Verify static evidence for a decorated provider callable.

    ``inspect.getfile`` for the live wrapped function points at the wrapper
    module, while ``inspect.getsource`` follows ``__wrapped__`` to the
    implementation method.  This entry point accepts both source texts and
    checks those two static facts without importing either module or claiming
    that a live callable or source owner is identical.
    """

    checks: dict[str, dict[str, Any]] = {}
    reasons: list[str] = []
    if not isinstance(descriptor, Mapping):
        return _wrapper_report(
            "rejected",
            {"descriptor": _check(None, "mapping", False)},
            reasons=["descriptor_not_mapping"],
        )

    implementation_basename = _is_strict_basename(implementation_filename)
    wrapper_basename = _is_strict_basename(wrapper_filename)
    checks["implementation_basename"] = _check(
        _basename_candidate(implementation_filename),
        implementation_filename,
        implementation_basename,
    )
    checks["wrapper_basename"] = _check(
        _basename_candidate(wrapper_filename),
        wrapper_filename,
        wrapper_basename,
    )
    if not implementation_basename:
        reasons.append("implementation_filename_not_strict_basename")
    if not wrapper_basename:
        reasons.append("wrapper_filename_not_strict_basename")

    expected_name_ok = expected_qualname == EXPECTED_QUALNAME
    checks["expected_qualname"] = _check(
        EXPECTED_QUALNAME,
        expected_qualname,
        expected_name_ok,
    )
    if not expected_name_ok:
        reasons.append("unexpected_expected_qualname")

    claimed: dict[str, Any] = {}
    present: dict[str, bool] = {}
    try:
        for key in _REQUIRED_DESCRIPTOR_FIELDS:
            claimed[key], present[key] = _descriptor_value(descriptor, key)
            if not present[key]:
                reasons.append(f"missing_{key}")
    except Exception:
        return _wrapper_report(
            "rejected",
            {"descriptor": _check(None, "readable mapping", False)},
            reasons=["descriptor_not_readable"],
        )

    wrapper_sha256: str | None = None
    if isinstance(wrapper_bytes, bytes):
        wrapper_sha256 = hashlib.sha256(wrapper_bytes).hexdigest()
    else:
        reasons.append("wrapper_bytes_not_bytes")
    checks["module"] = _check(
        wrapper_filename if isinstance(wrapper_filename, str) else None,
        claimed["module"],
        wrapper_basename
        and present["module"]
        and type(claimed["module"]) is str
        and claimed["module"] == wrapper_filename,
    )
    checks["module_sha256"] = _check(
        wrapper_sha256,
        claimed["module_sha256"],
        wrapper_sha256 is not None
        and present["module_sha256"]
        and type(claimed["module_sha256"]) is str
        and claimed["module_sha256"] == wrapper_sha256,
    )
    if not checks["module"]["match"]:
        reasons.append("wrapper_module_mismatch")
    if not checks["module_sha256"]["match"]:
        reasons.append("wrapper_module_sha256_mismatch")

    implementation_function: str | None = None
    implementation_ast_sha256: str | None = None
    implementation_class_count = 0
    implementation_method_count = 0
    implementation_decorator_ok = False
    implementation_syntax_ok = False
    if not isinstance(implementation_bytes, bytes):
        reasons.append("implementation_bytes_not_bytes")
    else:
        try:
            implementation_tree = ast.parse(implementation_bytes.decode("utf-8"))
            implementation_syntax_ok = True
            (
                implementation_class_count,
                implementation_methods,
                implementation_decorator_ok,
            ) = _implementation_target(implementation_tree)
            implementation_method_count = len(implementation_methods)
            if implementation_class_count == 1 and implementation_method_count == 1:
                implementation_function = EXPECTED_QUALNAME
                method_node = implementation_methods[0]
                method_module = ast.Module(body=[method_node], type_ignores=[])
                implementation_ast_sha256 = hashlib.sha256(
                    ast.dump(method_module).encode()
                ).hexdigest()
            if implementation_class_count == 0:
                reasons.append("target_class_not_found")
            elif implementation_class_count != 1:
                reasons.append("duplicate_target_class")
            elif implementation_method_count == 0:
                reasons.append("target_method_not_found")
            elif implementation_method_count != 1:
                reasons.append("duplicate_target_method")
            if not implementation_decorator_ok:
                reasons.append("implementation_decorator_mismatch")
        except Exception:
            reasons.append("implementation_not_valid_utf8_python")

    checks["implementation_target"] = _check(
        {
            "class": EXPECTED_CLASS if implementation_class_count == 1 else None,
            "method": EXPECTED_METHOD if implementation_method_count == 1 else None,
            "class_count": implementation_class_count,
            "method_count": implementation_method_count,
            "syntax": implementation_syntax_ok,
        },
        EXPECTED_QUALNAME,
        implementation_syntax_ok
        and implementation_class_count == 1
        and implementation_method_count == 1
        and expected_name_ok,
    )
    checks["function"] = _check(
        implementation_function,
        claimed["function"],
        expected_name_ok
        and implementation_function == EXPECTED_QUALNAME
        and present["function"]
        and type(claimed["function"]) is str
        and claimed["function"] == EXPECTED_QUALNAME,
    )
    checks["implementation_decorator"] = _check(
        "common.source_precision" if implementation_decorator_ok else None,
        "common.source_precision",
        implementation_decorator_ok,
    )
    checks["AST_sha256"] = _check(
        implementation_ast_sha256,
        claimed["AST_sha256"],
        implementation_ast_sha256 is not None
        and present["AST_sha256"]
        and type(claimed["AST_sha256"]) is str
        and claimed["AST_sha256"] == implementation_ast_sha256,
    )
    if not checks["function"]["match"]:
        reasons.append("function_mismatch")
    if not checks["AST_sha256"]["match"]:
        reasons.append("AST_sha256_mismatch")

    wrapper_syntax_ok = False
    wrapper_contract_ok = False
    wrapper_details: dict[str, Any] = {
        "source_precision_count": 0,
        "call_count": 0,
        "wraps_fn": False,
        "returns_call": False,
        "forwards_fn": False,
    }
    if isinstance(wrapper_bytes, bytes):
        try:
            wrapper_tree = ast.parse(wrapper_bytes.decode("utf-8"))
            wrapper_syntax_ok = True
            wrapper_contract_ok, wrapper_details = _wrapper_contract(wrapper_tree)
        except Exception:
            reasons.append("wrapper_not_valid_utf8_python")
    wrapper_details["syntax"] = wrapper_syntax_ok
    checks["wrapper_contract"] = _check(
        wrapper_details,
        {
            "source_precision": "@wraps(fn)",
            "returns": "call",
            "forwarding": "fn(*args, **kwargs)",
        },
        wrapper_syntax_ok and wrapper_contract_ok,
    )
    if not checks["wrapper_contract"]["match"]:
        reasons.append("wrapper_contract_mismatch")

    matched = all(item["match"] for item in checks.values())
    reasons = list(dict.fromkeys(reasons))
    if not matched:
        return _wrapper_report("rejected", checks, reasons=reasons)
    return _wrapper_report("wrapper_evidence_matched", checks)


__all__ = [
    "EXPECTED_CLASS",
    "EXPECTED_METHOD",
    "EXPECTED_QUALNAME",
    "SCHEMA",
    "verify",
    "verify_wrapped",
]
