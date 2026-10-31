"""Independent synthetic checks for the static Q4 provider provenance gate.

The fixture owns its source bytes and computes its descriptor from a parsed AST
without importing or executing that source.  ``run_fixture`` is intentionally
call-time only: importing this module performs no verification.
"""

from __future__ import annotations

import ast
import hashlib
from typing import Any, Callable


SCHEMA = "original_q4_provider_provenance_gate_fixture_v1"
QUALNAME = "CurrentOriginalN1CoreFixedPointFunctions.component_enclosure"
CLASS_NAME = "CurrentOriginalN1CoreFixedPointFunctions"
METHOD_NAME = "component_enclosure"
FILENAME = "synthetic_original_q4_provider.py"
FAMILY = {"name": "synthetic-q4-family", "revision": 1}
IMPLEMENTATION_FILENAME = "synthetic_original_q4_provider.py"
WRAPPER_FILENAME = "synthetic_switch_signed_integrals.py"

SOURCE = (
    b"class CurrentOriginalN1CoreFixedPointFunctions:\n"
    b"    def component_enclosure(self, quantity, rho, Z):\n"
    b"        return (quantity, rho, Z)\n"
)

IMPLEMENTATION = (
    b"import common\n"
    b"class CurrentOriginalN1CoreFixedPointFunctions:\n"
    b"    @common.source_precision\n"
    b"    def component_enclosure(self, quantity, rho, Z):\n"
    b"        return (quantity, rho, Z)\n"
)

WRAPPER = (
    b"from functools import wraps\n"
    b"def source_precision(fn):\n"
    b"    @wraps(fn)\n"
    b"    def call(*args, **kwargs):\n"
    b"        return fn(*args, **kwargs)\n"
    b"    return call\n"
)


def _method_node(source: bytes) -> ast.FunctionDef | ast.AsyncFunctionDef:
    """Select the synthetic method for descriptor construction only."""

    tree = ast.parse(source.decode("utf-8"))
    classes = [
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == CLASS_NAME
    ]
    methods = [
        node
        for class_node in classes
        for node in class_node.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == METHOD_NAME
    ]
    if len(classes) != 1 or len(methods) != 1:
        raise ValueError("synthetic descriptor source is not unique")
    return methods[0]


def _ast_sha256(source: bytes, *, allow_duplicate: bool = False) -> str:
    """Hash only the method node using the gate's documented AST shape."""

    tree = ast.parse(source.decode("utf-8"))
    classes = [
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == CLASS_NAME
    ]
    methods = [
        node
        for class_node in classes
        for node in class_node.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == METHOD_NAME
    ]
    if len(classes) != 1 or (not methods) or (len(methods) != 1 and not allow_duplicate):
        raise ValueError("synthetic descriptor source is not unique")
    serialized = ast.dump(ast.Module(body=[methods[0]], type_ignores=[])).encode()
    return hashlib.sha256(serialized).hexdigest()


def _descriptor(source: bytes = SOURCE, *, allow_duplicate: bool = False) -> dict[str, Any]:
    """Construct a descriptor from source bytes, never from a provider object."""

    # Keep the method selection independent from the verifier implementation;
    # this helper repeats the AST shape and hash recipe locally on purpose.
    method = _method_node(source) if not allow_duplicate else None
    if method is not None:
        serialized = ast.dump(ast.Module(body=[method], type_ignores=[])).encode()
        ast_digest = hashlib.sha256(serialized).hexdigest()
    else:
        ast_digest = _ast_sha256(source, allow_duplicate=True)
    return {
        "module": FILENAME,
        "module_sha256": hashlib.sha256(source).hexdigest(),
        "function": QUALNAME,
        "AST_sha256": ast_digest,
        "family": dict(FAMILY),
    }


def _wrapped_descriptor(
    implementation: bytes = IMPLEMENTATION,
    wrapper: bytes = WRAPPER,
) -> dict[str, Any]:
    """Build wrapped evidence from the implementation method AST and wrapper bytes."""

    method = _method_node(implementation)
    serialized = ast.dump(ast.Module(body=[method], type_ignores=[])).encode()
    return {
        "module": WRAPPER_FILENAME,
        "module_sha256": hashlib.sha256(wrapper).hexdigest(),
        "function": QUALNAME,
        "AST_sha256": hashlib.sha256(serialized).hexdigest(),
    }


def _call(
    verifier: Callable[..., Any],
    descriptor: dict[str, Any],
    source: bytes,
    filename: str = FILENAME,
) -> dict[str, Any]:
    return verifier(
        descriptor,
        source,
        filename,
        expected_qualname=QUALNAME,
        expected_family=FAMILY,
    )


def _call_wrapped(
    verifier: Callable[..., Any],
    descriptor: dict[str, Any],
    implementation: bytes = IMPLEMENTATION,
    wrapper: bytes = WRAPPER,
    implementation_filename: str = IMPLEMENTATION_FILENAME,
    wrapper_filename: str = WRAPPER_FILENAME,
) -> dict[str, Any]:
    return verifier(
        descriptor,
        implementation,
        implementation_filename,
        wrapper,
        wrapper_filename,
        expected_qualname=QUALNAME,
    )


def _rejected(report: Any) -> bool:
    return isinstance(report, dict) and report.get("status") == "rejected"


def run_fixture(verify: Callable[..., Any], verify_wrapped: Callable[..., Any]) -> dict[str, Any]:
    """Run positive and negative checks against an injected ``verify`` callable."""

    checks: dict[str, bool] = {}
    if not callable(verify) or not callable(verify_wrapped):
        return {"schema": SCHEMA, "passed": False, "checks": {"verify_callable": False}}

    base = _descriptor()
    valid = _call(verify, base, SOURCE)
    checks["consistent_verified"] = (
        isinstance(valid, dict)
        and valid.get("status") == "verified"
        and valid.get("passed") is True
        and all(
            isinstance(item, dict) and item.get("match") is True
            for item in valid.get("checks", {}).values()
        )
    )

    module_mismatch = dict(base)
    module_mismatch["module"] = "another_synthetic_provider.py"
    checks["module_mismatch_rejected"] = _rejected(
        _call(verify, module_mismatch, SOURCE)
    )

    sha_mismatch = dict(base)
    sha_mismatch["module_sha256"] = "f" * 64
    checks["module_sha_mismatch_rejected"] = _rejected(
        _call(verify, sha_mismatch, SOURCE)
    )

    ast_mismatch = dict(base)
    ast_mismatch["AST_sha256"] = "0" * 64
    checks["AST_mismatch_rejected"] = _rejected(
        _call(verify, ast_mismatch, SOURCE)
    )

    qualname_mismatch = dict(base)
    qualname_mismatch["function"] = "OtherProvider.component_enclosure"
    checks["qualname_mismatch_rejected"] = _rejected(
        _call(verify, qualname_mismatch, SOURCE)
    )

    duplicate_source = SOURCE.replace(
        b"        return (quantity, rho, Z)\n",
        b"        return (quantity, rho, Z)\n"
        b"    def component_enclosure(self, quantity, rho, Z):\n"
        b"        return (quantity, rho, Z)\n",
    )
    duplicate_descriptor = _descriptor(duplicate_source, allow_duplicate=True)
    checks["duplicate_method_rejected"] = _rejected(
        _call(verify, duplicate_descriptor, duplicate_source)
    )

    checks["path_basename_rejected"] = _rejected(
        _call(verify, base, SOURCE, "nested/" + FILENAME)
    )

    modified_source = SOURCE.replace(
        b"return (quantity, rho, Z)",
        b"return (quantity, rho, Z, 'modified')",
    )
    checks["modified_code_rejected"] = _rejected(
        _call(verify, base, modified_source)
    )

    wrapped_descriptor = _wrapped_descriptor()
    wrapped = _call_wrapped(verify_wrapped, wrapped_descriptor)
    checks["wrapped_positive_matched"] = (
        isinstance(wrapped, dict)
        and wrapped.get("status") == "wrapper_evidence_matched"
        and wrapped.get("passed") is True
        and all(
            isinstance(item, dict) and item.get("match") is True
            for item in wrapped.get("checks", {}).values()
        )
    )

    bad_decorator = IMPLEMENTATION.replace(
        b"@common.source_precision", b"@common.other_precision"
    )
    checks["decorator_error_rejected"] = _rejected(
        _call_wrapped(verify_wrapped, wrapped_descriptor, implementation=bad_decorator)
    )

    missing_wraps = WRAPPER.replace(b"@wraps(fn)", b"@missing_wraps(fn)")
    missing_wraps_descriptor = _wrapped_descriptor(wrapper=missing_wraps)
    checks["wrapper_missing_wraps_rejected"] = _rejected(
        _call_wrapped(verify_wrapped, missing_wraps_descriptor, wrapper=missing_wraps)
    )

    return {"schema": SCHEMA, "passed": all(checks.values()), "checks": checks}


__all__ = ["run_fixture"]
