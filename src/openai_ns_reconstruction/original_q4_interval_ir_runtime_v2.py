"""Interval runtime for the source-preserving Q4 v4 IR.

This module is an intentionally independent runtime boundary.  It does not
import the existing scalar evaluator, the interval kernel, a core provider,
or a graph/archive helper.  Every numeric operation is delegated to objects
injected by the caller.  In particular, a caller supplies the exact interval
class, arithmetic backend, component extractor, endpoint converter, and
range integrator.

The result of a successful call is a *conditional* interval enclosure.  The
runtime records what the injected callbacks returned, but it does not verify
source identity, endpoint conversion, binder provenance, or quadrature
correctness.  Those remain caller claims and ``certified`` is consequently
always ``False``.
"""

from __future__ import annotations

from collections.abc import Mapping
import inspect
from typing import Any


SCHEMA = "original_q4_interval_ir_runtime_v2"

_MAX_DEPTH = 128
_MAX_NODES = 100_000
_MAX_PARTITIONS = 100_000
_MISSING = object()

_ALLOWED_OPS = {
    "exact_rational",
    "sum",
    "product",
    "negative",
    "positive_quotient",
    "analytic_unary",
    "bound_variable",
    "definite_integral",
    "current_original_source_parameter",
    "current_original_unique_analytic_core_component",
}
_ALLOWED_ANALYTIC_UNARY = {"exp", "abs"}


class _Rejected(Exception):
    """Internal fail-closed signal."""

    def __init__(self, code: str, detail: Any = None):
        self.code = code
        self.detail = detail
        super().__init__(code)


def _member(value: Any, names: tuple[str, ...]) -> Any:
    """Read only explicitly named members from a mapping or object."""

    if isinstance(value, Mapping):
        for name in names:
            if name in value:
                return value[name]
        return _MISSING
    for name in names:
        try:
            return getattr(value, name)
        except Exception:
            continue
    return _MISSING


def _call_noarg(value: Any, error_code: str) -> Any:
    if not callable(value):
        return value
    try:
        return value()
    except _Rejected:
        raise
    except Exception as exc:
        raise _Rejected(error_code, type(exc).__name__) from exc


def _same_reference(left: Any, right: Any) -> bool:
    """Compare opaque source references without int/bool coercion."""

    if type(left) is not type(right):
        return False
    try:
        result = left == right
    except Exception:
        return False
    return type(result) is bool and result


def _copy_opaque(value: Any, depth: int = 0, active: set[int] | None = None) -> Any:
    """Copy JSON-like provenance while never stringifying arbitrary objects."""

    if depth > _MAX_DEPTH:
        raise _Rejected("metadata_depth_limit")
    if active is None:
        active = set()
    if isinstance(value, Mapping):
        identity = id(value)
        if identity in active:
            raise _Rejected("cyclic_metadata")
        active.add(identity)
        try:
            return {
                key: _copy_opaque(item, depth + 1, active)
                for key, item in value.items()
            }
        finally:
            active.remove(identity)
    if isinstance(value, list):
        identity = id(value)
        if identity in active:
            raise _Rejected("cyclic_metadata")
        active.add(identity)
        try:
            return [_copy_opaque(item, depth + 1, active) for item in value]
        finally:
            active.remove(identity)
    if isinstance(value, tuple):
        identity = id(value)
        if identity in active:
            raise _Rejected("cyclic_metadata")
        active.add(identity)
        try:
            return tuple(_copy_opaque(item, depth + 1, active) for item in value)
        finally:
            active.remove(identity)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise _Rejected("non_json_metadata")


def _require_interval(value: Any, interval_type: Any, code: str) -> Any:
    """Require the exact injected interval class and its enclosure tag."""

    if not isinstance(interval_type, type) or type(value) is not interval_type:
        raise _Rejected(code)
    if getattr(value, "kind", _MISSING) != "enclosure":
        raise _Rejected(code)
    return value


def _normalise_mode(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    return value.strip().lower()


def _valid_family(value: Any) -> bool:
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, Mapping):
        return bool(value)
    return False


def _valid_graph_hash(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(char in "0123456789abcdefABCDEF" for char in value)
    )


def _first(container: Any, names: tuple[str, ...], default: Any = _MISSING) -> Any:
    if not isinstance(container, Mapping):
        return default
    for name in names:
        if name in container:
            return container[name]
    return default


def _unwrap_expression(expression: Any) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    if not isinstance(expression, Mapping):
        raise _Rejected("expression_not_mapping")
    if "op" not in expression and "ir" in expression:
        if "status" in expression and expression.get("status") not in (None, "lowered"):
            raise _Rejected("expression_not_lowered")
        ir = expression.get("ir")
        if not isinstance(ir, Mapping):
            raise _Rejected("missing_ir")
        return ir, expression
    return expression, expression


def _oracle_report(oracle: Any) -> dict[str, Any]:
    if isinstance(oracle, (str, bytes, bytearray, list, tuple, set, frozenset)):
        raise _Rejected("invalid_oracle")
    report: dict[str, Any] = {}
    for name in ("reference_report", "mode_report", "admission_report"):
        member = _member(oracle, (name,))
        if member is _MISSING:
            continue
        value = _call_noarg(member, "oracle_report_error")
        if isinstance(value, Mapping):
            report.update(value)
    direct = _member(oracle, ("report",))
    if direct is not _MISSING:
        value = _call_noarg(direct, "oracle_report_error")
        if isinstance(value, Mapping):
            report.update(value)
    for name in (
        "mode",
        "reference_mode",
        "source_mode",
        "evaluation_mode",
        "source_family",
        "family",
        "graph_hash",
        "graph_sha256",
        "source_graph_hash",
        "source_graph_sha256",
    ):
        member = _member(oracle, (name,))
        if member is not _MISSING:
            report.setdefault(name, _call_noarg(member, "oracle_report_error"))
    return report


def _metadata_value(container: Any, names: tuple[str, ...]) -> Any:
    if not isinstance(container, Mapping):
        return _MISSING
    value = _first(container, names, _MISSING)
    if value is not _MISSING:
        return value
    report = container.get("report")
    if isinstance(report, Mapping):
        value = _first(report, names, _MISSING)
        if value is not _MISSING:
            return value
    return _MISSING


def _merge_values(left: Any, right: Any, code: str) -> Any:
    if left is not _MISSING and right is not _MISSING and left != right:
        raise _Rejected(code)
    return left if left is not _MISSING else right


def _admission_container(admission: Any) -> Mapping[str, Any]:
    if not isinstance(admission, Mapping):
        raise _Rejected("admission_required")
    nested = _first(
        admission,
        ("admission_certificate", "capture_safe_admission_certificate", "certificate"),
        _MISSING,
    )
    if nested is not _MISSING:
        if not isinstance(nested, Mapping):
            raise _Rejected("invalid_admission")
        # Keep top-level mode/family/hash available while using the explicit
        # certificate for the receipt records.
        merged = dict(nested)
        for key in (
            "mode",
            "reference_mode",
            "source_mode",
            "source_family",
            "family",
            "graph_hash",
            "graph_sha256",
            "source_graph_hash",
            "source_graph_sha256",
            "binder_receipts",
            "integral_binder_receipts",
            "binder_identity_receipts",
            "integral_binder_identity",
            "free_variable_namespace",
            "free_variable_namespaces",
            "free_variable_receipts",
            "namespace_qualification",
        ):
            if key in admission:
                merged[key] = admission[key]
        return merged
    return admission


def _receipt_list(certificate: Mapping[str, Any], names: tuple[str, ...], code: str) -> list[Any]:
    value = _first(certificate, names, _MISSING)
    if value is _MISSING or not isinstance(value, (list, tuple)):
        raise _Rejected(code)
    return list(value)


def _receipt_reference(candidate: Mapping[str, Any]) -> Any:
    return _first(candidate, ("source_ref", "source_node_reference"), _MISSING)


def _match_binder_receipts(
    certificate: Mapping[str, Any], expected: list[dict[str, Any]]
) -> None:
    receipts = _receipt_list(
        certificate,
        (
            "binder_receipts",
            "integral_binder_receipts",
            "binder_identity_receipts",
            "integral_binder_identity",
        ),
        "binder_receipt_required",
    )
    if len(receipts) != len(expected):
        raise _Rejected("binder_receipt_mismatch")
    unmatched = list(receipts)
    for item in expected:
        matched = False
        for index, candidate in enumerate(unmatched):
            if not isinstance(candidate, Mapping):
                continue
            ref = _receipt_reference(candidate)
            name = _first(candidate, ("source_name", "name"), _MISSING)
            if (
                name == item["source_name"]
                and candidate.get("token") == item["token"]
                and ref is not _MISSING
                and _same_reference(ref, item["source_ref"])
            ):
                unmatched.pop(index)
                matched = True
                break
        if not matched:
            raise _Rejected("binder_receipt_mismatch")


def _match_namespace_receipts(
    certificate: Mapping[str, Any], expected: list[dict[str, Any]]
) -> None:
    receipts = _receipt_list(
        certificate,
        (
            "free_variable_namespace",
            "free_variable_namespaces",
            "free_variable_receipts",
            "namespace_qualification",
        ),
        "free_variable_namespace_required",
    )
    if len(receipts) != len(expected):
        raise _Rejected("free_variable_namespace_mismatch")
    unmatched = list(receipts)
    for item in expected:
        matched = False
        for index, candidate in enumerate(unmatched):
            if not isinstance(candidate, Mapping):
                continue
            name = _first(candidate, ("name", "variable"), _MISSING)
            ref = _receipt_reference(candidate)
            qualified = candidate.get("qualified", candidate.get("admitted"))
            namespace = candidate.get("namespace", candidate.get("namespace_id"))
            has_namespace = isinstance(namespace, str) and bool(namespace.strip())
            if (
                name == item["name"]
                and ref is not _MISSING
                and _same_reference(ref, item["source_ref"])
                and (qualified is True or has_namespace)
            ):
                unmatched.pop(index)
                matched = True
                break
        if not matched:
            raise _Rejected("free_variable_namespace_mismatch")


def _inspect_ir(ir: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Validate the supported v4 shape and collect opaque admission records."""

    integrals: list[dict[str, Any]] = []
    free_variables: list[dict[str, Any]] = []
    integral_tokens: set[str] = set()
    seen_free: list[tuple[str, Any]] = []
    active: set[int] = set()
    count = [0]

    def add_free(name: str, source_ref: Any) -> None:
        for old_name, old_ref in seen_free:
            if old_name == name and _same_reference(old_ref, source_ref):
                return
        seen_free.append((name, source_ref))
        free_variables.append({"name": name, "source_ref": source_ref})

    def visit(node: Any, depth: int = 0) -> None:
        if depth > _MAX_DEPTH:
            raise _Rejected("expression_depth_limit")
        if not isinstance(node, Mapping):
            raise _Rejected("expression_node_not_mapping")
        identity = id(node)
        if identity in active:
            raise _Rejected("cyclic_expression")
        active.add(identity)
        count[0] += 1
        if count[0] > _MAX_NODES:
            active.remove(identity)
            raise _Rejected("expression_node_limit")
        try:
            if node.get("kind") == "fresh_binder":
                token = node.get("token")
                source_name = node.get("source_name")
                if not isinstance(token, str) or not token:
                    raise _Rejected("invalid_fresh_binder")
                if not isinstance(source_name, str) or not source_name:
                    raise _Rejected("invalid_fresh_binder_name")
                return
            op = node.get("op")
            if op == "function_substitution":
                raise _Rejected("function_substitution_not_runtime")
            if op not in _ALLOWED_OPS:
                raise _Rejected("unknown_op", op)
            if op == "exact_rational":
                numerator = node.get("numerator", node.get("value", _MISSING))
                denominator = node.get("denominator", 1)
                if type(numerator) is not int or type(denominator) is not int:
                    raise _Rejected("invalid_exact_rational")
                if denominator == 0:
                    raise _Rejected("zero_exact_rational_denominator")
                return
            if op == "bound_variable":
                name = node.get("name", node.get("variable", _MISSING))
                if not isinstance(name, str) or not name:
                    raise _Rejected("invalid_bound_name")
                source_ref = node.get(
                    "source_node_reference", node.get("source_ref", _MISSING)
                )
                add_free(name, source_ref)
                return
            if op in ("sum", "product"):
                field = "terms" if op == "sum" else "factors"
                children = node.get(field)
                if not isinstance(children, (list, tuple)) or not children:
                    raise _Rejected("malformed_operator_children")
                for child in children:
                    visit(child, depth + 1)
                return
            if op in ("negative", "analytic_unary"):
                argument = node.get("argument")
                if not isinstance(argument, Mapping):
                    raise _Rejected("malformed_unary_argument")
                if op == "analytic_unary":
                    name = node.get("function")
                    if not isinstance(name, str) or name.lower() not in _ALLOWED_ANALYTIC_UNARY:
                        raise _Rejected("unknown_analytic_unary")
                visit(argument, depth + 1)
                return
            if op == "positive_quotient":
                numerator = node.get("numerator")
                denominator = node.get("denominator")
                if not isinstance(numerator, Mapping) or not isinstance(denominator, Mapping):
                    raise _Rejected("malformed_quotient")
                visit(numerator, depth + 1)
                visit(denominator, depth + 1)
                return
            if op == "definite_integral":
                binder = node.get("binder")
                if not isinstance(binder, Mapping) or binder.get("kind") != "fresh_binder":
                    raise _Rejected("integral_binder_required")
                token = binder.get("token")
                source_name = binder.get("source_name")
                source_ref = node.get("source_ref", _MISSING)
                if not isinstance(token, str) or not token:
                    raise _Rejected("invalid_integral_binder")
                if not isinstance(source_name, str) or not source_name:
                    raise _Rejected("integral_binder_name_required")
                if token in integral_tokens:
                    raise _Rejected("duplicate_fresh_binder")
                integral_tokens.add(token)
                integrals.append(
                    {
                        "source_ref": source_ref,
                        "token": token,
                        "source_name": source_name,
                    }
                )
                for field in ("lower", "upper", "integrand"):
                    child = node.get(field)
                    if not isinstance(child, Mapping):
                        raise _Rejected("malformed_integral")
                    visit(child, depth + 1)
                return
            if op == "current_original_source_parameter":
                return
            if op == "current_original_unique_analytic_core_component":
                arguments = node.get("arguments")
                if not isinstance(arguments, Mapping):
                    raise _Rejected("malformed_component_arguments")
                for field in ("rho", "Z_variable"):
                    child = arguments.get(field)
                    if not isinstance(child, Mapping):
                        raise _Rejected("malformed_component_arguments")
                    visit(child, depth + 1)
                return
            raise _Rejected("unknown_op", op)
        finally:
            active.remove(identity)

    visit(ir)
    return integrals, free_variables


def _expression_claim(metadata: Mapping[str, Any], names: tuple[str, ...]) -> Any:
    return _metadata_value(metadata, names)


class _FreeEnvironment:
    """Explicit interval bindings; never treats an integer key as a node index."""

    _CONTAINER_KEYS = {
        "values",
        "bindings",
        "by_name",
        "by_source",
        "source_aware",
        "lexical_certificate",
        "capture_certificate",
        "certificate",
    }

    def __init__(self, raw: Any, interval_type: Any):
        if raw is None:
            raw = {}
        if not isinstance(raw, (Mapping, list, tuple)):
            raise _Rejected("free_variables_not_mapping")
        self.interval_type = interval_type
        self.bindings: list[dict[str, Any]] = []
        if isinstance(raw, Mapping):
            self._load_mapping(raw)
        else:
            for item in raw:
                self._load_record(item)

    def _add(self, name: Any, source_ref: Any, value: Any) -> None:
        if not isinstance(name, str) or not name:
            raise _Rejected("invalid_free_variable_name")
        value = _require_interval(value, self.interval_type, "free_variable_not_interval")
        for item in self.bindings:
            if item["name"] == name and _same_reference(item["source_ref"], source_ref):
                if item["value"] != value:
                    raise _Rejected("duplicate_free_variable")
                return
        self.bindings.append({"name": name, "source_ref": source_ref, "value": value})

    def _load_record(self, record: Any) -> None:
        if not isinstance(record, Mapping):
            raise _Rejected("invalid_free_variable_record")
        name = _first(record, ("name", "variable"), _MISSING)
        source_ref = _first(record, ("source_ref", "source_node_reference"), None)
        value = _first(record, ("value", "interval", "binding"), _MISSING)
        if value is _MISSING:
            raise _Rejected("missing_free_variable_value")
        self._add(name, source_ref, value)

    def _load_named_mapping(self, mapping: Any) -> None:
        if not isinstance(mapping, Mapping):
            raise _Rejected("free_variables_not_mapping")
        for key, value in mapping.items():
            if isinstance(key, str):
                self._add(key, None, value)
                continue
            if isinstance(key, tuple) and len(key) == 2:
                left, right = key
                if isinstance(left, str):
                    self._add(left, right, value)
                elif isinstance(right, str):
                    self._add(right, left, value)
                else:
                    raise _Rejected("invalid_free_variable_key")
                continue
            # A bare integer/object key is deliberately not interpreted as a
            # graph-node index or source identity.
            raise _Rejected("invalid_free_variable_key")

    def _load_source_mapping(self, mapping: Any) -> None:
        if isinstance(mapping, Mapping):
            for source_ref, values in mapping.items():
                if not isinstance(values, Mapping):
                    raise _Rejected("invalid_source_namespace")
                for name, value in values.items():
                    if not isinstance(name, str):
                        raise _Rejected("invalid_free_variable_name")
                    self._add(name, source_ref, value)
            return
        if isinstance(mapping, (list, tuple)):
            for item in mapping:
                self._load_record(item)
            return
        raise _Rejected("invalid_source_namespace")

    def _load_mapping(self, raw: Mapping[str, Any]) -> None:
        for key in ("values", "bindings", "by_name"):
            value = raw.get(key, _MISSING)
            if value is not _MISSING:
                self._load_named_mapping(value)
        for key in ("by_source", "source_aware"):
            value = raw.get(key, _MISSING)
            if value is not _MISSING:
                self._load_source_mapping(value)
        # A plain ``{"x": interval}`` remains a useful synthetic shape.
        for key, value in raw.items():
            if key in self._CONTAINER_KEYS:
                continue
            if isinstance(key, str):
                self._add(key, None, value)
                continue
            if isinstance(key, tuple) and len(key) == 2:
                left, right = key
                if isinstance(left, str):
                    self._add(left, right, value)
                elif isinstance(right, str):
                    self._add(right, left, value)
                else:
                    raise _Rejected("invalid_free_variable_key")
                continue
            raise _Rejected("invalid_free_variable_key")

    def resolve(self, name: Any, source_ref: Any) -> Any:
        if not isinstance(name, str) or not name:
            raise _Rejected("invalid_bound_name")
        exact = [
            item
            for item in self.bindings
            if item["name"] == name
            and source_ref is not _MISSING
            and source_ref is not None
            and item["source_ref"] is not None
            and _same_reference(item["source_ref"], source_ref)
        ]
        if len(exact) == 1:
            return exact[0]["value"]
        if len(exact) > 1:
            raise _Rejected("duplicate_free_variable")
        named = [item for item in self.bindings if item["name"] == name]
        if len(named) == 1 and named[0]["source_ref"] is None:
            return named[0]["value"]
        raise _Rejected("unbound_variable")


class _Runtime:
    def __init__(
        self,
        ir: Mapping[str, Any],
        *,
        oracle: Any,
        backend: Any,
        interval_type: Any,
        extract_component: Any,
        endpoint_converter: Any,
        partitions: int,
        free_variables: Any,
    ) -> None:
        self.ir = ir
        self.oracle = oracle
        self.backend = backend
        self.interval_type = interval_type
        self.extract_component = extract_component
        self.endpoint_converter = endpoint_converter
        self.partitions = partitions
        self.free = _FreeEnvironment(free_variables, interval_type)
        self.active: set[int] = set()
        self.count = 0
        self.component_receipts: list[dict[str, Any]] = []
        self.parameter_receipts: list[dict[str, Any]] = []
        self.integral_receipts: list[dict[str, Any]] = []

    def _backend_member(self, names: tuple[str, ...]) -> Any:
        member = _member(self.backend, names)
        if member is _MISSING:
            return _MISSING
        return member

    def _backend_call(
        self, names: tuple[str, ...], args: tuple[Any, ...], missing_code: str
    ) -> Any:
        method = self._backend_member(names)
        if callable(method):
            try:
                return _require_interval(method(*args), self.interval_type, "backend_non_interval")
            except _Rejected:
                raise
            except Exception as exc:
                raise _Rejected("backend_error", type(exc).__name__) from exc
        raise _Rejected(missing_code)

    def _rational(self, numerator: Any, denominator: Any) -> Any:
        if type(numerator) is not int or type(denominator) is not int:
            raise _Rejected("invalid_exact_rational")
        if denominator == 0:
            raise _Rejected("zero_exact_rational_denominator")
        return self._backend_call(("rational",), (numerator, denominator), "backend_rational_required")

    def _add(self, left: Any, right: Any) -> Any:
        method = self._backend_member(("add", "plus"))
        if callable(method):
            try:
                return _require_interval(method(left, right), self.interval_type, "backend_non_interval")
            except _Rejected:
                raise
            except Exception as exc:
                raise _Rejected("backend_error", type(exc).__name__) from exc
        try:
            return _require_interval(left + right, self.interval_type, "arithmetic_non_interval")
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("backend_arithmetic_required", type(exc).__name__) from exc

    def _multiply(self, left: Any, right: Any) -> Any:
        method = self._backend_member(("multiply", "mul", "times"))
        if callable(method):
            try:
                return _require_interval(method(left, right), self.interval_type, "backend_non_interval")
            except _Rejected:
                raise
            except Exception as exc:
                raise _Rejected("backend_error", type(exc).__name__) from exc
        try:
            return _require_interval(left * right, self.interval_type, "arithmetic_non_interval")
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("backend_arithmetic_required", type(exc).__name__) from exc

    def _negative(self, value: Any) -> Any:
        method = self._backend_member(("negative", "neg"))
        if callable(method):
            try:
                return _require_interval(method(value), self.interval_type, "backend_non_interval")
            except _Rejected:
                raise
            except Exception as exc:
                raise _Rejected("backend_error", type(exc).__name__) from exc
        try:
            return _require_interval(-value, self.interval_type, "arithmetic_non_interval")
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("backend_arithmetic_required", type(exc).__name__) from exc

    def _divide(self, numerator: Any, denominator: Any) -> Any:
        method = self._backend_member(("divide", "div", "quotient"))
        if callable(method):
            try:
                return _require_interval(method(numerator, denominator), self.interval_type, "backend_non_interval")
            except _Rejected:
                raise
            except Exception as exc:
                raise _Rejected("backend_error", type(exc).__name__) from exc
        try:
            return _require_interval(numerator / denominator, self.interval_type, "arithmetic_non_interval")
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("backend_arithmetic_required", type(exc).__name__) from exc

    def _strict_positive(self, value: Any) -> bool:
        method = self._backend_member(("is_positive", "positive"))
        if callable(method):
            try:
                result = method(value)
            except Exception as exc:
                raise _Rejected("backend_positive_error", type(exc).__name__) from exc
            if type(result) is not bool:
                raise _Rejected("backend_positive_check_not_bool")
            return result
        try:
            result = value > 0
        except Exception as exc:
            raise _Rejected("backend_positive_check_required", type(exc).__name__) from exc
        if type(result) is not bool:
            raise _Rejected("backend_positive_check_not_bool")
        return result

    def _unary(self, name: str, value: Any) -> Any:
        method = self._backend_member(("unary",))
        if callable(method):
            try:
                return _require_interval(method(name, value), self.interval_type, "backend_non_interval")
            except _Rejected:
                raise
            except Exception as exc:
                raise _Rejected("backend_error", type(exc).__name__) from exc
        if name == "abs":
            method = self._backend_member(("abs", "absolute"))
        else:
            method = self._backend_member(("exp",))
        if not callable(method):
            raise _Rejected("backend_unary_required")
        try:
            return _require_interval(method(value), self.interval_type, "backend_non_interval")
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("backend_error", type(exc).__name__) from exc

    def _source_record(self, node: Mapping[str, Any]) -> dict[str, Any]:
        source_row = node.get("source_row", _MISSING)
        if not isinstance(source_row, Mapping):
            raise _Rejected("missing_source_row")
        record = dict(source_row)
        record["op"] = node.get("op")
        if "source_ref" in node:
            record["source_ref"] = node.get("source_ref")
        record["source_row"] = dict(source_row)
        for name in (
            "provider",
            "provider_refs",
            "provenance",
            "quantity",
            "Z_derivative",
            "source_family",
        ):
            if name in node:
                record[name] = node[name]
        record["opaque_source_record"] = True
        return record

    def _oracle_method(self, names: tuple[str, ...], code: str) -> Any:
        method = _member(self.oracle, names)
        if not callable(method):
            raise _Rejected(code)
        return method

    def _component(self, node: Mapping[str, Any], bindings: dict[str, Any], depth: int) -> Any:
        arguments = node.get("arguments")
        if not isinstance(arguments, Mapping):
            raise _Rejected("malformed_component_arguments")
        rho = self._eval(arguments.get("rho"), bindings, depth + 1)
        z_value = self._eval(arguments.get("Z_variable"), bindings, depth + 1)
        method = self._oracle_method(("component",), "oracle_component_required")
        record = self._source_record(node)
        try:
            component = method(record, rho, z_value)
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("oracle_component_error", type(exc).__name__) from exc
        if _member(component, ("kind",)) != "enclosure":
            raise _Rejected("oracle_component_not_enclosure")
        if _member(component, ("actual_value_selected",)) is not False:
            raise _Rejected("oracle_component_selected_point")
        try:
            extraction = self.extract_component(component, self.endpoint_converter)
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("component_extraction_error", type(exc).__name__) from exc
        if _member(extraction, ("kind",)) != "conditional_enclosure":
            raise _Rejected("component_extraction_not_conditional")
        if _member(extraction, ("actual_value_selected",)) is not False:
            raise _Rejected("component_extraction_selected_point")
        interval = _member(extraction, ("interval",))
        interval = _require_interval(interval, self.interval_type, "component_interval_not_exact")
        extraction_receipt: dict[str, Any] = {
            "schema": SCHEMA,
            "source_ref": node.get("source_ref", None),
            "record": record,
            "rho": rho,
            "Z": z_value,
            "component_kind": "enclosure",
            "component_actual_value_selected": False,
            "extraction_kind": "conditional_enclosure",
            "extracted_interval": interval,
            "endpoint_conversion_verified": False,
            "source_identity_verified": False,
            "actual_value_selected": False,
            "certified": False,
        }
        component_receipt = _member(component, ("receipt", "last_receipt"))
        if component_receipt is not _MISSING:
            component_receipt = _call_noarg(component_receipt, "component_receipt_error")
            if isinstance(component_receipt, Mapping):
                extraction_receipt["component_receipt"] = component_receipt
        extraction_receipt_member = _member(extraction, ("receipt",))
        if extraction_receipt_member is not _MISSING:
            extraction_receipt_member = _call_noarg(
                extraction_receipt_member, "component_extraction_receipt_error"
            )
            if isinstance(extraction_receipt_member, Mapping):
                extraction_receipt["extraction_receipt"] = extraction_receipt_member
        self.component_receipts.append(extraction_receipt)
        return interval

    def _integrator_method(self) -> Any:
        for owner in (self.backend, self.oracle):
            method = _member(owner, ("integrate", "range_integrate"))
            if callable(method):
                return method
        raise _Rejected("range_integrator_required")

    def _call_integrator(self, method: Any, callback: Any, lower: Any, upper: Any) -> Any:
        """Invoke one of the explicit supported integrator signatures."""

        try:
            signature = inspect.signature(method)
        except (TypeError, ValueError):
            signature = None
        try:
            if signature is not None:
                params = list(signature.parameters.values())
                if "partitions" in signature.parameters or any(
                    item.kind is inspect.Parameter.VAR_KEYWORD for item in params
                ):
                    return method(callback, lower, upper, partitions=self.partitions)
                positional = [
                    item
                    for item in params
                    if item.kind
                    in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD)
                ]
                if any(item.kind is inspect.Parameter.VAR_POSITIONAL for item in params) or len(positional) >= 4:
                    return method(callback, lower, upper, self.partitions)
                return method(callback, lower, upper)
            return method(callback, lower, upper, self.partitions)
        except _Rejected:
            raise
        except Exception as exc:
            raise _Rejected("integrator_error", type(exc).__name__) from exc

    def _integral_result(self, result: Any, node: Mapping[str, Any], lower: Any, upper: Any, token: str) -> Any:
        if _member(result, ("kind",)) != "conditional_enclosure":
            raise _Rejected("integral_not_conditional_enclosure")
        if _member(result, ("actual_value_selected",)) is not False:
            raise _Rejected("integral_selected_point")
        interval = _member(result, ("interval",))
        interval = _require_interval(interval, self.interval_type, "integral_interval_not_exact")
        cells = _member(result, ("coverage_cells", "cells"))
        if not isinstance(cells, (list, tuple)):
            raise _Rejected("integral_coverage_required")
        checked_cells = []
        for cell in cells:
            checked_cells.append(_require_interval(cell, self.interval_type, "integral_cell_not_exact"))
        receipt: dict[str, Any] = {
            "schema": SCHEMA,
            "source_ref": node.get("source_ref", None),
            "token": token,
            "source_name": node.get("binder", {}).get("source_name"),
            "lower": lower,
            "upper": upper,
            "interval": interval,
            "coverage_cells": tuple(checked_cells),
            "coverage_count": len(checked_cells),
            "actual_value_selected": False,
            "certified": False,
            "coverage_verified": False,
        }
        receipt_member = _member(result, ("receipt",))
        if receipt_member is not _MISSING:
            receipt_member = _call_noarg(receipt_member, "integral_receipt_error")
            if isinstance(receipt_member, Mapping):
                receipt["integrator_receipt"] = receipt_member
        self.integral_receipts.append(receipt)
        return interval

    def _eval_integral(self, node: Mapping[str, Any], bindings: dict[str, Any], depth: int) -> Any:
        # Bounds are evaluated before introducing the fresh token, in the
        # captured outer environment.  No endpoint is converted to a point.
        lower = self._eval(node.get("lower"), bindings, depth + 1)
        upper = self._eval(node.get("upper"), bindings, depth + 1)
        binder = node.get("binder")
        if not isinstance(binder, Mapping) or binder.get("kind") != "fresh_binder":
            raise _Rejected("integral_binder_required")
        token = binder.get("token")
        if not isinstance(token, str) or not token:
            raise _Rejected("invalid_integral_binder")
        if token in bindings:
            raise _Rejected("binder_shadow")
        integrand = node.get("integrand")
        if not isinstance(integrand, Mapping):
            raise _Rejected("malformed_integral")
        captured = dict(bindings)

        def cell_callback(cell: Any) -> Any:
            cell = _require_interval(cell, self.interval_type, "integrator_cell_not_exact")
            inner = dict(captured)
            if token in inner:
                raise _Rejected("binder_shadow")
            inner[token] = cell
            return self._eval(integrand, inner, depth + 1)

        method = self._integrator_method()
        result = self._call_integrator(method, cell_callback, lower, upper)
        return self._integral_result(result, node, lower, upper, token)

    def _eval(self, node: Any, bindings: dict[str, Any], depth: int = 0) -> Any:
        if depth > _MAX_DEPTH:
            raise _Rejected("expression_depth_limit")
        if not isinstance(node, Mapping):
            raise _Rejected("expression_node_not_mapping")
        identity = id(node)
        if identity in self.active:
            raise _Rejected("cyclic_expression")
        self.active.add(identity)
        self.count += 1
        if self.count > _MAX_NODES:
            self.active.remove(identity)
            raise _Rejected("expression_node_limit")
        try:
            if node.get("kind") == "fresh_binder":
                token = node.get("token")
                if not isinstance(token, str) or not token:
                    raise _Rejected("invalid_fresh_binder")
                if token not in bindings:
                    raise _Rejected("unbound_fresh_binder")
                return _require_interval(bindings[token], self.interval_type, "invalid_bound_value")
            op = node.get("op")
            if op == "function_substitution":
                raise _Rejected("function_substitution_not_runtime")
            if op not in _ALLOWED_OPS:
                raise _Rejected("unknown_op", op)
            if op == "exact_rational":
                numerator = node.get("numerator", node.get("value", _MISSING))
                denominator = node.get("denominator", 1)
                if numerator is _MISSING:
                    raise _Rejected("invalid_exact_rational")
                return self._rational(numerator, denominator)
            if op == "bound_variable":
                name = node.get("name", node.get("variable", _MISSING))
                source_ref = node.get(
                    "source_node_reference", node.get("source_ref", _MISSING)
                )
                return self.free.resolve(name, source_ref)
            if op == "sum":
                terms = node.get("terms")
                if not isinstance(terms, (list, tuple)) or not terms:
                    raise _Rejected("malformed_sum")
                result = self._eval(terms[0], bindings, depth + 1)
                for term in terms[1:]:
                    result = self._add(result, self._eval(term, bindings, depth + 1))
                return result
            if op == "product":
                factors = node.get("factors")
                if not isinstance(factors, (list, tuple)) or not factors:
                    raise _Rejected("malformed_product")
                result = self._eval(factors[0], bindings, depth + 1)
                for factor in factors[1:]:
                    result = self._multiply(result, self._eval(factor, bindings, depth + 1))
                return result
            if op == "negative":
                return self._negative(self._eval(node.get("argument"), bindings, depth + 1))
            if op == "positive_quotient":
                numerator = self._eval(node.get("numerator"), bindings, depth + 1)
                denominator = self._eval(node.get("denominator"), bindings, depth + 1)
                if not self._strict_positive(denominator):
                    raise _Rejected("quotient_denominator_not_strictly_positive")
                return self._divide(numerator, denominator)
            if op == "analytic_unary":
                name = node.get("function")
                if not isinstance(name, str) or name.lower() not in _ALLOWED_ANALYTIC_UNARY:
                    raise _Rejected("unknown_analytic_unary")
                return self._unary(name.lower(), self._eval(node.get("argument"), bindings, depth + 1))
            if op == "definite_integral":
                return self._eval_integral(node, bindings, depth)
            if op == "current_original_source_parameter":
                record = self._source_record(node)
                method = self._oracle_method(("parameter",), "oracle_parameter_required")
                try:
                    value = method(record)
                except _Rejected:
                    raise
                except Exception as exc:
                    raise _Rejected("oracle_parameter_error", type(exc).__name__) from exc
                value = _require_interval(value, self.interval_type, "oracle_parameter_not_interval")
                self.parameter_receipts.append(
                    {
                        "source_ref": node.get("source_ref", None),
                        "interval": value,
                        "actual_value_selected": False,
                    }
                )
                return value
            if op == "current_original_unique_analytic_core_component":
                return self._component(node, bindings, depth)
            raise _Rejected("unknown_op", op)
        finally:
            self.active.remove(identity)

    def run(self) -> Any:
        return self._eval(self.ir, {}, 0)


def _base_result(
    status: str,
    *,
    mode: str | None,
    source_family: Any = None,
    graph_hash: Any = None,
    value: Any = None,
    error: dict[str, Any] | None = None,
    component_receipts: list[Any] | None = None,
    parameter_receipts: list[Any] | None = None,
    integral_receipts: list[Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": status,
        "ok": status == "conditional_enclosure",
        "value": None,
        "interval": value if status == "conditional_enclosure" else None,
        "mode": mode,
        "source_family": source_family,
        "graph_hash": graph_hash,
        "actual_value_selected": False,
        "certified": False,
        "provenance": {
            "caller_provenance": "not independently verified",
            "source_identity_verified": False,
            "endpoint_conversion_verified": False,
            "binder_identity_verified": False,
            "bounds_verified": False,
            "quadrature_verified": False,
            "certified": False,
        },
        "component_extractions": component_receipts or [],
        "component_extraction_receipts": component_receipts or [],
        "component_receipts": component_receipts or [],
        "parameter_receipts": parameter_receipts or [],
        "integral_receipts": integral_receipts or [],
        "integral_coverage_receipts": integral_receipts or [],
        "error": error,
        "limitations": [
            "Caller-supplied source, endpoint, binder, and quadrature claims are not independently verified.",
            "The interval is conditional on every injected callback returning a valid enclosure.",
            "No source graph node or integer graph index is resolved by this runtime.",
            "A point value, midpoint, tail cap, or raw scalar is rejected at interval boundaries.",
        ],
    }


def _admit(
    ir: Mapping[str, Any], metadata: Mapping[str, Any], report: Mapping[str, Any], admission: Any
) -> tuple[str, Any, Any, list[dict[str, Any]], list[dict[str, Any]]]:
    cert = _admission_container(admission) if admission is not None else {}
    mode_admission = _first(cert, ("mode", "reference_mode", "source_mode"), _MISSING)
    mode_oracle = _first(report, ("mode", "reference_mode", "source_mode", "evaluation_mode"), _MISSING)
    mode_expression = _metadata_value(metadata, ("mode", "reference_mode", "source_mode", "evaluation_mode"))
    mode_value = _MISSING
    for candidate in (mode_admission, mode_oracle, mode_expression):
        if candidate is _MISSING:
            continue
        mode_value = _merge_values(mode_value, candidate, "oracle_mode_mismatch")
    mode = _normalise_mode(mode_value)
    if mode not in ("synthetic_reference", "original_source"):
        raise _Rejected("oracle_mode_required")
    integrals, free_variables = _inspect_ir(ir)
    if mode == "synthetic_reference":
        return mode, None, None, integrals, free_variables

    expression_family = _expression_claim(metadata, ("source_family", "family"))
    expression_hash = _expression_claim(
        metadata,
        ("graph_hash", "graph_sha256", "source_graph_hash", "source_graph_sha256"),
    )
    cert_family = _first(cert, ("source_family", "family"), _MISSING)
    cert_hash = _first(
        cert,
        ("graph_hash", "graph_sha256", "source_graph_hash", "source_graph_sha256"),
        _MISSING,
    )
    oracle_family = _first(report, ("source_family", "family"), _MISSING)
    oracle_hash = _first(
        report,
        ("graph_hash", "graph_sha256", "source_graph_hash", "source_graph_sha256"),
        _MISSING,
    )
    if expression_family is _MISSING:
        raise _Rejected("expression_source_family_required")
    if expression_hash is _MISSING:
        raise _Rejected("expression_graph_hash_required")
    if cert_family is _MISSING or cert_hash is _MISSING:
        raise _Rejected("admission_source_identity_required")
    if mode_oracle is _MISSING or oracle_family is _MISSING or oracle_hash is _MISSING:
        raise _Rejected("oracle_source_identity_required")
    if not _valid_family(expression_family) or not _valid_graph_hash(expression_hash):
        raise _Rejected("source_identity_format_invalid")
    if cert_family != expression_family or cert_hash != expression_hash:
        raise _Rejected("admission_source_identity_mismatch")
    if oracle_family != expression_family:
        raise _Rejected("oracle_source_family_mismatch")
    if oracle_hash != expression_hash:
        raise _Rejected("oracle_graph_hash_mismatch")
    for record in integrals:
        if record["source_ref"] is _MISSING or record["source_ref"] is None:
            raise _Rejected("integral_source_ref_required")
    for record in free_variables:
        if record["source_ref"] is _MISSING or record["source_ref"] is None:
            raise _Rejected("free_variable_source_ref_required")
    _match_binder_receipts(cert, integrals)
    _match_namespace_receipts(cert, free_variables)
    return mode, expression_family, expression_hash, integrals, free_variables


def evaluate_interval(
    ir_or_lowering: Any,
    *,
    oracle: Any,
    backend: Any,
    interval_type: Any,
    extract_component: Any,
    endpoint_converter: Any,
    partitions: int,
    free_variables: Any = None,
    admission: Any = None,
) -> dict[str, Any]:
    """Evaluate one lowered v4 IR expression using exact injected intervals.

    ``backend`` must expose exact interval arithmetic and either it or the
    oracle must expose ``integrate(f, lo, hi, partitions)`` (keyword or
    positional ``partitions`` is accepted).  A successful result exposes its
    tagged enclosure under ``interval``; ``value`` remains None.
    """

    mode: str | None = None
    source_family: Any = None
    graph_hash: Any = None
    try:
        if oracle is None:
            raise _Rejected("oracle_required")
        if backend is None:
            raise _Rejected("backend_required")
        if not isinstance(interval_type, type):
            raise _Rejected("interval_type_required")
        if not callable(extract_component):
            raise _Rejected("extract_component_required")
        if not callable(endpoint_converter):
            raise _Rejected("endpoint_converter_required")
        if type(partitions) is not int or partitions <= 0 or partitions > _MAX_PARTITIONS:
            raise _Rejected("partitions_invalid")
        ir, metadata = _unwrap_expression(ir_or_lowering)
        report = _oracle_report(oracle)
        mode, source_family, graph_hash, _integrals, _free = _admit(
            ir, metadata, report, admission
        )
        runtime = _Runtime(
            ir,
            oracle=oracle,
            backend=backend,
            interval_type=interval_type,
            extract_component=extract_component,
            endpoint_converter=endpoint_converter,
            partitions=partitions,
            free_variables=free_variables,
        )
        value = runtime.run()
        value = _require_interval(value, interval_type, "runtime_result_not_interval")
        return _base_result(
            "conditional_enclosure",
            mode=mode,
            source_family=source_family,
            graph_hash=graph_hash,
            value=value,
            component_receipts=runtime.component_receipts,
            parameter_receipts=runtime.parameter_receipts,
            integral_receipts=runtime.integral_receipts,
        )
    except _Rejected as exc:
        error = {"code": exc.code}
        if exc.detail is not None:
            error["detail"] = exc.detail
        return _base_result(
            "rejected",
            mode=mode,
            source_family=source_family,
            graph_hash=graph_hash,
            error=error,
        )
    except Exception as exc:
        # Caller callbacks are outside this module's trust boundary.  Never
        # expose a partial accumulator or turn their failure into a value.
        return _base_result(
            "rejected",
            mode=mode,
            source_family=source_family,
            graph_hash=graph_hash,
            error={"code": "runtime_error", "type": type(exc).__name__},
        )


evaluate = evaluate_interval


__all__ = ["SCHEMA", "evaluate_interval", "evaluate"]
