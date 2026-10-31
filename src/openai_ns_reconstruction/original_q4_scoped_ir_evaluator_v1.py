"""Capture-safe runtime for the source-preserving Q4 scoped IR.

The v4 lowering module intentionally stops at an expression tree.  This
module is the small runtime boundary for that tree.  It does not open the
archived graph, inspect a graph handle, read a file, or infer a numerical
value from a source reference.  A caller must supply an oracle and a backend.

``evaluate`` accepts either the raw ``ir`` dictionary or a lowered-result
envelope containing an ``ir`` member.  The oracle must explicitly report one
of ``synthetic_reference`` and ``original_source``.  The latter additionally
requires source-family, graph-hash, and capture-safe binder admission
evidence in the expression/oracle metadata.  All successful and rejected
results deliberately say that caller provenance and numerical errors were
not independently verified.

The supported runtime nodes are the nodes emitted by
``original_q4_moment_lowering_v4``: exact rationals, sums, products,
negation, strictly-positive quotients, restricted analytic unary calls,
explicit bound variables, fresh binder tokens, definite integrals, and the
two opaque source records.  ``function_substitution`` is rejected because v4
has already lowered it.
"""

from collections.abc import Mapping


_MAX_DEPTH = 128
_MAX_NODES = 1024

_ALLOWED_OPS = {
    "sum",
    "product",
    "negative",
    "positive_quotient",
    "analytic_unary",
    "exact_rational",
    "bound_variable",
    "definite_integral",
    "current_original_source_parameter",
    "current_original_unique_analytic_core_component",
}

_ALLOWED_ANALYTIC_UNARY = {
    "abs",
    "cos",
    "cosh",
    "erf",
    "exp",
    "expm1",
    "log",
    "log1p",
    "sin",
    "sinh",
    "sqrt",
    "tan",
    "tanh",
}

_MISSING = object()


class _Rejected(Exception):
    """Internal fail-closed signal."""

    def __init__(self, code, detail=None):
        self.code = code
        self.detail = detail
        super().__init__(code)


def _member(obj, names):
    """Read one explicitly named member from an object or mapping."""

    if isinstance(obj, Mapping):
        for name in names:
            if name in obj:
                return obj[name]
        return _MISSING
    for name in names:
        try:
            return getattr(obj, name)
        except Exception:
            continue
    return _MISSING


def _clone_json(value, depth=0, active=None):
    """Clone metadata without importing or invoking an archive helper."""

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
                key: _clone_json(item, depth + 1, active)
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
            return [_clone_json(item, depth + 1, active) for item in value]
        finally:
            active.remove(identity)
    if isinstance(value, tuple):
        identity = id(value)
        if identity in active:
            raise _Rejected("cyclic_metadata")
        active.add(identity)
        try:
            return tuple(_clone_json(item, depth + 1, active) for item in value)
        finally:
            active.remove(identity)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    # Source rows from v4 are JSON-like.  Do not silently stringify an
    # arbitrary source/provider object, because that would make provenance
    # look numeric or authenticated when it is neither.
    raise _Rejected("non_json_source_record")


def _first(mapping, names, default=_MISSING):
    if not isinstance(mapping, Mapping):
        return default
    for name in names:
        if name in mapping:
            return mapping[name]
    return default


def _same_reference(left, right):
    """Compare source references without int/bool coercion."""

    if type(left) is not type(right):
        return False
    try:
        return bool(left == right)
    except Exception:
        return False


def _is_scalar(value):
    """Reject records, containers, and textual values at numeric boundaries."""

    if value is None or isinstance(value, bool):
        return False
    if isinstance(value, (Mapping, list, tuple, set, frozenset)):
        return False
    if isinstance(value, (str, bytes, bytearray)):
        return False
    return True


def _ensure_scalar(value, code):
    if not _is_scalar(value):
        raise _Rejected(code)
    return value


def _call_noarg(value, code):
    if not callable(value):
        return value
    try:
        return value()
    except _Rejected:
        raise
    except Exception:
        raise _Rejected(code)


def _oracle_report(oracle):
    """Read only the explicitly documented self-report members."""

    if isinstance(oracle, (str, bytes, bytearray, list, tuple, set, frozenset)):
        raise _Rejected("invalid_oracle")
    report = {}

    # A report method is useful for an object that keeps the admission data
    # private.  No arbitrary introspection is performed.
    for name in ("reference_report", "mode_report", "admission_report"):
        member = _member(oracle, (name,))
        if member is _MISSING:
            continue
        value = _call_noarg(member, "oracle_report_error")
        if isinstance(value, Mapping):
            report.update(value)

    direct_report = _member(oracle, ("report",))
    if direct_report is not _MISSING:
        value = _call_noarg(direct_report, "oracle_report_error")
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
        "admission_certificate",
        "capture_safe_admission_certificate",
        "capture_safe_certificate",
        "binder_identity_certificate",
    ):
        member = _member(oracle, (name,))
        if member is not _MISSING:
            report.setdefault(name, _call_noarg(member, "oracle_report_error"))
    return report


def _normalise_mode(value):
    if not isinstance(value, str):
        return None
    return value.strip().lower()


def _unwrap_expression(expression):
    """Return ``(ir, metadata)`` without interpreting metadata as an IR node."""

    if not isinstance(expression, Mapping):
        raise _Rejected("expression_not_mapping")
    if "op" not in expression and "ir" in expression:
        ir = expression.get("ir")
        if not isinstance(ir, Mapping):
            raise _Rejected("missing_ir")
        return ir, expression
    return expression, expression


def _metadata_value(metadata, names):
    if not isinstance(metadata, Mapping):
        return _MISSING
    for name in names:
        if name in metadata:
            return metadata[name]
    return _MISSING


def _report_value(report, names):
    return _metadata_value(report, names)


def _merge_explicit(metadata, report, names, code):
    left = _metadata_value(metadata, names)
    right = _report_value(report, names)
    if left is not _MISSING and right is not _MISSING:
        if left != right:
            raise _Rejected(code)
    value = left if left is not _MISSING else right
    return value


def _valid_family(value):
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, Mapping):
        return bool(value)
    return False


def _valid_graph_hash(value):
    # The runtime cannot independently recompute the archive hash.  It does
    # require a caller-supplied non-empty string and compares it exactly.
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def _certificate_from(metadata, report):
    names = (
        "admission_certificate",
        "capture_safe_admission_certificate",
        "capture_safe_certificate",
        "binder_identity_certificate",
    )
    candidates = []
    for container in (metadata, report):
        value = _metadata_value(container, names)
        if value is not _MISSING:
            candidates.append(value)
    if not candidates:
        # A self-report may put the certificate fields at its top level.
        if (
            _report_value(report, ("capture_safe", "capture_safe_for_this_ir"))
            is not _MISSING
            and _report_value(
                report,
                ("binder_identity", "binder_identities", "binder_identity_map"),
            )
            is not _MISSING
        ):
            return report
        return _MISSING
    certificate = candidates[0]
    if not isinstance(certificate, Mapping):
        raise _Rejected("invalid_admission_certificate")
    for other in candidates[1:]:
        if not isinstance(other, Mapping) or other != certificate:
            raise _Rejected("admission_certificate_mismatch")
    return certificate


def _fresh_tokens(ir):
    """Collect fresh tokens from expression structure, excluding opaque rows."""

    tokens = []
    active = set()
    count = [0]

    def visit(node, depth=0):
        if depth > _MAX_DEPTH:
            raise _Rejected("expression_depth_limit")
        if isinstance(node, Mapping):
            identity = id(node)
            if identity in active:
                raise _Rejected("cyclic_expression")
            active.add(identity)
            try:
                count[0] += 1
                if count[0] > _MAX_NODES:
                    raise _Rejected("expression_node_limit")
                if node.get("kind") == "fresh_binder":
                    token = node.get("token")
                    if not isinstance(token, str) or not token:
                        raise _Rejected("invalid_fresh_binder")
                    # A token may occur many times in one integrand.  It is
                    # the integral declaration, checked by
                    # ``_integral_receipts``, that must be unique.
                    if token not in tokens:
                        tokens.append(token)
                    return
                op = node.get("op")
                if op is None:
                    raise _Rejected("invalid_expression_node")
                if op not in _ALLOWED_OPS:
                    raise _Rejected("unknown_op")
                if op in ("sum", "product"):
                    field = "terms" if op == "sum" else "factors"
                    children = node.get(field)
                    if not isinstance(children, (list, tuple)) or not children:
                        raise _Rejected("malformed_operator_children")
                    for child in children:
                        visit(child, depth + 1)
                elif op in ("negative", "analytic_unary"):
                    child = node.get("argument")
                    if not isinstance(child, Mapping):
                        raise _Rejected("malformed_unary_argument")
                    visit(child, depth + 1)
                elif op == "positive_quotient":
                    for field in ("numerator", "denominator"):
                        child = node.get(field)
                        if not isinstance(child, Mapping):
                            raise _Rejected("malformed_quotient")
                        visit(child, depth + 1)
                elif op == "definite_integral":
                    binder = node.get("binder")
                    if not isinstance(binder, Mapping):
                        raise _Rejected("integral_binder_required")
                    visit(binder, depth + 1)
                    for field in ("lower", "upper", "integrand"):
                        child = node.get(field)
                        if not isinstance(child, Mapping):
                            raise _Rejected("malformed_integral")
                        visit(child, depth + 1)
                elif op == "current_original_unique_analytic_core_component":
                    arguments = node.get("arguments")
                    if not isinstance(arguments, Mapping):
                        raise _Rejected("malformed_component_arguments")
                    for field in ("rho", "Z_variable"):
                        child = arguments.get(field)
                        if not isinstance(child, Mapping):
                            raise _Rejected("malformed_component_arguments")
                        visit(child, depth + 1)
                # exact_rational, bound_variable, and source_parameter are
                # leaves.  Their source_row is opaque and is not traversed.
            finally:
                active.remove(identity)
        else:
            raise _Rejected("expression_node_not_mapping")

    visit(ir)
    return tokens


def _binder_identity_matches(certificate, tokens):
    capture = _first(
        certificate,
        ("capture_safe", "capture_safe_for_this_ir", "capture_safe_admission"),
        _MISSING,
    )
    if capture is not True:
        return False
    identity = _first(
        certificate,
        (
            "binder_identity",
            "binder_identities",
            "binder_identity_map",
            "binder_identity_certificate",
        ),
        _MISSING,
    )
    if identity is _MISSING:
        return False
    if isinstance(identity, Mapping):
        return all(token in identity or str(token) in identity for token in tokens)
    if isinstance(identity, (list, tuple, set, frozenset)):
        listed = set()
        for item in identity:
            if isinstance(item, str):
                listed.add(item)
            elif isinstance(item, Mapping) and isinstance(item.get("token"), str):
                listed.add(item["token"])
        return all(token in listed for token in tokens)
    return False


def _check_source_rows(ir, family):
    """Check explicit source-family annotations without reading a graph."""

    active = set()

    def visit(node, depth=0):
        if depth > _MAX_DEPTH:
            raise _Rejected("expression_depth_limit")
        if not isinstance(node, Mapping):
            return
        identity = id(node)
        if identity in active:
            raise _Rejected("cyclic_expression")
        active.add(identity)
        try:
            op = node.get("op")
            source_row = node.get("source_row")
            if isinstance(source_row, Mapping):
                row_family = _first(source_row, ("source_family", "family"), _MISSING)
                if row_family is not _MISSING and row_family != family:
                    raise _Rejected("source_family_mismatch")
            if node.get("kind") == "fresh_binder":
                return
            if op in ("sum", "product"):
                field = "terms" if op == "sum" else "factors"
                for child in node.get(field, ()):
                    visit(child, depth + 1)
            elif op in ("negative", "analytic_unary"):
                visit(node.get("argument"), depth + 1)
            elif op == "positive_quotient":
                visit(node.get("numerator"), depth + 1)
                visit(node.get("denominator"), depth + 1)
            elif op == "definite_integral":
                visit(node.get("lower"), depth + 1)
                visit(node.get("upper"), depth + 1)
                visit(node.get("integrand"), depth + 1)
            elif op == "current_original_unique_analytic_core_component":
                arguments = node.get("arguments", {})
                visit(arguments.get("rho"), depth + 1)
                visit(arguments.get("Z_variable"), depth + 1)
        finally:
            active.remove(identity)

    visit(ir)


def _integral_receipts(ir):
    """Collect integral binder relations without indexing by source refs."""

    records = []
    tokens = set()
    active = set()

    def visit(node, depth=0):
        if depth > _MAX_DEPTH:
            raise _Rejected("expression_depth_limit")
        if not isinstance(node, Mapping):
            return
        identity = id(node)
        if identity in active:
            raise _Rejected("cyclic_expression")
        active.add(identity)
        try:
            op = node.get("op")
            if op == "definite_integral":
                binder = node.get("binder")
                if not isinstance(binder, Mapping) or binder.get("kind") != "fresh_binder":
                    raise _Rejected("integral_binder_required")
                source_ref = node.get("source_ref", _MISSING)
                token = binder.get("token")
                source_name = binder.get("source_name")
                if source_ref is _MISSING or source_ref is None:
                    raise _Rejected("integral_source_ref_required")
                if not isinstance(token, str) or not token:
                    raise _Rejected("invalid_integral_binder")
                if not isinstance(source_name, str) or not source_name:
                    raise _Rejected("integral_binder_name_required")
                if token in tokens:
                    raise _Rejected("duplicate_fresh_binder")
                tokens.add(token)
                records.append(
                    {
                        "source_ref": _clone_json(source_ref),
                        "token": token,
                        "source_name": source_name,
                    }
                )
                visit(node.get("lower"), depth + 1)
                visit(node.get("upper"), depth + 1)
                visit(node.get("integrand"), depth + 1)
                return
            if node.get("kind") == "fresh_binder":
                return
            if op in ("sum", "product"):
                field = "terms" if op == "sum" else "factors"
                for child in node.get(field, ()):
                    visit(child, depth + 1)
            elif op in ("negative", "analytic_unary"):
                visit(node.get("argument"), depth + 1)
            elif op == "positive_quotient":
                visit(node.get("numerator"), depth + 1)
                visit(node.get("denominator"), depth + 1)
            elif op == "current_original_unique_analytic_core_component":
                arguments = node.get("arguments", {})
                visit(arguments.get("rho"), depth + 1)
                visit(arguments.get("Z_variable"), depth + 1)
        finally:
            active.remove(identity)

    visit(ir)
    return records


def _free_variable_records(ir):
    """Collect unresolved bound-variable leaves as opaque name/ref records."""

    records = []
    active = set()

    def add(record):
        for existing in records:
            if (
                existing["name"] == record["name"]
                and _same_reference(existing["source_ref"], record["source_ref"])
            ):
                return
        records.append(record)

    def visit(node, depth=0):
        if depth > _MAX_DEPTH:
            raise _Rejected("expression_depth_limit")
        if not isinstance(node, Mapping):
            return
        identity = id(node)
        if identity in active:
            raise _Rejected("cyclic_expression")
        active.add(identity)
        try:
            if node.get("kind") == "fresh_binder":
                return
            op = node.get("op")
            if op == "bound_variable":
                name = node.get("name", node.get("variable", _MISSING))
                source_ref = node.get(
                    "source_node_reference",
                    node.get("source_ref", _MISSING),
                )
                if not isinstance(name, str) or not name:
                    raise _Rejected("invalid_bound_name")
                if source_ref is _MISSING or source_ref is None:
                    raise _Rejected("free_variable_source_ref_required")
                add({"name": name, "source_ref": _clone_json(source_ref)})
                return
            if op in ("sum", "product"):
                field = "terms" if op == "sum" else "factors"
                for child in node.get(field, ()):
                    visit(child, depth + 1)
            elif op in ("negative", "analytic_unary"):
                visit(node.get("argument"), depth + 1)
            elif op == "positive_quotient":
                visit(node.get("numerator"), depth + 1)
                visit(node.get("denominator"), depth + 1)
            elif op == "definite_integral":
                visit(node.get("lower"), depth + 1)
                visit(node.get("upper"), depth + 1)
                visit(node.get("integrand"), depth + 1)
            elif op == "current_original_unique_analytic_core_component":
                arguments = node.get("arguments", {})
                visit(arguments.get("rho"), depth + 1)
                visit(arguments.get("Z_variable"), depth + 1)
        finally:
            active.remove(identity)

    visit(ir)
    return records


def _receipt_collection(certificate, names, error_code):
    value = _first(certificate, names, _MISSING)
    if value is _MISSING:
        raise _Rejected(error_code)
    if not isinstance(value, (list, tuple)):
        # A list of records forces opaque source refs to be compared as
        # values; using a source ref as a dictionary key would invite an
        # accidental graph-index interpretation.
        raise _Rejected(error_code)
    return value


def _receipt_matches(candidate, expected):
    if not isinstance(candidate, Mapping):
        return False
    if candidate.get("name", candidate.get("source_name")) != expected.get(
        "source_name", expected.get("name")
    ):
        return False
    candidate_token = candidate.get("token")
    if expected.get("token") is not None and candidate_token != expected.get("token"):
        return False
    candidate_ref = candidate.get("source_ref", candidate.get("source_node_reference", _MISSING))
    expected_ref = expected.get("source_ref", _MISSING)
    if candidate_ref is _MISSING or expected_ref is _MISSING:
        return False
    return _same_reference(candidate_ref, expected_ref)


def _binder_receipts_match(certificate, expected):
    receipts = _receipt_collection(
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
        return False
    unmatched = list(receipts)
    for record in expected:
        for index, candidate in enumerate(unmatched):
            if _receipt_matches(candidate, record):
                unmatched.pop(index)
                break
        else:
            return False
    return not unmatched


def _namespace_receipts_match(certificate, expected):
    receipts = _receipt_collection(
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
        return False
    unmatched = list(receipts)
    for record in expected:
        matched = False
        for index, candidate in enumerate(unmatched):
            if not isinstance(candidate, Mapping):
                continue
            candidate_name = candidate.get("name", candidate.get("variable"))
            candidate_ref = candidate.get(
                "source_ref",
                candidate.get("source_node_reference", _MISSING),
            )
            qualified = candidate.get("qualified", candidate.get("admitted"))
            namespace = candidate.get("namespace", candidate.get("namespace_id"))
            has_namespace = isinstance(namespace, str) and bool(namespace.strip())
            if (
                candidate_name == record["name"]
                and candidate_ref is not _MISSING
                and _same_reference(candidate_ref, record["source_ref"])
                and (qualified is True or has_namespace)
            ):
                unmatched.pop(index)
                matched = True
                break
        if not matched:
            return False
    return not unmatched


def _admit_mode(ir, metadata, report):
    # Mode is an oracle self-report.  Metadata on the expression may repeat
    # it for audit purposes, but may not invent an oracle mode.
    report_mode_value = _report_value(
        report,
        ("mode", "reference_mode", "source_mode", "evaluation_mode"),
    )
    if report_mode_value is _MISSING:
        raise _Rejected("oracle_mode_required")
    report_mode = _normalise_mode(report_mode_value)
    if report_mode not in ("synthetic_reference", "original_source"):
        raise _Rejected("oracle_mode_required")
    mode_value = _merge_explicit(
        metadata,
        report,
        ("mode", "reference_mode", "source_mode", "evaluation_mode"),
        "oracle_mode_mismatch",
    )
    mode = _normalise_mode(mode_value)
    if mode not in ("synthetic_reference", "original_source"):
        raise _Rejected("oracle_mode_required")
    if mode == "synthetic_reference":
        return mode, None, None, None

    report_family = _report_value(report, ("source_family", "family"))
    report_hash = _report_value(
        report,
        ("graph_hash", "graph_sha256", "source_graph_hash", "source_graph_sha256"),
    )
    if report_family is _MISSING:
        raise _Rejected("source_family_required")
    if report_hash is _MISSING:
        raise _Rejected("graph_hash_required")
    if _metadata_value(metadata, ("source_family", "family")) is _MISSING:
        raise _Rejected("expression_source_family_required")
    if _metadata_value(metadata, ("graph_hash", "graph_sha256", "source_graph_hash", "source_graph_sha256")) is _MISSING:
        raise _Rejected("expression_graph_hash_required")
    family = _merge_explicit(
        metadata,
        report,
        ("source_family", "family"),
        "source_family_mismatch",
    )
    graph_hash = _merge_explicit(
        metadata,
        report,
        ("graph_hash", "graph_sha256", "source_graph_hash", "source_graph_sha256"),
        "graph_hash_mismatch",
    )
    if not _valid_family(family):
        raise _Rejected("source_family_required")
    if not _valid_graph_hash(graph_hash):
        raise _Rejected("graph_hash_required")
    certificate = _certificate_from(metadata, report)
    if certificate is _MISSING:
        raise _Rejected("admission_certificate_required")
    cert_family = _first(certificate, ("source_family", "family"), _MISSING)
    if cert_family is _MISSING or cert_family != family:
        raise _Rejected("certificate_family_mismatch")
    cert_hash = _first(
        certificate,
        ("graph_hash", "graph_sha256", "source_graph_hash", "source_graph_sha256"),
        _MISSING,
    )
    if cert_hash is _MISSING or cert_hash != graph_hash:
        raise _Rejected("certificate_graph_hash_mismatch")
    tokens = _fresh_tokens(ir)
    if not _binder_identity_matches(certificate, tokens):
        raise _Rejected("binder_identity_certificate_required")
    integral_records = _integral_receipts(ir)
    if not _binder_receipts_match(certificate, integral_records):
        raise _Rejected("binder_receipt_mismatch")
    free_records = _free_variable_records(ir)
    if not _namespace_receipts_match(certificate, free_records):
        raise _Rejected("free_variable_namespace_mismatch")
    _check_source_rows(ir, family)
    return mode, family, graph_hash, certificate


class _FreeEnvironment:
    """Explicit free-variable environment with source-aware lookup first."""

    _META_KEYS = {
        "values",
        "bindings",
        "by_name",
        "by_source",
        "source_aware",
        "lexical_certificate",
        "__lexical_certificate__",
        "capture_certificate",
        "certificate",
    }

    def __init__(self, raw, *, allow_name):
        if raw is None:
            raw = {}
        if not isinstance(raw, Mapping):
            raise _Rejected("free_variables_not_mapping")
        self.raw = raw
        self.allow_name = allow_name
        self.maps = []
        self.source_maps = []
        self.certificates = []
        self._load(raw)

    def _load(self, raw):
        for key in (
            "lexical_certificate",
            "__lexical_certificate__",
            "capture_certificate",
            "certificate",
        ):
            if key in raw:
                self.certificates.append(raw[key])
        for key in ("values", "bindings", "by_name"):
            value = raw.get(key, _MISSING)
            if isinstance(value, Mapping):
                self.maps.append(value)
        self.maps.append(raw)
        for key in ("by_source", "source_aware"):
            value = raw.get(key, _MISSING)
            if isinstance(value, Mapping):
                self.source_maps.append(value)
        # A mapping under ``bindings`` may itself use source-aware tuple keys.
        for value in self.maps:
            for key, item in value.items():
                if isinstance(key, tuple) and len(key) == 2:
                    self.source_maps.append({key: item})

    def _source_lookup(self, name, source_ref):
        for mapping in self.source_maps:
            # Tuple keys are accepted in either human-readable order.  An
            # integer by itself is never interpreted as a graph node index.
            for key, value in mapping.items():
                if isinstance(key, tuple) and len(key) == 2:
                    left, right = key
                    if (
                        isinstance(left, str)
                        and left == name
                        and _same_reference(right, source_ref)
                    ) or (
                        isinstance(right, str)
                        and right == name
                        and _same_reference(left, source_ref)
                    ):
                        return value
                if isinstance(key, Mapping):
                    key_name = key.get("name", key.get("variable"))
                    key_ref = key.get("source_ref", key.get("source_node_reference"))
                    if key_name == name and _same_reference(key_ref, source_ref):
                        return value
                if _same_reference(key, source_ref) and isinstance(value, Mapping):
                    if name in value:
                        return value[name]
        return _MISSING

    def _lexical_certifies(self, name):
        for certificate in self.certificates:
            if certificate is True:
                return True
            if isinstance(certificate, Mapping):
                if certificate.get("name") == name:
                    return True
                for key in ("names", "variables", "bound_variables"):
                    values = certificate.get(key)
                    if isinstance(values, (list, tuple, set, frozenset)) and name in values:
                        return True
                if name in certificate and certificate[name] is not False:
                    return True
                if certificate.get("capture_safe") is True and certificate.get("all_names") is True:
                    return True
            elif isinstance(certificate, (list, tuple, set, frozenset)) and name in certificate:
                return True
        return False

    def resolve(self, name, source_ref):
        if not isinstance(name, str) or not name:
            raise _Rejected("invalid_bound_name")
        if source_ref is not _MISSING and source_ref is not None:
            value = self._source_lookup(name, source_ref)
            if value is not _MISSING:
                return _ensure_scalar(value, "invalid_bound_value")
        # Name-only lookup is intentionally certificate-gated.  This prevents
        # a graph handle such as integer 7 from becoming a variable value.
        if self.allow_name:
            for mapping in self.maps:
                if name in mapping and self._lexical_certifies(name):
                    return _ensure_scalar(mapping[name], "invalid_bound_value")
        raise _Rejected("unbound_variable")


class _Runtime:
    def __init__(self, ir, oracle, free_variables, backend, mode):
        self.ir = ir
        self.oracle = oracle
        self.free = _FreeEnvironment(free_variables, allow_name=(mode == "synthetic_reference"))
        self.backend = backend
        self.active = set()
        self.count = 0

    def _backend_member(self, names):
        member = _member(self.backend, names)
        if member is _MISSING:
            return _MISSING
        return member

    def _rational(self, numerator, denominator):
        if type(numerator) is not int or type(denominator) is not int:
            raise _Rejected("invalid_exact_rational")
        if denominator == 0:
            raise _Rejected("zero_exact_rational_denominator")
        method = self._backend_member(("rational",))
        if not callable(method):
            raise _Rejected("backend_rational_required")
        return _ensure_scalar(method(numerator, denominator), "backend_non_scalar")

    def _unary(self, name, value):
        method = self._backend_member(("unary",))
        if not callable(method):
            raise _Rejected("backend_unary_required")
        return _ensure_scalar(method(name, value), "backend_non_scalar")

    def _binary(self, names, left, right, fallback_code):
        method = self._backend_member(names)
        if callable(method):
            return _ensure_scalar(method(left, right), "backend_non_scalar")
        try:
            return _ensure_scalar(getattr(left, fallback_code)(right), "arithmetic_non_scalar")
        except (AttributeError, TypeError, ValueError):
            raise _Rejected("backend_arithmetic_required")

    def _add(self, left, right):
        method = self._backend_member(("add", "plus"))
        if callable(method):
            return _ensure_scalar(method(left, right), "backend_non_scalar")
        try:
            return _ensure_scalar(left + right, "arithmetic_non_scalar")
        except (TypeError, ValueError):
            raise _Rejected("backend_arithmetic_required")

    def _multiply(self, left, right):
        method = self._backend_member(("multiply", "mul", "times"))
        if callable(method):
            return _ensure_scalar(method(left, right), "backend_non_scalar")
        try:
            return _ensure_scalar(left * right, "arithmetic_non_scalar")
        except (TypeError, ValueError):
            raise _Rejected("backend_arithmetic_required")

    def _negative(self, value):
        method = self._backend_member(("negative", "neg"))
        if callable(method):
            return _ensure_scalar(method(value), "backend_non_scalar")
        try:
            return _ensure_scalar(-value, "arithmetic_non_scalar")
        except (TypeError, ValueError):
            raise _Rejected("backend_arithmetic_required")

    def _divide(self, left, right):
        method = self._backend_member(("divide", "div", "quotient"))
        if callable(method):
            return _ensure_scalar(method(left, right), "backend_non_scalar")
        try:
            return _ensure_scalar(left / right, "arithmetic_non_scalar")
        except (TypeError, ValueError, ZeroDivisionError):
            raise _Rejected("backend_arithmetic_required")

    def _strict_positive(self, value):
        method = self._backend_member(("is_positive", "positive"))
        if callable(method):
            result = method(value)
            if type(result) is not bool:
                raise _Rejected("backend_positive_check_not_bool")
            return result
        try:
            result = value > 0
        except (TypeError, ValueError):
            raise _Rejected("backend_positive_check_required")
        if type(result) is not bool:
            raise _Rejected("backend_positive_check_not_bool")
        return result

    def _oracle_method(self, names, code):
        method = _member(self.oracle, names)
        if not callable(method):
            raise _Rejected(code)
        return method

    def _source_record(self, node):
        source_row = node.get("source_row", _MISSING)
        if source_row is _MISSING:
            source_row = {}
        if not isinstance(source_row, Mapping):
            raise _Rejected("missing_source_row")
        record = _clone_json(source_row)
        if not isinstance(record, dict):
            raise _Rejected("invalid_source_row")
        # Preserve direct row fields for simple oracles, while also retaining
        # the exact v4 row and source/provider provenance as opaque data.
        record["op"] = node.get("op")
        if "source_ref" in node:
            record["source_ref"] = _clone_json(node.get("source_ref"))
        record["source_row"] = _clone_json(source_row)
        for name in (
            "provider",
            "provider_refs",
            "provenance",
            "quantity",
            "Z_derivative",
            "source_family",
        ):
            if name in node:
                record[name] = _clone_json(node[name])
        record["opaque_source_record"] = True
        return record

    def _eval(self, node, bindings, depth=0):
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
                return _ensure_scalar(bindings[token], "invalid_bound_value")

            op = node.get("op")
            if op not in _ALLOWED_OPS:
                if op == "function_substitution":
                    raise _Rejected("function_substitution_not_runtime")
                raise _Rejected("unknown_op")
            if op == "exact_rational":
                numerator = node.get("numerator", _MISSING)
                denominator = node.get("denominator", _MISSING)
                if numerator is _MISSING:
                    numerator = node.get("value", _MISSING)
                if denominator is _MISSING:
                    denominator = 1
                return self._rational(numerator, denominator)
            if op == "bound_variable":
                name = node.get("name", node.get("variable", _MISSING))
                source_ref = node.get(
                    "source_node_reference",
                    node.get("source_ref", _MISSING),
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
                argument = self._eval(node.get("argument"), bindings, depth + 1)
                return self._unary(name, argument)
            if op == "definite_integral":
                return self._eval_integral(node, bindings, depth)
            if op == "current_original_source_parameter":
                record = self._source_record(node)
                method = self._oracle_method(("parameter",), "oracle_parameter_required")
                return _ensure_scalar(method(record), "oracle_non_scalar_parameter")
            if op == "current_original_unique_analytic_core_component":
                arguments = node.get("arguments")
                if not isinstance(arguments, Mapping):
                    raise _Rejected("malformed_component_arguments")
                rho = self._eval(arguments.get("rho"), bindings, depth + 1)
                z_value = self._eval(arguments.get("Z_variable"), bindings, depth + 1)
                record = self._source_record(node)
                method = self._oracle_method(("component",), "oracle_component_required")
                return _ensure_scalar(
                    method(record, rho, z_value),
                    "oracle_non_scalar_component",
                )
            raise _Rejected("unknown_op")
        finally:
            self.active.remove(identity)

    def _eval_integral(self, node, bindings, depth):
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

        def sample(value):
            _ensure_scalar(value, "oracle_non_scalar_integral_sample")
            inner = dict(captured)
            if token in inner:
                raise _Rejected("binder_shadow")
            inner[token] = value
            return self._eval(integrand, inner, depth + 1)

        method = self._oracle_method(("integrate",), "oracle_integrate_required")
        return _ensure_scalar(method(sample, lower, upper), "oracle_non_scalar_integral")

    def run(self):
        return self._eval(self.ir, {}, 0)


def _result(status, *, value=None, mode=None, error=None, source_family=None, graph_hash=None):
    provenance = {
        "mode": mode,
        "caller_provenance": "not independently verified",
        "caller_provenance_independently_verified": False,
        "numeric_errors": "unverified",
        "numeric_errors_verified": False,
        "certification": False,
    }
    if source_family is not None:
        provenance["source_family"] = source_family
    if graph_hash is not None:
        provenance["graph_hash"] = graph_hash
    return {
        "status": status,
        "ok": status == "evaluated",
        "value": value if status == "evaluated" else None,
        "mode": mode,
        "provenance": provenance,
        "certification": False,
        "certified": False,
        "numeric_errors": "unverified",
        "numeric_errors_verified": False,
        "error": error,
        "limitations": [
            "Caller oracle provenance is not independently verified.",
            "Numeric error bounds are unverified; this result is not certified.",
            "Caller must not supply a bound, cap, or interval center as a source value; scalar provenance is not independently verified.",
            "function_substitution must already be lowered by the v4 compiler.",
        ],
    }


def evaluate(expression, *, oracle, free_variables, backend):
    """Evaluate one v4 IR expression through an explicit caller oracle.

    ``oracle`` must report ``synthetic_reference`` or ``original_source`` via
    ``mode``/``reference_mode`` (or a documented report method).  It must
    expose ``parameter(row)``, ``component(row, rho, Z)``, and
    ``integrate(callable, lower, upper)`` whenever the expression uses those
    source records.  ``backend`` must expose ``rational(n, d)`` and
    ``unary(name, value)``; ordinary arithmetic is used only when the returned
    scalar supports it directly or the backend supplies an arithmetic method.

    The return value is always a structured result.  Rejected inputs have
    ``status == 'rejected'`` and ``value is None``; they never receive an
    implicit zero or a saved archive value.
    """

    mode = None
    family = None
    graph_hash = None
    try:
        if oracle is None:
            raise _Rejected("oracle_required")
        if backend is None:
            raise _Rejected("backend_required")
        ir, metadata = _unwrap_expression(expression)
        report = _oracle_report(oracle)
        mode_value = _merge_explicit(
            metadata,
            report,
            ("mode", "reference_mode", "source_mode", "evaluation_mode"),
            "oracle_mode_mismatch",
        )
        mode = _normalise_mode(mode_value)
        mode, family, graph_hash, _certificate = _admit_mode(ir, metadata, report)
        value = _Runtime(ir, oracle, free_variables, backend, mode).run()
        return _result(
            "evaluated",
            value=value,
            mode=mode,
            source_family=family,
            graph_hash=graph_hash,
        )
    except _Rejected as exc:
        detail = exc.detail
        error = {"code": exc.code}
        if detail is not None:
            error["detail"] = detail
        return _result(
            "rejected",
            mode=mode,
            error=error,
            source_family=family,
            graph_hash=graph_hash,
        )
    except Exception as exc:
        # Oracle/backend implementations are caller code.  Keep their failure
        # fail-closed and do not turn an exception into a numeric value.
        return _result(
            "rejected",
            mode=mode,
            error={"code": "runtime_error", "type": type(exc).__name__},
            source_family=family,
            graph_hash=graph_hash,
        )


__all__ = ["evaluate"]
