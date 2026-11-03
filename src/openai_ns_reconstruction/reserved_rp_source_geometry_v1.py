"""Fail-closed source-graph geometry bridge for the reserved ``Rp`` path.

This module consumes already parsed archive reports and caller-supplied digest
strings.  It does not open files, hash raw bytes, decompress gzip, import a
producer, evaluate a source function, or compute numerical radius values.  The
outer runner is responsible for binding each raw archive's bytes to its
SHA-256 before calling this API; the bridge only checks the supplied digest
against the fixed digest pin and then inspects the supplied
``exact_graph_nodes`` rows.

The bridge has a deliberately narrow claim.  It proves the exact common graph
prefix and the symbolic scale recipe

``logRc_phys = original_Rc_offset + logminus``

``logRp_phys = logRc_phys + Tw - 2``

``Rp_phys / Rc_phys = exp(Tw - 2)``.

It also checks the source support endpoints and an explicit saved witness for
``Tw - 2 - log(2) > 0``.  This gives a *source-recipe support* statement.  It
does not authenticate a live owner, a complete ordinary-Z C1 matrix-row
identity, or the physical m15 function.  Those facts remain false in the
returned report even when every graph check passes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Final


SCHEMA: Final[str] = "reserved_rp_source_geometry_v1"

COMMON_PREFIX_LENGTH: Final[int] = 24952

RESERVED_ARCHIVE_RAW_SHA256: Final[str] = (
    "98fd45a48d72656ec91f19a5564577baed0dbb43fe85ce8c842336856dc17a17"
)
QUIET_C3_ARCHIVE_RAW_SHA256: Final[str] = (
    "cb2a69ca971965a6e487bd0fa818204216747094f37819a2cc3e42a3af79911a"
)
NATIVE_FRAME_ARCHIVE_RAW_SHA256: Final[str] = (
    "859b59d422104d184092bcf60f8afa95267335ce87e8b41f2c271c18130e9a28"
)

EXPECTED_PARAMETER_NODES: Final[dict[str, int]] = {
    "logP": 5,
    "logC": 6,
    "T": 7,
    "hbB": 8,
    "hbS": 9,
    "sc": 10,
    "mu": 17,
}

EXPECTED_RESERVED_LOG_RC_NODE: Final[int] = 24956
EXPECTED_RESERVED_LOG_2RC_NODE: Final[int] = 24957
EXPECTED_QUIET_OFFSET_NODE: Final[int] = 4893
EXPECTED_QUIET_RC_NODE: Final[int] = 25494
EXPECTED_QUIET_R2_NODE: Final[int] = 25495
EXPECTED_QUIET_RP_NODE: Final[int] = 25499
EXPECTED_QUIET_TW_NODE: Final[int] = 25007
EXPECTED_QUIET_LENGTH_NODE: Final[int] = 25010
EXPECTED_NATIVE_LOG_MINUS_NODE: Final[int] = 26684
EXPECTED_NATIVE_LOG_RC_NODE: Final[int] = 26685
EXPECTED_NATIVE_LOG_RP_NODE: Final[int] = 26687
EXPECTED_NATIVE_TW_NODE: Final[int] = 26013

EXPECTED_SUPPORT_ENDPOINTS: Final[tuple[tuple[Fraction, Fraction], ...]] = (
    (Fraction(119, 400), Fraction(121, 400)),
    (Fraction(129, 400), Fraction(131, 400)),
    (Fraction(139, 400), Fraction(141, 400)),
    (Fraction(149, 400), Fraction(151, 400)),
    (Fraction(159, 400), Fraction(161, 400)),
)
EXPECTED_ANGULAR_INDICES: Final[tuple[int, int, int]] = (1, 2, 3)
EXPECTED_AXIAL_INDICES: Final[tuple[int, int]] = (0, 4)
EXPECTED_QUIET_POSITIVE_EXPRESSION: Final[str] = "Tw-2-log(2)>0"


class SourceGeometryError(ValueError):
    """Raised when a caller-supplied source-graph precondition is absent."""


@dataclass(frozen=True, slots=True)
class AuthenticatedGraph:
    """A parsed archive report plus a caller-supplied raw digest.

    ``archive`` is expected to be the parsed JSON object, not compressed bytes.
    The outer runner must bind the digest to the raw archive bytes before
    constructing this value.  This consumer intentionally does not reopen,
    hash, or decompress an archive as part of a geometry check.
    """

    role: str
    archive: Mapping[str, Any]
    raw_sha256: str


@dataclass(frozen=True, slots=True)
class QuietLengthEvidence:
    """Saved source evidence that the quiet positive length was checked."""

    checked: bool
    expression: str = EXPECTED_QUIET_POSITIVE_EXPRESSION
    source: str = "caller-supplied quiet-length receipt"


@dataclass(frozen=True, slots=True)
class ReservedRpSourceGeometryReport:
    """JSON-safe result of the exact source-graph geometry audit."""

    schema: str
    raw_sha256: Mapping[str, str]
    checks: Mapping[str, bool]
    graphs_compared: tuple[str, ...]
    native_prefix_compared: bool
    common_prefix_length: int
    common_prefix_exact: bool
    node_ids: Mapping[str, int]
    expressions: Mapping[str, str]
    support_endpoints: tuple[tuple[str, str], ...]
    angular_support_endpoints: tuple[tuple[str, str], ...]
    all_support_upper_log2: str
    angular_support_upper_log2: str
    all_support_margin_to_two_rc_log2: str
    angular_support_margin_to_two_rc_log2: str
    quiet_positive_length_expression: str
    source_recipe_full_support: bool
    live_owner_proved: bool
    full_c1_matrix_proved: bool
    physical_m15_proved: bool
    limitations: tuple[str, ...]

    def report(self) -> dict[str, Any]:
        """Return a serialization-safe report without graph objects."""

        return {
            "schema": self.schema,
            "raw_sha256": dict(self.raw_sha256),
            "checks": dict(self.checks),
            "graphs_compared": list(self.graphs_compared),
            "native_prefix_compared": self.native_prefix_compared,
            "common_prefix": {
                "length": self.common_prefix_length,
                "exact": self.common_prefix_exact,
            },
            "node_ids": dict(self.node_ids),
            "expressions": dict(self.expressions),
            "support": {
                "all_endpoints_over_log2": [list(x) for x in self.support_endpoints],
                "angular_endpoints_over_log2": [
                    list(x) for x in self.angular_support_endpoints
                ],
                "all_support_upper_log2": self.all_support_upper_log2,
                "angular_support_upper_log2": self.angular_support_upper_log2,
                "all_support_margin_to_two_rc_log2": self.all_support_margin_to_two_rc_log2,
                "angular_support_margin_to_two_rc_log2": self.angular_support_margin_to_two_rc_log2,
                "quiet_positive_length": self.quiet_positive_length_expression,
                "source_recipe_full_support": self.source_recipe_full_support,
            },
            "owner_and_matrix_scope": {
                "live_owner_proved": self.live_owner_proved,
                "full_c1_matrix_proved": self.full_c1_matrix_proved,
                "physical_m15_proved": self.physical_m15_proved,
            },
            "limitations": list(self.limitations),
        }


Expr = tuple[Any, ...]


def _fail(message: str) -> None:
    raise SourceGeometryError(message)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _fail(f"{label} must be a mapping")
    return value


def _nodes(archive: Mapping[str, Any], label: str) -> Sequence[Mapping[str, Any]]:
    rows = archive.get("exact_graph_nodes")
    if not isinstance(rows, (list, tuple)) or not rows:
        _fail(f"{label}.exact_graph_nodes must be a nonempty list")
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            _fail(f"{label}.exact_graph_nodes[{index}] must be a mapping")
    return rows


def _digest_from_mapping(
    value: Mapping[str, Any], role: str, raw_sha256: Mapping[str, str] | None
) -> tuple[Mapping[str, Any], str]:
    nested = value.get("archive")
    archive = _mapping(nested if nested is not None else value, f"{role} archive")
    supplied = value.get("raw_sha256")
    if supplied is None:
        supplied = value.get("archive_raw_sha256")
    if supplied is None and raw_sha256 is not None:
        aliases = {
            "reserved": ("reserved", "reserved_linear_operator", "n1_reserved_linear_operator"),
            "quiet_c3": ("quiet_c3", "quiet", "c3_power", "C3_power_to_Rp"),
            "native_frame": ("native_frame", "c3_native_frame", "C3_Rp_native_frame"),
        }[role]
        for key in aliases:
            if key in raw_sha256:
                supplied = raw_sha256[key]
                break
    if not isinstance(supplied, str):
        _fail(f"{role} raw_sha256 is required; compressed data are not read here")
    return archive, supplied


def _coerce_graph(
    value: AuthenticatedGraph | Mapping[str, Any],
    role: str,
    raw_sha256: Mapping[str, str] | None,
) -> tuple[Mapping[str, Any], str]:
    if isinstance(value, AuthenticatedGraph):
        if value.role != role:
            _fail(f"caller-supplied graph role {value.role!r} is not {role!r}")
        return _mapping(value.archive, f"{role} archive"), value.raw_sha256
    if isinstance(value, Mapping):
        return _digest_from_mapping(value, role, raw_sha256)
    _fail(f"{role} must be AuthenticatedGraph or parsed archive mapping")


def _check_digest(role: str, digest: str, expected: str) -> str:
    if not isinstance(digest, str) or len(digest) != 64:
        _fail(f"{role} raw_sha256 must be a 64-character hexadecimal digest")
    lowered = digest.lower()
    if any(character not in "0123456789abcdef" for character in lowered):
        _fail(f"{role} raw_sha256 is not hexadecimal")
    if lowered != expected:
        _fail(f"{role} caller-supplied raw_sha256 does not match the pinned archive digest")
    return lowered


def _require_key(mapping: Mapping[str, Any], key: str, expected: Any, label: str) -> None:
    if mapping.get(key) != expected:
        _fail(f"{label}.{key} does not match the pinned source binding")


def _row(nodes: Sequence[Mapping[str, Any]], node: int, label: str) -> Mapping[str, Any]:
    if type(node) is not int or node < 0 or node >= len(nodes):
        _fail(f"{label} node {node!r} is outside exact_graph_nodes")
    return nodes[node]


def _require_row_fields(
    nodes: Sequence[Mapping[str, Any]], node: int, expected: Mapping[str, Any], label: str
) -> None:
    actual = _row(nodes, node, label)
    for key, value in expected.items():
        if actual.get(key) != value:
            _fail(f"{label} node {node} field {key!r} is not the expected graph edge")


def _fraction_literal(value: Any, label: str) -> Fraction:
    if not isinstance(value, str):
        _fail(f"{label} must be an exact rational string")
    try:
        return Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise SourceGeometryError(f"{label} is not an exact rational string") from exc


def _canonical_sum(children: Sequence[Expr]) -> Expr:
    flattened: list[Expr] = []
    for child in children:
        if child[0] == "sum":
            flattened.extend(child[1])
        else:
            flattened.append(child)
    return ("sum", tuple(sorted(flattened, key=repr)))


def _canonical_product(children: Sequence[Expr]) -> Expr:
    flattened: list[Expr] = []
    for child in children:
        if child[0] == "product":
            flattened.extend(child[1])
        else:
            flattened.append(child)
    return ("product", tuple(sorted(flattened, key=repr)))


def _graph_expr(
    nodes: Sequence[Mapping[str, Any]],
    node: int,
    stack: tuple[int, ...] = (),
    parameter_nodes: Mapping[int, str] | None = None,
) -> Expr:
    """Decode the small exact arithmetic language used by the target roots.

    Unknown operations fail closed.  In particular, this is not an evaluator:
    it returns an immutable expression tree and never calls ``exp``, ``log`` or
    a numerical/integral helper.
    """

    if node in stack:
        _fail(f"cycle in exact graph expression at node {node}")
    row = _row(nodes, node, "graph")
    if parameter_nodes is not None and node in parameter_nodes:
        return ("parameter", parameter_nodes[node])
    operation = row.get("operation")
    next_stack = stack + (node,)
    if operation == "exact_rational":
        numerator = row.get("numerator")
        denominator = row.get("denominator")
        if type(numerator) is not int or type(denominator) is not int or denominator == 0:
            _fail(f"node {node} has an invalid exact_rational")
        value = Fraction(numerator, denominator)
        return ("rational", value.numerator, value.denominator)
    if operation == "current_original_source_parameter":
        name = row.get("name")
        if not isinstance(name, str) or not name:
            _fail(f"node {node} has an invalid source parameter leaf")
        return ("parameter", name)
    if operation == "bound_variable":
        name = row.get("name")
        if not isinstance(name, str) or not name:
            _fail(f"node {node} has an invalid bound-variable leaf")
        return ("bound", name)
    if operation == "analytic_unary":
        name = row.get("name")
        argument = row.get("argument")
        if name not in ("log", "exp") or type(argument) is not int:
            _fail(f"node {node} has an unsupported analytic_unary")
        return (name, _graph_expr(nodes, argument, next_stack, parameter_nodes))
    if operation == "negative":
        argument = row.get("argument")
        if type(argument) is not int:
            _fail(f"node {node} has an invalid negative edge")
        child = _graph_expr(nodes, argument, next_stack, parameter_nodes)
        if child[0] == "rational":
            return _rational(-child[1], child[2])
        return ("negative", child)
    if operation in ("sum", "product"):
        arguments = row.get("arguments")
        if not isinstance(arguments, (list, tuple)) or not arguments:
            _fail(f"node {node} has an invalid {operation} edge list")
        if any(type(child) is not int for child in arguments):
            _fail(f"node {node} has a non-integer {operation} edge")
        children = [
            _graph_expr(nodes, child, next_stack, parameter_nodes) for child in arguments
        ]
        return _canonical_sum(children) if operation == "sum" else _canonical_product(children)
    _fail(f"node {node} uses unsupported exact graph operation {operation!r}")


def _expr_text(expr: Expr) -> str:
    tag = expr[0]
    if tag == "rational":
        numerator, denominator = expr[1], expr[2]
        if denominator == 1:
            return str(numerator)
        return f"{numerator}/{denominator}"
    if tag == "parameter" or tag == "bound":
        return str(expr[1])
    if tag == "log" or tag == "exp":
        return f"{tag}({_expr_text(expr[1])})"
    if tag == "negative":
        return f"-({_expr_text(expr[1])})"
    if tag == "sum":
        return " + ".join(_expr_text(child) for child in expr[1])
    if tag == "product":
        return "*".join(_expr_text(child) for child in expr[1])
    _fail(f"cannot format unknown expression tag {tag!r}")


def _rational(numerator: int, denominator: int = 1) -> Expr:
    value = Fraction(numerator, denominator)
    return ("rational", value.numerator, value.denominator)


def _negative(expr: Expr) -> Expr:
    return ("negative", expr)


def _sum(*children: Expr) -> Expr:
    return _canonical_sum(children)


def _product(*children: Expr) -> Expr:
    return _canonical_product(children)


def _log(expr: Expr) -> Expr:
    return ("log", expr)


def _exp(expr: Expr) -> Expr:
    return ("exp", expr)


def _parameter(name: str) -> Expr:
    return ("parameter", name)


def _expected_offset() -> Expr:
    log_p = _parameter("logP")
    log_c = _parameter("logC")
    hb_b = _parameter("hbB")
    sc = _parameter("sc")
    return _sum(
        _rational(1014),
        _product(_rational(10), log_c),
        _product(_rational(14), log_p),
        _product(_rational(-1, 2), hb_b, sc),
        _exp(_rational(40)),
        _log(_rational(25)),
        _log(_rational(11, 10)),
    )


def _expected_logminus() -> Expr:
    return _sum(
        _log(_rational(4)),
        _negative(_product(_rational(4), _parameter("logP"))),
        _rational(-1000),
        _product(_rational(1, 2), _parameter("hbB"), _parameter("sc")),
    )


def _expected_tw() -> Expr:
    return _product(_rational(-60), _log(_parameter("mu")))


def _parse_support_certificate(
    reserved_archive: Mapping[str, Any],
    explicit: Mapping[str, Any] | None,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    Mapping[str, bool],
]:
    certificate = explicit
    if certificate is None:
        certificate = _mapping(reserved_archive.get("support_certificate"), "support_certificate")
    expected = [
        [f"{left.numerator}/{left.denominator}", f"{right.numerator}/{right.denominator}"]
        for left, right in EXPECTED_SUPPORT_ENDPOINTS
    ]
    if certificate.get("new_log_supports_over_log2") != expected:
        _fail("reserved support endpoints are not the pinned source endpoints")
    if certificate.get("exact_reserved_gap_over_log2") != ["11/40", "17/40"]:
        _fail("reserved gap endpoints are not the pinned source endpoints")
    if certificate.get("axial_indices") != list(EXPECTED_AXIAL_INDICES):
        _fail("reserved axial support indices changed")
    if certificate.get("angular_indices") != list(EXPECTED_ANGULAR_INDICES):
        _fail("reserved angular support indices changed")
    if certificate.get("autonomous_support_choice") is not True:
        _fail("reserved support recipe is not marked as the saved autonomous choice")
    if certificate.get("old_corrections_and_all_local_derivatives_zero_on_new_support") is not True:
        _fail("old correction flatness on the new support was not checked")
    if certificate.get("nonlocal_histories_pressure_and_density_not_zeroed") is not True:
        _fail("nonlocal histories/pressure/density retention was not checked")
    return (
        EXPECTED_SUPPORT_ENDPOINTS,
        tuple(EXPECTED_SUPPORT_ENDPOINTS[index] for index in EXPECTED_ANGULAR_INDICES),
        {
            "support_certificate_exact": True,
            "support_indices_exact": True,
            "support_flatness_checked": True,
            "nonlocal_histories_retained": True,
        },
    )


def _coerce_quiet_length_evidence(
    evidence: QuietLengthEvidence | Mapping[str, Any] | None,
    source_bindings: Mapping[str, Any] | None,
) -> QuietLengthEvidence:
    candidate: Any = evidence
    if candidate is None and source_bindings is not None:
        candidate = source_bindings.get("quiet_positive_length")
        if candidate is None:
            checked = source_bindings.get("quiet_positive_length_checked")
            expression = source_bindings.get("quiet_positive_length_expression")
            if checked is not None or expression is not None:
                candidate = {
                    "checked": checked,
                    "expression": expression,
                    "source": source_bindings.get("quiet_positive_length_source"),
                }
    if isinstance(candidate, QuietLengthEvidence):
        result = candidate
    elif isinstance(candidate, Mapping):
        checked = candidate.get("checked")
        if checked is None:
            checked = candidate.get("positive_checked")
        expression = candidate.get("expression")
        if expression is None:
            expression = candidate.get("quiet_length_expression")
        source = candidate.get("source") or candidate.get("receipt") or "caller-supplied evidence"
        if type(checked) is not bool or not isinstance(expression, str):
            _fail("quiet positive-length evidence must include checked and expression")
        result = QuietLengthEvidence(checked=checked, expression=expression, source=str(source))
    else:
        _fail("explicit quiet positive-length source evidence is required")
    if result.checked is not True:
        _fail("quiet positive length was not checked as positive")
    if result.expression != EXPECTED_QUIET_POSITIVE_EXPRESSION:
        _fail("quiet positive-length expression is not Tw-2-log(2)>0")
    return result


def _validate_source_bindings(
    source_bindings: Mapping[str, Any] | None,
    reserved_archive: Mapping[str, Any],
    quiet_archive: Mapping[str, Any],
) -> None:
    if source_bindings is None:
        return
    _mapping(source_bindings, "source_bindings")
    for key in ("parameters", "parameter_nodes"):
        if key in source_bindings and source_bindings[key] != EXPECTED_PARAMETER_NODES:
            _fail(f"source_bindings.{key} does not match the common source parameter binding")
    if "source_family" in source_bindings:
        family = source_bindings["source_family"]
        if family != reserved_archive.get("source_family") or family != quiet_archive.get("source_family"):
            _fail("source_bindings.source_family is not the exact common source family")
    for role in ("reserved", "quiet_c3", "native_frame"):
        role_binding = source_bindings.get(role)
        if role_binding is not None and not isinstance(role_binding, Mapping):
            _fail(f"source_bindings.{role} must be a mapping when supplied")


def _validate_archive_reports(
    reserved_archive: Mapping[str, Any],
    quiet_archive: Mapping[str, Any],
    native_archive: Mapping[str, Any],
) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]:
    reserved_binding = _mapping(reserved_archive.get("actual_background_binding"), "reserved.actual_background_binding")
    reserved_functions = _mapping(
        reserved_archive.get("exact_normalized_basis_and_matrix_functions"),
        "reserved.exact_normalized_basis_and_matrix_functions",
    )
    quiet_functions = _mapping(
        quiet_archive.get("actual_C3_repaired_exit_power_and_pulse_functions"),
        "quiet.actual_C3_repaired_exit_power_and_pulse_functions",
    )
    native_functions = _mapping(
        native_archive.get("actual_Rp_native_frame_functions"),
        "native.actual_Rp_native_frame_functions",
    )
    for label, archive in (
        ("reserved", reserved_archive),
        ("quiet_c3", quiet_archive),
        ("native_frame", native_archive),
    ):
        family = archive.get("source_family")
        if not isinstance(family, Mapping):
            _fail(f"{label}.source_family is required for source identity")
    if not (
        reserved_archive["source_family"]
        == quiet_archive["source_family"]
        == native_archive["source_family"]
    ):
        _fail("reserved, quiet C3 and native frame source_family records differ")
    if reserved_archive.get("accepted_background_graph_prefix_length") != COMMON_PREFIX_LENGTH:
        _fail("reserved accepted graph prefix length is not 24952")
    if quiet_archive.get("accepted_original_C3_band_graph_prefix_length") != COMMON_PREFIX_LENGTH:
        _fail("quiet C3 accepted graph prefix length is not 24952")
    if native_archive.get("accepted_complete_Rc_graph_prefix_length") != 26011:
        _fail("native frame accepted graph prefix length is not 26011")
    if reserved_binding.get("parameters") != EXPECTED_PARAMETER_NODES:
        _fail("reserved parameter binding differs from the pinned source")
    if quiet_functions.get("parameters") != EXPECTED_PARAMETER_NODES:
        _fail("quiet C3 parameter binding differs from the pinned source")
    _require_key(reserved_binding, "logRc", EXPECTED_RESERVED_LOG_RC_NODE, "reserved.actual_background_binding")
    _require_key(reserved_binding, "log2Rc", EXPECTED_RESERVED_LOG_2RC_NODE, "reserved.actual_background_binding")
    for key, expected in {
        "original_Rc_offset": EXPECTED_QUIET_OFFSET_NODE,
        "Rc": EXPECTED_QUIET_RC_NODE,
        "R2": EXPECTED_QUIET_R2_NODE,
        "Rp": EXPECTED_QUIET_RP_NODE,
        "Tw": EXPECTED_QUIET_TW_NODE,
        "quiet_length": EXPECTED_QUIET_LENGTH_NODE,
    }.items():
        _require_key(quiet_functions, key, expected, "quiet C3 functions")
    _require_key(native_functions, "original_Rc_phase_offset", EXPECTED_QUIET_OFFSET_NODE, "native functions")
    _require_key(native_functions, "original_log_r_minus", EXPECTED_NATIVE_LOG_MINUS_NODE, "native functions")
    _require_key(native_functions, "original_logRc", EXPECTED_NATIVE_LOG_RC_NODE, "native functions")
    _require_key(native_functions, "original_logRp", EXPECTED_NATIVE_LOG_RP_NODE, "native functions")
    _require_key(native_functions, "Tw", EXPECTED_NATIVE_TW_NODE, "native functions")
    if quiet_functions.get("coordinate_identity") != "s=log(R/(2Rc)); 0<=s<=Tw-2-log2; Rp=Rc*exp(Tw-2)":
        _fail("quiet C3 coordinate identity is not the pinned source recipe")
    return reserved_binding, reserved_functions, quiet_functions, native_functions, reserved_archive["source_family"]


def _check_exact_graph_geometry(
    reserved_archive: Mapping[str, Any],
    quiet_archive: Mapping[str, Any],
    native_archive: Mapping[str, Any],
    reserved_binding: Mapping[str, Any],
    reserved_functions: Mapping[str, Any],
    quiet_functions: Mapping[str, Any],
    native_functions: Mapping[str, Any],
) -> tuple[dict[str, bool], dict[str, str]]:
    reserved_nodes = _nodes(reserved_archive, "reserved")
    quiet_nodes = _nodes(quiet_archive, "quiet_c3")
    native_nodes = _nodes(native_archive, "native_frame")
    parameter_nodes = {node: name for name, node in EXPECTED_PARAMETER_NODES.items()}
    if len(reserved_nodes) <= EXPECTED_RESERVED_LOG_2RC_NODE or len(quiet_nodes) <= EXPECTED_QUIET_RP_NODE:
        _fail("reserved or quiet graph is shorter than the pinned target nodes")
    if len(native_nodes) <= EXPECTED_NATIVE_LOG_RP_NODE:
        _fail("native frame graph is shorter than the pinned target nodes")
    if list(reserved_nodes[:COMMON_PREFIX_LENGTH]) != list(quiet_nodes[:COMMON_PREFIX_LENGTH]):
        _fail("reserved and quiet exact_graph_nodes do not share the pinned 24952-row prefix")

    # The prefix proves the actual shared source node; these root rows prove the
    # post-prefix constructions rather than relying on names or interval overlap.
    _require_row_fields(
        reserved_nodes,
        EXPECTED_QUIET_OFFSET_NODE,
        {"operation": "sum", "arguments": [4892, 2084, 2091, 649, 3, 645, 1507]},
        "reserved",
    )
    _require_row_fields(
        quiet_nodes,
        EXPECTED_QUIET_OFFSET_NODE,
        {"operation": "sum", "arguments": [4892, 2084, 2091, 649, 3, 645, 1507]},
        "quiet_c3",
    )
    _require_row_fields(
        native_nodes,
        EXPECTED_QUIET_OFFSET_NODE,
        {"operation": "sum", "arguments": [4892, 2084, 2091, 649, 3, 645, 1507]},
        "native_frame",
    )
    _require_row_fields(
        reserved_nodes,
        24955,
        {"operation": "sum", "arguments": [24952, 15, 24953, 24954]},
        "reserved",
    )
    _require_row_fields(
        reserved_nodes,
        EXPECTED_RESERVED_LOG_RC_NODE,
        {"operation": "sum", "arguments": [4893, 24955]},
        "reserved",
    )
    _require_row_fields(
        reserved_nodes,
        EXPECTED_RESERVED_LOG_2RC_NODE,
        {"operation": "sum", "arguments": [3184, 24956]},
        "reserved",
    )
    _require_row_fields(
        reserved_nodes,
        24978,
        {"operation": "bound_variable", "name": "n1_reserved_log_R_over_Rc"},
        "reserved",
    )
    _require_row_fields(
        quiet_nodes,
        EXPECTED_QUIET_RC_NODE,
        {"operation": "analytic_unary", "name": "exp", "argument": 4893},
        "quiet_c3",
    )
    _require_row_fields(
        quiet_nodes,
        EXPECTED_QUIET_R2_NODE,
        {"operation": "product", "arguments": [48, 25494]},
        "quiet_c3",
    )
    _require_row_fields(
        quiet_nodes,
        EXPECTED_QUIET_RP_NODE,
        {"operation": "product", "arguments": [25495, 25498]},
        "quiet_c3",
    )
    _require_row_fields(
        quiet_nodes,
        EXPECTED_QUIET_TW_NODE,
        {"operation": "product", "arguments": [25005, 25006]},
        "quiet_c3",
    )
    _require_row_fields(
        quiet_nodes,
        EXPECTED_QUIET_LENGTH_NODE,
        {"operation": "sum", "arguments": [25009, 12628]},
        "quiet_c3",
    )
    _require_row_fields(
        native_nodes,
        EXPECTED_NATIVE_LOG_MINUS_NODE,
        {"operation": "sum", "arguments": [26681, 15, 26682, 26683]},
        "native_frame",
    )
    _require_row_fields(
        native_nodes,
        EXPECTED_NATIVE_LOG_RC_NODE,
        {"operation": "sum", "arguments": [4893, 26684]},
        "native_frame",
    )
    _require_row_fields(
        native_nodes,
        26686,
        {"operation": "sum", "arguments": [26685, 26014]},
        "native_frame",
    )
    _require_row_fields(
        native_nodes,
        EXPECTED_NATIVE_LOG_RP_NODE,
        {"operation": "sum", "arguments": [26686, 26013]},
        "native_frame",
    )

    offset_expected = _expected_offset()
    logminus_expected = _expected_logminus()
    tw_expected = _expected_tw()
    offset_expressions = {
        "reserved": _graph_expr(reserved_nodes, EXPECTED_QUIET_OFFSET_NODE, parameter_nodes=parameter_nodes),
        "quiet_c3": _graph_expr(quiet_nodes, EXPECTED_QUIET_OFFSET_NODE, parameter_nodes=parameter_nodes),
        "native_frame": _graph_expr(native_nodes, EXPECTED_QUIET_OFFSET_NODE, parameter_nodes=parameter_nodes),
    }
    logminus_expressions = {
        "reserved": _graph_expr(reserved_nodes, 24955, parameter_nodes=parameter_nodes),
        "native_frame": _graph_expr(native_nodes, EXPECTED_NATIVE_LOG_MINUS_NODE, parameter_nodes=parameter_nodes),
    }
    if any(expression != offset_expected for expression in offset_expressions.values()):
        _fail("node 4893 does not decode to the pinned original_Rc_offset recipe")
    if logminus_expressions["reserved"] != logminus_expected or logminus_expressions["native_frame"] != logminus_expected:
        _fail("reserved/native logminus expressions do not match")
    if logminus_expressions["reserved"] != logminus_expressions["native_frame"]:
        _fail("reserved and native logminus expressions differ")

    quiet_rc = _graph_expr(quiet_nodes, EXPECTED_QUIET_RC_NODE, parameter_nodes=parameter_nodes)
    quiet_r2 = _graph_expr(quiet_nodes, EXPECTED_QUIET_R2_NODE, parameter_nodes=parameter_nodes)
    quiet_rp = _graph_expr(quiet_nodes, EXPECTED_QUIET_RP_NODE, parameter_nodes=parameter_nodes)
    quiet_tw = _graph_expr(quiet_nodes, EXPECTED_QUIET_TW_NODE, parameter_nodes=parameter_nodes)
    quiet_length = _graph_expr(quiet_nodes, EXPECTED_QUIET_LENGTH_NODE, parameter_nodes=parameter_nodes)
    native_tw = _graph_expr(native_nodes, EXPECTED_NATIVE_TW_NODE, parameter_nodes=parameter_nodes)
    native_logminus = _graph_expr(native_nodes, EXPECTED_NATIVE_LOG_MINUS_NODE, parameter_nodes=parameter_nodes)
    native_logrc = _graph_expr(native_nodes, EXPECTED_NATIVE_LOG_RC_NODE, parameter_nodes=parameter_nodes)
    native_logrp = _graph_expr(native_nodes, EXPECTED_NATIVE_LOG_RP_NODE, parameter_nodes=parameter_nodes)
    if quiet_rc != _exp(offset_expected):
        _fail("quiet Rc node 25494 is not exp(offset node 4893)")
    if quiet_tw != tw_expected or native_tw != tw_expected or quiet_tw != native_tw:
        _fail("quiet/native Tw expressions do not match the pinned recipe")
    quiet_length_expected = _sum(
        tw_expected,
        _rational(-2),
        _negative(_log(_rational(2))),
    )
    if quiet_length != quiet_length_expected:
        _fail("quiet_length node 25010 is not Tw-2-log(2)")
    if quiet_r2 != _product(_rational(2), _exp(offset_expected)):
        _fail("quiet R2 node 25495 is not 2*exp(offset)")
    if quiet_rp != _product(quiet_r2, _exp(quiet_length_expected)):
        _fail("quiet Rp node 25499 does not use the pinned quiet length")
    if native_logminus != logminus_expected:
        _fail("native node 26684 is not the pinned logminus expression")
    if native_logrc != _sum(offset_expected, logminus_expected):
        _fail("native node 26685 is not offset+logminus")
    if native_logrp != _sum(native_logrc, tw_expected, _rational(-2)):
        _fail("native node 26687 is not logRc+Tw-2")

    expressions = {
        "offset_4893": _expr_text(offset_expected),
        "logminus_24955": _expr_text(logminus_expected),
        "logRc_phys_24956": "offset_4893 + logminus_24955",
        "quiet_Rc_25494": "exp(offset_4893)",
        "quiet_Tw_25007": _expr_text(tw_expected),
        "quiet_length_25010": "Tw - 2 - log(2)",
        "quiet_Rp_25499": "2*exp(offset_4893)*exp(Tw-2-log(2))",
        "native_logminus_26684": _expr_text(logminus_expected),
        "native_logRc_phys_26685": "offset_4893 + logminus_24955",
        "native_logRp_phys_26687": "logRc_phys + Tw - 2",
        "Rp_over_Rc": "exp(Tw - 2)",
    }
    checks = {
        "reserved_quiet_prefix_24952_exact": True,
        "shared_offset_node_4893_exact": True,
        "reserved_logminus_24955_exact": True,
        "reserved_logRc_24956_exact": True,
        "reserved_log2Rc_24957_exact": True,
        "quiet_Rc_25494_exact": True,
        "quiet_Rp_25499_exact": True,
        "quiet_Tw_25007_exact": True,
        "quiet_length_25010_exact": True,
        "native_logminus_26684_exact": True,
        "native_logRc_26685_exact": True,
        "native_logRp_26687_exact": True,
        "quiet_native_Tw_expression_equal": True,
        "physical_log_radius_recipe_equal": True,
        "rp_over_rc_recipe_derived_without_numeric_overlap": True,
    }
    return checks, expressions


def audit_reserved_rp_source_geometry(
    *,
    reserved: AuthenticatedGraph | Mapping[str, Any],
    quiet_c3: AuthenticatedGraph | Mapping[str, Any],
    native_frame: AuthenticatedGraph | Mapping[str, Any],
    raw_sha256: Mapping[str, str] | None = None,
    source_bindings: Mapping[str, Any] | None = None,
    support_certificate: Mapping[str, Any] | None = None,
    quiet_length_evidence: QuietLengthEvidence | Mapping[str, Any] | None = None,
) -> ReservedRpSourceGeometryReport:
    """Audit caller-bound reserved/quiet/native source geometry.

    The three graph arguments may be :class:`AuthenticatedGraph` instances or
    parsed archive mappings containing ``exact_graph_nodes`` plus a raw digest
    (or a digest supplied through ``raw_sha256``).  The outer runner must bind
    each digest to the raw archive content before calling this API; this module
    does not perform that content-hash binding.  No archive path is accepted
    and no compressed bytes are read.  Missing evidence or any graph mismatch
    raises :class:`SourceGeometryError` before a report is returned.
    """

    reserved_archive, reserved_digest = _coerce_graph(reserved, "reserved", raw_sha256)
    quiet_archive, quiet_digest = _coerce_graph(quiet_c3, "quiet_c3", raw_sha256)
    native_archive, native_digest = _coerce_graph(native_frame, "native_frame", raw_sha256)
    digests = {
        "reserved": _check_digest("reserved", reserved_digest, RESERVED_ARCHIVE_RAW_SHA256),
        "quiet_c3": _check_digest("quiet_c3", quiet_digest, QUIET_C3_ARCHIVE_RAW_SHA256),
        "native_frame": _check_digest("native_frame", native_digest, NATIVE_FRAME_ARCHIVE_RAW_SHA256),
    }
    _validate_source_bindings(source_bindings, reserved_archive, quiet_archive)
    (
        reserved_binding,
        reserved_functions,
        quiet_functions,
        native_functions,
        _source_family,
    ) = _validate_archive_reports(reserved_archive, quiet_archive, native_archive)
    del reserved_binding, reserved_functions, quiet_functions, native_functions, _source_family
    checks, expressions = _check_exact_graph_geometry(
        reserved_archive,
        quiet_archive,
        native_archive,
        _mapping(reserved_archive["actual_background_binding"], "reserved.actual_background_binding"),
        _mapping(
            reserved_archive["exact_normalized_basis_and_matrix_functions"],
            "reserved.exact_normalized_basis_and_matrix_functions",
        ),
        _mapping(
            quiet_archive["actual_C3_repaired_exit_power_and_pulse_functions"],
            "quiet.actual_C3_repaired_exit_power_and_pulse_functions",
        ),
        _mapping(native_archive["actual_Rp_native_frame_functions"], "native.actual_Rp_native_frame_functions"),
    )
    all_support, angular_support, support_checks = _parse_support_certificate(
        reserved_archive, support_certificate
    )
    quiet_evidence = _coerce_quiet_length_evidence(quiet_length_evidence, source_bindings)
    checks.update(support_checks)
    checks["quiet_positive_length_checked"] = True
    checks["raw_sha256_pins_match_without_redecompression"] = True

    all_upper = all_support[-1][1]
    angular_upper = angular_support[-1][1]
    all_margin = Fraction(1) - all_upper
    angular_margin = Fraction(1) - angular_upper
    source_recipe_full_support = (
        all_upper < Fraction(1)
        and angular_upper < Fraction(1)
        and quiet_evidence.checked is True
        and quiet_evidence.expression == EXPECTED_QUIET_POSITIVE_EXPRESSION
    )
    if not source_recipe_full_support:
        _fail("source-recipe support conditions did not close")
    checks["all_reserved_support_before_two_rc_exact"] = True
    checks["angular_m15_support_before_two_rc_exact"] = True
    checks["source_recipe_support_before_rp"] = True

    node_ids = {
        "offset_4893": EXPECTED_QUIET_OFFSET_NODE,
        "reserved_logminus_24955": 24955,
        "reserved_logRc_phys_24956": EXPECTED_RESERVED_LOG_RC_NODE,
        "reserved_log2Rc_24957": EXPECTED_RESERVED_LOG_2RC_NODE,
        "quiet_Rc_25494": EXPECTED_QUIET_RC_NODE,
        "quiet_R2_25495": EXPECTED_QUIET_R2_NODE,
        "quiet_Rp_25499": EXPECTED_QUIET_RP_NODE,
        "quiet_Tw_25007": EXPECTED_QUIET_TW_NODE,
        "quiet_length_25010": EXPECTED_QUIET_LENGTH_NODE,
        "native_logminus_26684": EXPECTED_NATIVE_LOG_MINUS_NODE,
        "native_logRc_phys_26685": EXPECTED_NATIVE_LOG_RC_NODE,
        "native_logRp_phys_26687": EXPECTED_NATIVE_LOG_RP_NODE,
        "native_Tw_26013": EXPECTED_NATIVE_TW_NODE,
    }
    limitations = (
        "The graph/source recipe is checked against caller-supplied raw SHA pins and exact rows; raw content-to-digest binding is performed by the outer runner, and no live producer owner is authenticated.",
        "The reserved and quiet final target roots are in different post-prefix graph suffixes; this report does not promote final-node identity.",
        "The 24952-row common-prefix comparison is only reserved versus quiet_c3; native_frame is not compared row-by-row at that prefix.",
        "Full ordinary-Z C1 matrix-row identity and angular invertibility are not proved here.",
        "The physical m15 source function and its complete coefficient row are not proved here.",
        "Support is source-recipe support only; no numerical radius overlap or source value was used.",
    )
    return ReservedRpSourceGeometryReport(
        schema=SCHEMA,
        raw_sha256=digests,
        checks=checks,
        graphs_compared=("reserved", "quiet_c3"),
        native_prefix_compared=False,
        common_prefix_length=COMMON_PREFIX_LENGTH,
        common_prefix_exact=True,
        node_ids=node_ids,
        expressions=expressions,
        support_endpoints=tuple(
            (f"{left.numerator}/{left.denominator}", f"{right.numerator}/{right.denominator}")
            for left, right in all_support
        ),
        angular_support_endpoints=tuple(
            (f"{left.numerator}/{left.denominator}", f"{right.numerator}/{right.denominator}")
            for left, right in angular_support
        ),
        all_support_upper_log2=f"{all_upper.numerator}/{all_upper.denominator}",
        angular_support_upper_log2=f"{angular_upper.numerator}/{angular_upper.denominator}",
        all_support_margin_to_two_rc_log2=f"{all_margin.numerator}/{all_margin.denominator}",
        angular_support_margin_to_two_rc_log2=f"{angular_margin.numerator}/{angular_margin.denominator}",
        quiet_positive_length_expression=quiet_evidence.expression,
        source_recipe_full_support=True,
        live_owner_proved=False,
        full_c1_matrix_proved=False,
        physical_m15_proved=False,
        limitations=limitations,
    )


# Descriptive aliases keep the narrow consumer discoverable without creating a
# second implementation or a producer entry point.
verify_reserved_rp_source_geometry = audit_reserved_rp_source_geometry
build_reserved_rp_source_geometry_report = audit_reserved_rp_source_geometry


__all__ = [
    "AuthenticatedGraph",
    "COMMON_PREFIX_LENGTH",
    "EXPECTED_ANGULAR_INDICES",
    "EXPECTED_AXIAL_INDICES",
    "EXPECTED_PARAMETER_NODES",
    "EXPECTED_QUIET_POSITIVE_EXPRESSION",
    "EXPECTED_SUPPORT_ENDPOINTS",
    "NATIVE_FRAME_ARCHIVE_RAW_SHA256",
    "QUIET_C3_ARCHIVE_RAW_SHA256",
    "RESERVED_ARCHIVE_RAW_SHA256",
    "QuietLengthEvidence",
    "ReservedRpSourceGeometryReport",
    "SCHEMA",
    "SourceGeometryError",
    "audit_reserved_rp_source_geometry",
    "build_reserved_rp_source_geometry_report",
    "verify_reserved_rp_source_geometry",
]
