"""Raw-g5 lower-collar shear consumer for the two saved sigma intervals.

The runner supplies three already-decoded pieces of evidence:

* the selected raw-program rows for nodes ``1``, ``25085``, ``8``,
  ``35926``, ``35950``, ``35951`` and ``35954``;
* the saved exact ``Fraction`` intervals for nodes ``35949`` (``sigma``) and
  ``35953`` (``sigma_s``); and
* the geometry certificate for the source parameter ``hbB``.

This module deliberately does not read the graph or a receipt, decompress a
file, import a producer, or evaluate either sigma leaf.  It authenticates the
shape of the selected rows and performs only the four new local arithmetic
operations represented by the selected program rows:

``35926 = hbB - 1``
``35950 = 35926 * 35949``
``35951 = 1 + 35950``
``35954 = 35926 * 35953``

All endpoint arithmetic is exact closed-interval ``Fraction`` arithmetic.
The interval for ``hbB`` is an enclosure (the real certificate is
``[0, 1/2]``), never a midpoint or a finite-``N`` point estimate.  The result
is conditional on the caller's selected-row and saved-interval assertions;
it does not prove a live source owner, source-function equality, an
orientation/occurrence identity, the Kchi axial potential, or physical m15.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from typing import Any


SCHEMA = "n1_raw_shear_chi_interval_v1"
RESULT_SCHEMA = SCHEMA + "_result"
COMPONENT = "raw_g5_lower_collar_shear_chi_interval"
AUTHENTICATION_BASIS = "runner_asserted_selected_rows_and_saved_intervals"

GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
EXPECTED_GRAPH_SHA256 = GRAPH_SHA256

HB_NODE_ID = 8
CONSTANT_ONE_NODE_ID = 1
CONSTANT_MINUS_ONE_NODE_ID = 25085
HB_MINUS_ONE_NODE_ID = 35926
SIGMA_PRODUCT_NODE_ID = 35950
CHI_NODE_ID = 35951
SIGMA_S_PRODUCT_NODE_ID = 35954
SIGMA_NODE_ID = 35949
SIGMA_S_NODE_ID = 35953

REQUIRED_PROGRAM_NODE_IDS = (
    CONSTANT_ONE_NODE_ID,
    CONSTANT_MINUS_ONE_NODE_ID,
    HB_NODE_ID,
    HB_MINUS_ONE_NODE_ID,
    SIGMA_PRODUCT_NODE_ID,
    CHI_NODE_ID,
    SIGMA_S_PRODUCT_NODE_ID,
)
OPTIONAL_SAVED_LEAF_NODE_IDS = (SIGMA_NODE_ID, SIGMA_S_NODE_ID)
AUTHORIZED_NODE_IDS = REQUIRED_PROGRAM_NODE_IDS + OPTIONAL_SAVED_LEAF_NODE_IDS

HB_PARAMETER_NAME = "hbB"
REAL_HB_INTERVAL = (Fraction(0), Fraction(1, 2))
EXPECTED_ARITHMETIC_CALL_COUNT = 4

_MISSING = object()


class RawShearChiError(ValueError):
    """Base class for fail-closed selected-row and interval errors."""


class AuthenticationError(RawShearChiError):
    """Caller authentication metadata is missing or inconsistent."""


class GraphRecordError(RawShearChiError):
    """A selected record is not the required raw-program record."""


class SourceFamilyError(RawShearChiError):
    """The graph-node family and geometry-certificate family disagree."""


class IntervalInputError(RawShearChiError):
    """An interval is not an exact, ordered closed Fraction enclosure."""


class ZConditionError(RawShearChiError):
    """The fixed-s/hb ordinary-Z condition is missing or contradictory."""


class ArithmeticBudgetError(RawShearChiError):
    """The bounded new arithmetic budget was exceeded."""


@dataclass(frozen=True)
class FractionInterval:
    """A closed interval with exact rational endpoints."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if not isinstance(self.lower, Fraction) or not isinstance(self.upper, Fraction):
            raise TypeError("FractionInterval endpoints must be Fraction instances")
        if self.lower > self.upper:
            raise ValueError("FractionInterval endpoints must be ordered")

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower

    def contains_interval(self, other: "FractionInterval") -> bool:
        return self.lower <= other.lower and other.upper <= self.upper

    def as_dict(self) -> dict[str, Fraction]:
        return {"lower": self.lower, "upper": self.upper}


def _canonical(value: Any, label: str = "value") -> Any:
    """Return a structural key without relying on object identity."""

    if value is None:
        return ("none", None)
    if type(value) is bool:
        return ("bool", value)
    if type(value) is int:
        return ("int", value)
    if isinstance(value, Fraction):
        return ("fraction", value.numerator, value.denominator)
    if isinstance(value, FractionInterval):
        return ("interval", value.lower.numerator, value.lower.denominator, value.upper.numerator, value.upper.denominator)
    if type(value) is str:
        return ("str", value)
    if isinstance(value, Mapping):
        items: list[tuple[str, Any]] = []
        for key, item in value.items():
            if not isinstance(key, str):
                raise GraphRecordError(f"{label} has a non-string key")
            items.append((key, _canonical(item, f"{label}.{key}")))
        return ("mapping", tuple(sorted(items)))
    if isinstance(value, list):
        return ("list", tuple(_canonical(item, f"{label}[{i}]") for i, item in enumerate(value)))
    if isinstance(value, tuple):
        return ("tuple", tuple(_canonical(item, f"{label}[{i}]") for i, item in enumerate(value)))
    raise GraphRecordError(f"{label} is not a supported structural value")


def _same(left: Any, right: Any, label: str) -> bool:
    try:
        return _canonical(left, label + ".left") == _canonical(right, label + ".right")
    except RawShearChiError:
        raise
    except Exception as exc:  # pragma: no cover - hostile caller boundary
        raise AuthenticationError(f"{label} values are not comparable") from exc


def _alias(mapping: Mapping[str, Any], names: tuple[str, ...], label: str) -> Any:
    found = [(name, mapping[name]) for name in names if name in mapping]
    if not found:
        return _MISSING
    first_name, first = found[0]
    for name, value in found[1:]:
        if not _same(first, value, label):
            raise AuthenticationError(f"conflicting {label} aliases: {first_name}/{name}")
    return first


def _sha256(value: Any, label: str, *, synthetic: bool) -> str:
    if not isinstance(value, str) or not value:
        raise AuthenticationError(f"{label} must be a non-empty string")
    value = value.lower()
    if synthetic:
        return value
    if len(value) != 64 or any(ch not in "0123456789abcdef" for ch in value):
        raise AuthenticationError(f"{label} must be a SHA-256 hex string")
    if value != EXPECTED_GRAPH_SHA256:
        raise AuthenticationError("unexpected original raw-g5 graph SHA-256")
    return value


def _as_node_id(value: Any, label: str) -> int:
    if type(value) is not int:
        raise GraphRecordError(f"{label} must be an integer node id")
    return value


def _binary_mpf_fraction(value: Any, label: str) -> Fraction | object:
    """Decode the saved four-integer binary endpoint tuple, if present.

    The runner normally converts these tuples before calling this module.  A
    narrow decoder for the saved endpoint representation makes the boundary
    usable with the plain geometry metadata while still refusing floats.
    """

    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        return _MISSING
    if len(value) != 4 or any(type(item) is not int for item in value):
        return _MISSING
    sign, mantissa, exponent, precision = value
    if sign not in (0, 1) or mantissa < 0 or precision < 0:
        raise IntervalInputError(f"{label} is not a valid binary endpoint tuple")
    if mantissa == 0:
        return Fraction(0)
    if abs(exponent) > 1_000_000:
        raise IntervalInputError(f"{label} exponent exceeds the exact-rational bound")
    numerator = (-mantissa if sign else mantissa)
    if exponent >= 0:
        return Fraction(numerator * (2**exponent), 1)
    return Fraction(numerator, 2**(-exponent))


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    if isinstance(value, str):
        try:
            return Fraction(value.strip())
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise IntervalInputError(f"{label} is not an exact rational string") from exc
    if isinstance(value, Mapping):
        numerator = _alias(value, ("numerator", "num"), label + ".numerator")
        denominator = _alias(value, ("denominator", "den"), label + ".denominator")
        if numerator is not _MISSING and denominator is not _MISSING:
            if type(numerator) is not int or type(denominator) is not int or denominator == 0:
                raise IntervalInputError(f"{label} numerator/denominator are invalid")
            return Fraction(numerator, denominator)
        nested = _alias(value, ("value", "exact", "fraction"), label + ".nested")
        if nested is not _MISSING:
            return _fraction(nested, label)
    binary = _binary_mpf_fraction(value, label)
    if binary is not _MISSING:
        return binary  # type: ignore[return-value]
    if isinstance(value, float):
        raise IntervalInputError(f"{label} must not be a floating-point value")
    raise IntervalInputError(f"{label} is not an exact rational endpoint")


def _interval_metadata(value: Mapping[str, Any], label: str) -> None:
    for key in ("normalized", "is_normalized", "normalized_cutoff", "midpoint_used"):
        if key in value and value[key] not in (False, None):
            raise IntervalInputError(f"{label} contains forbidden point/normalized metadata: {key}")
    for key in ("finite_N", "finite_n", "finite_N_error_bound", "point_value"):
        if key in value and value[key] not in (False, None):
            raise IntervalInputError(f"{label} contains forbidden finite-N or point metadata: {key}")


def _fraction_interval(value: Any, label: str) -> FractionInterval:
    if isinstance(value, FractionInterval):
        return value
    if isinstance(value, Mapping):
        _interval_metadata(value, label)
        nested = _alias(value, ("interval", "enclosure", "bounds", "box"), label + ".nested")
        if nested is not _MISSING:
            return _fraction_interval(nested, label)
        lower = _alias(
            value,
            ("lower", "lo", "left", "lower_fraction", "lower_exact_mpf_tuple"),
            label + ".lower",
        )
        upper = _alias(
            value,
            ("upper", "hi", "right", "upper_fraction", "upper_exact_mpf_tuple"),
            label + ".upper",
        )
        if lower is _MISSING or upper is _MISSING:
            raise IntervalInputError(f"{label} requires lower and upper endpoints")
        return FractionInterval(_fraction(lower, label + ".lower"), _fraction(upper, label + ".upper"))
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        if len(value) != 2:
            raise IntervalInputError(f"{label} sequence must contain two endpoints")
        return FractionInterval(_fraction(value[0], label + ".lower"), _fraction(value[1], label + ".upper"))
    lower = getattr(value, "lower", _MISSING)
    upper = getattr(value, "upper", _MISSING)
    if lower is not _MISSING and upper is not _MISSING:
        return FractionInterval(_fraction(lower, label + ".lower"), _fraction(upper, label + ".upper"))
    # A scalar is accepted only as a degenerate exact enclosure.  The real
    # runner passes two endpoints; this convenience does not introduce a
    # midpoint or a numerical approximation.
    scalar = _fraction(value, label)
    return FractionInterval(scalar, scalar)


def _recursive_alias(value: Any, names: tuple[str, ...], label: str, *, depth: int = 4) -> Any:
    """Find a certificate field at the few saved geometry wrapper levels."""

    if depth < 0 or not isinstance(value, Mapping):
        return _MISSING
    direct = _alias(value, names, label)
    found: list[tuple[str, Any]] = []
    if direct is not _MISSING:
        found.append(("direct", direct))
    for key in (
        "geometry_certificate",
        "certificate",
        "geometry",
        "raw",
        "actual_local_j2_tail",
        "output",
        "result",
        "actual_result",
        "consumer_result",
        "direct_lower_collar",
        "selected",
        "scope",
    ):
        nested = value.get(key, _MISSING)
        if nested is _MISSING or not isinstance(nested, Mapping):
            continue
        item = _recursive_alias(nested, names, label, depth=depth - 1)
        if item is not _MISSING:
            found.append((key, item))
    if not found:
        return _MISSING
    first_name, first = found[0]
    for name, item in found[1:]:
        if not _same(first, item, label):
            raise AuthenticationError(f"conflicting {label} certificate fields: {first_name}/{name}")
    return first


def _contains_structural(value: Any, target: Any, *, depth: int = 8) -> bool:
    """Check a caller-supplied family receipt without treating IDs as owners."""

    if _same(value, target, "source family receipt"):
        return True
    if depth <= 0:
        return False
    if isinstance(value, Mapping):
        return any(_contains_structural(item, target, depth=depth - 1) for item in value.values())
    if isinstance(value, (list, tuple)):
        return any(_contains_structural(item, target, depth=depth - 1) for item in value)
    return False


def _row_payload(row: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(row, Mapping):
        raise GraphRecordError(f"{label} must be a mapping")
    raw = row.get("raw_record", _MISSING)
    if raw is _MISSING:
        return row
    if not isinstance(raw, Mapping):
        raise GraphRecordError(f"{label}.raw_record must be a mapping")
    # Keep selected-row metadata (notably source_family) available while
    # giving the operation/arguments fields precedence to the raw record.
    merged = dict(row)
    merged.update(raw)
    return merged


def _normalize_rows(rows: Any) -> dict[int, Mapping[str, Any]]:
    normalized: dict[int, Mapping[str, Any]] = {}
    if isinstance(rows, Mapping):
        iterable = rows.items()
        for raw_id, row in iterable:
            if isinstance(raw_id, str) and raw_id.isdigit():
                node_id = int(raw_id)
            elif type(raw_id) is int:
                node_id = raw_id
            else:
                # Non-numeric keys are packet metadata, not selected rows.
                continue
            if node_id in normalized:
                raise GraphRecordError(f"duplicate selected node {node_id}")
            normalized[node_id] = _row_payload(row, f"selected row {node_id}")
    elif isinstance(rows, Sequence) and not isinstance(rows, (str, bytes, bytearray)):
        for index, row in enumerate(rows):
            if not isinstance(row, Mapping):
                raise GraphRecordError(f"selected row {index} must be a mapping")
            node_value = _alias(row, ("node_id", "id", "node"), f"selected row {index}.node_id")
            if node_value is _MISSING:
                raise GraphRecordError(f"selected row {index} is missing node_id")
            node_id = _as_node_id(node_value, f"selected row {index}.node_id")
            if node_id in normalized:
                raise GraphRecordError(f"duplicate selected node {node_id}")
            normalized[node_id] = _row_payload(row, f"selected row {node_id}")
    else:
        raise GraphRecordError("selected rows must be a mapping or a sequence")
    missing = [node_id for node_id in REQUIRED_PROGRAM_NODE_IDS if node_id not in normalized]
    if missing:
        raise GraphRecordError(f"selected rows are missing required nodes: {missing}")
    extras = sorted(set(normalized) - set(AUTHORIZED_NODE_IDS))
    if extras:
        raise GraphRecordError(f"selected rows contain unauthorized nodes: {extras}")
    return normalized


def _operation(row: Mapping[str, Any], node_id: int) -> str:
    value = _alias(row, ("operation", "op"), f"node {node_id} operation")
    if not isinstance(value, str):
        raise GraphRecordError(f"node {node_id} operation must be a string")
    return value


def _arguments(row: Mapping[str, Any], node_id: int) -> list[int]:
    value = _alias(row, ("arguments", "operands", "terms", "factors"), f"node {node_id} arguments")
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise GraphRecordError(f"node {node_id} arguments must be a sequence")
    result = [_as_node_id(item, f"node {node_id} argument") for item in value]
    if len(result) != 2:
        raise GraphRecordError(f"node {node_id} must have exactly two arguments")
    return result


def _unary_argument(row: Mapping[str, Any], node_id: int) -> int:
    value = _alias(row, ("argument", "operand"), f"node {node_id} argument")
    if value is not _MISSING:
        return _as_node_id(value, f"node {node_id} argument")
    value = _alias(row, ("arguments", "operands"), f"node {node_id} arguments")
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)) or len(value) != 1:
        raise GraphRecordError(f"node {node_id} requires one unary argument")
    return _as_node_id(value[0], f"node {node_id} argument")


def _constant_value(row: Mapping[str, Any], node_id: int) -> Fraction:
    operation = _operation(row, node_id)
    if operation != "exact_rational":
        raise GraphRecordError(f"node {node_id} must use exact_rational, got {operation!r}")
    value = _alias(row, ("value", "number", "rational", "constant"), f"node {node_id} value")
    if value is _MISSING:
        numerator = _alias(row, ("numerator", "num"), f"node {node_id}.numerator")
        denominator = _alias(row, ("denominator", "den"), f"node {node_id}.denominator")
        if numerator is _MISSING or denominator is _MISSING:
            raise GraphRecordError(f"node {node_id} is missing exact rational value")
        value = {"numerator": numerator, "denominator": denominator}
    try:
        return _fraction(value, f"node {node_id} value")
    except IntervalInputError as exc:
        raise GraphRecordError(str(exc)) from exc


def _row_node_id(row: Mapping[str, Any], node_id: int) -> None:
    value = _alias(row, ("node_id", "id", "node"), f"node {node_id} identity")
    if value is not _MISSING and _as_node_id(value, f"node {node_id} identity") != node_id:
        raise GraphRecordError(f"node record identity disagrees with selected id {node_id}")


def _validate_program_rows(rows: dict[int, Mapping[str, Any]]) -> Any:
    for node_id, row in rows.items():
        _row_node_id(row, node_id)
    if _constant_value(rows[CONSTANT_ONE_NODE_ID], CONSTANT_ONE_NODE_ID) != Fraction(1):
        raise GraphRecordError("node 1 is not exact constant 1")
    minus_one = rows[CONSTANT_MINUS_ONE_NODE_ID]
    if _operation(minus_one, CONSTANT_MINUS_ONE_NODE_ID) != "negative":
        raise GraphRecordError("node 25085 must be the saved negative node")
    if _unary_argument(minus_one, CONSTANT_MINUS_ONE_NODE_ID) != CONSTANT_ONE_NODE_ID:
        raise GraphRecordError("node 25085 must negate node 1")

    hb_row = rows[HB_NODE_ID]
    hb_operation = _operation(hb_row, HB_NODE_ID)
    if hb_operation not in {
        "current_original_source_parameter",
        "original_source_parameter",
        "source_parameter",
    }:
        raise GraphRecordError(f"node 8 is not the original hbB source parameter: {hb_operation!r}")
    hb_name = _alias(hb_row, ("name", "parameter", "symbol"), "node 8 parameter name")
    if hb_name != HB_PARAMETER_NAME:
        raise GraphRecordError("node 8 source parameter is not hbB")

    expected = {
        HB_MINUS_ONE_NODE_ID: ("sum", [HB_NODE_ID, CONSTANT_MINUS_ONE_NODE_ID]),
        SIGMA_PRODUCT_NODE_ID: ("product", [HB_MINUS_ONE_NODE_ID, SIGMA_NODE_ID]),
        CHI_NODE_ID: ("sum", [CONSTANT_ONE_NODE_ID, SIGMA_PRODUCT_NODE_ID]),
        SIGMA_S_PRODUCT_NODE_ID: ("product", [HB_MINUS_ONE_NODE_ID, SIGMA_S_NODE_ID]),
    }
    for node_id, (operation, arguments) in expected.items():
        actual_operation = _operation(rows[node_id], node_id)
        if actual_operation != operation:
            raise GraphRecordError(f"node {node_id} operation must be {operation!r}")
        if _arguments(rows[node_id], node_id) != arguments:
            raise GraphRecordError(f"node {node_id} arguments do not match the saved raw program")

    family = _alias(hb_row, ("source_family", "family"), "node 8 source family")
    if family is _MISSING and isinstance(hb_row.get("source", _MISSING), Mapping):
        family = _alias(hb_row["source"], ("source_family", "family"), "node 8 nested source family")
    return family


def _graph_source_family(rows: dict[int, Mapping[str, Any]], metadata: Mapping[str, Any], supplied: Any) -> Any:
    row_family = _validate_program_rows(rows)
    metadata_family = _alias(metadata, ("source_family", "graph_source_family"), "graph source family")
    candidates = [item for item in (row_family, metadata_family, supplied) if item is not _MISSING and item is not None]
    if not candidates:
        raise SourceFamilyError("graph node 8 source_family is required")
    first = candidates[0]
    for candidate in candidates[1:]:
        if not _same(first, candidate, "graph source family"):
            raise SourceFamilyError("graph source family aliases disagree")
    return first


def _coalesce_field(current: Any, mapping: Mapping[str, Any], names: tuple[str, ...], label: str) -> Any:
    candidate = _alias(mapping, names, label)
    if current is _MISSING:
        return candidate
    if candidate is not _MISSING and not _same(current, candidate, label):
        raise AuthenticationError(f"conflicting {label} packet and argument values")
    return current


def _extract_packet(selected_rows: Any) -> tuple[Any, Mapping[str, Any]]:
    if not isinstance(selected_rows, Mapping):
        return selected_rows, {}
    rows = _alias(
        selected_rows,
        ("selected_rows", "program_rows", "graph_rows", "rows", "records"),
        "selected rows",
    )
    if rows is _MISSING:
        return selected_rows, selected_rows
    return rows, selected_rows


def _interval_for_saved_leaf(value: Any, node_id: int, quantity: str) -> FractionInterval:
    if isinstance(value, Mapping):
        marker = _alias(value, ("node_id", "node"), f"{quantity} node id")
        if marker is not _MISSING and _as_node_id(marker, f"{quantity} node id") != node_id:
            raise GraphRecordError(f"{quantity} interval is attached to the wrong node")
        kind = _alias(value, ("quantity", "kind"), f"{quantity} quantity")
        if kind is not _MISSING and kind not in (quantity, "sigma_s" if quantity == "sigma_s" else quantity):
            raise GraphRecordError(f"saved interval quantity is not {quantity}")
    return _fraction_interval(value, f"saved {quantity} interval")


def _validate_hb_evidence(evidence: Any, *, synthetic: bool) -> tuple[FractionInterval, Any, dict[str, Any]]:
    if not isinstance(evidence, Mapping):
        raise AuthenticationError("hb evidence must be a mapping")
    strict = _recursive_alias(evidence, ("strict_hb_admitted",), "strict_hb_admitted")
    independent = _recursive_alias(
        evidence,
        ("hb_sc_Z_independent", "hb_epsilon_Z_independent"),
        "hb_sc_Z_independent",
    )
    if strict is not True:
        raise AuthenticationError("hb evidence requires strict_hb_admitted=true")
    if independent is not True:
        raise AuthenticationError("hb evidence requires hb_sc_Z_independent=true")
    raw_interval = _recursive_alias(evidence, ("hb_interval",), "hb_interval")
    if raw_interval is _MISSING:
        raise IntervalInputError("hb evidence is missing hb_interval")
    hb_interval = _fraction_interval(raw_interval, "hb_interval")
    if not synthetic and (hb_interval.lower, hb_interval.upper) != REAL_HB_INTERVAL:
        raise IntervalInputError("real hb evidence must enclose exactly [0, 1/2]")
    family = _recursive_alias(evidence, ("source_family", "family"), "hb source family")
    normalized = {
        "strict_hb_admitted": True,
        "hb_sc_Z_independent": True,
        "hb_interval": hb_interval.as_dict(),
        "closed_hb_enclosure_allowed": True,
        "midpoint_used": False,
    }
    if family is not _MISSING:
        normalized["source_family"] = family
    return hb_interval, family, normalized


def _validate_z_condition(condition: Any, evidence: Any) -> dict[str, Any]:
    source = condition
    if source is _MISSING:
        source = evidence
    if not isinstance(source, Mapping):
        raise ZConditionError("fixed-s/hb Z condition must be a mapping")

    fixed_s = _recursive_alias(source, ("fixed_s", "s_fixed"), "fixed_s")
    fixed_hb = _recursive_alias(source, ("fixed_hb", "hb_fixed"), "fixed_hb")
    fixed_pair = _recursive_alias(source, ("fixed_s_hb", "fixed_s_and_hb", "ordinary_Z_geometry_fixed"), "fixed_s_hb")
    if fixed_pair is not _MISSING:
        if fixed_pair is not True:
            raise ZConditionError("the ordinary-Z geometry is not fixed in s and hb")
        fixed_s = True if fixed_s is _MISSING else fixed_s
        fixed_hb = True if fixed_hb is _MISSING else fixed_hb
    if fixed_s is not True or fixed_hb is not True:
        raise ZConditionError("fixed_s=true and fixed_hb=true are required")

    z1 = _recursive_alias(source, ("ordinary_Z1", "ordinary_z1", "Z1", "z1"), "ordinary_Z1")
    if z1 is _MISSING:
        # The saved geometry certificate records a fixed geometry and no
        # ordinary-Z differentiation; the explicit value remains conditional.
        z1 = Fraction(0)
    else:
        if _fraction(z1, "ordinary_Z1") != Fraction(0):
            raise ZConditionError("ordinary Z1 must be exactly zero under fixed s,hb")

    derivative = _recursive_alias(
        source,
        ("s_derivative_is_Z_derivative", "ordinary_s_derivative_is_Z_derivative"),
        "s/Z derivative relation",
    )
    not_derivative = _recursive_alias(
        source,
        ("s_derivative_not_Z_derivative",),
        "s/Z derivative separation",
    )
    differentiated = _recursive_alias(
        source,
        ("ordinary_Z_geometry_differentiated",),
        "ordinary-Z geometry differentiation",
    )
    if derivative is True:
        raise ZConditionError("the s derivative must not be treated as a Z derivative")
    if not_derivative is False:
        raise ZConditionError("the s derivative must not be treated as a Z derivative")
    if differentiated is True:
        raise ZConditionError("the saved geometry must not be Z-differentiated")
    q_relation = _recursive_alias(
        source,
        ("s_derivative_is_q_derivative", "q_derivative_is_s_derivative", "chi_s_is_q_derivative"),
        "s/q derivative relation",
    )
    if q_relation is True:
        raise ZConditionError("the s derivative must not be treated as a q derivative")
    if derivative is _MISSING and not_derivative is _MISSING and differentiated is _MISSING:
        raise ZConditionError("an explicit s-versus-Z derivative separation is required")
    if not_derivative is True or differentiated is False:
        derivative = False
    return {
        "fixed_s": True,
        "fixed_hb": True,
        "ordinary_Z1": Fraction(0),
        "s_derivative_is_Z_derivative": False,
        "conditional": True,
    }


class _IntervalBackend:
    """Private four-call exact interval backend used by the consumer."""

    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []

    def _record(self, operation: str, node_id: int, arguments: list[int], interval: FractionInterval) -> FractionInterval:
        if len(self.calls) >= EXPECTED_ARITHMETIC_CALL_COUNT:
            raise ArithmeticBudgetError("new arithmetic call budget exceeded")
        self.calls.append(
            {
                "node_id": node_id,
                "operation": operation,
                "arguments": list(arguments),
                "result_interval": interval.as_dict(),
            }
        )
        return interval

    def add(self, left: FractionInterval, right: FractionInterval, node_id: int, arguments: list[int]) -> FractionInterval:
        return self._record(
            "sum",
            node_id,
            arguments,
            FractionInterval(left.lower + right.lower, left.upper + right.upper),
        )

    def multiply(self, left: FractionInterval, right: FractionInterval, node_id: int, arguments: list[int]) -> FractionInterval:
        products = (
            left.lower * right.lower,
            left.lower * right.upper,
            left.upper * right.lower,
            left.upper * right.upper,
        )
        return self._record("product", node_id, arguments, FractionInterval(min(products), max(products)))


def compute(
    selected_rows: Any,
    sigma_interval: Any = _MISSING,
    sigma_s_interval: Any = _MISSING,
    hb_evidence: Any = _MISSING,
    *,
    graph_hash: Any = _MISSING,
    graph_sha256: Any = _MISSING,
    source_family: Any = _MISSING,
    z_condition: Any = _MISSING,
    rows_authenticated: Any = _MISSING,
    authenticated: Any = _MISSING,
    sigma_bounds_authenticated: Any = _MISSING,
    source_family_receipt: Any = _MISSING,
    synthetic: bool = False,
    **aliases: Any,
) -> dict[str, Any]:
    """Consume authenticated selected rows and return exact shear intervals.

    ``selected_rows`` may be the row mapping itself or a packet containing a
    ``selected_rows``/``program_rows`` field and the same metadata.  The
    runner must supply a graph hash and an authentication marker.  The two
    saved sigma intervals are accepted as prior results and are never
    recomputed here.
    """

    if type(synthetic) is not bool:
        raise AuthenticationError("synthetic must be a boolean")
    rows_raw, packet = _extract_packet(selected_rows)

    sigma_interval = _coalesce_field(sigma_interval, packet, ("sigma_interval", "sigma"), "sigma interval")
    sigma_s_interval = _coalesce_field(
        sigma_s_interval,
        packet,
        ("sigma_s_interval", "sigma_s", "d_sigma_ds", "derivative"),
        "sigma_s interval",
    )
    hb_evidence = _coalesce_field(
        hb_evidence,
        packet,
        ("hb_evidence", "geometry_certificate", "geometry", "certificate"),
        "hb evidence",
    )
    graph_hash = _coalesce_field(graph_hash, packet, ("graph_hash", "graph_sha256"), "graph hash")
    graph_hash = _coalesce_field(graph_hash, {"graph_hash": graph_sha256} if graph_sha256 is not _MISSING else {}, ("graph_hash",), "graph hash")
    source_family = _coalesce_field(source_family, packet, ("source_family", "graph_source_family"), "source family")
    source_family_receipt = _coalesce_field(
        source_family_receipt,
        packet,
        ("source_family_receipt", "family_receipt"),
        "source family receipt",
    )
    z_condition = _coalesce_field(z_condition, packet, ("z_condition", "Z_condition", "ordinary_Z_condition"), "Z condition")
    rows_authenticated = _coalesce_field(
        rows_authenticated,
        packet,
        ("rows_authenticated", "selected_rows_authenticated", "runner_authenticated", "authenticated"),
        "selected-row authentication",
    )
    if authenticated is not _MISSING:
        rows_authenticated = _coalesce_field(
            rows_authenticated,
            {"authenticated": authenticated},
            ("authenticated",),
            "selected-row authentication",
        )
    sigma_bounds_authenticated = _coalesce_field(
        sigma_bounds_authenticated,
        packet,
        ("sigma_bounds_authenticated", "saved_intervals_authenticated", "intervals_authenticated"),
        "saved-interval authentication",
    )
    for name in ("sigma", "sigma_s", "hb", "z", "Z"):
        if name in aliases:
            if name in {"sigma", "sigma_s"}:
                # The canonical positional arguments are preferred; aliases
                # are accepted only when the corresponding argument was absent.
                target = sigma_interval if name == "sigma" else sigma_s_interval
                if target is _MISSING:
                    if name == "sigma":
                        sigma_interval = aliases[name]
                    else:
                        sigma_s_interval = aliases[name]
                elif not _same(target, aliases[name], name):
                    raise AuthenticationError(f"conflicting {name} interval arguments")
            elif name == "hb":
                if hb_evidence is _MISSING:
                    hb_evidence = aliases[name]
                elif not _same(hb_evidence, aliases[name], "hb evidence"):
                    raise AuthenticationError("conflicting hb evidence arguments")
            else:
                if z_condition is _MISSING:
                    z_condition = aliases[name]
                elif not _same(z_condition, aliases[name], "Z condition"):
                    raise AuthenticationError("conflicting Z condition arguments")
    unknown_aliases = set(aliases) - {"sigma", "sigma_s", "hb", "z", "Z"}
    if unknown_aliases:
        raise AuthenticationError(f"unknown consumer arguments: {sorted(unknown_aliases)}")

    if rows_authenticated is not True:
        raise AuthenticationError("selected rows require an explicit runner authentication marker")
    if sigma_bounds_authenticated is not True:
        raise AuthenticationError("saved sigma intervals require an explicit runner authentication marker")
    graph_hash = _sha256(graph_hash, "graph_hash", synthetic=synthetic)
    rows = _normalize_rows(rows_raw)
    metadata = packet if isinstance(packet, Mapping) else {}
    graph_family = _graph_source_family(rows, metadata, source_family)
    if source_family_receipt is not _MISSING:
        if not isinstance(source_family_receipt, Mapping) or not _contains_structural(source_family_receipt, graph_family):
            raise SourceFamilyError("source family receipt does not contain graph node 8 source_family")

    if sigma_interval is _MISSING or sigma_s_interval is _MISSING or hb_evidence is _MISSING:
        raise AuthenticationError("sigma, sigma_s, and hb evidence are all required")
    sigma = _interval_for_saved_leaf(sigma_interval, SIGMA_NODE_ID, "sigma")
    sigma_s = _interval_for_saved_leaf(sigma_s_interval, SIGMA_S_NODE_ID, "sigma_s")
    if not synthetic and (sigma.lower < 0 or sigma_s.lower < 0):
        raise IntervalInputError("real saved sigma and sigma_s intervals must be nonnegative")
    hb, certificate_family, certificate_metadata = _validate_hb_evidence(hb_evidence, synthetic=synthetic)
    if not synthetic and certificate_family is _MISSING:
        raise SourceFamilyError("real hb evidence must carry source_family")
    if certificate_family is not _MISSING and not _same(graph_family, certificate_family, "source family"):
        raise SourceFamilyError("hb certificate source_family does not match graph node 8")
    z = _validate_z_condition(z_condition, hb_evidence)

    backend = _IntervalBackend()
    hb_minus_one = backend.add(hb, FractionInterval(Fraction(-1), Fraction(-1)), HB_MINUS_ONE_NODE_ID, [HB_NODE_ID, CONSTANT_MINUS_ONE_NODE_ID])
    sigma_product = backend.multiply(hb_minus_one, sigma, SIGMA_PRODUCT_NODE_ID, [HB_MINUS_ONE_NODE_ID, SIGMA_NODE_ID])
    chi = backend.add(FractionInterval(Fraction(1), Fraction(1)), sigma_product, CHI_NODE_ID, [CONSTANT_ONE_NODE_ID, SIGMA_PRODUCT_NODE_ID])
    sigma_s_product = backend.multiply(hb_minus_one, sigma_s, SIGMA_S_PRODUCT_NODE_ID, [HB_MINUS_ONE_NODE_ID, SIGMA_S_NODE_ID])
    if len(backend.calls) != EXPECTED_ARITHMETIC_CALL_COUNT:
        raise ArithmeticBudgetError("the consumer did not perform exactly the designated four arithmetic calls")

    source_family_out = graph_family
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": "conditional_raw_g5_lower_collar_shear_chi_interval",
        "graph_sha256": graph_hash,
        "selected_node_ids": list(REQUIRED_PROGRAM_NODE_IDS),
        "saved_sigma_node_ids": {"sigma": SIGMA_NODE_ID, "sigma_s": SIGMA_S_NODE_ID},
        "selected_rows_authenticated": True,
        "saved_intervals_authenticated": True,
        "authentication_basis": AUTHENTICATION_BASIS,
        "fraction_only_new_nodes": True,
        "reused_sigma_intervals": True,
        "hb_Z_independent": True,
        "source_family": source_family_out,
        "source_family_receipt_checked": source_family_receipt is not _MISSING,
        "source_identity_verified": False,
        "source_function_equality_proved": False,
        "owner_binding": False,
        "live_owner_verified": False,
        "orientation_preserved": True,
        "orientation_identity_proved": False,
        "hb_evidence": certificate_metadata,
        "hb_interval_is_parameter_enclosure": True,
        "z_condition": z,
        "ordinary_Z1": {"value": Fraction(0), "conditional": True, "fixed_s_hb": True},
        "coordinate": "raw_s",
        "s_derivative_is_Z_derivative": False,
        "s_derivative_is_q_derivative": False,
        "z_derivative": False,
        "q_derivative": False,
        "chi_s_is_Z_derivative": False,
        "chi_s_is_q_derivative": False,
        "ordinary_s_derivative": True,
        "lower_collar_shear": True,
        "quantity_role": "lower_collar_shear_chi",
        "not_Kchi_axial_potential": True,
        "not_physical_m15": True,
        "physical_m15": False,
        "full_g5": False,
        "sigma_recomputed": False,
        "sigma_s_recomputed": False,
        "prior_sigma_intervals_reused": True,
        "midpoint_used": False,
        "finite_N_error_bound_used": False,
        "arithmetic_call_count": len(backend.calls),
        "allowed_arithmetic_call_count": EXPECTED_ARITHMETIC_CALL_COUNT,
        "operation_records": backend.calls,
        "intervals": {
            "hb": hb.as_dict(),
            "hb_minus_one_35926": hb_minus_one.as_dict(),
            "sigma_35949": sigma.as_dict(),
            "sigma_s_35953": sigma_s.as_dict(),
            "chi_35951": chi.as_dict(),
            "chi_s_35954": sigma_s_product.as_dict(),
        },
        "chi_interval": chi.as_dict(),
        "chi_s_interval": sigma_s_product.as_dict(),
        "chi": chi.as_dict(),
        "chi_s": sigma_s_product.as_dict(),
        "35951": chi.as_dict(),
        "35954": sigma_s_product.as_dict(),
        "node_results": {
            str(CHI_NODE_ID): {"node_id": CHI_NODE_ID, "quantity": "chi", "interval": chi.as_dict()},
            str(SIGMA_S_PRODUCT_NODE_ID): {
                "node_id": SIGMA_S_PRODUCT_NODE_ID,
                "quantity": "chi_s",
                "interval": sigma_s_product.as_dict(),
            },
        },
        "synthetic_only": synthetic,
        "used_old_module": False,
        "used_old_fixture": False,
        "gzip_reads": False,
    }


consume = compute
evaluate = compute
build = compute
consume_selected_rows = compute


__all__ = [
    "SCHEMA",
    "RESULT_SCHEMA",
    "COMPONENT",
    "AUTHENTICATION_BASIS",
    "GRAPH_SHA256",
    "EXPECTED_GRAPH_SHA256",
    "HB_NODE_ID",
    "CONSTANT_ONE_NODE_ID",
    "CONSTANT_MINUS_ONE_NODE_ID",
    "HB_MINUS_ONE_NODE_ID",
    "SIGMA_PRODUCT_NODE_ID",
    "CHI_NODE_ID",
    "SIGMA_S_PRODUCT_NODE_ID",
    "SIGMA_NODE_ID",
    "SIGMA_S_NODE_ID",
    "REQUIRED_PROGRAM_NODE_IDS",
    "OPTIONAL_SAVED_LEAF_NODE_IDS",
    "AUTHORIZED_NODE_IDS",
    "HB_PARAMETER_NAME",
    "REAL_HB_INTERVAL",
    "EXPECTED_ARITHMETIC_CALL_COUNT",
    "RawShearChiError",
    "AuthenticationError",
    "GraphRecordError",
    "SourceFamilyError",
    "IntervalInputError",
    "ZConditionError",
    "ArithmeticBudgetError",
    "FractionInterval",
    "compute",
    "consume",
    "evaluate",
    "build",
    "consume_selected_rows",
]
