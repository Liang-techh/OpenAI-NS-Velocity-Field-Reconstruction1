"""Bounded source-only congruence for the saved Rm loop and density DAG.

This module is intentionally a consumer.  It accepts one already admitted
generic native packet, explicit chart bindings, and a source receipt for the
common phase.  It never imports a producer, checker, fixture, archive loader,
quadrature, or numerical evaluator.  The only values it forms are exact
rationals and symbolic source expressions used to compare the two saved
charts at their supplied native endpoints.

The public entry point is :func:`derive`.  It compares ordinary rows ``0``
and ``1`` of ``A``, ``B_over_Pstar``, ``A_phi``,
``B_over_Pstar_phi``, finite ``E_N`` and ``V_N``, and the five signed density
increments.  Leading
generic roots are admitted only through the accepted packet's explicit
``root_rows``.  A root row's saved canonical ``value.srepr`` is the canonical
atom; this deliberately makes equal values such as ``a`` and ``kappa`` share
an atom while retaining an exact zero as the exact rational zero.

Unknown reachable operations fail closed.  Compound source operations are
represented by hash atoms keyed by their canonical child expressions and all
mathematical descriptor inputs.  This permits exact arithmetic around an
analytic source operation without pretending to evaluate that operation.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

try:
    import sympy as _sp
except ImportError:  # pragma: no cover
    _sp = None


ORDERS = (0, 1)
GENERIC_ROOTS = (
    "E",
    "V",
    "a",
    "b",
    "p1",
    "p2",
    "t0",
    "kappa",
    "Delta",
)
LOOP_FIELDS = ("A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi")
FINITE_FIELDS = ("E_N", "V_N")
DENSITY_FIELDS = ("m", "h", "k", "e", "p")
TARGET_FIELDS = LOOP_FIELDS + FINITE_FIELDS + DENSITY_FIELDS
SIDES = ("postrestore", "patch")

_SIDE_ALIASES = {
    "postrestore": ("postrestore", "post_restore", "modified", "left"),
    "patch": ("patch", "finite", "right"),
}
_FUNCTION_GUARDS = (
    "Z_independent_substitution",
    "native_substitution_guard",
    "exact_native_chain_rule_substitution",
    "endpoint_substitution_guard",
    "guarded",
    "inverse_substitution_guard",
)
_ANALYTIC = {"exp", "expm1", "log", "log1p", "sqrt", "positive_sqrt", "sin", "cos"}
_OPAQUE_ANALYTIC = {
    "fractional_part",
    "gamma",
    "exprel",
    "log1prel",
    "original_flat_sigma",
    "original_flat_sigma_prime",
    "original_flat_sigma_second",
    "original_flat_sigma_third",
    "compact_raw_beta_derivative",
}
_COMPOUND_OPS = {
    "original_monotone_phase_inverse",
    "original_lazy_flat_branch",
    "exact_original_source_flat_sigma",
    "substitute_original_inverse_angle",
    "substitute_original_native_coordinate",
    "substitute_free_query_at_native_coordinate",
    "inverse_substitution",
    "native_substitution",
    "logical_or",
    "logical_and",
    "logical_not",
    "trim",
    "exact_real_comparison",
    "regular_even_radial_axis_continuation",
    "implicit_C1_Z_linear_system",
    "implicit_C2_ZZ_linear_system",
    "implicit_C3_ZZZ_linear_system",
    "source_bound_C1_Banach_limit",
    "C1_Banach_limit_component",
    "C2_Banach_limit_second_component",
    "C3_Banach_limit_third_component",
    "current_original_axis_pressure_function_derivative",
    "current_original_leading_function_recipe",
    "current_original_physical_radius",
    "current_original_physical_source_quantity",
    "current_original_physical_source_scale",
    "current_original_source_parameter",
    "current_original_unique_analytic_core_component",
    "current_original_unique_regular_n1_component",
    "original_Rm_unique_leading_implicit_control",
    "original_Rm_unique_leading_implicit_system",
    "original_compact_unique_Banach_limit_component",
    "original_compact_unique_Banach_limit_system",
    "proved_C1_zero_identity",
    "proved_C2_relative_terminal_zero_identity",
    "proved_C2_second_zero_identity",
    "proved_C3_relative_terminal_zero_identity",
    "proved_C3_third_zero_identity",
    "proved_original_compact_relative_C5_terminal_identity",
    "proved_original_n1_C3_transported_patch_integral_identity",
    "proved_original_quiet_linear_mass_source_identity",
    "compact_raw_beta",
    "compact_raw_beta_derivative",
}
_SOURCE_CONSTANT_OPS = {
    "exact_positive_integer_expression",
    "mathematical_pi",
    "mathematical_positive_infinity",
    "exact_original_source_heat_phi",
}
_SPECIAL_OPS = {
    "original_monotone_phase_inverse",
    "original_lazy_flat_branch",
    "exact_original_source_flat_sigma",
    "substitute_original_inverse_angle",
    "substitute_original_native_coordinate",
    "substitute_free_query_at_native_coordinate",
    "inverse_substitution",
    "native_substitution",
    "definite_integral",
    "positive_quotient",
}
_HANDLE_FIELDS = (
    "argument",
    "operand",
    "base",
    "body",
    "active_body",
    "flat_predicate",
    "phase_function",
    "target_phase",
    "inverse_angle",
    "query_coordinate",
    "free_query",
    "native_coordinate",
    "coordinate",
    "left",
    "right",
    "expression",
    "variable",
    "value",
    "integrand",
    "lower",
    "upper",
    "numerator",
    "denominator",
    "arguments",
    "operands",
)


class _Unresolved(ValueError):
    """A reachable source operation outside the bounded reducer."""


def _json_safe(value: Any) -> Any:
    if isinstance(value, Fraction):
        return {"__fraction__": [value.numerator, value.denominator]}
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [_json_safe(item) for item in value]
    return value


def canonical_digest(value: Any) -> str:
    """Return the stable source digest used by the consumer receipts."""

    try:
        raw = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError):
        raw = json.dumps(
            _json_safe(value),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _module_sha256() -> str:
    digest = hashlib.sha256()
    with Path(__file__).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _require_digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise ValueError(f"{label} must be a SHA-256 hex digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ValueError(f"{label} must be a SHA-256 hex digest") from exc
    return value.lower()


def _handle(value: Any, label: str, count: int) -> int:
    if isinstance(value, Mapping):
        value = value.get("node", value.get("handle", value.get("root")))
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer graph handle")
    if value < 0 or value >= count:
        raise ValueError(f"{label} is outside the source graph")
    return value


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{label} is not an exact rational")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, str):
        text = value.strip()
        try:
            pieces = text.split("/")
            if len(pieces) == 2:
                return Fraction(int(pieces[0]), int(pieces[1]))
            if len(pieces) == 1:
                return Fraction(int(pieces[0]), 1)
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{label} is not an exact rational") from exc
    if isinstance(value, Mapping):
        if "numerator" in value or "denominator" in value:
            return _fraction(value.get("numerator"), label) / _fraction(
                value.get("denominator", 1), label
            )
    raise ValueError(f"{label} is not an exact rational")


def _operation(row: Mapping[str, Any], index: int) -> str:
    operation = row.get("operation", row.get("op"))
    if not isinstance(operation, str) or not operation:
        raise _Unresolved(f"node {index} has no operation")
    return operation


def _field(row: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in row:
            return row[name]
    raise _Unresolved("missing " + "/".join(names))


def _args(row: Mapping[str, Any], index: int) -> list[Any]:
    raw = row.get("arguments", row.get("operands"))
    if not _is_sequence(raw):
        raise _Unresolved(f"node {index}.arguments is not a sequence")
    return list(raw)


def _exact_value(nodes: Sequence[Any], index: int) -> Fraction | None:
    try:
        row = nodes[index]
        if not isinstance(row, Mapping) or _operation(row, index) != "exact_rational":
            return None
        numerator = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
        denominator = _fraction(row.get("denominator", 1), f"node {index}.denominator")
        if denominator == 0:
            return None
        return numerator / denominator
    except (IndexError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return None


def _guarded(row: Mapping[str, Any]) -> bool:
    return any(row.get(name) is True for name in _FUNCTION_GUARDS)


def _simplify(value: Any) -> Any:
    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact source congruence")
    try:
        return _sp.cancel(_sp.together(value))
    except Exception:
        return _sp.simplify(value)


def _record(value: Any) -> dict[str, Any]:
    value = _simplify(value)
    return {
        "srepr": _sp.srepr(value),
        "text": str(value),
        "is_zero": bool(value == 0),
    }


def _expr_key(value: Any) -> str:
    return _sp.srepr(_simplify(value))


def _atom(label: str) -> Any:
    return _sp.Symbol(label, commutative=True)


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[str] = set()
    for value in values:
        marker = json.dumps(_json_safe(value), sort_keys=True, separators=(",", ":"), default=str)
        if marker not in seen:
            seen.add(marker)
            result.append(deepcopy(value))
    return result


def _proof_refs(value: Any) -> list[Any]:
    if value in (None, ""):
        return []
    if isinstance(value, Mapping):
        for key in ("proof_refs", "source_refs", "refs", "proof_ref", "source_ref"):
            if key in value:
                return _proof_refs(value[key])
        return [deepcopy(dict(value))]
    if _is_sequence(value):
        result: list[Any] = []
        for item in value:
            result.extend(_proof_refs(item))
        return result
    if isinstance(value, (str, int)):
        return [deepcopy(value)]
    raise ValueError("proof reference must be JSON-compatible data")


def _unwrap_chart(value: Mapping[str, Any]) -> tuple[Mapping[str, Any], list[Mapping[str, Any]]]:
    owners = [value]
    current = value
    for _ in range(5):
        if isinstance(current.get("chart"), Mapping):
            current = current["chart"]
            owners.append(current)
            continue
        if isinstance(current.get("bindings"), Mapping):
            current = current["bindings"]
            owners.append(current)
            continue
        if isinstance(current.get("charts"), Mapping):
            charts = current["charts"]
            for candidate in ("postrestore", "patch", "finite", "modified"):
                if isinstance(charts.get(candidate), Mapping):
                    current = charts[candidate]
                    owners.append(current)
                    break
            else:
                break
            continue
        break
    return current, owners


def _lookup_mapping(owners: Sequence[Mapping[str, Any]], *keys: str) -> Any:
    for owner in owners:
        for key in keys:
            value = owner.get(key)
            if value is not None:
                return value
    return None


def _rows(raw: Any, label: str, count: int) -> dict[int, int]:
    if isinstance(raw, Mapping):
        for wrapper in ("rows", "handles", "nodes", "values"):
            if wrapper in raw and isinstance(raw[wrapper], (Mapping, Sequence)) and not isinstance(raw[wrapper], (str, bytes, bytearray)):
                return _rows(raw[wrapper], label, count)
        result: dict[int, int] = {}
        for order in ORDERS:
            candidate = raw.get(order, raw.get(str(order)))
            if candidate is None:
                candidate = raw.get(f"Z{order}", raw.get(f"ordinary_Z{order}"))
            if candidate is None:
                raise ValueError(f"{label} lacks ordinary order {order}")
            result[order] = _handle(candidate, f"{label}[{order}]", count)
        return result
    if not _is_sequence(raw):
        raise ValueError(f"{label} must be a row sequence or mapping")
    values = list(raw)
    result = {}
    for order in ORDERS:
        if order >= len(values):
            raise ValueError(f"{label} lacks ordinary order {order}")
        result[order] = _handle(values[order], f"{label}[{order}]", count)
    return result


def _root_rows(chart: Mapping[str, Any], owners: Sequence[Mapping[str, Any]], count: int) -> dict[str, dict[int, int]]:
    roots = _lookup_mapping(owners, "roots", "generic_roots", "leading_roots")
    if not isinstance(roots, Mapping):
        raise ValueError("chart is missing explicit leading generic roots")
    out: dict[str, dict[int, int]] = {}
    for name in GENERIC_ROOTS:
        raw = roots.get(name)
        if raw is None:
            raw = roots.get(name.lower())
        if raw is None:
            raise ValueError(f"chart roots lacks {name}")
        out[name] = _rows(raw, f"chart.roots.{name}", count)
    return out


def _target_rows(chart: Mapping[str, Any], owners: Sequence[Mapping[str, Any]], count: int) -> dict[str, dict[int, int]]:
    loop = _lookup_mapping(owners, "loop", "loop_rows", "leading_loop")
    if not isinstance(loop, Mapping):
        raise ValueError("chart is missing explicit loop bindings")
    finite = _lookup_mapping(owners, "finite_N", "finite", "finite_rows")
    if not isinstance(finite, Mapping):
        raise ValueError("chart is missing explicit finite_N bindings")
    density = finite.get("five_signed_density_increments", finite.get("density_increments"))
    if not isinstance(density, Mapping):
        raise ValueError("chart is missing five_signed_density_increments")
    out: dict[str, dict[int, int]] = {}
    for name in LOOP_FIELDS:
        if name not in loop:
            raise ValueError(f"chart.loop lacks {name}")
        out[name] = _rows(loop[name], f"chart.loop.{name}", count)
    for name in FINITE_FIELDS:
        raw = finite.get(name)
        if raw is None:
            raise ValueError(f"chart.finite_N lacks exact {name} rows")
        out[name] = _rows(raw, f"chart.finite_N.{name}", count)
    for name in DENSITY_FIELDS:
        if name not in density:
            raise ValueError(f"chart.finite_N.five_signed_density_increments lacks {name}")
        out[name] = _rows(density[name], f"chart.finite_N.five_signed_density_increments.{name}", count)
    return out


def _parse_side_key(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    lowered = value.lower().replace("-", "_")
    for side, names in _SIDE_ALIASES.items():
        if lowered in names:
            return side
    return None


def _phase_side_values(raw: Any, count: int, label: str) -> dict[str, int]:
    if isinstance(raw, Mapping):
        result: dict[str, int] = {}
        for side in SIDES:
            candidate = None
            for key in _SIDE_ALIASES[side]:
                if key in raw:
                    candidate = raw[key]
                    break
            if candidate is None:
                raise ValueError(f"{label} lacks {side}")
            result[side] = _handle(candidate, f"{label}.{side}", count)
        return result
    if _is_sequence(raw) and len(raw) == 2:
        return {
            "postrestore": _handle(raw[0], f"{label}[0]", count),
            "patch": _handle(raw[1], f"{label}[1]", count),
        }
    raise ValueError(f"{label} must provide postrestore and patch handles")


def _validate_phase(nodes: Sequence[Mapping[str, Any]], proof: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(proof, Mapping):
        raise ValueError("phase_proof must be a mapping")
    # The phase receipt is an admission token for this exact saved loop
    # source.  These checks are deliberately source/static checks: the
    # consumer never reconstructs a phase value or differentiates the
    # discontinuous representative.
    loop_source_sha256 = proof.get("same_loop_source_sha256")
    if loop_source_sha256 != "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d":
        raise ValueError("phase_proof is bound to a different loop source")
    if proof.get("source_only") is not True:
        raise ValueError("phase_proof.source_only must be true")
    if proof.get("actual_phase_Z_independence_proved") is not True:
        raise ValueError("phase_proof.actual_phase_Z_independence_proved must be true")
    if proof.get("fractional_part_differentiated") is not False:
        raise ValueError("phase_proof.fractional_part_differentiated must be false")
    _require_digest(proof.get("radius_lemma_digest"), "phase_proof.radius_lemma_digest")
    _require_digest(proof.get("saved_lift_certificate_digest"), "phase_proof.saved_lift_certificate_digest")
    if proof.get("actual_common_phase_proved") is not True:
        raise ValueError("phase_proof.actual_common_phase_proved must be true")
    if proof.get("smooth_periodic_lift_bound") is not True:
        raise ValueError("phase_proof.smooth_periodic_lift_bound must be true")
    count = len(nodes)
    source_graph_digest = _require_digest(proof.get("source_graph_digest"), "phase_proof.source_graph_digest")
    if source_graph_digest != canonical_digest(list(nodes)):
        raise ValueError("phase_proof.source_graph_digest does not match the supplied source graph")
    descriptors = proof.get("exact_source_descriptors")
    if not isinstance(descriptors, Mapping) or not descriptors:
        raise ValueError("phase_proof.exact_source_descriptors is required")
    descriptor_rows: dict[str, Any] = {}
    for raw_index, descriptor in descriptors.items():
        try:
            index = int(raw_index)
        except (TypeError, ValueError) as exc:
            raise ValueError("phase_proof exact source descriptor has a noninteger handle") from exc
        if isinstance(raw_index, bool) or index < 0 or index >= count:
            raise ValueError("phase_proof exact source descriptor is outside the source graph")
        if str(index) in descriptor_rows:
            raise ValueError(f"phase_proof exact source descriptor duplicates node {index}")
        if not isinstance(descriptor, Mapping) or dict(nodes[index]) != dict(descriptor):
            raise ValueError(f"phase source descriptor changed at node {index}")
        descriptor_rows[str(index)] = deepcopy(dict(descriptor))
    # The phase constructor reports the complete closure used by its proof.
    # Preserve the optional explicit dependency list and reject foreign
    # handles when it is present.
    dependency_nodes = proof.get("Z_independent_dependency_nodes")
    if dependency_nodes is not None:
        if not _is_sequence(dependency_nodes):
            raise ValueError("phase_proof.Z_independent_dependency_nodes must be a sequence")
        dependency_nodes = sorted({_handle(item, "phase_proof.Z_independent_dependency_nodes", count) for item in dependency_nodes})
        missing_dependencies = [node for node in dependency_nodes if str(node) not in descriptor_rows]
        if missing_dependencies:
            raise ValueError(
                "phase_proof exact source descriptors omit declared phase dependencies: "
                + ",".join(str(node) for node in missing_dependencies)
            )
    endpoint_reduction = proof.get("endpoint_radius_reduction")
    if not isinstance(endpoint_reduction, Mapping):
        raise ValueError("phase_proof.endpoint_radius_reduction is required")
    for side in SIDES:
        side_reduction = next(
            (endpoint_reduction[key] for key in _SIDE_ALIASES[side] if key in endpoint_reduction),
            None,
        )
        if not isinstance(side_reduction, Mapping):
            raise ValueError(f"phase_proof.endpoint_radius_reduction lacks {side}")
        if "closed_node" not in side_reduction:
            raise ValueError(f"phase_proof.endpoint_radius_reduction.{side} lacks closed_node")
        _handle(side_reduction["closed_node"], f"phase_proof.endpoint_radius_reduction.{side}.closed_node", count)
        if "body_node" in side_reduction:
            _handle(side_reduction["body_node"], f"phase_proof.endpoint_radius_reduction.{side}.body_node", count)
    ordinary_orders = proof.get("ordinary_orders", [0, 1])
    if list(ordinary_orders) != [0, 1]:
        raise ValueError("phase_proof must cover ordinary Z orders [0, 1]")
    phase_nodes = proof.get("phase_nodes", proof.get("common_phase_nodes"))
    coords = proof.get("coordinates", proof.get("native_coordinates"))
    endpoints = proof.get("endpoints", proof.get("phase_endpoints"))
    if phase_nodes is None:
        phase_nodes = {
            "postrestore": proof.get("postrestore_phase_node"),
            "patch": proof.get("patch_phase_node"),
        }
    if coords is None:
        coords = {
            "postrestore": proof.get("postrestore_coordinate_node"),
            "patch": proof.get("patch_coordinate_node"),
        }
    phase = _phase_side_values(phase_nodes, count, "phase_proof.phase_nodes")
    coordinate = _phase_side_values(coords, count, "phase_proof.coordinates")
    if endpoints is None:
        endpoints = (1, 1)
    if isinstance(endpoints, Mapping):
        endpoint_values = []
        for side in SIDES:
            candidate = next((endpoints.get(key) for key in _SIDE_ALIASES[side] if key in endpoints), None)
            endpoint_values.append(candidate)
    else:
        endpoint_values = list(endpoints) if _is_sequence(endpoints) else []
    if endpoint_values != [1, 1]:
        raise ValueError("phase_proof endpoints must be exact native endpoints [1, 1]")
    kept = {
        "actual_common_phase_proved": True,
        "smooth_periodic_lift_bound": True,
        "source_only": True,
        "same_loop_source_sha256": loop_source_sha256,
        "actual_phase_Z_independence_proved": True,
        "fractional_part_differentiated": False,
        "radius_lemma_digest": proof["radius_lemma_digest"].lower(),
        "saved_lift_certificate_digest": proof["saved_lift_certificate_digest"].lower(),
        "source_graph_digest": source_graph_digest,
        "exact_source_descriptors": descriptor_rows,
        "Z_independent_dependency_nodes": dependency_nodes,
        "endpoint_radius_reduction": deepcopy(endpoint_reduction),
        "phase_nodes": phase,
        "coordinates": coordinate,
        "endpoints": [1, 1],
        "proof_refs": deepcopy(proof.get("proof_refs", proof.get("source_refs", proof.get("proof_ref")))),
    }
    return kept


def _native_port(
    chart: Mapping[str, Any],
    owners: Sequence[Mapping[str, Any]],
    side: str,
    phase_coordinate: int,
    count: int,
) -> dict[str, Any]:
    raw = _lookup_mapping(
        owners,
        "native_coordinate_to_one",
        "native_endpoint",
        "native_endpoint_map",
        "endpoint_substitution",
        "native_port",
    )
    coordinate: Any = None
    endpoint: Any = None
    if isinstance(raw, Mapping):
        side_raw = raw
        for key in _SIDE_ALIASES[side]:
            if isinstance(raw.get(key), Mapping):
                side_raw = raw[key]
                break
        coordinate = side_raw.get("coordinate_node", side_raw.get("native_coordinate", side_raw.get("node")))
        endpoint = side_raw.get("endpoint", side_raw.get("value", 1))
    if coordinate is None:
        coordinate = _lookup_mapping(owners, "native_coordinate_node", "native_coordinate", "coordinate_node", "coordinate")
    if coordinate is None:
        coordinate = phase_coordinate
    coordinate = _handle(coordinate, f"{side}.native_coordinate", count)
    if endpoint is None:
        endpoint = _lookup_mapping(owners, "native_endpoint_value", "native_endpoint_value_one")
    if endpoint is None:
        endpoint = 1
    if endpoint != 1:
        raise ValueError(f"{side}.native endpoint must be exact one")
    return {
        "coordinate_node": coordinate,
        "endpoint": 1,
        "source": "explicit native coordinate-to-one port or admitted common-phase coordinate",
    }


def _generic_rows(packet: Mapping[str, Any], count: int) -> tuple[dict[tuple[str, int], dict[str, Any]], dict[int, Any], list[Any]]:
    if not isinstance(packet, Mapping):
        raise ValueError("generic_packet must be a mapping")
    if packet.get("passed") is not True:
        raise ValueError("generic_packet.passed must be true")
    rows = packet.get("root_rows", packet.get("rows"))
    if not _is_sequence(rows):
        raise ValueError("generic_packet.root_rows must be an explicit sequence")
    expected = {(root, order) for root in GENERIC_ROOTS for order in ORDERS}
    seen: set[tuple[str, int]] = set()
    accepted: dict[tuple[str, int], dict[str, Any]] = {}
    aliases: dict[int, Any] = {}
    refs: list[Any] = []
    for offset, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"generic_packet.root_rows[{offset}] must be a mapping")
        root = row.get("root", row.get("name"))
        order = row.get("order", row.get("ordinary_Z_order"))
        if root not in GENERIC_ROOTS or isinstance(order, bool) or order not in ORDERS:
            raise ValueError(f"generic_packet.root_rows[{offset}] has an invalid root/order")
        key = (root, order)
        if key in seen:
            raise ValueError(f"duplicate generic root row {root}[{order}]")
        seen.add(key)
        if row.get("status") != "proved":
            raise ValueError(f"generic root row {root}[{order}] is not proved")
        side_rows: dict[str, Any] = {}
        for side in SIDES:
            value = row.get(side)
            if not isinstance(value, Mapping):
                raise ValueError(f"generic root row {root}[{order}] lacks {side} source row")
            if value.get("root_node", value.get("node")) is None:
                raise ValueError(f"generic root row {root}[{order}] lacks {side}.root_node")
            node = _handle(value.get("root_node", value.get("node")), f"generic {side}.{root}[{order}].root_node", count)
            saved = value.get("value")
            if not isinstance(saved, Mapping) or not isinstance(saved.get("srepr"), str):
                raise ValueError(f"generic root row {root}[{order}] lacks {side}.value.srepr")
            if not isinstance(saved.get("is_zero"), bool):
                raise ValueError(f"generic root row {root}[{order}] lacks boolean {side}.value.is_zero")
            side_rows[side] = {"node": node, "value": deepcopy(dict(saved))}
            refs.extend(_proof_refs(value.get("proof_refs", value.get("proof_ref"))))
        left_srepr = side_rows["postrestore"]["value"]["srepr"]
        right_srepr = side_rows["patch"]["value"]["srepr"]
        if left_srepr != right_srepr:
            raise ValueError(f"generic root row {root}[{order}] has unequal admitted canonical values")
        left_zero = side_rows["postrestore"]["value"]["is_zero"]
        right_zero = side_rows["patch"]["value"]["is_zero"]
        if left_zero != right_zero:
            raise ValueError(f"generic root row {root}[{order}] has unequal zero admission")
        atom = _sp.Integer(0) if left_zero else _atom("GENERIC_VALUE:" + hashlib.sha256(left_srepr.encode()).hexdigest())
        for side in SIDES:
            node = side_rows[side]["node"]
            previous = aliases.get(node)
            if previous is not None and _expr_key(previous) != _expr_key(atom):
                raise ValueError(f"generic source handle {node} has two admitted canonical values")
            aliases[node] = atom
        accepted[key] = {"root": root, "order": order, "sides": side_rows, "canonical_srepr": left_srepr, "is_zero": left_zero}
    if seen != expected or len(rows) != len(expected):
        raise ValueError("generic_packet.root_rows must contain exactly 18 explicit proved rows")
    return accepted, aliases, refs


def _explicit_aliases(chart: Mapping[str, Any], owners: Sequence[Mapping[str, Any]], side: str, count: int) -> list[dict[str, Any]]:
    raw = _lookup_mapping(owners, "root_aliases", "generic_root_aliases", "leading_root_aliases", "aliases")
    if raw is None:
        return []
    if isinstance(raw, Mapping):
        for key in _SIDE_ALIASES[side]:
            if isinstance(raw.get(key), (Mapping, Sequence)) and not isinstance(raw.get(key), (str, bytes, bytearray)):
                raw = raw[key]
                break
    entries: list[Any]
    if isinstance(raw, Mapping):
        entries = []
        for handle, atom in raw.items():
            entries.append({"node": handle, "canonical_atom": atom})
    elif _is_sequence(raw):
        entries = list(raw)
    else:
        raise ValueError(f"{side} root aliases must be a mapping or sequence")
    parsed: list[dict[str, Any]] = []
    for offset, item in enumerate(entries):
        if not isinstance(item, Mapping):
            raise ValueError(f"{side} root alias {offset} must be a mapping")
        node = _handle(item.get("node", item.get("handle", item.get("root_node"))), f"{side}.root_aliases[{offset}].node", count)
        parsed.append({"node": node, "canonical_atom": deepcopy(item.get("canonical_atom", item.get("atom", item.get("value"))))})
    return parsed


def _alias_marker(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, Mapping):
        return value.get("srepr", value.get("text", value.get("atom")))
    return str(value)


def _validate_explicit_aliases(
    chart: Mapping[str, Any],
    owners: Sequence[Mapping[str, Any]],
    side: str,
    aliases: Mapping[int, Any],
    count: int,
) -> list[dict[str, Any]]:
    explicit = _explicit_aliases(chart, owners, side, count)
    kept: list[dict[str, Any]] = []
    for item in explicit:
        node = item["node"]
        if node not in aliases:
            raise ValueError(f"{side}.root_aliases references a handle outside admitted generic root rows")
        marker = _alias_marker(item.get("canonical_atom"))
        if marker is not None:
            expected = aliases[node]
            expected_markers = {str(expected), _sp.srepr(expected)}
            if expected == 0:
                expected_markers.update(("0", "Integer(0)"))
            if marker not in expected_markers:
                raise ValueError(f"{side}.root_aliases canonical atom does not match the admitted generic value")
        kept.append({"node": node, "canonical_atom": marker, "validated_against_generic_packet": True})
    return kept


class _Evaluator:
    """Exact source reducer for the bounded reachable operation set."""

    def __init__(self, nodes: Sequence[Mapping[str, Any]], *, side: str, aliases: Mapping[int, Any], native_coordinate: int) -> None:
        self.nodes = nodes
        self.side = side
        self.aliases = dict(aliases)
        self.native_coordinate = native_coordinate
        self.memo: dict[int, Any] = {}
        self.used_nodes: set[int] = set()
        self.unresolved: list[dict[str, Any]] = []
        self.positive_quotients: list[dict[str, Any]] = []
        self.substitutions: list[dict[str, Any]] = []
        self.special_nodes: list[dict[str, Any]] = []
        self.endpoint_substitutions: list[dict[str, Any]] = []
        self._binder_depth = 0

    def _row(self, index: int) -> Mapping[str, Any]:
        if isinstance(index, bool) or not isinstance(index, int) or index < 0 or index >= len(self.nodes):
            raise _Unresolved(f"invalid source handle {index!r}")
        row = self.nodes[index]
        if not isinstance(row, Mapping):
            raise _Unresolved(f"node {index} is not a mapping")
        return row

    def _child(self, index: int, value: Any, label: str) -> int:
        child = _handle(value, f"node {index}.{label}", len(self.nodes))
        if index >= 0 and child >= index:
            raise _Unresolved(f"forward source reference from node {index}.{label} to node {child}")
        return child

    def _remember_special(self, index: int, row: Mapping[str, Any], operation: str) -> None:
        if operation in _SPECIAL_OPS:
            self.special_nodes.append({"node": index, "operation": operation, "descriptor": deepcopy(dict(row))})

    def _compound(self, operation: str, arguments: Sequence[Any], metadata: Mapping[str, Any] | None = None, *, index: int | None = None) -> Any:
        payload = {
            "operation": operation,
            "arguments": [_expr_key(value) for value in arguments],
            "metadata": _json_safe(metadata or {}),
        }
        marker = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()
        return _atom(f"COMPOUND:{operation}:{marker}")

    def _canonical_metadata(self, row: Mapping[str, Any], excluded: Sequence[str] = ()) -> dict[str, Any]:
        ignored = set(excluded) | {"operation", "op", "node", "handle", "root"}
        result: dict[str, Any] = {}
        for key, value in row.items():
            if key in ignored:
                continue
            if key == "angle_variable" or key.endswith("_angle_variable"):
                result[key] = "ALPHA_BINDER"
            else:
                result[key] = deepcopy(value)
        return result

    def _ref(self, row: Mapping[str, Any], *names: str) -> int:
        for name in names:
            if name in row:
                return self._child(-1, row[name], name)
        raise _Unresolved("missing " + "/".join(names))

    def _bound(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Any]) -> Any:
        if index in env:
            return env[index]
        name = row.get("name", row.get("variable"))
        if name in env:
            return env[name]
        if index == self.native_coordinate:
            self.used_nodes.add(index)
            self.endpoint_substitutions.append({"node": index, "endpoint": 1, "reason": "native coordinate to one"})
            return _sp.Integer(1)
        source_name = str(name) if isinstance(name, str) else f"node:{index}"
        return _atom(f"SOURCE_BOUND:node:{index}:" + source_name)

    def _aliases_allowed(self, env: Mapping[Any, Any]) -> bool:
        """Root/phase aliases are endpoint ports, never altered native fsubs."""

        if self.native_coordinate in env:
            try:
                return bool(_simplify(env[self.native_coordinate] - 1) == 0)
            except Exception:
                return False
        native_row = self._row(self.native_coordinate)
        native_name = native_row.get("name")
        if isinstance(native_name, str) and native_name in env:
            try:
                return bool(_simplify(env[native_name] - 1) == 0)
            except Exception:
                return False
        return True

    def _analytic(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Any]) -> Any:
        name = row.get("name", row.get("function", row.get("unary", row.get("kind"))))
        if not isinstance(name, str):
            raise _Unresolved(f"node {index} analytic_unary lacks a function name")
        child = self._child(index, _field(row, "argument", "operand", "value"), "argument")
        argument = self.rewrite(child, env)
        if name in _ANALYTIC:
            if name == "exp":
                return _simplify(_sp.exp(argument))
            if name == "expm1":
                return _simplify(_sp.exp(argument) - 1)
            if name == "log":
                return _simplify(_sp.log(argument))
            if name == "log1p":
                return _simplify(_sp.log(1 + argument))
            if name in {"sqrt", "positive_sqrt"}:
                return _simplify(_sp.sqrt(argument))
            if name == "sin":
                return _simplify(_sp.sin(argument))
            if name == "cos":
                return _simplify(_sp.cos(argument))
        if name in _OPAQUE_ANALYTIC:
            self._remember_special(index, row, "analytic_unary")
            return self._compound("analytic_unary:" + name, [argument], self._canonical_metadata(row, ("argument", "operand", "value")), index=index)
        raise _Unresolved(f"unsupported analytic unary {name!r}")

    def _positive_quotient(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Any]) -> Any:
        numerator = self._child(index, _field(row, "numerator"), "numerator")
        denominator = self._child(index, _field(row, "denominator"), "denominator")
        certificate = row.get("source_positive_certificate", row.get("positive_certificate", row.get("denominator_positive")))
        if certificate in (None, False, ""):
            raise _Unresolved(f"node {index} positive_quotient lacks a positive source guard")
        den = self.rewrite(denominator, env)
        if _simplify(den) == 0:
            raise _Unresolved(f"node {index} positive_quotient denominator is symbolically zero")
        num = self.rewrite(numerator, env)
        self.positive_quotients.append({"node": index, "numerator": numerator, "denominator": denominator, "certificate": deepcopy(certificate)})
        return _simplify(num / den)

    def _function_substitution(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Any]) -> Any:
        if not _guarded(row):
            raise _Unresolved(f"node {index} function_substitution lacks a guard")
        expression = self._child(index, _field(row, "expression", "body"), "expression")
        variable = self._child(index, _field(row, "variable", "native_variable"), "variable")
        value = self._child(index, _field(row, "value", "endpoint", "query"), "value")
        replacement = self.rewrite(value, env)
        scoped = dict(env)
        scoped[variable] = replacement
        variable_row = self._row(variable)
        if isinstance(variable_row.get("name"), str):
            scoped[variable_row["name"]] = replacement
        self.substitutions.append({"node": index, "expression": expression, "variable": variable, "value": value, "guarded": True, "endpoint": _record(replacement), "source_ref": deepcopy(row.get("source_ref", row.get("proof_ref")))})
        return self.rewrite(expression, scoped)

    def _special_substitution(self, index: int, row: Mapping[str, Any], operation: str, env: Mapping[Any, Any]) -> Any:
        if operation in {"substitute_original_inverse_angle", "inverse_substitution"}:
            if not isinstance(row.get("angle_variable"), str):
                raise _Unresolved(f"node {index} inverse substitution lacks angle_variable")
            body = self._child(index, _field(row, "body", "expression"), "body")
            inverse = self._child(index, _field(row, "inverse_angle", "inverse", "value"), "inverse_angle")
            # The angle is a lexical binder for the body.  Its spelling is
            # chart-local and must never enter the canonical key.
            scoped = dict(env)
            binder = _atom(f"BINDER:ANGLE:{self._binder_depth}")
            scoped[row["angle_variable"]] = binder
            if isinstance(row.get("angle_variable_node"), int) and not isinstance(row.get("angle_variable_node"), bool):
                scoped[row["angle_variable_node"]] = binder
            self._binder_depth += 1
            try:
                body_expr = self.rewrite(body, scoped)
            finally:
                self._binder_depth -= 1
            inverse_expr = self.rewrite(inverse, env)
            excluded = ["body", "expression", "inverse_angle", "inverse", "value", "angle_variable"]
            if isinstance(row.get("angle_variable_node"), int) and not isinstance(row.get("angle_variable_node"), bool):
                self._child(index, row["angle_variable_node"], "angle_variable_node")
                excluded.append("angle_variable_node")
            metadata = self._canonical_metadata(row, excluded)
            metadata["angle_variable"] = "ALPHA_BINDER"
            return self._compound(operation, (body_expr, inverse_expr), metadata, index=index)
        if operation in {"substitute_original_native_coordinate", "substitute_free_query_at_native_coordinate", "native_substitution"}:
            body = self._child(index, _field(row, "body", "expression"), "body")
            query = self._child(index, _field(row, "query_coordinate", "free_query", "value", "coordinate", "native_coordinate"), "query_coordinate")
            query_expr = self.rewrite(query, env)
            variable_raw = row.get("native_variable", row.get("variable", row.get("coordinate")))
            scoped = dict(env)
            if isinstance(variable_raw, int) and not isinstance(variable_raw, bool):
                variable = self._child(index, variable_raw, "native_variable")
                scoped[variable] = query_expr
                variable_row = self._row(variable)
                if isinstance(variable_row.get("name"), str):
                    scoped[variable_row["name"]] = query_expr
            elif isinstance(variable_raw, str):
                scoped[variable_raw] = query_expr
            else:
                raise _Unresolved(f"node {index} native substitution lacks a variable binder")
            body_expr = self.rewrite(body, scoped)
            metadata = self._canonical_metadata(row, ("body", "expression", "query_coordinate", "free_query", "value", "coordinate", "native_coordinate", "native_variable", "variable"))
            if "native_variable" in metadata:
                metadata["native_variable"] = "NATIVE_COORDINATE"
            if "variable" in metadata:
                metadata["variable"] = "NATIVE_COORDINATE"
            return self._compound(operation, (body_expr, query_expr), metadata, index=index)
        raise _Unresolved(f"unsupported substitution operation {operation!r}")

    def _phase_inverse(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Any]) -> Any:
        required = (
            "angle_variable",
            "derivative_positive_certificate",
            "flat_inverse",
            "phase_endpoints",
        )
        if not isinstance(row.get("angle_variable"), str) or any(row.get(key) in (None, "", False) for key in required[1:]):
            raise _Unresolved(f"node {index} original_monotone_phase_inverse lacks branch or positivity inputs")
        if row.get("branch_inverse_is_original_not_support_cap") is not True:
            raise _Unresolved(f"node {index} original inverse branch is not admitted")
        phase_function = self._child(index, _field(row, "phase_function"), "phase_function")
        target_phase = self._child(index, _field(row, "target_phase"), "target_phase")
        inverse_angle = _atom(f"BINDER:ANGLE:{self._binder_depth}")
        scoped = dict(env)
        scoped[row["angle_variable"]] = inverse_angle
        if isinstance(row.get("angle_variable_node"), int) and not isinstance(row.get("angle_variable_node"), bool):
            scoped[row["angle_variable_node"]] = inverse_angle
        self._binder_depth += 1
        try:
            phase_expr = self.rewrite(phase_function, scoped)
        finally:
            self._binder_depth -= 1
        target_expr = self.rewrite(target_phase, env)
        arguments = [phase_expr, target_expr]
        excluded = ["phase_function", "target_phase", "angle_variable"]
        if isinstance(row.get("angle_variable_node"), int) and not isinstance(row.get("angle_variable_node"), bool):
            self._child(index, row["angle_variable_node"], "angle_variable_node")
            excluded.append("angle_variable_node")
        # In the saved descriptor Delta, eta, and angle limits are source
        # handles, despite being serialized as integers.  Canonicalize them
        # as children so different graph offsets cannot create a false seam.
        for key in ("Delta", "delta", "eta", "eta21", "angle_lower", "angle_upper"):
            if key not in row or not isinstance(row[key], int) or isinstance(row[key], bool):
                continue
            child = self._child(index, row[key], key)
            arguments.append(self.rewrite(child, env))
            excluded.append(key)
        metadata = self._canonical_metadata(row, excluded)
        metadata["angle_variable"] = "ALPHA_BINDER"
        metadata["branch_inverse_is_original_not_support_cap"] = True
        metadata["derivative_positive_certificate"] = deepcopy(row["derivative_positive_certificate"])
        self._remember_special(index, row, "original_monotone_phase_inverse")
        return self._compound("original_monotone_phase_inverse", arguments, metadata, index=index)

    def _lazy_flat(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Any]) -> Any:
        if row.get("active_body_not_evaluated_on_flat_branch") is not True:
            raise _Unresolved(f"node {index} flat branch does not preserve lazy active body")
        predicate = self._child(index, _field(row, "flat_predicate"), "flat_predicate")
        active = self._child(index, _field(row, "active_body", "body"), "active_body")
        predicate_expr = self.rewrite(predicate, env)
        active_expr = self.rewrite(active, env)
        arguments = [predicate_expr, active_expr]
        excluded = ["flat_predicate", "active_body", "body"]
        for key in (
            "Delta",
            "delta",
            "eta",
            "eta21",
            "original_collar_coordinate",
            "original_collar_upper",
            "flat_value",
        ):
            if key not in row or not isinstance(row[key], int) or isinstance(row[key], bool):
                continue
            child = self._child(index, row[key], key)
            arguments.append(self.rewrite(child, env))
            excluded.append(key)
        metadata = self._canonical_metadata(row, excluded)
        metadata["flat_predicate_preserved"] = True
        metadata["active_body_not_evaluated_on_flat_branch"] = True
        self._remember_special(index, row, "original_lazy_flat_branch")
        return self._compound("original_lazy_flat_branch", arguments, metadata, index=index)

    def _integral(self, index: int, row: Mapping[str, Any], env: Mapping[Any, Any]) -> Any:
        integrand = self._child(index, _field(row, "integrand"), "integrand")
        lower = self._child(index, _field(row, "lower"), "lower")
        upper = self._child(index, _field(row, "upper"), "upper")
        lower_expr = self.rewrite(lower, env)
        upper_expr = self.rewrite(upper, env)
        variable = row.get("variable")
        if not isinstance(variable, (str, int)):
            raise _Unresolved(f"node {index} definite_integral lacks a lexical variable")
        binder = _atom(f"BINDER:ANGLE:{self._binder_depth}")
        scoped = dict(env)
        if isinstance(variable, int):
            scoped[variable] = binder
            binder_row = self._row(variable)
            if isinstance(binder_row.get("name"), str):
                scoped[binder_row["name"]] = binder
        else:
            scoped[variable] = binder
        self._binder_depth += 1
        try:
            integrand_expr = self.rewrite(integrand, scoped)
        finally:
            self._binder_depth -= 1
        metadata = self._canonical_metadata(row, ("integrand", "lower", "upper", "variable"))
        metadata["variable"] = "ALPHA_ANGLE_BINDER"
        metadata["lexical_alpha_renamed"] = True
        self._remember_special(index, row, "definite_integral")
        return self._compound("definite_integral", (integrand_expr, lower_expr, upper_expr), metadata, index=index)

    def _generic_compound(self, index: int, row: Mapping[str, Any], operation: str, env: Mapping[Any, Any]) -> Any:
        child_values: list[Any] = []
        excluded: list[str] = []
        for key in _HANDLE_FIELDS:
            if key not in row:
                continue
            value = row[key]
            if key in {"arguments", "operands"} and _is_sequence(value):
                for offset, child in enumerate(value):
                    child_index = self._child(index, child, f"{key}[{offset}]")
                    child_values.append(self.rewrite(child_index, env))
                excluded.append(key)
            elif isinstance(value, int) and not isinstance(value, bool):
                child_index = self._child(index, value, key)
                child_values.append(self.rewrite(child_index, env))
                excluded.append(key)
        metadata = self._canonical_metadata(row, excluded)
        if operation not in {"exact_real_comparison", "logical_or"}:
            # Implicit source functions and source parameters are opaque.
            # Matching names/descriptors do not identify distinct handles.
            metadata["source_node"] = index
        return self._compound(operation, child_values, metadata, index=index)

    def rewrite(self, index: int, env: Mapping[Any, Any] | None = None) -> Any:
        environment = dict(env or {})
        if index in environment:
            return environment[index]
        if index in self.aliases:
            if environment and not self._aliases_allowed(environment):
                raise _Unresolved(
                    f"endpoint alias node {index} reached under an altered native substitution"
                )
            self.used_nodes.add(index)
            return self.aliases[index]
        if index == self.native_coordinate:
            self.used_nodes.add(index)
            self.endpoint_substitutions.append({"node": index, "endpoint": 1, "reason": "native coordinate to one"})
            return _sp.Integer(1)
        if not environment and index in self.memo:
            return self.memo[index]
        row = self._row(index)
        operation = _operation(row, index)
        self.used_nodes.add(index)
        if operation == "exact_rational":
            numerator = _fraction(row.get("numerator", row.get("value")), f"node {index}.numerator")
            denominator = _fraction(row.get("denominator", 1), f"node {index}.denominator")
            if denominator == 0:
                raise _Unresolved(f"node {index} has a zero exact denominator")
            result = _sp.Rational(numerator.numerator, numerator.denominator) / _sp.Rational(denominator.numerator, denominator.denominator)
        elif operation == "bound_variable":
            result = self._bound(index, row, environment)
        elif operation in {"sum", "add"}:
            result = _simplify(sum((self.rewrite(self._child(index, value, f"arguments[{offset}]"), environment) for offset, value in enumerate(_args(row, index))), _sp.Integer(0)))
        elif operation in {"difference", "sub", "subtract"}:
            values = _args(row, index)
            if len(values) != 2:
                raise _Unresolved(f"node {index} difference requires two arguments")
            result = _simplify(self.rewrite(self._child(index, values[0], "arguments[0]"), environment) - self.rewrite(self._child(index, values[1], "arguments[1]"), environment))
        elif operation in {"product", "mul"}:
            result = _sp.Integer(1)
            for offset, value in enumerate(_args(row, index)):
                result = _simplify(result * self.rewrite(self._child(index, value, f"arguments[{offset}]"), environment))
        elif operation in {"negative", "neg"}:
            child = self._child(index, _field(row, "argument", "operand"), "argument")
            result = _simplify(-self.rewrite(child, environment))
        elif operation in {"positive_quotient", "quotient"}:
            if operation == "positive_quotient":
                result = self._positive_quotient(index, row, environment)
            else:
                numerator = self._child(index, _field(row, "numerator"), "numerator")
                denominator = self._child(index, _field(row, "denominator"), "denominator")
                den = self.rewrite(denominator, environment)
                if _simplify(den) == 0:
                    raise _Unresolved(f"node {index} quotient denominator is symbolically zero")
                result = _simplify(self.rewrite(numerator, environment) / den)
        elif operation == "function_substitution":
            result = self._function_substitution(index, row, environment)
        elif operation in {"analytic_unary", * _ANALYTIC}:
            result = self._analytic(index, row, environment)
        elif operation in {"integer_power", "power"}:
            child = self._child(index, _field(row, "argument", "operand", "base"), "argument")
            exponent = _fraction(row.get("exponent", row.get("power")), f"node {index}.exponent")
            if exponent.denominator != 1:
                raise _Unresolved(f"node {index} has a noninteger power")
            result = _simplify(self.rewrite(child, environment) ** exponent.numerator)
        elif operation == "definite_integral":
            result = self._integral(index, row, environment)
        elif operation == "original_monotone_phase_inverse":
            result = self._phase_inverse(index, row, environment)
        elif operation == "original_lazy_flat_branch":
            result = self._lazy_flat(index, row, environment)
        elif operation in {"exact_original_source_flat_sigma"}:
            argument = self._child(index, _field(row, "argument"), "argument")
            result = self._compound(operation, (self.rewrite(argument, environment),), self._canonical_metadata(row, ("argument",)), index=index)
            self._remember_special(index, row, operation)
        elif operation in {"substitute_original_inverse_angle", "substitute_original_native_coordinate", "substitute_free_query_at_native_coordinate", "inverse_substitution", "native_substitution"}:
            result = self._special_substitution(index, row, operation, environment)
            self._remember_special(index, row, operation)
        elif operation in _COMPOUND_OPS:
            result = self._generic_compound(index, row, operation, environment)
        elif operation in _SOURCE_CONSTANT_OPS:
            metadata = self._canonical_metadata(row)
            # A source constant is opaque and node-qualified.  Only an
            # explicit generic/phase port may identify two such handles.
            metadata["source_node"] = index
            result = self._compound(operation, (), metadata, index=index)
        else:
            raise _Unresolved(f"unsupported source operation at node {index}: {operation!r}")
        result = _simplify(result)
        if not environment:
            self.memo[index] = result
        return result

    def snapshot(self) -> dict[str, Any]:
        return {
            "used_nodes": sorted(self.used_nodes),
            "unresolved": deepcopy(self.unresolved),
            "positive_quotients": deepcopy(self.positive_quotients),
            "substitutions": deepcopy(self.substitutions),
            "endpoint_substitutions": deepcopy(self.endpoint_substitutions),
            "compound_nodes": deepcopy(self.special_nodes),
        }


def _side_form(evaluator: _Evaluator, node: int) -> tuple[Any | None, dict[str, Any]]:
    start = evaluator.snapshot()
    try:
        value = evaluator.rewrite(node)
    except (ValueError, _Unresolved) as exc:
        current = evaluator.snapshot()
        current["status"] = "unresolved"
        current["node"] = node
        current["error"] = str(exc)
        return None, current
    current = evaluator.snapshot()
    current["status"] = "source_form"
    current["node"] = node
    current["value"] = _record(value)
    current["delta_used_nodes"] = sorted(set(current["used_nodes"]) - set(start["used_nodes"]))
    return value, current


def derive(
    nodes: Sequence[Mapping[str, Any]],
    left_chart: Mapping[str, Any],
    right_chart: Mapping[str, Any],
    generic_packet: Mapping[str, Any],
    phase_proof: Mapping[str, Any],
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Return a source-only canonical congruence receipt for the actual DAG.

    ``left_chart`` and ``right_chart`` must expose explicit ``roots``,
    ``loop``, finite ``E_N``/``V_N``, and
    ``finite_N.five_signed_density_increments`` rows.  Their
    native coordinates are taken from an explicit native port when present,
    or from the common-phase coordinate handles in ``phase_proof``.  The
    function does not evaluate a source function numerically.
    """

    if _sp is None:  # pragma: no cover
        raise RuntimeError("SymPy is required for exact source congruence")
    if not _is_sequence(nodes):
        raise ValueError("nodes must be a sequence")
    expected = _require_digest(expected_nodes_digest, "expected_nodes_digest")
    actual = canonical_digest(list(nodes))
    if actual != expected:
        raise ValueError("nodes digest does not match expected_nodes_digest")
    if not isinstance(left_chart, Mapping) or not isinstance(right_chart, Mapping):
        raise ValueError("left_chart and right_chart must be mappings")
    if not isinstance(generic_packet, Mapping):
        raise ValueError("generic_packet must be a mapping")
    if not isinstance(phase_proof, Mapping):
        raise ValueError("phase_proof must be a mapping")
    original_nodes = deepcopy(list(nodes))
    generic, generic_aliases, generic_refs = _generic_rows(generic_packet, len(nodes))
    phase = _validate_phase(nodes, phase_proof)
    left_view, left_owners = _unwrap_chart(left_chart)
    right_view, right_owners = _unwrap_chart(right_chart)
    left_roots = _root_rows(left_view, left_owners, len(nodes))
    right_roots = _root_rows(right_view, right_owners, len(nodes))
    left_targets = _target_rows(left_view, left_owners, len(nodes))
    right_targets = _target_rows(right_view, right_owners, len(nodes))
    left_alias_rows = _validate_explicit_aliases(left_view, left_owners, "postrestore", generic_aliases, len(nodes))
    right_alias_rows = _validate_explicit_aliases(right_view, right_owners, "patch", generic_aliases, len(nodes))
    for name in GENERIC_ROOTS:
        for order in ORDERS:
            accepted = generic[(name, order)]["sides"]
            expected_left = accepted["postrestore"]["node"]
            expected_right = accepted["patch"]["node"]
            if left_roots[name][order] != expected_left:
                raise ValueError(
                    f"postrestore chart root {name}[{order}] does not match its accepted generic source row"
                )
            if right_roots[name][order] != expected_right:
                raise ValueError(
                    f"patch chart root {name}[{order}] does not match its accepted generic source row"
                )
    left_native = _native_port(left_view, left_owners, "postrestore", phase["coordinates"]["postrestore"], len(nodes))
    right_native = _native_port(right_view, right_owners, "patch", phase["coordinates"]["patch"], len(nodes))
    if left_native["coordinate_node"] != phase["coordinates"]["postrestore"]:
        raise ValueError("postrestore native coordinate is not the admitted common-phase coordinate")
    if right_native["coordinate_node"] != phase["coordinates"]["patch"]:
        raise ValueError("patch native coordinate is not the admitted common-phase coordinate")
    phase_atom = _atom("COMMON_PHASE:" + canonical_digest({"phase_nodes": phase["phase_nodes"], "coordinates": phase["coordinates"], "endpoints": phase["endpoints"]}))
    left_aliases = dict(generic_aliases)
    right_aliases = dict(generic_aliases)
    for side, aliases in (("postrestore", left_aliases), ("patch", right_aliases)):
        phase_node = phase["phase_nodes"][side]
        if phase_node in aliases and _expr_key(aliases[phase_node]) != _expr_key(phase_atom):
            raise ValueError(f"{side} common phase handle conflicts with an admitted generic root alias")
    if left_native["coordinate_node"] in left_aliases or right_native["coordinate_node"] in right_aliases:
        raise ValueError("native coordinate handle conflicts with an admitted generic root alias")
    left_aliases[phase["phase_nodes"]["postrestore"]] = phase_atom
    right_aliases[phase["phase_nodes"]["patch"]] = phase_atom
    left_eval = _Evaluator(nodes, side="postrestore", aliases=left_aliases, native_coordinate=left_native["coordinate_node"])
    right_eval = _Evaluator(nodes, side="patch", aliases=right_aliases, native_coordinate=right_native["coordinate_node"])
    rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for field in TARGET_FIELDS:
        for order in ORDERS:
            left_node = left_targets[field][order]
            right_node = right_targets[field][order]
            left_value, left_form = _side_form(left_eval, left_node)
            right_value, right_form = _side_form(right_eval, right_node)
            if left_value is None or right_value is None:
                status = "unresolved"
                unresolved.append({"gate": "source_reduction_unresolved", "field": field, "order": order, "postrestore": left_form.get("error"), "patch": right_form.get("error")})
            else:
                equal = bool(_simplify(left_value - right_value) == 0)
                status = "proved" if equal else "failed"
                if not equal:
                    unresolved.append({"gate": "canonical_congruence_mismatch", "field": field, "order": order, "postrestore": left_form.get("value"), "patch": right_form.get("value")})
            rows.append({
                "field": field,
                "order": order,
                "ordinary_Z_order": order,
                "status": status,
                "postrestore": left_form,
                "patch": right_form,
                "source_only": True,
            })
    phase_refs = _proof_refs(phase_proof.get("proof_refs", phase_proof.get("source_refs", phase_proof.get("proof_ref"))))
    all_refs = _unique(generic_refs + phase_refs)
    passed = len(rows) == len(TARGET_FIELDS) * len(ORDERS) and all(row["status"] == "proved" for row in rows) and not unresolved
    result = {
        "schema": "original_Rm_loop_density_identity_v1",
        "passed": bool(passed),
        "source_only": True,
        "sourceonly": True,
        "numeric_evaluation": False,
        "numeric_containment": False,
        "closure_claim": False,
        "source_digests": {
            "nodes": actual,
            "expected_nodes": expected,
            "left_chart": canonical_digest(left_chart),
            "right_chart": canonical_digest(right_chart),
            "generic_packet": canonical_digest(generic_packet),
            "phase_proof": canonical_digest(phase_proof),
            "module": _module_sha256(),
        },
        "nodes_digest": actual,
        "expected_nodes_digest": expected,
        "source_proof_refs": all_refs,
        "generic_packet": {
            "passed": True,
            "root_rows": len(generic),
            "canonical_values_from_saved_srepr": True,
            "explicit_source_rows_required": True,
        },
        "generic_aliases": [
            {
                "side": side,
                "root": root,
                "order": order,
                "node": generic[(root, order)]["sides"][side]["node"],
                "canonical_srepr": generic[(root, order)]["canonical_srepr"],
                "is_zero": generic[(root, order)]["is_zero"],
                "canonical_atom": "0" if generic[(root, order)]["is_zero"] else str(generic_aliases[generic[(root, order)]["sides"][side]["node"]]),
            }
            for side in SIDES
            for root in GENERIC_ROOTS
            for order in ORDERS
        ],
        "explicit_alias_validation": {"postrestore": left_alias_rows, "patch": right_alias_rows},
        "common_phase": {
            **deepcopy(phase),
            "canonical_atom": str(phase_atom),
            "aliases": [
                {"side": side, "node": phase["phase_nodes"][side], "coordinate_node": phase["coordinates"][side], "endpoint": 1}
                for side in SIDES
            ],
        },
        "native_ports": {"postrestore": left_native, "patch": right_native},
        "target_fields": list(TARGET_FIELDS),
        "ordinary_orders": list(ORDERS),
        "rows": rows,
        "row_passes": [row for row in rows if row["status"] == "proved"],
        "unresolved": unresolved,
        "unresolved_gates": unresolved,
        "used_nodes": {"postrestore": sorted(left_eval.used_nodes), "patch": sorted(right_eval.used_nodes)},
        "scope": {
            "source_graph_only": True,
            "source_graph_reduction_executed": True,
            "actual_saved_graph_run": True,
            "accepted_builder_replay": False,
            "accepted_checker_replay": False,
            "numeric_containment_claimed": False,
            "global_phase_claimed": False,
            "loop_closure_claimed": False,
            "ordinary_Z_orders": [0, 1],
            "fractional_part_differentiated": False,
            "definite_integrals_evaluated": False,
            "unknown_operations_fail_closed": True,
            "inputs_unmutated": deepcopy(list(nodes)) == original_nodes,
        },
        "AST": {
            "generic_packet_passed_and_explicit_rows": True,
            "loop_and_density_ports_explicit": True,
            "common_phase_proof_and_lift_bound": True,
            "native_coordinate_to_one": True,
            "original_inverse_branch_and_positive_guard_preserved": True,
            "flat_predicate_and_lazy_body_preserved": True,
        },
        "LF": {
            "exact_rationals": True,
            "AC_arithmetic_normalized": True,
            "positive_quotient_certificates_recorded": True,
            "alpha_renamed_angle_binders": True,
            "fractional_part_not_differentiated": True,
            "unhandled_operations_reported": True,
        },
    }
    json.dumps(result, allow_nan=False, default=str)
    return result


def prove(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility spelling for callers that name source consumers ``prove``."""

    return derive(*args, **kwargs)


__all__ = [
    "DENSITY_FIELDS",
    "FINITE_FIELDS",
    "GENERIC_ROOTS",
    "LOOP_FIELDS",
    "ORDERS",
    "SIDES",
    "TARGET_FIELDS",
    "canonical_digest",
    "derive",
    "prove",
]
