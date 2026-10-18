"""Append-only rebinding for the corrected reserved own-swirl source.

This module consumes an already materialised exact graph and a small set of
graph handles.  It never imports a producer, loads a report, evaluates a
coefficient, or constructs the original numerical provider.  The rebinder
walks only the supplied coefficient and field roots.  Rows whose operands
change are copied after the supplied graph; rows before the replacement
threshold are treated as immutable source ports.

``CorrectedOwnSwirlProvider`` is the data-only view used by the later physical
consumer.  Its point query returns a previously computed packet verbatim (by
deep copy) at the one supported exact point.  It deliberately has no
arithmetic or provider initialisation path.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from typing import Any, Callable


_ARITHMETIC_FIELDS: dict[str, tuple[str, ...]] = {
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "analytic_unary": ("argument",),
}
_CORE_BINDINGS = (
    "coefficients",
    "response_matrix",
    "linear_response_residuals",
    "repair_velocity",
)
_COEFFICIENT_ALIASES = {
    "coefficients": "coefficients",
    "coefficients5pairs": "coefficients",
}
_RESPONSE_ALIASES = {
    "response_matrix": "response_matrix",
    "response_matrix5x5pairs": "response_matrix",
}
_RESIDUAL_ALIASES = {
    "linear_response_residuals": "linear_response_residuals",
    "linear_response_residuals5pairs": "linear_response_residuals",
}


def canonical_digest(value: Any) -> str:
    """Return the canonical JSON SHA-256 digest for a JSON-compatible value."""

    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _handle(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer graph handle")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def _copy_nodes(value: Any) -> list[dict[str, Any]]:
    if not _is_sequence(value):
        raise TypeError("nodes must be a sequence of mappings")
    rows: list[dict[str, Any]] = []
    for index, row in enumerate(value):
        if not isinstance(row, Mapping):
            raise TypeError(f"nodes[{index}] must be a mapping")
        rows.append(deepcopy(dict(row)))
    if not rows:
        raise ValueError("nodes cannot be empty")
    return rows


def _pair(value: Any, name: str) -> tuple[int, int]:
    if not _is_sequence(value) or len(value) != 2:
        raise ValueError(f"{name} must contain exactly two graph handles")
    return (
        _handle(value[0], f"{name}[0]"),
        _handle(value[1], f"{name}[1]"),
    )


def _five_pairs(value: Any, name: str) -> Any:
    """Validate and return a private copy of a five-pair binding container."""

    if isinstance(value, Mapping):
        # The actual packet uses the five named coefficient roots.  A mapping
        # with other names is accepted only when it is unambiguously a
        # five-entry root map; unrelated metadata is not copied into the
        # returned source view.
        names = ("alpha1", "alpha2", "beta1", "beta2", "beta3")
        if all(name in value for name in names):
            return {
                key: list(_pair(value[key], f"{name}.{key}"))
                for key in names
            }
        if len(value) != 5:
            raise ValueError(f"{name} must contain five coefficient pairs")
        return {
            key: list(_pair(item, f"{name}.{key}"))
            for key, item in value.items()
        }
    if not _is_sequence(value) or len(value) != 5:
        raise ValueError(f"{name} must contain five coefficient pairs")
    return [list(_pair(item, f"{name}[{index}]")) for index, item in enumerate(value)]


def _matrix_pairs(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        if set(value) != {"pairs"}:
            raise ValueError(f"{name} mapping must contain only pairs")
        return {"pairs": _matrix_pairs(value["pairs"], f"{name}.pairs")}
    if not _is_sequence(value) or len(value) != 5:
        raise ValueError(f"{name} must be a five by five pair matrix")
    rows: list[list[list[int]]] = []
    for row_index, row in enumerate(value):
        if not _is_sequence(row) or len(row) != 5:
            raise ValueError(f"{name}[{row_index}] must contain five pairs")
        rows.append([
            list(_pair(item, f"{name}[{row_index}][{column}]"))
            for column, item in enumerate(row)
        ])
    return rows


def _rewrite_pair(value: Any, name: str, rewrite: Callable[[int], int]) -> list[int]:
    old = _pair(value, name)
    return [rewrite(old[0]), rewrite(old[1])]


def _rewrite_five_pairs(value: Any, name: str, rewrite: Callable[[int], int]) -> Any:
    checked = _five_pairs(value, name)
    if isinstance(checked, Mapping):
        return {
            key: _rewrite_pair(item, f"{name}.{key}", rewrite)
            for key, item in checked.items()
        }
    return [
        _rewrite_pair(item, f"{name}[{index}]", rewrite)
        for index, item in enumerate(checked)
    ]


def _rewrite_matrix(value: Any, name: str, rewrite: Callable[[int], int]) -> Any:
    checked = _matrix_pairs(value, name)
    if isinstance(checked, Mapping):
        return {"pairs": _rewrite_matrix(checked["pairs"], f"{name}.pairs", rewrite)}
    return [
        [
            _rewrite_pair(item, f"{name}[{row_index}][{column}]", rewrite)
            for column, item in enumerate(row_value)
        ]
        for row_index, row_value in enumerate(checked)
    ]


def _rewrite_residuals(value: Any, name: str, rewrite: Callable[[int], int]) -> Any:
    checked = _five_pairs(value, name)
    if isinstance(checked, Mapping):
        return {
            key: _rewrite_pair(item, f"{name}.{key}", rewrite)
            for key, item in checked.items()
        }
    return [
        _rewrite_pair(item, f"{name}[{index}]", rewrite)
        for index, item in enumerate(checked)
    ]


def _rewrite_repair_velocity(
    value: Any, name: str, rewrite: Callable[[int], int]
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    required = ("Uz1", "Utheta1", "F1", "physical_radius", "log_R_over_Rc")
    missing = [key for key in required if key not in value]
    if missing:
        raise ValueError(f"{name} is missing {', '.join(missing)}")
    return {
        "Uz1": _rewrite_pair(value["Uz1"], f"{name}.Uz1", rewrite),
        "Utheta1": _rewrite_pair(value["Utheta1"], f"{name}.Utheta1", rewrite),
        "F1": _rewrite_pair(value["F1"], f"{name}.F1", rewrite),
        "physical_radius": rewrite(
            _handle(value["physical_radius"], f"{name}.physical_radius")
        ),
        "log_R_over_Rc": rewrite(
            _handle(value["log_R_over_Rc"], f"{name}.log_R_over_Rc")
        ),
    }


def _rewrite_optional_m15(
    value: Any, path: str, rewrite: Callable[[int], int]
) -> list[int]:
    return _rewrite_pair(value, path, rewrite)


def _normalise_replacements(
    replacements: Mapping[Any, Any], node_count: int
) -> dict[int, int]:
    if not isinstance(replacements, Mapping) or not replacements:
        raise TypeError("replacements must be a nonempty mapping")
    result: dict[int, int] = {}
    for old, new in replacements.items():
        old_handle = _handle(old, "replacement source")
        new_handle = _handle(new, "replacement target")
        if old_handle >= node_count:
            raise ValueError("replacement source is outside nodes")
        # Replacement rows are already present in the immutable prefix.  They
        # are ports, so they are intentionally not traversed or copied.
        if new_handle >= node_count:
            raise ValueError("replacement target must be in the supplied prefix")
        if old_handle == new_handle:
            raise ValueError("replacement source and target must differ")
        result[old_handle] = new_handle
    return result


def _root_limit_guard(handle: int, path: str) -> None:
    if handle > 53397:
        raise ValueError(f"{path} source root exceeds 53397")


def _validate_root_pair(value: Any, path: str) -> tuple[int, int]:
    pair = _pair(value, path)
    _root_limit_guard(pair[0], path)
    _root_limit_guard(pair[1], path)
    return pair


def _validate_root_scalar(value: Any, path: str) -> int:
    handle = _handle(value, path)
    _root_limit_guard(handle, path)
    return handle


def _collect_core_roots(bindings: Mapping[str, Any]) -> dict[str, Any]:
    """Validate only the selected root fields and return a shape descriptor."""

    if not isinstance(bindings, Mapping):
        raise TypeError("bindings must be a mapping")
    selected: dict[str, Any] = {}

    coefficient_key = next((key for key in _COEFFICIENT_ALIASES if key in bindings), None)
    if coefficient_key is None:
        raise ValueError("bindings.coefficients is required")
    coefficients = _five_pairs(bindings[coefficient_key], coefficient_key)
    if isinstance(coefficients, Mapping):
        for key, pair in coefficients.items():
            _validate_root_pair(pair, f"{coefficient_key}.{key}")
    else:
        for index, pair in enumerate(coefficients):
            _validate_root_pair(pair, f"{coefficient_key}[{index}]")
    selected[coefficient_key] = coefficients

    response_key = next((key for key in _RESPONSE_ALIASES if key in bindings), None)
    if response_key is None:
        raise ValueError("bindings.response_matrix is required")
    response = _matrix_pairs(bindings[response_key], response_key)
    response_rows = response["pairs"] if isinstance(response, Mapping) else response
    for row_index, row in enumerate(response_rows):
        for column, pair in enumerate(row):
            _validate_root_pair(pair, f"{response_key}[{row_index}][{column}]")
    selected[response_key] = response

    residual_key = next((key for key in _RESIDUAL_ALIASES if key in bindings), None)
    if residual_key is None:
        raise ValueError("bindings.linear_response_residuals is required")
    residuals = _five_pairs(bindings[residual_key], residual_key)
    if isinstance(residuals, Mapping):
        for key, pair in residuals.items():
            _validate_root_pair(pair, f"{residual_key}.{key}")
    else:
        for index, pair in enumerate(residuals):
            _validate_root_pair(pair, f"{residual_key}[{index}]")
    selected[residual_key] = residuals

    if "repair_velocity" not in bindings:
        raise ValueError("bindings.repair_velocity is required")
    repair = bindings["repair_velocity"]
    if not isinstance(repair, Mapping):
        raise TypeError("bindings.repair_velocity must be a mapping")
    for key in ("Uz1", "Utheta1", "F1"):
        _validate_root_pair(repair.get(key), f"repair_velocity.{key}")
    for key in ("physical_radius", "log_R_over_Rc"):
        _validate_root_scalar(repair.get(key), f"repair_velocity.{key}")
    selected["repair_velocity"] = deepcopy(dict(repair))

    # These are optional caller-provided source ports.  They are deliberately
    # separate from the full-moment replacement and are never inferred from it.
    if "moments" in bindings:
        moments = bindings["moments"]
        if not isinstance(moments, Mapping) or "m15" not in moments:
            raise ValueError("bindings.moments.m15 is required when moments is supplied")
        _validate_root_pair(moments["m15"], "moments.m15")
        selected["moments"] = {"m15": list(_pair(moments["m15"], "moments.m15"))}
    if "segments" in bindings:
        segments = bindings["segments"]
        if not isinstance(segments, Mapping):
            raise TypeError("bindings.segments must be a mapping")
        pre_rp = segments.get("preRp")
        if not isinstance(pre_rp, Mapping) or "m15" not in pre_rp:
            raise ValueError(
                "bindings.segments.preRp.m15 is required when segments is supplied"
            )
        _validate_root_pair(pre_rp["m15"], "segments.preRp.m15")
        selected["segments"] = {"preRp": {"m15": list(_pair(pre_rp["m15"], "segments.preRp.m15"))}}
    return selected


def _operand_references(row: Mapping[str, Any], index: int) -> list[tuple[str, int]]:
    operation = row.get("operation")
    fields = _ARITHMETIC_FIELDS.get(operation)
    if fields is None:
        raise ValueError(
            f"visited node {index} uses unsupported operation {operation!r}; "
            "only the declared arithmetic operand schemas are allowed"
        )
    references: list[tuple[str, int]] = []
    for field in fields:
        if field not in row:
            raise ValueError(f"visited node {index} ({operation}) is missing {field}")
        value = row[field]
        if field == "arguments":
            if not _is_sequence(value) or not value:
                raise ValueError(f"visited node {index} arguments must be nonempty")
            values = value
        else:
            values = (value,)
        for offset, child in enumerate(values):
            references.append(
                (
                    f"{field}[{offset}]" if field == "arguments" else field,
                    _handle(child, f"node {index}.{field}"),
                )
            )
    return references


def _replace_operands(
    row: dict[str, Any], references: list[tuple[str, int]], mapped: list[int]
) -> None:
    operation = row["operation"]
    fields = _ARITHMETIC_FIELDS[operation]
    cursor = 0
    for field in fields:
        if field == "arguments":
            original = row[field]
            values = list(original)
            for offset in range(len(values)):
                values[offset] = mapped[cursor]
                cursor += 1
            row[field] = tuple(values) if isinstance(original, tuple) else values
        else:
            row[field] = mapped[cursor]
            cursor += 1


def _updated_binding_roots(
    selected: Mapping[str, Any], rewrite: Callable[[int], int]
) -> dict[str, Any]:
    updated: dict[str, Any] = {}
    for key, value in selected.items():
        if key in _COEFFICIENT_ALIASES:
            updated[key] = _rewrite_five_pairs(value, key, rewrite)
        elif key in _RESPONSE_ALIASES:
            updated[key] = _rewrite_matrix(value, key, rewrite)
        elif key in _RESIDUAL_ALIASES:
            updated[key] = _rewrite_residuals(value, key, rewrite)
        elif key == "repair_velocity":
            updated[key] = _rewrite_repair_velocity(value, key, rewrite)
        elif key == "moments":
            updated[key] = {"m15": _rewrite_optional_m15(value["m15"], "moments.m15", rewrite)}
        elif key == "segments":
            updated[key] = {
                "preRp": {
                    "m15": _rewrite_optional_m15(
                        value["preRp"]["m15"], "segments.preRp.m15", rewrite
                    )
                }
            }
    return updated


def rebind(
    nodes: Sequence[Mapping[str, Any]],
    bindings: Mapping[str, Any],
    replacements: Mapping[int, int],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Rebind selected source roots into an immutable append-only graph.

    The function copies no existing rows into ``delta_nodes``.  It appends a
    copied row only when a visited post-threshold arithmetic operand changes.
    The returned ``updated_bindings`` contains only the four core source
    fields, plus optional ``moments.m15`` and ``segments.preRp.m15`` supplied
    by the caller.  It never infers a segment root from a full-moment
    replacement.
    """

    graph = _copy_nodes(nodes)
    actual_digest = canonical_digest(graph)
    if not isinstance(expected_nodes_digest, str) or expected_nodes_digest != actual_digest:
        raise ValueError("nodes canonical digest does not match expected_nodes_digest")
    replacement_map = _normalise_replacements(replacements, len(graph))
    selected = _collect_core_roots(bindings)
    threshold = min(replacement_map)

    cache: dict[int, int] = {}
    visiting: set[int] = set()
    delta: list[dict[str, Any]] = []

    def rebind_handle(handle: int) -> int:
        handle = _handle(handle, "source handle")
        if handle >= len(graph):
            raise ValueError(f"source handle {handle} is outside nodes")
        if handle in cache:
            return cache[handle]
        if handle in replacement_map:
            mapped = replacement_map[handle]
            cache[handle] = mapped
            return mapped
        if handle < threshold:
            cache[handle] = handle
            return handle
        if handle in visiting:
            raise ValueError(f"cycle detected at visited node {handle}")
        visiting.add(handle)
        row = graph[handle]
        references = _operand_references(row, handle)
        mapped_operands: list[int] = []
        changed = False
        for path, child in references:
            if child >= handle:
                raise ValueError(
                    f"forward reference from node {handle} {path} to node {child}"
                )
            mapped_child = rebind_handle(child)
            mapped_operands.append(mapped_child)
            changed = changed or mapped_child != child
        if not changed:
            mapped = handle
        else:
            copied = deepcopy(row)
            _replace_operands(copied, references, mapped_operands)
            mapped = len(graph) + len(delta)
            for _, child in _operand_references(copied, mapped):
                if child >= mapped:
                    raise ValueError(
                        f"appended node {mapped} has a forward operand {child}"
                    )
            delta.append(copied)
        visiting.remove(handle)
        cache[handle] = mapped
        return mapped

    updated = _updated_binding_roots(selected, rebind_handle)

    # The corrected RHS changes only the beta/swirl branch.  These ports are
    # source identities and must remain byte-for-byte equal to their old roots.
    response_name = next(key for key in _RESPONSE_ALIASES if key in selected)
    if updated[response_name] != selected[response_name]:
        raise ValueError("response_matrix roots changed during own-swirl rebinding")
    coefficient_name = next(key for key in _COEFFICIENT_ALIASES if key in selected)
    original_coefficients = selected[coefficient_name]
    rebound_coefficients = updated[coefficient_name]
    if isinstance(original_coefficients, Mapping):
        for alpha in ("alpha1", "alpha2"):
            if alpha in original_coefficients and rebound_coefficients[alpha] != original_coefficients[alpha]:
                raise ValueError(f"{alpha} root changed during own-swirl rebinding")
    else:
        if rebound_coefficients[:2] != original_coefficients[:2]:
            raise ValueError("alpha1/alpha2 roots changed during own-swirl rebinding")
    original_repair = selected["repair_velocity"]
    if updated["repair_velocity"]["Uz1"] != list(_pair(original_repair["Uz1"], "repair_velocity.Uz1")):
        raise ValueError("Uz1 roots changed during own-swirl rebinding")

    joined = graph + deepcopy(delta)
    new_digest = canonical_digest(joined)
    mapping = dict(cache)
    return {
        "passed": True,
        "source_only": True,
        "source_executed": False,
        "accepted_ancestors_replayed": False,
        "delta_nodes": deepcopy(delta),
        "delta_node_count": len(delta),
        "mapping": mapping,
        "updated_bindings": deepcopy(updated),
        "nodes_digest": actual_digest,
        "new_graph_digest": new_digest,
        "new_graph_node_count": len(joined),
        "new_graph_count": len(joined),
        "replacement_count": len(replacement_map),
        "replacement_handles_in_prefix": True,
        "graph_immutable": True,
        "new_source_terms": False,
        "old_numerics_transferred": False,
        "numeric_Utheta1_evaluated": False,
        "numeric_F1_evaluated": False,
        "global_closure": False,
        "provider_default_chain_reconfigured": False,
        "response_matrix_roots_unchanged": True,
        "alpha1_alpha2_roots_unchanged": True,
        "Uz1_roots_unchanged": True,
    }


def _split_binding_packet(binding_packet: Mapping[str, Any]) -> tuple[Any, Any]:
    if not isinstance(binding_packet, Mapping):
        raise TypeError("binding_packet must be a mapping")
    explicit_fields = binding_packet.get("field_source_bindings")
    explicit_coefficients = binding_packet.get("coefficient_source_bindings")
    if explicit_fields is not None or explicit_coefficients is not None:
        if explicit_fields is None or explicit_coefficients is None:
            raise ValueError(
                "binding_packet must provide both field_source_bindings and coefficient_source_bindings"
            )
        return deepcopy(explicit_fields), deepcopy(explicit_coefficients)
    roots = binding_packet.get("updated_bindings", binding_packet)
    if not isinstance(roots, Mapping):
        raise TypeError("binding_packet.updated_bindings must be a mapping")
    coefficients = {
        key: deepcopy(roots[key])
        for key in (*_COEFFICIENT_ALIASES, *_RESPONSE_ALIASES, *_RESIDUAL_ALIASES)
        if key in roots
    }
    fields = {
        key: deepcopy(roots[key])
        for key in ("repair_velocity", "moments", "segments")
        if key in roots
    }
    if not coefficients and not fields:
        raise ValueError("binding_packet contains no rebound source roots")
    return fields, coefficients


class CorrectedOwnSwirlProvider:
    """Data-only corrected own-swirl source/provider view."""

    def __init__(self, binding_packet: Mapping[str, Any], cached_point_packet: Mapping[str, Any]):
        fields, coefficients = _split_binding_packet(binding_packet)
        if not isinstance(cached_point_packet, Mapping):
            raise TypeError("cached_point_packet must be a mapping")
        self._field_bindings = deepcopy(fields)
        self._coefficient_bindings = deepcopy(coefficients)
        self._cached_point_packet = deepcopy(dict(cached_point_packet))

    def field_source_bindings(self) -> dict[str, Any]:
        """Return rebound field roots without exposing provider-owned state."""

        return deepcopy(self._field_bindings)

    def coefficient_source_bindings(self) -> dict[str, Any]:
        """Return rebound coefficient roots without exposing provider-owned state."""

        return deepcopy(self._coefficient_bindings)

    def evaluate(self, Z: Any, x: Any) -> dict[str, Any]:
        """Return the cached corrected m 15 packet at the one supported point."""

        if type(Z) is not Fraction or Z != Fraction(371, 1000):
            raise ValueError("only exact Z=Fraction(371, 1000) is available")
        if type(x) is not int or x != 2:
            raise ValueError("only x=2 is available")
        return deepcopy(self._cached_point_packet)


__all__ = ["CorrectedOwnSwirlProvider", "canonical_digest", "rebind"]
