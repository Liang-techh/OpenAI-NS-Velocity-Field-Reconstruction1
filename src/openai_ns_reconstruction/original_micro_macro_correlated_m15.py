"""Join the accepted correlated micro source to a new macro C1 source.

The accepted reports are consumed as receipts and saved graph patches.  This
module never imports or executes an accepted producer.  It relocates the
small RQy suffix into the correlated graph namespace, then appends one new
macro six-integrand pressure source.  The returned packet is a source graph;
it deliberately contains no signed numerical evaluations.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

from .original_signed_graph_pairs import GraphPairs, Pair


PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_micro_macro_RQy_seam.json"
)
PARENT_SHA256 = "1ee3ff1234cad6bf79fb24a0a558be50a3c10058a7746c176ea5d4b54d260c3d"
CORRELATED = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_micro_correlated_m15.json"
)
CORRELATED_SHA256 = "7868f1899398a54aaeba1b799b8b1015465718445444966648295e66ae6c97aa"
Q4_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
)
Q4_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
)
Q4_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"
)
Q4_SHA256 = "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
Q4_SOURCE_SHA256 = "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"
Q4_CHECK_SHA256 = "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f"
GRAPH_HELPER = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"
GRAPH_HELPER_SHA256 = "1e9699cd84a1c82dc7b6577f922f1da1362cbedcfef148e7cec1d146a99aaaff"
JOIN_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join_check.json"
)
JOIN_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.py"
)
JOIN_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.json.gz"
)
JOIN_PINS = {
    Path(JOIN_CHECK).name: "49e7d4f0c566abcbfe9833c5cc3c63875b1c80fa6e6025601eab8aaa83956225",
    Path(JOIN_SOURCE).name: "3f0bb8e0363990ab4d851020e988dd463061086d3ad25e7f0e59c047fa3c5bc0",
    Path(JOIN_ARCHIVE).name: "0c673c0979a69f0a6f6fd48e2ee72adae2223f6e384d537d6e39801c15af79e2",
}

BASE_NODE_COUNT = 183700
MACRO_COORDINATE = 70265
MACRO_RADIUS = 71448
MACRO_LENGTH = 66987
MACRO_U = (140685, 140686, 140687, 140688, 140689, 140690)
MACRO_Q = (140780, 140784, 140790, 140796, 140802)
MACRO_Q_Y = (140868, 140872, 140878, 140884, 140890)
MACRO_LOWER_RQY = (141027, 141028, 141029, 141030)
MACRO_UPPER_RQY = (141114, 141115, 141116, 141117)
NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")

_REF_FIELDS: dict[str, tuple[str, ...]] = {
    "sum": ("arguments",),
    "product": ("arguments",),
    "negative": ("argument",),
    "analytic_unary": ("argument",),
    "positive_quotient": ("numerator", "denominator"),
    "function_substitution": ("expression", "variable", "value"),
    "definite_integral": ("integrand", "variable", "lower", "upper"),
}


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def _gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path}")
    return value


def _digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def _refs(row: Mapping[str, Any]) -> list[tuple[str, int]]:
    operation = row.get("operation")
    fields = _REF_FIELDS.get(operation, ())
    result: list[tuple[str, int]] = []
    for field in fields:
        if field not in row:
            raise ValueError(f"{operation} row is missing {field}")
        value = row[field]
        if operation == "definite_integral" and field == "variable" and isinstance(value, str):
            continue
        values = value if field == "arguments" else (value,)
        if not isinstance(values, (list, tuple)):
            raise TypeError(f"{operation}.{field} must be a sequence")
        for offset, child in enumerate(values):
            if isinstance(child, bool) or not isinstance(child, int) or child < 0:
                raise TypeError(f"{operation}.{field}[{offset}] must be a nonnegative handle")
            result.append((f"{field}[{offset}]" if field == "arguments" else field, child))
    return result


def _append_saved(graph: GraphPairs, rows: Sequence[Mapping[str, Any]], start: int) -> list[int]:
    handles: list[int] = []
    for offset, source_row in enumerate(rows):
        row = copy.deepcopy(dict(source_row))
        operation = row.pop("operation", None)
        _require(isinstance(operation, str), f"saved row {offset} has no operation")
        handle = graph.node(operation, **row)
        _require(handle == start + offset, "saved graph suffix changed append order")
        handles.append(handle)
    return handles


def _relocate_row(
    source_row: Mapping[str, Any],
    row_offset: int,
    old_base: int,
    new_base: int,
    suffix_count: int,
) -> dict[str, Any]:
    row = copy.deepcopy(dict(source_row))
    operation = row.get("operation")
    _require(isinstance(operation, str), "relocated graph row has no operation")
    for field in _REF_FIELDS.get(operation, ()):
        if operation == "definite_integral" and field == "variable" and isinstance(row.get(field), str):
            continue
        values = row[field] if field == "arguments" else [row[field]]
        _require(isinstance(values, (list, tuple)), f"{operation}.{field} is not a sequence")
        mapped: list[int] = []
        for child in values:
            if isinstance(child, bool) or not isinstance(child, int) or child < 0:
                raise TypeError("graph reference is not a nonnegative integer")
            if child < old_base:
                mapped.append(child)
            else:
                _require(child < old_base + row_offset, "relocated suffix references a future row")
                mapped.append(new_base + child - old_base)
        row[field] = mapped if field == "arguments" else mapped[0]
    return row


def _relocate_suffix(
    graph: GraphPairs,
    rows: Sequence[Mapping[str, Any]],
    old_base: int,
    new_base: int,
) -> tuple[list[int], list[dict[str, Any]]]:
    handles: list[int] = []
    receipt: list[dict[str, Any]] = []
    for offset, source_row in enumerate(rows):
        relocated = _relocate_row(source_row, offset, old_base, new_base, len(rows))
        before = copy.deepcopy(dict(source_row))
        operation = relocated.pop("operation", None)
        _require(isinstance(operation, str), "relocated row operation missing")
        handle = graph.node(operation, **relocated)
        expected = new_base + offset
        _require(handle == expected, "relocated graph suffix changed append order")
        handles.append(handle)
        ref_names = [field for field, _ in _refs(source_row)]
        relocated_full = {"operation": operation, **relocated}
        before_refs = {field: value for field, value in _refs(source_row)}
        after_refs = {field: value for field, value in _refs(relocated_full)}
        receipt.append(
            {
                "old_handle": old_base + offset,
                "new_handle": handle,
                "operation": operation,
                "reference_fields": ref_names,
                "references_before": before_refs,
                "references_after": after_refs,
                "literal_metadata_preserved": {
                    key: before[key] == relocated_full[key]
                    for key in before
                    if key not in _REF_FIELDS.get(operation, ())
                },
            }
        )
    return handles, receipt


def _map_handle(handle: int, old_base: int, new_base: int, suffix_count: int) -> int:
    if isinstance(handle, bool) or not isinstance(handle, int) or handle < 0:
        raise TypeError("graph handle must be a nonnegative integer")
    if handle < old_base:
        return handle
    _require(handle < old_base + suffix_count, "saved suffix endpoint is outside relocation")
    return new_base + handle - old_base


def _map_handles(handles: Sequence[int], old_base: int, new_base: int, suffix_count: int) -> list[int]:
    return [_map_handle(value, old_base, new_base, suffix_count) for value in handles]


def _poly_add(
    left: Mapping[tuple[str, ...], Fraction],
    right: Mapping[tuple[str, ...], Fraction],
) -> dict[tuple[str, ...], Fraction]:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, Fraction(0)) + coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def _poly_mul(
    left: Mapping[tuple[str, ...], Fraction],
    right: Mapping[tuple[str, ...], Fraction],
) -> dict[tuple[str, ...], Fraction]:
    result: dict[tuple[str, ...], Fraction] = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(sorted(left_monomial + right_monomial))
            result[monomial] = result.get(monomial, Fraction(0)) + left_coefficient * right_coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def _poly_neg(value: Mapping[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    return {monomial: -coefficient for monomial, coefficient in value.items() if coefficient}


def _parse_endpoint_polynomial(
    nodes: Sequence[Mapping[str, Any]],
    index: int,
    aliases: Mapping[int, str],
    memo: dict[int, dict[tuple[str, ...], Fraction]] | None = None,
) -> dict[tuple[str, ...], Fraction]:
    """Parse only the endpoint arithmetic admitted by the source receipt."""

    if memo is None:
        memo = {}
    if index in memo:
        return memo[index]
    if index in aliases:
        result = {(aliases[index],): Fraction(1)}
        memo[index] = result
        return result
    row = nodes[index]
    operation = row.get("operation")
    if operation == "exact_rational":
        numerator = row.get("numerator")
        denominator = row.get("denominator")
        _require(isinstance(numerator, int) and isinstance(denominator, int) and denominator != 0,
                 "endpoint exact rational is malformed")
        result = {(): Fraction(numerator, denominator)}
    elif operation == "sum":
        result = {}
        for child in row.get("arguments", ()):
            result = _poly_add(result, _parse_endpoint_polynomial(nodes, child, aliases, memo))
    elif operation == "product":
        result = {(): Fraction(1)}
        for child in row.get("arguments", ()):
            result = _poly_mul(result, _parse_endpoint_polynomial(nodes, child, aliases, memo))
    elif operation == "negative":
        result = _poly_neg(_parse_endpoint_polynomial(nodes, row["argument"], aliases, memo))
    else:
        raise ValueError(f"unexpected endpoint arithmetic operation at {index}: {operation}")
    memo[index] = result
    return result


def _endpoint_polynomial(
    nodes: Sequence[Mapping[str, Any]],
    handle: int,
    side: str,
    order: int,
    physical: Mapping[str, Sequence[int]],
) -> tuple[dict[tuple[str, ...], Fraction], dict[str, Any]]:
    variable = 66988 if side == "micro" else MACRO_COORDINATE
    value = 48 if side == "micro" else 0
    wrapper = nodes[handle]
    _require(wrapper.get("operation") == "function_substitution", f"{side} endpoint wrapper operation changed")
    _require(wrapper.get("variable") == variable and wrapper.get("value") == value,
             f"{side} endpoint cutpoint changed at ordinary Z{order}")
    _require(wrapper.get("Z_independent_substitution") is True,
             f"{side} endpoint substitution is not Z independent at ordinary Z{order}")
    aliases: dict[int, str] = {
        int(physical["radius"]): "R",
    }
    for index, name in zip(physical["Q"][:2], ("Q0", "Q1")):
        aliases[int(index)] = name
    for index, name in zip(physical["U"][:2], ("U0", "U1")):
        aliases[int(index)] = name
    polynomial = _parse_endpoint_polynomial(nodes, wrapper["expression"], aliases)
    return polynomial, {
        "handle": handle,
        "wrapper": copy.deepcopy(dict(wrapper)),
        "side": side,
        "ordinary_Z_order": order,
        "canonical_polynomial": {
            "monomials": [list(monomial) for monomial in sorted(polynomial)],
            "coefficients": [str(polynomial[monomial]) for monomial in sorted(polynomial)],
        },
    }


def _verify_join_receipts(root: Path, q4: Mapping[str, Any]) -> dict[str, Any]:
    check = _json(root / JOIN_CHECK)
    archive = _gzip_json(root / JOIN_ARCHIVE)
    archive_nodes = archive.get("exact_graph_nodes")
    _require(isinstance(archive_nodes, list) and len(archive_nodes) >= BASE_NODE_COUNT
             and archive_nodes[:BASE_NODE_COUNT] == q4.get("exact_graph_nodes"),
             "join archive Q4 graph prefix differs from pinned Q4 archive")
    checks = check.get("checks", {})
    smooth = checks.get("actual_smooth_source_seam", {})
    leading = checks.get("actual_leading_recovery_and_common_N", {})
    _require(smooth.get("passed") is True, "join smooth source receipt failed")
    _require(smooth.get("actual_comparison_phi_U_physical_y1_to_y5_flat_by_original_alpha_Leibniz") is True,
             "join physical phi/U source receipt missing")
    _require(smooth.get("comparison_six_history_jets_from_common_inlet_and_source_ODE") is True,
             "join history source receipt missing")
    _require(smooth.get("actual_comparison_direction_jets_are_common_smooth_Z_compositions") is True,
             "join common smooth Z source receipt missing")
    _require(smooth.get("actual_phi_U_six_histories_and_generic_background_physical_y0_to5_match_by_ODE_induction") is True,
             "join ODE physical source receipt missing")
    _require(leading.get("passed") is True, "join leading recovery receipt failed")
    _require(leading.get("common_U_and_mean_physical_jets_imply_common_Q_then_common_S_and_N_y0_to2") is True,
             "join common physical U/Q source receipt missing")
    _require(leading.get("physical_R_y_equals_R_and_both_native_Jacobians_retained") is True,
             "join physical radius/Jacobian receipt missing")
    bindings = archive.get("micro_macro_mixed_source_join_bindings", {})
    _require(isinstance(bindings, Mapping), "join archive bindings missing")
    _require(bindings.get("coordinate") == MACRO_COORDINATE and bindings.get("radius") == MACRO_RADIUS
             and bindings.get("jacobian") == MACRO_LENGTH and bindings.get("N") == 18,
             "join common cutpoint bindings changed")
    loop = bindings.get("loop", {})
    _require(loop.get("second_micro", {}).get("coordinate") == 66988
             and loop.get("second_micro", {}).get("radius") == 68644
             and loop.get("frozen_macro", {}).get("coordinate") == MACRO_COORDINATE
             and loop.get("frozen_macro", {}).get("radius") == MACRO_RADIUS,
             "join chart cutpoint/radius bindings changed")
    micro_physical = bindings.get("modified", {}).get("second_micro", {}).get("physical", {})
    macro_physical = bindings.get("modified", {}).get("frozen_macro", {}).get("physical", {})
    for fields in (micro_physical, macro_physical):
        _require(isinstance(fields, Mapping), "join physical source fields missing")
        for name in ("U", "Q"):
            _require(isinstance(fields.get(name), list) and len(fields[name]) >= 2,
                     f"join physical {name} source rows missing")
    receipt = {
        "passed": True,
        "check_target": {
            "actual_smooth_source_seam": copy.deepcopy(smooth),
            "actual_leading_recovery_and_common_N": copy.deepcopy(leading),
        },
        "common_cutpoints": {
            "micro_coordinate": loop["second_micro"]["coordinate"],
            "micro_upper_value": 48,
            "macro_coordinate": loop["frozen_macro"]["coordinate"],
            "macro_lower_value": 0,
            "common_N": bindings["N"],
        },
        "physical_source_rows": {
            "micro": {
                "radius": loop["second_micro"]["radius"],
                "U": list(micro_physical["U"][:2]),
                "Q": list(micro_physical["Q"][:2]),
            },
            "macro": {
                "radius": loop["frozen_macro"]["radius"],
                "U": list(macro_physical["U"][:2]),
                "Q": list(macro_physical["Q"][:2]),
            },
        },
        "archive_prefix_equal_to_q4": True,
    }
    return receipt


def _prove_internal_B_seam(
    nodes: Sequence[Mapping[str, Any]],
    micro_fields: Mapping[str, Sequence[int]],
    macro_fields: Mapping[str, Sequence[int]],
    join_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    """Prove the six endpoint B-field rows through common physical leaves."""

    physical = join_receipt["physical_source_rows"]
    micro_physical = dict(physical["micro"])
    macro_physical = dict(physical["macro"])
    fields = ("R_squared_q", "R_squared_q_squared", "R_squared_U_squared")
    rows: list[dict[str, Any]] = []
    expected_micro = {
        "R_squared_q": ("Q0", "R", "R"),
        "R_squared_q_squared": ("Q0", "Q0", "R", "R"),
        "R_squared_U_squared": ("R", "R", "U0", "U0"),
    }
    expected_micro_z1 = {
        "R_squared_q": ("Q1", "R", "R"),
        "R_squared_q_squared": ("Q0", "Q1", "R", "R"),
        "R_squared_U_squared": ("R", "R", "U0", "U1"),
    }
    for order in range(2):
        for field in fields:
            micro_polynomial, micro_row = _endpoint_polynomial(
                nodes, int(micro_fields[field][order]), "micro", order, micro_physical,
            )
            macro_polynomial, macro_row = _endpoint_polynomial(
                nodes, int(macro_fields[field][order]), "macro", order, macro_physical,
            )
            _require(micro_polynomial == macro_polynomial,
                     f"common physical endpoint polynomial mismatch for {field} at ordinary Z{order}")
            expected_monomial = (expected_micro if order == 0 else expected_micro_z1)[field]
            expected_coefficient = Fraction(2 if order == 1 and field != "R_squared_q" else 1)
            _require(micro_polynomial == {tuple(sorted(expected_monomial)): expected_coefficient},
                     f"unexpected endpoint polynomial for {field} at ordinary Z{order}")
            rows.append({
                "ordinary_Z_order": order,
                "field": field,
                "micro": micro_row,
                "macro": macro_row,
                "canonical_difference": {},
                "identity": "common physical U/Q and radius leaves give identical R^2 endpoint row",
            })
    return {
        "passed": True,
        "ordinary_Z_orders": [0, 1],
        "fields": list(fields),
        "rows": rows,
        "row_count": len(rows),
        "canonical_digest": _digest(rows),
        "source_receipt": {
            "join_archive_physical_source_rows": copy.deepcopy(dict(physical)),
            "join_archive_common_cutpoints": copy.deepcopy(dict(join_receipt["common_cutpoints"])),
            "accepted_join_common_physical_source": True,
            "accepted_join_common_Q_y_rows": True,
        },
    }


def _verify_increment_graph_rows(
    nodes: Sequence[Mapping[str, Any]],
    increment: Sequence[int],
    S1: Sequence[int],
    upper: Sequence[int],
    lower: Sequence[int],
    label: str,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    _require(len(increment) == 2 and len(S1) == 2 and len(upper) == 2 and len(lower) == 2,
             f"{label} increment must have two ordinary Z rows")
    for order in range(2):
        row = nodes[increment[order]]
        _require(row.get("operation") == "sum" and len(row.get("arguments", ())) == 2,
                 f"{label} increment row changed at ordinary Z{order}")
        half_handle, negative_handle = row["arguments"]
        _require(nodes[half_handle] == {
            "operation": "product",
            "arguments": [S1[order], 28],
        }, f"{label} half-S1 row changed at ordinary Z{order}")
        negative_row = nodes[negative_handle]
        _require(negative_row.get("operation") == "negative",
                 f"{label} endpoint sign changed at ordinary Z{order}")
        boundary_handle = negative_row["argument"]
        boundary_row = nodes[boundary_handle]
        _require(boundary_row.get("operation") == "sum" and len(boundary_row.get("arguments", ())) == 2,
                 f"{label} boundary row changed at ordinary Z{order}")
        boundary_negative = boundary_row["arguments"][1]
        _require(boundary_row["arguments"][0] == upper[order],
                 f"{label} upper endpoint changed at ordinary Z{order}")
        _require(nodes[boundary_negative] == {
            "operation": "negative",
            "argument": lower[order],
        }, f"{label} lower endpoint sign changed at ordinary Z{order}")
        rows.append({
            "ordinary_Z_order": order,
            "increment": increment[order],
            "half_S1": half_handle,
            "boundary": boundary_handle,
            "upper": upper[order],
            "lower": lower[order],
            "row": copy.deepcopy(dict(row)),
        })
    return {"passed": True, "label": label, "rows": rows}


def _verify_endpoint_wrappers(
    nodes: Sequence[Mapping[str, Any]],
    fields: Mapping[str, Sequence[int]],
    variable: int,
    value: int,
    label: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for field in ("R_squared_q", "R_squared_q_squared", "R_squared_U_squared"):
        for order, handle in enumerate(fields[field]):
            row = nodes[handle]
            _require(row.get("operation") == "function_substitution"
                     and row.get("variable") == variable
                     and row.get("value") == value
                     and row.get("Z_independent_substitution") is True,
                     f"{label} endpoint wrapper changed for {field} at ordinary Z{order}")
            rows.append({
                "field": field,
                "ordinary_Z_order": order,
                "handle": handle,
                "row": copy.deepcopy(dict(row)),
            })
    return rows


def _depends_on_z(nodes: Sequence[Mapping[str, Any]], index: int, memo: dict[int, bool]) -> bool:
    if index in memo:
        return memo[index]
    row = nodes[index]
    if row.get("operation") == "bound_variable" and row.get("name") == "Z":
        memo[index] = True
        return True
    result = any(_depends_on_z(nodes, child, memo) for _, child in _refs(row))
    memo[index] = result
    return result


def _verify_parent(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    path = root / PARENT
    _require(_sha(path) == PARENT_SHA256, "RQy seam parent report changed")
    on_disk = _json(path)
    _require(on_disk == parent, "build parent does not equal the pinned RQy seam report")
    seam = parent.get("actual_micro_macro_RQy_seam")
    _require(isinstance(seam, Mapping) and seam.get("passed") is True, "accepted RQy seam packet missing")
    _require(seam.get("graph_patch", {}).get("base_node_count") == BASE_NODE_COUNT, "RQy base prefix changed")
    _require(seam.get("graph_patch", {}).get("appended_node_count") == 102, "RQy suffix size changed")
    _require(seam.get("graph_patch", {}).get("base_prefix_unchanged") is True, "RQy prefix receipt failed")
    _require(seam.get("contract", {}).get("cancellation_only_in_adjacent_joined_expression") is True,
             "RQy cancellation scope changed")
    _require(seam.get("contract", {}).get("no_outer_cap_deletion") is True,
             "RQy outer source was deleted")
    _require(seam.get("contract", {}).get("whole_first_to_macro_pair_not_claimed") is True,
             "RQy parent scope unexpectedly claims an absolute first-to-macro pair")
    guards = seam.get("source_guards", {})
    _require(guards.get("passed") is True, "accepted source join guards failed")
    _require(guards.get("accepted_join_common_physical_source") is True,
             "common physical source receipt missing")
    full_qy = guards.get("accepted_join_full_Q_y_rows")
    if not isinstance(full_qy, Mapping):
        full_qy = guards.get("q4_rows", {}).get("physical_Q_y_rows_verified")
    _require(isinstance(full_qy, Mapping), "full Q_y source rows missing")
    _require(guards.get("join_archive_full_Q_y_bindings") is True,
             "accepted full Q_y graph binding receipt missing")
    _require(guards.get("full_Q4_physical_Q_y_rows_bound") is True,
             "accepted full Q_y physical-row receipt missing")
    _require(guards.get("accepted_join_physical_R_and_native_J_distinct") is True,
             "physical radius/native Jacobian distinction missing")
    smooth = guards.get("accepted_smooth_source_premises", {})
    for name in (
        "common_global_log_radius_phase_at_micro2_macro0_and_one_N",
        "original_monotone_inverse_uniqueness_and_periodic_lift_transfer_all_y4_Z5_loop_rows",
        "actual_micro_complete_correction_is_macro_incoming_so_mean_values_match",
        "original_mean_ODE_and_one_minus_delta_recovery_give_common_y4_Z4_Q",
    ):
        _require(smooth.get(name) is True, f"accepted smooth premise missing: {name}")
    proof_rows = seam.get("proof_rows")
    _require(isinstance(proof_rows, list) and len(proof_rows) == 4,
             "top-level RQy proof rows missing")
    expected_micro_upper = [133745, 133746, 133747, 133748]
    expected_macro_lower = [183704, 183705, 183706, 183707]
    for order, proof in enumerate(proof_rows):
        _require(isinstance(proof, Mapping) and proof.get("ordinary_Z_order") == order,
                 f"RQy proof row order changed at ordinary Z{order}")
        second = proof.get("second_micro", {})
        macro = proof.get("frozen_macro", {})
        _require(second.get("R_Q_y_upper_at_s_2") == expected_micro_upper[order],
                 f"RQy proof micro endpoint changed at ordinary Z{order}")
        _require(macro.get("R_Q_y_lower_at_m_0") == expected_macro_lower[order],
                 f"RQy proof macro endpoint changed at ordinary Z{order}")
        micro_identity = proof.get("factorized_micro_endpoint_identity", {})
        macro_identity = proof.get("factorized_macro_lower_identity", {})
        conditional = proof.get("conditional_boundary_identity", {})
        _require(micro_identity.get("passed") is True
                 and micro_identity.get("actual_wrapper") == expected_micro_upper[order]
                 and micro_identity.get("variable") == 66988
                 and micro_identity.get("value") == 48
                 and micro_identity.get("Z_independent_substitution") is True,
                 f"RQy micro factorized identity missing at ordinary Z{order}")
        _require(macro_identity.get("passed") is True
                 and macro_identity.get("actual_wrapper") == expected_macro_lower[order]
                 and macro_identity.get("variable") == MACRO_COORDINATE
                 and macro_identity.get("value") == 0
                 and macro_identity.get("Z_independent_substitution") is True,
                 f"RQy macro factorized identity missing at ordinary Z{order}")
        for name in (
            "common_Q_y_source_rows",
            "common_physical_radius_cutpoint",
            "one_common_N_and_global_phase",
            "mean_ODE_and_complete_correction_premises",
        ):
            _require(conditional.get(name) is True,
                     f"RQy conditional identity missing {name} at ordinary Z{order}")
    result = dict(seam)
    result["accepted_full_Q_y_rows"] = copy.deepcopy(dict(full_qy))
    result["proof_rows"] = copy.deepcopy(proof_rows)
    return result


def _verify_correlated(
    root: Path,
    latest_parent_packet: Mapping[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    path = root / CORRELATED
    _require(_sha(path) == CORRELATED_SHA256, "correlated micro report changed")
    report = _json(path)
    packet = report.get("actual_micro_correlated_m15")
    _require(isinstance(packet, Mapping) and packet.get("passed") is True, "correlated packet missing")
    _require(
        isinstance(latest_parent_packet, Mapping),
        "latest RQy parent lacks the correlated micro packet identity",
    )
    _require(
        dict(packet) == dict(latest_parent_packet),
        "standalone correlated packet differs from latest RQy parent's correlated packet",
    )
    _require(packet.get("accepted_ancestors_replayed") is False, "accepted correlated producer was replayed")
    _require(packet.get("accepted_seam_fingerprints_replayed") is False, "accepted seam checker was replayed")
    _require(packet.get("actual_full_Q4_six_I1_source_implemented") is True, "correlated six-I1 source missing")
    _require(packet.get("physical_I1_measure") == "R^2*hb*ds exactly once", "correlated measure receipt changed")
    exact = packet.get("exact_algebra", {})
    _require(exact.get("passed") is True and exact.get("no_Z_factor_on_U_QZ_volume_coefficient") is True,
             "accepted exact pressure receipt missing")
    _require(all(value is True for value in exact.get("checks", {}).values()), "accepted pressure identity failed")
    graph_patch = packet.get("graph_patch", {})
    _require(graph_patch.get("base_graph_node_count") == BASE_NODE_COUNT, "correlated Q4 prefix changed")
    _require(graph_patch.get("base_sha256") == Q4_SHA256, "correlated Q4 archive pin changed")
    rows = graph_patch.get("appended_nodes")
    _require(isinstance(rows, list) and len(rows) == 624, "correlated suffix size changed")
    _require(graph_patch.get("old_graph_unchanged") is True, "correlated Q4 prefix mutation receipt failed")
    return dict(packet), dict(graph_patch)


def _verify_source_pins(root: Path, seam: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    expected = {
        Q4_ARCHIVE: Q4_SHA256,
        Q4_SOURCE: Q4_SOURCE_SHA256,
        Q4_CHECK: Q4_CHECK_SHA256,
        GRAPH_HELPER: GRAPH_HELPER_SHA256,
        JOIN_CHECK: JOIN_PINS[Path(JOIN_CHECK).name],
        JOIN_SOURCE: JOIN_PINS[Path(JOIN_SOURCE).name],
        JOIN_ARCHIVE: JOIN_PINS[Path(JOIN_ARCHIVE).name],
    }
    pins: dict[str, dict[str, Any]] = {}
    for path, digest in expected.items():
        actual = _sha(root / path)
        _require(actual == digest, f"source pin changed: {path}")
        pins[path] = {"path": path, "expected_sha256": digest, "passed": True}
    report_pins = seam.get("source_pins", {})
    for name, digest in JOIN_PINS.items():
        _require(report_pins.get(name) == digest, f"RQy receipt pin changed: {name}")
    return pins


def _macro_source(
    graph: GraphPairs,
    seam_ports: Mapping[str, Sequence[int]],
    q4: Mapping[str, Any],
    rqy_start: int,
) -> dict[str, Any]:
    bindings = q4["modified_Q4_pressure_C3_bindings"]["charts"]["frozen_macro"]
    physical = bindings["physical"]
    _require(tuple(physical["U"]) == MACRO_U, "macro U source rows changed")
    _require(tuple(physical["Q"]) == MACRO_Q, "macro Q source rows changed")
    _require(tuple(physical["Q_y"]) == MACRO_Q_Y, "macro Q_y source rows changed")
    _require(tuple(bindings["finite_boundary"]["R_Q_y"]) == MACRO_LOWER_RQY,
             "macro lower R_Q_y source rows changed")
    for index, expression in zip(MACRO_LOWER_RQY, MACRO_Q_Y[:4]):
        row = q4["exact_graph_nodes"][index]
        _require(row == {"operation": "product", "arguments": [expression, MACRO_RADIUS]},
                 f"macro R_Q_y factor changed at {index}")
    nodes = graph.nodes
    memo: dict[int, bool] = {}
    _require(not _depends_on_z(nodes, MACRO_RADIUS, memo), "macro radius depends on free Z")
    _require(not _depends_on_z(nodes, MACRO_LENGTH, memo), "macro length depends on free Z")
    _require(nodes[MACRO_RADIUS] == {"operation": "product", "arguments": [65604, 71447]},
             "macro radius source row changed")
    _require(nodes[MACRO_LENGTH] == {"operation": "sum", "arguments": [66984, 66986]},
             "macro length source row changed")
    _require(tuple(seam_ports["frozen_macro_lower"]) == tuple(rqy_start + offset for offset in range(4, 8)),
             "RQy lower target relocation mismatch")
    for offset, handle in enumerate(seam_ports["frozen_macro_lower"]):
        _require(
            nodes[handle]
            == {
                "operation": "function_substitution",
                "expression": MACRO_LOWER_RQY[offset],
                "variable": MACRO_COORDINATE,
                "value": 0,
                "Z_independent_substitution": True,
                "source_endpoint_trace": True,
            },
            "relocated macro lower endpoint changed",
        )
    for offset, handle in enumerate(MACRO_UPPER_RQY):
        _require(
            nodes[handle]
            == {
                "operation": "function_substitution",
                "expression": MACRO_LOWER_RQY[offset],
                "variable": MACRO_COORDINATE,
                "value": 1,
                "Z_independent_substitution": True,
            },
            "macro upper endpoint changed",
        )

    suffix_start_node = len(graph.nodes)
    suffix_start_append = len(graph.appended)
    z: Pair = (20, graph.constant(1))
    delta = graph.fixed(25220)
    one = graph.fixed(graph.constant(1))
    two_id = graph.constant(2)
    d = graph.padd(one, graph.pneg(graph.pmul(z, z)))
    h = graph.padd(one, graph.pneg(delta))
    L = graph.padd(one, graph.pneg(graph.pmul(delta, z, z)))
    L2 = graph.pmul(L, L)
    divide = lambda pair, denominator=L: graph.pdiv(pair, denominator, "same positive physical L=1-delta*Z^2")
    half = lambda pair: graph.pscale(pair, graph.constant("1/2"))
    U: Pair = tuple(physical["U"][:2])  # type: ignore[assignment]
    U_Z: Pair = tuple(physical["U"][1:3])  # type: ignore[assignment]
    Q: Pair = tuple(physical["Q"][:2])  # type: ignore[assignment]
    Q_Z: Pair = tuple(physical["Q"][1:3])  # type: ignore[assignment]
    source = {
        "q": Q,
        "q_Z": Q_Z,
        "q_squared": graph.pmul(Q, Q),
        "U_q_Z": graph.pmul(U, Q_Z),
        "U_squared": graph.pmul(U, U),
        "U_U_Z": graph.pmul(U, U_Z),
    }
    radius_squared = graph.mul(MACRO_RADIUS, MACRO_RADIUS)
    weight_scalar = graph.mul(radius_squared, MACRO_LENGTH)
    weight = graph.fixed(weight_scalar)
    integrals: dict[str, Pair] = {}
    for name in NAMES:
        integrals[name] = graph.integral(
            graph.pmul(weight, source[name]),
            MACRO_COORDINATE,
            0,
            1,
            f"actual_macro_correlated_I1_{name}",
            {
                "chart": "frozen_macro",
                "physical_radius_power": 2,
                "physical_length_measure_applied_once": True,
                "no_total_y_derivative_in_volume": True,
                "ordinary_Z_rows": [0, 1],
            },
        )
    lower = {
        "R_squared_q": graph.at(graph.pscale(Q, radius_squared), MACRO_COORDINATE, 0),
        "R_squared_q_squared": graph.at(graph.pscale(source["q_squared"], radius_squared), MACRO_COORDINATE, 0),
        "R_squared_U_squared": graph.at(graph.pscale(source["U_squared"], radius_squared), MACRO_COORDINATE, 0),
        "R_Q_y": tuple(seam_ports["frozen_macro_lower"][:2]),
    }
    upper = {
        "R_squared_q": graph.at(graph.pscale(Q, radius_squared), MACRO_COORDINATE, 1),
        "R_squared_q_squared": graph.at(graph.pscale(source["q_squared"], radius_squared), MACRO_COORDINATE, 1),
        "R_squared_U_squared": graph.at(graph.pscale(source["U_squared"], radius_squared), MACRO_COORDINATE, 1),
        "R_Q_y": tuple(MACRO_UPPER_RQY[:2]),
    }
    boundary = lambda key: graph.padd(upper[key], graph.pneg(lower[key]))
    B = graph.padd(
        divide(boundary("R_squared_q")),
        half(boundary("R_squared_q_squared")),
        graph.pneg(graph.pscale(divide(graph.pmul(z, z, boundary("R_squared_U_squared")), L2), two_id)),
    )
    D = graph.padd(
        graph.pneg(divide(integrals["q"])),
        half(divide(graph.pmul(h, z, integrals["q_Z"]))),
        graph.pneg(half(integrals["q_squared"])),
        divide(graph.pmul(d, integrals["U_q_Z"])),
        graph.pscale(divide(graph.pmul(h, z, z, integrals["U_squared"]), L2), two_id),
        graph.pscale(divide(graph.pmul(z, d, integrals["U_U_Z"]), L2), two_id),
    )
    S1 = graph.padd(B, D)
    increment = graph.padd(half(S1), graph.pneg(boundary("R_Q_y")))
    return {
        "passed": True,
        "ordinary_Z_orders": [0, 1],
        "six_I1_names": list(NAMES),
        "physical_fields": {"U": list(physical["U"][:3]), "Q": list(physical["Q"][:3]),
                             "U_Z": list(physical["U"][1:3]), "Q_Z": list(physical["Q"][1:3])},
        "physical_measure": {
            "radius": MACRO_RADIUS,
            "radius_power": 2,
            "length": MACRO_LENGTH,
            "weight": weight_scalar,
            "R_squared_length_applied_once": True,
            "Z_independent": True,
        },
        "coordinate": MACRO_COORDINATE,
        "domain": [0, 1],
        "source": source,
        "integrals": integrals,
        "endpoint_fields": {"lower": lower, "upper": upper},
        "formula_context": {
            "z": list(z),
            "delta": list(delta),
            "L": list(L),
            "L2": list(L2),
            "two_id": two_id,
        },
        "reduced_boundary_B": B,
        "reduced_volume_D": D,
        "reduced_S1": S1,
        "macro_increment": increment,
        "outer_RQy": {"lower": list(seam_ports["frozen_macro_lower"]), "upper": list(MACRO_UPPER_RQY)},
        "graph_suffix_start": suffix_start_node,
        "graph_suffix_append_index": suffix_start_append,
        "formula": {
            "B": "R^2 Q/L + R^2 Q^2/2 - 2 Z^2 R^2 U^2/L^2",
            "D": "-q/L + h Z q_Z/(2L) - q^2/2 + d U q_Z/L + 2 h Z^2 U^2/L^2 + 2 Z d U U_Z/L^2",
            "increment": "(B+D)/2 - R_Q_y(upper) + R_Q_y(lower)",
            "L": "1-delta*Z^2",
            "coefficients": {
                "U_squared": "2*h*Z^2/L^2",
                "U_U_Z": "2*Z*d/L^2",
                "U_q_Z": "d/L",
            },
        },
        "exact_pressure_source_receipt": {
            "source_formula": "physical reduced_pressure_S with R^2*length measure",
            "no_extra_Z_on_U_q_Z": True,
            "Pstar_not_added_to_I1_weight": True,
        },
        "appended_nodes": copy.deepcopy(graph.appended[suffix_start_append:]),
    }


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Return one relocated, source-linked graph namespace for micro plus macro."""

    root_path = Path(root)
    seam = _verify_parent(root_path, parent)
    source_pins = _verify_source_pins(root_path, seam)
    correlated, correlated_patch = _verify_correlated(
        root_path,
        parent.get("actual_micro_correlated_m15"),
    )

    q4 = _gzip_json(root_path / Q4_ARCHIVE)
    base_nodes = q4.get("exact_graph_nodes")
    _require(isinstance(base_nodes, list) and len(base_nodes) == BASE_NODE_COUNT, "Q4 graph prefix length changed")
    join_receipt = _verify_join_receipts(root_path, q4)
    _require(correlated_patch["base_archive"] == Q4_ARCHIVE, "correlated base archive changed")
    _require(seam["graph_patch"]["base_archive"].endswith(Path(Q4_ARCHIVE).name), "RQy base archive changed")
    graph = GraphPairs(base_nodes)
    base_snapshot = copy.deepcopy(base_nodes)
    correlated_suffix_start = len(graph.nodes)
    correlated_ids = _append_saved(graph, correlated_patch["appended_nodes"], correlated_suffix_start)
    _require(graph.nodes[:BASE_NODE_COUNT] == base_snapshot, "Q4 graph prefix was mutated")
    _require(correlated_ids[-1] == BASE_NODE_COUNT + 623, "correlated suffix endpoint changed")

    rqy_patch = seam["graph_patch"]
    rqy_rows = rqy_patch["appended_nodes"]
    _require(isinstance(rqy_rows, list) and len(rqy_rows) == 102, "RQy appended graph missing")
    rqy_start = len(graph.nodes)
    rqy_ids, relocation_rows = _relocate_suffix(graph, rqy_rows, BASE_NODE_COUNT, rqy_start)
    _require(len(rqy_ids) == 102 and rqy_ids[0] == rqy_start, "RQy relocation failed")
    seam_endpoints = rqy_patch["endpoint_handles"]
    remapped_ports = {
        name: _map_handles(values, BASE_NODE_COUNT, rqy_start, len(rqy_rows))
        for name, values in seam_endpoints.items()
    }
    _require(remapped_ports["second_micro_upper"] == list(seam_endpoints["second_micro_upper"]),
             "Q4 second-micro endpoint was incorrectly relocated")
    _require(remapped_ports["frozen_macro_lower"] == rqy_ids[4:8], "RQy macro lower target changed")
    _require(remapped_ports["signed_cancellation"] == [rqy_start + i for i in (41, 59, 77, 95)],
             "RQy signed cancellation targets changed")

    macro = _macro_source(graph, remapped_ports, q4, rqy_start)
    macro_suffix_start = macro["graph_suffix_append_index"]
    macro_suffix = copy.deepcopy(graph.appended[macro_suffix_start:])
    _require(macro_suffix == macro["appended_nodes"], "macro suffix receipt changed")
    join_receipt["endpoint_proof_source"] = {
        "micro_second_upper": copy.deepcopy(correlated_patch["endpoint_fields"]["second_micro"]["upper"]),
        "macro_lower": copy.deepcopy(macro["endpoint_fields"]["lower"]),
    }

    micro_patch = correlated_patch
    micro_pair = tuple(micro_patch["actual_micro_pair_increment"])
    _require(len(micro_pair) == 2 and all(isinstance(value, int) for value in micro_pair),
             "accepted micro pair source missing")
    micro_upper = tuple(micro_patch["endpoint_fields"]["second_micro"]["upper"]["R_Q_y"][:2])
    seam_upper_full = tuple(remapped_ports["second_micro_upper"])
    seam_upper = seam_upper_full[:2]
    endpoint_wrapper_identity: list[dict[str, Any]] = []
    for order, (correlated_handle, seam_handle) in enumerate(zip(micro_upper, seam_upper)):
        _require(graph.nodes[correlated_handle] == graph.nodes[seam_handle],
                 f"second-micro R_Q_y wrapper mismatch at ordinary Z{order}")
        endpoint_wrapper_identity.append({
            "ordinary_Z_order": order,
            "correlated_handle": correlated_handle,
            "seam_handle": seam_handle,
            "row": copy.deepcopy(graph.nodes[seam_handle]),
            "same_expression_variable_value_and_Z_independence": True,
        })
    micro_B = tuple(micro_patch["reduced_boundary"])
    micro_D = tuple(micro_patch["reduced_volume"])
    micro_S1 = tuple(micro_patch["reduced_S1"])
    for order, handle in enumerate(micro_S1):
        _require(graph.nodes[handle] == {
            "operation": "sum",
            "arguments": [micro_B[order], micro_D[order]],
        }, f"accepted micro B+D seam row changed at ordinary Z{order}")
    macro_B = tuple(macro["reduced_boundary_B"])
    macro_D = tuple(macro["reduced_volume_D"])
    macro_S1 = tuple(macro["reduced_S1"])
    for order, handle in enumerate(macro_S1):
        _require(graph.nodes[handle] == {
            "operation": "sum",
            "arguments": [macro_B[order], macro_D[order]],
        }, f"macro B+D seam row changed at ordinary Z{order}")
    first_lower = tuple(micro_patch["endpoint_fields"]["first_micro"]["lower"]["R_Q_y"][:2])
    macro_upper = tuple(macro["outer_RQy"]["upper"][:2])
    micro_increment_proof = _verify_increment_graph_rows(
        graph.nodes,
        micro_pair,
        micro_S1,
        tuple(micro_patch["endpoint_fields"]["second_micro"]["upper"]["R_Q_y"][:2]),
        first_lower,
        "accepted_micro_pair",
    )
    macro_increment_proof = _verify_increment_graph_rows(
        graph.nodes,
        tuple(macro["macro_increment"]),
        macro_S1,
        macro_upper,
        tuple(macro["endpoint_fields"]["lower"]["R_Q_y"][:2]),
        "new_macro_increment",
    )
    internal_B_seam = _prove_internal_B_seam(
        graph.nodes,
        micro_patch["endpoint_fields"]["second_micro"]["upper"],
        macro["endpoint_fields"]["lower"],
        join_receipt,
    )
    joined_append_start = len(graph.appended)
    first_lower_B_fields = micro_patch["endpoint_fields"]["first_micro"]["lower"]
    macro_upper_B_fields = macro["endpoint_fields"]["upper"]
    first_lower_B_wrappers = _verify_endpoint_wrappers(
        graph.nodes, first_lower_B_fields, 66988, 97295, "first-micro lower B",
    )
    macro_upper_B_wrappers = _verify_endpoint_wrappers(
        graph.nodes, macro_upper_B_fields, MACRO_COORDINATE, 1, "macro upper B",
    )
    outer_B_boundary = {
        key: graph.padd(
            tuple(macro_upper_B_fields[key]),
            graph.pneg(tuple(first_lower_B_fields[key])),
        )
        for key in ("R_squared_q", "R_squared_q_squared", "R_squared_U_squared")
    }
    formula_context = macro["formula_context"]
    z = tuple(formula_context["z"])
    L = tuple(formula_context["L"])
    L2 = tuple(formula_context["L2"])
    two_id = int(formula_context["two_id"])
    divide_outer = lambda pair, denominator=L: graph.pdiv(pair, denominator, "same positive physical L=1-delta*Z^2")
    half_outer_field = lambda pair: graph.pscale(pair, graph.constant("1/2"))
    outer_B = graph.padd(
        divide_outer(outer_B_boundary["R_squared_q"]),
        half_outer_field(outer_B_boundary["R_squared_q_squared"]),
        graph.pneg(
            graph.pscale(
                divide_outer(graph.pmul(z, z, outer_B_boundary["R_squared_U_squared"]), L2),
                two_id,
            )
        ),
    )
    outer_D = graph.padd(micro_D, macro_D)
    outer_S1 = graph.padd(outer_B, outer_D)
    half_outer = graph.pscale(outer_S1, graph.constant("1/2"))
    uncontracted_joined = graph.padd(micro_pair, macro["macro_increment"])
    shifted_tail = graph.padd(half_outer, graph.pneg(macro_upper))
    no_seam_joined = graph.padd(shifted_tail, first_lower)
    joined_suffix = copy.deepcopy(graph.appended[joined_append_start:])
    contraction = remapped_ports["signed_cancellation"]
    contraction_rows: list[dict[str, Any]] = []
    for order, handle in enumerate(contraction):
        negative_handle = rqy_start + 40 + 18 * order
        lower_handle = rqy_start + 4 + order
        _require(graph.nodes[negative_handle] == {
            "operation": "negative",
            "argument": seam_upper_full[order],
        }, f"RQy negative endpoint row changed at ordinary Z{order}")
        _require(graph.nodes[handle] == {
            "operation": "sum",
            "arguments": [negative_handle, lower_handle],
        }, f"RQy signed contraction row changed at ordinary Z{order}")
        contraction_rows.append({
            "ordinary_Z_order": order,
            "handle": handle,
            "negative_second_micro_upper": seam_upper_full[order],
            "frozen_macro_lower": lower_handle,
            "row": copy.deepcopy(graph.nodes[handle]),
        })
    rqy_identity_rows: list[dict[str, Any]] = []
    for order, proof in enumerate(seam["proof_rows"]):
        micro_identity = proof["factorized_micro_endpoint_identity"]
        macro_identity = proof["factorized_macro_lower_identity"]
        conditional = proof["conditional_boundary_identity"]
        rqy_identity_rows.append({
            "ordinary_Z_order": order,
            "micro_wrapper": micro_identity["actual_wrapper"],
            "macro_wrapper_before_relocation": macro_identity["actual_wrapper"],
            "macro_wrapper_after_relocation": remapped_ports["frozen_macro_lower"][order],
            "signed_contraction_row": contraction_rows[order],
            "factorized_micro_endpoint_identity": copy.deepcopy(micro_identity),
            "factorized_macro_lower_identity": copy.deepcopy(macro_identity),
            "conditional_boundary_identity": copy.deepcopy(conditional),
        })
    rqy_identity = {
        "passed": len(rqy_identity_rows) == 4,
        "rows": rqy_identity_rows,
        "micro_wrappers_are_exact_base_ports": [row["micro_wrapper"] for row in rqy_identity_rows]
        == [133745, 133746, 133747, 133748],
        "macro_wrappers_are_exact_relocated_ports": [row["macro_wrapper_after_relocation"] for row in rqy_identity_rows]
        == list(remapped_ports["frozen_macro_lower"]),
    }
    rqy_identity["passed"] = bool(
        rqy_identity["passed"]
        and rqy_identity["micro_wrappers_are_exact_base_ports"]
        and rqy_identity["macro_wrappers_are_exact_relocated_ports"]
    )
    joined_expression = {
        "ordinary_Z_orders": [0, 1],
        "micro_whole_pair": list(micro_pair),
        "macro_increment": list(macro["macro_increment"]),
        "uncontracted_joined": list(uncontracted_joined),
        "joined": list(no_seam_joined),
        "no_seam_joined": list(no_seam_joined),
        "formula": "(B_outer+D_micro+D_macro)/2 - R_Q_y(m=1) + R_Q_y(first_micro lower)",
        "contracted_formula": "(outer_B + D_micro + D_macro)/2 - R_Q_y(m=1) + R_Q_y(first_micro lower)",
        "reduced_boundary_B": {
            "micro": list(micro_B),
            "macro": list(macro_B),
            "outer": list(outer_B),
            "outer_boundary_fields": {key: list(value) for key, value in outer_B_boundary.items()},
            "outer_endpoint_wrappers": {
                "first_micro_lower": first_lower_B_wrappers,
                "macro_upper": macro_upper_B_wrappers,
            },
        },
        "reduced_volume_D": {"micro": list(micro_D), "macro": list(macro_D), "outer": list(outer_D)},
        "shared_U_Q_B_seam": {
            "micro_B_plus_D": list(micro_S1),
            "macro_B_plus_D": list(macro_S1),
            "outer_B_plus_D": list(outer_S1),
            "same_source_physical_U_Q": internal_B_seam["passed"],
            "internal_endpoint_identity": internal_B_seam,
            "ordinary_Z_orders": [0, 1],
        },
        "uncontracted_equals_contracted": {
            "passed": bool(
                micro_increment_proof["passed"]
                and macro_increment_proof["passed"]
                and internal_B_seam["passed"]
                and rqy_identity["passed"]
            ),
            "uncontracted": list(uncontracted_joined),
            "contracted": list(no_seam_joined),
            "accepted_micro_pair_decomposition": micro_increment_proof,
            "new_macro_increment_decomposition": macro_increment_proof,
            "RQy_identity": rqy_identity,
            "B_internal_identity": internal_B_seam["canonical_digest"],
            "algebra": "micro_pair + macro_increment = half(B_outer+D_micro+D_macro)-R_Q_y(m=1)+R_Q_y(first_micro lower)",
        },
        "rqy_boundary_contraction": contraction,
        "rqy_boundary_contraction_rows": contraction_rows,
        "rqy_signed_terms": "-R_Q_y(second_micro,s=2)+R_Q_y(frozen_macro,m=0)",
        "contraction_only_in_joined_expression": True,
        "latest_second_micro_increment_not_added": True,
        "first_micro_lower_endpoint_retained": True,
        "macro_outer_RQy_retained": True,
        "joined_shifted_tail": list(shifted_tail),
        "shifted_tail_formula": ".5*(B_outer+D_micro+D_macro)-R_Q_y(m=1)",
        "pressure_increment": list(no_seam_joined),
        "absolute_M15_formula": "M15(m=1)=accepted_H+joined_shifted_tail",
        "pressure_increment_scope": "first_micro_sc/2_to_macro_m=1",
        "pressure_increment_source_defined": True,
        "pressure_increment_numerical_enclosure": False,
        "original_first_lower_M15_retained": True,
        "shifted_H_anchor_retained": True,
        "shifted_H_anchor_source": "accepted_correlated.shifted_inlet_source_digest",
        "shifted_H_anchor_digest": correlated["shifted_inlet_source_digest"],
        "endpoint_wrapper_identity": endpoint_wrapper_identity,
    }
    relocated_rows = graph.appended[
        len(correlated_patch["appended_nodes"]):
        len(correlated_patch["appended_nodes"]) + len(rqy_rows)
    ]
    accepted_rows = copy.deepcopy(correlated_patch["appended_nodes"])
    relocated_rows_copy = copy.deepcopy(relocated_rows)
    macro_rows = copy.deepcopy(macro_suffix)
    joined_rows = copy.deepcopy(joined_suffix)
    merged_suffix = accepted_rows + relocated_rows_copy + macro_rows + joined_rows
    _require(len(base_nodes) + len(merged_suffix) == len(graph.nodes),
             "serialized graph suffix does not cover the merged namespace")
    _require(graph.nodes[BASE_NODE_COUNT:] == merged_suffix,
             "serialized graph suffix differs from the constructed namespace")
    segment_starts = {
        "accepted_correlated": BASE_NODE_COUNT,
        "relocated_rqy": BASE_NODE_COUNT + len(accepted_rows),
        "new_macro": BASE_NODE_COUNT + len(accepted_rows) + len(relocated_rows_copy),
        "joined_expression": BASE_NODE_COUNT + len(accepted_rows) + len(relocated_rows_copy) + len(macro_rows),
    }
    graph_start_digest = _digest({"base_sha256": Q4_SHA256, "base_node_count": BASE_NODE_COUNT,
                                  "correlated_suffix": correlated_patch["appended_nodes"],
                                  "relocated_rqy_suffix": relocated_rows_copy, "macro_suffix": macro_rows,
                                  "joined_suffix": joined_rows})
    source_pins[PARENT] = {"path": PARENT, "expected_sha256": PARENT_SHA256, "passed": True}
    source_pins[CORRELATED] = {"path": CORRELATED, "expected_sha256": CORRELATED_SHA256, "passed": True}
    source_guards = {
        "accepted_correlated_source": {
            "passed": True,
            "standalone_equals_latest_parent_packet": True,
            "packet_digest": _digest(correlated),
            "exact_algebra_receipt": copy.deepcopy(correlated["exact_algebra"]),
            "physical_measure": correlated["physical_I1_measure"],
            "accepted_producer_executed": False,
        },
        "accepted_rqy_join": copy.deepcopy(seam["source_guards"]),
        "accepted_rqy_proof_rows": copy.deepcopy(seam["proof_rows"]),
        "join_physical_source_binding": copy.deepcopy(join_receipt),
        "operation_specific_relocation": {
            "passed": True,
            "reference_fields": {key: list(value) for key, value in _REF_FIELDS.items()},
            "exact_rational_numerator_denominator_are_literals": True,
            "definite_integral_variable_strings_preserved": True,
            "graph_metadata_integers_not_relocated": True,
            "source_graph_prefix_unchanged": graph.nodes[:BASE_NODE_COUNT] == base_snapshot,
        },
        "macro_radius_and_length": {
            "passed": True,
            "radius_Z_independent": True,
            "length_Z_independent": True,
            "physical_R_squared_length_applied_once": True,
        },
        "second_micro_endpoint_wrapper_identity": {
            "passed": True,
            "rows": endpoint_wrapper_identity,
            "same_expression_variable_value_and_Z_independence": True,
        },
        "rqy_boundary_contraction": {
            "passed": True,
            "rows": contraction_rows,
            "opposite_FTC_signs": True,
        },
        "reduced_B_seam": {
            "passed": internal_B_seam["passed"],
            "micro_B_plus_D_rows": list(micro_S1),
            "macro_B_plus_D_rows": list(macro_S1),
            "shared_physical_U_Q_source": internal_B_seam["passed"],
            "internal_endpoint_identity": internal_B_seam,
            "outer_B_boundary_fields": {key: list(value) for key, value in outer_B_boundary.items()},
            "outer_endpoint_wrappers": {
                "first_micro_lower": first_lower_B_wrappers,
                "macro_upper": macro_upper_B_wrappers,
            },
            "ordinary_Z_orders": [0, 1],
        },
    }
    safe_graphs = {
        "base_q4": {"archive": Q4_ARCHIVE, "sha256": Q4_SHA256, "node_count": BASE_NODE_COUNT},
        "accepted_correlated": {
            "base_node_count": BASE_NODE_COUNT,
            "appended_node_count": len(correlated_patch["appended_nodes"]),
            "appended_nodes": copy.deepcopy(correlated_patch["appended_nodes"]),
            "source_handles": {"actual_micro_pair_increment": list(micro_pair)},
        },
        "relocated_rqy": {
            "old_base_node_count": BASE_NODE_COUNT,
            "new_base_node_count": correlated_suffix_start + len(correlated_patch["appended_nodes"]),
            "appended_node_count": len(relocation_rows),
            "appended_nodes": copy.deepcopy(
                [dict(row) for row in relocated_rows]
            ),
            "relocation_receipt": relocation_rows,
            "target_ports": remapped_ports,
        },
        "new_macro": {
            "base_node_count": rqy_start + len(rqy_rows),
            "appended_node_count": len(macro_suffix),
            "appended_nodes": macro_suffix,
            "source_handles": {"six_I1": macro["integrals"], "increment": macro["macro_increment"]},
        },
        "joined_expression": {
            "base_node_count": joined_append_start + BASE_NODE_COUNT,
            "appended_node_count": len(joined_suffix),
            "appended_nodes": joined_suffix,
            "source_handles": {
                "uncontracted": list(uncontracted_joined),
                "no_seam": list(no_seam_joined),
                "shifted_tail": list(shifted_tail),
                "pressure_increment": list(no_seam_joined),
            },
        },
    }
    return {
        "passed": True,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_pins": source_pins,
        "input_pins": {
            CORRELATED: {"path": CORRELATED, "expected_sha256": CORRELATED_SHA256, "passed": True},
            PARENT: {"path": PARENT, "expected_sha256": PARENT_SHA256, "passed": True},
        },
        "source_guards": source_guards,
        "graph_namespace": {
            "base_archive": Q4_ARCHIVE,
            "base_node_count": BASE_NODE_COUNT,
            "accepted_correlated_suffix_count": len(correlated_patch["appended_nodes"]),
            "relocated_rqy_suffix_count": len(rqy_rows),
            "new_macro_suffix_count": len(macro_suffix),
            "joined_expression_suffix_count": len(joined_suffix),
            "merged_node_count": len(graph.nodes),
            "old_graph_unchanged": graph.nodes[:BASE_NODE_COUNT] == base_snapshot,
            "namespace_digest": graph_start_digest,
            "serialized_suffix_count": len(merged_suffix),
            "segment_starts": segment_starts,
            "serialized_suffix_covers_every_new_node": True,
        },
        "safe_graphs": safe_graphs,
        "merged_namespace": {
            "base_archive": Q4_ARCHIVE,
            "base_node_count": BASE_NODE_COUNT,
            "appended_nodes": merged_suffix,
            "segment_starts": segment_starts,
            "merged_node_count": len(graph.nodes),
            "roundtrip_verified": True,
        },
        "remapped_accepted_seam_target_ports": remapped_ports,
        "macro_six_I1_source": macro,
        "actual_macro_six_I1_source": macro,
        "remaining_outer_macro_RQy": macro["outer_RQy"],
        "actual_joined_pressure_function": joined_expression,
        "joined_pressure_function": joined_expression,
        "accepted_rqy_seam_identity_consumed": {
            "signed_cancellation_rows": contraction,
            "signed_cancellation_row_checks": contraction_rows,
            "accepted_full_Q_y_rows": copy.deepcopy(seam["accepted_full_Q_y_rows"]),
            "parent_proof_rows": copy.deepcopy(seam["proof_rows"]),
            "second_micro_wrapper_identity": endpoint_wrapper_identity,
            "common_physical_source": join_receipt["passed"],
            "R_Z_independent": seam["source_guards"]["radius_Z_independent_transitive"],
            "join_physical_source_binding": copy.deepcopy(join_receipt),
            "identity_used_only_for_joined_expression": True,
        },
        "old_212_terms_and_errors_unchanged": True,
        "old_six_numeric_terms_and_all_error_caps_unchanged": True,
        "signed_function_values_evaluated": False,
        "signed_six_integrals_evaluated": False,
        "global_goal_complete": False,
        "accepted_ancestors_replayed": False,
        "numerical_phase_or_quadrature_evaluated": False,
    }


__all__ = ["build"]
