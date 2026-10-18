"""Source-bound VALUE seam for physical ``Q_y`` at ``R_m``.

This consumer is deliberately smaller than the older slow native seam.  It
does not run a provider, replay a saved operator, differentiate an opaque
VALUE atom, or evaluate a numerical field.  It validates the pinned physical
ports and the source formula for ordinary ``Z0`` and then rewrites that
bounded DAG with formal VALUE atoms supplied by the admitted loop, modified
moment, and modified radial receipts.

The only physical pressure claim made here is the source identity
``Q_y[Z0]``.  The ``Q_y[Z1]`` targets are retained as open handles and are
never accepted as a result of this module.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from copy import deepcopy
import hashlib
import json
from typing import Any

try:
    import sympy as _sp
except ImportError:  # pragma: no cover - callers receive a closed failure
    _sp = None


ORDERS = (0, 1)
SIDES = ("postrestore", "patch")
SCHEMA = "original_Rm_physical_qy0_value_seam_v1"
MOMENT_FIELDS = ("m", "h", "k", "e", "p")
LOOP_FIELDS = (
    "A",
    "B_over_Pstar",
    "A_phi",
    "B_over_Pstar_phi",
    "E_N",
    "V_N",
    "m",
    "h",
    "k",
    "e",
    "p",
)
RADIAL_FIELDS = ("E_y", "V_y")

EXPECTED_GRAPH_NODE_COUNT = 205300
EXPECTED_GRAPH_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_FINITE_PREFIX_COUNT = 205147
EXPECTED_FINITE_PREFIX_DIGEST = "c7ab7da8cc5c51cd782a5794e3076a6cf021a45235aae1c75bd246b01491b639"
EXPECTED_LOOP_GRAPH_COUNT = 205200
EXPECTED_LOOP_GRAPH_DIGEST = "6b6f1a9461bbc4bb33f25f8a1153bbda96b360fb09af19bb57dc6e521d078960"
EXPECTED_MOMENT_GRAPH_COUNT = 205240
EXPECTED_MOMENT_GRAPH_DIGEST = "ff5c2c5407a1a38ed5c20e9a3aec338946ebae9533f813bdc1215885c0df3234"
FINITE_ARCHIVE_SHA256 = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
LOOP_PORTS_SHA256 = "0331c5fc37bd3ec6c0d386c38b07d9863ede622e14f1c6dc65bf8cdf675a16b9"
ACCEPTED_REPORT_SHA256 = "e649b4f1556c3c04eadbf9469d79e14795d6be2da8a8db52669a9e5fdebc05ea"
PRODUCER_MODULE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
)
PRODUCER_MODULE_SHA256 = "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"

EXTERNAL_MEAN_PREMISE_SHA256 = "3f936cbb5f91fd05677f0c292ed783ad72788b0971ed869298f43be9dfc0850d"
EXTERNAL_MEAN_CERTIFICATE_DIGEST = "fc315269a1cca8bda01f7a47b876ef1f7abef4d25fe85c9ce6cf8df64061530e"

PSTAR_NODE = 33572
Z_NODE = 20
POST_COORDINATE = 84528
PATCH_COORDINATE = 185676
POST_RADIUS = 84533
PATCH_RADIUS = 187214
POST_J = 1
PATCH_J = 187215
COMMON_L_NODE = 25272
COMMON_DELTA_ONE_NODE = 25222
COMMON_D_NODE = 25251


EXPECTED_LOOP_PORTS: dict[str, dict[tuple[str, int], tuple[int, int]]] = {
    "postrestore": {
        ("A", 0): (121070, 197890),
        ("A", 1): (121071, 197891),
        ("B_over_Pstar", 0): (121108, 197928),
        ("B_over_Pstar", 1): (121109, 197929),
        ("A_phi", 0): (121189, 198009),
        ("A_phi", 1): (121190, 198010),
        ("B_over_Pstar_phi", 0): (121291, 198123),
        ("B_over_Pstar_phi", 1): (121292, 198124),
        ("E_N", 0): (121560, 198412),
        ("E_N", 1): (121562, 198414),
        ("V_N", 0): (121571, 198423),
        ("V_N", 1): (121572, 198424),
        ("m", 0): (121350, 198182),
        ("m", 1): (121351, 198183),
        ("h", 0): (121324, 198156),
        ("h", 1): (121327, 198159),
        ("k", 0): (121424, 198266),
        ("k", 1): (121425, 198267),
        ("e", 0): (121554, 198406),
        ("e", 1): (121555, 198407),
        ("p", 0): (121542, 198394),
        ("p", 1): (121543, 198395),
    }
}
EXPECTED_LOOP_PORTS["patch"] = {
    key: pair for key, pair in EXPECTED_LOOP_PORTS["postrestore"].items()
}

# The saved rows use the first element for postrestore and the second for
# patch.  Keeping the table explicit makes accidental chart interchange a
# hard failure rather than an ancestry or containment match.
EXPECTED_MOMENT_PORTS: dict[str, dict[tuple[str, int], int]] = {
    "postrestore": {
        ("m", 0): 121720,
        ("m", 1): 121721,
        ("h", 0): 121726,
        ("h", 1): 121727,
        ("k", 0): 121732,
        ("k", 1): 121733,
        ("e", 0): 121738,
        ("e", 1): 121739,
        ("p", 0): 121744,
        ("p", 1): 121745,
    },
    "patch": {
        ("m", 0): 198591,
        ("m", 1): 198592,
        ("h", 0): 198597,
        ("h", 1): 198598,
        ("k", 0): 198603,
        ("k", 1): 198604,
        ("e", 0): 198609,
        ("e", 1): 198610,
        ("p", 0): 198615,
        ("p", 1): 198616,
    },
}

EXPECTED_RADIAL_PORTS: dict[str, dict[str, tuple[int, int]]] = {
    "postrestore": {
        "E_y": (182725, 182728),
        "V_y": (182757, 182758),
    },
    "patch": {
        "E_y": (203774, 203777),
        "V_y": (203806, 203807),
    },
}

EXPECTED_PHYSICAL_PORTS: dict[str, dict[str, tuple[int, int]]] = {
    "postrestore": {
        "U": (183268, 183269),
        "U_y": (183274, 183275),
        "m": (183280, 183281),
        "m_y": (183292, 183293),
        "Q_y": (183451, 183455),
    },
    "patch": {
        "U": (204317, 204318),
        "U_y": (204323, 204324),
        "m": (204329, 204330),
        "m_y": (204341, 204342),
        "Q_y": (204500, 204504),
    },
}

EXPECTED_SOURCE_BINDINGS: dict[str, Any] = {
    "postrestore": {
        "chart": {
            "coordinate": POST_COORDINATE,
            "endpoint": 1,
            "radius": POST_RADIUS,
            "physical_measure": POST_J,
            "Pstar": PSTAR_NODE,
        },
        "radial": {key: list(value) for key, value in EXPECTED_RADIAL_PORTS["postrestore"].items()},
        "physical": {key: list(value) for key, value in EXPECTED_PHYSICAL_PORTS["postrestore"].items()},
        "total_V_y": list(EXPECTED_RADIAL_PORTS["postrestore"]["V_y"]),
    },
    "patch": {
        "chart": {
            "coordinate": PATCH_COORDINATE,
            "endpoint": 1,
            "radius": PATCH_RADIUS,
            "physical_measure": PATCH_J,
            "Pstar": PSTAR_NODE,
        },
        "radial": {key: list(value) for key, value in EXPECTED_RADIAL_PORTS["patch"].items()},
        "physical": {key: list(value) for key, value in EXPECTED_PHYSICAL_PORTS["patch"].items()},
        "total_V_y": list(EXPECTED_RADIAL_PORTS["patch"]["V_y"]),
    },
}

_COMMON_SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


class PhysicalQyValueSeamError(ValueError):
    """Malformed admission data or an unresolved source VALUE closure."""


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise PhysicalQyValueSeamError(message)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def canonical_digest(value: Any) -> str:
    """Return the stable JSON digest used by saved source receipts."""

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _digest_hex(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, f"{label} must be a SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise PhysicalQyValueSeamError(f"{label} must be a SHA-256 digest") from exc
    return value.lower()


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), f"{label} must be a mapping")
    return value


def _handle(value: Any, label: str, count: int) -> int:
    if isinstance(value, Mapping):
        value = value.get("node", value.get("handle", value.get("root")))
    _require(isinstance(value, int) and not isinstance(value, bool), f"{label} must be an integer graph handle")
    _require(0 <= value < count, f"{label} is outside the source graph")
    return value


def _row(nodes: Sequence[Mapping[str, Any]], index: int, label: str = "node") -> Mapping[str, Any]:
    i = _handle(index, label, len(nodes))
    value = nodes[i]
    _require(isinstance(value, Mapping), f"{label} {i} is not a mapping")
    return value


def _operation(nodes: Sequence[Mapping[str, Any]], index: int, label: str = "node") -> str:
    value = _row(nodes, index, label).get("operation", _row(nodes, index, label).get("op"))
    _require(isinstance(value, str) and value, f"{label} {index} has no operation")
    return value


def _exact(node: Mapping[str, Any], numerator: int, denominator: int = 1) -> bool:
    return (
        node.get("operation", node.get("op")) == "exact_rational"
        and node.get("numerator", node.get("value")) == numerator
        and node.get("denominator", 1) == denominator
    )


def _form(value: Any, label: str) -> dict[str, Any]:
    value = _mapping(value, label)
    srepr = value.get("srepr")
    text = value.get("text", srepr)
    is_zero = value.get("is_zero")
    _require(isinstance(srepr, str) and srepr, f"{label}.srepr must be nonempty text")
    _require(isinstance(text, str) and text, f"{label}.text must be nonempty text")
    _require(isinstance(is_zero, bool), f"{label}.is_zero must be boolean")
    if is_zero:
        _require(srepr in {"Integer(0)", "0"}, f"{label} claims a non-exact zero")
    else:
        _require(srepr not in {"Integer(0)", "0"}, f"{label} has a zero form without zero admission")
    return {"srepr": srepr, "text": text, "is_zero": is_zero}


def _side(raw: Any, side: str, label: str) -> Mapping[str, Any]:
    raw = _mapping(raw, label)
    if side in raw:
        return _mapping(raw[side], f"{label}.{side}")
    aliases = ("left", "post") if side == "postrestore" else ("right",)
    for name in aliases:
        if name in raw:
            return _mapping(raw[name], f"{label}.{name}")
    raise PhysicalQyValueSeamError(f"{label} lacks {side}")


def _pair(value: Any, label: str) -> tuple[Any, Any]:
    if isinstance(value, Mapping):
        if "postrestore" in value or "patch" in value:
            return value.get("postrestore"), value.get("patch")
        if "0" in value or "1" in value:
            return value.get("0"), value.get("1")
    if _is_sequence(value) and len(value) == 2:
        return value[0], value[1]
    raise PhysicalQyValueSeamError(f"{label} must contain two side values")


def _port_pair(value: Any, label: str) -> tuple[Any, Any]:
    """Read the ordinary ``Z0``/``Z1`` prefix of a source port list.

    Archive bindings carry six chart entries for most fields (and five for
    physical ``Q``/``Q_y``), while this seam consumes only the two ordinary
    orders.  A side-local binding is therefore allowed to contain a longer
    source list, but its first two entries remain mandatory and are never
    confused with the paired ``postrestore``/``patch`` receipt shape used by
    admission packets.
    """

    if isinstance(value, Mapping):
        if "ordinary_order_0" in value or "ordinary_order_1" in value:
            return value.get("ordinary_order_0"), value.get("ordinary_order_1")
        if "0" in value or "1" in value:
            return value.get("0"), value.get("1")
    if _is_sequence(value) and len(value) >= 2:
        return value[0], value[1]
    raise PhysicalQyValueSeamError(f"{label} must contain ordinary Z0 and Z1 ports")


def _row_order(row: Mapping[str, Any], label: str) -> tuple[str, int]:
    field = row.get("field")
    order = row.get("order", row.get("ordinary_order", row.get("ordinary_Z_order")))
    _require(isinstance(field, str) and field, f"{label}.field is missing")
    _require(isinstance(order, int) and not isinstance(order, bool) and order in ORDERS, f"{label}.order is invalid")
    return field, order


def _owner_node(owner: Mapping[str, Any], label: str, count: int) -> int:
    return _handle(owner.get("node", owner.get("root_node", owner.get("target_handle"))), label, count)


def _validate_graph_digest(nodes: Sequence[Mapping[str, Any]], expected_nodes_digest: str) -> dict[str, str]:
    _require(_is_sequence(nodes), "nodes must be a sequence")
    _require(len(nodes) == EXPECTED_GRAPH_NODE_COUNT, f"source graph must contain {EXPECTED_GRAPH_NODE_COUNT} nodes")
    _require(all(isinstance(value, Mapping) for value in nodes), "nodes must contain mappings")
    actual = canonical_digest(list(nodes))
    _require(actual == expected_nodes_digest, "nodes digest does not match expected_nodes_digest")
    _require(actual == EXPECTED_GRAPH_DIGEST, "source graph digest is not the accepted Rm graph")
    prefixes = {
        "finite_prefix": (EXPECTED_FINITE_PREFIX_COUNT, EXPECTED_FINITE_PREFIX_DIGEST),
        "loop_prefix": (EXPECTED_LOOP_GRAPH_COUNT, EXPECTED_LOOP_GRAPH_DIGEST),
        "moment_prefix": (EXPECTED_MOMENT_GRAPH_COUNT, EXPECTED_MOMENT_GRAPH_DIGEST),
    }
    actual_prefixes: dict[str, str] = {}
    for name, (count, digest) in prefixes.items():
        prefix = canonical_digest(list(nodes[:count]))
        _require(prefix == digest, f"{name} digest is not the accepted source prefix")
        actual_prefixes[name] = prefix
    actual_prefixes["full_graph"] = actual
    return actual_prefixes


def _validate_common_phase(packet: Mapping[str, Any], count: int) -> dict[str, Any]:
    common = packet.get("common_phase_source")
    _require(isinstance(common, Mapping), "bounded_loop_ports_admission.common_phase_source is missing")
    _require(common.get("schema") == "actual_Rm_common_phase_source_v1", "common phase schema changed")
    for key in ("passed", "source_only", "actual_common_phase_proved", "smooth_periodic_lift_bound", "actual_phase_Z_independence_proved"):
        _require(common.get(key) is True, f"common phase {key} must be true")
    for key in ("fractional_part_differentiated", "numerical_phase_value_computed", "native_loop_derivative_proved", "density_seam_proved", "inverse_conditioning_proved"):
        _require(common.get(key) is False, f"common phase {key} must be false")
    _require(common.get("phase_rows") == [
        {
            "coordinate": POST_COORDINATE,
            "endpoint": 1,
            "radius": 84533,
            "global_offset": 119414,
            "actual_phase": 119416,
            "common_N": 18,
            "common_Ra": 65604,
            "common_phase_origin": 24954,
        },
        {
            "coordinate": PATCH_COORDINATE,
            "endpoint": 1,
            "radius": PATCH_RADIUS,
            "global_offset": 196044,
            "actual_phase": 196046,
            "common_N": 18,
            "common_Ra": 65604,
            "common_phase_origin": 24954,
        },
    ], "common phase rows changed")
    _require(common.get("phase_nodes") == [119416, 196046], "common phase handles changed")
    _require(common.get("coordinates") == [POST_COORDINATE, PATCH_COORDINATE], "common phase coordinates changed")
    _require(common.get("endpoints") == [1, 1], "common phase endpoints must be exact one")
    return {
        "phase_nodes": [119416, 196046],
        "coordinates": [POST_COORDINATE, PATCH_COORDINATE],
        "endpoints": [1, 1],
    }


def _validate_loop_ports(packet: Mapping[str, Any], count: int) -> dict[tuple[str, int], dict[str, Any]]:
    _require(packet.get("schema") == "admitted_Rm_loop_source_ports_v1", "loop_ports_admission schema changed")
    _require(packet.get("passed") is True, "loop_ports_admission.passed must be true")
    _require(packet.get("source_only") is True, "loop_ports_admission must be source-only")
    _require(packet.get("mathematical_operator_executed") is False, "loop_ports_admission executed an operator")
    _require(packet.get("numerical_values_computed") is False, "loop_ports_admission computed numerical values")
    _require(packet.get("native_loop_seam_proved") is False, "loop_ports_admission overclaims the native loop seam")
    _require(packet.get("full_goal_complete") is False, "loop_ports_admission overclaims the full goal")
    _require(packet.get("graph_node_count") == EXPECTED_LOOP_GRAPH_COUNT, "loop_ports_admission graph count changed")
    _require(packet.get("ordinary_orders") == [0, 1], "loop_ports_admission ordinary orders changed")
    _validate_common_phase(packet, count)
    rows = packet.get("rows")
    _require(_is_sequence(rows) and len(rows) == len(LOOP_FIELDS) * len(ORDERS), "loop_ports_admission must contain 22 rows")
    accepted: dict[tuple[str, int], dict[str, Any]] = {}
    seen: set[tuple[str, int]] = set()
    for offset, raw in enumerate(rows):
        row = _mapping(raw, f"loop row {offset}")
        field, order = _row_order(row, f"loop row {offset}")
        _require(field in LOOP_FIELDS, f"loop row {offset} has an unsupported field")
        key = (field, order)
        _require(key not in seen, f"duplicate loop row {field}[{order}]")
        seen.add(key)
        _require(row.get("status") == "proved", f"loop row {field}[{order}] is not proved")
        side_data: dict[str, dict[str, Any]] = {}
        for side_index, side in enumerate(SIDES):
            owner = _side(row, side, f"loop row {field}[{order}]")
            node = _owner_node(owner, f"loop {side}.{field}[{order}].node", count)
            expected = EXPECTED_LOOP_PORTS[side][key][side_index]
            _require(node == expected, f"loop {side}.{field}[{order}] handle changed")
            form = _form(owner.get("value"), f"loop {side}.{field}[{order}].value")
            side_data[side] = {"node": node, "value": form}
        _require(side_data["postrestore"]["value"] == side_data["patch"]["value"], f"loop VALUE forms differ for {field}[{order}]")
        accepted[key] = {"field": field, "order": order, **side_data}
    _require(seen == {(field, order) for field in LOOP_FIELDS for order in ORDERS}, "loop rows are incomplete")
    return accepted


def _validate_modified_moment(packet: Mapping[str, Any], count: int) -> dict[tuple[str, int], dict[str, Any]]:
    _require(packet.get("schema") == "actual_Rm_modified_moment_value_native_germ_v1", "modified_moment_germ schema changed")
    for key, expected in (
        ("passed", True),
        ("source_only", True),
        ("actual_saved_graph_run", False),
        ("actual_source_graph_inspected", True),
        ("numeric_values_computed", False),
        ("old_mathematical_operator_replayed", False),
        ("full_goal_complete", False),
        ("moment_value_native_seam_proved", True),
        ("E_V_native_seam_proved", False),
        ("physical_RQy_seam_proved", False),
    ):
        _require(packet.get(key) is expected, f"modified_moment_germ.{key} has an invalid admission flag")
    _require(_digest_hex(packet.get("nodes_digest"), "modified_moment_germ.nodes_digest") == EXPECTED_MOMENT_GRAPH_DIGEST,
             "modified_moment_germ graph digest changed")
    _require(packet.get("ordinary_orders") == [0, 1], "modified_moment_germ ordinary orders changed")
    _require(packet.get("fields") == list(MOMENT_FIELDS), "modified_moment_germ fields changed")
    delta = _mapping(packet.get("graph_delta"), "modified_moment_germ.graph_delta")
    _require(delta.get("base_graph_node_count") == EXPECTED_MOMENT_GRAPH_COUNT, "modified_moment_germ graph delta base count changed")
    _require(_digest_hex(delta.get("base_nodes_sha256"), "modified_moment_germ.graph_delta.base_nodes_sha256") == EXPECTED_MOMENT_GRAPH_DIGEST,
             "modified_moment_germ graph delta base digest changed")
    _require(delta.get("appended_node_count") == 60 and delta.get("full_node_count") == EXPECTED_GRAPH_NODE_COUNT,
             "modified_moment_germ graph delta size changed")
    _require(_digest_hex(delta.get("full_nodes_sha256"), "modified_moment_germ.graph_delta.full_nodes_sha256") == EXPECTED_GRAPH_DIGEST,
             "modified_moment_germ graph delta full digest changed")
    rows = packet.get("rows")
    _require(_is_sequence(rows) and len(rows) == len(MOMENT_FIELDS) * len(ORDERS), "modified_moment_germ must contain 10 rows")
    accepted: dict[tuple[str, int], dict[str, Any]] = {}
    seen: set[tuple[str, int]] = set()
    for offset, raw in enumerate(rows):
        row = _mapping(raw, f"modified moment row {offset}")
        field, order = _row_order(row, f"modified moment row {offset}")
        _require(field in MOMENT_FIELDS, f"modified moment row {offset} has an unsupported field")
        key = (field, order)
        _require(key not in seen, f"duplicate modified moment row {field}[{order}]")
        seen.add(key)
        _require(row.get("status") == "proved", f"modified moment row {field}[{order}] is not proved")
        _require(row.get("value_identity_proved") is True and row.get("native_identity_proved") is True,
                 f"modified moment row {field}[{order}] lacks its identity receipts")
        side_data: dict[str, dict[str, Any]] = {}
        for side in SIDES:
            owner = _side(row, side, f"modified moment row {field}[{order}]")
            body = _handle(owner.get("actual_modified_body"), f"modified moment {side}.{field}[{order}].actual_modified_body", count)
            _require(body == EXPECTED_MOMENT_PORTS[side][key], f"modified moment {side}.{field}[{order}] body changed")
            side_data[side] = {"actual_modified_body": body}
        if field == "m":
            _require(row.get("value_status") == "proved_from_saved_mean_premise", f"modified moment m[{order}] must retain its saved mean premise")
            ref = row.get("leading_mean_proof_ref")
            _require(isinstance(ref, Mapping), f"modified moment m[{order}] mean premise reference is missing")
            _require(ref.get("report_sha256") == EXTERNAL_MEAN_PREMISE_SHA256,
                     f"modified moment m[{order}] mean premise report changed")
            _require(ref.get("certificate_digest") == EXTERNAL_MEAN_CERTIFICATE_DIGEST,
                     f"modified moment m[{order}] mean premise certificate changed")
        else:
            _require(row.get("leading_value_status_preserved") == "proved", f"modified moment {field}[{order}] leading value status changed")
            _require(row.get("leading_mean_proof_ref") is None, f"modified moment {field}[{order}] imported an unrelated mean premise")
        accepted[key] = {
            "field": field,
            "order": order,
            "sides": side_data,
            "value_identity_proved": True,
            "native_identity_proved": True,
            "mean_premise": {
                "label": "external_premise" if field == "m" else None,
                "report_sha256": EXTERNAL_MEAN_PREMISE_SHA256 if field == "m" else None,
                "certificate_digest": EXTERNAL_MEAN_CERTIFICATE_DIGEST if field == "m" else None,
            },
        }
    _require(seen == {(field, order) for field in MOMENT_FIELDS for order in ORDERS}, "modified moment rows are incomplete")
    return accepted


def _value_pair_handle(owner: Mapping[str, Any], label: str, count: int) -> int:
    value = owner.get("target_handle", owner.get("node", owner.get("root_node")))
    return _handle(value, label, count)


def _validate_modified_radial(packet: Mapping[str, Any], count: int) -> dict[tuple[str, int], dict[str, Any]]:
    schema = packet.get("schema")
    schema_text = schema.lower() if isinstance(schema, str) else ""
    _require(
        isinstance(schema, str)
        and "value" in schema_text
        and "seam" in schema_text
        and any(token in schema_text for token in ("radial", "total", "e_y", "v_y")),
             "modified_radial_value_seam schema is missing or unsupported")
    _require(packet.get("passed") is True, "modified_radial_value_seam.passed must be true")
    _require(packet.get("source_only") is True, "modified_radial_value_seam must be source-only")
    for key in (
        "actual_saved_graph_run",
        "numeric_values_computed",
        "numeric_values_evaluated",
        "numeric_field_evaluated",
        "old_mathematical_operator_replayed",
        "old_operators_replayed",
        "new_native_derivative_generated",
        "source_executed",
    ):
        if key in packet:
            _require(packet.get(key) is False, f"modified_radial_value_seam.{key} must be false")
    for key in ("nodes_digest", "graph_digest", "graphdigest", "expected_nodes_digest"):
        if key in packet:
            _require(_digest_hex(packet.get(key), f"modified_radial_value_seam.{key}") == EXPECTED_GRAPH_DIGEST,
                     f"modified_radial_value_seam.{key} changed")
    rows = packet.get("canonical_value_pairs")
    _require(_is_sequence(rows) and len(rows) == len(RADIAL_FIELDS) * len(ORDERS),
             "modified_radial_value_seam must contain four canonical VALUE pairs")
    accepted: dict[tuple[str, int], dict[str, Any]] = {}
    seen: set[tuple[str, int]] = set()
    for offset, raw in enumerate(rows):
        row = _mapping(raw, f"modified radial row {offset}")
        field, order = _row_order(row, f"modified radial row {offset}")
        _require(field in RADIAL_FIELDS, f"modified radial row {offset} has an unsupported field")
        key = (field, order)
        _require(key not in seen, f"duplicate modified radial row {field}[{order}]")
        seen.add(key)
        _require(row.get("equal") is True and row.get("status") in {"proved", "source_identity_proved"},
                 f"modified radial row {field}[{order}] is not passed")
        side_data: dict[str, dict[str, Any]] = {}
        for side_index, side in enumerate(SIDES):
            # Producers may put the handle beside the side VALUE form or at
            # the row level as ``postrestore_target_handle`` /
            # ``patch_target_handle``.  Accept both authenticated shapes;
            # never infer a handle from a VALUE form or from a chart-local
            # unused port.
            owner_value = row.get(side)
            owner = owner_value if isinstance(owner_value, Mapping) else {}
            handle_value = owner.get("target_handle", owner.get("node", owner.get("root_node")))
            if handle_value is None:
                handle_value = row.get(f"{side}_target_handle")
            handle = _handle(handle_value, f"modified radial {side}.{field}[{order}].target_handle", count)
            _require(handle == EXPECTED_RADIAL_PORTS[side][field][order], f"modified radial {side}.{field}[{order}] target changed")
            form_value = owner.get("value", owner.get("form"))
            if form_value is None and "srepr" in owner:
                form_value = owner
            if form_value is None:
                form_value = row.get(f"{side}_value", row.get(f"{side}_value_form"))
            form = _form(form_value, f"modified radial {side}.{field}[{order}].value")
            side_data[side] = {"target_handle": handle, "value": form}
        _require(side_data["postrestore"]["value"] == side_data["patch"]["value"], f"modified radial VALUE forms differ for {field}[{order}]")
        accepted[key] = {"field": field, "order": order, **side_data}
    _require(seen == {(field, order) for field in RADIAL_FIELDS for order in ORDERS}, "modified radial rows are incomplete")
    # A future producer must not smuggle an already-closed Q_y pair through
    # the radial receipt.  Its contract is only E_y and total V_y.
    for row in rows:
        if isinstance(row, Mapping) and row.get("field") == "Q_y":
            raise PhysicalQyValueSeamError("modified radial receipt contains a forbidden Q_y row")
    return accepted


def _lookup_side(raw: Mapping[str, Any], side: str) -> Mapping[str, Any]:
    payload: Any = raw.get("ports", raw)
    payload = _mapping(payload, "actual_source_bindings")
    if side in payload:
        return _mapping(payload[side], f"actual_source_bindings.{side}")
    charts = payload.get("charts")
    if isinstance(charts, Mapping) and side in charts:
        return _mapping(charts[side], f"actual_source_bindings.charts.{side}")
    raise PhysicalQyValueSeamError(f"actual_source_bindings lacks {side}")


def _extract_group(side_data: Mapping[str, Any], name: str, aliases: Sequence[str] = ()) -> Mapping[str, Any]:
    for key in (name, *aliases):
        value = side_data.get(key)
        if isinstance(value, Mapping):
            return value
    raise PhysicalQyValueSeamError(f"actual_source_bindings side lacks {name}")


def _validate_source_bindings(raw: Mapping[str, Any] | None, count: int) -> dict[str, Any]:
    _require(raw is not None, "actual_source_bindings is required for source-bound Q_y congruence")
    supplied = True
    payload = _mapping(raw, "actual_source_bindings")
    if "passed" in payload:
        _require(payload.get("passed") is True, "actual_source_bindings.passed must be true")
    if "source_only" in payload:
        _require(payload.get("source_only") is True, "actual_source_bindings must be source-only")
    for key in ("changed_Z", "changed_z", "changed_native_environment", "native_environment_changed"):
        if key in payload:
            _require(payload.get(key) is False, f"actual_source_bindings.{key} must be false")
    checked: dict[str, Any] = {}
    for side in SIDES:
        side_data = _lookup_side(payload, side)
        for key in ("changed_Z", "changed_z", "changed_native_environment", "native_environment_changed"):
            if key in side_data:
                _require(side_data.get(key) is False, f"{side}.{key} must be false")
        # The modified postrestore archive exposes its chart groups directly
        # under ``charts.postrestore`` and has no nested ``chart`` metadata;
        # the finite patch archive carries that metadata under ``chart``.
        # Use the pinned postrestore constants only for metadata absent from
        # the archive, while still checking every supplied field.
        chart = side_data.get("chart")
        if not isinstance(chart, Mapping):
            chart = {}
        radial = _extract_group(side_data, "radial")
        physical = _extract_group(side_data, "physical")
        for key in ("total_fast_phase_terms_included", "physical_y_equals_native_partial_divided_by_Jacobian"):
            if key in radial:
                _require(radial.get(key) is True, f"{side}.radial.{key} guard changed")
        coordinate_value = chart.get("coordinate", POST_COORDINATE if side == "postrestore" else PATCH_COORDINATE)
        coordinate = _handle(coordinate_value, f"{side}.chart.coordinate", count)
        _require(coordinate == (POST_COORDINATE if side == "postrestore" else PATCH_COORDINATE), f"{side} coordinate changed")
        endpoint = chart.get("endpoint", chart.get("endpoints", 1))
        _require(endpoint == 1, f"{side} endpoint is not exact one")
        radius_value = chart.get("radius", POST_RADIUS if side == "postrestore" else PATCH_RADIUS)
        radius = _handle(radius_value, f"{side}.chart.radius", count)
        _require(radius == (POST_RADIUS if side == "postrestore" else PATCH_RADIUS), f"{side} radius changed")
        measure_value = chart.get("physical_measure", chart.get("J", chart.get("physical_J")))
        if measure_value is None:
            measure_value = POST_J if side == "postrestore" else PATCH_J
        measure = _handle(measure_value, f"{side}.chart.physical_measure", count)
        _require(measure == (POST_J if side == "postrestore" else PATCH_J), f"{side} physical J changed")
        pstar = _handle(chart.get("Pstar", chart.get("pstar", PSTAR_NODE)), f"{side}.chart.Pstar", count)
        _require(pstar == PSTAR_NODE, f"{side} Pstar changed")
        checked_side: dict[str, Any] = {
            "chart": {"coordinate": coordinate, "endpoint": 1, "radius": radius, "physical_measure": measure, "Pstar": pstar},
            "radial": {},
            "physical": {},
        }
        for field in RADIAL_FIELDS:
            values = radial.get(field)
            if values is None and field == "V_y":
                values = side_data.get("total_V_y")
            pair = _port_pair(values, f"{side}.radial.{field}")
            expected = EXPECTED_RADIAL_PORTS[side][field]
            # A side-local source map may give [Z0,Z1], while a full map can
            # give {postrestore,patch}; after selecting the side both must be
            # the exact ordinary-order pair.
            selected = tuple(_handle(value, f"{side}.radial.{field}[{i}]", count) for i, value in enumerate(pair))
            _require(selected == expected, f"{side}.radial.{field} ports changed")
            checked_side["radial"][field] = list(selected)
        for field in EXPECTED_PHYSICAL_PORTS[side]:
            pair = _port_pair(physical.get(field), f"{side}.physical.{field}")
            selected = tuple(_handle(value, f"{side}.physical.{field}[{i}]", count) for i, value in enumerate(pair))
            _require(selected == EXPECTED_PHYSICAL_PORTS[side][field], f"{side}.physical.{field} ports changed")
            checked_side["physical"][field] = list(selected)
        if "m_y" in radial:
            radial_mean_y = tuple(_handle(value, f"{side}.radial.m_y[{i}]", count)
                                  for i, value in enumerate(_port_pair(radial["m_y"], f"{side}.radial.m_y")))
            _require(tuple(checked_side["physical"]["m_y"]) != radial_mean_y,
                     f"{side} physical mean_y must not reuse radial m_y")
        checked[side] = checked_side
    return {"supplied": supplied, "ports": checked, "digest": canonical_digest(checked)}


def _assert_row(nodes: Sequence[Mapping[str, Any]], index: int, operation: str, label: str, **expected: Any) -> Mapping[str, Any]:
    row = _row(nodes, index, label)
    _require(_operation(nodes, index, label) == operation, f"{label} operation changed")
    for key, value in expected.items():
        _require(row.get(key) == value, f"{label}.{key} changed")
    return row


def _validate_common_source_geometry(nodes: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    _require(_exact(_row(nodes, 0), 0), "graph node 0 must remain exact zero")
    _assert_row(nodes, PSTAR_NODE, "analytic_unary", "Pstar", name="exp", argument=5)
    _assert_row(nodes, 5, "sum", "Pstar source", arguments=[3, 4])
    _assert_row(nodes, 3, "analytic_unary", "Pstar exp", name="exp", argument=2)
    _require(_exact(_row(nodes, 2), 40), "Pstar N source changed")
    _require(_exact(_row(nodes, 4), 11), "Pstar source offset changed")
    _assert_row(nodes, COMMON_DELTA_ONE_NODE, "sum", "1-delta", arguments=[1, 25221])
    _assert_row(nodes, 25221, "negative", "delta complement", argument=25220)
    _assert_row(nodes, COMMON_D_NODE, "sum", "1-Z^2", arguments=[1, 25250])
    _assert_row(nodes, 25250, "negative", "Z^2 complement", argument=25249)
    _assert_row(nodes, 25249, "product", "Z^2", arguments=[Z_NODE, Z_NODE])
    _assert_row(nodes, COMMON_L_NODE, "sum", "L", arguments=[1, 25271])
    _assert_row(nodes, 25271, "negative", "L complement", argument=25270)
    _assert_row(nodes, 25270, "product", "delta Z^2", arguments=[25220, 25249])
    _assert_row(nodes, 25220, "analytic_unary", "delta source", name="exp", argument=25219)
    _assert_row(nodes, 25219, "negative", "delta exponent", argument=25218)
    _assert_row(nodes, 25218, "sum", "delta source exponent", arguments=[14, 4194])
    _require(_exact(_row(nodes, 1), 1), "graph node 1 must remain exact one")
    _require(_exact(_row(nodes, 48), 2), "graph node 48 must remain exact two")
    _require(_exact(_row(nodes, 13), 4), "graph node 13 must remain exact four")
    _assert_row(nodes, Z_NODE, "bound_variable", "Z", name="Z")
    return {
        "Pstar": PSTAR_NODE,
        "Z": Z_NODE,
        "delta_one_minus": COMMON_DELTA_ONE_NODE,
        "d_one_minus_Z2": COMMON_D_NODE,
        "L": COMMON_L_NODE,
        "physical_J_applied_once": True,
    }


_QY_Z0_LAYOUT: dict[str, dict[str, int]] = {
    "postrestore": {
        "target": 183451,
        "numerator": 183446,
        "two_Z_Uy": 183398,
        "Uy_term": 183386,
        "Uy": 183274,
        "negative_mean_terms": 183441,
        "mean_terms": 183436,
        "one_delta_Z_mean": 183415,
        "Z_mean": 183403,
        "mean0": 183292,
        "d_shift_mean": 183420,
        "mean1": 183293,
    },
    "patch": {
        "target": 204500,
        "numerator": 204495,
        "two_Z_Uy": 204447,
        "Uy_term": 204435,
        "Uy": 204323,
        "negative_mean_terms": 204490,
        "mean_terms": 204485,
        "one_delta_Z_mean": 204464,
        "Z_mean": 204452,
        "mean0": 204341,
        "d_shift_mean": 204469,
        "mean1": 204342,
    },
}


def _validate_physical_formula(nodes: Sequence[Mapping[str, Any]], source: Mapping[str, Any]) -> dict[str, Any]:
    _validate_common_source_geometry(nodes)
    for side in SIDES:
        physical = source["ports"][side]["physical"]
        radial = source["ports"][side]["radial"]
        pstar = source["ports"][side]["chart"]["Pstar"]
        _require(pstar == PSTAR_NODE, f"{side} Pstar source changed")
        for order in ORDERS:
            u = physical["U"][order]
            uy = physical["U_y"][order]
            mean = physical["m"][order]
            mean_y = physical["m_y"][order]
            v = EXPECTED_LOOP_PORTS[side][("V_N", order)][0 if side == "postrestore" else 1]
            vy = radial["V_y"][order]
            m = EXPECTED_MOMENT_PORTS[side][("m", order)]
            _assert_row(nodes, u, "product", f"{side}.U[{order}]", arguments=[v, pstar])
            _assert_row(nodes, uy, "product", f"{side}.U_y[{order}]", arguments=[vy, pstar])
            _assert_row(nodes, physical["m"][order], "product", f"{side}.m[{order}]", arguments=[m, pstar])
            neg = 183286 + order if side == "postrestore" else 204335 + order
            _assert_row(nodes, neg, "negative", f"{side}.negative_mean[{order}]", argument=mean)
            _assert_row(nodes, mean_y, "sum", f"{side}.mean_y[{order}]", arguments=[u, neg])
        layout = _QY_Z0_LAYOUT[side]
        target = physical["Q_y"][0]
        _require(target == layout["target"], f"{side}.Q_y Z0 target changed")
        _assert_row(nodes, target, "positive_quotient", f"{side}.Q_y[Z0]", numerator=layout["numerator"], denominator=COMMON_L_NODE)
        target_row = _row(nodes, target)
        _require(target_row.get("source_positive_certificate") == "same original L=1-delta*Z^2",
                 f"{side}.Q_y[Z0] positivity guard changed")
        _assert_row(nodes, layout["numerator"], "sum", f"{side}.Q_y numerator", arguments=[layout["two_Z_Uy"], layout["negative_mean_terms"]])
        _assert_row(nodes, layout["two_Z_Uy"], "product", f"{side}.2 Z Uy", arguments=[layout["Uy_term"], 48])
        _assert_row(nodes, layout["Uy_term"], "product", f"{side}.Z Uy", arguments=[Z_NODE, layout["Uy"]])
        _assert_row(nodes, layout["negative_mean_terms"], "negative", f"{side}.negative mean numerator", argument=layout["mean_terms"])
        _assert_row(nodes, layout["mean_terms"], "sum", f"{side}.mean numerator terms", arguments=[layout["one_delta_Z_mean"], layout["d_shift_mean"]])
        _assert_row(nodes, layout["one_delta_Z_mean"], "product", f"{side}.(1-delta) Z mean", arguments=[COMMON_DELTA_ONE_NODE, layout["Z_mean"]])
        _assert_row(nodes, layout["Z_mean"], "product", f"{side}.Z mean", arguments=[Z_NODE, layout["mean0"]])
        _assert_row(nodes, layout["d_shift_mean"], "product", f"{side}.d shift(mean)", arguments=[COMMON_D_NODE, layout["mean1"]])
        _require(physical["Q_y"][1] == EXPECTED_PHYSICAL_PORTS[side]["Q_y"][1], f"{side}.Q_y Z1 target changed")
    return {
        "source_formula": "Q_y=div(2*Z*trim(U_y)-((1-delta)*Z*trim(U-m)+(1-Z^2)*shift(U-m)))",
        "expanded_formula": "Q_y[Z0]=(2*Z*U_y[Z0]-((1-delta)*Z*mean_y[Z0]+(1-Z^2)*mean_y[Z1]))/(1-delta*Z^2)",
        "mean_y_source": "U-mean",
        "radial_mean_y_used": False,
        "Pstar_applied_once": True,
        "J_applied_once": True,
        "radius_preserved_once": True,
        "Z1_targets_open": True,
    }


class _BoundedEvaluator:
    """Exact symbolic rewrite for the bounded physical Q_y Z0 DAG."""

    _ANALYTIC = {"exp", "expm1", "log", "log1p", "sqrt", "sin", "cos"}
    _COMMON_ATOMS = {PSTAR_NODE: "PSTAR", 25220: "DELTA_EXP", COMMON_DELTA_ONE_NODE: "ONE_MINUS_DELTA", COMMON_D_NODE: "ONE_MINUS_Z2", COMMON_L_NODE: "L"}

    def __init__(self, nodes: Sequence[Mapping[str, Any]], aliases: Mapping[int, Any]) -> None:
        _require(_sp is not None, "SymPy is required for exact source normalization")
        self.nodes = nodes
        self.aliases = dict(aliases)
        self.cache: dict[int, Any] = {}
        self.active: set[int] = set()
        self.used_nodes: set[int] = set()
        self.operation_counts: Counter[str] = Counter()
        self.alias_hits: list[int] = []

    def _atom(self, label: str) -> Any:
        return _sp.Symbol(label, commutative=True)

    def _child(self, parent: int, value: Any, label: str) -> int:
        child = _handle(value, f"{label} child", len(self.nodes))
        _require(child < parent, f"{label} has a forward source reference")
        return child

    def rewrite(self, index: int) -> Any:
        index = _handle(index, "rewrite handle", len(self.nodes))
        if index in self.cache:
            return self.cache[index]
        _require(index not in self.active, f"cyclic source DAG at node {index}")
        self.used_nodes.add(index)
        if index in self.aliases:
            self.alias_hits.append(index)
            value = self.aliases[index]
            self.cache[index] = value
            return value
        if index in self._COMMON_ATOMS:
            value = self._atom(f"SOURCE:{self._COMMON_ATOMS[index]}")
            self.cache[index] = value
            self.alias_hits.append(index)
            return value
        row = _row(self.nodes, index, "rewrite node")
        operation = _operation(self.nodes, index, "rewrite node")
        self.operation_counts[operation] += 1
        self.active.add(index)
        try:
            if operation == "exact_rational":
                numerator = row.get("numerator", row.get("value"))
                denominator = row.get("denominator", 1)
                _require(isinstance(numerator, int) and isinstance(denominator, int) and denominator != 0,
                         f"exact rational {index} is malformed")
                value = _sp.Rational(numerator, denominator)
            elif operation == "bound_variable":
                name = row.get("name")
                _require(index == Z_NODE and name == "Z", f"changed or unadmitted bound variable at node {index}")
                value = _sp.Symbol("Z", commutative=True)
            elif operation in {"sum", "add"}:
                args = row.get("arguments", row.get("operands"))
                _require(_is_sequence(args) and args, f"sum {index} has no arguments")
                value = sum((self.rewrite(self._child(index, child, f"sum {index}")) for child in args), _sp.Integer(0))
            elif operation in {"product", "mul"}:
                args = row.get("arguments", row.get("operands"))
                _require(_is_sequence(args) and args, f"product {index} has no arguments")
                value = _sp.Integer(1)
                for child in args:
                    value *= self.rewrite(self._child(index, child, f"product {index}"))
            elif operation in {"negative", "neg"}:
                value = -self.rewrite(self._child(index, row.get("argument"), f"negative {index}"))
            elif operation in {"difference", "sub", "subtract"}:
                args = row.get("arguments", row.get("operands"))
                _require(_is_sequence(args) and len(args) == 2, f"difference {index} must have two arguments")
                value = self.rewrite(self._child(index, args[0], f"difference {index} lhs")) - self.rewrite(self._child(index, args[1], f"difference {index} rhs"))
            elif operation in {"quotient", "div"}:
                args = row.get("arguments", row.get("operands"))
                _require(_is_sequence(args) and len(args) == 2, f"quotient {index} must have two arguments")
                numerator = self.rewrite(self._child(index, args[0], f"quotient {index} numerator"))
                denominator = self.rewrite(self._child(index, args[1], f"quotient {index} denominator"))
                _require(_sp.simplify(denominator) != 0, f"quotient {index} has an exact zero denominator")
                value = numerator / denominator
            elif operation == "positive_quotient":
                _require(row.get("source_positive_certificate") not in (None, "", False), f"positive quotient {index} lacks its source guard")
                numerator = self.rewrite(self._child(index, row.get("numerator"), f"positive quotient {index} numerator"))
                denominator = self.rewrite(self._child(index, row.get("denominator"), f"positive quotient {index} denominator"))
                _require(_sp.simplify(denominator) != 0, f"positive quotient {index} has an exact zero denominator")
                value = numerator / denominator
            elif operation == "analytic_unary":
                name = row.get("name")
                _require(name in self._ANALYTIC, f"unsupported analytic operation {name!r} at node {index}")
                argument = self.rewrite(self._child(index, row.get("argument"), f"analytic {index}"))
                value = getattr(_sp, name)(argument)
            elif operation in {
                "function_substitution",
                "substitute_original_native_coordinate",
                "substitute_original_inverse_angle",
                "inverse_substitution",
                "definite_integral",
                "original_lazy_flat_branch",
                "compact_raw_beta",
                "exact_original_source_flat_sigma",
                "current_original_source_parameter",
                "current_original_leading_function_recipe",
                "current_original_physical_source_quantity",
                "current_original_physical_source_scale",
                "current_original_physical_radius",
                "current_original_unique_analytic_core_component",
                "current_original_axis_pressure_function_derivative",
                "original_Rm_unique_leading_implicit_control",
                "exact_positive_integer_expression",
                "mathematical_pi",
            }:
                raise PhysicalQyValueSeamError(
                    f"source opaque operation {operation!r} reached outside an authenticated common atom at node {index}"
                )
            else:
                raise PhysicalQyValueSeamError(f"unknown source operation at node {index}: {operation!r}")
            self.cache[index] = value
            return value
        finally:
            self.active.remove(index)


def _formal_atom(prefix: str, token: str) -> Any:
    _require(_sp is not None, "SymPy is required for formal VALUE atoms")
    return _sp.Symbol(f"{prefix}:{token}", commutative=True)


def _record(value: Any) -> dict[str, Any]:
    _require(_sp is not None, "SymPy is required for formal VALUE records")
    zero = _sp.simplify(value) == 0
    return {
        "srepr": _sp.srepr(value),
        "text": str(value),
        "is_zero": bool(zero) if isinstance(zero, bool) else False,
        "formal": True,
        "numeric_value": False,
    }


def _alias_forms(loop: Mapping[tuple[str, int], Mapping[str, Any]], moments: Mapping[tuple[str, int], Mapping[str, Any]], radial: Mapping[tuple[str, int], Mapping[str, Any]]) -> tuple[dict[int, Any], dict[str, Any]]:
    aliases: dict[int, Any] = {}
    roles: dict[str, Any] = {"V_N": {}, "m": {}, "V_y": {}}
    for order in ORDERS:
        guard = {
            "ordinary_Z_order": order,
            "Z_handle": Z_NODE,
            "endpoint": 1,
            "postrestore_coordinate": POST_COORDINATE,
            "patch_coordinate": PATCH_COORDINATE,
            "changed_Z": False,
            "changed_native_environment": False,
        }
        loop_row = loop[("V_N", order)]
        form = loop_row["postrestore"]["value"]
        token = canonical_digest(form)
        atom = _sp.Integer(0) if form["is_zero"] else _formal_atom("VALUE_ALIAS:V_N", f"{order}:{token}")
        for side in SIDES:
            node = loop_row[side]["node"]
            aliases[node] = atom
        roles["V_N"][str(order)] = {
            "atom": _record(atom),
            "ordinary_order": order,
            "differentiated": False,
            "postrestore_handle": loop_row["postrestore"]["node"],
            "patch_handle": loop_row["patch"]["node"],
            "value_form": deepcopy(form),
            "zero": form["is_zero"],
            "topology": "loop.V_N",
            "guard": deepcopy(guard),
        }
        moment_row = moments[("m", order)]
        token = canonical_digest({"field": "m", "order": order, "body_pair": [moment_row[side]["actual_modified_body"] for side in SIDES], "status": "proved"})
        atom = _formal_atom("VALUE_ALIAS:modified_m", f"{order}:{token}")
        for side in SIDES:
            aliases[moment_row[side]["actual_modified_body"]] = atom
        roles["m"][str(order)] = {
            "atom": _record(atom),
            "ordinary_order": order,
            "differentiated": False,
            "postrestore_handle": moment_row["postrestore"]["actual_modified_body"],
            "patch_handle": moment_row["patch"]["actual_modified_body"],
            "value_identity": True,
            "native_identity": True,
            "mean_premise": deepcopy(moment_row["mean_premise"]),
            "zero": None,
            "topology": "modified_moment.m",
            "guard": deepcopy(guard),
        }
        radial_row = radial[("V_y", order)]
        form = radial_row["postrestore"]["value"]
        token = canonical_digest(form)
        atom = _sp.Integer(0) if form["is_zero"] else _formal_atom("VALUE_ALIAS:V_y", f"{order}:{token}")
        for side in SIDES:
            aliases[radial_row[side]["target_handle"]] = atom
        roles["V_y"][str(order)] = {
            "atom": _record(atom),
            "ordinary_order": order,
            "differentiated": False,
            "postrestore_handle": radial_row["postrestore"]["target_handle"],
            "patch_handle": radial_row["patch"]["target_handle"],
            "value_form": deepcopy(form),
            "zero": form["is_zero"],
            "topology": "modified_radial.V_y",
            "guard": deepcopy(guard),
        }
    return aliases, roles


def _canonical_equal(left: Any, right: Any) -> bool:
    _require(_sp is not None, "SymPy is required for exact source comparison")
    try:
        return bool(_sp.simplify(left - right) == 0)
    except Exception:
        return False


def static_source_proof() -> dict[str, Any]:
    """Return static provenance and scope metadata without touching a graph."""

    return {
        "schema": "original_Rm_physical_qy0_value_seam_static_v1",
        "source_module": PRODUCER_MODULE,
        "source_module_sha256": PRODUCER_MODULE_SHA256,
        "source_family": deepcopy(_COMMON_SOURCE_FAMILY),
        "source_function": "physical_pressure_fields",
        "source_lines": [193, 202],
        "finite_archive_sha256": FINITE_ARCHIVE_SHA256,
        "loop_ports_sha256": LOOP_PORTS_SHA256,
        "accepted_report_sha256": ACCEPTED_REPORT_SHA256,
        "physical_formula": {
            "U": "Pstar*fields[V]",
            "U_y": "Pstar*radial[V_y]",
            "mean": "Pstar*fields[m]",
            "mean_y": "U-mean",
            "Q_y": "div(2*Z*trim(U_y)-((1-delta)*Z*trim(mean_y)+(1-Z^2)*shift(mean_y)))",
            "L": "1-delta*Z^2",
        },
        "scope": {
            "ordinary_Q_y_Z0_only": True,
            "ordinary_Q_y_Z1_proved": False,
            "full_RQy_germ_proved": False,
            "numeric_values_computed": False,
            "pressure_quadrature_executed": False,
            "old_operator_replayed": False,
            "differentiated_opaque_alias": False,
        },
    }


def derive(
    nodes: Sequence[Mapping[str, Any]],
    bounded_loop_ports_admission: Mapping[str, Any],
    modified_moment_germ: Mapping[str, Any],
    modified_radial_value_seam: Mapping[str, Any],
    actual_source_bindings: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Consume the admitted source packets and prove physical ``Q_y[Z0]``.

    ``actual_source_bindings`` is mandatory: the result is source-bound only
    when the caller supplies the postrestore and finite patch archive slices.
    All four VALUE receipts remain mandatory, and the radial receipt must pass
    both total ``V_y`` orders.
    """

    _require(_sp is not None, "SymPy is required for exact source normalization")
    expected = _digest_hex(expected_nodes_digest, "expected_nodes_digest")
    graph_digests = _validate_graph_digest(nodes, expected)
    loop_packet = _mapping(bounded_loop_ports_admission, "bounded_loop_ports_admission")
    moment_packet = _mapping(modified_moment_germ, "modified_moment_germ")
    radial_packet = _mapping(modified_radial_value_seam, "modified_radial_value_seam")
    source = _validate_source_bindings(actual_source_bindings, len(nodes))
    loop = _validate_loop_ports(loop_packet, len(nodes))
    moments = _validate_modified_moment(moment_packet, len(nodes))
    radial = _validate_modified_radial(radial_packet, len(nodes))
    formula = _validate_physical_formula(nodes, source)

    aliases, alias_roles = _alias_forms(loop, moments, radial)
    # Do not let a paired VALUE atom enter a changed Z/native environment.  The
    # bounded target grammar has no substitution node; encountering one in a
    # reachable target therefore fails closed inside _BoundedEvaluator.
    post_eval = _BoundedEvaluator(nodes, aliases)
    patch_eval = _BoundedEvaluator(nodes, aliases)
    post_expr = post_eval.rewrite(EXPECTED_PHYSICAL_PORTS["postrestore"]["Q_y"][0])
    patch_expr = patch_eval.rewrite(EXPECTED_PHYSICAL_PORTS["patch"]["Q_y"][0])
    equal = _canonical_equal(post_expr, patch_expr)
    pair = {
        "field": "Q_y",
        "ordinary_order": 0,
        "postrestore_target_handle": EXPECTED_PHYSICAL_PORTS["postrestore"]["Q_y"][0],
        "patch_target_handle": EXPECTED_PHYSICAL_PORTS["patch"]["Q_y"][0],
        "postrestore": _record(post_expr),
        "patch": _record(patch_expr),
        "equal": equal,
        "status": "proved" if equal else "failed",
        "source_formula": formula["expanded_formula"],
        "numeric_value": False,
    }
    z1 = {
        "field": "Q_y",
        "ordinary_order": 1,
        "postrestore_target_handle": EXPECTED_PHYSICAL_PORTS["postrestore"]["Q_y"][1],
        "patch_target_handle": EXPECTED_PHYSICAL_PORTS["patch"]["Q_y"][1],
        "status": "open",
        "equal": False,
        "proved": False,
        "reason": "Q_y Z1 remains outside this bounded consumer, even when its graph is source-shaped",
    }
    passed = bool(equal)
    return {
        "schema": SCHEMA,
        "passed": passed,
        "source_only": True,
        "actual_source_graph_inspected": True,
        "actual_saved_graph_run": False,
        "numeric_values_computed": False,
        "pressure_quadrature_executed": False,
        "old_mathematical_operator_replayed": False,
        "accepted_graph_replayed": False,
        "differentiated_opaque_alias": False,
        "expected_nodes_digest": expected,
        "nodes_digest": graph_digests["full_graph"],
        "graphdigest": graph_digests["full_graph"],
        "graph_node_count": len(nodes),
        "graph_prefix_digests": graph_digests,
        "canonical_value_pairs": [pair],
        "value_pairs": [pair],
        "q_y_z1": z1,
        "full_RQy_germ_Z1_proved": False,
        "physical_RQy_Z1_proved": False,
        "conditional_on": {
            "modified_radial_value_seam_passed": True,
            "required_radial_fields": list(RADIAL_FIELDS),
            "required_V_y_orders": list(ORDERS),
            "E_y_used_in_Qy0": False,
            "V_y_used_in_Qy0": True,
        },
        "physical_formula": formula,
        "physical_ports": deepcopy(source["ports"]),
        "alias_roles": alias_roles,
        "alias_policy": {
            "source_bound_endpoint_alias_only": True,
            "endpoint": 1,
            "postrestore_coordinate": POST_COORDINATE,
            "patch_coordinate": PATCH_COORDINATE,
            "ordinary_Z1_alias_is_independent": True,
            "ordinary_Z1_was_not_differentiated": True,
            "aliases_disabled_under_changed_Z_or_native_environment": True,
        },
        "walker": {
            "postrestore_used_node_count": len(post_eval.used_nodes),
            "patch_used_node_count": len(patch_eval.used_nodes),
            "postrestore_used_node_digest": canonical_digest(sorted(post_eval.used_nodes)),
            "patch_used_node_digest": canonical_digest(sorted(patch_eval.used_nodes)),
            "postrestore_operation_counts": dict(sorted(post_eval.operation_counts.items())),
            "patch_operation_counts": dict(sorted(patch_eval.operation_counts.items())),
            "postrestore_alias_hits": sorted(set(post_eval.alias_hits)),
            "patch_alias_hits": sorted(set(patch_eval.alias_hits)),
        },
        "input_data_digests": {
            "bounded_loop_ports_admission": canonical_digest(loop_packet),
            "modified_moment_germ": canonical_digest(moment_packet),
            "modified_radial_value_seam": canonical_digest(radial_packet),
            "actual_source_bindings": source["digest"],
        },
        "provenance": {
            "producer_module": PRODUCER_MODULE,
            "producer_module_sha256": PRODUCER_MODULE_SHA256,
            "source_family": deepcopy(_COMMON_SOURCE_FAMILY),
            "finite_archive_sha256": FINITE_ARCHIVE_SHA256,
            "loop_ports_sha256": LOOP_PORTS_SHA256,
            "accepted_report_sha256": ACCEPTED_REPORT_SHA256,
            "modified_mean_external_premise": {
                "label": "external_premise",
                "report_sha256": EXTERNAL_MEAN_PREMISE_SHA256,
                "certificate_digest": EXTERNAL_MEAN_CERTIFICATE_DIGEST,
            },
        },
        "preserved_ledger": {
            "all_212_forward_terms_and_errors_unchanged": True,
            "structural_channels": 6,
            "signed_channels": 5,
            "core_swirl_channel_is_error_only": True,
            "reserved_channels": 1,
            "raw_channels": 213,
            "reduced_channels": 183,
            "P1_equals_m15": False,
            "physical_Q_equals_loop_q": False,
        },
        "scope": {
            "ordinary_Q_y_Z0_value_congruence": passed,
            "ordinary_Q_y_Z1_value_congruence": False,
            "full_physical_join": False,
            "global_five_moment_closure": False,
            "heat_layer": False,
            "stress_layer": False,
            "n_recursion": False,
            "corrections": False,
            "Cartesian_NS": False,
            "full_goal_complete": False,
        },
    }


__all__ = [
    "EXPECTED_GRAPH_DIGEST",
    "EXPECTED_GRAPH_NODE_COUNT",
    "EXPECTED_LOOP_PORTS",
    "EXPECTED_MOMENT_PORTS",
    "EXPECTED_PHYSICAL_PORTS",
    "EXPECTED_RADIAL_PORTS",
    "EXPECTED_SOURCE_BINDINGS",
    "PhysicalQyValueSeamError",
    "SCHEMA",
    "canonical_digest",
    "derive",
    "static_source_proof",
]
