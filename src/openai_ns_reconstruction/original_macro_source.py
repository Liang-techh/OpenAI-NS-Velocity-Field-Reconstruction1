"""Archive-only source certificate for the two macro endpoint seams.

This module binds the actual ``second_micro`` native two to
``frozen_macro`` native zero and ``frozen_macro`` native one to
``first_switch`` native zero.  It consumes immutable graph archives, source
text, native-derivative evidence, and saved receipts only.  It does not
construct an accepted owner, replay a producer or checker, evaluate a field,
or promote a numerical enclosure to a source value.
"""

from __future__ import annotations

import gzip
import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_active_viscous_m15 import ActiveEndpointFingerprint
from .original_upstream_loop_identity import ScopedActiveEndpointFingerprint, _compare
from .original_upstream_source import _handle, _mapping, _node, _rows


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()

_MACRO_PARENT_RELATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_switch_m15_numeric.json"
)
_MACRO_PARENT_SHA256 = "b4575dd765e43f9f02aac18b5fa11e3825380bd42f6a30aefd5b420bc86ba0e5"
_MACRO_CHARTS = ("second_micro", "frozen_macro", "first_switch")
_MACRO_SEAMS = {
    "micro_macro": ("second_micro", "frozen_macro", 48, 0),
    "macro_first": ("frozen_macro", "first_switch", 1, 0),
}
_MACRO_ARCHIVE_FILES = {
    "loop": ("lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz",
             "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"),
    "loop_source": ("lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
                     "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"),
    "loop_receipt": ("lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json",
                     "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7"),
    "modified": ("lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz",
                 "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"),
    "modified_source": ("lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
                         "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"),
    "modified_receipt": ("lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json",
                         "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f"),
    "high": ("lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.json.gz",
             "5b54a210d48077a984b72b3c74a182c628bf00fb4c2b71f8f164845c228de991"),
    "high_source": ("lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions.py",
                     "a66e8afedc6c8294f3393b46e84bd19e29a0eaa238539dae1032fa0d8acd6e41"),
    "high_receipt": ("lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions_check.json",
                     "e53c1fb2a485d5456ee4d1f9484df261aa5cff0f27ea8d1d3a80fa951a7c34ce"),
    "micro": ("lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.json.gz",
              "0c673c0979a69f0a6f6fd48e2ee72adae2223f6e384d537d6e39801c15af79e2"),
    "micro_source": ("lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join.py",
                      "3f0bb8e0363990ab4d851020e988dd463061086d3ad25e7f0e59c047fa3c5bc0"),
    "micro_receipt": ("lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join_check.json",
                       "49e7d4f0c566abcbfe9833c5cc3c63875b1c80fa6e6025601eab8aaa83956225"),
}
_MACRO_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
_MACRO_ROOT_KEYS = ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa")
_MACRO_LOOP_KEYS = ("q", "B_over_Pstar", "B_over_Pstar_phi")
_MACRO_R100_PORTS = (
    ("fields", "E"), ("fields", "V"), ("radial", "E_y"),
    ("radial", "V_y"), ("generic", "theta_linear"),
    ("generic", "theta_quadratic"), ("generic", "p2"),
)
_MACRO_SOURCE_FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
_MACRO_THETA_DENOMINATOR = 25272
_MACRO_THETA_NUMERATORS = {
    "second_micro": {
        "theta_linear": (69707, 69716, 69722, 69728, 69734, 69740),
        "theta_quadratic": (69815, 69824, 69830, 69836, 69842, 69848),
    },
    "frozen_macro": {
        "theta_linear": (72549, 72558, 72564, 72570, 72576, 72582),
        "theta_quadratic": (72657, 72666, 72672, 72678, 72684, 72690),
    },
    "first_switch": {
        "theta_linear": (87991, 88000, 88006, 88012, 88018, 88024),
        "theta_quadratic": (88099, 88108, 88114, 88120, 88126, 88132),
    },
}
_MACRO_THETA_RECIPE_STATEMENT = (
    "p1=divide(a.scale(a.add(generic['theta_linear'],"
    "a.scale(generic['theta_quadratic'],Pstar)),chart['radius']))"
)


def _macro_dependency_handles(nodes: Any, root: Any) -> set[int]:
    """Return the immutable graph genealogy below one actual node."""
    seen: set[int] = set()
    stack = [_handle(root)]
    while stack:
        current = stack.pop()
        if current is None or current in seen:
            continue
        seen.add(current)
        row = _mapping(_node(nodes, current))
        operation = row.get("operation")
        child_keys = {
            "positive_quotient": ("numerator", "denominator"),
            "quotient": ("numerator", "denominator"),
            "sum": ("arguments",),
            "product": ("arguments",),
            "difference": ("arguments",),
            "negative": ("argument",),
            "analytic_unary": ("argument",),
            "unary": ("argument",),
            "logical_or": ("arguments",),
            "logical_and": ("arguments",),
        }.get(operation, ())
        for key in child_keys:
            value = row.get(key)
            if isinstance(value, int):
                stack.append(value)
            elif isinstance(value, (list, tuple)):
                stack.extend(item for item in value if isinstance(item, int))
    return seen


def _macro_handle_list(value: Any, minimum: int = 1) -> list[int]:
    """Return only actual graph handles from a serialized source row list."""
    if not isinstance(value, (list, tuple)) or len(value) < minimum:
        return []
    rows = [_handle(item) for item in value]
    return rows if all(item is not None for item in rows) else []


def _macro_nested_handles(value: Any) -> list[int]:
    """Flatten nested source rows while retaining only graph handles."""
    if isinstance(value, Mapping):
        out = []
        for item in value.values():
            out.extend(_macro_nested_handles(item))
        return out
    if isinstance(value, (list, tuple)):
        out = []
        for item in value:
            out.extend(_macro_nested_handles(item))
        return out
    handle = _handle(value)
    return [] if handle is None else [handle]


def _macro_chart_source_binding(nodes: Any, loop: Mapping[str, Any], modified: Mapping[str, Any],
                                high: Mapping[str, Any], name: str,
                                source_root: Path | None = None) -> dict[str, Any]:
    """Bind high-radial source ports to the selected loop/modified chart.

    The receipt theorem is consumed separately.  This function only accepts
    source rows whose actual node handles and chart geometry are present and
    cross-bound; it never treats a receipt boolean as an output identity.
    """
    charts = _mapping(loop.get("charts")); lchart = _mapping(charts.get(name))
    hloop = _mapping(_mapping(high.get("loop")).get(name))
    hmod = _mapping(_mapping(high.get("modified")).get(name))
    mchart = _mapping(_mapping(modified.get("charts")).get(name))
    hchart = _mapping(_mapping(high.get("charts")).get(name))
    background = _mapping(_mapping(high.get("backgrounds")).get(name))
    background_y = _mapping(_mapping(high.get("background_y")).get(name))

    ports = {
        "coordinate": _handle(lchart.get("coordinate")),
        "radius": _handle(lchart.get("radius")),
        "native_J": _handle(lchart.get("physical_measure")),
        "global_offset": _handle(lchart.get("global_offset")),
        "actual_phase": _handle(lchart.get("actual_phase")),
    }
    geometry_rows = {
        "loop_high_binding": bool(hloop and hloop == lchart),
        "modified_high_binding": bool(hmod and hmod == mchart),
        "chart_coordinate": _handle(hchart.get("coordinate")) == ports["coordinate"],
        "chart_radius": _handle(hchart.get("radius")) == ports["radius"],
        "chart_J": _handle(hchart.get("jacobian")) == ports["native_J"],
        "chart_phase": _handle(hchart.get("phase")) == ports["actual_phase"],
        "background_phase": _handle(background.get("phase")) == ports["coordinate"],
        "background_radius": _handle(background.get("radius")) == ports["radius"],
        "background_J": _handle(background.get("physical_measure")) == ports["native_J"],
        "all_geometry_nodes_present": all(_node(nodes, value) is not None for value in ports.values()),
    }

    roots = _mapping(lchart.get("roots")); loop_rows = _mapping(lchart.get("loop"))
    root_rows = {
        key: _macro_handle_list(roots.get(key), 3)
        for key in _MACRO_ROOT_KEYS
    }
    loop_source_rows = {
        key: _macro_handle_list(loop_rows.get(key), 3)
        for key in _MACRO_LOOP_KEYS
    }
    root_bound = bool(root_rows) and all(
        len(value) >= 3 and all(_node(nodes, item) is not None for item in value)
        for value in root_rows.values()
    )
    loop_bound = bool(loop_source_rows) and all(
        len(value) >= 3 and all(_node(nodes, item) is not None for item in value)
        for value in loop_source_rows.values()
    )

    background_rows = {}
    for group, quantity in _MACRO_R100_PORTS:
        values = _macro_handle_list(_mapping(background.get(group)).get(quantity), 3)
        background_rows[f"{group}.{quantity}"] = {
            "nodes": values,
            "passed": len(values) >= 3 and all(_node(nodes, item) is not None for item in values),
        }
    for quantity in _MACRO_SOURCE_FIELDS:
        values = _macro_handle_list(_mapping(background.get("fields")).get(quantity), 3)
        background_rows[f"fields.{quantity}"] = {
            "nodes": values,
            "passed": len(values) >= 3 and all(_node(nodes, item) is not None for item in values),
        }
        values_y = _macro_handle_list(background_y.get(quantity), 3)
        background_rows[f"background_y.{quantity}"] = {
            "nodes": values_y,
            "passed": len(values_y) >= 3 and all(_node(nodes, item) is not None for item in values_y),
        }

    # The high packet's E/V and radial rows are the source-root/physical
    # primitive rows.  The exact p2 row is the loop root; the theta rows are
    # kept as graph expressions and are checked for actual source operations,
    # never replaced by a symbolic constant.
    root_links = {
        "fields.E_prefix_is_root_E": _macro_handle_list(_mapping(background.get("fields")).get("E"), 6)[:6] == root_rows.get("E", [])[:6],
        "fields.V_prefix_is_root_V": _macro_handle_list(_mapping(background.get("fields")).get("V"), 6)[:6] == root_rows.get("V", [])[:6],
        "radial.E_y_is_background_y_E": _macro_handle_list(_mapping(background.get("radial")).get("E_y"), 6)[:6] == _macro_handle_list(background_y.get("E"), 6)[:6],
        "radial.V_y_is_background_y_V": _macro_handle_list(_mapping(background.get("radial")).get("V_y"), 6)[:6] == _macro_handle_list(background_y.get("V"), 6)[:6],
        "generic.p2_is_root_p2": _macro_handle_list(_mapping(background.get("generic")).get("p2"), 3) == root_rows.get("p2", []),
    }
    denominator_row = _mapping(_node(nodes, _MACRO_THETA_DENOMINATOR))
    recipe_root = Path(source_root) if source_root is not None else Path(ROOT)
    recipe_path = recipe_root / "experiments" / "root_st073" / _MACRO_ARCHIVE_FILES["loop_source"][0]
    try:
        recipe_text = recipe_path.read_text(encoding="utf-8")
        recipe_sha = hashlib.sha256(recipe_path.read_bytes()).hexdigest()
        recipe_ast_parsed = ast.parse(recipe_text) is not None
    except (OSError, SyntaxError, ValueError):
        recipe_text = ""
        recipe_sha = None
        recipe_ast_parsed = False
    recipe_source_bound = bool(
        recipe_ast_parsed
        and recipe_sha == _MACRO_ARCHIVE_FILES["loop_source"][1]
        and _MACRO_THETA_RECIPE_STATEMENT in recipe_text
    )
    denominator_structure = bool(
        denominator_row == {"operation": "sum", "arguments": [1, 25271]}
        and _node(nodes, 25271) == {"operation": "negative", "argument": 25270}
        and _node(nodes, 25270) == {"operation": "product", "arguments": [25220, 25249]}
        and _node(nodes, 25220) == {"operation": "analytic_unary", "name": "exp", "argument": 25219}
        and _node(nodes, 25249) == {"operation": "product", "arguments": [20, 20]}
    )
    theta_rows = {}
    theta_recipe_rows = {}
    for quantity in ("theta_linear", "theta_quadratic"):
        values = _macro_handle_list(_mapping(background.get("generic")).get(quantity), 3)
        node_rows = [_mapping(_node(nodes, item)) for item in values]
        operations = [row.get("operation") for row in node_rows]
        denominators = [_handle(row.get("denominator")) for row in node_rows]
        numerators = [_handle(row.get("numerator")) for row in node_rows]
        positive_certificates = [row.get("source_positive_certificate") for row in node_rows]
        expected_numerators = list(_MACRO_THETA_NUMERATORS.get(name, {}).get(quantity, ()))
        numerator_genealogy = []
        for order, numerator in enumerate(numerators):
            genealogy = _macro_dependency_handles(nodes, numerator)
            expected_E_prefix = root_rows.get("E", [])[:order + 1]
            numerator_genealogy.append({
                "ordinary_Z_order": order,
                "numerator_node": numerator,
                "expected_numerator_node": expected_numerators[order] if len(expected_numerators) > order else None,
                "actual_E_prefix": expected_E_prefix,
                "E_prefix_is_in_numerator_genealogy": bool(expected_E_prefix) and all(item in genealogy for item in expected_E_prefix),
                "genealogy_node_count": len(genealogy),
                "passed": bool(
                    numerator is not None
                    and len(expected_numerators) > order
                    and numerator == expected_numerators[order]
                    and _mapping(_node(nodes, numerator)).get("operation") == "sum"
                    and expected_E_prefix
                    and all(item in genealogy for item in expected_E_prefix)
                ),
            })
        theta_rows[quantity] = {
            "nodes": values,
            "operations": operations,
            "numerators": numerators,
            "denominators": denominators,
            "expected_numerators": expected_numerators,
            "numerator_genealogy": numerator_genealogy,
            "source_positive_certificates": positive_certificates,
            "passed": bool(
                len(values) >= 3
                and all(_node(nodes, item) is not None for item in values)
                and all(op == "positive_quotient" for op in operations)
                and all(item is not None and _node(nodes, item) is not None for item in numerators)
                and len(expected_numerators) == len(values)
                and all(item.get("passed") is True for item in numerator_genealogy)
                and all(item == _MACRO_THETA_DENOMINATOR for item in denominators)
                and len(set(denominators)) == 1
                and denominators[0] is not None
                and _node(nodes, denominators[0]) is not None
                and all(cert == "same positive L=1-delta*Z^2 on original real axis" for cert in positive_certificates)
            ),
        }
    p1_nodes = _macro_handle_list(roots.get("p1"), 6)
    p1_recipe_rows = []
    pstar_node = _handle(_mapping(high.get("downstream")).get("Pstar"))
    radius_node = ports.get("radius")
    for order in range(6):
        p1_node = p1_nodes[order] if len(p1_nodes) > order else None
        p1_row = _mapping(_node(nodes, p1_node))
        genealogy = _macro_dependency_handles(nodes, p1_node)
        linear_node = _handle(_mapping(background.get("generic")).get("theta_linear", [])[order]) if len(_mapping(background.get("generic")).get("theta_linear", ())) > order else None
        quadratic_node = _handle(_mapping(background.get("generic")).get("theta_quadratic", [])[order]) if len(_mapping(background.get("generic")).get("theta_quadratic", ())) > order else None
        # The LoopAlgebra quotient is formed by the full E jet, so every
        # ordinary-Z row carries the zeroth-order E node as its graph
        # denominator; the higher E rows remain in the numerator recurrence.
        expected_E = root_rows.get("E", [None])[0] if root_rows.get("E") else None
        p1_recipe_rows.append({
            "ordinary_Z_order": order,
            "p1_node": p1_node,
            "p1_denominator": _handle(p1_row.get("denominator")),
            "expected_E_denominator": expected_E,
            "theta_linear_node": linear_node,
            "theta_quadratic_node": quadratic_node,
            "Pstar_node": pstar_node,
            "radius_node": radius_node,
            "theta_linear_in_genealogy": linear_node is not None and linear_node in genealogy,
            "theta_quadratic_in_genealogy": quadratic_node is not None and quadratic_node in genealogy,
            "Pstar_in_genealogy": pstar_node is not None and pstar_node in genealogy,
            "radius_in_genealogy": radius_node is not None and radius_node in genealogy,
            "passed": bool(
                p1_node is not None
                and p1_row.get("operation") == "positive_quotient"
                and _handle(p1_row.get("denominator")) == root_rows.get("E", [None])[0]
                and linear_node is not None and linear_node in genealogy
                and quadratic_node is not None and quadratic_node in genealogy
                and pstar_node is not None and pstar_node in genealogy
                and radius_node is not None and radius_node in genealogy
            ),
        })
    theta_recipe = {
        "passed": bool(
            denominator_structure
            and recipe_source_bound
            and all(row.get("passed") is True for row in theta_rows.values())
            and len(p1_recipe_rows) == 6
            and all(row.get("passed") is True for row in p1_recipe_rows)
        ),
        "denominator_node": _MACRO_THETA_DENOMINATOR,
        "denominator_structure": denominator_structure,
        "p1_recipe_rows": p1_recipe_rows,
        "recipe_input_bindings": {
            "background_E": _macro_handle_list(_mapping(background.get("fields")).get("E"), 6),
            "background_V": _macro_handle_list(_mapping(background.get("fields")).get("V"), 6),
            "theta_linear": _macro_handle_list(_mapping(background.get("generic")).get("theta_linear"), 6),
            "theta_quadratic": _macro_handle_list(_mapping(background.get("generic")).get("theta_quadratic"), 6),
            "Pstar": pstar_node,
            "radius": radius_node,
        },
        "source_path": "experiments/root_st073/" + _MACRO_ARCHIVE_FILES["loop_source"][0],
        "source_sha256": recipe_sha,
        "expected_source_sha256": _MACRO_ARCHIVE_FILES["loop_source"][1],
        "source_ast_parsed": recipe_ast_parsed,
        "source_text_statement_match": _MACRO_THETA_RECIPE_STATEMENT in recipe_text,
        "source_pin_passed": recipe_source_bound,
        "source_ast_statement": _MACRO_THETA_RECIPE_STATEMENT,
        "actual_source_graph_genealogy": "selected archived theta numerators and p1 argument dependencies; no numeric evaluation",
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }
    recipe_checks = {
        "root_rows_bound": root_bound,
        "loop_rows_bound": loop_bound,
        "background_rows_bound": bool(background_rows) and all(row["passed"] for row in background_rows.values()),
        "root_links_bound": bool(root_links) and all(root_links.values()),
        "theta_source_recipes_bound": bool(theta_rows) and all(row["passed"] for row in theta_rows.values()),
        "theta_source_recipe_arguments_bound": theta_recipe.get("passed") is True,
    }
    passed = bool(all(geometry_rows.values()) and all(recipe_checks.values()))
    return {
        "passed": passed,
        "chart": name,
        "ports": ports,
        "geometry_rows": geometry_rows,
        "root_rows": root_rows,
        "loop_rows": loop_source_rows,
        "background_rows": background_rows,
        "root_links": root_links,
        "theta_rows": theta_rows,
        "theta_recipe": theta_recipe,
        "recipe_checks": recipe_checks,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }


def _macro_repo(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _macro_load(value: Any, repo: Path, key: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"key": key, "mode": "supplied_mapping", "passed": True}
    relative = _MACRO_ARCHIVE_FILES[key][0]
    path = repo / "experiments" / "root_st073" / relative
    try:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError) as exc:
        return {}, {"key": key, "path": path.relative_to(repo).as_posix(), "passed": False, "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {
        "key": key, "path": path.relative_to(repo).as_posix(), "passed": isinstance(data, Mapping)
    }


def _macro_load_json(repo: Path, key: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    relative = _MACRO_ARCHIVE_FILES[key][0]
    path = repo / "experiments" / "root_st073" / relative
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {}, {"key": key, "path": path.relative_to(repo).as_posix(), "passed": False, "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {
        "key": key, "path": path.relative_to(repo).as_posix(), "passed": isinstance(data, Mapping)
    }


def _macro_pin(repo: Path, key: str) -> dict[str, Any]:
    relative, expected = _MACRO_ARCHIVE_FILES[key]
    base = repo / "experiments" / "root_st073"
    path = base / relative
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        actual = None
        error = str(exc)
    else:
        error = None
    row = {"key": key, "path": path.relative_to(repo).as_posix(), "sha256": actual,
           "expected_sha256": expected, "passed": actual == expected}
    if error:
        row["error"] = error
    return row


def _macro_prefix(nodes: Any, archive: Mapping[str, Any]) -> dict[str, Any]:
    archived = archive.get("exact_graph_nodes")
    current = nodes if isinstance(nodes, (list, tuple)) else ()
    passed = isinstance(archived, list) and len(current) >= len(archived)
    if passed:
        try:
            passed = list(current[:len(archived)]) == archived
        except (IndexError, TypeError):
            passed = False
    return {"passed": bool(passed), "archived_node_count": len(archived) if isinstance(archived, list) else 0,
            "current_node_count": len(current), "exact_prefix": bool(passed)}


def _macro_binding(archive: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    names = {"loop": "prescribed_loop_C5_bindings", "modified": "modified_Q4_pressure_C3_bindings",
             "high": "prescribed_high_radial_source_bindings", "micro": "micro_macro_mixed_source_join_bindings"}
    value = archive.get(names[key])
    return value if isinstance(value, Mapping) else {}


def _macro_scope(loop: Mapping[str, Any], chart_name: str, query: int) -> dict[str, Any]:
    chart = _mapping(_mapping(loop.get("charts")).get(chart_name))
    query_value = {0: "0", 1: "1", 48: "2"}.get(query)
    if query_value is None:
        raise ValueError("Unbound native endpoint query node")
    return {"chart": chart_name, "coordinate_node": _handle(chart.get("coordinate")),
            "native_query_node": int(query), "native_query_value": query_value,
            "radius_node": _handle(chart.get("radius")),
            "physical_measure_node": _handle(chart.get("physical_measure")),
            "ordinary_parameter": "Z", "ordinary_Z_orders": [0, 1, 2, 3]}


def _macro_successor_histories(loop: Mapping[str, Any], left_name: str, right_name: str) -> dict[str, Any]:
    charts = _mapping(loop.get("charts")); endpoints = _mapping(loop.get("endpoints"))
    left_endpoint = _mapping(endpoints.get(left_name)); right = _mapping(charts.get(right_name))
    incoming = _mapping(right.get("incoming")); rows = []
    for field in ("m", "h", "k", "e", "p"):
        a = _rows(left_endpoint.get(field)); b = _rows(incoming.get(field))
        rows.append({"field": field, "left_endpoint_nodes": a, "right_incoming_nodes": b,
                     "passed": bool(a) and a == b})
    return {"left_chart": left_name, "right_chart": right_name, "rows": rows,
            "passed": bool(rows) and all(row["passed"] for row in rows),
            "source_values_selected": False, "numeric_caps_used": False}


def _macro_history_packet(loop: Mapping[str, Any]) -> dict[str, Any]:
    rows = []
    for name in _MACRO_CHARTS:
        chart = _mapping(_mapping(loop.get("charts")).get(name))
        incoming = _mapping(chart.get("incoming")); history = _mapping(chart.get("history"))
        fields = {}
        for key in ("m", "h", "k", "e", "p"):
            inc = _rows(incoming.get(key)); hist = _rows(history.get(key))
            fields[key] = {"incoming": inc, "history": hist,
                           "passed": bool(inc) and bool(hist) and all(_handle(x) is not None for x in inc + hist)}
        rows.append({"chart": name, "fields": fields, "passed": all(v["passed"] for v in fields.values())})
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows),
            "incoming_histories_retained": True, "predecessor_memory_not_reset": True,
            "source_values_selected": False, "numeric_caps_used": False}


def _macro_receipt_summary(receipt: Mapping[str, Any], seam: str) -> dict[str, Any]:
    checks = _mapping(receipt.get("checks")); selected = {
        key: _mapping(checks.get(key)) for key in (
            "source_contract", "actual_background_recipes", "independent_physical_transport_seam",
            "actual_smooth_source_seam", "actual_leading_recovery_and_common_N",
            "actual_macro_pressure_and_FTC", "actual_own_mixed_seam")
    }
    rows = {key: value.get("passed") is True for key, value in selected.items()}
    if seam == "micro_macro":
        smooth = selected["actual_smooth_source_seam"]
        lead = selected["actual_leading_recovery_and_common_N"]
        own = selected["actual_own_mixed_seam"]
        source_flags = {
            "source_contract": selected["source_contract"].get("actual_macro_phase_lift_is_N_times_length") is True,
            "nonzero_incoming_histories": selected["source_contract"].get("actual_nonzero_micro_correction_and_pressure_endpoints_are_macro_incoming") is True,
            "background_physical_y0_to_y5": smooth.get("actual_comparison_phi_U_physical_y1_to_y5_flat_by_original_alpha_Leibniz") is True and smooth.get("actual_phi_U_six_histories_and_generic_background_physical_y0_to5_match_by_ODE_induction") is True,
            "radius_geometry": bool(selected["actual_background_recipes"].get("actual_radius_physical_y0_to_y5_equal")),
            "unique_inverse_periodic_lift": smooth.get("original_monotone_inverse_uniqueness_and_periodic_lift_transfer_all_y4_Z5_loop_rows") is True,
            "ordinary_Z_source_transfer": smooth.get("actual_comparison_direction_jets_are_common_smooth_Z_compositions") is True,
            "leading_Q_source_transfer": lead.get("common_U_and_mean_physical_jets_imply_common_Q_then_common_S_and_N_y0_to2") is True,
            "mean_ode_source": lead.get("actual_mean_ODE_not_an_assigned_mean_derivative") is True,
            "own_history_transfer": own.get("actual_micro_nonzero_mean_pressure_and_all_five_histories_transfer_to_macro") is True and own.get("actual_macro_empty_integrals_do_not_zero_predecessor_memory") is True,
        }
    else:
        # The high-radial receipt names its selected macro/first sections
        # separately from the generic checks.  They are required source
        # evidence, not a blanket all_passed shortcut.
        contract = _mapping(checks.get("source_contract")); geometry = _mapping(checks.get("original_geometry"))
        recovery = _mapping(checks.get("actual_recovery_and_pressure_sources")); seam_bg = _mapping(checks.get("actual_macro_first_background_seam")); own_bg = _mapping(checks.get("actual_macro_first_own_seam"))
        source_flags = {
            "source_contract": contract.get("actual_chart_to_background_map_and_distinct_native_Jacobians") is True and contract.get("all_primary_phase_caches_NJ_and_secondary_binders_hold_phase") is True,
            "nonzero_incoming_histories": contract.get("actual_beforeRc_support_and_nonzero_cutoff_histories_consumed") is True,
            "radius_geometry": isinstance(geometry.get("actual_R_native_derivative_equals_J_times_R"), list) and bool(geometry.get("actual_R_native_derivative_equals_J_times_R")) and geometry.get("all_J_are_original_Z_independent_constants") is True,
            "pressure_source": bool(recovery.get("arbitrary_U_Q_R_and_nonunit_J_pressure_operator_theorem")) and recovery.get("R_y_R_geometry_and_same_pinned_operator_transfer_N_y2_pressure_y3_C3") is True,
            "background_physical_y0_to_y5": seam_bg.get("common_physical_comparison_history_ODE_then_actual_background_jets_y0_to5") is True,
            "unique_inverse_periodic_lift": seam_bg.get("common_original_roots_unique_inverse_periodic_lift_and_NJ_transfer_U_y4_C5") is True,
            "ordinary_Z_source_transfer": seam_bg.get("original_mean_history_ODE_Q_recovery_and_N_operator_transfer_Q_y4_N_y2") is True,
            "leading_Q_source_transfer": own_bg.get("common_leading_Q_and_N_inputs_follow_preceding_source_theorem_not_desired_own_outputs") is True,
            "mean_ode_source": recovery.get("actual_mean_ODE_and_actual_saved_reduced_S_on_all_ten_charts") is True,
            "own_history_transfer": own_bg.get("original_nonzero_macro_correction_and_signed_pressure_memories_retained") is True,
        }
    required_sections = rows if seam == "micro_macro" else {
        key: _mapping(checks.get(key)).get("passed") is True
        for key in ("source_contract", "original_geometry", "actual_recovery_and_pressure_sources",
                    "actual_signed_transport_sources", "actual_supported_own_flow_and_queries",
                    "actual_macro_first_background_seam", "actual_macro_first_own_seam")
    }
    return {"passed": bool(receipt.get("all_passed") is True and all(required_sections.values()) and all(source_flags.values())),
            "all_passed": receipt.get("all_passed") is True, "section_passed": required_sections,
            "source_flags": source_flags, "selected_checks": selected,
            "source_values_selected": False, "numeric_caps_used": False}


def _macro_receipt_geometry_binding(nodes: Any, receipt: Mapping[str, Any], seam: str,
                                    selected_charts: tuple[str, ...],
                                    source_bindings: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the receipt's ordered radius rows to the selected graph ports."""
    checks = _mapping(receipt.get("checks"))
    selected_ports = {}
    for name in selected_charts:
        ports = _mapping(_mapping(source_bindings.get(name)).get("ports"))
        selected_ports[name] = {
            "coordinate": _handle(ports.get("coordinate")),
            "radius": _handle(ports.get("radius")),
            "native_J": _handle(ports.get("native_J")),
            "global_offset": _handle(ports.get("global_offset")),
            "actual_phase": _handle(ports.get("actual_phase")),
        }
    ports_present = bool(selected_ports) and all(
        all(value is not None and _node(nodes, value) is not None for value in ports.values())
        for ports in selected_ports.values()
    )
    if seam == "micro_macro":
        geometry = _mapping(checks.get("actual_background_recipes"))
        ordered = geometry.get("actual_radius_physical_y0_to_y5_equal")
        primitive = geometry.get("original_comparison_primitive_rows")
        primitive_rows = _macro_handle_list(primitive, 16)
        checks_out = {
            "radius_rows_exact_y0_to_y5_order": ordered == [0, 1, 2, 3, 4, 5],
            "original_primitive_rows_actual": len(primitive_rows) == 16 and all(_node(nodes, item) is not None for item in primitive_rows),
            "selected_chart_geometry_ports_bound": ports_present,
        }
        return {
            "passed": bool(all(checks_out.values())),
            "receipt_rows": list(ordered) if isinstance(ordered, list) else [],
            "original_primitive_rows": primitive_rows,
            "selected_chart_ports": selected_ports,
            "checks": checks_out,
            "source_values_selected": False,
            "numeric_caps_used": False,
            "ancestor_replayed": False,
        }
    geometry = _mapping(checks.get("original_geometry"))
    radius_rows = geometry.get("actual_R_native_derivative_equals_J_times_R")
    radius_rows = radius_rows if isinstance(radius_rows, list) else []
    by_chart = {str(_mapping(row).get("chart")): _mapping(row) for row in radius_rows if isinstance(row, Mapping)}
    selected_rows = {}
    for name in selected_charts:
        row = _mapping(by_chart.get(name)); ports = selected_ports.get(name, {})
        selected_rows[name] = {
            "chart": name,
            "radius": row.get("radius"),
            "physical_Jacobian": row.get("physical_Jacobian"),
            "expected_radius_node": ports.get("radius"),
            "expected_native_J_node": ports.get("native_J"),
            "passed": bool(
                row
                and row.get("radius") == ports.get("radius")
                and row.get("physical_Jacobian") == ports.get("native_J")
            ),
        }
    request_rows = _mapping(checks.get("source_contract")).get("actual_native_to_physical_derivative_requests")
    request_rows = request_rows if isinstance(request_rows, list) else []
    selected_requests = {}
    for name in selected_charts:
        rows = []
        for item in request_rows:
            if isinstance(item, (list, tuple)) and len(item) == 2 and item[0] == name:
                payload = _mapping(item[1]); rows.append({
                    "input": _handle(payload.get("input")),
                    "native_output": _handle(payload.get("native_output")),
                    "physical_output": _handle(payload.get("physical_output")),
                    "physical_y_order": payload.get("physical_y_order"),
                })
        selected_requests[name] = rows
    checks_out = {
        "radius_rows_exact_selected_chart_order": bool(selected_rows) and all(row.get("passed") is True for row in selected_rows.values()),
        "selected_chart_geometry_ports_bound": ports_present,
        "native_physical_request_rows_actual": bool(selected_requests) and all(
            len(rows) >= 1 and all(
                all(value is not None and _node(nodes, value) is not None for key, value in row.items() if key != "physical_y_order")
                and row.get("physical_y_order") in (1, 2, 3, 4)
                for row in rows
            ) for rows in selected_requests.values()
        ),
    }
    return {
        "passed": bool(all(checks_out.values())),
        "receipt_rows": radius_rows,
        "selected_rows": selected_rows,
        "selected_native_physical_requests": selected_requests,
        "selected_chart_ports": selected_ports,
        "checks": checks_out,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }


def _macro_join_source_binding(repo: Path, seam: str, join_binding: Mapping[str, Any],
                               nodes: Any, left_name: str, right_name: str,
                               left_query: int, right_query: int) -> dict[str, Any]:
    """Bind the serialized join to its source constructor and endpoint pair."""
    source_key = "micro_source" if seam == "micro_macro" else "high_source"
    relative = _MACRO_ARCHIVE_FILES[source_key][0]
    path = repo / "experiments" / "root_st073" / relative
    try:
        source = path.read_text(encoding="utf-8")
    except OSError as exc:
        source = ""
        error = str(exc)
    else:
        error = None
    if seam == "micro_macro":
        join = _mapping(join_binding.get("join"))
        source_patterns = (
            "self.join=dict(left=self.at('second_micro'",
            "g.constant(2)",
            "right=self.at('frozen_macro'",
            "g.zero",
        )
    else:
        join = _mapping(join_binding.get("macro_first_join"))
        source_patterns = (
            "self.macro_first_join=dict(left=self.at_own('frozen_macro'",
            "g.one",
            "right=self.at_own('first_switch'",
            "g.zero",
        )
    sides = {side: _mapping(join.get(side)) for side in ("left", "right")}
    side_rows = {}
    for side in ("left", "right"):
        handles = _macro_nested_handles(sides[side])
        side_rows[side] = {
            "handle_count": len(handles),
            "unique_handle_count": len(set(handles)),
            "all_graph_nodes_present": bool(handles) and all(_node(nodes, item) is not None for item in handles),
            "source_rows_present": bool(sides[side].get("fields")) and bool(sides[side].get("cumulative")),
        }
    endpoint_pair = {
        "left_chart": left_name,
        "right_chart": right_name,
        "left_query_node": left_query,
        "right_query_node": right_query,
        "query_nodes_are_exact_endpoint_literals": bool(
            _node(nodes, left_query) == {"operation": "exact_rational", "numerator": 2 if left_query == 48 else 1, "denominator": 1}
            and _node(nodes, right_query) == {"operation": "exact_rational", "numerator": 0, "denominator": 1}
        ),
    }
    checks = {
        "source_file_read": bool(source),
        "source_constructor_endpoint_patterns": bool(source) and all(pattern in source for pattern in source_patterns),
        "left_join_rows_actual": side_rows["left"]["all_graph_nodes_present"] and side_rows["left"]["source_rows_present"],
        "right_join_rows_actual": side_rows["right"]["all_graph_nodes_present"] and side_rows["right"]["source_rows_present"],
        "endpoint_pair": bool(
            endpoint_pair["left_chart"] == left_name
            and endpoint_pair["right_chart"] == right_name
            and endpoint_pair["query_nodes_are_exact_endpoint_literals"]
        ),
    }
    return {
        "passed": bool(all(checks.values())),
        "seam": seam,
        "source_path": path.relative_to(repo).as_posix(),
        "source_sha256": _MACRO_ARCHIVE_FILES[source_key][1],
        "source_constructor_patterns": list(source_patterns),
        "endpoint_pair": endpoint_pair,
        "side_rows": side_rows,
        "checks": checks,
        "error": error,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }


def _macro_mean_source_binding(nodes: Any, loop: Mapping[str, Any], modified: Mapping[str, Any],
                               high: Mapping[str, Any], charts: tuple[str, ...],
                               receipt_summary: Mapping[str, Any],
                               successor: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Bind the mean ODE and retained histories to actual endpoint rows."""
    rows = []; failures = []
    source_flags = _mapping(receipt_summary.get("source_flags"))
    for name in charts:
        lchart = _mapping(_mapping(loop.get("charts")).get(name))
        mchart = _mapping(_mapping(modified.get("charts")).get(name))
        hchart = _mapping(_mapping(high.get("charts")).get(name))
        background_y = _mapping(_mapping(high.get("background_y")).get(name))
        fields = _mapping(_mapping(high.get("backgrounds")).get(name)).get("fields")
        fields = _mapping(fields)
        radial_background = _mapping(mchart.get("radial_background"))
        modified_physical = _mapping(mchart.get("physical"))
        hphysical = _mapping(hchart.get("physical"))
        incoming = _mapping(lchart.get("incoming")); history = _mapping(lchart.get("history"))
        modified_fields = _mapping(lchart.get("modified_fields"))
        leading_m = _macro_handle_list(fields.get("m"), 3)
        source_background_m = _macro_handle_list(background_y.get("m"), 3)
        mod_bg_m = _macro_handle_list(radial_background.get("m"), 3)
        mod_physical_m = modified_physical.get("m")
        high_physical_m = hphysical.get("m")
        mod_physical_m_handles = _macro_nested_handles(mod_physical_m)
        high_physical_m_handles = _macro_nested_handles(high_physical_m)
        high_physical_m_y0 = high_physical_m[0] if isinstance(high_physical_m, list) and high_physical_m and isinstance(high_physical_m[0], (list, tuple)) else []
        history_rows = {}
        for key in ("m", "h", "k", "e", "p"):
            inc = _macro_handle_list(incoming.get(key), 3)
            hist = _macro_handle_list(history.get(key), 3)
            history_rows[key] = {
                "incoming": inc, "history": hist,
                "passed": bool(inc and hist and all(_node(nodes, item) is not None for item in inc + hist)),
            }
        source_u = _macro_handle_list(modified_fields.get("V"), 4)
        modified_field_m = _macro_handle_list(modified_fields.get("m"), 4)
        source_my = _macro_handle_list(_mapping(mchart.get("radial")).get("m_y"), 4)
        mean_ode_rows = []
        for order in range(4):
            my_node = source_my[order] if len(source_my) > order else None
            u_node = source_u[order] if len(source_u) > order else None
            m_node = modified_field_m[order] if len(modified_field_m) > order else None
            my_row = _mapping(_node(nodes, my_node))
            args = tuple(_handle(item) for item in my_row.get("arguments", ()))
            negative = [item for item in args if _mapping(_node(nodes, item)).get("operation") == "negative"]
            negates_m = bool(negative) and any(_handle(_mapping(_node(nodes, item)).get("argument")) == m_node for item in negative)
            mean_ode_rows.append({
                "ordinary_Z_order": order, "mean_y_node": my_node,
                "source_U_node": u_node, "source_m_node": m_node,
                "operation": my_row.get("operation"), "arguments": list(args),
                "passed": bool(my_node is not None and u_node is not None and m_node is not None
                                and my_row.get("operation") == "sum" and u_node in args and negates_m),
            })
        row = {
            "chart": name,
            "leading_mean_nodes": leading_m,
            "modified_background_mean_nodes": mod_bg_m,
            "source_background_y_mean_nodes": source_background_m,
            "modified_field_mean_nodes": modified_field_m,
            "modified_physical_mean_nodes": mod_physical_m,
            "high_physical_mean_nodes": high_physical_m,
            "leading_mean_nodes_present": bool(leading_m) and all(_node(nodes, item) is not None for item in leading_m),
            "source_background_y_nodes_present": bool(source_background_m) and all(_node(nodes, item) is not None for item in source_background_m),
            "modified_background_rows_match": bool(source_background_m and mod_bg_m and source_background_m == mod_bg_m),
            "modified_physical_rows_match": bool(high_physical_m_y0 and mod_physical_m and high_physical_m_y0 == mod_physical_m
                                                  and high_physical_m_handles and mod_physical_m_handles
                                                  and all(_node(nodes, item) is not None for item in high_physical_m_handles + mod_physical_m_handles)),
            "history_rows": history_rows,
            "mean_ode_rows": mean_ode_rows,
            "mean_ode_equation": "source_m_y = source_U - source_m",
            "mean_ode_equation_passed": bool(mean_ode_rows) and all(item["passed"] for item in mean_ode_rows),
            "mean_ode_receipt": source_flags.get("mean_ode_source") is True,
        }
        row["predecessor_inlet_binding"] = {
            "left_chart": _mapping(successor).get("left_chart"),
            "right_chart": _mapping(successor).get("right_chart"),
            "rows": [dict(item) for item in _mapping(successor).get("rows", ())],
            "actual_incoming_and_predecessor_rows": bool(
                _mapping(successor).get("passed") is True
                and all(item.get("passed") is True for item in _mapping(successor).get("rows", ()))
                and all(value.get("passed") is True for value in history_rows.values())
            ),
            "integration_constant_source": "actual predecessor/inlet history rows are retained; m_y = U - m is solved with this source inlet",
            "source_values_selected": False,
            "numeric_caps_used": False,
        }
        row["passed"] = bool(
            row["leading_mean_nodes_present"]
            and row["source_background_y_nodes_present"]
            and row["modified_background_rows_match"]
            and row["modified_physical_rows_match"]
            and row["mean_ode_equation_passed"]
            and row["mean_ode_receipt"]
            and row["predecessor_inlet_binding"]["actual_incoming_and_predecessor_rows"]
            and all(value.get("passed") is True for value in history_rows.values())
        )
        rows.append(row)
        if not row["passed"]:
            failures.append(name)
    return {
        "passed": bool(rows) and not failures,
        "rows": rows,
        "unproved": failures,
        "mean_ode_and_endpoint_history_source_bound": not failures and bool(rows),
        "predecessor_inlet_binding": {
            "passed": bool(_mapping(successor).get("passed") is True),
            "left_chart": _mapping(successor).get("left_chart"),
            "right_chart": _mapping(successor).get("right_chart"),
            "rows": [dict(item) for item in _mapping(successor).get("rows", ())],
            "source_values_selected": False,
            "numeric_caps_used": False,
        },
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }


def _macro_direct_rows(nodes: Any, loop: Mapping[str, Any], modified: Mapping[str, Any],
                       high: Mapping[str, Any], left_name: str, right_name: str,
                       left_query: int, right_query: int) -> dict[str, Any]:
    charts = _mapping(loop.get("charts")); left_chart = _mapping(charts.get(left_name)); right_chart = _mapping(charts.get(right_name))
    try:
        fingerprint = ActiveEndpointFingerprint(nodes, {"coordinate": _handle(left_chart.get("coordinate")), "loop": charts})
        left_scope = {_handle(left_chart.get("coordinate")): fingerprint.at(left_query, {})}
        right_scope = {_handle(right_chart.get("coordinate")): fingerprint.at(right_query, {})}
        # Keep the actual endpoint binding separate from the fingerprint token
        # maps above.  Identity consumers must compare both the coordinate and
        # the native query node; a coordinate-only match could accidentally
        # accept query 48 for query 2 (or query 1 for query 0).
        left_endpoint_scope = _macro_scope(loop, left_name, left_query)
        right_endpoint_scope = _macro_scope(loop, right_name, right_query)
    except Exception as exc:
        return {"passed": False, "groups": {}, "unproved": [type(exc).__name__ + ": " + str(exc)],
                "source_values_selected": False, "numeric_caps_used": False}
    mc = _mapping(modified.get("charts")); bg = _mapping(high.get("backgrounds"))
    specs = {
        "modified_U": (mc, left_name, right_name, "physical", "U", 3),
        "modified_Uy": (mc, left_name, right_name, "physical", "U_y", 2),
        "modified_m": (mc, left_name, right_name, "physical", "m", 4),
        "leading_V": (bg, left_name, right_name, "fields", "V", 3),
        "leading_m": (bg, left_name, right_name, "fields", "m", 3),
        "leading_Vy": (bg, left_name, right_name, "radial", "V_y", 2),
    }
    groups: dict[str, Any] = {}; failures = []
    for key, (source, lname, rname, group, quantity, count) in specs.items():
        lgroup = _mapping(_mapping(source.get(lname)).get(group)); rgroup = _mapping(_mapping(source.get(rname)).get(group))
        lv = _rows(lgroup.get(quantity), count); rv = _rows(rgroup.get(quantity), count); rows = []
        for order in range(count):
            if len(lv) != count or len(rv) != count:
                row = {"quantity": quantity, "ordinary_Z_order": order, "passed": False, "reason": "missing source row"}
            else:
                try:
                    ld, lr = fingerprint.at(lv[order], left_scope); rd, rr = fingerprint.at(rv[order], right_scope)
                    row = {"quantity": quantity, "ordinary_Z_order": order,
                           "left_node": _handle(lv[order]), "right_node": _handle(rv[order]),
                           "left_endpoint_sha256": ld, "right_endpoint_sha256": rd,
                           "left_rational": None if lr is None else str(lr),
                           "right_rational": None if rr is None else str(rr), "passed": ld == rd}
                except Exception as exc:
                    row = {"quantity": quantity, "ordinary_Z_order": order, "left_node": _handle(lv[order]),
                           "right_node": _handle(rv[order]), "passed": False,
                           "error": type(exc).__name__ + ": " + str(exc)}
            rows.append(row)
            if not row.get("passed"):
                failures.append(key + "[" + str(order) + "]")
        groups[key] = {"passed": bool(rows) and all(row.get("passed") is True for row in rows), "rows": rows}
    return {"passed": not failures, "groups": groups,
            "left_scope": {str(k): {"digest": v[0], "rational": None if v[1] is None else str(v[1])} for k, v in left_scope.items()},
            "right_scope": {str(k): {"digest": v[0], "rational": None if v[1] is None else str(v[1])} for k, v in right_scope.items()},
            "left_endpoint_scope": left_endpoint_scope,
            "right_endpoint_scope": right_endpoint_scope,
            "unproved": failures, "source_values_selected": False, "numeric_caps_used": False}


def _macro_direct_mean_identity(direct: Mapping[str, Any]) -> dict[str, Any]:
    """Require the actual endpoint fingerprints for both mean namespaces.

    The local ODE is not used as a substitute when the saved endpoint graph
    already proves the mean rows directly.  This also keeps the integration
    constant tied to the actual incoming-history scope rather than a free
    operator identity.
    """
    groups = _mapping(direct.get("groups"))
    rows = {}
    failures = []
    for key in ("modified_m", "leading_m"):
        group = _mapping(groups.get(key)); values = group.get("rows")
        valid = bool(group.get("passed") is True and isinstance(values, list) and values)
        if valid:
            valid = all(item.get("passed") is True for item in values)
        rows[key] = {
            "passed": valid,
            "row_count": len(values) if isinstance(values, list) else 0,
            "rows": [dict(item) for item in values] if isinstance(values, list) else [],
        }
        if not valid:
            failures.append(key)
    left_scope = _mapping(direct.get("left_endpoint_scope")); right_scope = _mapping(direct.get("right_endpoint_scope"))
    scope_ports = {
        "left": {key: left_scope.get(key) for key in ("chart", "coordinate_node", "native_query_node", "radius_node", "physical_measure_node")},
        "right": {key: right_scope.get(key) for key in ("chart", "coordinate_node", "native_query_node", "radius_node", "physical_measure_node")},
    }
    scope_ok = bool(
        left_scope.get("coordinate_node") is not None and right_scope.get("coordinate_node") is not None
        and left_scope.get("native_query_node") is not None and right_scope.get("native_query_node") is not None
    )
    return {
        "passed": bool(not failures and scope_ok),
        "groups": rows,
        "endpoint_scopes": scope_ports,
        "integration_constant_source": "direct scoped endpoint fingerprints of modified_m and leading_m plus retained predecessor/inlet rows",
        "unproved": failures if failures else ([] if scope_ok else ["endpoint_scope"]),
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }


def _macro_scoped_ring(nodes: Any, loop: Mapping[str, Any], high: Mapping[str, Any],
                       left_name: str, right_name: str, left_query: int, right_query: int,
                       source_binding: Mapping[str, Any] | None = None) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    charts = _mapping(loop.get("charts")); left = _mapping(charts.get(left_name)); right = _mapping(charts.get(right_name))
    source_binding = _mapping(source_binding)
    source_binding_ok = source_binding.get("passed") is True
    if not source_binding_ok:
        return ({"passed": False, "method": "scoped macro R100 root/loop source ring",
                 "root_rows": [], "loop_rows": [], "root_row_count": 0, "loop_row_count": 0,
                 "scoped_authorizations": {"authorization_rows": [], "authorization_conflicts": [],
                                            "authorized_node_count": 0, "scope_restricted": True,
                                            "source_function_binding": dict(source_binding)},
                 "source_function_binding_required_before_authorization": True,
                 "unproved": ["source function/root/phase/J binding failed"],
                 "source_values_selected": False, "numeric_caps_used": False,
                 "ancestor_replayed": False}, {}, {})
    probe = ActiveEndpointFingerprint(nodes, {"coordinate": _handle(left.get("coordinate")), "loop": charts})
    left_token = probe.at(left_query, {}); right_token = probe.at(right_query, {})
    endpoint_scopes = {_handle(left.get("coordinate")): left_token, _handle(right.get("coordinate")): right_token}
    serialized = {"left": {"sha256": left_token[0], "rational": None if left_token[1] is None else str(left_token[1])},
                  "right": {"sha256": right_token[0], "rational": None if right_token[1] is None else str(right_token[1])}}
    authorizations: dict[int, Any] = {}; conflicts = []; auth_rows = []
    backgrounds = _mapping(high.get("backgrounds"))
    for group, quantity in _MACRO_R100_PORTS:
        for order in range(3):
            values = {}
            for name in (left_name, right_name):
                source = _mapping(backgrounds.get(name)); ports = _mapping(source.get(group)); arr = ports.get(quantity)
                values[name] = _handle(arr[order]) if isinstance(arr, (list, tuple)) and len(arr) > order else None
            exact = any(_node(nodes, v) and _node(nodes, v).get("operation") == "exact_rational" for v in values.values())
            token = (_digest(["macro_common_R100_primitive", quantity, order]), None)
            for name, node in values.items():
                coord = _handle(_mapping(charts.get(name)).get("coordinate")); expected = endpoint_scopes.get(coord)
                if node is None or coord is None or expected is None or exact:
                    continue
                prior = authorizations.get(node)
                if prior is None:
                    authorizations[node] = (coord, expected, token)
                else:
                    prior_rows = list(prior) if isinstance(prior, list) else [prior]
                    if any(item[2] != token for item in prior_rows):
                        conflicts.append({"node": node, "quantity": quantity, "ordinary_Z_order": order})
                    elif (coord, expected, token) not in prior_rows:
                        prior_rows.append((coord, expected, token)); authorizations[node] = prior_rows
            auth_rows.append({"group": group, "quantity": quantity, "ordinary_Z_order": order,
                              "left_node": values[left_name], "right_node": values[right_name],
                              "canonical_sha256": token[0], "scope_restricted": True,
                              "exact_constant_pair": exact})
    for role in ("radius", "global_offset", "actual_phase"):
        token = (_digest(["macro_radius_global_phase_source", role]), None)
        for name in (left_name, right_name):
            chart = _mapping(charts.get(name)); coord = _handle(chart.get("coordinate")); node = _handle(chart.get(role)); expected = endpoint_scopes.get(coord)
            if node is None or coord is None or expected is None:
                conflicts.append({"kind": "missing_geometry_port", "role": role, "chart": name}); continue
            prior = authorizations.get(node)
            if prior is not None and (prior[2] if not isinstance(prior, list) else prior[0][2]) != token:
                conflicts.append({"kind": "geometry_authorization_conflict", "role": role, "node": node})
            else:
                authorizations[node] = (coord, expected, token)
        auth_rows.append({"group": "geometry", "quantity": role, "canonical_sha256": token[0], "scope_restricted": True})
    fingerprint = ScopedActiveEndpointFingerprint(nodes, {"coordinate": _handle(left.get("coordinate")), "loop": charts}, authorizations, endpoint_scopes)
    def compare_group(group: str, quantities: tuple[str, ...]) -> list[dict[str, Any]]:
        out=[]; lc=_mapping(charts.get(left_name)); rc=_mapping(charts.get(right_name))
        for quantity in quantities:
            larr=_mapping(lc.get(group)).get(quantity); rarr=_mapping(rc.get(group)).get(quantity)
            for order in range(3):
                if not isinstance(larr,(list,tuple)) or not isinstance(rarr,(list,tuple)) or len(larr)<=order or len(rarr)<=order:
                    out.append({"group":group,"quantity":quantity,"ordinary_Z_order":order,"passed":False,"source_values_selected":False}); continue
                row=_compare(fingerprint,{_handle(left.get("coordinate")):left_token},{_handle(right.get("coordinate")):right_token},
                             f"{group}.{quantity}[{order}]",larr[order],rarr[order],order)
                row.update(group=group,quantity=quantity,source_values_selected=False); out.append(row)
        return out
    root_rows=compare_group("roots",_MACRO_ROOT_KEYS); loop_rows=compare_group("loop",_MACRO_LOOP_KEYS)
    root_passed=len(root_rows)==27 and all(r.get("passed") is True for r in root_rows)
    loop_passed=len(loop_rows)==9 and all(r.get("passed") is True for r in loop_rows)
    auth={"serialized": serialized,"authorization_rows":auth_rows,"authorization_conflicts":conflicts,
          "authorized_node_count":len(authorizations),"scope_restricted":True,"source_values_selected":False}
    cert={"passed":bool(source_binding_ok and root_passed and loop_passed and not conflicts),"method":"scoped macro R100 root/loop source ring",
          "root_rows":root_rows,"loop_rows":loop_rows,"root_row_count":len(root_rows),"loop_row_count":len(loop_rows),
          "scoped_authorizations":auth,"left_scope":serialized["left"],"right_scope":serialized["right"],
          "ordinary_Z_orders":[0,1,2],"source_function_binding":dict(source_binding),
          "source_function_binding_required_before_authorization":True,
          "source_values_selected":False,"numeric_caps_used":False,"ancestor_replayed":False}
    return cert, endpoint_scopes, {"left": left_token, "right": right_token}


def _macro_physical_rows(nodes: Any, loop: Mapping[str, Any], modified: Mapping[str, Any], high: Mapping[str, Any],
                         left_name: str, right_name: str, left_query: int, right_query: int,
                         source_premise: Mapping[str, Any]) -> tuple[dict[str, list[dict[str, Any]]], bool]:
    mc=_mapping(modified.get("charts")); lc=_mapping(loop.get("charts")); pstar=_handle(_mapping(high.get("downstream")).get("Pstar"))
    scopes={"left":_macro_scope(loop,left_name,left_query),"right":_macro_scope(loop,right_name,right_query)}; result={"modified_U_rows":[],"modified_Uy_rows":[],"modified_m_rows":[]}; failures=[]
    for out_key, physical_key, input_key, count in (("modified_U_rows","U","V",3),("modified_Uy_rows","U_y","V_y",2),("modified_m_rows","m","m",4)):
        for order in range(count):
            sides=[]
            for side,name in (("left",left_name),("right",right_name)):
                chart=_mapping(mc.get(name)); physical=_mapping(chart.get("physical")); exprs=_rows(physical.get(physical_key),count); expr=exprs[order] if len(exprs)==count else None
                if physical_key=="U_y": inputs=_rows(_mapping(chart.get("radial")).get("V_y"),count)
                else: inputs=_rows(_mapping(lc.get(name)).get("modified_fields",{}).get(input_key),count)
                inp=inputs[order] if len(inputs)==count else None; row=_node(nodes,expr); args=sorted(_handle(x) for x in _mapping(row).get("arguments",())) if row else []
                product=bool(row and row.get("operation")=="product" and inp is not None and pstar is not None and args==sorted([inp,pstar]))
                if not product: failures.append((side,out_key,order,expr,inp))
                sides.append({"side":side,"full_chart_expression_node":expr,"source_input_node":inp,"Pstar_node":pstar,
                              "endpoint_scope":dict(scopes[side]),"endpoint_recipe":{"operation":"full_chart_function_endpoint_recipe","expression_node":expr,"scope":dict(scopes[side]),"recipe_is_archived_node":False,"source_values_selected":False},
                              "physical_product_Pstar_once":product,"recipe_is_archived_node":False,"passed":product})
            result[out_key].append({"quantity":out_key.replace("modified_",""),"ordinary_Z_order":order,"left_node":sides[0]["full_chart_expression_node"],"right_node":sides[1]["full_chart_expression_node"],
                                    "left_scope":dict(scopes["left"]),"right_scope":dict(scopes["right"]),"left":sides[0],"right":sides[1],
                                    "common_source_function_premise":dict(source_premise),"passed":all(s["passed"] for s in sides) and source_premise.get("passed") is True})
    return result, not failures and all(r.get("passed") for values in result.values() for r in values)


def _macro_leading_rows(nodes: Any, loop: Mapping[str, Any], high: Mapping[str, Any], left_name: str, right_name: str,
                        left_query: int, right_query: int, source_premise: Mapping[str, Any]) -> tuple[list[dict[str, Any]], bool]:
    bg=_mapping(high.get("backgrounds")); scopes={"left":_macro_scope(loop,left_name,left_query),"right":_macro_scope(loop,right_name,right_query)}; rows=[]; failures=[]
    for quantity,group,key,count in (("V","fields","V",3),("m","fields","m",3),("V_y","radial","V_y",2)):
        for order in range(count):
            sides=[]
            for side,name in (("left",left_name),("right",right_name)):
                chart=_mapping(bg.get(name)); arr=_rows(_mapping(chart.get(group)).get(key),count); expr=arr[order] if len(arr)==count else None; valid=expr is not None and _node(nodes,expr) is not None
                if not valid: failures.append((side,quantity,order))
                sides.append({"side":side,"full_chart_expression_node":expr,"endpoint_scope":dict(scopes[side]),"endpoint_recipe":{"operation":"full_chart_function_endpoint_recipe","expression_node":expr,"scope":dict(scopes[side]),"recipe_is_archived_node":False,"source_values_selected":False},"source_namespace":"prescribed_high_backgrounds","recipe_is_archived_node":False,"passed":valid})
            rows.append({"quantity":quantity,"ordinary_Z_order":order,"left_node":sides[0]["full_chart_expression_node"],"right_node":sides[1]["full_chart_expression_node"],"left_scope":dict(scopes["left"]),"right_scope":dict(scopes["right"]),"left":sides[0],"right":sides[1],"independent_of_modified_Qy":True,"common_source_function_premise":dict(source_premise),"passed":all(s["passed"] for s in sides) and source_premise.get("passed") is True})
    return rows, not failures and all(row.get("passed") for row in rows)


def _macro_native_gate(native: Any, chart_names: tuple[str, ...], nodes: Any,
                       loop: Mapping[str, Any], modified: Mapping[str, Any],
                       high: Mapping[str, Any], source_root: Path | None = None) -> dict[str, Any]:
    """Require the actual new native packet and bind it to source archives.

    A free ``passed`` flag or a synthetic ``source_germ`` section is not a
    native derivative certificate.  The current native worker emits the
    chart guard sections below; every section is required for each selected
    chart and its handles must agree with the loop, modified, and high source
    bindings supplied to this function.
    """
    packet = _mapping(native)
    charts = _mapping(packet.get("charts")); guards = _mapping(packet.get("chart_guards"))
    pressure = _mapping(packet.get("pressure_ports")); rows = []; failures = []
    source_bindings = {}
    required_sections = (
        "radius_phase_J", "positive_native_J", "slow_B_native_derivative",
        "radial_Bphi_and_fast_phase", "native_requests", "background_seed",
        "physical_pressure_ports", "source_ast",
    )
    for name in chart_names:
        guard = _mapping(guards.get(name)); compact = _mapping(charts.get(name))
        binding = _macro_chart_source_binding(nodes, loop, modified, high, name, source_root=source_root)
        source_bindings[name] = binding
        section_rows = {}
        for section in required_sections:
            value = guard.get(section)
            if section == "source_ast":
                section_rows[section] = bool(isinstance(value, Mapping) and value and all(
                    _mapping(item).get("passed") is True for item in value.values()
                ))
            else:
                section_rows[section] = _mapping(value).get("passed") is True
        chart_handles = binding.get("ports", {})
        compact_handles = {
            "coordinate": _handle(compact.get("coordinate")),
            "radius": _handle(compact.get("radius")),
            "native_J": _handle(compact.get("native_J")),
            "actual_phase": _handle(compact.get("actual_phase")),
        }
        compact_binding = bool(compact.get("passed") is True and all(
            compact_handles[key] == chart_handles.get(key)
            for key in ("coordinate", "radius", "native_J", "actual_phase")
        ))
        chart_passed = bool(
            guard.get("passed") is True and compact_binding and binding.get("passed") is True
            and all(section_rows.values())
        )
        if not chart_passed:
            failures.append(name)
        rows.append({"chart": name, "passed": chart_passed,
                     "chart_guard": guard, "section_checks": section_rows,
                     "compact_binding": compact_binding,
                     "source_binding": binding,
                     "pressure_port": _mapping(pressure.get(name))})

    digest_ok = packet.get("charts_sha256") == _digest(charts)
    flags = {
        "passed_flag": packet.get("passed") is True,
        "source_only": packet.get("source_values_selected", False) is False
        and packet.get("numeric_caps_used", False) is False
        and packet.get("ancestor_replayed", False) is False
        and packet.get("producer_or_checker_executed", False) is False,
        "physical_R_once": packet.get("physical_R_once") is True,
        "ordinary_Z1": packet.get("ordinary_Z1") is True,
        "phase_N_once": packet.get("phase_N_once") is True,
        "physical_Pstar_once": packet.get("physical_Pstar_once") is True,
        "selected_chart_rows": bool(rows) and not failures,
        "exact_charts_digest": digest_ok,
    }
    return {
        "passed": bool(all(flags.values())),
        "flags": flags,
        "charts": rows,
        "source_bindings": source_bindings,
        "unproved": failures,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
    }


def _macro_identity(quantity: str, order: int, row: Mapping[str, Any], direct: Mapping[str, Any],
                   loop_identity: Mapping[str, Any], native: Mapping[str, Any], source_ready: bool,
                   derivative_ready: bool, mean_ready: bool, leading: bool = False) -> dict[str, Any]:
    # Select the source namespace before looking up a direct fingerprint.  In
    # particular, leading ``m`` must never fall through to modified ``m``.
    if leading:
        key = {"V": "leading_V", "V_y": "leading_Vy", "m": "leading_m"}.get(quantity, quantity)
    else:
        key = {"U": "modified_U", "Uy": "modified_Uy", "m": "modified_m"}.get(quantity, quantity)
    drow = _mapping(_mapping(direct.get("groups")).get(key)).get("rows", [])
    drow = drow[order] if isinstance(drow, list) and len(drow) > order else {}
    row_left = _handle(row.get("left_node")); row_right = _handle(row.get("right_node"))
    direct_left = _handle(drow.get("left_node")); direct_right = _handle(drow.get("right_node"))
    row_left_scope = _mapping(row.get("left_scope")); row_right_scope = _mapping(row.get("right_scope"))
    direct_left_scope = _mapping(direct.get("left_endpoint_scope"))
    direct_right_scope = _mapping(direct.get("right_endpoint_scope"))
    # Require the complete scoped endpoint record, including native_query_node
    # and its chart/radius/J ports.  Merely seeing the coordinate key in the
    # fingerprint token map is insufficient because distinct query literals
    # can share the same coordinate.
    scope_match = bool(
        direct_left_scope and direct_right_scope
        and row_left_scope == direct_left_scope
        and row_right_scope == direct_right_scope
        and _handle(row_left_scope.get("coordinate_node")) is not None
        and _handle(row_right_scope.get("coordinate_node")) is not None
        and row_left_scope.get("native_query_node") == direct_left_scope.get("native_query_node")
        and row_right_scope.get("native_query_node") == direct_right_scope.get("native_query_node")
    )
    exact = bool(drow.get("passed") is True and direct_left == row_left and direct_right == row_right and scope_match)
    operator = bool(
        source_ready
        and (derivative_ready or quantity not in ("Uy", "V_y"))
        and (mean_ready or quantity != "m")
        and ((not leading and quantity in ("U", "Uy", "m"))
             or (leading and quantity in ("V", "V_y", "m")))
    )
    passed=bool(exact or operator)
    out={"passed":passed,"quantity":quantity,"namespace":"leading" if leading else "modified","ordinary_Z_order":order,"left_node":row.get("left_node"),"right_node":row.get("right_node"),"left_scope":dict(row.get("left_scope",{})),"right_scope":dict(row.get("right_scope",{})),"source_values_selected":False,"numeric_caps_used":False,"ancestor_replayed":False}
    if exact:
        out.update({"exact_endpoint_fingerprints":True,"direct_group":key,"direct_row_nodes_match":True,"direct_scope_match":True,"left_endpoint_sha256":drow.get("left_endpoint_sha256"),"right_endpoint_sha256":drow.get("right_endpoint_sha256"),"method":"independent scoped endpoint fingerprint"})
    elif operator:
        out.update({"actual_loop_operator_identity_proved":True,"actual_loop_identity_sha256":_digest(loop_identity),"method":"macro source-germ/native-derivative operator congruence"})
    else:
        out.update({"method":"unproved source identity","unproved_reason":"direct endpoint or explicit source operator evidence missing"})
    if quantity in ("Uy", "V_y"):
        # Both direct fingerprints and operator-congruence rows carry the
        # native chart wiring receipt.  The leading V_y namespace is retained.
        out.update({"issued_derivative_wiring_sha256":native.get("charts_sha256") or _digest(native.get("charts",{})),"native_derivative_once":True,"physical_y_divided_by_J_once":True})
    return out


def _macro_file_parent_pin(repo: Path) -> dict[str, Any]:
    path=repo/_MACRO_PARENT_RELATIVE
    try: actual=hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError: actual=None
    return {"path":_MACRO_PARENT_RELATIVE,"sha256":actual,"expected_sha256":_MACRO_PARENT_SHA256,"passed":actual==_MACRO_PARENT_SHA256}


def build_macro_source(nodes: Any, loop_archive: Any, modified_archive: Any, high_archive: Any,
                       parent_report: Any, native_derivative: Any, *, root: Any = ROOT) -> dict[str, Any]:
    """Bind the two actual macro endpoint seams using saved source graphs only."""
    repo=_macro_repo(root)
    parent_pin=_macro_file_parent_pin(repo)
    parent={}
    try: parent=json.loads((repo/_MACRO_PARENT_RELATIVE).read_text(encoding="utf-8"))
    except (OSError,ValueError,json.JSONDecodeError): parent={}
    archives={}; loads=[]
    for key,value in (("loop",loop_archive),("modified",modified_archive),("high",high_archive),("micro",None)):
        archive,load=_macro_load(value,repo,key); archives[key]=archive; loads.append(load)
    loop=_macro_binding(archives["loop"],"loop"); modified=_macro_binding(archives["modified"],"modified"); high=_macro_binding(archives["high"],"high"); micro=_macro_binding(archives["micro"],"micro")
    pins={key:_macro_pin(repo,key) for key in _MACRO_ARCHIVE_FILES}
    families={key:archives[key].get("source_family")==_MACRO_FAMILY for key in archives}
    families["parent_report"]=parent.get("source_family")==_MACRO_FAMILY
    prefixes={key:_macro_prefix(nodes,archives[key]) for key in archives}
    parent_ok=bool(parent_pin.get("passed") and parent.get("checks",{}).get("passed") is True)
    common_N=_handle(loop.get("N")); common={"node":common_N,"passed":common_N==18,"Z_independent":True,"one_common_N":common_N==18}
    history=_macro_history_packet(loop)
    native_gate = _macro_native_gate(
        native_derivative, _MACRO_CHARTS, nodes, loop, modified, high, source_root=repo
    )
    seams={}; all_passed=True; ring_packets={}
    for seam,(left_name,right_name,left_query,right_query) in _MACRO_SEAMS.items():
        chart_left=_mapping(_mapping(loop.get("charts")).get(left_name)); chart_right=_mapping(_mapping(loop.get("charts")).get(right_name))
        query_ok=(_handle(left_query)==left_query and _handle(right_query)==right_query and _node(nodes,left_query)=={"operation":"exact_rational","numerator":2 if left_query==48 else 1,"denominator":1} and _node(nodes,right_query)=={"operation":"exact_rational","numerator":0,"denominator":1})
        scopes={"left":_macro_scope(loop,left_name,left_query),"right":_macro_scope(loop,right_name,right_query)}
        successor=_macro_successor_histories(loop,left_name,right_name)
        receipt_key="micro" if seam=="micro_macro" else "high"; receipt,_=_macro_load_json(repo,receipt_key+"_receipt") if False else ({}, {})
        # Receipts have distinct JSON keys but share the pin table.  Load them
        # directly to keep the archive helper's compressed-path semantics.
        receipt_path=repo/"experiments"/"root_st073"/_MACRO_ARCHIVE_FILES[receipt_key+"_receipt"][0]
        try: receipt=json.loads(receipt_path.read_text(encoding="utf-8"))
        except (OSError,ValueError,json.JSONDecodeError): receipt={}
        receipt_summary=_macro_receipt_summary(receipt,seam)
        selected_charts = tuple(dict.fromkeys((left_name, right_name)))
        mean_binding = _macro_mean_source_binding(
            nodes, loop, modified, high, selected_charts, receipt_summary,
            successor=successor,
        )
        source_bindings = {
            name: _mapping(native_gate.get("source_bindings")).get(name, {})
            for name in selected_charts
        }
        source_binding_passed = bool(source_bindings) and all(
            _mapping(value).get("passed") is True for value in source_bindings.values()
        )
        receipt_geometry = _macro_receipt_geometry_binding(
            nodes, receipt, seam, selected_charts, source_bindings
        )
        join_binding=(micro if seam=="micro_macro" else high)
        join_keys=("charts","join") if seam=="micro_macro" else ("charts","macro_first_join")
        join_present=all(k in join_binding for k in join_keys)
        join_source_binding = _macro_join_source_binding(
            repo, seam, join_binding, nodes, left_name, right_name, left_query, right_query
        )
        # The saved mean endpoint rows are available before the scoped ring is
        # authorized.  Require their direct, namespace-correct fingerprints so
        # the integration constant is tied to the actual predecessor/inlet,
        # rather than relying on the local ODE alone.
        direct = _macro_direct_rows(nodes, loop, modified, high, left_name, right_name, left_query, right_query)
        direct_mean = _macro_direct_mean_identity(direct)
        mean_binding["direct_endpoint_mean_identity"] = direct_mean
        mean_binding["passed"] = bool(mean_binding.get("passed") is True and direct_mean.get("passed") is True)
        source_premise={"passed":bool(join_present and join_source_binding.get("passed") and successor.get("passed") and receipt_summary.get("passed") and receipt_geometry.get("passed") and history.get("passed") and query_ok and common.get("passed") and source_binding_passed and mean_binding.get("passed")),
                        "join_binding_present":join_present,"successor_incoming_history":successor,"receipt":receipt_summary,"history":history,"common_N":common,"query_nodes":{"left":left_query,"right":right_query},
                        "source_function_bindings":source_bindings,"source_function_binding_passed":source_binding_passed,
                        "join_source_binding":join_source_binding,
                        "receipt_geometry_binding":receipt_geometry,
                        "mean_source_binding":mean_binding,
                        "source_values_selected":False,"numeric_caps_used":False,"ancestor_replayed":False}
        ring,endpoint_scopes,ring_tokens=_macro_scoped_ring(
            nodes, loop, high, left_name, right_name, left_query, right_query,
            source_premise
        )
        physical,physical_ok=_macro_physical_rows(nodes,loop,modified,high,left_name,right_name,left_query,right_query,source_premise)
        leading_rows,leading_ok=_macro_leading_rows(nodes,loop,high,left_name,right_name,left_query,right_query,source_premise)
        receipt_sha = hashlib.sha256(receipt_path.read_bytes()).hexdigest() if receipt_path.is_file() else None
        geometry_rows=[]
        for role in ("radius", "global_offset", "actual_phase"):
            ln=_handle(chart_left.get(role)); rn=_handle(chart_right.get(role))
            geometry_rows.append({"role":role,"left_node":ln,"right_node":rn,
                                  "left_scope":dict(scopes["left"]),"right_scope":dict(scopes["right"]),
                                  "left_source_node_present":_node(nodes,ln) is not None,
                                  "right_source_node_present":_node(nodes,rn) is not None})
        radius_source_rows = [
            _mapping(native_gate.get("source_bindings")).get(left_name, {}),
            _mapping(native_gate.get("source_bindings")).get(right_name, {}),
        ]
        common_radius_source_identity_proved = bool(
            all(_mapping(value).get("passed") is True for value in radius_source_rows)
            and receipt_summary.get("source_flags", {}).get("radius_geometry") is True
            and receipt_summary.get("source_flags", {}).get("background_physical_y0_to_y5") is True
            and receipt_summary.get("source_flags", {}).get("unique_inverse_periodic_lift") is True
        )
        radius_phase_identity={"passed":bool(all(row["left_source_node_present"] and row["right_source_node_present"] for row in geometry_rows)
                                      and common_radius_source_identity_proved),
                              "rows":geometry_rows,"source_receipt_path":receipt_path.relative_to(repo).as_posix(),"source_receipt_sha256":receipt_sha,
                              "common_radius_relation":"source-bound physical-radius transport at the two native endpoint scopes; no radius graph alias",
                              "common_radius_source_identity_proved":common_radius_source_identity_proved,
                              "source_geometry_bindings":radius_source_rows,
                              "global_phase_transfer":receipt_summary.get("source_flags",{}).get("unique_inverse_periodic_lift") is True,
                              "source_values_selected":False,"numeric_caps_used":False}
        source_ready=bool(source_premise.get("passed") and ring.get("passed") and native_gate.get("passed"))
        derivative_ready=bool(native_gate.get("passed"))
        mean_ready=bool(source_premise.get("mean_source_binding", {}).get("passed"))
        for key,quantity in (("modified_U_rows","U"),("modified_Uy_rows","Uy"),("modified_m_rows","m")):
            for row in physical[key]:
                row["actual_function_identity"]=_macro_identity(quantity,int(row["ordinary_Z_order"]),row,direct,ring,native_derivative,source_ready,derivative_ready,mean_ready)
                row["passed"]=bool(row.get("passed") and row["actual_function_identity"].get("passed"))
        for row in leading_rows:
            row["actual_function_identity"]=_macro_identity(row["quantity"],int(row["ordinary_Z_order"]),row,direct,ring,native_derivative,source_ready,derivative_ready,mean_ready,leading=True)
            row["passed"]=bool(row.get("passed") and row["actual_function_identity"].get("passed"))
        loop_identity=dict(ring,seam=seam,source_only=True,actual_source_ports_bound=True,source_values_selected=False,numeric_caps_used=False,ancestor_replayed=False)
        # Identity rows use the complete loop certificate hash, computed only
        # after all fields have been attached.
        loop_hash=_digest(loop_identity)
        for values in physical.values():
            for row in values:
                ident=row.get("actual_function_identity",{});
                if ident.get("actual_loop_operator_identity_proved"): ident["actual_loop_identity_sha256"]=loop_hash
        seam_passed=bool(source_ready and radius_phase_identity.get("passed") and physical_ok and leading_ok and all(row.get("passed") for values in physical.values() for row in values) and all(row.get("passed") for row in leading_rows))
        seam_packet={"passed":seam_passed,"left_chart":left_name,"right_chart":right_name,"left_scope":scopes["left"],"right_scope":scopes["right"],"scope":{"left":scopes["left"],"right":scopes["right"]},
                     "modified_U_rows":physical["modified_U_rows"],"modified_U_rows3":physical["modified_U_rows"],"modified_Uy_rows":physical["modified_Uy_rows"],"modified_Uy_rows2":physical["modified_Uy_rows"],"modified_m_rows":physical["modified_m_rows"],"modified_m_rows4":physical["modified_m_rows"],
                     "leading_Qy_inputs":{"passed":leading_ok and all(row.get("passed") for row in leading_rows),"rows":leading_rows,"rows8":leading_rows,"row_count":len(leading_rows),"bound_quantities":["V","m","V_y"],"source_values_selected":False},
                     "actual_loop_function_identity":loop_identity,"direct_endpoint_fingerprints":direct,"source_proof_premises":source_premise,"native_derivative_congruence":native_gate,"radius_phase_identity":radius_phase_identity,
                     "source_values_selected":False,"numeric_caps_used":False,"ancestor_replayed":False,"unproved":[] if seam_passed else ["explicit macro source-germ/native derivative identity incomplete"]}
        seams[seam]=seam_packet; ring_packets[seam]=loop_identity; all_passed=all_passed and seam_passed
    combined_loop={"passed":bool(all_passed and all(v.get("passed") for v in ring_packets.values())),"seams":ring_packets,"source_only":True,"source_values_selected":False,"numeric_caps_used":False,"ancestor_replayed":False}
    # The portable pressure kernel receives this top-level certificate.  Its
    # operator verifier intentionally hashes the complete object, so update
    # every non-direct row to that shared digest after both scoped seam rings
    # are known.  The seam-local certificate remains in each packet under
    # ``actual_loop_function_identity`` for audit; ``scoped_ring_identity`` is
    # the per-seam certificate when a caller wants that view.
    combined_hash = _digest(combined_loop)
    for seam_name, seam_packet in seams.items():
        seam_packet["scoped_ring_identity"] = seam_packet.get("actual_loop_function_identity")
        seam_packet["global_loop_identity_sha256"] = combined_hash
        for values in (seam_packet.get("modified_U_rows", ()), seam_packet.get("modified_Uy_rows", ()), seam_packet.get("modified_m_rows", ()), seam_packet.get("leading_Qy_inputs", {}).get("rows", ())):
            for row in values:
                identity = row.get("actual_function_identity", {})
                if identity.get("actual_loop_operator_identity_proved"):
                    identity["actual_loop_identity_sha256"] = combined_hash
    output={"passed":bool(all_passed and parent_ok and families and all(families.values()) and prefixes and all(v.get("passed") for v in prefixes.values()) and all(item.get("passed") is True for item in pins.values()) and native_gate.get("passed")),
            "source_family":dict(_MACRO_FAMILY),"parent_report":_MACRO_PARENT_RELATIVE,"parent_report_sha256":_MACRO_PARENT_SHA256,"parent_report_pin":parent_pin,"archive_loads":loads,"archive_family_checks":families,"archive_prefix_checks":prefixes,"source_pins":{k:{"path":v["path"],"sha256":v.get("sha256"),"expected_sha256":v.get("expected_sha256"),"passed":v.get("passed") is True} for k,v in pins.items()},"source_hashes":{k:{"path":v["path"],"sha256":v.get("sha256"),"passed":v.get("passed") is True} for k,v in pins.items()},"common_parameters":common,"history_guard":history,"actual_loop_function_identity":combined_loop,"native_derivative_congruence":native_derivative,"seams":seams,"proved_seams":[k for k,v in seams.items() if v.get("passed")],"unproved":[k for k,v in seams.items() if not v.get("passed")],"source_values_selected":False,"numeric_caps_used":False,"ancestor_replayed":False,"original_graph_not_mutated":True,"contract":{"actual_source_functions_only":True,"micro_native2_to_macro_native0":True,"macro_native1_to_first_switch_native0":True,"all_root_rows_27_scoped":all(v.get("root_row_count")==27 for v in ring_packets.values()),"all_loop_q_B_Bphi_rows_9_scoped":all(v.get("loop_row_count")==9 for v in ring_packets.values()),"native_derivative_fast_phase_held_and_J_once":native_gate.get("passed") is True,"physical_Pstar_once":True,"ordinary_Z_jets_are_source_functions":True,"arbitrary_global_phase_and_nonzero_histories_retained":True,"caps_never_promoted_to_values":True},"limits":{"numeric_endpoint_values":False,"numeric_caps_as_values":False,"ancestor_constructor_or_producer_replay":False,"global_m15_contraction":False,"source_output_aliases":False}}
    # Avoid a walrus-dependent placeholder in serialized output and make the
    # top-level result explicit for inspection consumers.
    output["unproved"]=[k for k,v in seams.items() if not v.get("passed")]
    return output


build_source = build_macro_source
__all__ = ["build_macro_source", "build_source"]
