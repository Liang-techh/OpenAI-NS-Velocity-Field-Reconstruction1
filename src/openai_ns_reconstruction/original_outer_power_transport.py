"""Source-only O3 power transport certificate.

The O3 power chart is quiet because its *source* roots satisfy ``V=0`` and
the original loop takes the lazy flat branch when ``Delta >= eta``.  The
finite-N graph still contains the incoming correction history, so this file
does not look for literal zero rows in the modified packet.  It proves the
zero velocity consequence through the graph products and records the exact
mass ODE which remains in the retained history.

Only serialized archives, pinned source text, and graph metadata are read.
No owner, producer, checker, numerical source, interval cap, or ancestor
calculation is executed.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT


_OUTER_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)
def _py_path(archive: str) -> str:
    return archive[:-8] + ".py" if archive.endswith(".json.gz") else archive + ".py"


def _check_path(source: str, suffix: str) -> str:
    return source[:-3] + suffix if source.endswith(".py") else source + suffix


_OUTER_SOURCE = _py_path(_OUTER_ARCHIVE)
_OUTER_CHECK_SOURCE = _check_path(_OUTER_SOURCE, "_check.py")
_OUTER_RECEIPT = _check_path(_OUTER_SOURCE, "_check.json")
_LOOP_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
)
_LOOP_SOURCE = _py_path(_LOOP_ARCHIVE)
_LOOP_CHECK_SOURCE = _check_path(_LOOP_SOURCE, "_check.py")
_LOOP_RECEIPT = _check_path(_LOOP_SOURCE, "_check.json")
_LEADING_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity.json.gz"
)
_LEADING_SOURCE = _py_path(_LEADING_ARCHIVE)
_LEADING_CHECK_SOURCE = _check_path(_LEADING_SOURCE, "_check.py")
_LEADING_RECEIPT = _check_path(_LEADING_SOURCE, "_check.json")
_RC_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.json.gz"
)
_RC_SOURCE = _py_path(_RC_ARCHIVE)
_RC_CHECK_SOURCE = _check_path(_RC_SOURCE, "_check.py")
_RC_RECEIPT = _check_path(_RC_SOURCE, "_check.json")
_O3_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.json.gz"
)
_O3_SOURCE = _py_path(_O3_ARCHIVE)
_O3_CHECK_SOURCE = _check_path(_O3_SOURCE, "_check.py")
_O3_RECEIPT = _check_path(_O3_SOURCE, "_check.json")
_POWER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp.py"
)


_EXPECTED_HASHES = {
    _OUTER_ARCHIVE: "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
    _OUTER_SOURCE: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    _OUTER_CHECK_SOURCE: "855c1b31c42c1db22ce503e2a31e5aa20d13dbf7d5d1feb2869c048c64560a00",
    _OUTER_RECEIPT: "ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba",
    _LOOP_ARCHIVE: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    _LOOP_SOURCE: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    _LOOP_CHECK_SOURCE: "066e6e09ccbb4dc28d7b1881197a426beb6bbc7f92d393e7f3e646f27287914b",
    _LOOP_RECEIPT: "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7",
    _LEADING_ARCHIVE: "e2c5b0d363df8c24d5c0aeaeb2d919474d7a7db2db3a5cb8ddacc623d6249f11",
    _LEADING_SOURCE: "94bf63d7ea54b7e3a29d8c97efefa7e242614f9fb79f5d3d88f8fc6402fced47",
    _LEADING_CHECK_SOURCE: "cc2da7bff8a9030d438e4afeafb331178a5f3fc3a51451771260bf81191c3b36",
    _LEADING_RECEIPT: "988b3438c94c6f222dded2dc69d631aa6f99673110db9534a1ac2b80099c53b9",
    _RC_ARCHIVE: "93364a44e86b79dfbe0a942f9d6e6157a0eb0b2e149dae58833cfc50a41d72d0",
    _RC_SOURCE: "9d73659a1c76163dddee6653af53f56ebffd7045986d871dc559ba1e053620b7",
    _RC_CHECK_SOURCE: "9fc3a4921ab1a7c633f3dd60acd175b65803871ab0cc7963e1b964d90942fda9",
    _RC_RECEIPT: "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158",
    _O3_ARCHIVE: "8ad3463e8c1dd88da2dbd38a102c8aeb9892d854a9a245f795ea5da7cda6459d",
    _O3_SOURCE: "83e5a8bdae75d20a261d24dd22d30cc7cb5e603e8fa2f0133958ca880b5e76db",
    _O3_CHECK_SOURCE: "66625aba667bbca404a5cb2a2ebbc2cc7fb1f0fa128b55e2af8e3dccd6c1023d",
    _O3_RECEIPT: "30bd8cd0d502a4e2e9a8a45b9ad37fa55054ea3dd3e1dd690c5fc3c7700a3270",
    _POWER_SOURCE: "f1630fb589f9e5dc74d2f4c5ac3bff4f3e396cc91be141c7dfe6fce5bd95ec94",
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    node = getattr(value, "node", None)
    return node if isinstance(node, int) else None


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None or index < 0:
        return None
    try:
        item = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return item if isinstance(item, Mapping) else None


def _operands(row: Mapping[str, Any]) -> tuple[int, ...]:
    values: list[Any] = []
    for key in ("arguments", "argument", "numerator", "denominator"):
        value = row.get(key)
        if isinstance(value, (list, tuple)):
            values.extend(value)
        elif value is not None:
            values.append(value)
    return tuple(index for index in (_handle(value) for value in values) if index is not None)


def _reachable(nodes: Any, root: Any, target: Any, *, limit: int = 512) -> bool:
    start, wanted = _handle(root), _handle(target)
    if start is None or wanted is None:
        return False
    stack = [start]
    seen: set[int] = set()
    while stack and len(seen) < limit:
        index = stack.pop()
        if index == wanted:
            return True
        if index in seen:
            continue
        seen.add(index)
        stack.extend(_operands(_row(nodes, index) or {}))
    return False


def _exact_zero(nodes: Any, value: Any) -> bool:
    row = _row(nodes, value)
    return bool(row and row.get("operation") == "exact_rational" and row.get("numerator") == 0 and row.get("denominator") == 1)


def _direct_product(nodes: Any, value: Any, left: Any, right: Any) -> bool:
    row = _row(nodes, value)
    args = tuple(_handle(item) for item in row.get("arguments", ())) if row else ()
    return bool(row and row.get("operation") == "product" and len(args) == 2 and set(args) == {_handle(left), _handle(right)})


def _direct_sum_with_negative(nodes: Any, value: Any, positive: Any, negative: Any) -> bool:
    row = _row(nodes, value)
    if not row or row.get("operation") != "sum" or len(row.get("arguments", ())) != 2:
        return False
    args = tuple(_handle(item) for item in row["arguments"])
    neg = _row(nodes, args[1]) if len(args) == 2 else None
    return bool(args[0] == _handle(positive) and neg and neg.get("operation") == "negative" and _handle(neg.get("argument")) == _handle(negative))


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _load_json(repo: Path, relative: str) -> dict[str, Any]:
    try:
        value = json.loads((repo / relative).read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, json.JSONDecodeError):
        return {}


def _load_gzip(repo: Path, relative: str) -> dict[str, Any]:
    try:
        with gzip.open(repo / relative, "rt", encoding="utf-8") as stream:
            value = json.load(stream)
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, EOFError, json.JSONDecodeError):
        return {}


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {"path": relative, "exists": path.is_file(), "sha256": actual,
            "expected_sha256": expected, "passed": actual == expected}


def _find_packet(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    if key in value and isinstance(value[key], Mapping):
        return value
    for child_key in ("base", "archive", "packet", "outer", "leading", "true_outer", "outer_packet"):
        child = value.get(child_key)
        if isinstance(child, Mapping):
            found = _find_packet(child, key)
            if found:
                return found
    return {}


def _source_text(repo: Path, relative: str) -> str:
    try:
        return (repo / relative).read_text(encoding="utf-8")
    except OSError:
        return ""


def _ast_defs(text: str, wanted: tuple[str, ...]) -> bool:
    if not text:
        return False
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return False
    names = {node.name for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
    return all(name in names for name in wanted)


def _source_guard(repo: Path) -> dict[str, Any]:
    outer = _source_text(repo, _OUTER_SOURCE)
    loop = _source_text(repo, _LOOP_SOURCE)
    leading = _source_text(repo, _LEADING_SOURCE)
    power = _source_text(repo, _POWER_SOURCE)
    o3 = _source_text(repo, _O3_SOURCE)
    rows = {
        "outer_ast_and_pressure_geometry": _ast_defs(outer, ("geometry", "pressure_geometry")) and all(
            text in outer for text in (
                "O3_power", "physical_Jacobian_applied_once=True",
                "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,R))",
                "original_normalized_history_ODE", "original_power",
            )
        ),
        "power_ast_and_zero_velocity_profile": _ast_defs(power, ("power_flow",)) and all(
            text in power for text in (
                "alpha=g.add(g.constant('1/2'),mu)",
                "hist=dict(m=alg.scale(own['m'],decay(1))",
                "profiles=[dict(E=alg.scale(E,g.mul(g.constant((-1)**j),*[alpha]*j)),V=zero)",
                "complete_histories_y0_y1_y2_y3_y4_C3",
            )
        ),
        "leading_power_port": _ast_defs(leading, ("own_first_y", "build_source")) and all(
            text in leading for text in (
                "D=dict(m=V", "flow=quiet.power_flow", "packet('O3_power'",
                "profiles_y0_y1_y2_y3_y4_C3",
            )
        ),
        "loop_lazy_flat_source": _ast_defs(loop, ("lazy", "sigma")) and all(
            text in loop for text in (
                "g.node('original_lazy_flat_branch'",
                "active_body_not_evaluated_on_flat_branch=True",
                "same_original_flat_extension_through_order5=True",
            )
        ),
        "mu_admission_source": _ast_defs(o3, ("source_mu_admission",)) and all(
            text in o3 for text in (
                "logmu=c.ln(c.mpf('.001'))-4*(c.exp(40)+11)",
                "margin=logmu+c.ln(2)-eta_log",
                "ep(margin)[0]<0",
            )
        ),
    }
    return {"passed": all(rows.values()), "rows": rows, "no_source_execution": True}


def _receipt_guard(receipt: Mapping[str, Any], gates: tuple[str, ...], checks: tuple[str, ...] = ()) -> dict[str, Any]:
    check_map = _mapping(receipt.get("checks"))
    gate_rows = {name: receipt.get(name) is True for name in gates}
    check_rows = {name: _mapping(check_map.get(name)).get("passed") is True for name in checks}
    # Some accepted packets (the O3/Rc owner and the C3 identity owner) put
    # their source gates at the receipt top level rather than under checks.
    direct_rows = {name: receipt.get(name) is True for name in checks if name not in check_map}
    passed = receipt.get("all_passed") is True and all(gate_rows.values()) and all(check_rows.values()) and all(direct_rows.values())
    return {"passed": passed, "all_passed": receipt.get("all_passed") is True,
            "gates": gate_rows, "checks": {**check_rows, **direct_rows}}


def _family(*packets: Mapping[str, Any]) -> tuple[bool, Any]:
    values = [packet.get("source_family") for packet in packets if packet.get("source_family") is not None]
    return bool(values) and all(value == values[0] for value in values), values[0] if values else None


def _flat_rows(nodes: Any, rows: Any, delta: Any, eta: Any, predicate: Any, *, native_derivative: bool = False) -> dict[str, Any]:
    rows = list(rows) if isinstance(rows, (list, tuple)) else []
    details: list[dict[str, Any]] = []
    for order, value in enumerate(rows):
        row = _row(nodes, value) or {}
        expected = {
            "operation": row.get("operation") == "original_lazy_flat_branch",
            "Delta": _handle(row.get("Delta")) == _handle(delta),
            "eta": _handle(row.get("eta")) == _handle(eta),
            "flat_value": _handle(row.get("flat_value")) == 0,
            "flat_predicate": _handle(row.get("flat_predicate")) == _handle(predicate),
            "active_body_not_evaluated": row.get("active_body_not_evaluated_on_flat_branch") is True,
            "ordinary_Z_order": row.get("ordinary_slow_Z_derivative_order") == order,
            "same_flat_extension": row.get("same_original_flat_extension_through_order5") is True,
        }
        if native_derivative:
            expected["same_native_derivative"] = row.get("same_original_lazy_native_derivative") is True
        details.append({"node": _handle(value), "checks": expected, "passed": all(expected.values())})
    return {"passed": len(rows) == 6 and all(item["passed"] for item in details), "rows": details}


def _flat_dependencies(nodes: Any, rows: Any, delta: Any, eta: Any, predicate: Any) -> dict[str, Any]:
    """Check products/sums whose descendants contain the lazy flat branch."""
    rows = list(rows) if isinstance(rows, (list, tuple)) else []
    details: list[dict[str, Any]] = []
    for order, value in enumerate(rows):
        stack = [_handle(value)]
        seen: set[int] = set()
        candidates: list[int] = []
        while stack and len(seen) < 512:
            index = stack.pop()
            if index is None or index in seen:
                continue
            seen.add(index)
            row = _row(nodes, index) or {}
            if row.get("operation") == "original_lazy_flat_branch":
                candidates.append(index)
            stack.extend(_operands(row))
        matched = []
        for index in candidates:
            row = _row(nodes, index) or {}
            matched.append({
                "node": index,
                "Delta": _handle(row.get("Delta")) == _handle(delta),
                "eta": _handle(row.get("eta")) == _handle(eta),
                "flat_value": _handle(row.get("flat_value")) == 0,
                "flat_predicate": _handle(row.get("flat_predicate")) == _handle(predicate),
                "active_body_not_evaluated": row.get("active_body_not_evaluated_on_flat_branch") is True,
                "ordinary_Z_order": row.get("ordinary_slow_Z_derivative_order") == order,
                "same_flat_extension": row.get("same_original_flat_extension_through_order5") is True,
                "same_native_derivative": row.get("same_original_lazy_native_derivative") is True,
            })
        passed = bool(matched) and any(all(value for key, value in item.items() if key != "node") for item in matched)
        details.append({"node": _handle(value), "lazy_descendants": matched, "passed": passed})
    return {"passed": len(rows) == 6 and all(item["passed"] for item in details), "rows": details}


def _extract_outer(value: Any, fallback: Mapping[str, Any]) -> Mapping[str, Any]:
    packet = _find_packet(value, "true_outer_C6_C5_bindings")
    return packet if packet else fallback


def build_outer_power_transport(
    nodes: Any,
    joined: Mapping[str, Any] | None,
    outer: Mapping[str, Any] | None,
    leading: Mapping[str, Any] | None,
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Build the source certificate for the whole O3 power chart.

    ``nodes`` is the caller's joined graph.  The other arguments may be raw
    archive packets or runner wrappers; pinned archives supply immutable
    fallback bindings when a wrapper omits a packet.
    """
    repo = _repo_root(root)
    loaded_outer = _load_gzip(repo, _OUTER_ARCHIVE)
    loaded_loop = _load_gzip(repo, _LOOP_ARCHIVE)
    loaded_leading = _load_gzip(repo, _LEADING_ARCHIVE)
    loaded_rc = _load_gzip(repo, _RC_ARCHIVE)
    loaded_o3 = _load_gzip(repo, _O3_ARCHIVE)
    outer_packet = _find_packet(outer, "true_outer_C6_C5_bindings") or loaded_outer
    loop_packet = _find_packet(joined, "prescribed_loop_C5_bindings") or loaded_loop
    leading_packet = _find_packet(leading, "actual_original_leading_C3_source_functions") or loaded_leading
    rc_packet = loaded_rc
    o3_packet = loaded_o3
    outer_binding = _mapping(outer_packet.get("true_outer_C6_C5_bindings"))
    loop_binding = _mapping(loop_packet.get("prescribed_loop_C5_bindings"))
    leading_binding = _mapping(leading_packet.get("actual_original_leading_C3_source_functions"))
    rc_binding = _mapping(rc_packet.get("original_outer_source_binding"))
    chart = _mapping(_mapping(outer_binding.get("charts")).get("O3_power"))
    roots = _mapping(chart.get("roots"))
    loop = _mapping(chart.get("loop"))
    finite = _mapping(chart.get("finite_N"))
    radial = _mapping(chart.get("radial"))
    physical = _mapping(chart.get("physical"))
    boundary = _mapping(chart.get("boundary"))
    source_family_ok, source_family = _family(outer_packet, loop_packet, leading_packet, rc_packet, o3_packet)

    pins = {name: _pin(repo, path, expected) for name, (path, expected) in {
        "outer_archive": (_OUTER_ARCHIVE, _EXPECTED_HASHES[_OUTER_ARCHIVE]),
        "outer_source": (_OUTER_SOURCE, _EXPECTED_HASHES[_OUTER_SOURCE]),
        "outer_checker": (_OUTER_CHECK_SOURCE, _EXPECTED_HASHES[_OUTER_CHECK_SOURCE]),
        "outer_receipt": (_OUTER_RECEIPT, _EXPECTED_HASHES[_OUTER_RECEIPT]),
        "loop_archive": (_LOOP_ARCHIVE, _EXPECTED_HASHES[_LOOP_ARCHIVE]),
        "loop_source": (_LOOP_SOURCE, _EXPECTED_HASHES[_LOOP_SOURCE]),
        "loop_checker": (_LOOP_CHECK_SOURCE, _EXPECTED_HASHES[_LOOP_CHECK_SOURCE]),
        "loop_receipt": (_LOOP_RECEIPT, _EXPECTED_HASHES[_LOOP_RECEIPT]),
        "leading_archive": (_LEADING_ARCHIVE, _EXPECTED_HASHES[_LEADING_ARCHIVE]),
        "leading_source": (_LEADING_SOURCE, _EXPECTED_HASHES[_LEADING_SOURCE]),
        "leading_checker": (_LEADING_CHECK_SOURCE, _EXPECTED_HASHES[_LEADING_CHECK_SOURCE]),
        "leading_receipt": (_LEADING_RECEIPT, _EXPECTED_HASHES[_LEADING_RECEIPT]),
        "rc_archive": (_RC_ARCHIVE, _EXPECTED_HASHES[_RC_ARCHIVE]),
        "rc_source": (_RC_SOURCE, _EXPECTED_HASHES[_RC_SOURCE]),
        "rc_checker": (_RC_CHECK_SOURCE, _EXPECTED_HASHES[_RC_CHECK_SOURCE]),
        "rc_receipt": (_RC_RECEIPT, _EXPECTED_HASHES[_RC_RECEIPT]),
        "o3_archive": (_O3_ARCHIVE, _EXPECTED_HASHES[_O3_ARCHIVE]),
        "o3_source": (_O3_SOURCE, _EXPECTED_HASHES[_O3_SOURCE]),
        "o3_checker": (_O3_CHECK_SOURCE, _EXPECTED_HASHES[_O3_CHECK_SOURCE]),
        "o3_receipt": (_O3_RECEIPT, _EXPECTED_HASHES[_O3_RECEIPT]),
        "power_source": (_POWER_SOURCE, _EXPECTED_HASHES[_POWER_SOURCE]),
    }.items()}

    outer_receipt = _load_json(repo, _OUTER_RECEIPT)
    loop_receipt = _load_json(repo, _LOOP_RECEIPT)
    leading_receipt = _load_json(repo, _LEADING_RECEIPT)
    rc_receipt = _load_json(repo, _RC_RECEIPT)
    o3_receipt = _load_json(repo, _O3_RECEIPT)
    receipt_guards = {
        "outer": _receipt_guard(outer_receipt,
            ("current_original_true_six_outer_prescribed_C6_generic_C5_functions_installed",
             "current_original_all_17_chart_axial6_generic_source_installed",
             "current_original_true_six_outer_charts_finite_N_C5_functions_installed",
             "current_original_all_17_chart_finite_N_C5_functions_installed",
             "current_original_true_all17_modified_Q4_total_radial_C5_functions_installed",
             "current_original_true_rminus_to_Rc_signed_pressure_C3_functions_installed"),
            ("actual_source", "geometry", "original_seed_theorem", "actual_profiles", "actual_quiet_power", "exponential_transport", "exponential_pressure")),
        "loop": _receipt_guard(loop_receipt,
            ("current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed",
             "current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed",
             "current_original_true_prescribed_rminus_to_Rm_correction_history_C5_functions_installed"),
            ("actual_source", "lazy_flat")),
        "leading": _receipt_guard(leading_receipt, ("current_original_outer_leading_C3_source_ODE_and_seam_functions_installed",)),
        "rc": _receipt_guard(rc_receipt, (
            "current_original_whole_Z_full_phase_all_N_outer_C1_functions_normalized_Rc_and_relative_targets_installed",
            "all_four_Z_live_P0_Rm_Rh_Rc_positive_mu_amplitude_incoming_and_pressure_memory_checked",
            "original_frontends_pure_parent_recipes_physical_measures_and_quiet_incoming_retained",
            "actual_uniform_N_scaled_Rc_targets_installed",
            "actual_all_regions_all_N_adapter_installed")),
        "o3": _receipt_guard(o3_receipt, (
            "current_original_whole_Z_same_N_actual_O3_source_and_correction_through_Rc_enclosed",
            "original_whole_Z_O3_history_functions_positive_mu_critical_quotients_and_signed_C0_Z_drivers_checked",
            "genuine_Rd_correction_only_incoming_true_memory_and_Rc_exit_checked",
            "positive_mu_source_Delta_active_transition_and_exact_quiet_power_retained",
            "original_kernel_cache_reuses_only_exact_same_context_t_mu_partition",
            "full_closed_real_axial_and_O3_Rc_cover_checked",
            "conservative_full_phase_covers_not_phase_average_cancellation")),
    }

    source_guard = _source_guard(repo)
    family_guard = {
        "passed": source_family_ok,
        "source_family": source_family,
        "N": outer_binding.get("N"),
        "eta": outer_binding.get("eta"),
        "dstar": outer_binding.get("dstar"),
        "loop_N": loop_binding.get("N"),
    }
    graph_prefix = False
    outer_nodes = outer_packet.get("exact_graph_nodes", ())
    try:
        graph_prefix = isinstance(nodes, (list, tuple)) and isinstance(outer_nodes, (list, tuple)) and len(nodes) >= len(outer_nodes) and list(nodes[:len(outer_nodes)]) == list(outer_nodes)
    except (TypeError, IndexError):
        graph_prefix = False

    eta = outer_binding.get("eta")
    delta_rows = roots.get("Delta", ())
    delta = delta_rows[0] if isinstance(delta_rows, (list, tuple)) and delta_rows else None
    predicate = loop.get("flat_predicate")
    flat_checks = {name: _flat_rows(nodes, loop.get(name), delta, eta, predicate) for name in ("A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi")}
    radial_flat = _flat_dependencies(nodes, radial.get("slow_B_y"), delta, eta, predicate)
    root_checks = {
        "shape": all(isinstance(roots.get(name), (list, tuple)) and len(roots[name]) == 6 for name in ("E", "V", "a", "b", "Delta", "kappa", "t0")),
        "leading_V_exact_zero": all(_exact_zero(nodes, value) for value in roots.get("V", ())),
        "root_b_t0_exact_zero": all(_exact_zero(nodes, value) for name in ("b", "t0") for value in roots.get(name, ())),
        "positive_a_rows": all((_row(nodes, value) or {}).get("operation") == "positive_quotient" for value in roots.get("a", ())),
        "delta_source_identity": _mapping(_mapping(o3_packet.get("original_source_bindings")).get("original_source")).get("exact_collected_Delta_power") == "2*mu",
        "delta_predicate": (_row(nodes, predicate) or {}).get("operation") == "exact_real_comparison" and (_row(nodes, predicate) or {}).get("operator") == ">=",
    }
    mu_proof = _mapping(outer_binding.get("mu_proof"))
    o3_source = _mapping(_mapping(o3_packet.get("original_source_bindings")).get("original_source"))
    mu_checks = {
        "positive_mu_and_2mu_ge_eta": mu_proof.get("positive_mu_and_power_2mu_ge_eta_and_Rp_reservation") is True,
        "bounds_not_field_values": mu_proof.get("bounds_not_field_values") is True and _mapping(mu_proof.get("eta_proof")).get("parameter_ranges_not_source_field_values") is True,
        "exact_mu_recipe": mu_proof.get("exact_original_mu_recipe") == "mu=1/1000*exp(-4*(exp(40)+11))",
        "exact_o3_delta_source": o3_source.get("exact_collected_Delta_power") == "2*mu",
        "quiet_source_gate": o3_receipt.get("positive_mu_source_Delta_active_transition_and_exact_quiet_power_retained") is True,
    }

    leading_chart = _mapping(_mapping(leading_binding.get("charts")).get("O3_power"))
    leading_v = leading_chart.get("V", ())
    leading_checks = {
        "chart_present": bool(leading_chart),
        "V_rows_exact_zero": isinstance(leading_v, (list, tuple)) and len(leading_v) == 4 and all(_handle(value) == 0 for value in leading_v),
        "E_rows_present": isinstance(leading_chart.get("E"), (list, tuple)) and len(leading_chart["E"]) == 4,
        "power_flow_source_bound": _mapping(_mapping(outer_binding.get("bindings")).get("original_power")).get("function") == "power_flow",
        "leading_ODE_source_bound": _mapping(_mapping(outer_binding.get("bindings")).get("original_normalized_history_ODE")).get("function") == "own_first_y",
        "leading_source_graph_prefix_recorded": leading_packet.get("original_outer_graph_prefix_length") == 24752,
    }

    pstar_candidates: set[int] = set()
    u_rows, uy_rows = physical.get("U", ()), physical.get("U_y", ())
    for value in list(u_rows) + list(uy_rows):
        row = _row(nodes, value)
        if row and row.get("operation") == "product":
            pstar_candidates.update(_handle(item) for item in row.get("arguments", ()) if _handle(item) is not None)
    pstar = None
    if len(u_rows) == 6 and len(uy_rows) == 6:
        factors = []
        for value in list(u_rows) + list(uy_rows):
            row = _row(nodes, value) or {}
            factors.append({_handle(item) for item in row.get("arguments", ()) if _handle(item) is not None})
        common = set.intersection(*factors) if factors else set()
        if len(common) == 1:
            pstar = next(iter(common))
    finite_v = finite.get("V_N", ())
    b_rows = loop.get("B_over_Pstar", ())
    vy_rows = radial.get("V_y", ())
    bphi_rows = loop.get("B_over_Pstar_phi", ())
    velocity_checks = {
        "finite_V_N_from_flat_B": len(finite_v) == 6 and len(b_rows) == 6 and all(_direct_product(nodes, finite_v[i], b_rows[i], 99411) for i in range(6)),
        "physical_U_one_Pstar": pstar is not None and len(u_rows) == 6 and all(_direct_product(nodes, u_rows[i], finite_v[i], pstar) for i in range(6)),
        "physical_Uy_one_Pstar": pstar is not None and len(uy_rows) == 6 and len(vy_rows) == 6 and all(_direct_product(nodes, uy_rows[i], vy_rows[i], pstar) for i in range(6)),
        "radial_Vy_slow_and_fast_flat": len(vy_rows) == 6 and len(radial.get("slow_B_y", ())) == 6 and len(bphi_rows) == 6 and all(_reachable(nodes, vy_rows[i], radial["slow_B_y"][i]) and _reachable(nodes, vy_rows[i], bphi_rows[i]) for i in range(6)),
        "radial_chain_rule": radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True and radial.get("total_fast_phase_terms_included") is True,
    }
    mass_rows = physical.get("m", ())
    mass_y_rows = physical.get("m_y", ())
    mass_checks = {
        "modified_mass_ODE_U_minus_m": len(mass_rows) == 6 and len(mass_y_rows) == 6 and all(_direct_sum_with_negative(nodes, mass_y_rows[i], u_rows[i], mass_rows[i]) for i in range(6)),
        "leading_mass_ODE_source": leading_checks["leading_ODE_source_bound"] and source_guard["rows"]["leading_power_port"],
        "radius_Ry_equals_R_source": "physical_Jacobian_applied_once=True" in _source_text(repo, _OUTER_SOURCE) and _mapping(rc_binding).get("original_global_phase") == "frac(N*(logRm+original_outer_offset-logRa-hb*s_c/2))",
        "retained_incoming_m": isinstance(chart.get("incoming", {}).get("m"), (list, tuple)) and len(chart["incoming"]["m"]) == 6 and any(not _exact_zero(nodes, value) for value in chart["incoming"]["m"]),
        "memory_not_reset": _mapping(outer_binding.get("definition")).get("O3_power_quiet_local_density_does_not_reset_correction_memory") is True and o3_receipt.get("genuine_Rd_correction_only_incoming_true_memory_and_Rc_exit_checked") is True,
    }
    qy_rows = physical.get("Q_y", ())
    rqy_rows = boundary.get("R_Q_y", ())
    pressure_checks = {
        "Qy_rows": isinstance(qy_rows, (list, tuple)) and len(qy_rows) >= 4,
        "R_Qy_rows": isinstance(rqy_rows, (list, tuple)) and len(rqy_rows) == 4,
        "R_Qy_radius_once": len(rqy_rows) == 4 and len(qy_rows) >= 4 and all(_direct_product(nodes, rqy_rows[i], qy_rows[i], chart.get("radius")) for i in range(4)),
        "physical_measure_once": _handle(chart.get("physical_measure")) == 1,
    }
    proof_passed = bool(
        source_guard["passed"] and all(row["passed"] for row in pins.values())
        and all(guard["passed"] for guard in receipt_guards.values())
        and source_family_ok and graph_prefix
        and all(root_checks.values()) and all(mu_checks.values())
        and all(leading_checks.values()) and all(flat["passed"] for flat in flat_checks.values())
        and radial_flat["passed"] and all(velocity_checks.values()) and all(mass_checks.values()) and all(pressure_checks.values())
    )
    power_certificate = {
        "passed": proof_passed,
        "method": "actual_O3_power_whole_chart_flat_branch_and_mass_transport_source_guard",
        "source_only": True,
        "whole_chart_U_zero_by_flat_semantics": all(flat["passed"] for flat in flat_checks.values()) and velocity_checks["finite_V_N_from_flat_B"] and velocity_checks["physical_U_one_Pstar"],
        "whole_chart_Uy_zero_by_flat_semantics": radial_flat["passed"] and velocity_checks["radial_Vy_slow_and_fast_flat"] and velocity_checks["physical_Uy_one_Pstar"],
        "literal_modified_graph_zero_not_required": True,
        "mean_not_zeroed": True,
        "mean_ODE_m_y_equals_U_minus_m": mass_checks["modified_mass_ODE_U_minus_m"],
        "leading_ODE_m_y_equals_minus_m": mass_checks["leading_mass_ODE_source"],
        "delta_M_y_equals_minus_delta_M": mass_checks["modified_mass_ODE_U_minus_m"] and mass_checks["leading_mass_ODE_source"],
        "R_y_equals_R": mass_checks["radius_Ry_equals_R_source"],
        "D_definition": "D=R*(Q_y_modified-Q_y_leading), constant in physical y and ordinary Z1",
        "constant_D": mass_checks["modified_mass_ODE_U_minus_m"] and mass_checks["leading_mass_ODE_source"] and mass_checks["radius_Ry_equals_R_source"] and pressure_checks["R_Qy_radius_once"],
        "lower_upper_D_equal": mass_checks["modified_mass_ODE_U_minus_m"] and mass_checks["leading_mass_ODE_source"] and mass_checks["radius_Ry_equals_R_source"] and pressure_checks["R_Qy_radius_once"],
        "D_constant_in_physical_y": mass_checks["modified_mass_ODE_U_minus_m"] and mass_checks["leading_mass_ODE_source"] and mass_checks["radius_Ry_equals_R_source"],
        "physical_Q_y_rows": list(qy_rows) if isinstance(qy_rows, (list, tuple)) else [],
        "physical_R_Q_y_rows": list(rqy_rows) if isinstance(rqy_rows, (list, tuple)) else [],
        "Pstar_node": pstar,
        "radius_node": _handle(chart.get("radius")),
        "N_node": _handle(outer_binding.get("N")),
        "eta_node": _handle(eta),
        "mu_node": _handle(outer_binding.get("mu")),
        "Delta_nodes": list(delta_rows) if isinstance(delta_rows, (list, tuple)) else [],
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "producer_or_checker_executed": False,
        "source_pins": pins,
        "receipt_guards": receipt_guards,
        "checks": {
            "source_guard": source_guard,
            "source_family": family_guard,
            "graph_prefix": graph_prefix,
            "root_rows": root_checks,
            "mu_proof": mu_checks,
            "leading_ports": leading_checks,
            "flat_loop": flat_checks,
            "radial_flat": radial_flat,
            "velocity_ports": velocity_checks,
            "mass_transport": mass_checks,
            "pressure_ports": pressure_checks,
        },
        "unproved": [] if proof_passed else [
            name for name, ok in {
                "source_guard": source_guard["passed"], "source_pins": all(row["passed"] for row in pins.values()),
                "accepted_receipts": all(guard["passed"] for guard in receipt_guards.values()), "source_family": source_family_ok,
                "graph_prefix": graph_prefix, "root_rows": all(root_checks.values()), "mu_proof": all(mu_checks.values()),
                "leading_ports": all(leading_checks.values()), "flat_loop": all(flat["passed"] for flat in flat_checks.values()) and radial_flat["passed"],
                "velocity_ports": all(velocity_checks.values()), "mass_transport": all(mass_checks.values()),
                "pressure_ports": all(pressure_checks.values()),
            }.items() if not ok
        ],
    }
    # Keep this alias explicit: downstream numeric code consumes the nested
    # certificate, while source audits often look for the top-level key.
    return {
        "passed": power_certificate["passed"],
        "method": "actual_O3_power_whole_chart_transport_source_guard",
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "producer_or_checker_executed": False,
        "source_pins": pins,
        "sourcepins": pins,
        "power_transport_certificate": power_certificate,
        "certificate": power_certificate,
        "checks": power_certificate["checks"],
        "unproved": power_certificate["unproved"],
        "limits": {
            "numeric_endpoint_values": False,
            "numeric_caps_as_values": False,
            "whole_chart_velocity_literal_zero": False,
            "global_m15_contraction": False,
            "full_Navier_Stokes_completion": False,
        },
    }


__all__ = ["build_outer_power_transport"]
