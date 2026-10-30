"""One-entry source-bound uniform first-word consumer.

The runner reads the frozen selected-word receipt only after checking its raw
SHA-256.  It checks the three source modules with AST guards, injects the
saved common ``selected_Z_interval`` into the unique ``q = 1`` row, and calls
the new certificate leaf exactly once.  No old producer, archive, graph, or
ancestor helper is imported or executed here.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import json
import math
import sys
import traceback
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "n1_uniform_first_word_v1"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

REPORT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"
NUMERIC = "src/openai_ns_reconstruction/original_n1_selected_word_numeric.py"
MEASURE = "src/openai_ns_reconstruction/original_n1_cutoff_measure.py"
TAIL = "src/openai_ns_reconstruction/original_n1_selected_word_tail.py"
LEAF = "src/openai_ns_reconstruction/n1_uniform_first_word_certificate.py"
COUPLED_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_coupled_first_step.py"
RUNNER = "scripts/run_n1_uniform_first_word_once.py"
REQUIRED = frozenset((REPORT, NUMERIC, MEASURE, TAIL, LEAF, COUPLED_SOURCE, RUNNER))
FIXED = {
    REPORT: "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e",
    COUPLED_SOURCE: "97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7",
}

PRECISION_BITS = 1800
FIELDS = ("F1", "Uz1")
ORDERS = (0, 1)
ALL_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(value: Any, message: str) -> None:
    if not value:
        raise ValueError(message)


def mapping(value: Any, name: str) -> Mapping[str, Any]:
    require(isinstance(value, Mapping), f"{name} must be a mapping")
    return value


def exclusive(path: Path, value: Any) -> str:
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return sha(raw)


def preflight() -> tuple[dict[str, str], str, str]:
    """Validate the pending capsule without parsing an input before its hash."""

    require(not ENTRY.exists() and not RESULT.exists(), "already entered")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") == COMPONENT + "_pending"
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True
        and brief.get("one_entry_only") is True,
        "pending guard",
    )
    checkpoint = mapping(brief.get("source_checkpoint"), "brief.source_checkpoint")
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "checkpoint path")

    source_raw = SOURCE.read_bytes()
    require(sha(source_raw) == checkpoint.get("sha256"), "checkpoint hash")
    source = json.loads(source_raw)
    require(source.get("component") == COMPONENT and source.get("one_entry_only") is True, "scope")

    hashes: dict[str, str] = {}
    for item in source.get("files", ()):  # raw bytes first; no JSON input is parsed here
        item_map = mapping(item, "source file entry")
        path = item_map.get("path")
        require(path in REQUIRED and path not in hashes, "unexpected input")
        raw = (ROOT / path).read_bytes()
        actual = sha(raw)
        require(actual == item_map.get("sha256"), "changed " + str(path))
        if path in FIXED:
            require(actual == FIXED[path], "fixed input SHA changed " + str(path))
        hashes[path] = actual
    require(set(hashes) == REQUIRED, "missing input")
    return hashes, sha(source_raw), sha(brief_raw)


def _function(tree: ast.AST, name: str) -> ast.AST:
    matches = [node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]
    require(len(matches) == 1, f"expected one function {name}")
    return matches[0]


def _statement_dump(source: str) -> str:
    parsed = ast.parse(source, mode="exec")
    require(len(parsed.body) == 1, "AST fragment must contain one statement")
    return ast.dump(parsed.body[0], include_attributes=False)


def _expression_dump(source: str) -> str:
    return ast.dump(ast.parse(source, mode="eval").body, include_attributes=False)


def _has_statement(node: ast.AST, source: str) -> bool:
    wanted = _statement_dump(source)
    return any(ast.dump(child, include_attributes=False) == wanted for child in ast.walk(node))


def _has_expression(node: ast.AST, source: str) -> bool:
    wanted = _expression_dump(source)
    return any(ast.dump(child, include_attributes=False) == wanted for child in ast.walk(node))


def _has_constant(node: ast.AST, value: Any) -> bool:
    return any(isinstance(child, ast.Constant) and child.value == value for child in ast.walk(node))


def _has_constant_containing(node: ast.AST, fragment: str) -> bool:
    return any(
        isinstance(child, ast.Constant) and isinstance(child.value, str) and fragment in child.value
        for child in ast.walk(node)
    )


def _has_call_keyword(node: ast.AST, key: str, value_source: str) -> bool:
    wanted = _expression_dump(value_source)
    for child in ast.walk(node):
        if not isinstance(child, ast.Call):
            continue
        for keyword in child.keywords:
            if keyword.arg == key and ast.dump(keyword.value, include_attributes=False) == wanted:
                return True
    return False


def _source_ast_guards(raw: Mapping[str, bytes]) -> dict[str, Any]:
    numeric_text = raw[NUMERIC].decode("utf-8")
    measure_text = raw[MEASURE].decode("utf-8")
    tail_text = raw[TAIL].decode("utf-8")
    coupled_text = raw[COUPLED_SOURCE].decode("utf-8")
    leaf_text = raw[LEAF].decode("utf-8")
    numeric_tree = ast.parse(numeric_text, filename=NUMERIC)
    measure_tree = ast.parse(measure_text, filename=MEASURE)
    tail_tree = ast.parse(tail_text, filename=TAIL)
    coupled_tree = ast.parse(coupled_text, filename=COUPLED_SOURCE)
    leaf_tree = ast.parse(leaf_text, filename=LEAF)

    transport = _function(numeric_tree, "transport")
    local_radii = _function(numeric_tree, "local_radii")
    numeric_build = _function(numeric_tree, "build")
    measure_build = _function(measure_tree, "build")
    tail_rows = _function(tail_tree, "_raw_tail_rows")
    tail_hydrate = _function(tail_tree, "hydrate_selected_word_tail")
    coupled_tail = _function(coupled_tree, "preconditioned_tail")
    coupled_owner_tail = _function(coupled_tree, "tail")
    coupled_evaluate = _function(coupled_tree, "evaluate")

    checks = {
        "leaf_assemble_function": any(
            isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == "assemble"
            for node in ast.walk(leaf_tree)
        ),
        "transport_y_is_box_of_y_upper": _has_statement(transport, "y = box(c, y_upper)"),
        "transport_is_Z_independent": not any(
            isinstance(node, ast.Name) and node.id == "Z" for node in ast.walk(transport)
        ),
        "build_qsc_is_sc_times_q": _has_statement(numeric_build, "qsc = sc * q.numerator / q.denominator"),
        "build_y_upper_is_qsc_over_8": _has_statement(numeric_build, "y_upper = endpoint(c, qsc / 8)"),
        "build_prefix_is_W0_plus_W1": _has_statement(
            numeric_build,
            "prefix = {f: linear_combination(c, (c.mpf(1), w0[f]), (c.mpf(1), w1[f])) for f in FIELDS}",
        ),
        "local_dx_is_qsc_over_4": _has_statement(local_radii, "dx = qsc / 4"),
        "local_G_is_positive_linear_qsc": _has_statement(
            local_radii,
            "G = {i: qsc * positive_bound(c, Ai[str(i)]) / 2 for i in (4, 5, 6)}",
        ),
        "local_W0_radius_structure": _has_statement(
            local_radii,
            "r0 = [dx * G[5], dx * G[6], dx * G[6], G[4], G[5], G[6]]",
        ),
        "local_W1_radius_structure": _has_statement(
            local_radii,
            "r1 = [dx**2 * H, dx**2 * H, dx**2 * H, dx * H, dx * H, dx * H]",
        ),
        "local_radius_boxes_are_positive": _has_expression(local_radii, "box(c, v)"),
        "cutoff_eta_formula": _has_constant(measure_build, "1-sigma(4*q-2)"),
        "cutoff_eta_support": _has_call_keyword(measure_build, "eta_support_q", '["0", "3/4"]'),
        "cutoff_q_interval": _has_call_keyword(measure_build, "integration_q_interval", '["0", "1"]'),
        "cutoff_sc_exact_tuple": _has_call_keyword(measure_build, "source_sc_exact_mpf_tuple", "list(bits)"),
        "cutoff_sc_positive_margin": _has_expression(measure_build, "bits[1] > 0"),
        "cutoff_sc_small_margin": _has_expression(measure_build, "bits[2] + bits[3] <= -10"),
        "cutoff_hb_positive_less_half": _has_constant_containing(measure_build, "0<hb<1/2"),
        "cutoff_measure_upper_bound": _has_expression(measure_build, "c.mpf(3) * sc / 32"),
        "cutoff_eta_nonnegative_proof": _has_constant_containing(measure_build, "eta>=0"),
        "tail_hydrate_function": isinstance(tail_hydrate, ast.FunctionDef),
        "tail_N1_majorant_formula": _has_statement(
            tail_rows,
            "recomputed_iv = mbase_iv * beta_iv ** (WORD_DEPTH + 1) / (1 - beta_iv)",
        ),
        "tail_raw_majorant_scope_source": _has_constant(
            tail_hydrate, "uniform N=1 raw core/first-collar tail majorant"
        ),
        "tail_majorant_not_field_source": _has_constant_containing(
            tail_hydrate, "does not evaluate a signed first-chart field"
        ),
        "coupled_preconditioned_tail_formula": _has_expression(
            coupled_tail, "Mbase*beta**(N+1)/(1-beta)"
        ) and _has_constant(coupled_tail, "Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!"),
        "coupled_all_derivative_words": _has_call_keyword(
            coupled_tail, "all_k_derivative_words_allowed", "True"
        ) and _has_call_keyword(coupled_tail, "old_unpreconditioned_tail_reused", "False"),
        "coupled_tail_uses_N1": _has_statement(
            coupled_owner_tail,
            "tail=preconditioned_tail(c,view['initial_Gg_common_tube_bound'],c.sqrt(rho_max),rows['Cp'],rows['C1'],eta/8,1)",
        ),
        "coupled_tail_broadcasts_one_error_to_each_field": _has_statement(
            coupled_evaluate,
            "errors=[c.mpf([-ep(v)[1],ep(v)[1]]) for v in (error,derivative_error)]",
        ) and _has_statement(
            coupled_evaluate,
            "series={name:total[name]+algebra.lift(op.Jet(c,errors)) for name in FIELDS}",
        ),
    }
    require(all(checks.values()), "source AST guard failed")
    return checks


def _tuple4(value: Any, label: str) -> tuple[int, int, int, int]:
    require(isinstance(value, (list, tuple)) and len(value) == 4, f"{label} must be a precise four-entry MPF tuple")
    result = tuple(int(item) for item in value)
    require(result[0] in (0, 1) and result[1] >= 0 and result[3] >= 0, f"invalid {label}")
    return result  # type: ignore[return-value]


def _short_scientific(value: Sequence[int]) -> str:
    """Short display only; never expands a huge binary exponent."""

    sign, mantissa, exponent, _bits = _tuple4(value, "display endpoint")
    if mantissa == 0:
        return "0"
    log10_value = math.log10(mantissa) + exponent * math.log10(2.0)
    decimal_exponent = math.floor(log10_value)
    leading = 10.0 ** (log10_value - decimal_exponent)
    return f"{'-' if sign else ''}{leading:.8g}e{decimal_exponent:+d}"


def _record(value: Any, label: str) -> dict[str, Any]:
    require(hasattr(value, "_mpi_"), f"{label} is not an interval")
    lower = _tuple4(value._mpi_[0], label + ".lower")
    upper = _tuple4(value._mpi_[1], label + ".upper")
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "lower_scientific_display": _short_scientific(lower),
        "upper_scientific_display": _short_scientific(upper),
        "directed": True,
    }


def _mp_tuple_record(value: Any, side: str, label: str) -> tuple[int, int, int, int]:
    packet = mapping(value, label)
    return _tuple4(packet.get(f"{side}_exact_mpf_tuple"), f"{label}.{side}")


def _interval_context() -> MPIntervalContext:
    context = MPIntervalContext()
    context.prec = PRECISION_BITS
    return context


def _tail_contract(tail: Mapping[str, Any], family: Any, source_ast: Mapping[str, Any]) -> dict[str, Any]:
    normalized = mapping(tail.get("normalized_error_mapping", {}), "actual_selected_raw_word_tail.normalized_error_mapping")
    rows = mapping(tail.get("raw_tail_enclosures"), "actual_selected_raw_word_tail.raw_tail_enclosures")
    checks = mapping(tail.get("checks"), "actual_selected_raw_word_tail.checks")
    require(
        (checks.get("passed") is True if "passed" in checks else bool(checks) and all(value is True for value in checks.values())),
        "saved raw tail checks did not pass",
    )
    required = {
        "passed": tail.get("passed") is True,
        "source_only": tail.get("source_only") is True,
        "source_family": tail.get("source_family") == family,
        "N1_bound_hydrated": tail.get("actual_N1_tail_bound_hydrated") is True,
        "natural_word_depth_N1": tail.get("natural_word_depth") == 1,
        "majorant_not_value": tail.get("tail_is_a_directed_majorant_not_a_field_value") is True,
        "no_numeric_field": tail.get("numeric_actual_first_chart_field_evaluated") is False,
        "no_signed_point": tail.get("signed_point_value_evaluated") is False,
        "raw_tail_only": normalized.get("raw_tail_only") is True,
        "normalization_not_applied": normalized.get("applied") is False,
        "ordinary_Z_orders": tail.get("ordinary_Z_orders") == [0, 1],
    }
    require(all(required.values()), "saved raw source-bound N=1 tail contract failed")
    for order in ORDERS:
        row = mapping(rows.get(str(order), rows.get(order)), f"raw_tail_enclosures[{order}]")
        require(row.get("directed_upper_available") is not False, f"raw tail Z{order} is not directed")
        _mp_tuple_record(row.get("tail_upper"), "lower", f"raw tail Z{order}.tail_upper")
        _mp_tuple_record(row.get("tail_upper"), "upper", f"raw tail Z{order}.tail_upper")
    arithmetic = mapping(tail.get("arithmetic"), "actual_selected_raw_word_tail.arithmetic")
    require(
        tail.get("scope") == "uniform N=1 raw core/first-collar tail majorant"
        and tail.get("tail_source_sample") == "actual_first_keep tail constants only; saved field point excluded"
        and arithmetic.get("natural_word_depth") == 1
        and arithmetic.get("retained_preconditioned_indices") == [0, 1]
        and arithmetic.get("formula") == "Mbase_upper * beta_N_upper^(N+1) / (1-beta_N_upper)",
        "saved raw tail scope/norm semantics changed",
    )
    require(
        source_ast.get("tail_hydrate_function") is True
        and source_ast.get("tail_N1_majorant_formula") is True
        and source_ast.get("tail_raw_majorant_scope_source") is True
        and source_ast.get("tail_majorant_not_field_source") is True,
        "selected raw tail source AST contract changed",
    )
    pins = mapping(tail.get("source_pins"), "actual_selected_raw_word_tail.source_pins")
    coupled_pin = mapping(pins.get("coupled_source"), "actual_selected_raw_word_tail.source_pins.coupled_source")
    require(
        coupled_pin.get("path") == COUPLED_SOURCE
        and coupled_pin.get("sha256_match") is True
        and coupled_pin.get("expected_sha256") == FIXED[COUPLED_SOURCE],
        "saved raw tail coupled-source pin is not accepted",
    )
    require(
        source_ast.get("coupled_preconditioned_tail_formula") is True
        and source_ast.get("coupled_all_derivative_words") is True
        and source_ast.get("coupled_tail_uses_N1") is True
        and source_ast.get("coupled_tail_broadcasts_one_error_to_each_field") is True,
        "coupled first-step source does not prove componentwise raw tail domination",
    )
    required["source_AST_N1_majorant"] = True
    required["coupled_source_pin"] = True
    required["raw_vector_tail_dominates_each_field"] = True
    return required


def _prefix_receipt_contract(report: Mapping[str, Any], enclosure: Mapping[str, Any], source_ast: Mapping[str, Any]) -> bool:
    rows = enclosure.get("selected_coordinate_rows")
    require(isinstance(rows, list) and rows, "saved prefix rows are missing")
    q1_rows = [mapping(row, "selected row") for row in rows if _parse_q(mapping(row, "selected row").get("q"), "selected row.q") == 1]
    require(len(q1_rows) == 1, "saved prefix receipt does not have one q=1 row")
    ordinary = mapping(q1_rows[0].get("ordinary_Z"), "saved q1 ordinary_Z")
    for order in ORDERS:
        order_row = mapping(ordinary.get(str(order), ordinary.get(order)), f"saved q1 ordinary_Z[{order}]")
        for field in FIELDS:
            sector_sets: dict[str, set[tuple[Any, ...]]] = {}
            for component in ("W0", "W1", "prefix"):
                component_map = mapping(order_row.get(component), f"saved q1 Z{order}.{component}")
                terms = component_map.get(field)
                require(isinstance(terms, list) and terms, f"saved q1 Z{order}.{component}.{field} is empty")
                keys: set[tuple[Any, ...]] = set()
                for term in terms:
                    term_map = mapping(term, f"saved q1 Z{order}.{component}.{field} term")
                    exponents = term_map.get("source_exponents")
                    require(isinstance(exponents, list) and len(exponents) == 4, "saved prefix source sector changed")
                    keys.add(tuple(exponents))
                sector_sets[component] = keys
            require(sector_sets["prefix"] == sector_sets["W0"] | sector_sets["W1"], "saved prefix sectors are not W0 union W1")
    require(source_ast.get("build_prefix_is_W0_plus_W1") is True,
            "selected-word build AST does not prove prefix=W0+W1")
    return True


def _section_checks(report: Mapping[str, Any], section: Mapping[str, Any], name: str) -> dict[str, Any]:
    root_checks = mapping(report.get("checks"), "saved selected-word report.checks")
    require(root_checks.get("passed") is True, "saved selected-word report checks did not pass")
    section_checks = section.get("checks")
    if section_checks is None:
        # The source producer stores its direct ``passed`` result in the
        # section; the enclosing report still supplies the saved checks.
        require(section.get("passed") is True, f"{name} lacks passed source receipt")
        return {"enclosing_checks_passed": True, "section_passed": True}
    checks = mapping(section_checks, name + ".checks")
    require(
        (checks.get("passed") is True if "passed" in checks else bool(checks) and all(value is True for value in checks.values())),
        name + " checks did not pass",
    )
    return {"enclosing_checks_passed": True, "section_checks_passed": True}


def _parse_q(value: Any, label: str) -> Fraction:
    require(isinstance(value, str), label + " must be an exact rational string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(label + " is not an exact rational") from exc
    return result


def _selected_rows(report: Mapping[str, Any], enclosure: Mapping[str, Any]) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]], Any]:
    rows = enclosure.get("selected_coordinate_rows")
    require(isinstance(rows, list) and rows, "selected first-word rows missing")
    q_values = [_parse_q(mapping(row, f"selected row {index}").get("q"), f"selected row {index}.q") for index, row in enumerate(rows)]
    require(min(q_values) == 0 and max(q_values) == 1, "saved selected q rows do not span exact [0,1]")
    q1_indices = [index for index, q in enumerate(q_values) if q == 1]
    require(len(q1_indices) == 1, "saved selected rows do not contain a unique q=1 row")
    selected_z = enclosure.get("selected_Z_interval")
    require(selected_z is not None, "actual_selected_first_word_enclosure.selected_Z_interval is missing")
    prepared: list[Mapping[str, Any]] = []
    for index, row in enumerate(rows):
        row_map = dict(mapping(row, f"selected row {index}"))
        if q_values[index] == 1:
            row_map["selected_Z_interval"] = selected_z
        prepared.append(row_map)
    q1_rows = [prepared[q1_indices[0]]]
    return [mapping(row, f"prepared selected row {index}") for index, row in enumerate(prepared)], q1_rows, selected_z


def _check_enclosure(
    report: Mapping[str, Any], enclosure: Mapping[str, Any], family: Any, q1_rows: Sequence[Mapping[str, Any]], selected_z: Any,
    context: MPIntervalContext, source_ast: Mapping[str, Any],
) -> dict[str, Any]:
    require(enclosure.get("passed") is True, "saved selected first-word enclosure did not pass")
    require(enclosure.get("source_family") == family, "selected first-word source family changed")
    required = {
        "actual_finite_words": enclosure.get("actual_selected_finite_word_enclosures_evaluated") is True,
        "actual_sector_radii": enclosure.get("actual_directed_sector_coefficients_and_local_radii_evaluated") is True,
        "ordinary_Z_0_1": enclosure.get("ordinary_Z_orders") == [0, 1],
        "all_Z_not_claimed": enclosure.get("all_Z_functional_coverage") is False,
        "positive_exponential_not_materialized": enclosure.get("physical_exponential_source_factors_materialized") is False,
        "tail_separate": enclosure.get("raw_allword_tail_kept_separate") is True,
        "source_interval_errors_retained": enclosure.get("core_quadrature_and_source_interval_errors_retained") is True,
        "old_point_values_not_used": enclosure.get("old_actual_first_keep_point_values_used") is False,
        "old_core_not_promoted": enclosure.get("old_core_polynomial_promoted_to_word_prefix") is False,
        "source_not_executed": enclosure.get("source_executed") is False,
        "ancestors_not_replayed": enclosure.get("accepted_ancestors_replayed") is False,
        "local_integral_not_evaluated": enclosure.get("signed_local_collar_integral_evaluated") is False,
        "C_inner_not_evaluated": enclosure.get("signed_C_inner_m15_value_evaluated") is False,
    }
    require(all(required.values()), "saved selected first-word scope flags failed")
    required["prefix_receipt_is_W0_plus_W1"] = _prefix_receipt_contract(report, enclosure, source_ast)

    tube_equal = True
    local_nonnegative = True
    for row_index, row in enumerate(q1_rows):
        require(row.get("q") == "1", "prepared q1 row changed")
        ordinary = mapping(row.get("ordinary_Z"), "q1 ordinary_Z")
        for order in ORDERS:
            order_row = mapping(ordinary.get(str(order), ordinary.get(order)), f"q1 ordinary_Z[{order}]")
            require(order_row.get("allword_tail_included") is False, f"q1 Z{order} already includes all-word tail")
            radii = mapping(order_row.get("local_raw_error_radii"), f"q1 Z{order}.local_raw_error_radii")
            for component in ("W0", "W1", "prefix"):
                component_map = mapping(radii.get(component), f"q1 Z{order}.local_raw_error_radii.{component}")
                for field in ALL_FIELDS:
                    value = component_map.get(field)
                    if value is None:
                        local_nonnegative = False
                        continue
                    lower = _mp_tuple_record(value, "lower", f"q1 Z{order} {component}.{field}")
                    iv = context.make_mpf((lower, lower))
                    if not bool(iv.a >= 0):
                        local_nonnegative = False
        row_tube = row.get("selected_Z_interval")
        tube_equal = tube_equal and row_tube == selected_z
    require(tube_equal and local_nonnegative, "q1 saved Z tube or local-radius positivity changed")

    proof = enclosure.get("raw_local_error_proof")
    if proof is None:
        proof = report.get("raw_local_error_proof")
    proof = mapping(proof, "raw_local_error_proof")
    require(
        proof.get("local_W0") == "[dx*G5,dx*G6,dx*G6,G4,G5,G6]"
        and proof.get("local_W1") == "[dx^2*H,dx^2*H,dx^2*H,dx*H,dx*H,dx*H]"
        and proof.get("H") == "Cp*Mbase+C1*Mbase/loss on Z radius eta/8",
        "saved local-radius proof structure changed",
    )
    return required | {
        "q1_unique": True,
        "q1_allword_tail_excluded": True,
        "selected_Z_interval_injected": True,
        "local_radii_nonnegative": True,
    }


def _measure_view(report: Mapping[str, Any]) -> Mapping[str, Any]:
    for key in (
        "actual_cutoff_complete_support_measure",
        "actual_selected_cutoff_measure",
        "actual_cutoff_measure",
        "cutoff_measure",
        "actual_n1_cutoff_measure",
    ):
        value = report.get(key)
        if isinstance(value, Mapping):
            return value
    # Some accepted packets publish the measure fields at the report root.
    require("eta" in report and "source_sc_exact_mpf_tuple" in report, "saved cutoff measure section missing")
    return report


def _cutoff_measure_contract(report: Mapping[str, Any], family: Any, context: MPIntervalContext) -> tuple[dict[str, Any], Any]:
    measure = _measure_view(report)
    passed = measure.get("checks", {}).get("passed") if isinstance(measure.get("checks"), Mapping) else measure.get("passed")
    require(passed is True, "saved cutoff measure did not pass")
    require(measure.get("source_family") == family, "cutoff measure source family changed")
    require(measure.get("eta") == "1-sigma(4*q-2)", "cutoff eta definition changed")
    require(measure.get("integration_q_interval") == ["0", "1"], "cutoff q integration interval changed")
    require(measure.get("eta_support_q") == ["0", "3/4"], "cutoff eta support changed")
    require(measure.get("source_symmetry") == "sigma(t)+sigma(1-t)=1, including flat support endpoints", "cutoff reflection proof missing")
    require(measure.get("measure_is_Z_independent") is True, "cutoff measure is not Z independent")
    require(measure.get("complete_first_collar_support_included") is True, "complete eta support not included")
    require(measure.get("source_executed") is False and measure.get("accepted_ancestors_replayed") is False, "cutoff measure replayed a source")
    proof = str(measure.get("inequality_proof", ""))
    require("0<=q<=3/4" in proof and "0<hb<1/2" in proof and "eta>=0" in proof, "cutoff sc/hb/eta proof missing")
    bits = _tuple4(measure.get("source_sc_exact_mpf_tuple"), "source_sc_exact_mpf_tuple")
    sc = context.make_mpf((bits, bits))
    require(bool(sc.a > 0) and bool(sc.b < context.mpf(1) / 1024), "saved sc is not positive microscopic width")
    contract = {
        "passed": True,
        "eta_formula": True,
        "eta_support_q_0_to_3_over_4": True,
        "eta_in_unit_interval_from_reflection_and_nonnegative_proof": True,
        "sc_positive_and_small": True,
        "hb_positive_less_half_from_saved_proof": True,
        "measure_Z_independent": True,
        "complete_support": True,
        "signed_pressure_not_zeroed": measure.get("scalar_weight_bound_is_not_signed_m15_value") is True,
    }
    require(contract["signed_pressure_not_zeroed"], "cutoff measure was conflated with signed pressure/m15")
    return contract, sc


def _source_binding(
    report: Mapping[str, Any], enclosure: Mapping[str, Any], tail: Mapping[str, Any], family: Any,
    enclosure_flags: Mapping[str, Any], tail_flags: Mapping[str, Any], selected_z: Any,
) -> tuple[dict[str, Any], dict[str, str]]:
    rows = enclosure.get("selected_coordinate_rows")
    same_tube = isinstance(rows, list) and all(mapping(row, "selected row").get("selected_Z_interval", selected_z) == selected_z for row in rows if _parse_q(mapping(row, "selected row").get("q"), "selected row.q") == 1)
    values = {
        "same_source_family": enclosure.get("source_family") == family and tail.get("source_family") == family,
        "same_raw_first_function": enclosure_flags.get("actual_finite_words") is True and enclosure_flags.get("source_not_executed") is True,
        "same_raw_core_entry": enclosure_flags.get("source_interval_errors_retained") is True and enclosure_flags.get("old_point_values_not_used") is True,
        "same_raw_all_word_tail": enclosure_flags.get("tail_separate") is True and tail_flags.get("raw_tail_only") is True,
        "raw_vector_tail_dominates_each_field": (
            tail_flags.get("raw_vector_tail_dominates_each_field") is True
            and tail_flags.get("coupled_source_pin") is True
            and tail_flags.get("source_AST_N1_majorant") is True
        ),
        "ordinary_Z_rows_same_source": enclosure_flags.get("ordinary_Z_0_1") is True and enclosure.get("source_family") == family,
        "same_saved_Z_tube": same_tube and selected_z is not None,
        "prefix_is_homogeneous_not_function_value": enclosure_flags.get("positive_exponential_not_materialized") is True and enclosure_flags.get("old_core_not_promoted") is True,
        "tail_is_majorant_not_function_value": tail_flags.get("majorant_not_value") is True and tail_flags.get("no_numeric_field") is True,
        "analysis_checks_not_function_values": enclosure_flags.get("local_integral_not_evaluated") is True and enclosure_flags.get("C_inner_not_evaluated") is True and tail_flags.get("no_signed_point") is True,
        "core_extension_not_used_as_collar": enclosure_flags.get("old_core_not_promoted") is True and enclosure_flags.get("old_point_values_not_used") is True,
    }
    evidence = {
        "same_source_family": "actual_selected_first_word_enclosure.source_family == actual_selected_raw_word_tail.source_family == report.source_family",
        "same_raw_first_function": "actual_selected_first_word_enclosure.actual_selected_finite_word_enclosures_evaluated and source_executed is False",
        "same_raw_core_entry": "core_quadrature_and_source_interval_errors_retained and old_actual_first_keep_point_values_used is False",
        "same_raw_all_word_tail": "raw_allword_tail_kept_separate and actual_selected_raw_word_tail.normalized_error_mapping.raw_tail_only",
        "raw_vector_tail_dominates_each_field": "saved coupled_source pin plus coupled preconditioned-tail AST broadcasts one norm error to every FIELDS component",
        "ordinary_Z_rows_same_source": "ordinary_Z_orders == [0, 1] on selected enclosure and source family",
        "same_saved_Z_tube": "q=1 selected_Z_interval equals the injected enclosure-level selected_Z_interval",
        "prefix_is_homogeneous_not_function_value": "physical_exponential_source_factors_materialized is False and old_core_polynomial_promoted_to_word_prefix is False",
        "tail_is_majorant_not_function_value": "tail_is_a_directed_majorant_not_a_field_value and numeric_actual_first_chart_field_evaluated is False",
        "analysis_checks_not_function_values": "signed local collar/C_inner are explicitly not evaluated and signed point is false",
        "core_extension_not_used_as_collar": "old core polynomial and old point values are explicitly excluded",
    }
    require(all(values.values()), "source-binding evidence is incomplete")
    return values, evidence


def _phase_monotonicity(
    report: Mapping[str, Any], enclosure: Mapping[str, Any], measure_flags: Mapping[str, Any], source_ast: Mapping[str, Any],
    q_rows: Sequence[Mapping[str, Any]], q1_rows: Sequence[Mapping[str, Any]], sc: Any, context: MPIntervalContext,
    enclosure_flags: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    q_values = [_parse_q(mapping(row, "selected row").get("q"), "selected row.q") for row in q_rows]
    q_bounds = min(q_values) == 0 and max(q_values) == 1
    q1 = q1_rows[0]
    q1_y = q1.get("y_upper")
    require(q1_y is not None, "q=1 transport row lacks y_upper")
    q1_y_iv = _interval_from_record(context, q1_y, "q1.y_upper")
    sc_over_8 = sc / 8
    y_exact = q1_y_iv.a == sc_over_8.a and q1_y_iv.b == sc_over_8.b
    proof = enclosure.get("raw_local_error_proof") or report.get("raw_local_error_proof")
    proof = mapping(proof, "raw_local_error_proof")
    eta_values = mapping(measure_flags.get("exact_selected_eta", {}), "exact_selected_eta") if isinstance(measure_flags.get("exact_selected_eta"), Mapping) else {}
    eta_samples = eta_values or mapping(_measure_view(report).get("exact_selected_eta", {}), "exact_selected_eta")
    eta_sample_ok = all(str(eta_samples.get(key)) in {"0", "1/2", "1"} for key in ("0", "1/2", "5/8", "3/4", "1"))
    local_ast = source_ast["local_dx_is_qsc_over_4"] and source_ast["local_G_is_positive_linear_qsc"] and source_ast["local_W0_radius_structure"] and source_ast["local_W1_radius_structure"] and source_ast["local_radius_boxes_are_positive"]
    values = {
        "q_interval_is_exact_0_to_1": q_bounds and _measure_view(report).get("integration_q_interval") == ["0", "1"],
        "eta_in_unit_interval": measure_flags["eta_in_unit_interval_from_reflection_and_nonnegative_proof"] and eta_sample_ok,
        "eta_support_through_3_over_4": measure_flags["eta_support_q_0_to_3_over_4"],
        "eta_Z_zero": _measure_view(report).get("measure_is_Z_independent") is True,
        "sc_positive": bool(sc.a > 0),
        "hb_positive_less_half": measure_flags["hb_positive_less_half_from_saved_proof"] and source_ast["cutoff_hb_positive_less_half"],
        "transport_q1_y_box_is_0_to_sc_over_8": source_ast["transport_y_is_box_of_y_upper"] and source_ast["build_y_upper_is_qsc_over_8"] and y_exact,
        "all_q_y_is_subset_of_q1_box": q_bounds and bool(sc.a > 0) and bool(sc.b > 0) and y_exact,
        "transport_is_Z_independent": source_ast["transport_is_Z_independent"],
        "local_radii_monotone_in_qsc": local_ast and proof.get("dx_upper") == "dx=2*expm1(hb*q*sc/8)<=q*sc/4" and proof.get("local_W0") == "[dx*G5,dx*G6,dx*G6,G4,G5,G6]" and proof.get("local_W1") == "[dx^2*H,dx^2*H,dx^2*H,dx*H,dx*H,dx*H]",
        "local_radii_nonnegative": enclosure_flags.get("local_radii_nonnegative") is True,
    }
    require(all(values.values()), "phase monotonicity evidence is incomplete")
    derivation = {
        "exact_selected_q_values": [str(value) for value in sorted(set(q_values))],
        "q_lower_bound": "0",
        "q_upper_bound": "1",
        "qsc_definition": "qsc=sc*q",
        "qsc_lower_bound": "0 because sc>0 and q>=0",
        "qsc_upper_bound": "sc because q<=1",
        "y_definition": "y_upper=qsc/8",
        "all_q_y_upper_bound": "0<=y_upper<=sc/8; q=1 attains the upper bound",
        "caller_premises_retained": ["source_sc_exact_mpf_tuple is exact and positive", "source AST proves 0<hb<1/2", "eta support is the saved [0,3/4] source support"],
        "not_a_function_value": True,
    }
    return values, derivation


def _interval_from_record(context: MPIntervalContext, value: Any, label: str) -> Any:
    lower = _mp_tuple_record(value, "lower", label)
    upper = _mp_tuple_record(value, "upper", label)
    interval = context.make_mpf((lower, upper))
    require(interval.a <= interval.b, f"reversed saved interval {label}")
    return interval


def _serialize(value: Any, label: str = "result") -> Any:
    if hasattr(value, "_mpi_"):
        return _record(value, label)
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            lower = _tuple4(value["lower_exact_mpf_tuple"], label + ".lower")
            upper = _tuple4(value["upper_exact_mpf_tuple"], label + ".upper")
            projected = dict(value)
            projected["lower_exact_mpf_tuple"] = list(lower)
            projected["upper_exact_mpf_tuple"] = list(upper)
            projected["lower_scientific_display"] = _short_scientific(lower)
            projected["upper_scientific_display"] = _short_scientific(upper)
            projected["directed"] = True
            return {str(key): _serialize(item, f"{label}.{key}") for key, item in projected.items()}
        return {str(key): _serialize(item, f"{label}.{key}") for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize(item, f"{label}[{index}]") for index, item in enumerate(value)]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"cannot serialize {label}: {type(value).__name__}")


def _consume(hashes: Mapping[str, str]) -> dict[str, Any]:
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry " + path)
        raw[path] = data

    # The frozen report is parsed only after the raw-byte hash check above.
    report = mapping(json.loads(raw[REPORT]), "selected-word report")
    enclosure = mapping(report.get("actual_selected_first_word_enclosure"), "actual_selected_first_word_enclosure")
    tail = mapping(report.get("actual_selected_raw_word_tail"), "actual_selected_raw_word_tail")
    family = report.get("source_family")
    require(isinstance(family, Mapping) and bool(family), "saved selected-word source_family missing")
    source_ast = _source_ast_guards(raw)
    _section_checks(report, enclosure, "actual_selected_first_word_enclosure")
    tail_flags = _tail_contract(tail, family, source_ast)

    context = _interval_context()
    prepared_rows, q1_rows, selected_z = _selected_rows(report, enclosure)
    enclosure_flags = _check_enclosure(report, enclosure, family, q1_rows, selected_z, context, source_ast)
    measure_contract, sc = _cutoff_measure_contract(report, family, context)
    # Add the exact source-side eta samples only after checking the measure map.
    measure_for_samples = _measure_view(report)
    measure_flags = dict(measure_contract)
    measure_flags["exact_selected_eta"] = measure_for_samples.get("exact_selected_eta", {})
    phase_flags, phase_derivation = _phase_monotonicity(
        report, enclosure, measure_flags, source_ast, prepared_rows, q1_rows, sc, context, enclosure_flags
    )
    source_binding, source_binding_evidence = _source_binding(
        report, enclosure, tail, family, enclosure_flags, tail_flags, selected_z
    )

    def interval(value: Any) -> Any:
        return _interval_from_record(context, value, "leaf input")

    def hull(values: Sequence[Any]) -> Any:
        values = list(values)
        require(values, "empty directed hull")
        lower = values[0]._mpi_[0]
        upper = values[0]._mpi_[1]
        for value in values[1:]:
            candidate_lower = value._mpi_[0]
            candidate_upper = value._mpi_[1]
            if context.make_mpf((candidate_lower, candidate_lower)).a < context.make_mpf((lower, lower)).a:
                lower = candidate_lower
            if context.make_mpf((candidate_upper, candidate_upper)).a > context.make_mpf((upper, upper)).a:
                upper = candidate_upper
        result = context.make_mpf((lower, upper))
        require(result.a <= result.b, "directed hull reversed")
        return result

    def nonnegative(value: Any) -> bool:
        return bool(value.a >= 0)

    spec = importlib.util.spec_from_file_location("n1_uniform_first_word_certificate_leaf", ROOT / LEAF)
    require(spec is not None and spec.loader is not None, "leaf import")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    # This is the sole consumer call.  In particular, no selected-word/tail
    # producer and no ancestor receipt is imported or replayed.
    certificate = module.assemble(
        q1_rows,
        tail,
        interval=interval,
        hull=hull,
        nonnegative=nonnegative,
        source_binding=source_binding,
        phase_monotonicity=phase_flags,
        fields=FIELDS,
        orders=ORDERS,
    )
    serialized = _serialize(certificate, "certificate")
    return {
        "status": "conditional_uniform_certificate",
        "saved_selected_words_report_sha256": hashes[REPORT],
        "source_family": dict(family),
        "source_AST": source_ast,
        "source_binding": source_binding,
        "source_binding_evidence": source_binding_evidence,
        "phase_monotonicity": phase_flags,
        "phase_monotonicity_derivation": phase_derivation,
        "cutoff_measure_contract": measure_contract,
        "selected_Z_interval": _serialize(selected_z, "selected_Z_interval"),
        "selected_q1_row_count": 1,
        "fields": list(FIELDS),
        "ordinary_Z_orders": list(ORDERS),
        "zero_after_cutoff_scope": {
            "fields": list(FIELDS),
            "eta_zero_region": "q in [3/4,1]",
            "pressure_zeroed": False,
            "pressure_reason": "The leaf is called only for F1 and Uz1; the cutoff measure is not a signed pressure or m15 value.",
        },
        "q_bounds_are_caller_conditioned": all(phase_flags.values()),
        "all_Z_functional_coverage": False,
        "terminal_value_claim": False,
        "signed_sectors_preserved": serialized.get("signed_sectors_are_not_function_values") is True,
        "positive_radii_separate_from_signed_sectors": serialized.get("positive_radii_are_not_signed_increments") is True,
        "raw_all_word_tail_separate": source_binding.get("same_raw_all_word_tail") is True,
        "certificate": serialized,
    }


def main() -> int:
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(
        ENTRY,
        {
            "component": COMPONENT,
            "entered": True,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "consumer_calls": 0,
        },
    )
    out: dict[str, Any] = {
        "component": COMPONENT,
        "source_sha256": source_sha,
        "entry_sha256": entry_sha,
        "status": "failed",
        "consumer_calls": 0,
    }
    try:
        out["consumer_calls"] = 1
        out["proof"] = _consume(hashes)
        out["status"] = "returned"
        out["error"] = None
    except BaseException as exc:  # preserve the one-entry result even on a failed guard
        out["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha = exclusive(RESULT, out)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "status": out["status"],
                "source_sha256": source_sha,
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
            }
        )
    )
    return 0 if out["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
