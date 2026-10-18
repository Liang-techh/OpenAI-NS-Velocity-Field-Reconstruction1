"""Compose the reviewed native-control directed repairs.

This coordinator is deliberately a small integration boundary.  It reads the
accepted directed-repair report as data, binds the new common inputs once, and
then passes those inputs through the new long-reference and restoration
providers.  It never calls an accepted stage builder, inspector, self-check,
quadrature, or fixture.  The only fixture allowed here is the new closed-form
fixture owned by ``original_native_directed_common``.

The providers return raw interval trees.  Before the packet leaves this
module every interval is converted recursively with ``core.record`` so that a
runner can inspect the saved report without importing any mathematical
provider.
"""
from __future__ import annotations

import ast
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp


BASE = Path("experiments/root_st073")
PARENT_NAME = (
    "lei_ren_part1_paper_compliant_current_original_native_control_directed_repair.json"
)
PARENT_PATH = (BASE / PARENT_NAME).as_posix()
PARENT_SHA256 = "2ba16bd7c1bfd2fae90a5f13c69736d2764cccc049f8737e653fd10f9f8dbf3c"
PARENT_BYTES = 39425302

KEY = "actual_native_directed_controls_v2"
SCHEMA = "original_native_directed_controls_v2_v1"
MODULE = "src/openai_ns_reconstruction/original_native_directed_controls_v2.py"
RUNNER = "scripts/run_original_native_directed_controls_v2.py"
COMMON = "src/openai_ns_reconstruction/original_native_directed_common.py"
LONG = "src/openai_ns_reconstruction/original_long_reference_directed_controls.py"
RESTORATION = "src/openai_ns_reconstruction/original_restoration_directed_controls.py"
INTERVAL_PRECISION_BITS = 1600
ORDERS = tuple(range(5))

# These are intentionally pending until the source review is complete.  A
# pending value is a hard failure in ``build`` and in the runner, so a draft
# cannot be mistaken for a reviewed publication.
REVIEWED_SOURCE_HASHES = {
    COMMON: "5aa14af159faba82da6e9a73c21c49f83ec9003bf56cd40474301908c4882cb6",
    LONG: "29b5d50c4ea558a08974eb42e3a5f2746b121f9d2565685f80b5d726d7861116",
    RESTORATION: "6255fd90689e2557f5c1493087cdcc643d2ce5e1fec4ad61c128f989d94e4d2e",
}


def require(value: Any, message: str) -> None:
    if not value:
        raise ValueError(message)


def sha(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _compact(value: str) -> str:
    return "".join(value.split())


def _safe_relative(path: str) -> str:
    candidate = Path(path)
    require(not candidate.is_absolute(), "Absolute source pin path")
    require(".." not in candidate.parts, "Parent traversal in source pin path")
    return candidate.as_posix()


def _read_parent(root: Path, parent: Mapping[str, Any]) -> None:
    path = root / PARENT_PATH
    require(path.is_file(), "Directed-controls parent report is missing")
    raw = path.read_bytes()
    require(len(raw) == PARENT_BYTES, "Directed-controls parent byte count changed")
    require(hashlib.sha256(raw).hexdigest() == PARENT_SHA256,
            "Directed-controls parent bytes changed")
    require(json.loads(raw) == parent, "Build parent mapping is not canonical")


def _pin_path_for_input(name: str) -> str:
    candidate = Path(name)
    if candidate.parts and candidate.parts[0] == "experiments":
        return candidate.as_posix()
    # One inherited safe-helper pin is stored by the accepted producer as a
    # traversal from ``experiments/root_st073``.  Preserve that exact key so
    # the new packet can freeze the accepted input map without rewriting it.
    if ".." in candidate.parts:
        return candidate.as_posix()
    return (BASE / candidate.name).as_posix()


def _parent_pin_records(root: Path, parent: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Return every accepted input and implementation pin in canonical form."""
    records: dict[str, dict[str, Any]] = {}
    for name, expected in parent.get("input_hashes", {}).items():
        path = _pin_path_for_input(name)
        source = root / BASE / name if ".." in Path(name).parts else root / path
        actual = sha(source)
        require(actual == expected, "Changed pinned input: " + name)
        records[path] = {
            "path": path,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "passed": True,
            "source_executed": False,
        }
    for path, expected in parent.get("implementation_hashes", {}).items():
        records[_safe_relative(path)] = _pin(root, _safe_relative(path), expected)
    return records


def _check_parent_pins(root: Path, parent: Mapping[str, Any]) -> int:
    count = 0
    for name, expected in parent.get("input_hashes", {}).items():
        path = root / BASE / name if ".." in Path(name).parts else root / _pin_path_for_input(name)
        require(path.is_file(), "Pinned input is missing: " + str(path))
        require(sha(path) == expected, "Pinned input changed: " + str(path))
        count += 1
    for relative, expected in parent.get("implementation_hashes", {}).items():
        path = root / _safe_relative(relative)
        require(path.is_file(), "Pinned implementation is missing: " + relative)
        require(sha(path) == expected, "Pinned implementation changed: " + relative)
        count += 1
    require(count == len(parent.get("input_hashes", {}))
            + len(parent.get("implementation_hashes", {})),
            "Not all accepted source pins were checked")
    return count


def _pin(root: Path, relative: str, expected: str | None = None) -> dict[str, Any]:
    relative = _safe_relative(relative)
    actual = sha(root / relative)
    if expected is not None:
        require(actual == expected, "Changed pinned source: " + relative)
    return {
        "path": relative,
        "expected_sha256": actual if expected is None else expected,
        "actual_sha256": actual,
        "passed": True,
        "source_executed": False,
    }


def _reviewed_pins(root: Path) -> dict[str, dict[str, Any]]:
    pins: dict[str, dict[str, Any]] = {}
    for relative, expected in REVIEWED_SOURCE_HASHES.items():
        require(expected != "REVIEW_PENDING",
                "Independent source review is pending: " + relative)
        pins[relative] = _pin(root, relative, expected)
    return pins


def _raw_endpoint(value: Any) -> Any:
    """Read one exact MPF tuple without converting through decimal text."""
    require(isinstance(value, (list, tuple)) and len(value) == 4,
            "Malformed MPF endpoint tuple")
    require(all(type(part) is int for part in value),
            "MPF endpoint components must be integers")
    sign, mantissa, exponent, bit_count = value
    require(sign in (0, 1) and mantissa >= 0,
            "Malformed finite MPF sign or mantissa")
    require(bit_count == mantissa.bit_length(),
            "Malformed finite MPF bit count")
    if mantissa == 0:
        require(tuple(value) == (0, 0, 0, 0),
                "Noncanonical MPF zero tuple")
    else:
        require(mantissa & 1, "Noncanonical MPF mantissa")
    return tuple(value)


def _record(value: Any, label: str = "interval") -> int:
    """Validate one serialized finite directed interval and return one count."""
    require(isinstance(value, Mapping) and value.get("directed") is True,
            "Directed interval record required: " + label)
    lower = _raw_endpoint(value.get("lower_exact_mpf_tuple"))
    upper = _raw_endpoint(value.get("upper_exact_mpf_tuple"))
    require(mpf_cmp(lower, upper) <= 0, "Reversed directed interval: " + label)
    return 1


def _record_count(value: Any, label: str = "packet") -> int:
    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            return _record(value, label)
        return sum(_record_count(item, label + "." + str(key))
                   for key, item in value.items())
    if isinstance(value, (list, tuple)):
        return sum(_record_count(item, label + "[" + str(index) + "]")
                   for index, item in enumerate(value))
    return 0


def _serialize(value: Any, core: Any, label: str = "raw") -> Any:
    """Serialize an interval tree using ``core.record`` at every interval."""
    if hasattr(value, "_mpi_"):
        return core.record(value)
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            require(isinstance(key, (str, int, float, bool)) or key is None,
                    "Non-JSON mapping key in " + label)
            result[str(key) if not isinstance(key, str) else key] = _serialize(
                item, core, label + "." + str(key)
            )
        return result
    if isinstance(value, (list, tuple)):
        return [_serialize(item, core, label + "[]") for item in value]
    require(value is None or isinstance(value, (str, int, float, bool)),
            "Unsupported raw interval value in " + label)
    return value


def _validate_ordered_lists(value: Any, label: str = "packet") -> None:
    """Require every explicit ordinary-Z row list to be ordered and unique."""
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            child = label + "." + key_text
            if key_text in {"ordinary_Z_orders", "ordinary_Z_order_list"}:
                require(isinstance(item, list)
                        and item == list(range(len(item))),
                        "Unordered ordinary-Z orders at " + child)
            if key_text == "p2_ordinary_Z_orders":
                require(item == [0, 1, 2, 3],
                        "p2 ordinary-Z orders must be four rows at " + child)
            if isinstance(item, list) and item and all(
                    isinstance(row, Mapping) and "ordinary_Z_order" in row
                    for row in item
            ):
                orders = [row["ordinary_Z_order"] for row in item]
                require(orders == list(range(len(item))),
                        "Unordered ordinary-Z records at " + child)
            _validate_ordered_lists(item, child)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _validate_ordered_lists(item, label + "[" + str(index) + "]")


def _ast_guard(root: Path, parent: Mapping[str, Any], section: str,
               name: str) -> dict[str, Any]:
    guard_group = parent.get(section, {}).get("source_guards", {})
    if isinstance(guard_group, Mapping):
        old = guard_group.get(name)
    elif isinstance(guard_group, list):
        matches = [row for row in guard_group
                   if isinstance(row, Mapping) and row.get("function") == name]
        old = matches[0] if len(matches) == 1 else None
    else:
        old = None
    require(isinstance(old, Mapping), "Missing parent source guard: " + section + "." + name)
    relative = old.get("path")
    function = old.get("function")
    snippets = old.get("required_text")
    require(isinstance(relative, str) and isinstance(function, str)
            and isinstance(snippets, list),
            "Malformed parent source guard: " + name)
    text = (root / _safe_relative(relative)).read_text(encoding="utf-8")
    tree = ast.parse(text, filename=relative)
    candidates = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        if node.name != function:
            continue
        segment = ast.get_source_segment(text, node) or ""
        if all(_compact(snippet) in _compact(segment) for snippet in snippets):
            candidates.append(node)
    require(len(candidates) == 1, "Source AST guard is ambiguous: " + name)
    node = candidates[0]
    ast_hash = hashlib.sha256(
        ast.dump(node, include_attributes=False).encode()
    ).hexdigest()
    source_hash = sha(root / _safe_relative(relative))
    require(ast_hash == old.get("AST_sha256"), "Changed source AST: " + name)
    require(source_hash == old.get("source_sha256", parent.get("implementation_hashes", {})
                                  .get(relative, source_hash)),
            "Changed guarded source bytes: " + name)
    return {
        "passed": True,
        "path": relative,
        "function": function,
        "required_text": list(snippets),
        "AST_sha256": ast_hash,
        "source_sha256": source_hash,
        "source_executed": False,
        "parent_guard_reused_as_data": True,
    }


def _source_guards(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    # The parent guard records are evidence.  The coordinator independently
    # parses the reviewed source and never imports or executes these modules.
    for section in ("actual_long_reference_native_controls",
                    "actual_restoration_native_controls"):
        inherited = parent.get(section, {})
        require(isinstance(inherited, Mapping)
                and inherited.get("source_executed") is not True
                and inherited.get("accepted_ancestors_replayed") is False,
                "Inherited source execution claim changed: " + section)
    specs = (
        ("full_p2", "actual_long_reference_native_controls", "full_p2"),
        ("functions_source", "actual_long_reference_native_controls", "actual_constructor"),
        ("bridge_fields", "actual_long_reference_native_controls", "generic_fields"),
        ("roots_source", "actual_long_reference_native_controls", "roots"),
        ("loop_source", "actual_long_reference_native_controls", "loop"),
        ("anchored_gradient", "actual_long_reference_native_controls", "anchored_gradient"),
        ("history_densities", "actual_long_reference_native_controls", "history_densities"),
        ("complete_histories", "actual_long_reference_native_controls", "complete_histories"),
        ("sigma_formula", "actual_long_reference_native_controls", "sigma_formula"),
        ("physical_geometry", "actual_long_reference_background_I1_numeric", "__init__"),
        ("restoration_functions", "actual_restoration_native_controls", "c6_chartmaps"),
        ("restoration_p2", "actual_restoration_native_controls", "q4_full_formula"),
    )
    result = {
        output: _ast_guard(root, parent, section, source_name)
        for output, section, source_name in specs
    }
    result["source_proof"] = {
        "passed": all(row.get("passed") is True for row in result.values()),
        "parent_guards_reused_as_data": True,
        "source_executed": False,
        "accepted_source_builders_called": False,
        "guarded_functions": [output for output, _, _ in specs],
    }
    require(result["source_proof"]["passed"], "Source proof guards failed")
    return result


def _source_premises(parent: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the saved source premises consumed by the corrected long worker."""
    section = parent.get("actual_long_reference_native_controls", {})
    identities = section.get("source_function_identities", {})
    sigma = identities.get("sigma_derivative_bound", {})
    require(isinstance(identities, Mapping) and isinstance(sigma, Mapping)
            and identities.get("passed") is True
            and sigma.get("passed") is True
            and sigma.get("value_range") == [0, 1]
            and sigma.get("derivative_range") == [0, 8]
            and sigma.get("argument") ==
                "t=(1-2q)^2, x=8*sqrt(t)/(1-t)^2, sigma'=8*(1+3t)/(1-t)^3/cosh(x)^2"
            and sigma.get("inequality") ==
                "cosh(x)^2>=1+x^2; 64-(6-3t+t^2)*(1-t)=58+9t-4t^2+t^3>0 for 0<=t<1"
            and sigma.get("positive_polynomial") == "t**3 - 4*t**2 + 9*t + 58"
            and sigma.get("flat_ends_by_definition") is True,
            "Saved sigma derivative source premise changed")
    history = section.get("history_rate_source", {})
    rates = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
    require(isinstance(history, Mapping)
            and history.get("passed") is True
            and history.get("source_executed") is False
            and history.get("source") ==
                "experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
            and history.get("bridge_alias") == "RATES=bridge.RATES"
            and history.get("rates") == rates,
            "Saved history-rate source premise changed")
    gradient = section.get("anchored_gradient", {})
    require(isinstance(gradient, Mapping)
            and gradient.get("ordinary_not_taylor") is True
            and len(gradient.get("G_Z1_through_Z4", [])) == 4,
            "Saved anchored-gradient source premise changed")
    sigma_receipt = {
        "passed": True,
        "value_range": list(sigma["value_range"]),
        "derivative_range": list(sigma["derivative_range"]),
        "argument": sigma.get("argument"),
        "inequality": sigma.get("inequality"),
        "positive_polynomial": sigma.get("positive_polynomial"),
        "flat_ends_by_definition": True,
        "source_executed": False,
    }
    history_receipt = {
        "passed": True,
        "source": history.get("source"),
        "bridge_alias": history.get("bridge_alias"),
        "rates": dict(history["rates"]),
        "active_native_rates": {name: history["rates"][name]
                                 for name in ("A", "B", "C")},
        "source_executed": False,
    }
    return {
        "passed": True,
        "anchored_gradient": {
            "passed": True,
            "ordinary_not_taylor": True,
            "row_count": 4,
            "source_executed": False,
        },
        "sigma_derivative_bound": sigma_receipt,
        "history_rate_source": history_receipt,
        "source_executed": False,
    }


def _validate_execution_flags(value: Any, label: str = "packet") -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key).lower()
            if key_text in {
                "source_executed", "accepted_ancestors_replayed",
                "accepted_builders_called", "accepted_builders_or_owners_replayed",
                "accepted_checkers_called", "accepted_checkers_replayed",
                "accepted_self_checks_called", "accepted_generic_self_checks_called",
                "accepted_recurrences_replayed", "accepted_quadratures_replayed",
                "old_stage_calls", "old_stage_called", "old_builders_called",
            }:
                require(item is False, "Forbidden execution flag at " + label + "." + str(key))
            _validate_execution_flags(item, label + "." + str(key))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _validate_execution_flags(item, label + "[" + str(index) + "]")


def _validate_component(value: Any, name: str) -> None:
    require(isinstance(value, Mapping), name + " component must be a mapping")
    # The long provider intentionally returns a raw chart dictionary without
    # a top-level pass flag.  Its sourceproof and limits are the admission
    # fields; restoration also carries the explicit pass flag.
    if value.get("passed") is not True:
        require(isinstance(value.get("charts"), Mapping)
                and isinstance(value.get("sourceproof"), Mapping)
                and value.get("sourceproof", {}).get("source_executed") is False,
                name + " component did not publish a verified raw chart packet")
    _validate_execution_flags(value, name)
    _validate_ordered_lists(value, name)
    # These are intentionally checked when present and never inferred from a
    # missing field.  A worker that has not published the field is rejected by
    # the fail-closed packet checks below.
    if "source_executed" in value:
        require(value["source_executed"] is False, name + " source was executed")


def _strict_negative(value: Any, label: str) -> bool:
    if hasattr(value, "_mpi_"):
        endpoints = value._mpi_
        require(isinstance(endpoints, (list, tuple)) and len(endpoints) == 2,
                "Malformed raw directed interval: " + label)
        lower = _raw_endpoint(endpoints[0])
        upper = _raw_endpoint(endpoints[1])
        require(mpf_cmp(lower, upper) <= 0,
                "Reversed raw directed interval: " + label)
    elif isinstance(value, Mapping) and "upper_exact_mpf_tuple" in value:
        _record(value, label)
        upper = tuple(value["upper_exact_mpf_tuple"])
    else:
        return False
    zero = (0, 0, 0, 0)
    return mpf_cmp(upper, zero) < 0


def _strict_controls(long_packet: Mapping[str, Any], restoration_packet: Mapping[str, Any]) -> dict[str, Any]:
    # These are the provider's fixed schema paths.  Avoid a broad name search
    # that could accidentally certify an unrelated diagnostic field.
    slope = long_packet.get("slope", {}).get("margin")
    require(_strict_negative(slope, "slope.margin"),
            "Actual strict slope margin is missing or nonnegative")
    restoration_chart = restoration_packet.get("charts", {}).get("restoration", {})
    bmargin = restoration_chart.get("bmargin", {}).get("log_abs_b_minus_log_0_1")
    require(_strict_negative(bmargin, "charts.restoration.bmargin.log_abs_b_minus_log_0_1"),
            "Actual negative b margin is missing or nonnegative")
    return {
        "passed": True,
        "strict_slope_records": ["raw_interval_tree.long_reference.slope.margin"],
        "negative_b_margin_records": [
            "raw_interval_tree.restoration.charts.restoration.bmargin."
            "log_abs_b_minus_log_0_1"
        ],
        "exact_tuple_comparisons": True,
    }


def _native_cap_rows(rows: Any, label: str, serialized: bool = False,
                     count: int = 5) -> str:
    """Validate one ordered raw or explicitly packed cap sequence."""
    require(isinstance(rows, (list, tuple)) and len(rows) == count,
            label + " must contain " + str(count) + " rows")
    packed = all(isinstance(row, Mapping) for row in rows)
    raw = all(row is None or hasattr(row, "_mpi_") for row in rows)
    require((packed and not raw) or (raw and not packed),
            label + " must be uniformly raw or explicitly packed")
    if raw:
        for index, row in enumerate(rows):
            if row is None:
                continue
            endpoints = row._mpi_
            require(isinstance(endpoints, (list, tuple)) and len(endpoints) == 2,
                    "Malformed raw cap row: " + label + "[" + str(index) + "]")
            lower = _raw_endpoint(endpoints[0])
            upper = _raw_endpoint(endpoints[1])
            require(mpf_cmp(lower, upper) <= 0,
                    "Reversed raw cap row: " + label + "[" + str(index) + "]")
        return "raw"
    required = {"ordinary_Z_order", "exact_zero", "log_absolute_upper_bound"}
    for index, row in enumerate(rows):
        require(set(row) == required,
                "Packed cap row keys are not explicit: " + label + "[" + str(index) + "]")
        require(row["ordinary_Z_order"] == index
                and row["exact_zero"] in (True, False),
                "Packed cap row order changed: " + label + "[" + str(index) + "]")
        bound = row["log_absolute_upper_bound"]
        if row["exact_zero"]:
            require(bound is None,
                    "Exact-zero packed cap has a bound: " + label + "[" + str(index) + "]")
        elif serialized:
            _record(bound, label + "[" + str(index) + "]")
        else:
            require(hasattr(bound, "_mpi_"),
                    "Packed cap bound is not raw: " + label + "[" + str(index) + "]")
    return "packed"


def _finite_interval_or_record(value: Any, label: str) -> None:
    if hasattr(value, "_mpi_"):
        endpoints = value._mpi_
        require(isinstance(endpoints, (list, tuple)) and len(endpoints) == 2,
                "Malformed finite interval: " + label)
        lower = _raw_endpoint(endpoints[0])
        upper = _raw_endpoint(endpoints[1])
        require(mpf_cmp(lower, upper) <= 0,
                "Reversed finite interval: " + label)
        return
    _record(value, label)


def _lower_positive(value: Any, label: str) -> None:
    if hasattr(value, "_mpi_"):
        endpoints = value._mpi_
        require(isinstance(endpoints, (list, tuple)) and len(endpoints) == 2,
                "Malformed positive interval: " + label)
        lower = _raw_endpoint(endpoints[0])
    else:
        _record(value, label)
        lower = tuple(value["lower_exact_mpf_tuple"])
    require(mpf_cmp(lower, (0, 0, 0, 0)) > 0,
            "Positive lower endpoint is not strict: " + label)


def _packed_upper_negative(row: Any, label: str) -> None:
    required = {"ordinary_Z_order", "exact_zero", "log_absolute_upper_bound"}
    require(isinstance(row, Mapping) and set(row) == required
            and row.get("exact_zero") is False,
            "Nonzero packed row required: " + label)
    bound = row["log_absolute_upper_bound"]
    if hasattr(bound, "_mpi_"):
        upper = _raw_endpoint(bound._mpi_[1])
    else:
        _record(bound, label)
        upper = tuple(bound["upper_exact_mpf_tuple"])
    require(mpf_cmp(upper, (0, 0, 0, 0)) < 0,
            "Packed log cap is not negative: " + label)


def _cut_numeric_receipt(context: MPIntervalContext, shared: Mapping[str, Any],
                         long_packet: Mapping[str, Any],
                         restoration_packet: Mapping[str, Any],
                         core: Any, common_module: Any) -> dict[str, Any]:
    """Bind a numeric lower bound for the active ``cut=1-Delta/eta`` branch."""
    rows = []
    for section_name, packet in (("long_reference", long_packet),
                                 ("restoration", restoration_packet)):
        charts = packet.get("charts", {})
        for chart_name, chart in charts.items():
            require(isinstance(chart, Mapping), "Cut chart is not a mapping")
            conditioning = chart.get("conditioning")
            roots = chart.get("root_caps")
            require(isinstance(conditioning, Mapping)
                    and isinstance(roots, Mapping)
                    and conditioning.get("Delta_negative") is True,
                    "Numeric cut source fields are missing: " + chart_name)
            eta_log = conditioning.get("eta_log", shared.get("logeta"))
            require(hasattr(eta_log, "_mpi_"),
                    "Numeric cut eta source is not raw: " + chart_name)
            kappa_rows = roots.get("kappa")
            require(isinstance(kappa_rows, (list, tuple)) and len(kappa_rows) == 4,
                    "Numeric cut kappa source is missing: " + chart_name)
            kappa_row = kappa_rows[0]
            require(isinstance(kappa_row, Mapping)
                    and kappa_row.get("exact_zero") is False,
                    "Numeric cut kappa source is zero: " + chart_name)
            kappa_log = kappa_row.get("log_absolute_upper_bound")
            require(hasattr(kappa_log, "_mpi_"),
                    "Numeric cut kappa cap is not raw: " + chart_name)
            # eta is positive by its exact source definition.  Its finite
            # logarithm is kept logarithmic even when eta is too small to
            # materialize.  eta<=1 and kappa<1 imply cut>2 algebraically.
            require(common_module.upper(eta_log) <= 0,
                    "Cut proof does not establish eta<=1")
            require(common_module.upper(kappa_log) < 0,
                    "Cut proof does not establish kappa<1")
            cut_lower = context.mpf(2)
            require(common_module.lower(cut_lower) > 1,
                    "Numeric cut lower bound is not above one: " + chart_name)
            rows.append({
                "section": section_name,
                "chart": chart_name,
                "formula": "cut=1-Delta/eta",
                "numeric_lower_bound": core.record(cut_lower),
                "eta_source_log": core.record(eta_log),
                "kappa_log_cap": core.record(kappa_log),
                "stored_logarithms_exponentiated": False,
                "derived_from": "kappa_upper<1, eta_upper<=1, Delta=-(2-kappa)",
                "Delta_negative": True,
                "source_executed": False,
            })
    require(len(rows) == 4, "Numeric cut receipt does not cover four charts")
    return {
        "passed": True,
        "formula": "cut=1-Delta/eta",
        "numeric_lower_bound_above_one": True,
        "rows": rows,
        "source_executed": False,
    }


def _canonical_public_rows(value: Any, count: int, label: str, common: Any) -> Any:
    """Pack raw provider rows once while preserving already packed rows."""
    require(isinstance(value, (list, tuple)) and len(value) == count,
            label + " must contain " + str(count) + " rows")
    raw = all(row is None or hasattr(row, "_mpi_") for row in value)
    required = {"ordinary_Z_order", "exact_zero", "log_absolute_upper_bound"}
    packed = all(isinstance(row, Mapping) and set(row) == required for row in value)
    require(raw or packed, label + " is neither raw nor explicitly packed")
    if raw:
        return common.pack(list(value))
    for index, row in enumerate(value):
        require(row["ordinary_Z_order"] == index
                and row["exact_zero"] in (True, False),
                "Packed row order changed: " + label + "[" + str(index) + "]")
    return value


def _canonicalize_public_packet(packet: Mapping[str, Any], common: Any) -> None:
    """Canonicalize public arrays after both workers have consumed raw rows."""
    charts = packet.get("charts")
    require(isinstance(charts, Mapping), "Provider charts are missing")
    if "ordinary_Z_orders" not in packet:
        require(all(chart.get("ordinary_Z_orders") == [0, 1, 2, 3, 4]
                    for chart in charts.values()), "Provider chart orders changed")
        packet["ordinary_Z_orders"] = [0, 1, 2, 3, 4]
    for chart_name, chart in charts.items():
        require(isinstance(chart, Mapping), "Provider chart is not a mapping")
        for field in ("phi_relative", "E_relative", "inverse_Utheta", "p2"):
            if field in chart:
                count = 4 if field == "p2" else 5
                chart[field] = _canonical_public_rows(
                    chart[field], count, str(chart_name) + "." + field, common
                )
        native = chart.get("native_caps")
        require(isinstance(native, Mapping),
                "Provider native caps are missing: " + str(chart_name))
        for field, rows in native.items():
            native[field] = _canonical_public_rows(
                rows, 5, str(chart_name) + ".native_caps." + str(field), common
            )
        root_caps = chart.get("root_caps")
        require(isinstance(root_caps, Mapping),
                "Provider root caps are missing: " + str(chart_name))
        for field, rows in root_caps.items():
            root_caps[field] = _canonical_public_rows(
                rows, 4, str(chart_name) + ".root_caps." + str(field), common
            )
        if "p2" in chart and "p2_ordinary_Z_orders" not in chart:
            chart["p2_ordinary_Z_orders"] = [0, 1, 2, 3]


def _shape_receipt(shared: Mapping[str, Any], long_packet: Mapping[str, Any],
                   restoration_packet: Mapping[str, Any]) -> dict[str, Any]:
    # The common provider publishes the independent five-row sources.  The
    # two chart providers must publish full Z0..Z4 rows and four-row p2 rows.
    for field in ("F2_logs5", "P0_logs5", "q_log_values5", "inverse_L_logs5"):
        rows = shared.get(field)
        require(isinstance(rows, list) and len(rows) == 5,
                "Common five-row field is missing: " + field)
    encoding_receipt = {}
    for name, packet in (("long", long_packet), ("restoration", restoration_packet)):
        orders = packet.get("ordinary_Z_orders")
        require(orders == [0, 1, 2, 3, 4],
                name + " full ordinary-Z orders are not five rows")
        charts = packet.get("charts")
        expected_charts = (
            {"long_reshape", "reference"}
            if name == "long" else {"restoration", "post_restore"}
        )
        require(isinstance(charts, Mapping) and set(charts) >= expected_charts,
                name + " chart set changed")
        for chart_name, chart in charts.items():
            require(chart.get("ordinary_Z_orders") == [0, 1, 2, 3, 4],
                    name + "." + str(chart_name) + " full rows are not five")
            p2_orders = chart.get("p2_ordinary_Z_orders")
            if p2_orders is not None:
                require(p2_orders == [0, 1, 2, 3],
                        name + "." + str(chart_name) + " p2 rows are not four")
            p2_rows = chart.get("p2")
            require(isinstance(p2_rows, (list, tuple)) and len(p2_rows) == 4,
                    name + "." + str(chart_name) + " p2 field is not four rows")
            for field in ("phi_relative", "E_relative", "inverse_Utheta"):
                rows = chart.get(field)
                require(isinstance(rows, (list, tuple)) and len(rows) == 5,
                        name + "." + str(chart_name) + " " + field
                        + " field is not five rows")
                _native_cap_rows(
                    rows, name + "." + str(chart_name) + "." + field
                )
            _native_cap_rows(
                p2_rows, name + "." + str(chart_name) + ".p2", count=4
            )
            native_caps = chart.get("native_caps")
            require(isinstance(native_caps, Mapping)
                    and set(native_caps) == {"phi", "U", "M", "A", "B", "C"},
                    name + "." + str(chart_name) + " native cap names changed")
            encodings = {}
            for field, rows in native_caps.items():
                encodings[field] = _native_cap_rows(
                    rows, name + "." + str(chart_name) + "." + field
                )
            encoding_receipt[name + "." + str(chart_name)] = encodings
            root_caps = chart.get("root_caps")
            root_names = {"q_squared", "q", "u", "h", "r", "alpha", "K"}
            require(isinstance(root_caps, Mapping) and set(root_caps) >= root_names,
                    name + "." + str(chart_name) + " root caps are incomplete")
            require("kappa" in root_caps,
                    name + "." + str(chart_name) + " kappa cap is missing")
            for root_name, rows in root_caps.items():
                _native_cap_rows(
                    rows, name + "." + str(chart_name) + ".root_caps." + str(root_name),
                    count=4
                )
            _packed_upper_negative(
                root_caps["kappa"][0],
                name + "." + str(chart_name) + ".root_caps.kappa[0]",
            )
            conditioning = chart.get("conditioning")
            required_conditioning = {
                "log_D_lower", "log_lambda_lower", "positive_a_lower",
                "positive_gamma_lower", "positive_q_squared_lower_log",
                "Delta_negative",
            }
            require(isinstance(conditioning, Mapping)
                    and set(conditioning) >= required_conditioning
                    and conditioning.get("Delta_negative") is True,
                    name + "." + str(chart_name) + " conditioning is incomplete")
            for field in ("log_D_lower", "log_lambda_lower", "positive_a_lower",
                          "positive_gamma_lower", "positive_q_squared_lower_log"):
                _finite_interval_or_record(
                    conditioning[field],
                    name + "." + str(chart_name) + ".conditioning." + field,
                )
            for field in ("positive_a_lower", "positive_gamma_lower"):
                _lower_positive(
                    conditioning[field],
                    name + "." + str(chart_name) + ".conditioning." + field,
                )
            if str(chart_name) == "post_restore":
                for field in ("b", "t0"):
                    rows = root_caps[field]
                    require(all(
                        isinstance(row, Mapping)
                        and row.get("exact_zero") is True
                        and row.get("log_absolute_upper_bound") is None
                        for row in rows
                    ), "Post-restore " + field + " rows are not exact zero")
                detail = chart.get("source_detail")
                require(isinstance(detail, Mapping)
                        and detail.get("physical_U_y") is True
                        and detail.get("b_formula") == "b=2*U_y/Utheta"
                        and detail.get("b_logs") == [None] * 5,
                        "Post-restore U_y=0 source receipt is missing")
    return {
        "passed": True,
        "common_five_row_fields": ["F2_logs5", "P0_logs5", "q_log_values5", "inverse_L_logs5"],
        "chart_full_rows": True,
        "chart_p2_four_rows": True,
        "native_cap_rows_raw_or_packed": True,
        "native_cap_encodings": encoding_receipt,
        "root_names_complete": True,
        "conditioning_checked": True,
    }


def _same_raw(left: Any, right: Any) -> bool:
    if hasattr(left, "_mpi_") or hasattr(right, "_mpi_"):
        return (hasattr(left, "_mpi_") and hasattr(right, "_mpi_")
                and tuple(left._mpi_[0]) == tuple(right._mpi_[0])
                and tuple(left._mpi_[1]) == tuple(right._mpi_[1]))
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        return (set(left) == set(right)
                and all(_same_raw(left[key], right[key]) for key in left))
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        return len(left) == len(right) and all(
            _same_raw(a, b) for a, b in zip(left, right)
        )
    return left == right


def _memory_rows(value: Any, label: str) -> Any:
    """Decode only the common packed five-row form for numeric comparison."""
    if isinstance(value, list) and value and all(
            isinstance(row, Mapping) and "log_absolute_upper_bound" in row
            for row in value
    ):
        orders = [row.get("ordinary_Z_order") for row in value]
        require(orders == list(range(len(value))),
                "Unordered packed memory rows: " + label)
        result = []
        for index, row in enumerate(value):
            require(row.get("exact_zero") in (True, False),
                    "Malformed packed memory row: " + label + "[" + str(index) + "]")
            bound = row.get("log_absolute_upper_bound")
            if row["exact_zero"]:
                require(bound is None, "Nonempty exact-zero memory row: " + label)
                result.append(None)
            else:
                require(hasattr(bound, "_mpi_"),
                        "Packed memory row is not a raw interval: " + label)
                result.append(bound)
        return result
    return value


def _memory_receipt(long_packet: Mapping[str, Any], restoration_packet: Mapping[str, Any]) -> dict[str, Any]:
    long_reference = long_packet.get("charts", {}).get("reference", {})
    restoration_chart = restoration_packet.get("charts", {}).get("restoration", {})
    post_chart = restoration_packet.get("charts", {}).get("post_restore", {})
    expected = long_reference.get("native_caps")
    actual = restoration_chart.get("incoming_memory")
    post_expected = restoration_chart.get("native_caps")
    post_actual = post_chart.get("incoming_memory")
    require(isinstance(expected, Mapping) and isinstance(actual, Mapping),
            "Restoration source-memory input rows are missing")
    require(set(expected) >= {"A", "B", "C"}
            and set(actual) >= {"A", "B", "C"},
            "Restoration source-memory fields are incomplete")
    require(all(_same_raw(_memory_rows(expected[key], "long_reference." + key),
                         _memory_rows(actual[key], "restoration." + key))
                for key in ("A", "B", "C")),
            "Restoration memory rows do not equal corrected long reference rows")
    require(isinstance(post_expected, Mapping) and isinstance(post_actual, Mapping)
            and all(_same_raw(_memory_rows(post_expected[key], "restoration." + key),
                              _memory_rows(post_actual[key], "post_restore." + key))
                    for key in ("A", "B", "C")),
            "Post-restore memory rows do not equal restoration rows")
    return {
        "passed": True,
        "source_packet": "long_reference_directed_controls",
        "restoration_input_receipt_present": True,
        "numeric_actual_rows_equal": True,
        "fields": ["A", "B", "C"],
        "post_restore_rows_equal": True,
        "source_executed": False,
    }


def _same_serialized_record(value: Any, expected: Any, label: str) -> bool:
    """Compare one raw provider interval with one saved directed record.

    The background radii are physical geometry data inherited from the
    accepted source.  They are compared through their exact MPF endpoint
    tuples; decimal display text and recomputation are deliberately excluded.
    """
    require(hasattr(value, "_mpi_"), "Physical radius is not a raw interval: " + label)
    endpoints = value._mpi_
    require(isinstance(endpoints, (list, tuple)) and len(endpoints) == 2,
            "Malformed physical radius interval: " + label)
    lower = _raw_endpoint(endpoints[0])
    upper = _raw_endpoint(endpoints[1])
    _record(expected, label)
    return (lower == tuple(expected["lower_exact_mpf_tuple"])
            and upper == tuple(expected["upper_exact_mpf_tuple"]))


def _physical_geometry_receipt(parent: Mapping[str, Any],
                               long_packet: Mapping[str, Any],
                               restoration_packet: Mapping[str, Any]) -> dict[str, Any]:
    """Bind unchanged physical radii to the saved background geometry.

    Only the physical ``log_R`` scales are reused here.  Native absolute-log
    caps remain newly computed controls and are never admitted as a source of
    the radius values.
    """
    background = parent.get("actual_long_reference_background_I1_numeric", {})
    old_restoration = parent.get("actual_restoration_native_controls", {})
    old_background_charts = background.get("charts", {})
    old_restoration_charts = old_restoration.get("charts", {})
    comparisons = (
        ("long_reshape", long_packet, ("charts", "long_reshape", "scale", "logR"),
         old_background_charts, ("long_reshape", "log_R_upper")),
        ("reference", long_packet, ("charts", "reference", "scale", "logR"),
         old_background_charts, ("reference", "log_R_upper")),
        ("restoration", restoration_packet,
         ("charts", "restoration", "scales", "log_R_upper"),
         old_restoration_charts, ("restoration", "log_R_upper")),
        ("post_restore", restoration_packet,
         ("charts", "post_restore", "scales", "log_R_upper"),
         old_restoration_charts, ("post_restore", "log_R_upper")),
    )
    receipt_rows = []
    for name, packet, provider_path, old_charts, old_path in comparisons:
        current: Any = packet
        for key in provider_path:
            require(isinstance(current, Mapping) and key in current,
                    "Physical radius provider path is missing: " + name)
            current = current[key]
        old_chart_name, old_field = old_path
        require(isinstance(old_charts, Mapping)
                and isinstance(old_charts.get(old_chart_name), Mapping),
                "Saved physical radius chart is missing: " + name)
        expected = old_charts[old_chart_name].get(old_field)
        require(_same_serialized_record(current, expected, "physical_geometry." + name),
                "Physical radius changed: " + name)
        receipt_rows.append({
            "chart": name,
            "provider_path": ".".join(
                (("long_reference",) if packet is long_packet else ("restoration",))
                + provider_path
            ),
            "saved_path": (
                ("actual_long_reference_background_I1_numeric.charts." +
                 old_chart_name + "." + old_field)
                if old_charts is old_background_charts else
                ("actual_restoration_native_controls.charts." +
                 old_chart_name + "." + old_field)
            ),
            "exact_endpoint_tuples_equal": True,
        })
    return {
        "passed": True,
        "unaffected_physical_geometry_reused": True,
        "saved_radii_exactly_equal": True,
        "unsafe_abs_native_cap_reused": False,
        "source_guard": "physical_geometry",
        "formula_source": {
            "long_reshape": "logR=ln(110)+T",
            "reference": "logR=ln(110)+logref-8",
        },
        "comparisons": receipt_rows,
        "source_executed": False,
    }


def _preservation_receipt(parent: Mapping[str, Any], long_packet: Mapping[str, Any],
                          restoration_packet: Mapping[str, Any]) -> dict[str, Any]:
    required = (
        "old_212_terms_and_errors_unchanged",
        "six_signed_integrals_retained",
        "incoming_histories_not_reset",
    )
    found: dict[str, bool] = {}
    old_restoration = parent.get("actual_restoration_native_controls", {})
    old_preservation = old_restoration.get("preservation", {})
    require(old_preservation.get("old_212_terms_and_errors_unchanged") is True,
            "Historical 212 terms and errors were not retained")
    require(old_preservation.get("six_signed_integrals_retained") is True,
            "Historical signed integrals were not retained")
    require((long_packet.get("sourceproof", {}).get("nonreset_histories") is True)
            and (restoration_packet.get("source_inputs", {})
                 .get("corrected_reference_A_B_C_memory") is True),
            "Corrected source histories were reset")
    found.update({
        "old_212_terms_and_errors_unchanged": True,
        "six_signed_integrals_retained": True,
        "incoming_histories_not_reset": True,
    })
    # The parent packet must retain its six signed terms and 212 cumulative
    # terms; report-level deep equality is enforced by the runner as well.
    require(parent.get("checks", {}).get("directed_records") == 256,
            "Directed-repair predecessor evidence changed")
    return {"passed": True, "flags": found, "parent_math_data_retained": True}


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Build the new standalone corrected native-control packet."""
    root = Path(root)
    require(isinstance(parent, Mapping), "Parent report mapping required")
    _read_parent(root, parent)
    parent_pin_count = _check_parent_pins(root, parent)
    integrity = parent.get("actual_abs_log_integrity", {})
    gate = integrity.get("self_hash_fixture_receipt", {}).get("gate_fixture", {})
    require(integrity.get("inherited_log_cap_certification_verified") is False
            and integrity.get("revalidation_required") is True
            and gate.get("no_revalidation_receipt_admitted_yet") is True,
            "Historical abs-log integrity guard was weakened")

    # Import only the three new providers.  The parent source and old stage
    # modules remain serialized evidence and are never imported here.
    from . import original_core_signed_direction as core
    from . import original_native_directed_common as common
    from . import original_long_reference_directed_controls as long_controls
    from . import original_restoration_directed_controls as restoration_controls

    reviewed_pins = _reviewed_pins(root)
    source_pins = _parent_pin_records(root, parent)
    source_pins[PARENT_PATH] = _pin(root, PARENT_PATH, PARENT_SHA256)
    source_pins.update(reviewed_pins)
    # The coordinator and runner pin themselves dynamically.  Their hashes
    # cannot be constants in the coordinator's review table because changing
    # the table would change the coordinator hash (and the runner is still
    # writing the report).  The runner independently verifies its own bytes.
    source_pins[MODULE] = _pin(root, MODULE)
    source_pins[RUNNER] = _pin(root, RUNNER)
    guards = _source_guards(root, parent)
    premises = _source_premises(parent)

    context = MPIntervalContext()
    context.prec = INTERVAL_PRECISION_BITS
    shared = common.build_inputs(context, parent)
    require(isinstance(shared, Mapping), "Common directed inputs must be a mapping")
    fixture = common.fresh_fixture() if hasattr(common, "fresh_fixture") else None
    require(isinstance(fixture, Mapping) and fixture.get("passed") is True,
            "New closed-form common fixture did not pass")
    require(fixture.get("closed_log_fixtures") == 2
            and fixture.get("independent_pressure_sectors") == 6
            and fixture.get("nonzero_pressure_rows") == 24
            and fixture.get("historical_inputs_used") is False
            and fixture.get("accepted_self_checks_called") is False,
            "New common fixture receipt changed")
    long_packet = long_controls.compute(context, parent, shared)
    restoration_packet = restoration_controls.compute(
        context, parent, shared, long_packet
    )
    _canonicalize_public_packet(long_packet, common)
    _canonicalize_public_packet(restoration_packet, common)
    _validate_component(long_packet, "long_reference")
    _validate_component(restoration_packet, "restoration")
    _validate_ordered_lists(shared, "common")
    shape = _shape_receipt(shared, long_packet, restoration_packet)
    cut_numeric = _cut_numeric_receipt(
        context, shared, long_packet, restoration_packet, core, common
    )
    strict = _strict_controls(long_packet, restoration_packet)
    memory = _memory_receipt(long_packet, restoration_packet)
    physical_geometry = _physical_geometry_receipt(
        parent, long_packet, restoration_packet
    )
    preservation = _preservation_receipt(parent, long_packet, restoration_packet)

    # Providers may contribute additional source pins.  They are merged only
    # after the canonical parent and reviewed self pins have been checked.
    for packet_name, packet in (("long_reference", long_packet),
                                ("restoration", restoration_packet)):
        packet_pins = packet.get("source_pins", packet.get("sourcepins", {}))
        require(isinstance(packet_pins, Mapping),
                packet_name + " source pins are missing")
        for relative, pin in packet_pins.items():
            require(isinstance(pin, Mapping), "Malformed provider source pin")
            path = _safe_relative(str(relative))
            require(pin.get("path") == path
                    and pin.get("passed") is True
                    and pin.get("source_executed") is False,
                    "Noncanonical provider source pin: " + path)
            actual = sha(root / path)
            require(actual == pin.get("expected_sha256") == pin.get("actual_sha256"),
                    "Provider source pin changed: " + path)
            if path in source_pins:
                require(source_pins[path]["actual_sha256"] == actual,
                        "Provider pin conflicts with frozen parent pin: " + path)
            else:
                source_pins[path] = _pin(root, path, actual)

    raw_interval_tree = {
        "shared": _serialize(shared, core, "shared"),
        "long_reference": _serialize(long_packet, core, "long_reference"),
        "restoration": _serialize(restoration_packet, core, "restoration"),
    }
    records = _record_count(raw_interval_tree)
    require(records > 0, "No canonical directed interval records were produced")
    _validate_ordered_lists(raw_interval_tree, "raw_interval_tree")
    json.dumps(raw_interval_tree, allow_nan=False)

    components = {
        "common_inputs": {"passed": True, "build_calls": 1,
                          "fixture_calls": 1, "source_executed": False},
        "long_reference": {"passed": True, "compute_calls": 1,
                            "source_executed": False},
        "restoration": {"passed": True, "compute_calls": 1,
                         "source_executed": False},
    }
    packet = {
        "schema": SCHEMA,
        "passed": True,
        "parent_scope": {
            "report": PARENT_PATH,
            "report_sha256": PARENT_SHA256,
            "report_bytes": PARENT_BYTES,
            "exact_parent_mapping_checked": True,
            "all_parent_pins_checked": True,
            "parent_pin_count": parent_pin_count,
        },
        "source_family": parent["source_family"],
        "source_pins": source_pins,
        "source_guards": guards,
        "source_premises": premises,
        "raw_interval_tree": raw_interval_tree,
        "fixture": dict(fixture),
        "components": components,
        "shape_receipt": shape,
        "cut_numeric_receipt": cut_numeric,
        "strict_control_receipt": strict,
        "memory_receipt": memory,
        "physical_geometry": physical_geometry,
        "preservation": preservation,
        "scope": {
            "full_corrected_native_controls_standalone": True,
            "corrected_native_controls_verified": True,
            "inherited_loop_revalidated": False,
            "current_N_revalidated": False,
            "full_pressure_revalidated": False,
            "global_five_moment_closure": False,
            "full_goal_complete": False,
            "inherited_log_caps_revalidated": False,
            "historical_abs_log_integrity_guard_preserved": True,
        },
        "checks": {
            "passed": True,
            "parent_bytes_checked": True,
            "parent_mapping_checked": True,
            "all_parent_pins_checked": True,
            "source_pins_checked": True,
            "source_guards_checked": True,
            "source_premises_checked": True,
            "context_precision_bits": INTERVAL_PRECISION_BITS,
            "common_built_once": True,
            "long_computed_after_common": True,
            "restoration_computed_after_long": True,
            "raw_interval_tree_serialized": True,
            "canonical_records": records,
            "full_five_row_shapes": True,
            "four_row_p2_shapes": True,
            "strict_slope_exact_tuple_check": True,
            "negative_b_margin_exact_tuple_check": True,
            "source_memory_input_checked": True,
            "physical_geometry_checked": True,
            "numeric_cut_checked": True,
            "old_errors_and_212_retained": True,
            "historical_directed_records": parent["checks"]["directed_records"],
            "accepted_ancestors_replayed": False,
            "source_executed": False,
        },
        "execution_audit": {
            "new_builder": "original_native_directed_controls_v2.build",
            "accepted_builders_or_owners_replayed": False,
            "accepted_checkers_replayed": False,
            "accepted_self_checks_replayed": False,
            "accepted_recurrences_replayed": False,
            "accepted_quadratures_replayed": False,
            "old_stage_calls": False,
            "source_executed": False,
            "common_build_inputs_calls": 1,
            "common_fresh_fixture_calls": 1,
            "long_compute_calls": 1,
            "restoration_compute_calls": 1,
        },
    }
    _validate_execution_flags(packet)
    json.dumps(packet, allow_nan=False)
    return packet


__all__ = [
    "BASE", "PARENT_NAME", "PARENT_PATH", "PARENT_SHA256", "PARENT_BYTES",
    "KEY", "SCHEMA", "MODULE", "RUNNER", "COMMON", "LONG", "RESTORATION",
    "INTERVAL_PRECISION_BITS", "REVIEWED_SOURCE_HASHES", "build", "require", "sha",
]
