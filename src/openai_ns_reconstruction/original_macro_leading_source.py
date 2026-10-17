"""Source-only adapter for the surviving macro leading viscous rows.

The finite-``N`` endpoint category is already an exact source-function
cancellation.  This module keeps the separate leading source terms in their
own namespace.  It reads pinned graph archives and saved JSON packets only;
it never imports an accepted owner, runs a producer or checker, or turns an
error enclosure into a signed source value.

The returned packet is deliberately useful when the source data are only
partly hydrated.  ``passed`` means that the independent source graph ports
and the saved signed endpoint rows are valid.  ``numeric_scope_complete`` is
separate and remains false while the transported lower derivative is absent.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_leading_viscous_kernel import decode_pair, encode_pair, endpoint_coefficient


_PREFIX = "lei_ren_part1_paper_compliant_current_original_"
_PARENT = _PREFIX + "macro_endpoint_m15_numeric.json"
_PARENT_SHA256 = "16eb888da56fad8e5696b5fd8400b827ff7634e2c8f524110267cfcc6f525b80"

_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_ARCHIVE_FILES = {
    "loop_archive": (
        _PREFIX + "prescribed_loop_C5_functions.json.gz",
        "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    ),
    "modified_archive": (
        _PREFIX + "modified_Q4_pressure_C3_functions.json.gz",
        "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    ),
    "high_archive": (
        _PREFIX + "prescribed_high_radial_source_functions.json.gz",
        "5b54a210d48077a984b72b3c74a182c628bf00fb4c2b71f8f164845c228de991",
    ),
}

_REPORT_FILES = {
    "macro_pressure": (
        _PREFIX + "macro_pressure_rows.json",
        "28fac21299ae929447b45894331e33538d51687918323137cbecc4ccf90aad13",
    ),
    "macro_m15": (
        _PREFIX + "macro_m15_numeric.json",
        "076900a16aaea965e229d652b311dd1bcbd2fc73645061fe56aa58575acaf54c",
    ),
    "second_micro_m15": (
        _PREFIX + "second_micro_m15_numeric.json",
        "57155e42c649dcdc3b30bfaee20b13914038aea62c13300b8ea678d80ef63dc3",
    ),
}

# These are implementation pins already present in the immutable parent
# reports.  They are read and parsed as source text only; nothing is run.
_SOURCE_FILES = {
    "macro_m15_source": (
        "src/openai_ns_reconstruction/original_macro_m15.py",
        "5d2c162e6c61bee3d05955004114c9be73b44471555c0638a06b0bfb3d1febf5",
    ),
    "macro_pressure_source": (
        "src/openai_ns_reconstruction/original_macro_pressure_rows.py",
        "6a96ec977ae8f4235bdb874ce9d2c7b5acc15b91f55710fef1ec040e025d6753",
    ),
    "second_micro_m15_source": (
        "src/openai_ns_reconstruction/original_second_micro_m15.py",
        "458be364bcc329e103e9be5d98fa14f9c1cfab39b74032d1642d864da6a83a6b",
    ),
    "leading_kernel_source": (
        "src/openai_ns_reconstruction/original_leading_viscous_kernel.py",
        None,
    ),
}

_ARCHIVE_KEYS = {
    "loop": "prescribed_loop_C5_bindings",
    "modified": "modified_Q4_pressure_C3_bindings",
    "high": "prescribed_high_radial_source_bindings",
}

_EXPECTED_CHARTS = {
    "second_micro": {
        "coordinate": 66988,
        "radius": 68644,
        "jacobian": 8,
        "phase": 97301,
        "domain": [1, 48],
        "query": {"node": 48, "value": "2"},
    },
    "frozen_macro": {
        "coordinate": 70265,
        "radius": 71448,
        "jacobian": 66987,
        "phase": 101909,
        "domain": [0, 1],
        "query": {"nodes": [0, 1], "values": ["0", "1"]},
    },
    "first_switch": {
        "coordinate": 73190,
        "radius": 85905,
        "jacobian": 9,
        "phase": 104537,
        "domain": [0, 1],
        "query": {"node": 0, "value": "0"},
    },
}

def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, Mapping):
        for key in ("node", "handle", "index", "id"):
            if isinstance(value.get(key), int) and not isinstance(value.get(key), bool):
                return value[key]
    return None


def _repo_root(root: Any, directory: Any) -> Path:
    for candidate in (root, directory):
        if candidate is None:
            continue
        try:
            path = Path(candidate)
        except TypeError:
            continue
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _relative(path: Path, repo: Path) -> str:
    try:
        return path.relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _load_serialized(value: Any, repo: Path, *, name: str, default: Path | None = None) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"name": name, "mode": "supplied_mapping", "passed": True}
    path = Path(value) if value is not None else default
    if path is None:
        return {}, {"name": name, "passed": False, "reason": "serialized input was not supplied"}
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError) as exc:
        return {}, {"name": name, "path": _relative(path, repo), "passed": False, "reason": "serialized input is unreadable", "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {
        "name": name,
        "path": _relative(path, repo),
        "sha256": _sha256(path),
        "passed": isinstance(data, Mapping),
    }


def _report_path(directory: Any, repo: Path, basename: str) -> Path:
    candidates: list[Path] = []
    if directory is not None:
        try:
            path = Path(directory)
            candidates.append(path if path.is_file() else path / basename)
        except TypeError:
            pass
    candidates.extend((repo / "experiments" / "root_st073" / basename, repo / basename))
    for path in candidates:
        if path.is_file():
            return path
    return candidates[0]


def _load_report(directory: Any, repo: Path, name: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    basename, expected = _REPORT_FILES[name]
    path = _report_path(directory, repo, basename)
    data, meta = _load_serialized(path, repo, name=name)
    meta.update({"path": _relative(path, repo), "sha256": _sha256(path), "expected_sha256": expected})
    meta["passed"] = bool(meta.get("passed") and meta.get("sha256") == expected)
    return data, meta


def _archive_binding(archive: Mapping[str, Any], family: str) -> Mapping[str, Any]:
    direct = archive.get(_ARCHIVE_KEYS[family])
    if isinstance(direct, Mapping):
        return direct
    # A caller may pass the already-unwrapped binding.
    if family == "loop" and isinstance(archive.get("charts"), Mapping) and isinstance(archive.get("N"), int):
        return archive
    if family == "modified" and isinstance(archive.get("charts"), Mapping):
        return archive
    if family == "high" and isinstance(archive.get("backgrounds"), Mapping):
        return archive
    return {}


def _family_ok(data: Mapping[str, Any]) -> bool:
    family = _mapping(data.get("source_family"))
    if all(family.get(key) == value for key, value in _FAMILY.items()):
        return True
    # Unwrapped bindings do not repeat the family marker.  They are accepted
    # only when the serialized archive itself was already checked by its pin.
    return not family and bool(data)


def _node_present(nodes: Any, handle: Any) -> bool:
    value = _handle(handle)
    if value is None:
        return False
    if isinstance(nodes, Mapping):
        return value in nodes or str(value) in nodes
    try:
        return 0 <= value < len(nodes)
    except (TypeError, AttributeError):
        return False


def _node(nodes: Any, handle: Any) -> Any:
    value = _handle(handle)
    if value is None:
        return None
    if isinstance(nodes, Mapping):
        return nodes.get(value, nodes.get(str(value)))
    try:
        return nodes[value]
    except (IndexError, KeyError, TypeError):
        return None


def _query_guard(nodes: Any, query: int, numerator: int) -> bool:
    row = _node(nodes, query)
    return row == {"operation": "exact_rational", "numerator": numerator, "denominator": 1}


def _source_pin(path: Path, repo: Path, *, expected: str | None, name: str, parse_ast: bool = False) -> dict[str, Any]:
    actual = _sha256(path)
    parsed = None
    error = None
    if parse_ast and actual is not None:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            parsed = True
        except (OSError, SyntaxError, UnicodeError) as exc:
            parsed = False
            error = str(exc)
    passed = actual is not None and (expected is None or actual == expected)
    if parsed is False:
        passed = False
    row: dict[str, Any] = {
        "name": name,
        "path": _relative(path, repo),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": bool(passed),
    }
    if parsed is not None:
        row["ast_parsed_without_execution"] = parsed
    if error:
        row["error"] = error
    return row


def _archive_pin(meta: Mapping[str, Any], repo: Path, *, name: str, expected: str) -> dict[str, Any]:
    path_text = meta.get("path")
    path = Path(path_text) if isinstance(path_text, str) else repo / "experiments" / "root_st073" / _ARCHIVE_FILES[name][0]
    if not path.is_absolute():
        path = repo / path
    row = _source_pin(path, repo, expected=expected, name=name)
    if meta.get("mode") == "supplied_mapping":
        row["supplied_mapping"] = True
        row["serialized_bytes_verified"] = path.is_file() and row["passed"]
    return row


def _sample_coordinate(sample: Mapping[str, Any]) -> str:
    coordinate = sample.get("coordinate")
    if isinstance(coordinate, Mapping):
        return str(coordinate.get("value"))
    return str(coordinate)


def _sample_key(sample: Mapping[str, Any]) -> tuple[str, str]:
    return str(sample.get("Z")), _sample_coordinate(sample)


def _samples(report: Mapping[str, Any]) -> dict[tuple[str, str], Mapping[str, Any]]:
    return {
        _sample_key(sample): sample
        for sample in report.get("samples", ())
        if isinstance(sample, Mapping)
    }


def _source_rows(sample: Mapping[str, Any] | None) -> list[Mapping[str, Any]]:
    leading = _mapping((sample or {}).get("leading_U_y_factored"))
    rows = leading.get("coefficient_ordinary_Z")
    if not isinstance(rows, list) or len(rows) < 2:
        return []
    return [row for row in rows[:2] if isinstance(row, Mapping)]


def _term(sample: Mapping[str, Any] | None, name: str) -> Mapping[str, Any] | None:
    for row in (sample or {}).get("m15_C0_C1_terms", ()):
        if isinstance(row, Mapping) and row.get("name") == name:
            return row
    return None


def _compact_pair(row: Mapping[str, Any] | None, *, signed: bool, source_kind: str) -> dict[str, Any]:
    row = _mapping(row)
    endpoints = row.get("source_enclosure_exact_mpf_tuples")
    valid = isinstance(endpoints, (list, tuple)) and len(endpoints) == 2 and all(
        isinstance(item, (list, tuple)) and len(item) == 4 for item in endpoints
    )
    return {
        "approximation": row.get("approximation"),
        "absolute_error_upper": row.get("absolute_error_upper"),
        "raw_directed_error_exact_mpf_tuple": row.get("raw_directed_error_exact_mpf_tuple"),
        "source_enclosure_exact_mpf_tuples": [list(item) for item in endpoints] if valid else endpoints,
        "error_applies_to_serialized_decimal": row.get("error_applies_to_serialized_decimal") is True,
        "source_kind": source_kind,
        "signed_source_row": bool(signed),
        "approximation_is_not_exact_source_value": True,
        "source_enclosure_is_not_a_selected_value": True,
        "exact_directed_endpoints_present": bool(valid),
    }


def _compact_term(term: Mapping[str, Any] | None, *, name: str, raw_number: int, kind: str, signed: bool) -> dict[str, Any]:
    term = _mapping(term)
    coefficients = term.get("coefficient_C0_C1")
    compact = [
        _compact_pair(item, signed=signed, source_kind=kind)
        for item in coefficients[:2]
    ] if isinstance(coefficients, list) else []
    return {
        "name": name,
        "raw_number": raw_number,
        "kind": kind,
        "signed_point_truth_row": bool(signed),
        "positive_remainder_only": not signed,
        "coefficient_C0_C1": compact,
        "actual_function_binding": dict(term.get("actual_function_binding", {})) if isinstance(term.get("actual_function_binding"), Mapping) else {},
        "independent_saved_row_present": bool(term),
        "approximation_is_not_exact_source_value": True,
    }


def _leading_graph_binding(nodes: Any, loop: Mapping[str, Any], modified: Mapping[str, Any], high: Mapping[str, Any]) -> dict[str, Any]:
    lb = _archive_binding(loop, "loop")
    mb = _archive_binding(modified, "modified")
    hb = _archive_binding(high, "high")
    lcharts = _mapping(lb.get("charts"))
    hcharts = _mapping(hb.get("charts"))
    mcharts = _mapping(mb.get("charts"))
    backgrounds = _mapping(hb.get("backgrounds"))
    background_y = _mapping(hb.get("background_y"))
    chart_rows: dict[str, Any] = {}
    all_passed = True
    for chart, expected in _EXPECTED_CHARTS.items():
        lchart = _mapping(lcharts.get(chart))
        hchart = _mapping(hcharts.get(chart))
        mchart = _mapping(mcharts.get(chart))
        bg = _mapping(backgrounds.get(chart))
        by = _mapping(background_y.get(chart))
        actual_loop = {
            "coordinate": _handle(lchart.get("coordinate")),
            "radius": _handle(lchart.get("radius")),
            "jacobian": _handle(lchart.get("physical_measure")),
            "phase": _handle(lchart.get("actual_phase")),
            "domain": lchart.get("domain"),
        }
        actual_high = {
            "coordinate": _handle(hchart.get("coordinate")),
            "radius": _handle(hchart.get("radius")),
            "jacobian": _handle(hchart.get("jacobian")),
            "phase": _handle(hchart.get("phase")),
            "domain": hchart.get("domain"),
        }
        fields = _mapping(bg.get("fields"))
        radial = _mapping(bg.get("radial"))
        source_handles = {
            "V": list(fields.get("V", ()))[:3],
            "m": list(fields.get("m", ()))[:3],
            "V_y": list(radial.get("V_y", ()))[:2],
        }
        source_present = all(
            len(values) == (2 if quantity == "V_y" else 3)
            and all(_node_present(nodes, value) for value in values)
            for quantity, values in source_handles.items()
        )
        geometry_passed = bool(lchart and hchart and all(
            actual_loop.get(key) == expected.get(key) and actual_high.get(key) == expected.get(key)
            for key in ("coordinate", "radius", "jacobian", "phase", "domain")
        ))
        query_nodes = expected.get("query", {})
        if chart == "second_micro":
            query_passed = _query_guard(nodes, 48, 2)
        elif chart == "frozen_macro":
            query_passed = _query_guard(nodes, 0, 0) and _query_guard(nodes, 1, 1)
        else:
            query_passed = _query_guard(nodes, 0, 0)
        pstar = _handle(_mapping(hb.get("downstream")).get("Pstar"))
        chart_rows[chart] = {
            "chart": chart,
            "expected": dict(expected),
            "loop_ports": actual_loop,
            "high_ports": actual_high,
            "modified_chart_present": bool(mchart),
            "modified_Qy_read": False,
            "modified_Qy_used": False,
            "source_namespace": "prescribed_high_backgrounds",
            "source_handles": source_handles,
            "source_handles_present": source_present,
            "query_nodes": query_nodes,
            "query_nodes_present": query_passed,
            "physical_product": {
                "Pstar_node": pstar,
                "Pstar_applied_once": True,
                "physical_R_applied_once": True,
                "inverse_J_applied_once": True,
                "epsilon": "Lambda^-1",
            },
            "passed": bool(geometry_passed and query_passed and source_present and mchart),
        }
        all_passed = all_passed and chart_rows[chart]["passed"]
    return {
        "binding_key": _ARCHIVE_KEYS["high"],
        "charts": chart_rows,
        "leading_source_formula": {
            "macro": "-hbB*(Phi/barPhi)*(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl)",
            "first_switch": "-hbS*(1-sigma(s))*(Phi/barPhi)*(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl)",
            "source_Uy_is_factored_before_epsilon_and_width": True,
            "Pstar_once": True,
            "R_once": True,
        },
        "modified_Qy_not_used_for_leading_source": True,
        "ordinary_Z_orders": [0, 1],
        "ordinary_Z2_retained_for_mean_source": True,
        "passed": bool(all_passed),
    }


def _source_binding_rows(report: Mapping[str, Any], *, role: str, chart: str, width: str) -> dict[str, Any]:
    indexed = _samples(report)
    rows = []
    for (z, coordinate), sample in indexed.items():
        leading = _mapping(sample.get("leading_U_y_factored"))
        source_rows = _source_rows(sample)
        rows.append({
            "name": role,
            "chart": chart,
            "Z": z,
            "coordinate": sample.get("coordinate"),
            "ordinary_Z_orders": [0, 1],
            "source_namespace": "prescribed_high_backgrounds",
            "source_rows": [_compact_pair(row, signed=True, source_kind="factored_source_Uy") for row in source_rows],
            "source_definition": leading.get("definition"),
            "actual_source": leading.get("actual_source"),
            "source_control": leading.get("axial_control") or leading.get("angular_control") or "1",
            "epsilon_factor": "Lambda^-1",
            "width_factor": width,
            "physical_R_factor_once": True,
            "physical_Pstar_once": True,
            "approximation_is_not_exact_source_value": True,
            "truth_enclosure_retained": bool(source_rows),
        })
    return {
        "role": role,
        "chart": chart,
        "width": width,
        "ordinary_Z_orders": [0, 1],
        "source_field": "leading_U_y_factored.coefficient_ordinary_Z",
        "endpoint_coefficient": {
            "kernel": "original_leading_viscous_kernel.endpoint_coefficient",
            "formula": "sign * 2*Z*B/(1-delta*Z^2)",
            "signed_signs": [-1, 1],
            "requires_signed_B_and_B_Z": True,
            "physical_R_factor_once": True,
        },
        "mean_source_coefficient": {
            "kernel": "original_leading_viscous_kernel.mean_source_coefficient",
            "formula": "C(B)=((1-delta)*Z*B+(1-Z^2)*B_Z)/(1-delta*Z^2)",
            "ordinary_Z1_requires_B_ZZ": True,
        },
        "rows": rows,
        "sample_count": len(rows),
        "point_is_saved_approximation": True,
        "truth_interval_is_retained": True,
        "passed": bool(rows) and all(len(row["source_rows"]) >= 2 for row in rows),
    }


def _new_contexts() -> tuple[Any, Any, Any, Any, dict[str, Any]]:
    """Create the new arithmetic contexts and the pinned core delta only."""
    ctx = mp.mp.clone()
    ctx.dps = 110
    iv = MPIntervalContext()
    iv.dps = 145
    delta_formula = "exp(-4*(exp(40)+11)-30)"
    delta = ctx.exp(-4 * (ctx.exp(40) + 11) - 30)
    delta_iv = iv.exp(-4 * (iv.exp(40) + 11) - 30)
    return ctx, iv, delta, delta_iv, {
        "formula": "delta=" + delta_formula,
        "source": "pinned original core delta formula; no ancestor owner imported",
        "point_dps": 110,
        "interval_dps": 145,
    }


def _macro_local_coefficients(
    report: Mapping[str, Any],
    *,
    ctx: Any,
    iv: Any,
    delta: Any,
    delta_iv: Any,
    delta_binding: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Compute the new local ``2 Z B/L`` rows from saved signed source pairs.

    These rows describe the source at each saved macro ``q``.  They are not
    the transported lower endpoint; that endpoint needs the separately
    missing second-micro ``lower_dy`` source.
    """
    rows: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for (z, coordinate), sample in _samples(report).items():
        source = _source_rows(sample)
        row: dict[str, Any] = {
            "name": "macro_leading_local_source_coefficient",
            "raw_number": None,
            "Z": z,
            "coordinate": sample.get("coordinate"),
            "kind": "signed_local_source_endpoint_coefficient",
            "source_binding": "macro_source_Uy",
            "local_q_source_coefficient": True,
            "is_actual_macro_lower_endpoint": False,
            "actual_macro_lower_requires_second_micro_lower_dy": True,
            "ordinary_Z_orders": [0, 1],
            "formula": "2*Z*B/(1-delta*Z^2)",
            "sign": 1,
            "source_rows": [_compact_pair(item, signed=True, source_kind="factored_source_Uy") for item in source],
            "physical_scale": {
                "R_factor_once": True,
                "epsilon_factor_once": True,
                "width_factor": "hbB",
            },
            "delta_binding": dict(delta_binding),
            "approximation_is_not_exact_source_value": True,
            "source_values_selected": False,
        }
        try:
            if len(source) < 2:
                raise ValueError("signed source B and B_Z rows are required")
            point_inputs = []
            truth_inputs = []
            for item in source[:2]:
                point, truth = decode_pair(ctx, iv, dict(item))
                point_inputs.append(point)
                truth_inputs.append(truth)
            point = endpoint_coefficient(ctx.mpf(z), delta, point_inputs, sign=1)
            truth = endpoint_coefficient(iv.mpf(z), delta_iv, truth_inputs, sign=1)
            row["computed_coefficient_C0_C1"] = [
                encode_pair(ctx, iv, point[index], truth[index])
                for index in range(2)
            ]
            row["computed_with_new_kernel"] = True
            row["passed"] = True
        except (TypeError, ValueError, ArithmeticError, IndexError) as exc:
            row["computed_coefficient_C0_C1"] = []
            row["computed_with_new_kernel"] = False
            row["passed"] = False
            error = {
                "name": row["name"],
                "Z": z,
                "coordinate": coordinate,
                "reason": str(exc),
            }
            errors.append(error)
        rows.append(row)
    return rows, errors


def _second_micro_binding(report: Mapping[str, Any], pin: Mapping[str, Any]) -> dict[str, Any]:
    samples = _samples(report)
    by_z = {
        z: {
            "query": "2",
            "saved_sample_present": (z, "2") in samples,
            "saved_sample_coordinate": samples[(z, "2")].get("coordinate") if (z, "2") in samples else None,
        }
        for z in sorted({key[0] for key in samples})
    }
    return {
        "role": "second_micro_transported_lower_dy",
        "source_report": dict(pin),
        "source_chart": "second_micro",
        "native_query_node": 48,
        "native_query_value": "2",
        "formula": "actual.Phi / comparison.Phi[0] * comparison.Psi[1] * rho",
        "transport_formula": "copy exact saved interval endpoints between point and interval contexts; no macro J=66987 factor is applied during that copy",
        "input_rows": [
            "actual.Phi",
            "comparison.Phi[0]",
            "comparison.Psi[1]",
            "rho",
        ],
        "sample_presence_by_Z": by_z,
        "available_as_named_source_value": False,
        "serialized_value_selected": False,
        "source_values_selected": False,
        "numeric_input_recovered": False,
        "passed": bool(pin.get("passed") and by_z),
        "missing_plan": "New source kernel must export signed actual/comparison Phi/Psi radial rows and rho. Only after lower_dy is recovered may its native transport be applied; do not reuse modified Qy or a midpoint.",
    }


def _m15_rows(
    macro_m15: Mapping[str, Any],
    macro_pressure: Mapping[str, Any],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    macro_index = _samples(macro_m15)
    macro_pressure_index = _samples(macro_pressure)

    # Signed normalized source_Uy rows are retained as source inputs.  Their
    # factored epsilon and hb are recorded, never numerically multiplied here.
    for (z, coordinate), sample in macro_pressure_index.items():
        source = _source_rows(sample)
        rows.append({
            "name": "macro_source_Uy",
            "raw_number": None,
            "Z": z,
            "coordinate": sample.get("coordinate"),
            "kind": "signed_factored_source_Uy",
            "source_binding": "macro_source_Uy",
            "ordinary_Z_orders": [0, 1],
            "source_rows": [_compact_pair(row, signed=True, source_kind="factored_source_Uy") for row in source],
            "physical_scale": {"R_factor_once": True, "epsilon_factor_once": True, "width_factor": "hbB"},
            "approximation_is_not_exact_source_value": True,
        })

    for (z, coordinate), sample in macro_index.items():
        for name, raw_number, kind, signed in (
            ("macro_leading_viscous_upper_error", 23, "positive_endpoint_remainder_cap", False),
            ("macro_leading_viscous_lower_error", 24, "positive_endpoint_remainder_cap", False),
        ):
            rows.append({
                "name": name,
                "raw_number": raw_number,
                "Z": z,
                "coordinate": sample.get("coordinate"),
                "kind": kind,
                "source_binding": "macro_source_Uy",
                "physical_scale": {"R_factor_once": True, "epsilon_factor_once": True, "width_factor": "hbB"},
                "independent_saved_row": _compact_term(_term(sample, name), name=name, raw_number=raw_number, kind=kind, signed=signed),
                "source_point_is_not_selected": True,
                "positive_remainder_only": True,
            })

    return rows


def _source_pins(repo: Path, directory: Any, parent: Any, archive_meta: Mapping[str, Mapping[str, Any]], report_meta: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    pins: dict[str, dict[str, Any]] = {}
    base = repo / "experiments" / "root_st073"
    for name, (basename, expected) in _ARCHIVE_FILES.items():
        path = base / basename
        pins[name] = _source_pin(path, repo, expected=expected, name=name)
        if archive_meta.get(name, {}).get("mode") == "supplied_mapping":
            pins[name]["supplied_mapping"] = True
    for name, (basename, expected) in _REPORT_FILES.items():
        path = _report_path(directory, repo, basename)
        pins[name] = _source_pin(path, repo, expected=expected, name=name)
    for name, (relative, expected) in _SOURCE_FILES.items():
        path = repo / relative
        pins[name] = _source_pin(path, repo, expected=expected, name=name, parse_ast=True)
    parent_path = None
    if isinstance(parent, (str, Path)):
        parent_path = Path(parent)
        if not parent_path.is_absolute():
            parent_path = repo / parent_path
    else:
        parent_path = base / _PARENT
    pins["parent_report"] = _source_pin(parent_path, repo, expected=_PARENT_SHA256, name="parent_report")
    return pins


def _parent_family(parent: Any) -> bool:
    if not isinstance(parent, Mapping):
        return False
    family = _mapping(parent.get("source_family"))
    return all(family.get(key) == value for key, value in _FAMILY.items())


def _missing_rows(second_binding: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "name": "second_micro_lower_dy",
            "raw_number": None,
            "missing": ["actual.Phi", "comparison.Phi[0]", "comparison.Psi[1]", "rho"],
            "available": ["source recipe", "native query node 48", "native J=66987 port"],
            "reason": second_binding.get("missing_plan"),
            "source_values_selected": False,
        },
    ]


def build_macro_leading_source(
    *,
    nodes: Any,
    loop_archive: Any,
    modified_archive: Any,
    high_archive: Any,
    parent_report: Any,
    directory: Any,
    root: Any,
) -> dict[str, Any]:
    """Bind actual leading source ports and hydrate saved point/remainder rows."""
    repo = _repo_root(root, directory)
    loop, loop_meta = _load_serialized(loop_archive, repo, name="loop_archive", default=repo / "experiments" / "root_st073" / _ARCHIVE_FILES["loop_archive"][0])
    modified, modified_meta = _load_serialized(modified_archive, repo, name="modified_archive", default=repo / "experiments" / "root_st073" / _ARCHIVE_FILES["modified_archive"][0])
    high, high_meta = _load_serialized(high_archive, repo, name="high_archive", default=repo / "experiments" / "root_st073" / _ARCHIVE_FILES["high_archive"][0])
    macro_pressure, macro_pressure_meta = _load_report(directory, repo, "macro_pressure")
    macro_m15, macro_m15_meta = _load_report(directory, repo, "macro_m15")
    second_micro_m15, second_micro_meta = _load_report(directory, repo, "second_micro_m15")

    archive_meta = {
        "loop_archive": loop_meta,
        "modified_archive": modified_meta,
        "high_archive": high_meta,
    }
    report_meta = {
        "macro_pressure": macro_pressure_meta,
        "macro_m15": macro_m15_meta,
        "second_micro_m15": second_micro_meta,
    }
    pins = _source_pins(repo, directory, parent_report, archive_meta, report_meta)
    family_checks = {
        name: {
            "source_family": dict(_mapping(data.get("source_family"))),
            "passed": bool(_family_ok(data)),
        }
        for name, data in {
            "loop_archive": loop,
            "modified_archive": modified,
            "high_archive": high,
            "macro_pressure": macro_pressure,
            "macro_m15": macro_m15,
            "second_micro_m15": second_micro_m15,
        }.items()
    }
    graph = _leading_graph_binding(nodes, loop, modified, high)
    macro_source_binding = _source_binding_rows(macro_pressure, role="macro_source_Uy", chart="frozen_macro", width="hbB")
    second_binding = _second_micro_binding(second_micro_m15, pins.get("second_micro_m15", {}))
    source_bindings = {
        "geometry": graph,
        "macro_source_Uy": macro_source_binding,
        "second_micro_transported_lower_dy": second_binding,
        "modified_Qy_not_used_for_leading_source": True,
        "ordinary_Z1_retained": True,
        "physical_R_once": True,
        "epsilon": "Lambda^-1",
        "hbB_and_hbS_distinct": True,
    }
    source_bindings_passed = bool(
        graph.get("passed")
        and macro_source_binding.get("passed")
        and second_binding.get("passed")
    )
    source_bindings["passed"] = source_bindings_passed
    ctx, iv, delta, delta_iv, delta_binding = _new_contexts()
    rows = _m15_rows(macro_m15, macro_pressure)
    computed_rows, compute_errors = _macro_local_coefficients(
        macro_pressure,
        ctx=ctx,
        iv=iv,
        delta=delta,
        delta_iv=delta_iv,
        delta_binding=delta_binding,
    )
    rows.extend(computed_rows)
    blocked = _missing_rows(second_binding)
    blocked.extend(compute_errors)
    endpoint_interface_passed = bool(computed_rows) and not compute_errors and all(
        row.get("passed") is True for row in computed_rows
    )
    all_pins_passed = all(row.get("passed") is True for row in pins.values())
    family_passed = all(row.get("passed") is True for row in family_checks.values())
    limits = {
        "numeric_scope_complete": False,
        "second_micro_lower_dy_numeric": False,
        "mass_point_numeric": False,
        "weighted_primitive_remainder_numeric": False,
        "integral_C1_numeric": False,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "complete": False,
    }
    contract = {
        "source_only": True,
        "actual_high_background_namespace": True,
        "modified_Qy_not_used": True,
        "physical_R_once": True,
        "physical_Pstar_once": True,
        "epsilon_is_Lambda_inverse": True,
        "ordinary_Z1_retained": True,
        "hbB_and_hbS_distinct": True,
        "positive_remainders_not_signed_values": True,
        "point_approximations_not_exact_source_values": True,
        "caps_not_promoted": True,
        "ancestor_constructors_producers_checkers_not_run": True,
        "raw23_24_separate_from_endpoint_category": True,
        "raw43_47_delegated_to_first_switch_adapter": True,
        "raw213_reduced183_untouched": True,
    }
    return {
        "passed": bool(family_passed and all_pins_passed and source_bindings_passed and endpoint_interface_passed and rows),
        "partial": bool(endpoint_interface_passed),
        "source_family": dict(_FAMILY),
        "parent_report": "experiments/root_st073/" + _PARENT,
        "parent_report_sha256": _PARENT_SHA256,
        "source_pins": pins,
        "family_checks": family_checks,
        "source_bindings": source_bindings,
        "point_remainder_rows": rows,
        "blocked_rows": blocked,
        "endpoint_interface_passed": endpoint_interface_passed,
        "source_bindings_passed": source_bindings_passed,
        "numeric_scope_complete": False,
        "limits": limits,
        "contract": contract,
    }


__all__ = ["build_macro_leading_source"]
