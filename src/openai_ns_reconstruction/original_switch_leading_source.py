"""Source-only adapter for the first-switch leading viscous rows.

The finite-``N`` endpoint category is already an exact source-function
cancellation.  This module handles the separate leading source terms 43--47
from their saved pressure and m15 packets.  It deliberately does not import
an owner, construct a numerical producer, or turn an enclosure into a
function value.  Endpoint coefficients are recomputed from the saved signed
``U_y`` source jets and are compared with the independent m15 rows.  The
mean-source rows remain blocked until the source supplies the missing second
ordinary derivative and weighted primitive remainder; the exact phase-one
pulse mass branch is bound separately without executing the ancestor.
"""

from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_leading_viscous_kernel import (
    decode_pair,
    encode_pair,
    endpoint_coefficient,
    mean_source_value,
)


_PREFIX = "lei_ren_part1_paper_compliant_current_original_"
_REPORTS = {
    "first_switch_pressure": _PREFIX + "first_switch_pressure_rows.json",
    "first_switch_m15": _PREFIX + "first_switch_m15_numeric.json",
    "macro_pressure": _PREFIX + "macro_pressure_rows.json",
    "macro_m15": _PREFIX + "macro_m15_numeric.json",
    "second_micro_m15": _PREFIX + "second_micro_m15_numeric.json",
}
_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
_ARCHIVE_KEYS = {
    "loop": "prescribed_loop_C5_bindings",
    "modified": "modified_Q4_pressure_C3_bindings",
    "high": "prescribed_high_radial_source_bindings",
}
_EXPECTED_CHARTS = {
    "second_micro": {"coordinate": 66988, "radius": 68644, "jacobian": 8, "phase": 97301, "domain": [1, 48]},
    "frozen_macro": {"coordinate": 70265, "radius": 71448, "jacobian": 66987, "phase": 101909, "domain": [0, 1]},
    "first_switch": {"coordinate": 73190, "radius": 85905, "jacobian": 9, "phase": 104537, "domain": [0, 1]},
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _repo_root(root: Any, directory: Any) -> Path:
    for candidate in (root, directory):
        if candidate is None:
            continue
        path = Path(candidate)
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


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()


def _load_serialized(value: Any, repo: Path, *, name: str, default: Path | None = None) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"name": name, "mode": "supplied_mapping", "passed": True}
    path = Path(value) if value is not None else default
    if path is None:
        return {}, {"name": name, "passed": False, "reason": "no serialized input supplied"}
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
        path = Path(directory)
        candidates.append(path if path.is_file() else path / basename)
    candidates.extend((repo / "experiments" / "root_st073" / basename, repo / basename))
    for path in candidates:
        if path.is_file():
            return path
    return candidates[0]


def _load_report(directory: Any, repo: Path, name: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = _report_path(directory, repo, _REPORTS[name])
    data, meta = _load_serialized(path, repo, name=name)
    meta["path"] = _relative(path, repo)
    meta["sha256"] = _sha256(path)
    return data, meta


def _family_ok(data: Mapping[str, Any]) -> bool:
    family = _mapping(data.get("source_family"))
    return all(family.get(key) == value for key, value in _FAMILY.items())


def _archive_binding(archive: Mapping[str, Any], family: str) -> Mapping[str, Any]:
    return _mapping(archive.get(_ARCHIVE_KEYS[family]))


def _handle_presence(nodes: Any, handle: Any) -> bool:
    if not isinstance(handle, int):
        return False
    if isinstance(nodes, Mapping):
        return handle in nodes or str(handle) in nodes
    try:
        return 0 <= handle < len(nodes)
    except TypeError:
        return False


def _source_handle_guard(nodes: Any, archive: Mapping[str, Any], family: str) -> dict[str, Any]:
    binding = _archive_binding(archive, family)
    charts = _mapping(binding.get("charts"))
    rows: dict[str, Any] = {}
    for chart, expected in _EXPECTED_CHARTS.items():
        item = _mapping(charts.get(chart))
        if family == "loop":
            # The loop receipt names the physical Jacobian ``physical_measure``
            # and stores the global phase as ``actual_phase``.  The high-radial
            # receipt uses the shorter ``jacobian``/``phase`` names.
            actual = {
                "coordinate": item.get("coordinate"),
                "radius": item.get("radius"),
                "jacobian": item.get("physical_measure"),
                "phase": item.get("actual_phase"),
                "domain": item.get("domain"),
            }
        else:
            actual = {key: item.get(key) for key in expected}
        row = {
            "chart": chart,
            "expected": expected,
            "actual": actual,
            "chart_present": bool(item),
            "passed": bool(item) and all(actual.get(key) == value for key, value in expected.items()),
        }
        handles = {key: actual.get(key) for key in ("coordinate", "radius", "jacobian", "phase") if actual.get(key) is not None}
        row["handles_present_in_nodes"] = all(_handle_presence(nodes, value) for value in handles.values()) if handles else False
        rows[chart] = row
    return {
        "family": family,
        "binding_key": _ARCHIVE_KEYS[family],
        "rows": rows,
        "passed": bool(rows) and all(
            row["passed"] and row["handles_present_in_nodes"]
            for row in rows.values()
        ),
    }


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


def _source_rows(sample: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    leading = _mapping(sample.get("leading_U_y_factored"))
    rows = leading.get("coefficient_ordinary_Z")
    if not isinstance(rows, list) or len(rows) < 2 or not all(isinstance(row, Mapping) for row in rows[:2]):
        return []
    return [copy.deepcopy(row) for row in rows[:2]]


def _m15_term(sample: Mapping[str, Any], name: str) -> Mapping[str, Any] | None:
    for term in sample.get("m15_C0_C1_terms", ()):
        if isinstance(term, Mapping) and term.get("name") == name:
            return term
    return None


def _delta_source_guard(repo: Path) -> dict[str, Any]:
    """Pin and AST-check the actual core definition of ``delta``."""
    relative = "src/openai_ns_reconstruction/original_core_numeric.py"
    path = repo / relative
    digest = _sha256(path)
    required = {
        ("self.logP", "c.exp(40)+11"),
        ("self.delta", "c.exp(-4*self.logP-30)"),
    }
    found: set[tuple[str, str]] = set()
    error = None
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            expression = ast.unparse(node.value).replace(" ", "")
            for target in node.targets:
                target_text = ast.unparse(target).replace(" ", "")
                pair = (target_text, expression)
                if pair in required:
                    found.add(pair)
    except (OSError, SyntaxError, UnicodeError) as exc:
        error = str(exc)
    return {
        "path": relative,
        "sha256": digest,
        "required_ast": sorted([list(pair) for pair in required]),
        "found_ast": sorted([list(pair) for pair in found]),
        "passed": bool(digest and error is None and found == required),
        "error": error,
    }


def _mass_source_guard(repo: Path) -> dict[str, Any]:
    """Pin the exact phase-zero/phase-one mass branches without executing them."""
    relative = "experiments/root_st073/lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
    path = repo / relative
    digest = _sha256(path)
    expected = {
        "pulse_control_masses": {
            "phase_zero": ("hi==0", "return(c.mpf(0),c.mpf(0))"),
            "phase_one": ("lo==hi==1", "return(c.mpf('.5'),c.mpf('.5'))"),
        },
        "CompliantActualBridgeIntegrals.pulse_masses": {
            "phase_zero": ("hi==0", "return(c.mpf(0),c.mpf(0))"),
            "phase_one": ("lo==hi==1", "return(c.mpf('.5'),c.mpf('.5'))"),
        },
    }
    found: dict[str, dict[str, bool]] = {
        name: {branch: False for branch in rows}
        for name, rows in expected.items()
    }
    error = None

    def branch_present(function: ast.FunctionDef, test_text: str, return_text: str) -> bool:
        for statement in function.body:
            if not isinstance(statement, ast.If):
                continue
            actual_test = ast.unparse(statement.test).replace(" ", "")
            if actual_test != test_text:
                continue
            for child in statement.body:
                if isinstance(child, ast.Return) and child.value is not None:
                    actual_return = ast.unparse(child).replace(" ", "")
                    if actual_return == return_text:
                        return True
        return False

    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        functions: dict[str, ast.FunctionDef] = {}
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                functions[node.name] = node
            elif isinstance(node, ast.ClassDef):
                for child in node.body:
                    if isinstance(child, ast.FunctionDef):
                        functions[f"{node.name}.{child.name}"] = child
        for name, branches in expected.items():
            function = functions.get(name)
            if function is None:
                continue
            for branch, (test_text, return_text) in branches.items():
                found[name][branch] = branch_present(function, test_text, return_text)
    except (OSError, SyntaxError, UnicodeError) as exc:
        error = str(exc)
    return {
        "path": relative,
        "sha256": digest,
        "required_ast": {
            name: {branch: [test_text, return_text] for branch, (test_text, return_text) in branches.items()}
            for name, branches in expected.items()
        },
        "found_ast": found,
        "phase_zero_exact_zero": all(rows["phase_zero"] for rows in found.values()),
        "phase_one_exact_half": all(rows["phase_one"] for rows in found.values()),
        "passed": bool(
            digest
            and error is None
            and all(all(branches.values()) for branches in found.values())
        ),
        "error": error,
    }


def _contexts(repo: Path) -> tuple[Any, Any, Any, Any, dict[str, Any]]:
    ctx = mp.mp.clone()
    ctx.dps = 110
    iv = MPIntervalContext()
    iv.dps = 145
    source_guard = _delta_source_guard(repo)
    delta_binding = {
        "formula": "delta=exp(-4*(exp(40)+11)-30)",
        "source": "src/openai_ns_reconstruction/original_core_numeric.py:144-148",
        "source_pin": source_guard,
        "source_constant_defined_analytically": True,
        "numerical_delta_has_directed_enclosure": True,
        "source_value_selected": False,
    }
    delta = ctx.exp(-4 * (ctx.exp(40) + 11) - 30)
    delta_iv = iv.exp(-4 * (iv.exp(40) + 11) - 30)
    return ctx, iv, delta, delta_iv, delta_binding


def _interval_or_none(value: Any) -> Any | None:
    try:
        return mp.mpf(tuple(value))
    except (TypeError, ValueError):
        return None


def _source_pair(ctx: Any, iv: Any, row: Mapping[str, Any]) -> tuple[Any, Any]:
    return decode_pair(ctx, iv, dict(row))


def _source_summary(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "approximation": row.get("approximation"),
        "absolute_error_upper": row.get("absolute_error_upper"),
        "source_enclosure_exact_mpf_tuples": copy.deepcopy(row.get("source_enclosure_exact_mpf_tuples")),
        "approximation_is_not_exact_source_value": True,
    }


def _computed_endpoint(ctx: Any, iv: Any, z: str, delta: Any, delta_iv: Any, source_rows: Sequence[Mapping[str, Any]], sign: int) -> dict[str, Any]:
    point_pairs = [_source_pair(ctx, iv, row)[0] for row in source_rows[:2]]
    truth_pairs = [_source_pair(ctx, iv, row)[1] for row in source_rows[:2]]
    point = endpoint_coefficient(ctx.mpf(z), delta, point_pairs, sign=sign)
    truth = endpoint_coefficient(iv.mpf(z), delta_iv, truth_pairs, sign=sign)
    return {
        "point": [encode_pair(ctx, iv, point[index], truth[index]) for index in range(2)],
        "truth": [encode_pair(ctx, iv, point[index], truth[index]) for index in range(2)],
    }


def _numeric_agreement(ctx: Any, computed: Mapping[str, Any], raw: Mapping[str, Any]) -> dict[str, Any]:
    """Compare signed point decimals without using an enclosure cap.

    The saved ``absolute_error_upper`` is a source enclosure radius.  It is
    intentionally excluded from this audit: a large cap cannot certify a
    recomputed signed coefficient.  The comparison uses the 85-digit saved
    decimal precision with a fixed relative tolerance.
    """
    tolerance = ctx.mpf("1e-70")
    try:
        cp = ctx.mpf(str(computed["approximation"]))
        rp = ctx.mpf(str(raw["approximation"]))
    except (ValueError, TypeError):
        return {"passed": False, "error": "serialized point is not numeric"}
    difference = abs(cp - rp)
    scale = max(abs(cp), abs(rp))
    if scale == 0:
        relative = ctx.mpf(0)
        passed = difference == 0
    else:
        relative = difference / scale
        passed = relative <= tolerance
    return {
        "passed": bool(passed),
        "absolute_difference": ctx.nstr(difference, 50),
        "relative_difference": ctx.nstr(relative, 50),
        "relative_tolerance": "1e-70",
        "zero_scale": bool(scale == 0),
    }


def _interval_bounds(value: Any) -> tuple[Any, Any] | None:
    try:
        rows = value._mpi_
        return mp.mpf(tuple(rows[0])), mp.mpf(tuple(rows[1]))
    except (AttributeError, TypeError, ValueError):
        return None


def _truth_overlaps(ctx: Any, computed: Mapping[str, Any], raw: Mapping[str, Any]) -> bool:
    try:
        crows = computed["source_enclosure_exact_mpf_tuples"]
        rrows = raw["source_enclosure_exact_mpf_tuples"]
        cends = [ctx.mpf(tuple(row)) for row in crows]
        rends = [ctx.mpf(tuple(row)) for row in rrows]
        cl, cu = min(cends), max(cends)
        rl, ru = min(rends), max(rends)
        return cl <= ru and rl <= cu
    except (KeyError, TypeError, ValueError):
        return False


def _compare_endpoint(ctx: Any, iv: Any, computed: Mapping[str, Any], raw_rows: Sequence[Mapping[str, Any]], *, empty: bool) -> dict[str, Any]:
    if empty:
        return {
            "passed": True,
            "status": "not_compared_empty_increment",
            "reason": "The accepted m15 producer serializes the coordinate 0 increment as zero; this is not a function-zero proof.",
        }
    point = [_numeric_agreement(ctx, computed["point"][index], raw_rows[index]) for index in range(2)]
    truth = [_truth_overlaps(ctx, computed["truth"][index], raw_rows[index]) for index in range(2)]
    return {
        "passed": bool(all(row.get("passed") is True for row in point)),
        "status": "compared_independent_saved_row",
        "point_agreement": point,
        "truth_intervals_overlap_informational": truth,
        "truth_overlap_is_not_a_pass_gate": True,
    }


def _compare_scalar(ctx: Any, computed: Mapping[str, Any], raw: Mapping[str, Any]) -> dict[str, Any]:
    """Compare one signed C0 point while retaining truth overlap as information."""
    point = _numeric_agreement(ctx, computed, raw)
    return {
        "passed": point.get("passed") is True,
        "status": "compared_independent_saved_row",
        "point_agreement": point,
        "truth_intervals_overlap_informational": _truth_overlaps(ctx, computed, raw),
        "truth_overlap_is_not_a_pass_gate": True,
    }


def _exact_source_pair(ctx: Any, iv: Any, value: Any, *, reason: str) -> dict[str, Any]:
    """Serialize a value fixed by a source branch, with its exactness explicit."""
    pair = encode_pair(ctx, iv, ctx.mpf(value), iv.mpf(value))
    pair.update({
        "approximation_is_exact_source_value": True,
        "source_enclosure_is_not_a_selected_value": False,
        "exact_source_branch": reason,
    })
    return pair


def _geometry_bindings(nodes: Any, loop: Mapping[str, Any], modified: Mapping[str, Any], high: Mapping[str, Any]) -> dict[str, Any]:
    loop_guard = _source_handle_guard(nodes, loop, "loop")
    high_guard = _source_handle_guard(nodes, high, "high")
    modified_binding = _archive_binding(modified, "modified")
    modified_charts = _mapping(modified_binding.get("charts"))
    modified_rows = {
        chart: {"chart_present": chart in modified_charts, "keys": sorted(_mapping(modified_charts.get(chart)).keys())}
        for chart in _EXPECTED_CHARTS
    }
    return {
        "loop": loop_guard,
        "high": high_guard,
        "modified": {"binding_key": _ARCHIVE_KEYS["modified"], "rows": modified_rows, "passed": bool(modified_rows) and all(row["chart_present"] for row in modified_rows.values())},
        "physical_geometry": {
            "second_micro": {"coordinate": 66988, "radius": 68644, "J": 8, "phase": 97301},
            "frozen_macro": {"coordinate": 70265, "radius": 71448, "J": 66987, "phase": 101909},
            "first_switch": {"coordinate": 73190, "radius": 85905, "J": 9, "phase": 104537},
            "physical_measure": "dy=hbS ds on first_switch; physical R appears once",
        },
        "passed": bool(loop_guard["passed"] and high_guard["passed"] and modified_rows and all(row["chart_present"] for row in modified_rows.values())),
    }


def _source_bindings(
    *,
    nodes: Any,
    loop: Mapping[str, Any],
    modified: Mapping[str, Any],
    high: Mapping[str, Any],
    first_pressure: Mapping[str, Any],
    macro_pressure: Mapping[str, Any],
    second_micro_m15: Mapping[str, Any],
    pins: Mapping[str, Any],
    delta_binding: Mapping[str, Any],
    mass_source: Mapping[str, Any],
) -> dict[str, Any]:
    geometry = _geometry_bindings(nodes, loop, modified, high)
    first_samples = _samples(first_pressure)
    macro_samples = _samples(macro_pressure)
    micro_samples = _samples(second_micro_m15)
    macro_endpoint_rows = {
        z: _source_rows(sample)
        for (z, coordinate), sample in macro_samples.items()
        if coordinate == "1"
    }
    first_sample_count = len(first_samples)
    micro_queries = {
        z: {
            "query": "2",
            "saved_sample_present": (z, "2") in micro_samples,
            "saved_sample_coordinate": _mapping(micro_samples.get((z, "2")).get("coordinate")) if (z, "2") in micro_samples else None,
        }
        for z in sorted({key[0] for key in first_samples})
    }
    return {
        "geometry": geometry,
        "source_family": dict(_FAMILY),
        "macro_source_Uy": {
            "source_report": pins.get("macro_pressure"),
            "query": "frozen_macro coordinate 1",
            "ordinary_Z_orders": [0, 1],
            "rows_by_Z": {z: [_source_summary(row) for row in rows] for z, rows in macro_endpoint_rows.items()},
            "point_is_legacy_midpoint_approximation": True,
            "truth_interval_is_retained": True,
        },
        "first_switch_upper_source_Uy": {
            "source_report": pins.get("first_switch_pressure"),
            "ordinary_Z_orders": [0, 1],
            "sample_count": first_sample_count,
            "source_field": "leading_U_y_factored.coefficient_ordinary_Z",
            "point_is_legacy_midpoint_approximation": True,
            "truth_interval_is_retained": True,
            "exact_s1_control_is_separate": True,
        },
        "first_switch_lower_source_Uy": {
            "source_report": pins.get("macro_pressure"),
            "source_field": "macro source_Uy at frozen_macro coordinate 1",
            "alias_to_modified_Qy": False,
            "point_is_legacy_midpoint_approximation": True,
            "truth_interval_is_retained": True,
        },
        "second_micro_transported_lower_dy": {
            "source_report": pins.get("second_micro_m15"),
            "source_chart": "second_micro",
            "query": "2",
            "sample_presence_by_Z": micro_queries,
            "formula": "lower.actual.Phi/lower.comparison.Phi[0] * lower.comparison.Psi[1] * lower.rho",
            "transfer_formula": "ordinary first_switch J interval transfer before source combination",
            "serialized_value_selected": False,
            "available_as_named_source_value": False,
            "reason": "The immutable m15 packet does not serialize lower_dy; the source recipe is bound for a later bounded kernel.",
        },
        "mass": {
            "source": "original_first_switch_m15.source_norms",
            "source_file": "src/openai_ns_reconstruction/original_first_switch_m15.py:100",
            "formula": "mass, _ = bridge.pulse_masses(phase); mass = transfer(mass)",
            "serialized_value_selected": False,
            "available_as_named_source_value": False,
            "s1_prefix_evidence_is_not_pulse_mass": True,
            "exact_branch_source": mass_source,
            "phase_zero_exact_mass": mass_source.get("phase_zero_exact_zero") is True,
            "phase_one_exact_mass": mass_source.get("phase_one_exact_half") is True,
            "interior_mass_source_unresolved": True,
        },
        "weighted_primitive_remainder_per_hb": {
            "source": "original_first_switch_m15.source_norms",
            "source_file": "src/openai_ns_reconstruction/original_first_switch_m15.py:119-127",
            "formula": "s*mass*weighted_derivative + hcap*angular*exp(ell)*mass*weighted",
            "serialized_value_selected": False,
            "available_as_named_source_value": False,
        },
        "delta": dict(delta_binding),
        "modified_Qy_not_used_for_leading_source": True,
        "ordinary_Z1_retained": True,
        "passed": bool(geometry.get("passed") and macro_endpoint_rows and first_sample_count),
    }


def _mapping_fingerprint(value: Mapping[str, Any], binding_key: str | None = None) -> str:
    binding = _mapping(value.get(binding_key)) if binding_key else {}
    return _digest({
        "source_family": _mapping(value.get("source_family")),
        "top_level_keys": sorted(str(key) for key in value),
        "binding_keys": sorted(str(key) for key in binding),
    })


def _source_pins(
    repo: Path,
    report_meta: Mapping[str, Mapping[str, Any]],
    archive_meta: Mapping[str, Mapping[str, Any]],
    archive_values: Mapping[str, Mapping[str, Any]],
    parent_report: Any,
    delta_source: Mapping[str, Any],
    mass_source: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return filesystem pins separately from owner-supplied mappings."""
    pins: dict[str, Any] = {}
    archive_inputs: dict[str, Any] = {}

    for name, meta in report_meta.items():
        path, digest = meta.get("path"), meta.get("sha256")
        pins[name] = {"path": path, "sha256": digest, "passed": bool(path and digest and meta.get("passed") is True)}
    for name, meta in archive_meta.items():
        path, digest = meta.get("path"), meta.get("sha256")
        if meta.get("mode") == "supplied_mapping":
            value = archive_values.get(name, {})
            archive_inputs[name] = {
                "mode": "supplied_mapping",
                "digest": _mapping_fingerprint(value, _ARCHIVE_KEYS.get(name.removesuffix("_archive"))),
                "passed": meta.get("passed") is True,
            }
        else:
            pins[name] = {"path": path, "sha256": digest, "passed": bool(path and digest and meta.get("passed") is True)}

    if isinstance(parent_report, (str, Path)):
        path = Path(parent_report)
        if not path.is_absolute():
            path = repo / path
        digest = _sha256(path)
        pins["parent_report"] = {"path": _relative(path, repo), "sha256": digest, "passed": digest is not None}
    elif isinstance(parent_report, Mapping):
        archive_inputs["parent_report"] = {
            "mode": "supplied_mapping",
            "digest": _digest({key: parent_report.get(key) for key in ("source_family", "checks", "source_derivative_theorem", "source_seam_proofs")} ),
            "passed": True,
        }
    else:
        archive_inputs["parent_report"] = {"mode": "missing", "passed": False}

    pins["delta_source"] = {
        "path": delta_source.get("path"),
        "sha256": delta_source.get("sha256"),
        "passed": delta_source.get("passed") is True,
    }
    pins["mass_source"] = {
        "path": mass_source.get("path"),
        "sha256": mass_source.get("sha256"),
        "passed": mass_source.get("passed") is True,
    }
    return pins, archive_inputs


def build_switch_leading_source(
    *,
    nodes: Any,
    loop_archive: Any,
    modified_archive: Any,
    high_archive: Any,
    parent_report: Any,
    directory: Any,
    root: Any,
) -> dict[str, Any]:
    """Build the source-bound first-switch leading point/remainder view.

    The function accepts mappings so a caller can pass already-loaded pinned
    archives.  Paths are loaded as JSON/gzip only.  No accepted producer,
    owner, checker, fixture, or numerical ancestor is imported or executed.
    """
    repo = _repo_root(root, directory)
    loop, loop_meta = _load_serialized(loop_archive, repo, name="loop_archive", default=repo / "experiments" / "root_st073" / (_PREFIX + "prescribed_loop_C5_functions.json.gz"))
    modified, modified_meta = _load_serialized(modified_archive, repo, name="modified_archive", default=repo / "experiments" / "root_st073" / (_PREFIX + "modified_Q4_pressure_C3_functions.json.gz"))
    high, high_meta = _load_serialized(high_archive, repo, name="high_archive", default=repo / "experiments" / "root_st073" / (_PREFIX + "prescribed_high_radial_source_functions.json.gz"))
    first_pressure, first_pressure_meta = _load_report(directory, repo, "first_switch_pressure")
    first_m15, first_m15_meta = _load_report(directory, repo, "first_switch_m15")
    macro_pressure, macro_pressure_meta = _load_report(directory, repo, "macro_pressure")
    macro_m15, macro_m15_meta = _load_report(directory, repo, "macro_m15")
    second_micro_m15, second_micro_meta = _load_report(directory, repo, "second_micro_m15")

    report_meta = {
        "first_switch_pressure": first_pressure_meta,
        "first_switch_m15": first_m15_meta,
        "macro_pressure": macro_pressure_meta,
        "macro_m15": macro_m15_meta,
        "second_micro_m15": second_micro_meta,
    }
    archive_meta = {"loop_archive": loop_meta, "modified_archive": modified_meta, "high_archive": high_meta}
    archive_values = {"loop_archive": loop, "modified_archive": modified, "high_archive": high}
    delta_source = _delta_source_guard(repo)
    mass_source = _mass_source_guard(repo)
    pins, archive_inputs = _source_pins(
        repo,
        report_meta,
        archive_meta,
        archive_values,
        parent_report,
        delta_source,
        mass_source,
    )
    family_rows = {
        name: {"passed": _family_ok(data), "source_family": _mapping(data.get("source_family"))}
        for name, data in {
            "loop_archive": loop,
            "modified_archive": modified,
            "high_archive": high,
            "first_switch_pressure": first_pressure,
            "first_switch_m15": first_m15,
            "macro_pressure": macro_pressure,
            "macro_m15": macro_m15,
            "second_micro_m15": second_micro_m15,
        }.items()
    }
    ctx, iv, delta, delta_iv, delta_binding = _contexts(repo)
    bindings = _source_bindings(
        nodes=nodes,
        loop=loop,
        modified=modified,
        high=high,
        first_pressure=first_pressure,
        macro_pressure=macro_pressure,
        second_micro_m15=second_micro_m15,
        pins=pins,
        delta_binding=delta_binding,
        mass_source=mass_source,
    )

    first_index = _samples(first_pressure)
    first_m15_index = _samples(first_m15)
    macro_index = _samples(macro_pressure)
    rows: list[dict[str, Any]] = []
    blocked: list[dict[str, Any]] = []
    comparison_count = 0
    comparison_passed = 0
    computed_count = 0
    endpoint_empty_count = 0
    endpoint_failure_count = 0
    integral_c0_comparison_count = 0
    integral_c0_comparison_passed = 0
    for key, pressure_sample in first_index.items():
        z, coordinate = key
        m15_sample = first_m15_index.get(key)
        macro_sample = macro_index.get((z, "1"))
        upper_source = _source_rows(pressure_sample)
        lower_source = _source_rows(macro_sample) if macro_sample is not None else []
        raw_upper = _m15_term(m15_sample or {}, "first_switch_viscous_upper")
        raw_lower = _m15_term(m15_sample or {}, "first_switch_viscous_lower")
        empty = coordinate in {"0", "0.0", "+0"}
        for endpoint, source_name, source_rows, raw_term, sign, hb in (
            ("upper", "first_switch_upper_source_Uy", upper_source, raw_upper, -1, "hbS"),
            ("lower", "first_switch_lower_source_Uy", lower_source, raw_lower, 1, "hbB"),
        ):
            row: dict[str, Any] = {
                "name": "first_switch_viscous_" + endpoint,
                "raw_number": 43 if endpoint == "upper" else 44,
                "Z": z,
                "coordinate": copy.deepcopy(pressure_sample.get("coordinate")),
                "source_binding": source_name,
                "sign": sign,
                "formula": "sign * 2*Z*B/(1-delta*Z^2), ordinary Z0/Z1",
                "physical_scale": {
                    "R_factor_once": True,
                    "reference_radius": 100,
                    "reference_radius_source": "physical(1,1)=log(100)-logLambda",
                    "actual_radius_node": 85905 if endpoint == "upper" else 100,
                    "actual_radius_difference_retained_in_raw45": endpoint == "upper",
                    "radius_binding": (
                        "raw43 uses external R100; actual R(s)-100 is retained in raw45"
                        if endpoint == "upper"
                        else "raw44 uses frozen_macro R100 exactly"
                    ),
                    "functional_identity": "E=2Z*R(s)*Uy/L is a separate root source identity",
                    "epsilon_factor_once": True,
                    "width_factor": hb,
                    "formula": f"R100*epsilon*coefficient({endpoint}_source_Uy)*{hb}",
                },
                "source_point": [_source_summary(item) for item in source_rows] if source_rows else [],
                "source_truth": [_source_summary(item) for item in source_rows] if source_rows else [],
                "independent_saved_row": copy.deepcopy(raw_term.get("coefficient_C0_C1")) if raw_term else None,
                "approximation_is_not_exact_source_value": True,
            }
            if len(source_rows) >= 2 and raw_term is not None:
                try:
                    computed = _computed_endpoint(ctx, iv, z, delta, delta_iv, source_rows, sign)
                    row["recomputed"] = computed["point"]
                    row["recomputed_truth"] = computed["truth"]
                    row["kernel"] = "original_leading_viscous_kernel.endpoint_coefficient"
                    comparison = _compare_endpoint(ctx, iv, computed, raw_term.get("coefficient_C0_C1", ()), empty=empty)
                    row["comparison"] = comparison
                    computed_count += 1
                    if comparison.get("status") == "compared_independent_saved_row":
                        comparison_count += 1
                        if comparison.get("passed"):
                            comparison_passed += 1
                    elif comparison.get("status") == "not_compared_empty_increment":
                        endpoint_empty_count += 1
                    else:
                        endpoint_failure_count += 1
                except (TypeError, ValueError, ArithmeticError) as exc:
                    row["recomputed"] = None
                    row["recomputed_truth"] = None
                    row["comparison"] = {"passed": False, "status": "source_arithmetic_failed", "error": str(exc)}
                    endpoint_failure_count += 1
            else:
                row["recomputed"] = None
                row["recomputed_truth"] = None
                row["comparison"] = {"passed": False, "status": "missing_source_or_independent_row"}
                endpoint_failure_count += 1
            rows.append(row)

        raw_radius = _m15_term(m15_sample or {}, "first_switch_viscous_upper_radius_error")
        rows.append({
            "name": "first_switch_viscous_upper_radius_error",
            "raw_number": 45,
            "Z": z,
            "coordinate": copy.deepcopy(pressure_sample.get("coordinate")),
            "kind": "positive_radius_remainder_cap",
            "source_binding": "first_switch_upper_source_Uy.truth",
            "physical_scale": {
                "R_factor_once": True,
                "reference_radius": 100,
                "actual_radius_node": 85905,
                "radius_difference_source": "actual R(s)-100 retained in raw45",
                "epsilon_factor_once": True,
                "width_factors": ["hbS", "hbS"],
                "formula": "(actual R(s)-100)*epsilon*hbS^2*source_truth",
            },
            "source_point_is_not_selected": True,
            "independent_saved_row": copy.deepcopy(raw_radius.get("coefficient_C0_C1")) if raw_radius else None,
            "independent_saved_log_scale": copy.deepcopy(raw_radius.get("log_scale_terms")) if raw_radius else None,
            "comparison": {"passed": bool(raw_radius), "status": "retained_positive_or_zero_cap" if raw_radius else "missing_independent_row"},
        })

        raw_integral = _m15_term(m15_sample or {}, "first_switch_viscous_source_integral")
        integral_row: dict[str, Any] = {
            "name": "first_switch_viscous_source_integral",
            "raw_number": 46,
            "Z": z,
            "coordinate": copy.deepcopy(pressure_sample.get("coordinate")),
            "kind": "signed_point_plus_remainder_source_integral",
            "source_binding": "first_switch_lower_source_Uy",
            "physical_scale": {
                "R_factor_once": True,
                "reference_radius": 100,
                "R_source": "serialized physical(1,1)=log(100)-logLambda; root functional E uses actual R(s) separately",
                "epsilon_factor_once": True,
                "width_factors": ["hbS", "hbS"],
                "formula": "R100*epsilon*hbS^2*C(B)*mass",
            },
            "independent_saved_row": copy.deepcopy(raw_integral.get("coefficient_C0_C1")) if raw_integral else None,
            "independent_saved_log_scale": copy.deepcopy(raw_integral.get("log_scale_terms")) if raw_integral else None,
            "approximation_is_not_exact_source_value": True,
        }
        phase_zero = empty
        phase_one = coordinate in {"1", "1.0", "+1"}
        if phase_zero:
            zero_pair = _exact_source_pair(ctx, iv, 0, reason="pulse_masses hi==0 exact branch")
            integral_row.update({
                "mass": zero_pair,
                "signed_integral_C0": zero_pair,
                "signed_integral_C1": zero_pair,
                "status": "empty_exact_zero_scope",
                "empty_scope_is_not_a_H_or_mean_zero": True,
                "comparison": {
                    "passed": bool(raw_integral),
                    "status": "not_compared_empty_integral",
                    "reason": "The saved s0 increment is an exact zero branch; it is not a function-zero proof for H or the mean.",
                },
            })
        elif len(lower_source) >= 2:
            try:
                lower_point = [_source_pair(ctx, iv, item)[0] for item in lower_source[:2]]
                lower_truth = [_source_pair(ctx, iv, item)[1] for item in lower_source[:2]]
                c0_point = mean_source_value(ctx.mpf(z), delta, lower_point)
                c0_truth = mean_source_value(iv.mpf(z), delta_iv, lower_truth)
                c0_without_mass = encode_pair(ctx, iv, c0_point, c0_truth)
                integral_row["C0_recomputed_without_mass"] = c0_without_mass
                if phase_one and mass_source.get("phase_one_exact_half") is True:
                    mass_pair = _exact_source_pair(ctx, iv, "0.5", reason="pulse_masses lo==hi==1 exact branch")
                    integral_point = c0_point * ctx.mpf("0.5")
                    integral_truth = c0_truth * iv.mpf("0.5")
                    integral_c0 = encode_pair(ctx, iv, integral_point, integral_truth)
                    integral_row.update({
                        "mass": mass_pair,
                        "signed_integral_C0": integral_c0,
                        "signed_integral_C1": None,
                        "status": "signed_C0_computed_C1_blocked",
                        "mass_scope": "exact_phase_one_half_only",
                        "approximation_is_not_exact_source_value": True,
                    })
                    raw_rows = raw_integral.get("coefficient_C0_C1", ()) if raw_integral else ()
                    if isinstance(raw_rows, (list, tuple)) and raw_rows and isinstance(raw_rows[0], Mapping):
                        comparison = _compare_scalar(ctx, integral_c0, raw_rows[0])
                        integral_row["comparison"] = comparison
                        integral_c0_comparison_count += 1
                        if comparison.get("passed") is True:
                            integral_c0_comparison_passed += 1
                    else:
                        integral_row["comparison"] = {
                            "passed": False,
                            "status": "missing_independent_saved_C0_row",
                        }
                else:
                    integral_row["status"] = "blocked_missing_source_inputs"
            except (TypeError, ValueError, ArithmeticError) as exc:
                integral_row["C0_recomputed_without_mass"] = {"error": str(exc)}
                integral_row["status"] = "blocked_source_arithmetic_failed"
        else:
            integral_row["status"] = "blocked_missing_source_inputs"
        rows.append(integral_row)
        if not phase_zero:
            missing = ["B_ZZ2 for ordinary C1"]
            if not phase_one or mass_source.get("phase_one_exact_half") is not True:
                missing.extend(["pulse mass", "signed truth mass product"])
            blocked.append({
                "name": "first_switch_viscous_source_integral",
                "raw_number": 46,
                "Z": z,
                "coordinate": coordinate,
                "missing": missing,
                "available": ["C(B) C0 recipe from signed lower B/B_Z", "saved signed raw46 point/truth row"],
                "reason": (
                    "C(B) ordinary Z1 requires explicit B_ZZ2; phase-one mass is the exact source half branch."
                    if phase_one and mass_source.get("phase_one_exact_half") is True
                    else "C(B) ordinary Z1 requires explicit B_ZZ2; interior pulse mass and its source enclosure are not serialized."
                ),
            })

        raw_remainder = _m15_term(m15_sample or {}, "first_switch_viscous_source_integral_error")
        rows.append({
            "name": "first_switch_viscous_source_integral_error",
            "raw_number": 47,
            "Z": z,
            "coordinate": copy.deepcopy(pressure_sample.get("coordinate")),
            "kind": "positive_source_remainder_cap",
            "source_binding": "weighted_primitive_remainder_per_hb",
            "physical_scale": {
                "R_factor_once": True,
                "reference_radius": 100,
                "R_source": "serialized physical(1,1)=log(100)-logLambda; actual R(s) belongs to root functional binding",
                "epsilon_factor_once": True,
                "width_factors": ["hbS", "hbS", "hbS"],
                "formula": "R100*epsilon*hbS^3*symmetric(C(weighted_primitive_remainder_per_hb))",
            },
            "source_point_is_not_selected": True,
            "independent_saved_row": copy.deepcopy(raw_remainder.get("coefficient_C0_C1")) if raw_remainder else None,
            "independent_saved_log_scale": copy.deepcopy(raw_remainder.get("log_scale_terms")) if raw_remainder else None,
            "comparison": {"passed": bool(raw_remainder), "status": "retained_positive_remainder_cap" if raw_remainder else "missing_independent_row"},
        })
        blocked.append({
            "name": "first_switch_viscous_source_integral_error",
            "raw_number": 47,
            "Z": z,
            "coordinate": coordinate,
            "missing": ["weighted_primitive_remainder_per_hb", "B_ZZ2 for C1", "pulse mass"],
            "available": ["saved positive raw47 enclosure"],
            "reason": "Raw47 is an error enclosure only. Its positive value cannot be promoted to a source function or used to infer raw46.",
        })

    first_s1 = [
        sample
        for (z, coordinate), sample in first_index.items()
        if coordinate in {"1", "1.0", "+1"}
    ]
    s1_flags = [
        _mapping(sample.get("leading_U_y_factored")).get("actual_unmodified_U_y_exact_zero_at_s1") is True
        for sample in first_s1
    ]
    bindings["first_switch_exact_control_zero_at_s1"] = {
        "passed": bool(first_s1) and all(s1_flags),
        "sample_count": len(first_s1),
        "theorem_field": "leading_U_y_factored.actual_unmodified_U_y_exact_zero_at_s1",
        "not_inferred_from_serialized_raw43": True,
    }
    family_passed = all(row["passed"] for row in family_rows.values())
    source_pins_passed = bool(pins) and all(row.get("passed") is True for row in pins.values())
    source_bindings_passed = bool(
        bindings.get("passed")
        and bindings["first_switch_exact_control_zero_at_s1"]["passed"]
        and bindings["delta"].get("source_pin", {}).get("passed") is True
        and bindings["mass"].get("exact_branch_source", {}).get("passed") is True
    )
    endpoint_interface_passed = bool(
        computed_count > 0
        and comparison_count > 0
        and comparison_count == comparison_passed
        and endpoint_failure_count == 0
        and computed_count == comparison_count + endpoint_empty_count
    )
    integral_c0_interface_passed = bool(
        integral_c0_comparison_count > 0
        and integral_c0_comparison_count == integral_c0_comparison_passed
    )
    limits = {
        "source_family_passed": family_passed,
        "source_pins_passed": source_pins_passed,
        "source_bindings_passed": source_bindings_passed,
        "endpoint_interface_passed": endpoint_interface_passed,
        "integral_C0_interface_passed": integral_c0_interface_passed,
        "endpoint_rows_computed": computed_count,
        "endpoint_rows_compared": comparison_count,
        "endpoint_rows_comparison_passed": comparison_passed,
        "endpoint_rows_not_compared_empty": endpoint_empty_count,
        "endpoint_rows_failed": endpoint_failure_count,
        "endpoint_kernel_complete": endpoint_interface_passed,
        "integral_C0_rows_compared": integral_c0_comparison_count,
        "integral_C0_rows_comparison_passed": integral_c0_comparison_passed,
        "integral_C0_interface_passed": integral_c0_interface_passed,
        "integral_C0_partial_only": True,
        "integral_C1_requires_B_ZZ2": True,
        "mass_source_serialized": False,
        "mass_exact_phase_zero_branch_bound": bindings["mass"]["phase_zero_exact_mass"],
        "mass_exact_phase_one_branch_bound": bindings["mass"]["phase_one_exact_mass"],
        "mass_interior_unresolved": True,
        "weighted_remainder_source_serialized": False,
        "physical_R_factor_once": True,
        "hbB_and_hbS_distinct": True,
        "ordinary_Z1_retained": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "complete": False,
    }
    return {
        "passed": bool(family_passed and source_pins_passed and source_bindings_passed and endpoint_interface_passed and not blocked),
        "partial": endpoint_interface_passed,
        "source_bindings_passed": source_bindings_passed,
        "endpoint_interface_passed": endpoint_interface_passed,
        "integral_C0_interface_passed": integral_c0_interface_passed,
        "complete": False,
        "source_family": dict(_FAMILY),
        "source_pins": pins,
        "archive_inputs": archive_inputs,
        "family_checks": family_rows,
        "source_bindings": bindings,
        "point_remainder_rows": rows,
        "blocked_rows": blocked,
        "limits": limits,
    }


__all__ = ["build_switch_leading_source", "endpoint_coefficient"]
