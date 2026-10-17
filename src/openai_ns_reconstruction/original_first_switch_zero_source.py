"""Source-only proof for the first-switch axial source at ``s = 1``.

The first-switch pressure packet contains a useful endpoint flag, but that
flag is not the proof.  This module binds the flag to the source bodies and
to the exact endpoint branch of the original flat cutoff.  It reads saved
JSON and source text only; it does not import an owner, producer, fixture,
field evaluator, quadrature routine, or numerical ancestor.

The conclusion is deliberately narrow: the actual first-switch ``U_y``
source atom ``E(B_up)`` is zero at ``s = 1`` as a function of the valid
ordinary ``Z`` variable.  It says nothing about the actual mean or an ``H``
field.  Positive radius/error caps remain enclosures and are never used as
the zero proof.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


_PARENT_RELATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_macro_coupled_m15_numeric.json"
)
_PARENT_SHA256 = "b63476443d5de5c1fca24ed5979ccb8a71273c42ff7a0b9b9094d8f80782f177"
_PREFIX = "lei_ren_part1_paper_compliant_current_original_"
_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_REPORTS = {
    "first_switch_pressure_rows": _PREFIX + "first_switch_pressure_rows.json",
    "first_switch_numeric": _PREFIX + "first_switch_numeric.json",
    "first_switch_m15": _PREFIX + "first_switch_m15_numeric.json",
}
_IMPLEMENTATION_SOURCES = {
    "pressure_rows": "src/openai_ns_reconstruction/original_first_switch_pressure_rows.py",
    "first_switch_numeric": "src/openai_ns_reconstruction/original_first_switch_numeric.py",
    "first_switch_m15": "src/openai_ns_reconstruction/original_first_switch_m15.py",
    "macro_pressure_rows": "src/openai_ns_reconstruction/original_macro_pressure_rows.py",
    "first_collar_transport": "src/openai_ns_reconstruction/original_first_collar_transport.py",
}
_INPUT_SOURCES = {
    "flat_sigma": "experiments/root_st073/lei_ren_part1_paper_compliant_flat_pulse_derivatives.py",
    "prescribed_first_switch": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
    ),
    "physical_ftc": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions.py"
    ),
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
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


def _work_directory(repo: Path, directory: Any) -> Path:
    if directory is None:
        return repo / "experiments" / "root_st073"
    path = Path(directory)
    return path.parent if path.is_file() else path


def _relative(path: Path, repo: Path) -> str:
    try:
        return path.relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _resolve(repo: Path, path: Any, *, default: Path | None = None) -> Path:
    if path is None:
        if default is None:
            return repo
        return default
    value = Path(path)
    if value.is_absolute():
        return value
    # Callers commonly pass either a repository-relative path or a path
    # relative to experiments/root_st073.  Prefer an existing repository path.
    candidate = repo / value
    if candidate.exists() or not default:
        return candidate
    return default / value


def _load_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _empty_failure(
    gaps: list[str], *, repo: Path, directory: Path, parent_pin: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    return {
        "passed": False,
        "source_family": dict(_FAMILY),
        "compact_proofs_passed": False,
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": dict(parent_pin or {"passed": False}),
        "parent_mapping_hash_bind": {
            "passed": False,
            "supplied_mapping_sha256": None,
            "latest_mapping_sha256": None,
        },
        "input_hashes": {},
        "implementation_hashes": {},
        "source_pins": {},
        "AST": {},
        "ast_function_hashes": {},
        "inherited_fields_digests": {},
        "ordinary_Z_orders": [0, 1],
        "ordinary_orders": {"required": [0, 1], "available": [], "passed": False},
        "exact_zero_identity": {
            "passed": False,
            "source_atom": "E(B_up)",
            "kappa": "1-sigma(s)",
            "truth_source_multiplied_before_ordinary_Z": False,
            "valid_Z_domain": "-1 < Z < 1",
        },
        "mapping": {},
        "zero_scope": {
            "E_source_zero": False,
            "actual_mean_zero_claim": False,
            "H_field_zero_claim": False,
            "own_mean_retained": True,
        },
        "sourceflagtruth": {
            "passed": False,
            "truth_source_multiplied_before_ordinary_Z": False,
            "normalization_error_zero_required": False,
        },
        "sigmaendpointbranch": {
            "passed": False,
            "identity": "sigma(1)=1",
        },
        "ordinaryZsourcezeroorders": {
            "passed": False,
            "orders": [0, 1],
            "valid_Z_domain": "-1 < Z < 1",
        },
        "scopephase1": {
            "passed": False,
            "phase": "s=1",
            "source_atom": "E(B_up)",
        },
        "capnotselected": {
            "passed": False,
            "positive_cap_value_zero_claim": False,
        },
        "finiteN_meannotzeroed": {
            "passed": False,
            "finite_N_corrections_retained": True,
            "actual_mean_zero_claim": False,
            "H_field_zero_claim": False,
        },
        "source_executed": False,
        "numeric_field_execution": False,
        "root": str(repo),
        "directory": str(directory),
        "gaps": gaps,
    }


def _function_node(tree: ast.AST, function: str, class_name: str | None) -> ast.AST:
    scope: ast.AST = tree
    if class_name is not None:
        classes = [node for node in ast.walk(tree) if isinstance(node, ast.ClassDef) and node.name == class_name]
        if len(classes) != 1:
            raise ValueError(f"unique class required: {class_name}")
        scope = classes[0]
    nodes = [
        node
        for node in ast.walk(scope)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == function
    ]
    if len(nodes) != 1:
        raise ValueError(f"unique function required: {function}")
    return nodes[0]


def _parsed_fragment(fragment: str) -> ast.AST:
    tree = ast.parse(fragment)
    if len(tree.body) != 1:
        raise ValueError(f"one AST fragment required: {fragment}")
    return tree.body[0]


def _contains(node: ast.AST, wanted: ast.AST) -> bool:
    target = ast.dump(wanted, include_attributes=False)
    return any(ast.dump(part, include_attributes=False) == target for part in ast.walk(node))


def _ast_guard(
    path: Path,
    repo: Path,
    function: str,
    fragments: tuple[str, ...],
    *,
    class_name: str | None = None,
    expressions: tuple[str, ...] = (),
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "function": function,
        "class": class_name,
        "required_AST": {fragment: False for fragment in fragments},
        "required_expressions": {expression: False for expression in expressions},
        "function_AST_sha256": None,
        "source_path": _relative(path, repo),
        "source_executed": False,
        "passed": False,
    }
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        function_node = _function_node(tree, function, class_name)
        for fragment in fragments:
            row["required_AST"][fragment] = _contains(function_node, _parsed_fragment(fragment))
        for expression in expressions:
            row["required_expressions"][expression] = _contains(function_node, ast.parse(expression, mode="eval").body)
        row["function_AST_sha256"] = hashlib.sha256(
            ast.dump(function_node, include_attributes=False).encode("utf-8")
        ).hexdigest()
        row["passed"] = all(row["required_AST"].values()) and all(row["required_expressions"].values())
    except (OSError, SyntaxError, UnicodeError, ValueError, StopIteration) as exc:
        row["error"] = str(exc)
    return row


def _parent_pin(path: Path, repo: Path, parent_group: Mapping[str, Any], *, kind: str, name: str) -> dict[str, Any]:
    relative = _relative(path, repo)
    key = relative if kind == "implementation" else path.name
    expected = parent_group.get(key)
    actual = _sha256(path)
    return {
        "name": name,
        "kind": kind,
        "key": key,
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "parent_pin_present": isinstance(expected, str),
        "passed": isinstance(expected, str) and actual == expected,
        "source_executed": False,
    }


def _sample_coordinate(sample: Mapping[str, Any]) -> str:
    value = sample.get("coordinate")
    if isinstance(value, Mapping):
        return str(value.get("value"))
    return str(value)


def _samples_for_z(report: Mapping[str, Any], z: str, coordinate: str | None = None) -> list[Mapping[str, Any]]:
    result = []
    for sample in report.get("samples", ()):
        if not isinstance(sample, Mapping) or str(sample.get("Z")) != str(z):
            continue
        if coordinate is not None and _sample_coordinate(sample) not in {coordinate, "1.0", "+1"}:
            continue
        result.append(sample)
    return result


def _family_ok(report: Mapping[str, Any]) -> bool:
    return _mapping(report.get("source_family")) == _FAMILY


def _json_pin_for_report(
    report: Mapping[str, Any], name: str, path: Path, expected: str | None
) -> dict[str, Any]:
    actual = _sha256(path)
    return {
        "name": name,
        "kind": "input",
        "key": path.name,
        "path": path.as_posix(),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": isinstance(report, Mapping) and actual == expected,
        "source_executed": False,
    }


def _digest_fields(parent: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    gaps: list[str] = []
    fields = (
        "inherited_m15_packet",
        "signed_m15_query_views",
        "macro_weighted_source_rows",
        "native_macro_seam",
        "signed_ledger_binding",
    )
    result: dict[str, Any] = {}
    for name in fields:
        value = parent.get(name)
        if value is None:
            gaps.append(f"parent field missing for digest: {name}")
            continue
        calculated = _digest(value)
        row: dict[str, Any] = {"sha256": calculated, "present": True}
        stored_name = name + "_digest"
        if stored_name in parent:
            row["stored_sha256"] = parent.get(stored_name)
            row["stored_matches"] = parent.get(stored_name) == calculated
            if not row["stored_matches"]:
                gaps.append(f"parent stored digest mismatch: {stored_name}")
        result[name] = row
    return result, gaps


def _load_parent_report(parent_report: Any, repo: Path, directory: Path) -> tuple[Mapping[str, Any], dict[str, Any], list[str]]:
    gaps: list[str] = []
    latest_path = repo / _PARENT_RELATIVE
    latest_sha = _sha256(latest_path)
    pin = {
        "path": _relative(latest_path, repo),
        "sha256": latest_sha,
        "expected_sha256": _PARENT_SHA256,
        "passed": latest_sha == _PARENT_SHA256,
        "source_executed": False,
    }
    if not pin["passed"]:
        gaps.append(f"latest parent report pin mismatch: expected {_PARENT_SHA256}, got {latest_sha}")
    try:
        latest = _load_json(latest_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {}, pin, gaps + [f"cannot load latest parent report: {exc}"]

    supplied: Mapping[str, Any]
    if isinstance(parent_report, Mapping):
        supplied = parent_report
    elif isinstance(parent_report, (str, Path)):
        candidate = _resolve(repo, parent_report, default=directory)
        try:
            supplied = _load_json(candidate)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            return {}, pin, gaps + [f"cannot load supplied parent report: {exc}"]
    else:
        return {}, pin, gaps + ["latest macro-coupled parent mapping is required"]

    latest_mapping_sha256 = _digest(latest)
    supplied_mapping_sha256 = _digest(supplied)
    pin["latest_mapping_sha256"] = latest_mapping_sha256
    pin["supplied_mapping_sha256"] = supplied_mapping_sha256
    pin["mapping_match"] = supplied_mapping_sha256 == latest_mapping_sha256
    if not pin["mapping_match"]:
        gaps.append("supplied parent mapping does not match the pinned latest parent report")
    return supplied, pin, gaps


def build_first_switch_zero_source(parent_report: Any, *, root: Any, directory: Any) -> dict[str, Any]:
    """Return a compact source-bound proof of ``U_y(1, Z) == 0``.

    All source and report reads are fail-closed.  The endpoint rows are used
    as coverage and provenance evidence only; no serialized approximation,
    cap, midpoint, underflow result, or interval overlap is used as an
    identity.
    """

    repo = _repo_root(root, directory)
    work = _work_directory(repo, directory)
    parent, parent_pin, gaps = _load_parent_report(parent_report, repo, work)
    if not parent:
        return _empty_failure(gaps, repo=repo, directory=work, parent_pin=parent_pin)

    source_family = _mapping(parent.get("source_family"))
    if source_family != _FAMILY:
        gaps.append("parent source_family does not match the active original family")
    checks = _mapping(parent.get("checks"))
    if checks.get("passed") is not True:
        gaps.append("parent checks.passed is not true")
    if _mapping(parent.get("contract")).get("raw43_raw45_and_all_other_errors_retained") is not True:
        gaps.append("parent raw43/raw45 retention contract is missing")

    parent_inputs = _mapping(parent.get("input_hashes"))
    parent_implementations = _mapping(parent.get("implementation_hashes"))
    if not parent_inputs:
        gaps.append("parent input_hashes mapping is missing")
    if not parent_implementations:
        gaps.append("parent implementation_hashes mapping is missing")

    source_pins: dict[str, Any] = {}
    input_hashes: dict[str, str] = {}
    implementation_hashes: dict[str, str] = {}
    paths: dict[str, Path] = {}
    for semantic, relative in _IMPLEMENTATION_SOURCES.items():
        path = repo / relative
        pin = _parent_pin(path, repo, parent_implementations, kind="implementation", name=semantic)
        source_pins[semantic] = pin
        paths[semantic] = path
        if pin["passed"]:
            implementation_hashes[pin["key"]] = pin["sha256"]
        else:
            gaps.append(f"implementation source pin failed: {relative}")
    for semantic, relative in _INPUT_SOURCES.items():
        path = repo / relative
        pin = _parent_pin(path, repo, parent_inputs, kind="input", name=semantic)
        source_pins[semantic] = pin
        paths[semantic] = path
        if pin["passed"]:
            input_hashes[pin["key"]] = pin["sha256"]
        else:
            gaps.append(f"input source pin failed: {relative}")

    # The inherited endpoint packet is needed only to verify the two
    # unchanged slots.  It remains a read-only source pin from the parent.
    inherited_packet = _mapping(parent.get("inherited_m15_packet"))
    inherited_name = inherited_packet.get("report")
    inherited_path = work / inherited_name if isinstance(inherited_name, str) else work / "missing-inherited-endpoint.json"
    inherited_expected = inherited_packet.get("report_sha256")
    inherited_report: Mapping[str, Any] = {}
    if isinstance(inherited_name, str):
        inherited_pin = _parent_pin(
            inherited_path, repo, parent_inputs, kind="input", name="inherited_endpoint_report"
        )
        source_pins["inherited_endpoint_report"] = inherited_pin
        if inherited_pin["passed"]:
            input_hashes[inherited_pin["key"]] = inherited_pin["sha256"]
        else:
            gaps.append("inherited endpoint report pin failed")
        if inherited_pin["sha256"] != inherited_expected:
            gaps.append("inherited endpoint report SHA disagrees with parent packet")
        try:
            inherited_report = _load_json(inherited_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            gaps.append(f"cannot load inherited endpoint report: {exc}")
    else:
        gaps.append("parent inherited_m15_packet report name is missing")

    report_values: dict[str, Mapping[str, Any]] = {}
    for semantic, basename in _REPORTS.items():
        path = work / basename
        expected = parent_inputs.get(path.name)
        try:
            report = _load_json(path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            report = {}
            gaps.append(f"cannot load saved {semantic} report: {exc}")
        report_values[semantic] = report
        pin = _json_pin_for_report(report, semantic, path, expected)
        source_pins[semantic + "_report"] = pin
        if pin["passed"]:
            input_hashes[path.name] = pin["sha256"]
        else:
            gaps.append(f"saved {semantic} report pin failed")

    inherited_digests, digest_gaps = _digest_fields(parent)
    gaps.extend(digest_gaps)

    # Source bodies are parsed only.  These fragments establish the factor,
    # endpoint branch, valid Z domain, and the extraction order into C0/C1.
    ast_rows: dict[str, Any] = {
        "pressure_rows.inputs": _ast_guard(
            paths["pressure_rows"],
            repo,
            "inputs",
            (
                "incoming = self.macro_rows.inputs(Z, iv.mpf(1))",
                "source = self.source_Uy(o._prepared[Z], geo, 5, phase=s, ell_switch=ell)",
                "point['Psi_y'], truth['Psi_y'] = source['normalized_value_point'], source['normalized_value_truth']",
                "out = dict(point=point, truth=truth, geometry=geo, source_Uy=source, incoming=incoming, switch_evidence=evidence, switch_Psi_change_error=psi_error, switch_mean_change_error=mean_error)",
            ),
            class_name="OriginalFirstSwitchPressureRows",
        ),
        "pressure_rows.source_Uy": _ast_guard(
            paths["pressure_rows"],
            repo,
            "source_Uy",
            (
                "complement = 1 - self.t.sigma_jets(iv, phase)[0]",
                "truth = (raw['truth'] + o._norm_error_jet(extra, order)) * complement",
            ),
            class_name="OriginalFirstSwitchPressureRows",
            expressions=(
                "{name: row * complement for name, row in raw['parts'].items()}",
                "J.constant(c, 0, order)",
            ),
        ),
        "pressure_rows.evaluate": _ast_guard(
            paths["pressure_rows"],
            repo,
            "evaluate",
            (),
            class_name="OriginalFirstSwitchPressureRows",
            expressions=(
                "s._mpi_ == self.iv.mpf(1)._mpi_",
                "'100*exp(hbS*s)'",
                "'dy=hbS ds; R*dR=R^2*hbS ds'",
            ),
        ),
        "first_switch_numeric.coordinate": _ast_guard(
            paths["first_switch_numeric"],
            repo,
            "_coordinate",
            (
                "if lo < 0 or hi > 1: raise ValueError('First-switch phase in[0,1] required')",
            ),
            class_name="OriginalFirstSwitchNumeric",
        ),
        "first_switch_m15.leading": _ast_guard(
            paths["first_switch_m15"],
            repo,
            "leading",
            (
                "low = macro['source_Uy']",
                "up = v['source_Uy']",
                "endpoints = [(coefficient(c, up['point'], o.core.delta) * (-1), coefficient(iv, up['truth'], o.core.delta_iv) * (-1)), (coefficient(c, low['point'], o.core.delta), coefficient(iv, low['truth'], o.core.delta_iv))]",
                "radius_error = self.symmetric(coefficient(iv, up['truth'], o.core.delta_iv)) * s * iv.exp(o.hcap * s)",
            ),
            class_name="OriginalFirstSwitchM15",
        ),
        "first_switch_m15.evaluate": _ast_guard(
            paths["first_switch_m15"],
            repo,
            "evaluate",
            (
                "append('first_switch_viscous_' + endpoint, p, t, [physical(1, 1), hbS if k == 0 else hbB], actual_Uy_source_endpoint=True, endpoint=endpoint)",
                "append('first_switch_viscous_upper_radius_error', zero, zeroi if empty else lead['radius_error'], [physical(1, 1), hbS, hbS], error_enclosure_only=True)",
            ),
            class_name="OriginalFirstSwitchM15",
        ),
        "macro_pressure_rows.valid_Z": _ast_guard(
            paths["macro_pressure_rows"],
            repo,
            "evaluate",
            (),
            class_name="OriginalMacroPressureRows",
            expressions=("-1 < c.mpf(Z) < 1",),
        ),
        "flat_sigma.sigma_jets": _ast_guard(
            paths["flat_sigma"],
            repo,
            "sigma_jets",
            (
                "if lo <= 0: rows.append(IntervalTaylor.constant(c, 0, 4))",
                "if hi >= 1: rows.append(IntervalTaylor.constant(c, 1, 4))",
                "x = c.mpf(x)",
                "return hull_jets(c, rows)",
            ),
            expressions=(
                "hi > 0 and lo <= mp.mpf('.5')",
                "lo < 1 and hi >= mp.mpf('.5')",
            ),
        ),
        "prescribed_source.first_switch_fields": _ast_guard(
            paths["prescribed_first_switch"],
            repo,
            "first_switch_fields",
            (
                "drive = full_drive(a, known, radius, F0, Pstar)",
                "chi_axial = g.mul(hb, g.sub(g.one, sigma(g, phase)))",
                "Uy = a.scale(a.mul(ratio, drive), g.neg(chi_axial))",
            ),
        ),
        "prescribed_source.first_switch_radius": _ast_guard(
            paths["prescribed_first_switch"],
            repo,
            "__init__",
            (),
            class_name="CurrentOriginalPrescribedR100RmC6Functions",
            expressions=(
                "self.native_parameter('hbS', 'original switch_owner.flow.h; not merged with hbB by name')",
                "g.mul(c(100), g.unary('exp', g.mul(hb, s)))",
            ),
        ),
        "prescribed_source.positive_hbS": _ast_guard(
            paths["prescribed_first_switch"],
            repo,
            "parameter_geometry_guard",
            (),
            class_name="CurrentOriginalPrescribedR100RmC6Functions",
            expressions=(
                "ep(shared.logh)",
                "ep(shared.logh)[1] < ep(c.ln(hcap))[0]",
            ),
        ),
        "cutoff_selector.transport_sigma": _ast_guard(
            repo / "src/openai_ns_reconstruction/original_first_collar_transport.py",
            repo,
            "__init__",
            (
                "self.sigma_jets = importlib.import_module(PREFIX + 'compliant_flat_pulse_derivatives').sigma_jets",
            ),
            class_name="OriginalFirstCollarTransport",
        ),
        "physical_ftc.pressure_density": _ast_guard(
            paths["physical_ftc"],
            repo,
            "pressure_density",
            (),
            expressions=("a.scale(Py, g.neg(R))",),
        ),
        "physical_ftc.pressure_FTC": _ast_guard(
            paths["physical_ftc"],
            repo,
            "pressure_FTC",
            (),
            expressions=("a.scale(a.add(Qy, Qyy), R)",),
        ),
    }
    for ast_name, ast_row in ast_rows.items():
        ast_path = ast_row.get("source_path")
        matches = [
            (pin_name, pin)
            for pin_name, pin in source_pins.items()
            if isinstance(pin, Mapping) and pin.get("path") == ast_path
        ]
        ast_row["source_pin_key"] = matches[0][0] if matches else None
        ast_row["source_pin_passed"] = bool(matches) and all(pin.get("passed") is True for _, pin in matches)
        if ast_row["source_pin_passed"] is not True:
            gaps.append(f"AST source pin failed: {ast_name}")
            ast_row["passed"] = False
    ast_hashes = {
        name: row.get("function_AST_sha256")
        for name, row in ast_rows.items()
        if isinstance(row, Mapping)
    }
    for name, row in ast_rows.items():
        if row.get("passed") is not True:
            gaps.append(f"AST source guard failed: {name}")

    # The parent views are executable 182-slot views.  We inspect only the
    # immutable names/digests around raw 43 and raw 45; no term value is read.
    views = parent.get("signed_m15_query_views")
    view_rows: list[dict[str, Any]] = []
    view_zs: list[str] = []
    if not isinstance(views, list) or not views:
        gaps.append("parent signed_m15_query_views are missing")
        views = []
    for view in views:
        if not isinstance(view, Mapping):
            gaps.append("malformed signed m15 view")
            continue
        z = str(view.get("Z"))
        view_zs.append(z)
        entries = view.get("entries")
        expected_indices = [index for index in range(183) if index != 39]
        index_list = [entry.get("parent_reduced_index") for entry in entries] if isinstance(entries, list) else []
        row: dict[str, Any] = {"Z": z, "view_terms": len(entries) if isinstance(entries, list) else 0}
        if index_list != expected_indices:
            gaps.append(f"signed view index chain changed at Z={z}")
        for raw_number, reduced_index, expected_name in (
            (43, 38, "first_switch_viscous_upper"),
            (45, 40, "first_switch_viscous_upper_radius_error"),
        ):
            matching = [entry for entry in entries if isinstance(entry, Mapping) and entry.get("parent_reduced_index") == reduced_index] if isinstance(entries, list) else []
            ok = bool(matching) and matching[0].get("name") == expected_name and matching[0].get("operation") == "unchanged_parent_reference"
            if not ok:
                gaps.append(f"raw {raw_number} mapping to reduced {reduced_index} is not retained at Z={z}")
            row[f"raw{raw_number}"] = {
                "raw_index": raw_number,
                "reduced_index": reduced_index,
                "name": expected_name,
                "active": ok,
                "operation": matching[0].get("operation") if matching else None,
                "parent_term_digest": matching[0].get("parent_term_digest") if matching else None,
            }
        view_rows.append(row)

    # Verify those parent term digests against the inherited 183-term packet.
    inherited_samples = {str(sample.get("Z")): sample for sample in inherited_report.get("samples", ()) if isinstance(sample, Mapping)}
    for row in view_rows:
        sample = inherited_samples.get(row["Z"])
        terms = sample.get("correlated_reduced_m15_C0_C1_terms") if isinstance(sample, Mapping) else None
        for raw_number, reduced_index in ((43, 38), (45, 40)):
            item = row.get(f"raw{raw_number}", {})
            actual_digest = None
            if isinstance(terms, list) and len(terms) == 183:
                actual_digest = _digest(terms[reduced_index])
            item["inherited_term_digest"] = actual_digest
            item["inherited_term_digest_matches"] = actual_digest is not None and actual_digest == item.get("parent_term_digest")
            if not item["inherited_term_digest_matches"]:
                gaps.append(f"inherited digest mismatch for raw {raw_number} at Z={row['Z']}")

    # Endpoint report coverage is retained as provenance.  The identity below
    # is source algebra, so no approximation or cap is inspected.
    pressure = report_values["first_switch_pressure_rows"]
    numeric = report_values["first_switch_numeric"]
    m15 = report_values["first_switch_m15"]
    endpoint_flags: list[dict[str, Any]] = []
    endpoint_source_structure = True
    ordinary_available: set[int] = set()
    for z in view_zs:
        p_rows = _samples_for_z(pressure, z, "1")
        n_rows = _samples_for_z(numeric, z, "1")
        m_rows = _samples_for_z(m15, z, "1")
        if len(p_rows) != 1:
            gaps.append(f"expected one first-switch pressure endpoint sample at Z={z}")
            endpoint_source_structure = False
        if len(n_rows) != 1:
            gaps.append(f"expected one first-switch numeric endpoint sample at Z={z}")
            endpoint_source_structure = False
        if len(m_rows) != 1:
            gaps.append(f"expected one first-switch m15 endpoint sample at Z={z}")
            endpoint_source_structure = False
        if not p_rows:
            continue
        sample = p_rows[0]
        leading = _mapping(sample.get("leading_U_y_factored"))
        source_cover = _mapping(sample.get("source_cover"))
        geometry = _mapping(sample.get("geometry"))
        contract = _mapping(sample.get("contract"))
        coefficient_rows = leading.get("coefficient_ordinary_Z")
        if isinstance(coefficient_rows, list):
            ordinary_available.update(index for index in range(min(3, len(coefficient_rows))))
        flag = leading.get("actual_unmodified_U_y_exact_zero_at_s1") is True
        structural = (
            flag
            and leading.get("actual_source") == "-hbS*(1-sigma(s))*(Phi/barPhi)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)"
            and leading.get("axial_control") == "1-sigma(s)"
            and leading.get("definition") == "U_y=epsilon*hbS*(1-sigma(s))*(original normalized current-radius drive with full angular error)"
            and source_cover.get("native_source_domain") == ["0", "1"]
            and contract.get("exact_control_zero_at_s1_not_misused_for_finite_N") is True
            and isinstance(coefficient_rows, list)
            and len(coefficient_rows) >= 2
            and geometry.get("R_definition") == "100*exp(hbS*s)"
            and geometry.get("native_measure") == "dy=hbS ds; R*dR=R^2*hbS ds"
            and geometry.get("original_hbS_is_not_hb_to_a_power") is True
            and geometry.get("defining_width_not_selected_from_cap") is True
        )
        if not flag:
            gaps.append(f"endpoint flag is not exact at Z={z}")
        if not structural:
            endpoint_source_structure = False
            endpoint_flags.append({"Z": z, "passed": False, "counterevidence": "endpoint metadata or ordinary source rows are incomplete"})
        else:
            endpoint_flags.append({"Z": z, "passed": True, "ordinary_rows": len(coefficient_rows), "source_values_selected": False, "numeric_caps_used": False})

    # The old m15 row names and FTC orientation are checked without treating
    # their coefficient payloads as source values.
    m15_endpoint_structure = True
    m15_contract_flags: list[dict[str, Any]] = []
    for z in view_zs:
        rows = _samples_for_z(m15, z, "1")
        if not rows:
            m15_endpoint_structure = False
            continue
        m15_contract = _mapping(rows[0].get("contract"))
        m15_contract_flags.append(
            {
                "Z": z,
                "both_finite_N_pressure_endpoints": m15_contract.get("both_finite_N_pressure_endpoints") is True,
                "complete_prefix_and_quadratic_errors": m15_contract.get("complete_prefix_and_quadratic_errors") is True,
                "original_selected_N_unchanged": m15_contract.get("original_selected_N_unchanged") is True,
                "exact_mean_ODE_used_for_viscous_increment": m15_contract.get("exact_mean_ODE_used_for_viscous_increment") is True,
                "inherited_exp_minus_Y_mean_mode_not_rounded_away": m15_contract.get("inherited_exp_minus_Y_mean_mode_not_rounded_away") is True,
            }
        )
        names = {str(term.get("name")): term for term in rows[0].get("m15_C0_C1_terms", ()) if isinstance(term, Mapping)}
        upper = _mapping(names.get("first_switch_viscous_upper"))
        radius = _mapping(names.get("first_switch_viscous_upper_radius_error"))
        upper_binding = _mapping(upper.get("actual_function_binding"))
        radius_binding = _mapping(radius.get("actual_function_binding"))
        good = (
            upper_binding.get("actual_Uy_source_endpoint") is True
            and upper_binding.get("endpoint") == "upper"
            and radius_binding.get("error_enclosure_only") is True
            and isinstance(upper.get("coefficient_C0_C1"), list)
            and len(upper.get("coefficient_C0_C1")) >= 2
            and isinstance(radius.get("coefficient_C0_C1"), list)
            and len(radius.get("coefficient_C0_C1")) >= 2
        )
        if not good:
            m15_endpoint_structure = False
            gaps.append(f"past m15 raw 43/raw 45 source categories are incomplete at Z={z}")

    ftc = _mapping(parent.get("signed_ledger_binding"))
    ftc_path = paths["physical_ftc"]
    ftc_evidence = {
        "passed": (
            ftc.get("passed") is True
            and ftc.get("independent_leading_FTC_passed") is True
            and ftc.get("mean_ODE_partition") == "Hlower-Hupper=Elower-Eupper+integral R*C(Uy) dy"
            and ftc.get("absolute_mean_retained_in_H") is True
            and ftc.get("actual_minus_proxy_category_reconstructed_from_original_physical_FTC") is True
            and ftc.get("signed_endpoint_orientation") == {"native_lower": 1, "macro_upper": -1, "first_switch_lower_E": 1}
            and ftc.get("old_raw_indices") == [23, 24, 44]
            and ftc.get("old_reduced_indices") == [19, 20, 39]
        ),
        "physical_density_FTC_source": ftc.get("physical_density_FTC_source"),
        "independent_leading_FTC_digest": ftc.get("independent_leading_FTC_digest"),
        "absolute_mean_retained_in_H": ftc.get("absolute_mean_retained_in_H") is True,
        "source_values_selected": False,
        "numeric_caps_used": False,
    }
    if not ftc_evidence["passed"]:
        gaps.append("inherited physical FTC source mapping is incomplete")

    ast_passed = bool(ast_rows) and all(row.get("passed") is True for row in ast_rows.values())
    ordinary_passed = endpoint_source_structure and m15_endpoint_structure and {0, 1}.issubset(ordinary_available)
    finite_n_not_zeroed = bool(
        m15_contract_flags
        and all(
            row["both_finite_N_pressure_endpoints"]
            and row["complete_prefix_and_quadratic_errors"]
            and row["original_selected_N_unchanged"]
            for row in m15_contract_flags
        )
    )
    mean_not_zeroed = bool(
        ftc_evidence["absolute_mean_retained_in_H"]
        and m15_contract_flags
        and all(
            row["exact_mean_ODE_used_for_viscous_increment"]
            and row["inherited_exp_minus_Y_mean_mode_not_rounded_away"]
            for row in m15_contract_flags
        )
    )
    source_family_passed = _family_ok(pressure) and _family_ok(numeric) and _family_ok(m15)
    if not source_family_passed:
        gaps.append("first-switch saved reports do not retain the active source family")

    truth_factor = ast_rows["pressure_rows.source_Uy"].get("passed") is True
    sigma_endpoint = ast_rows["flat_sigma.sigma_jets"].get("passed") is True
    sigma_interior_guards = bool(
        sigma_endpoint
        and ast_rows["flat_sigma.sigma_jets"].get("required_expressions", {}).get(
            "hi > 0 and lo <= mp.mpf('.5')"
        ) is True
        and ast_rows["flat_sigma.sigma_jets"].get("required_expressions", {}).get(
            "lo < 1 and hi >= mp.mpf('.5')"
        ) is True
    )
    sigma_selector = ast_rows["cutoff_selector.transport_sigma"].get("passed") is True
    exact_flag = ast_rows["pressure_rows.evaluate"].get("passed") is True
    source_body = ast_rows["prescribed_source.first_switch_fields"].get("passed") is True
    radius_geometry = ast_rows["prescribed_source.first_switch_radius"].get("passed") is True
    positive_hbS = ast_rows["prescribed_source.positive_hbS"].get("passed") is True
    z_domain = ast_rows["macro_pressure_rows.valid_Z"].get("passed") is True
    ordinary_extraction = ast_rows["first_switch_m15.leading"].get("passed") is True
    exact_source = bool(
        truth_factor
        and sigma_endpoint
        and sigma_selector
        and exact_flag
        and source_body
        and radius_geometry
        and positive_hbS
        and z_domain
        and ordinary_extraction
        and ordinary_passed
    )
    counterevidence: list[str] = []
    if not truth_factor:
        counterevidence.append("source_Uy truth jet is not visibly multiplied by kappa(s) before extraction")
    if not sigma_endpoint:
        counterevidence.append("sigma_jets does not expose the exact hi >= 1 -> sigma = 1 endpoint branch")
    if not sigma_selector:
        counterevidence.append("actual t.sigma_jets is not source-bound to the pinned flat sigma_jets implementation")
    if not sigma_interior_guards:
        counterevidence.append("sigma_jets interior branches are not source-bound for the exact x=1 exclusion")
    if not exact_flag:
        counterevidence.append("evaluate does not retain the exact s._mpi_ == 1 endpoint flag")
    if not source_body:
        counterevidence.append("original first-switch source body does not bind chi_axial = hbS*(1-sigma)")
    if not radius_geometry:
        counterevidence.append("prescribed first-switch source does not bind R(s)=100*exp(hbS*s)")
    if not positive_hbS:
        counterevidence.append("prescribed first-switch geometry does not expose the positive hbS guard")
    if not z_domain:
        counterevidence.append("valid Z domain -1 < Z < 1 is not source-bound")
    if not ordinary_extraction:
        counterevidence.append("m15 ordinary extraction does not consume the actual truth source")
    counterevidence.extend(
        f"saved endpoint counterevidence at Z={row['Z']}: {row.get('counterevidence')}"
        for row in endpoint_flags
        if row.get("passed") is not True
    )
    if counterevidence:
        gaps.extend(f"counterevidence: {item}" for item in counterevidence)

    mapping = {
        "raw43": {
            "raw_index": 43,
            "reduced_index": 38,
            "name": "first_switch_viscous_upper",
            "formula": "-100*epsilon*hbS*E(Bup)",
            "source_endpoint": "upper actual U_y(s=1,Z)",
            "active_in_parent_views": all(row.get("raw43", {}).get("active") for row in view_rows),
            "retirement_reason": "actual_same_Uy_source_zero_only",
            "positive_cap_value_zero_claim": False,
            "source_values_selected": False,
            "numeric_caps_used": False,
            "per_view": view_rows,
        },
        "raw45": {
            "raw_index": 45,
            "reduced_index": 40,
            "name": "first_switch_viscous_upper_radius_error",
            "kind": "actual_radius_correction",
            "formula": "-(R(s)-100)*epsilon*hbS*E(Bup)",
            "equivalent_formula": "100*epsilon*hbS^2*[-s*exprel(hbS*s)*E(Bup)]",
            "exprel": "exprel(t)=(exp(t)-1)/t",
            "radius_definition": "R(s)=100*exp(hbS*s)",
            "old_positive_cap_may_be_nonzero": True,
            "source_endpoint": "same upper actual U_y(s=1,Z) source atom",
            "active_in_parent_views": all(row.get("raw45", {}).get("active") for row in view_rows),
            "retirement_reason": "actual_same_Uy_source_zero_only",
            "positive_cap_value_zero_claim": False,
            "source_values_selected": False,
            "numeric_caps_used": False,
            "per_view": view_rows,
        },
    }

    sourceflagtruth = bool(truth_factor and ordinary_extraction and ordinary_passed)
    sigmaendpointbranch = bool(sigma_endpoint and sigma_selector and sigma_interior_guards and exact_flag)
    ordinaryZsourcezeroorders = bool(ordinary_passed and ordinary_extraction and z_domain)
    scopephase1 = bool(
        exact_source
        and ftc_evidence["passed"]
        and all(row.get("passed") is True for row in endpoint_flags)
    )
    capnotselected = bool(
        all(
            row.get("source_values_selected") is False and row.get("numeric_caps_used") is False
            for row in endpoint_flags
        )
        and ftc_evidence["source_values_selected"] is False
        and ftc_evidence["numeric_caps_used"] is False
        and all(
            item.get("source_values_selected") is False and item.get("numeric_caps_used") is False
            for item in mapping.values()
        )
    )
    finiteN_meannotzeroed = bool(finite_n_not_zeroed and mean_not_zeroed)
    compact_proofs_passed = all(
        (
            sourceflagtruth,
            sigmaendpointbranch,
            ordinaryZsourcezeroorders,
            scopephase1,
            capnotselected,
            finiteN_meannotzeroed,
        )
    )

    result = {
        "passed": bool(
            not gaps
            and parent_pin.get("passed") is True
            and source_family == _FAMILY
            and checks.get("passed") is True
            and ast_passed
            and source_family_passed
            and ftc_evidence["passed"]
            and exact_source
            and all(row.get("passed") is True for row in endpoint_flags)
            and compact_proofs_passed
        ),
        "source_family": dict(_FAMILY),
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": parent_pin,
        "parent_mapping_hash_bind": {
            "passed": parent_pin.get("passed") is True and parent_pin.get("mapping_match") is True,
            "supplied_mapping_sha256": parent_pin.get("supplied_mapping_sha256"),
            "latest_mapping_sha256": parent_pin.get("latest_mapping_sha256"),
        },
        "input_hashes": input_hashes,
        "implementation_hashes": implementation_hashes,
        "source_pins": source_pins,
        "AST": ast_rows,
        "ast_function_hashes": ast_hashes,
        "inherited_fields_digests": inherited_digests,
        "compact_proofs_passed": compact_proofs_passed,
        "source_family_passed": source_family_passed,
        "ordinary_Z_orders": [0, 1],
        "ordinary_orders": {
            "required": [0, 1],
            "available": sorted(ordinary_available),
            "passed": ordinary_passed,
            "source_function_order": 5,
            "source_values_selected": False,
            "numeric_caps_used": False,
        },
        "valid_Z_domain": {
            "expression": "-1 < Z < 1",
            "source_bound": z_domain,
            "sampled_Z": view_zs,
            "samples_are_coverage_only": True,
        },
        "endpoint_coverage": endpoint_flags,
        "exact_zero_identity": {
            "passed": exact_source,
            "source_atom": "E(B_up)",
            "identity": "U_y(1,Z) = epsilon*hbS*(1-sigma(1))*D(1,Z) = 0",
            "kappa": "1-sigma(s)",
            "kappa_endpoint": "sigma(1)=1 exactly from sigma_jets endpoint branch",
            "truth_source_multiplied_before_ordinary_Z": truth_factor and ordinary_extraction,
            "ordinary_Z_orders": [0, 1],
            "valid_Z_domain": "-1 < Z < 1",
            "epsilon_materialized": False,
            "width_materialized": False,
            "proof_is_source_algebra": True,
            "not_from_numeric_rows_caps_midpoints_underflow_or_interval_overlap": True,
            "exact_flag_only_point_branch": not truth_factor,
            "counterevidence": counterevidence,
        },
        "physical_ftc": ftc_evidence,
        "mapping": mapping,
        # These compact keys are the strict packet boundary.  They carry
        # source identities and negative claims explicitly, without exposing
        # any numerical field value as the proof.
        "sourceflagtruth": {
            "passed": sourceflagtruth,
            "truth_source_multiplied_before_ordinary_Z": sourceflagtruth,
            "kappa": "1-sigma(s)",
            "actual_unmodified_truth_source": True,
            "nonlinear_error_norm_zero_required": False,
            "parts_are_multiplied_by_kappa": truth_factor,
        },
        "sigmaendpointbranch": {
            "passed": sigmaendpointbranch,
            "identity": "sigma(1)=1",
            "exact_branch": "if hi >= 1: sigma=1",
            "source_function": "sigma_jets",
            "actual_t_sigma_jets_bound_to_same_source": sigma_selector,
            "exact_x1_point_branch": exact_flag,
            "nonconstant_interior_guards_source_bound": sigma_interior_guards,
            "nonconstant_row_enters_hull_at_x1": False,
            "interval_straddle_used": False,
        },
        "ordinaryZsourcezeroorders": {
            "passed": ordinaryZsourcezeroorders,
            "orders": [0, 1],
            "source_function_order": 5,
            "valid_Z_domain": "-1 < Z < 1",
            "extraction_uses_actual_truth_source": ordinary_extraction,
            "zero_is_source_atom_identity": True,
        },
        "scopephase1": {
            "passed": scopephase1,
            "phase": "s=1",
            "source_atom": "E(B_up)",
            "scope": "Efirst only; H and actual mean retained",
            "radius_geometry": "R(s)=100*exp(hbS*s)",
            "positive_hbS_unmaterialized": positive_hbS,
        },
        "capnotselected": {
            "passed": capnotselected,
            "raw43_positive_cap_value_zero_claim": False,
            "raw45_old_cap_may_be_nonzero": True,
            "source_values_selected": False,
            "numeric_caps_used": False,
        },
        "finiteN_meannotzeroed": {
            "passed": finiteN_meannotzeroed,
            "finiteN": finite_n_not_zeroed,
            "meannotzeroed": mean_not_zeroed,
            "finite_N_corrections_retained": finite_n_not_zeroed,
            "actual_mean_zero_claim": False,
            "H_field_zero_claim": False,
            "own_mean_retained": True,
            "m15_contracts": m15_contract_flags,
        },
        "zero_scope": {
            "E_source_zero": exact_source,
            "actual_mean_zero_claim": False,
            "H_field_zero_claim": False,
            "own_mean_retained": True,
            "statement": "Only the source E(B_up) atom is zeroed; no actual mean or H field is claimed zero.",
        },
        "source_executed": False,
        "numeric_field_execution": False,
        "accepted_fixture_replayed": False,
        "owners_or_producers_imported": False,
        "gaps": gaps,
    }
    return result


__all__ = ["build_first_switch_zero_source"]
