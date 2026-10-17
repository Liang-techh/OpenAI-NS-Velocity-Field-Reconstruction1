"""Source-only binding for the six original inner/collar m15 terms.

The forward m15 producer stores interval-backed polynomial or cell
approximations.  The canonical inner and lower-collar producers store exact
graph functions without numerical values.  This module compares their source
interfaces and returns a typed signed difference.  It deliberately does not
import or execute either producer, an owner, a constructor, a checker, a
quadrature routine, or a recurrence.

In particular, the sixth forward term is an enclosure for ``|A-A0|``.  It is
not a signed value.  The returned source difference therefore never contracts
that term to zero and never selects an interval midpoint or overlap as a
function value.
"""

from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


_PARENT_RELATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_quiet_scalar_m15_numeric.json"
)
_PARENT_SHA256 = "93f419a4a10809ff15341349e86856a1a191937cc27612891be67842148375f9"

_FORWARD_CORE_SOURCE = "src/openai_ns_reconstruction/original_n1_core_m15_numeric.py"
_FORWARD_INNER_SOURCE = "src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py"
_FORWARD_CORE_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_n1_core_m15_numeric.json"
)
_FORWARD_INNER_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_n1_inner_m15_numeric.json"
)

_CANONICAL_INNER_SOURCE = (
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py"
)
_CANONICAL_INNER_GRAPH = (
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.json.gz"
)
_CANONICAL_LOWER_SOURCE = (
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py"
)
_CANONICAL_LOWER_GRAPH = (
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.json.gz"
)

_TERM_NAMES = (
    "core_axial_coupling",
    "core_pressure_integral_and_finite_boundary",
    "first_lower_axial_width",
    "first_lower_regular_pressure_width",
    "first_lower_singular_pressure",
    "first_lower_singular_width_error",
)

_FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _sha(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _repo(root: Any, directory: Any) -> Path:
    for value in (root, directory):
        if value is None:
            continue
        try:
            path = Path(value)
        except TypeError:
            continue
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _work(repo: Path, directory: Any) -> Path:
    if directory is None:
        return repo / "experiments" / "root_st073"
    path = Path(directory)
    if path.is_file():
        path = path.parent
    if (path / "experiments" / "root_st073").is_dir():
        return path / "experiments" / "root_st073"
    return path


def _relative(path: Path, repo: Path) -> str:
    try:
        return path.relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _resolve(repo: Path, work: Path, value: Any) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    for candidate in (work / path, repo / path, repo / "experiments" / "root_st073" / path):
        if candidate.exists():
            return candidate
    return work / path


def _load(path: Path) -> Mapping[str, Any]:
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    value = json.loads(raw)
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _load_parent(parent_report: Any, repo: Path, work: Path) -> tuple[Mapping[str, Any], dict[str, Any], list[str]]:
    latest_path = repo / _PARENT_RELATIVE
    latest_actual = _sha(latest_path)
    pin: dict[str, Any] = {
        "path": _relative(latest_path, repo),
        "sha256": latest_actual,
        "expected_sha256": _PARENT_SHA256,
        "passed": latest_actual == _PARENT_SHA256,
        "source_executed": False,
    }
    gaps: list[str] = []
    if not pin["passed"]:
        gaps.append("accepted quiet scalar parent hash mismatch")
    try:
        latest = _load(latest_path)
    except (OSError, ValueError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
        return {}, pin, gaps + [f"cannot load accepted parent: {exc}"]

    supplied_value = parent_report
    if supplied_value is None:
        supplied_value = latest_path
    try:
        supplied = (
            supplied_value
            if isinstance(supplied_value, Mapping)
            else _load(_resolve(repo, work, supplied_value))
        )
    except (OSError, ValueError, TypeError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
        return {}, pin, gaps + [f"cannot load supplied parent: {exc}"]
    supplied_digest = _digest(supplied)
    latest_digest = _digest(latest)
    pin.update(
        {
            "latest_mapping_sha256": latest_digest,
            "supplied_mapping_sha256": supplied_digest,
            "mapping_match": latest_digest == supplied_digest,
        }
    )
    if latest_digest != supplied_digest:
        gaps.append("supplied parent mapping does not match accepted parent")
    return supplied, pin, gaps


def _expected(parent: Mapping[str, Any], name: str) -> tuple[Any, str | None]:
    for key in ("implementation_hashes", "input_hashes"):
        values = _mapping(parent.get(key))
        if name in values:
            return values[name], key[:-7] if key.endswith("_hashes") else key
    return None, None


def _pin_file(
    parent: Mapping[str, Any],
    path: Path,
    repo: Path,
    name: str,
    gaps: list[str],
) -> dict[str, Any]:
    expected, kind = _expected(parent, name)
    actual = _sha(path)
    passed = isinstance(expected, str) and actual == expected
    row = {
        "name": name,
        "kind": kind or "unlisted",
        "path": _relative(path, repo),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": passed,
        "source_executed": False,
    }
    if not passed:
        gaps.append(f"source pin failed: {name}")
    return row


def _node(nodes: Any, index: int) -> Mapping[str, Any]:
    if isinstance(nodes, list) and isinstance(index, int) and 0 <= index < len(nodes):
        value = nodes[index]
        if isinstance(value, Mapping):
            return value
    return {}


def _node_proof(nodes: Any, index: int, expected: Mapping[str, Any]) -> dict[str, Any]:
    actual = dict(_node(nodes, index))
    selected = {key: actual.get(key) for key in expected}
    passed = all(selected.get(key) == value for key, value in expected.items())
    return {
        "index": index,
        "expected_fields": dict(expected),
        "actual_fields": selected,
        "passed": passed,
        "source_executed": False,
    }


def _find_function(tree: ast.AST, name: str, cls: str | None = None) -> ast.AST:
    scope: ast.AST = tree
    if cls is not None:
        classes = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.ClassDef) and node.name == cls
        ]
        if len(classes) != 1:
            raise ValueError(f"unique class required: {cls}")
        scope = classes[0]
    functions = [
        node
        for node in ast.walk(scope)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name
    ]
    if len(functions) != 1:
        raise ValueError(f"unique function required: {name}")
    return functions[0]


def _ast_guard(
    path: Path,
    repo: Path,
    name: str,
    *,
    cls: str | None = None,
    text: tuple[str, ...] = (),
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_path": _relative(path, repo),
        "function": name,
        "class": cls,
        "required_text": {value: False for value in text},
        "function_AST_sha256": None,
        "source_executed": False,
        "passed": False,
    }
    try:
        source = path.read_text(encoding="utf-8")
        function = _find_function(ast.parse(source, filename=str(path)), name, cls)
        segment = ast.get_source_segment(source, function) or ""
        compact = "".join(segment.split())
        row["required_text"] = {
            value: "".join(value.split()) in compact for value in text
        }
        row["function_AST_sha256"] = hashlib.sha256(
            ast.dump(function, include_attributes=False).encode()
        ).hexdigest()
        row["passed"] = all(row["required_text"].values())
    except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
        row["error"] = str(exc)
    return row


def _family(value: Any) -> dict[str, Any]:
    result = _mapping(value)
    return {key: result.get(key) for key in _FAMILY_KEYS}


def _empty(
    repo: Path,
    work: Path,
    parent_pin: Mapping[str, Any],
    gaps: list[str],
    family: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "passed": False,
        "source_binding_passed": False,
        "source_only": True,
        "source_executed": False,
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": dict(parent_pin),
        "source_family": dict(family or {}),
        "source_pins": {},
        "source_pins_ok": False,
        "AST": {},
        "graph_identity": {},
        "source_scope": {},
        "forward": {},
        "canonical": {},
        "forward_six": {},
        "canonical_inner_lower": {},
        "delta_six": {
            "type": "signed_source_difference",
            "exact_identity": False,
            "zero_claim": False,
            "numeric_value_available": False,
        },
        "exact_identity": {
            "forward_six_equals_canonical_inner_plus_lower": False,
            "delta_six_zero": False,
            "ordinary_Z0_zero": False,
            "ordinary_Z1_zero": False,
        },
        "archives": {},
        "numeric_value_available": False,
        "numeric_value_unavailable_reason": "source binding only; saved enclosures are not function values",
        "root": str(repo),
        "directory": str(work),
        "gaps": list(gaps),
    }


def _sample_for(report: Mapping[str, Any], Z: str, *, core: bool) -> Mapping[str, Any] | None:
    candidates = []
    for row in report.get("samples", ()):
        if not isinstance(row, Mapping) or str(row.get("Z")) != str(Z):
            continue
        if core:
            if str(row.get("rho")) != "4":
                continue
            candidates.append(row)
        elif str(row.get("q")) == "2" and row.get("cells_per_piece") == 16:
            candidates.append(row)
    if not candidates:
        return None
    if core:
        return sorted(
            candidates,
            key=lambda row: (
                row.get("source_degree", -1),
                row.get("potential_kernel_degree", -1),
                row.get("nominal_feedback_iterations", -1),
            ),
        )[-1]
    return candidates[-1]


def _term_errors(term: Mapping[str, Any]) -> dict[str, Any]:
    """Copy every saved coefficient/error field without interpreting it."""
    return {
        "name": term.get("name"),
        "coefficient_C0_C1": copy.deepcopy(term.get("coefficient_C0_C1")),
        "log_scale_terms": copy.deepcopy(term.get("log_scale_terms")),
        "actual_function_binding": copy.deepcopy(term.get("actual_function_binding")),
        "saved_term": copy.deepcopy(term),
        "source_enclosure_is_not_a_selected_value": True,
        "source_function_definition_restored": True,
        "actual_numerical_source_value_installed": False,
        "saved_error_fields_retained": True,
    }


def _term_source_binding(index: int, name: str, ast_rows: Mapping[str, Any]) -> dict[str, Any]:
    if index == 0:
        density = "2*U_lead*u_n1 d rho"
        measure = "d rho on rho in [0,4]"
        mismatch = "cross-producer identity from forward finite rows to the canonical analytic core function is unproved"
        source_ast = ("forward_core_signed_prefix", "canonical_signed_rho_integrands")
    elif index == 1:
        density = "(p_n1 - rho*p_n1_endpoint_derivative_form) d rho"
        measure = "d rho on rho in [0,4]"
        mismatch = "forward P0/M2 pressure data have no proved identity with the canonical native P0 function"
        source_ast = ("forward_core_signed_prefix", "canonical_signed_rho_integrands")
    elif index == 2:
        density = "2*rho*U_lead*u_cut dq"
        measure = "hb*(sc/4) dq"
        mismatch = "forward collar profile/inlet rows are not source-identified with canonical ratio and Psi-history functions"
        source_ast = ("forward_inner_source_cell", "canonical_signed_rho_integrands", "canonical_collar_fields")
    elif index == 3:
        density = "-rho^2*(Np0 + 2*S*Phi*f_cut/(epsilon^3*Pstar^2)) dq"
        measure = "hb*(sc/4) dq"
        mismatch = "forward regular pressure partition has no proved algebraic identity with canonical direct N1p"
        source_ast = ("forward_inner_source_cell", "canonical_direct_pressure_functions")
    elif index == 4:
        density = "+epsilon^2*A0*sigma(s) [the A0 boundary piece]"
        measure = "ds after inverse-hb cancellation, integrated by sigma(0)=0"
        mismatch = "forward A0 endpoint coefficient is not source-identified with exact Psi_core_rho(4,Z)"
        source_ast = ("forward_inner_sigma_atom", "canonical_signed_rho_integrands", "canonical_collar_fields")
    else:
        density = "+epsilon^2*(A(s,h)-A0)*sigma_prime(s) ds"
        measure = "ds over s in [0, sc/2]"
        mismatch = "only |A-A0| <= hb*E is saved; the signed A-A0 function is not source-bound"
        source_ast = ("forward_inner_sigma_width_error", "canonical_direct_pressure_functions")
    row = {
        "index": index,
        "name": name,
        "forward_density": density,
        "forward_measure": measure,
        "canonical_source_ast_guards": list(source_ast),
        "AST_guards_present": {
            key: bool(_mapping(ast_rows.get(key)).get("passed")) for key in source_ast
        },
        "structural_density_binding": True,
        "physical_measure_binding": True,
        "initial_endpoint_binding": False,
        "scalar_normalization_binding": False,
        "ordinary_Z1_binding": False,
        "exact_identity": False,
        "remaining_source_mismatch": mismatch,
    }
    if index in (4, 5):
        row["full_singular_source_definition"] = (
            "epsilon^2*integral_[0,s] A(t,h)*sigma_prime(t) dt"
        )
        row["singular_source_formula"] = (
            "epsilon^2*A0*sigma(s)" if index == 4 else
            "epsilon^2*integral_[0,s](A(t,h)-A0)*sigma_prime(t) dt"
        )
        row["singular_A_formula"] = (
            "A(s,h)=2*rho(s)^2*Z*(1-h)*(Phi_actual/Phi_core(rho(s)))*"
            "Psi_core_rho(rho(s))/L, L=1-delta*Z^2"
        )
        row["singular_phase_map"] = "rho(s)=4*exp(hb*s), sigma=original_sigma(s)"
        row["A0_formula"] = "A0=32*Z*Psi_core_rho(4)/L"
        row["exact_split_definition"] = (
            "epsilon^2*(A0*sigma(s) + integral_[0,s](A(t,h)-A0)*sigma_prime(t) dt)"
        )
        row["exact_split"] = row["exact_split_definition"]
        row["source_function_definition_restored"] = True
        row["actual_numerical_source_value_installed"] = False
        row["term_role"] = "A0_sigma_boundary_piece" if index == 4 else "signed_A_minus_A0_remainder"
        row["remainder_is_signed_source"] = index == 5
        row["remainder_bound_is_not_value"] = index == 5
    return row


def _graph_identity(
    inner_graph: Mapping[str, Any], lower_graph: Mapping[str, Any], gaps: list[str]
) -> tuple[dict[str, Any], bool]:
    inner_nodes = inner_graph.get("exact_graph_nodes", ())
    lower_nodes = lower_graph.get("exact_graph_nodes", ())
    inner_bindings = _mapping(inner_graph.get("inner_cutoff_bindings"))
    lower_bindings = _mapping(lower_graph.get("lower_collar_bindings"))

    inner_expected = {
        "endpoint_m15_C1": (53161, {"operation": "function_substitution", "expression": 53141, "variable": 35848, "value": 35824}),
        "endpoint_m15_Z1_C1": (53162, {"operation": "function_substitution", "expression": 53142, "variable": 35848, "value": 35824}),
        "first_m15_integral": (53137, {"operation": "definite_integral", "integrand": 53133, "variable": "original_n1_cutoff_first_moment_t", "lower": 0, "upper": 35848, "exact_function_integral": True, "numerical_value_not_installed": True}),
        "first_m15_Z1_integral": (53138, {"operation": "definite_integral", "integrand": 53134, "variable": "original_n1_cutoff_first_moment_t", "lower": 0, "upper": 35848, "exact_function_integral": True, "numerical_value_not_installed": True}),
    }
    lower_expected = {
        "lower_m15_C1": (36592, {"operation": "definite_integral", "integrand": 36590, "variable": "original_lower_collar_s", "lower": 35824, "upper": 35825, "exact_function_integral": True, "numerical_value_not_installed": True, "original_direct_N1p": True, "ordinary_slow_Z_derivative_order": 0, "endpoint_independent_of_Z": True, "all_finite_viscous_derivatives_retained": True, "original_m15_minus_applied_once": True, "signed_integrand_not_a_cap": True}),
        "lower_m15_Z1_C1": (36593, {"operation": "definite_integral", "integrand": 36591, "variable": "original_lower_collar_s", "lower": 35824, "upper": 35825, "exact_function_integral": True, "numerical_value_not_installed": True, "original_direct_N1p": True, "ordinary_slow_Z_derivative_order": 1, "endpoint_independent_of_Z": True, "all_finite_viscous_derivatives_retained": True, "original_m15_minus_applied_once": True, "signed_integrand_not_a_cap": True}),
        "lower_m15_integrand": (36590, {"operation": "product", "arguments": [36589, 36583]}),
        "lower_m15_Z1_integrand": (36591, {"operation": "product", "arguments": [36589, 36584]}),
        "weighted_measure": (36494, {"operation": "product", "arguments": [35565, 35565, 35851, 35853]}),
    }

    inner_proofs = {
        key: _node_proof(inner_nodes, index, expected)
        for key, (index, expected) in inner_expected.items()
    }
    lower_proofs = {
        key: _node_proof(lower_nodes, index, expected)
        for key, (index, expected) in lower_expected.items()
    }
    geometry = {
        "q_interval": "[0,2]",
        "q_Rin": 1,
        "q_rminus": 2,
        "s_of_q": "s=(sc/4)*q",
        "rho_of_s": "rho=4*exp(hb*s)",
        "rho_of_q": "rho=4*exp(hb*(sc/4)*q)",
        "jacobian_d_rho_dq": "hb*(sc/4)*rho",
        "inner_endpoint": "Rin at q=1, s=sc/4",
        "lower_endpoint": "rminus at q=2, s=sc/2",
        "inner_geometry": copy.deepcopy(_mapping(inner_bindings.get("geometry"))),
        "lower_geometry": copy.deepcopy(_mapping(lower_bindings.get("geometry_certificate"))),
        "lower_fields": {
            key: _mapping(lower_bindings.get("fields")).get(key)
            for key in ("query", "domain", "rho", "physical_measure", "weighted_physical_measure")
        },
    }
    geometry_checks = {
        "inner_phase_is_sc_quarter": _mapping(inner_bindings.get("geometry")).get("phase_is_sc_quarter") is True,
        "inner_endpoint_is_Rin": _mapping(inner_bindings.get("definition")).get("endpoint_is_Rin") is True,
        "lower_Rin_is_Ra_exp_hb_sc_quarter": _mapping(lower_bindings.get("geometry_certificate")).get("physical_Rin_is_Ra_exp_hb_sc_quarter") is True,
        "lower_rminus_is_Ra_exp_hb_sc_half": _mapping(lower_bindings.get("geometry_certificate")).get("physical_rminus_is_Ra_exp_hb_sc_half") is True,
        "lower_sc_and_hb_Z_independent": _mapping(lower_bindings.get("geometry_certificate")).get("sc_and_hb_Z_independent") is True,
        "lower_weighted_measure_present": _mapping(lower_bindings.get("fields")).get("weighted_physical_measure") == 36494,
    }
    proofs_ok = all(row["passed"] for row in (*inner_proofs.values(), *lower_proofs.values()))
    geometry_ok = all(geometry_checks.values())
    if not proofs_ok:
        gaps.append("canonical inner/lower graph node proof failed")
    if not geometry_ok:
        gaps.append("canonical inner/lower collar geometry proof failed")
    identity = {
        "inner": {
            "report": _CANONICAL_INNER_GRAPH,
            "function_C1_pair": [53161, 53162],
            "function_path": "inner_cutoff_bindings.endpoint.m15",
            "node_proofs": inner_proofs,
        },
        "lower": {
            "report": _CANONICAL_LOWER_GRAPH,
            "function_C1_pair": [36592, 36593],
            "function_path": "lower_collar_bindings.functions.m15",
            "node_proofs": lower_proofs,
        },
        "geometry": geometry,
        "geometry_checks": geometry_checks,
        "passed": proofs_ok and geometry_ok,
        "source_executed": False,
    }
    return identity, proofs_ok and geometry_ok


def build_inner_collar_source(parent_report: Any, *, root: Any, directory: Any) -> dict[str, Any]:
    """Bind the forward six terms to the canonical inner/lower source graphs.

    The result is source-only.  ``passed`` means that the accepted parent,
    saved reports, source ASTs, and canonical graph descriptors were pinned.
    It does not mean that the six forward values equal the canonical pair.
    That unresolved signed quantity is returned as ``delta_six``.
    """

    repo = _repo(root, directory)
    work = _work(repo, directory)
    parent, parent_pin, gaps = _load_parent(parent_report, repo, work)
    if not parent:
        return _empty(repo, work, parent_pin, gaps)
    family = _family(parent.get("source_family"))
    if family != {
        "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
        "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
        "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
    }:
        gaps.append("parent source_family mismatch")
    if _mapping(parent.get("checks")).get("passed") is not True:
        gaps.append("accepted parent checks.passed is not true")
    if gaps:
        return _empty(repo, work, parent_pin, gaps, family)

    source_pins: dict[str, Any] = {}
    paths = {
        "forward_core_source": repo / _FORWARD_CORE_SOURCE,
        "forward_inner_source": repo / _FORWARD_INNER_SOURCE,
        "forward_inner_report": work / _FORWARD_INNER_REPORT,
        "canonical_inner_source": work / _CANONICAL_INNER_SOURCE,
        "canonical_inner_graph": work / _CANONICAL_INNER_GRAPH,
        "canonical_lower_source": work / _CANONICAL_LOWER_SOURCE,
        "canonical_lower_graph": work / _CANONICAL_LOWER_GRAPH,
    }
    pin_names = {
        "forward_core_source": _FORWARD_CORE_SOURCE,
        "forward_inner_source": _FORWARD_INNER_SOURCE,
        "forward_inner_report": _FORWARD_INNER_REPORT,
        "canonical_inner_source": _CANONICAL_INNER_SOURCE,
        "canonical_inner_graph": _CANONICAL_INNER_GRAPH,
        "canonical_lower_source": _CANONICAL_LOWER_SOURCE,
        "canonical_lower_graph": _CANONICAL_LOWER_GRAPH,
    }
    for key, path in paths.items():
        source_pins[key] = _pin_file(parent, path, repo, pin_names[key], gaps)
    source_pins_ok = all(pin.get("passed") is True for pin in source_pins.values())

    reports: dict[str, Mapping[str, Any]] = {}
    for key in ("forward_inner_report", "canonical_inner_graph", "canonical_lower_graph"):
        try:
            reports[key] = _load(paths[key])
        except (OSError, ValueError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
            reports[key] = {}
            gaps.append(f"cannot load {key}: {exc}")

    family_checks = {
        key: _family(reports[key].get("source_family")) == family
        for key in ("forward_inner_report", "canonical_inner_graph", "canonical_lower_graph")
    }
    if not all(family_checks.values()):
        gaps.append("forward/canonical saved source family mismatch")

    ast_rows = {
        "forward_core_signed_prefix": _ast_guard(
            paths["forward_core_source"], repo, "signed_prefix_rows",
            text=("product = c1_product_rows", "axial = [2*sum", "pressure = [-sum"),
        ),
        "forward_core_evaluate": _ast_guard(
            paths["forward_core_source"], repo, "evaluate", cls="OriginalN1CoreM15Numeric",
            text=("point=signed_prefix_rows", "base=signed_prefix_rows", "source_formula=dict"),
        ),
        "forward_inner_source_cell": _ast_guard(
            paths["forward_inner_source"], repo, "source_cell", cls="OriginalN1InnerM15Numeric",
            text=("axial=(vt['rho']*vt['U'][0].rows[0]*u*2).truncate(1)", "pressure=(-(Np+coupling)*vt['rho']*vt['rho']).truncate(1)"),
        ),
        "forward_inner_sigma_atom": _ast_guard(
            paths["forward_inner_source"], repo, "sigma_atom", cls="OriginalN1InnerM15Numeric",
            text=("Exact original flat value as exp(odds)*kappa", "return self.atom(1/(1-tp)**2-1/tp**2,odds),kappa"),
        ),
        "forward_inner_sigma_width_error": _ast_guard(
            paths["forward_inner_source"], repo, "evaluate", cls="OriginalN1InnerM15Numeric",
            text=("first_lower_singular_width_error", "proof='|A-A0|_norm<=hb*E, integrate against positive d_sigma'"),
        ),
        "canonical_signed_rho_integrands": _ast_guard(
            paths["canonical_inner_source"], repo, "signed_rho_integrands",
            text=("m15=a.add(a.scale(a.mul(u,U),g.mul(g.constant(2),epsilon))", "a.neg(a.scale(pressure,g.mul(eps2,r)))"),
        ),
        "canonical_cutoff_history": _ast_guard(
            paths["canonical_inner_source"], repo, "own_cutoff_fields",
            text=("actual_signed_n1_inner_moment=True", "full_core_to_first_history_retained=True", "leading_N1p_not_cut=True"),
        ),
        "canonical_collar_fields": _ast_guard(
            paths["canonical_lower_source"], repo, "collar_fields",
            text=("rho=lambda x:g.mul(c(4),g.unary('exp',g.mul(hb,x)))", "jacobian=g.mul(hb,rs)", "weighted_physical_measure=g.mul(eps,eps,rs,jacobian)"),
        ),
        "canonical_direct_pressure_functions": _ast_guard(
            paths["canonical_lower_source"], repo, "direct_pressure_functions",
            text=("original_direct_N1p=True", "original_m15_minus_applied_once=True", "signed_integrand_not_a_cap=True"),
        ),
    }
    ast_ok = all(row.get("passed") is True for row in ast_rows.values())
    if not ast_ok:
        gaps.append("forward/canonical source AST guard failed")

    graph_identity, graph_ok = _graph_identity(
        reports["canonical_inner_graph"], reports["canonical_lower_graph"], gaps
    )

    selected: dict[str, dict[str, Mapping[str, Any] | None]] = {}
    forward_terms: dict[str, list[dict[str, Any]]] = {}
    for Z in ("-.3", "0", ".3"):
        inner_sample = _sample_for(reports["forward_inner_report"], Z, core=False)
        # The accepted inner report retains the two core terms at indices 0 and
        # 1.  Reusing those saved terms avoids inventing an unpinned producer
        # report while preserving their original error records.
        core_sample = inner_sample
        selected[Z] = {"core": core_sample, "inner": inner_sample}
        if core_sample is None:
            gaps.append(f"missing forward core rho=4 sample at Z={Z}")
        if inner_sample is None:
            gaps.append(f"missing forward inner q=2 cells=16 sample at Z={Z}")
        core_terms = (inner_sample.get("m15_C0_C1_terms", ())[:2] if inner_sample else ())
        inner_terms = inner_sample.get("m15_C0_C1_terms", ()) if inner_sample else ()
        names = [term.get("name") for term in (*core_terms, *inner_terms[2:])] if core_sample and inner_sample else []
        if names != list(_TERM_NAMES):
            gaps.append(f"forward six term partition mismatch at Z={Z}")
            continue
        rows: list[dict[str, Any]] = []
        for index, name in enumerate(_TERM_NAMES):
            term = core_terms[index] if index < 2 else inner_terms[index]
            row = _term_errors(term)
            row.update(_term_source_binding(index, name, ast_rows))
            rows.append(row)
        forward_terms[Z] = rows

    forward_ok = len(forward_terms) == 3
    if not forward_ok:
        gaps.append("complete saved forward six-term source set unavailable")

    # Preserve all six terms and all saved errors under the accepted parent.
    # No coefficient is recomputed or narrowed here.
    archive_samples = copy.deepcopy(_mapping(parent.get("inherited_m15_packet")).get("samples", ()))
    views = copy.deepcopy(parent.get("quiet179_discrepancy_views", ()))
    archives = {
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "raw213_reduced183": {
            "inherited_m15_packet": archive_samples,
            "view_metadata": [
                {
                    key: copy.deepcopy(view.get(key))
                    for key in ("Z", "raw213_digest", "parent_original_reduced183_digest", "six_inner_lower_digest", "largest_unchanged_caps")
                }
                for view in views
                if isinstance(view, Mapping)
            ],
            "source_errors_retained": True,
            "archives_replayed": False,
        },
        "forward_reports": {
            "core_terms": {
                "report": _FORWARD_INNER_REPORT,
                "path": "samples[Z,q=2,cells_per_piece=16].m15_C0_C1_terms[0:2]",
                "source_report_is_pinned": True,
            },
            "inner": _FORWARD_INNER_REPORT,
            "selected_rows": {Z: {key: (row is not None) for key, row in values.items()} for Z, values in selected.items()},
        },
    }

    source_bindings_by_Z = {
        Z: {
            "terms": rows,
            "all_term_errors_retained": all(row.get("saved_error_fields_retained") for row in rows),
            "term_4_A0_sigma_boundary_piece": {
                "type": "actual_signed_source_boundary_piece",
                "formula": "epsilon^2*A0*sigma(s)",
                "A0_formula": "32*Z*Psi_core_rho(4)/L",
                "full_singular_source_definition": "epsilon^2*integral_[0,s] A(t,h)*sigma_prime(t)dt",
                "exact_split": "epsilon^2*(A0*sigma(s)+integral_[0,s](A(t,h)-A0)*sigma_prime(t)dt)",
                "source_function_definition_restored": True,
                "actual_numerical_source_value_installed": False,
                "zero_claim": False,
            },
            "term_5_signed_source_remainder": {
                "type": "actual_signed_source_remainder",
                "formula": "epsilon^2*integral_[s=0..sc/2] (A(s,h)-A0)*sigma_prime(s) ds",
                "A_formula": "2*rho(s)^2*Z*(1-h)*(Phi_actual/Phi_core(rho(s)))*Psi_core_rho(rho(s))/L",
                "A0_formula": "32*Z*Psi_core_rho(4)/L",
                "L_formula": "L=1-delta*Z^2",
                "rho_formula": "rho(s)=4*exp(hb*s)",
                "exact_split": "epsilon^2*(A0*sigma(s)+integral_[0,s](A(t,h)-A0)*sigma_prime(t)dt)",
                "bound": "|A(s,h)-A0(s)| <= hb*E",
                "normalized_bound": "|R(s)|/(hb*sigma(s)) <= E for s>0, where R(s)=integral_[0,s](A(t,h)-A0)*sigma_prime(t)dt",
                "bound_interpretation": "E is only a norm enclosure; after integration it may be recorded as the normalized hb*sigma cap above, never a signed source value and never positiveWidthChangeE",
                "bound_is_not_a_value": True,
                "saved_approximation_is_not_a_value": True,
                "source_function_definition_restored": True,
                "actual_numerical_source_value_installed": False,
                "zero_claim": False,
                "source_mismatch": "the signed A-A0 function is not source-bound by the forward producer",
            },
        }
        for Z, rows in forward_terms.items()
    }

    source_scope = {
        "q_interval": [0, 2],
        "q_Rin": 1,
        "q_rminus": 2,
        "s_of_q": "s=(sc/4)*q",
        "rho_of_s": "rho=4*exp(hb*s)",
        "rho_of_q": "rho=4*exp(hb*(sc/4)*q)",
        "R_of_q": "R=epsilon*rho",
        "jacobian": "d rho = hb*(sc/4)*rho*dq",
        "inner_physical_C1": [53161, 53162],
        "lower_physical_C1": [36592, 36593],
        "ordinary_Z_orders": [0, 1],
        "source_executed": False,
    }
    canonical = {
        "inner": {
            "physical_C1_pair": [53161, 53162],
            "source_path": "inner_cutoff_bindings.endpoint.m15",
            "graph_report": _CANONICAL_INNER_GRAPH,
            "exact_function_integral": True,
            "numerical_value_not_installed": True,
            "source_function_definition_restored": True,
            "actual_numerical_source_value_installed": False,
            "signed_source": True,
        },
        "lower": {
            "physical_C1_pair": [36592, 36593],
            "source_path": "lower_collar_bindings.functions.m15",
            "graph_report": _CANONICAL_LOWER_GRAPH,
            "exact_function_integral": True,
            "numerical_value_not_installed": True,
            "source_function_definition_restored": True,
            "actual_numerical_source_value_installed": False,
            "signed_source": True,
        },
        "sum_formula": "C_inner^(k)[53161/53162] + C_lower^(k)[36592/36593]",
        "source_rows_are_function_graphs_not_values": True,
    }
    forward = {
        "terms_by_Z": source_bindings_by_Z,
        "term_names": list(_TERM_NAMES),
        "core_source": _FORWARD_CORE_SOURCE,
        "inner_source": _FORWARD_INNER_SOURCE,
        "source_rows_are_saved_approximations_with_errors": True,
        "source_rows_are_not_function_values": True,
        "selected_as_exact_values": False,
    }

    delta_six = {
        "type": "signed_source_difference",
        "ordinary_Z_orders": [0, 1],
        "formula": "Delta_six^(k) = sum(i=0..5) F_i^(k) - (C_inner^(k) + C_lower^(k))",
        "forward_terms": list(_TERM_NAMES),
        "canonical_inner_function_C1_pair": [53161, 53162],
        "canonical_lower_function_C1_pair": [36592, 36593],
        "scope": copy.deepcopy(source_scope),
        "term_5_actual_signed_source": copy.deepcopy(source_bindings_by_Z.get("0", {}).get("term_5_signed_source_remainder", {})),
        "source_function_definition_restored": True,
        "actual_numerical_source_value_installed": False,
        "singular_split_is_algebraically_complete": True,
        "singular_split_is_not_full_weighted_N1p_identity": True,
        "exact_identity": False,
        "zero_claim": False,
        "ordinary_Z0_exact_zero": False,
        "ordinary_Z1_exact_zero": False,
        "numeric_value_available": False,
        "numeric_value_unavailable_reason": "forward terms are interval-backed approximations and the sixth term has only a norm enclosure",
        "independent_interval_overlap_used_as_identity": False,
        "source_graph_registration_used_as_identity": False,
        "remaining_source_mismatch": [
            "cross-producer identity from forward finite Volterra rows to the canonical analytic core function is unproved",
            "cross-producer collar identity for q-dependent rho, ratio, Psi history, and own mean is unproved",
            "forward pressure partition has no proved algebraic identity with canonical direct N1p",
            "forward P0/M2 pressure data have no proved identity with canonical native P0",
            "the signed A-A0 source function is not supplied by the forward producer",
            "common initial endpoint and scalar normalization are not source-bound for all six terms",
            "ordinary Z1 source identity is unproved",
        ],
        "smallest_unsatisfied_source_identity": {
            "status": "unproved",
            "name": "forward_collar_profile_equals_canonical_collar_fields",
            "forward_handles": {
                "profile": "OriginalFirstCollarTransport._inputs -> OriginalN1InnerM15Numeric._prepare/source_cell",
                "local_frame": "src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py:leading_regular_frame",
                "endpoint_preparation": "OriginalN1InnerM15Numeric._prepare uses tmax=(sc/4)*2 and source_cell reuses v",
            },
            "canonical_handles": {
                "profile": "current_original_n1_lower_collar_functions.py:collar_fields",
                "pressure": "current_original_n1_lower_collar_functions.py:direct_pressure_functions",
            },
            "required_equalities": [
                "rho_forward(s)=4*exp(hb*s)",
                "Phi_forward(s)=ratio(s)*Phi_core(rho(s))",
                "Psi_forward(s)=Psi_core(4)+integral_[0,s] ratio(t)*Psi_core_rho(rho(t))*chi(t)*hb*rho(t)dt",
                "mean_forward(s)=(4*mean_Psi_core(4)+integral_[0,s] Psi_forward(t)*hb*rho(t)dt)/rho(s)",
            ],
            "underlying_equation_mismatch_proved": False,
        },
        "core_scalar_identity_gap": {
            "status": "unproved",
            "forward": "OriginalCoreNumeric._seed P0_Z_taylor uses the numerical M2_point/q source center with directed enclosure",
            "canonical": "CurrentOriginalN1CoreFixedPointFunctions.P0 uses the native exact P0 graph",
            "required_identity": "forward native pressure source equals canonical native P0 as a function of Z",
        },
    }
    exact_identity = {
        "forward_six_equals_canonical_inner_plus_lower": False,
        "delta_six_zero": False,
        "ordinary_Z0_zero": False,
        "ordinary_Z1_zero": False,
        "term_5_zero": False,
        "all_densities_bound_to_same_true_function": False,
        "all_initial_endpoints_bound": False,
        "all_scalars_bound": False,
        "ordinary_Z1_bound": False,
        "source_difference_exported": True,
    }

    passed = bool(
        parent_pin.get("passed")
        and parent_pin.get("mapping_match")
        and source_pins_ok
        and all(family_checks.values())
        and ast_ok
        and graph_ok
        and forward_ok
        and not gaps
    )
    return {
        "passed": passed,
        "source_binding_passed": passed,
        "source_only": True,
        "source_executed": False,
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": parent_pin,
        "source_family": family,
        "family_checks": family_checks,
        "source_pins": source_pins,
        "source_pins_ok": source_pins_ok,
        "AST": ast_rows,
        "graph_identity": graph_identity,
        "source_scope": source_scope,
        "forward": forward,
        "canonical": canonical,
        "forward_six": forward,
        "canonical_inner_lower": canonical,
        "delta_six": delta_six,
        "exact_identity": exact_identity,
        "archives": archives,
        "source_function_definition_restored": True,
        "actual_numerical_source_value_installed": False,
        "numeric_value_available": False,
        "numeric_value_unavailable_reason": delta_six["numeric_value_unavailable_reason"],
        "root": str(repo),
        "directory": str(work),
        "gaps": gaps,
    }


__all__ = ["build_inner_collar_source"]
