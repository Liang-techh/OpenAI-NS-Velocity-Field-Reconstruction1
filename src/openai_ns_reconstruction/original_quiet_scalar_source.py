from __future__ import annotations

"""Source-only binding for the original quiet scalar increment.

This module reads pinned JSON or gzip receipts and parses pinned source ASTs.
It deliberately does not import or execute a producer, owner, constructor,
quadrature, recurrence, or numerical checker.  Saved decimal enclosures are
returned as provenance only; the caller must keep them as enclosures.
"""

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


_PARENT_RELATIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_joint_reserved_m15_numeric.json"
)
_PARENT_SHA256 = "a46413fd755aa96c385aa7b81e1919fe0d22838868b8b5a7d4348708dd16dad6"
_PREFIX = "lei_ren_part1_paper_compliant_current_original_"
_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
_REPORTS = {
    "axial_m15": _PREFIX + "axial_m15_numeric.json",
    "compact_quiet_m15": _PREFIX + "compact_quiet_m15_numeric.json",
    "preRp_graph": _PREFIX + "n1_preRp_moment_functions.json.gz",
    "native_frame": _PREFIX + "C3_Rp_native_frame.json.gz",
    "true_quiet": _PREFIX + "true_quiet_Rp_C5_pressure_C3_functions.json.gz",
    "quiet_source": _PREFIX + "n1_quiet_pressure_source.json.gz",
    # These saved bridge receipts carry the accepted physical-radius recipe.
    # They are read as provenance only; no checker or bridge producer is run.
    "preRp_bridge": _PREFIX + "n1_preRp_numeric_bridge.json.gz",
    "preRp_bridge_check": _PREFIX + "n1_preRp_numeric_bridge_check.json",
}
_SOURCES = {
    "quiet_functions": _PREFIX + "n1_preRp_moment_functions.py",
    "canonical_constructor": _PREFIX + "n1_quiet_pressure_source.py",
    "flat_sigma": "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py",
    "axial_pressure_binding": "src/openai_ns_reconstruction/original_axial_pressure_rows.py",
    "axial_memory_log_beta": "src/openai_ns_reconstruction/original_axial_memory.py",
    "axial_rate_evaluate": "src/openai_ns_reconstruction/original_axial_rate_memory.py",
    "compact_transport": "src/openai_ns_reconstruction/original_compact_pressure_rows.py",
    "forward_evaluate": "src/openai_ns_reconstruction/original_compact_quiet_m15.py",
    "mean_selection": "src/openai_ns_reconstruction/original_post_axial_pressure_rows.py",
    "post_axial_anchor": "src/openai_ns_reconstruction/original_post_axial_m15.py",
    "axial_rate_memory": "src/openai_ns_reconstruction/original_axial_rate_memory.py",
    "axial_geometry": "src/openai_ns_reconstruction/original_axial_pressure_rows.py",
    "core_log_amplitude": "src/openai_ns_reconstruction/original_core_log_amplitude.py",
    "reference_numeric": "src/openai_ns_reconstruction/original_reference_numeric.py",
    "preRp_bridge": _PREFIX + "n1_preRp_numeric_bridge.py",
    "preRp_bridge_check": _PREFIX + "n1_preRp_numeric_bridge_check.py",
    "preRc_endpoint": _PREFIX + "n1_preRc_pressure_endpoint.py",
}
_IMPLEMENTATIONS = tuple(_SOURCES.values())


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


def _find(repo: Path, work: Path, key: str) -> Path:
    raw = Path(key)
    candidates = [raw] if raw.is_absolute() else [work / raw, repo / raw, repo / "experiments" / "root_st073" / raw]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _load(path: Path) -> Mapping[str, Any]:
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    value = json.loads(raw)
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _function(tree: ast.AST, name: str, cls: str | None = None) -> ast.AST:
    scope: ast.AST = tree
    if cls is not None:
        classes = [n for n in ast.walk(tree) if isinstance(n, ast.ClassDef) and n.name == cls]
        if len(classes) != 1:
            raise ValueError(f"unique class required: {cls}")
        scope = classes[0]
    functions = [n for n in ast.walk(scope) if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == name]
    if len(functions) != 1:
        raise ValueError(f"unique function required: {name}")
    return functions[0]


def _statement(text: str) -> ast.AST:
    tree = ast.parse(text)
    if len(tree.body) != 1:
        raise ValueError(text)
    return tree.body[0]


def _contains(node: ast.AST, wanted: ast.AST) -> bool:
    target = ast.dump(wanted, include_attributes=False)
    return any(ast.dump(part, include_attributes=False) == target for part in ast.walk(node))


def _ast_guard(
    path: Path,
    repo: Path,
    name: str,
    *,
    cls: str | None = None,
    statements: tuple[str, ...] = (),
    text: tuple[str, ...] = (),
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_path": _relative(path, repo),
        "function": name,
        "class": cls,
        "required_AST": {value: False for value in statements},
        "required_text": {value: False for value in text},
        "function_AST_sha256": None,
        "source_executed": False,
        "passed": False,
    }
    try:
        source = path.read_text(encoding="utf-8")
        function = _function(ast.parse(source, filename=str(path)), name, cls)
        row["required_AST"] = {value: _contains(function, _statement(value)) for value in statements}
        segment = ast.get_source_segment(source, function) or ""
        compact = "".join(segment.split())
        row["required_text"] = {value: "".join(value.split()) in compact for value in text}
        row["function_AST_sha256"] = hashlib.sha256(ast.dump(function, include_attributes=False).encode()).hexdigest()
        row["passed"] = all(row["required_AST"].values()) and all(row["required_text"].values())
    except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
        row["error"] = str(exc)
    return row


def _source_module_ast_sha(path: Path, name: str, cls: str | None = None) -> str | None:
    """Hash a parsed source segment exactly as the saved graph fingerprint does."""
    try:
        source = path.read_text(encoding="utf-8")
        function = _function(ast.parse(source, filename=str(path)), name, cls)
        segment = ast.get_source_segment(source, function) or ""
        return hashlib.sha256(ast.dump(ast.parse(segment), include_attributes=False).encode()).hexdigest()
    except (OSError, SyntaxError, UnicodeError, ValueError):
        return None


def _node(nodes: Any, index: int) -> Mapping[str, Any]:
    if isinstance(nodes, list) and 0 <= index < len(nodes) and isinstance(nodes[index], Mapping):
        return nodes[index]
    return {}


def _node_proofs(nodes: Any, expected: Mapping[int, Mapping[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for index, wanted in expected.items():
        actual = dict(_node(nodes, index))
        result[str(index)] = {"expected": dict(wanted), "actual": actual, "passed": actual == dict(wanted)}
    return result


def _node_field_proofs(nodes: Any, expected: Mapping[int, Mapping[str, Any]]) -> dict[str, Any]:
    """Check selected saved node fields without replaying the graph."""
    result: dict[str, Any] = {}
    for index, wanted in expected.items():
        actual = dict(_node(nodes, index))
        passed = all(actual.get(key) == value for key, value in wanted.items())
        result[str(index)] = {"expected_fields": dict(wanted),
                              "actual_fields": {key: actual.get(key) for key in wanted},
                              "passed": passed, "source_executed": False}
    return result


def _pin(path: Path, repo: Path, expected: Any, name: str, kind: str) -> dict[str, Any]:
    actual = _sha(path)
    return {
        "name": name,
        "kind": kind,
        "key": path.name if kind == "input" else _relative(path, repo),
        "path": _relative(path, repo),
        "sha256": actual,
        "expected_sha256": expected,
        "parent_pin_present": isinstance(expected, str),
        "passed": isinstance(expected, str) and actual == expected,
        "source_executed": False,
    }


def _verify_map(
    hashes: Mapping[str, Any],
    repo: Path,
    work: Path,
    *,
    label: str,
    source_pins: dict[str, Any],
    selected: dict[str, str],
    gaps: list[str],
) -> tuple[int, int]:
    passed = total = 0
    for key, expected in hashes.items():
        if not isinstance(key, str):
            continue
        total += 1
        path = _find(repo, work, key)
        actual = _sha(path)
        ok = isinstance(expected, str) and actual == expected
        source_pins[f"transitive:{label}:{key}"] = {
            "name": f"{label}:{key}",
            "kind": f"transitive_{label}",
            "key": key,
            "path": _relative(path, repo),
            "sha256": actual,
            "expected_sha256": expected,
            "parent_expected_sha256": expected,
            "passed": ok,
            "source_executed": False,
        }
        if ok:
            passed += 1
            selected[key] = actual
        else:
            gaps.append(f"{label} hash failed: {key}")
    return passed, total


def _parent_path(repo: Path, work: Path, value: Any) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    candidate = repo / path
    return candidate if candidate.exists() else work / path


def _load_parent(value: Any, repo: Path, work: Path) -> tuple[Mapping[str, Any], dict[str, Any], list[str]]:
    latest_path = repo / _PARENT_RELATIVE
    latest_pin = {
        "path": _relative(latest_path, repo),
        "sha256": _sha(latest_path),
        "expected_sha256": _PARENT_SHA256,
        "passed": _sha(latest_path) == _PARENT_SHA256,
        "source_executed": False,
    }
    gaps: list[str] = []
    if not latest_pin["passed"]:
        gaps.append("latest parent report pin mismatch")
    try:
        latest = _load(latest_path)
    except (OSError, ValueError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
        return {}, latest_pin, gaps + [f"cannot load latest parent: {exc}"]
    try:
        supplied = value if isinstance(value, Mapping) else _load(_parent_path(repo, work, value))
    except (OSError, ValueError, TypeError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
        return {}, latest_pin, gaps + [f"cannot load supplied parent: {exc}"]
    latest_digest, supplied_digest = _digest(latest), _digest(supplied)
    latest_pin.update({"latest_mapping_sha256": latest_digest, "supplied_mapping_sha256": supplied_digest,
                       "mapping_match": latest_digest == supplied_digest})
    if latest_digest != supplied_digest:
        gaps.append("supplied parent mapping does not match pinned latest report")
    return supplied, latest_pin, gaps


def _empty(gaps: list[str], repo: Path, work: Path, pin: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "passed": False,
        "source_only": True,
        "source_executed": False,
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": dict(pin),
        "source_family": dict(_FAMILY),
        "input_hashes": {},
        "implementation_hashes": {},
        "source_pins": {},
        "source_pins_ok": False,
        "added_pins": {},
        "transitive_hash_validation": {"passed": False, "input_passed": 0, "input_total": 0,
                                        "implementation_passed": 0, "implementation_total": 0},
        "AST": {},
        "graph_identity": {},
        "forward_quiet": {},
        "forward": {},
        "canonical_quiet": {},
        "canonical": {},
        "mean_source": {},
        "mean_source_summary": {},
        "quiet_identity": {},
        "quiet_scalar_identity": {},
        "physical_scale_guards": {},
        "numeric_value_available": False,
        "numeric_value_unavailable_reason": "source binding only; saved enclosures are not function values",
        "root": str(repo),
        "directory": str(work),
        "gaps": list(gaps),
    }


def _expected(mapping: Mapping[str, Any], key: str, path: Path) -> Any:
    return mapping.get(key, mapping.get(path.name))


def _rows_for_coordinate(report: Mapping[str, Any], chart: str | None = None) -> list[tuple[int, Mapping[str, Any]]]:
    result = []
    for index, row in enumerate(report.get("samples", ())):
        if not isinstance(row, Mapping) or row.get("native_coordinate") != "1":
            continue
        if chart is not None and row.get("chart") != chart:
            continue
        result.append((index, row))
    return result


def build_quiet_scalar_source(parent_report: Any, *, root: Any, directory: Any) -> dict[str, Any]:
    """Return source-bound forward and canonical quiet scalar handles."""
    repo = _repo(root, directory)
    work = _work(repo, directory)
    parent, parent_pin, gaps = _load_parent(parent_report, repo, work)
    if not parent:
        return _empty(gaps, repo, work, parent_pin)
    if _mapping(parent.get("source_family")) != _FAMILY:
        gaps.append("parent source_family mismatch")
    if _mapping(parent.get("checks")).get("passed") is not True:
        gaps.append("parent checks.passed is not true")

    parent_inputs = _mapping(parent.get("input_hashes"))
    parent_impl = _mapping(parent.get("implementation_hashes"))
    source_pins: dict[str, Any] = {}
    selected_inputs: dict[str, str] = {}
    selected_impl: dict[str, str] = {}
    input_passed, input_total = _verify_map(parent_inputs, repo, work, label="input",
                                            source_pins=source_pins, selected=selected_inputs, gaps=gaps)
    impl_passed, impl_total = _verify_map(parent_impl, repo, work, label="implementation",
                                          source_pins=source_pins, selected=selected_impl, gaps=gaps)

    paths: dict[str, Path] = {}
    reports: dict[str, Mapping[str, Any]] = {}
    for key, basename in _REPORTS.items():
        path = work / basename
        paths[key] = path
        pin = _pin(path, repo, _expected(parent_inputs, basename, path), key + "_report", "input")
        source_pins[key + "_report"] = pin
        if not pin["passed"]:
            gaps.append("report pin failed: " + basename)
        try:
            reports[key] = _load(path)
        except (OSError, ValueError, json.JSONDecodeError, gzip.BadGzipFile) as exc:
            reports[key] = {}
            gaps.append(f"cannot load {key}: {exc}")

    source_paths: dict[str, Path] = {}
    for key, name in _SOURCES.items():
        path = repo / name if name.startswith("src/") else work / name
        source_paths[key] = path
        if name.startswith("src/"):
            expected, kind = _expected(parent_impl, name, path), "implementation"
        else:
            expected, kind = _expected(parent_inputs, name, path), "input"
        pin = _pin(path, repo, expected, key + "_source", kind)
        source_pins[key + "_source"] = pin
        if not pin["passed"]:
            gaps.append("source pin failed: " + name)
    source_pins_ok = all(pin.get("passed") is True for pin in source_pins.values())

    family_checks = {key: _mapping(report.get("source_family")) == _FAMILY for key, report in reports.items()}
    if not all(family_checks.values()):
        gaps.append("quiet source report family mismatch")

    ast_rows = {
        "quiet_functions": _ast_guard(
            source_paths["quiet_functions"], repo, "quiet_functions",
            statements=(
                "original=linear.realization.restore_graph(raw['exact_graph_nodes'])",
                "kappa=ref(f['actual_identified_Rp_histories_C3']['m'][1])",
                "K=g.mul(Rp,scale,kappa)",
                "m15=g.mul(g.constant('-1/4'),K,K,logratio)",
            ),
            text=("actual_identified_Rp_histories_C3", "constant_Z_identity=quiet_owner.actual_source_identity()",
                  "source.C1Function(m15,g.zero)"),
        ),
        "canonical_constructor": _ast_guard(
            source_paths["canonical_constructor"], repo, "__init__",
            cls="CurrentOriginalN1QuietPressureSource",
            statements=(
                "Y=s.exp(40)",
                "kernel=s.cancel(kappa*s.exp(Tw+logP+12)/4-s.exp(1-Y))",
                "self.kernel_node=int(str(kernel.func).rsplit('_',1)[1])",
            ),
            text=(
                "expected=s.exp(y-u)*s.Function('original_sigma')(1-s.log(y)/40)",
                "node['operation']!='definite_integral'",
                "if s.simplify(actual-expected)!=0 or _graph_expression(self.nodes,node['lower'],{})!=1",
                "if s.simplify(M-s.exp(logP)*kappa)!=0",
            ),
        ),
        "flat_sigma": _ast_guard(
            source_paths["flat_sigma"], repo, "sigma_jets",
            statements=("x=c.mpf(x)", "rows=[]", "return hull_jets(c,rows)"),
            text=("Ordinary Taylor coefficients at every point of a real interval",
                  "x=c.mpf(x);lo,hi=endpoints(x);rows=[]"),
        ),
        "axial_rate_evaluate": _ast_guard(
            source_paths["axial_rate_evaluate"], repo, "evaluate", cls="OriginalAxialRateMemory",
            statements=(
                "dp=c.expm1(40*sp)",
                "U=c.exp(40*sp)",
                "r0=1-sp",
                "if power not in (1,2):raise ValueError('Original beta power1/2 required')",
                "rvalue=initial-ctx.ln(1-(qorigin+scale*z)/(ctx.mpf(U) if ctx is c else Ui))/40",
                "u=(ctx.mpf(U) if ctx is c else Ui)-qorigin-scale*z",
                "result=iv.exp(-ratei*dt)+iv.exp(iv.mpf(anchor))*iv.mpf(scale)*integral",
            ),
            text=(
                "physical_measure='dq=d(d-v); original exponential Jacobian already absorbed'",
                "local_axial_zero_does_not_zero_memory=True",
            ),
        ),
        "axial_memory_log_beta": _ast_guard(
            source_paths["axial_memory_log_beta"], repo, "_log_beta", cls="OriginalAxialMemory",
            statements=("sj=self.sigma_jets(ctx,r[0])", "odds=(1-r)**-2-r**-2"),
            text=("return odds-self._ln(1+odds.exp())", "if ctx is self.iv and lo>=self.ctx.mpf('.25')"),
        ),
        "axial_pressure_binding": _ast_guard(
            source_paths["axial_pressure_binding"], repo, "__init__", cls="OriginalAxialPressureRows",
            statements=("self.memory=OriginalAxialRateMemory(self.ctx,self.iv,self.J)",),
        ),
        "axial_rate_default": _ast_guard(
            source_paths["axial_rate_memory"], repo, "__init__", cls="OriginalAxialRateMemory",
            statements=(
                "if kernel_rate not in (1,2):raise ValueError('Original pressure kernel rate one or two required')",
                "self.rate,self.rate_iv=ctx.mpf(kernel_rate),iv.mpf(kernel_rate)",
            ),
            text=("kernel_rate=1",),
        ),
        "axial_pressure_memory": _ast_guard(
            source_paths["axial_pressure_binding"], repo, "scalars", cls="OriginalAxialPressureRows",
            statements=("memory=self.memory.evaluate(value,1)",),
        ),
        "axial_pressure_fields": _ast_guard(
            source_paths["axial_pressure_binding"], repo, "fields", cls="OriginalAxialPressureRows",
            statements=("M=4*z*mean",),
            text=("Exact analytic cancellation in the original mean recovery",),
        ),
        "compact_transport": _ast_guard(
            source_paths["compact_transport"], repo, "scalars", cls="OriginalCompactPressureRows",
            statements=("_,p,t=self.coordinate(value)",),
            text=("mean=self.mean_point*self.ctx.exp(-14)/p",),
        ),
        "quiet_transport": _ast_guard(
            source_paths["compact_transport"], repo, "scalars", cls="OriginalQuietPressureRows",
            statements=("g=self.geometry(value)",),
            text=("mean=self.mean_point*(self.ctx.exp(-14)/2)*self.ctx.exp(-g['local_point'])",),
        ),
        "compact_geometry": _ast_guard(
            source_paths["compact_transport"], repo, "geometry", cls="OriginalCompactPressureRows",
            text=("terms=list(self.axial.geometry('1')['atomic_log_radius_terms'])+[(c.mpf(14),iv.mpf(14)),(sp,st)]",
                  "physical_jacobian_point=1/p"),
        ),
        "quiet_geometry": _ast_guard(
            source_paths["compact_transport"], repo, "geometry", cls="OriginalQuietPressureRows",
            text=("terms=list(self.axial.geometry('1')['atomic_log_radius_terms'])+[(c.mpf(14),iv.mpf(14)),(c.log(2),iv.log(2)),(sp,st)]",),
        ),
        "axial_geometry": _ast_guard(
            source_paths["axial_geometry"], repo, "geometry", cls="OriginalAxialPressureRows",
            statements=("dp,dt=c.expm1(40*sp),iv.expm1(40*st)",),
            text=("terms=[(c.log(110),iv.log(110)),(10*ref.logCp,10*ref.logC),(10*ref.logPp,10*ref.logP),(c.mpf(1),iv.mpf(1)),(dp,dt)]",),
        ),
        "quiet_length_positive": _ast_guard(
            source_paths["compact_transport"], repo, "__init__", cls="OriginalQuietPressureRows",
            statements=(
                "self.length=-60*self.logmu-2-self.ctx.log(2)",
                "self.lengthi=-60*self.logmui-2-self.iv.log(2)",
            ),
            text=("Original positive quiet length required",),
        ),
        "core_log_amplitude": _ast_guard(
            source_paths["core_log_amplitude"], repo, "__init__", cls="OriginalCoreLogAmplitude",
            statements=("self.lp,self.lp_iv=c.exp(40)+11,iv.exp(40)+11",),
            text=("selected_logCstar", "must be an exact singleton"),
        ),
        "reference_numeric": _ast_guard(
            source_paths["reference_numeric"], repo, "__init__", cls="OriginalReferenceNumeric",
            statements=(
                "a=self.owner.amplitude",
                "self.logC,self.logCp=iv.mpf(a.logC_iv),c.mpf(a.logC)",
                "self.logP,self.logPp=iv.mpf(a.lp_iv),c.mpf(a.lp)",
            ),
        ),
        "preRp_bridge_identity": _ast_guard(
            source_paths["preRp_bridge"], repo, "assert_source_identity",
            text=(
                "same_original_radius_units_shared_Rc_cancelled_before_boxing=True",
                "checked_C1_input_and_endpoint_receipt_chains_are_authoritative=True",
                "original_physical_units=copy.deepcopy(self.numeric.parents['old']['original_physical_units'])",
            ),
        ),
        "preRp_bridge_check": _ast_guard(
            source_paths["preRp_bridge_check"], repo, "run",
            text=("owner.source_identity['passed'] and raw['source_identity']==owner.source_identity",),
        ),
        "preRc_geometry_recipe": _ast_guard(
            source_paths["preRc_endpoint"], repo, "_symbolic_geometry",
            statements=(
                "expected_log_rc=source_c.ln(110)+10*log_c+11*source_c.exp(40)+124",
                "expected_log_rp=expected_log_rc-2-60*source_c.ln(mu_box)",
            ),
            text=(
                "source_logRc_recipe",
                "log(110)+10*selected_logCstar+11*exp(40)+124",
                "source_logRp_recipe",
                "logRc-2-60*log(mu)",
            ),
        ),
        "forward_evaluate": _ast_guard(
            source_paths["forward_evaluate"], repo, "evaluate", cls="OriginalCompactQuietM15",
            statements=(
                "anchor=self.rows.axial.geometry('1')",
                "coefficientp=J.constant(c,-4*self.rows.mean_point**2*geo['local_point'],1)",
                "coefficientt=J.constant(iv,-4*self.rows.mean_truth**2*geo['local'],1)",
            ),
            text=("exact_original_mean_square_integral", "original_positive_axial_terminal_mean",
                  "global_radius_and_global_mean_paired=True", "half_S1_already_included=True"),
        ),
        "mean_selection": _ast_guard(
            source_paths["mean_selection"], repo, "__init__", cls="OriginalPostAxialPressureRows",
            statements=(
                "endpoints=[p for p in data['pressure_row_samples'] if p['native_coordinate']=='1']",
                "self.mean_source_summary=copy.deepcopy(endpoints[0]['leading_scalars']['mean'])",
                "self.mean_point=self.ctx.mpf(self.mean_source_summary['approximation'])",
                "self.mean_truth=self.iv.mpf([self.ctx.make_mpf(tuple(v)) for v in self.mean_source_summary['source_enclosure_exact_mpf_tuples']])",
            ),
            text=("Original scalar mean must be independent of Z", "Original nonzero terminal memory must be preserved"),
        ),
        "post_axial_anchor": _ast_guard(
            source_paths["post_axial_anchor"], repo, "evaluate", cls="OriginalPostAxialM15",
            statements=("anchor=self.rows.axial.geometry('1')", "mp,mt=self.rows.mean_point,self.rows.mean_truth"),
            text=("exact_mean_square_integral", "global_mean_and_radius_anchor=True"),
        ),
    }
    ast_passed = all(row.get("passed") is True for row in ast_rows.values())
    if not ast_passed:
        gaps.append("quiet source AST guard failed")

    # The post-axial source selects one saved positive mean summary from the
    # axial coordinate 1 rows and checks equality across Z.  Keep all of the
    # summary, including its exact enclosure tuples, but never call it a value.
    axial_endpoints = [
        (index, row) for index, row in enumerate(reports["axial_m15"].get("pressure_row_samples", ()))
        if isinstance(row, Mapping) and row.get("native_coordinate") == "1"
    ]
    mean_summary: Mapping[str, Any] = {}
    mean_proof: Mapping[str, Any] = {}
    mean_equal = bool(axial_endpoints)
    if axial_endpoints:
        first = axial_endpoints[0][1]
        mean_summary = _mapping(_mapping(first).get("leading_scalars")).get("mean", {})
        mean_proof = _mapping(first.get("positive_mean_memory_receipt"))
        mean_equal = all(_mapping(_mapping(row).get("leading_scalars")).get("mean") == mean_summary
                         for _, row in axial_endpoints)
    if not axial_endpoints:
        gaps.append("no saved axial pressure row at native coordinate 1")
    elif not mean_equal:
        gaps.append("saved axial mean summaries differ across coordinate 1 rows")
    if not isinstance(mean_summary, Mapping) or not mean_summary.get("source_enclosure_exact_mpf_tuples"):
        gaps.append("saved axial mean enclosure is missing")

    mean_source = {
        "summary": dict(mean_summary),
        "positive_memory_receipt": dict(mean_proof),
        "source_report": _relative(paths["axial_m15"], repo),
        "source_path": "pressure_row_samples[native_coordinate=1].leading_scalars.mean",
        "endpoint_rows": [index for index, _ in axial_endpoints],
        "same_summary_for_all_Z": mean_equal,
        "selected_as_exact_function_value": False,
        "decimal_approximation_is_exact": False,
        "source_enclosure_is_function_value": False,
    }

    quiet_rows = _rows_for_coordinate(reports["compact_quiet_m15"], "quiet")
    forward_receipts: dict[str, Any] = {}
    forward_ok = bool(quiet_rows)
    for index, row in quiet_rows:
        terms = row.get("m15_C0_C1_terms")
        incoming = row.get("incoming_m15_terms")
        exact = terms[-1] if isinstance(terms, list) and terms else {}
        binding = _mapping(_mapping(exact).get("actual_function_binding"))
        binding_ok = (
            isinstance(terms, list) and len(terms) == 212
            and isinstance(incoming, list) and len(incoming) == 211
            and str(binding.get("exact_primitive", "")).startswith("-4*R_axial_end")
            and binding.get("common_source_function") == "original_positive_axial_terminal_mean"
            and binding.get("half_S1_already_included") is True
            and binding.get("global_radius_and_global_mean_paired") is True
        )
        if not binding_ok:
            forward_ok = False
        forward_receipts[str(row.get("Z"))] = {
            "sample_index": index,
            "chart": row.get("chart"),
            "native_coordinate": row.get("native_coordinate"),
            "incoming_term_count": len(incoming) if isinstance(incoming, list) else None,
            "saved_output_term_count": len(terms) if isinstance(terms, list) else None,
            "exact_increment_binding": dict(binding),
            "log_scale_terms": list(exact.get("log_scale_terms", ())) if isinstance(exact, Mapping) else [],
            "atomic_log_radius_terms": list(_mapping(row.get("geometry")).get("atomic_log_radius_terms", ())),
            "physical_local_log_distance": _mapping(row.get("geometry")).get("physical_local_log_distance"),
            "geometry_path": "samples[chart=quiet,native_coordinate=1].geometry.physical_local_log_distance",
            "scale_rows_path": "samples[chart=quiet,native_coordinate=1].m15_C0_C1_terms[-1].log_scale_terms",
            "selected_as_exact_value": False,
            "enclosure_only": True,
        }
    if not forward_ok:
        gaps.append("saved forward quiet increment receipt failed")

    quiet_source = _mapping(reports["quiet_source"].get("actual_quiet_source"))
    quiet_source_identity = _mapping(reports["quiet_source"].get("source_identity"))
    quiet_source_guards = {
        "source_scaled_Z1_exact_zero": quiet_source.get("source_scaled_Z1_exact_zero") is True,
        "physical_Rc_is_absolute": quiet_source_identity.get("physical_Rc_is_absolute_not_original_phase_offset") is True,
        "signed_coefficients": quiet_source.get("signed_coefficients_not_magnitude_caps") is True,
        "quiet_length_node_present": isinstance(quiet_source_identity.get("exact_original_quiet_length_node"), int),
        "source_identity_available": bool(quiet_source_identity),
    }
    if not all(quiet_source_guards.values()):
        gaps.append("saved quiet source identity guard failed")

    pre = reports["preRp_graph"]
    defining = _mapping(pre.get("defining_bindings"))
    quiet_binding = _mapping(defining.get("quiet"))
    signed_pieces = _mapping(pre.get("signed_pieces"))
    signed_functions = _mapping(pre.get("signed_pressure_C1_functions"))
    graph_expected = {
        31204: {"operation": "product", "arguments": [31199, 31203]},
        31205: {"operation": "product", "arguments": [31162, 25203, 31204]},
        31206: {"operation": "product", "arguments": [48, 30828]},
        31213: {"operation": "sum", "arguments": [31161, 31212]},
        31215: {"operation": "product", "arguments": [31214, 31205, 31205, 31213]},
        31219: {"operation": "sum", "arguments": [30707, 31155]},
        31220: {"operation": "sum", "arguments": [30708, 31156]},
        31221: {"operation": "sum", "arguments": [31219, 31215]},
    }
    graph_proofs = _node_proofs(pre.get("exact_graph_nodes"), graph_expected)
    graph_ok = all(row.get("passed") is True for row in graph_proofs.values())
    graph_checks = {
        "quiet_m15_pair": signed_pieces.get("quiet", {}).get("m15") == [31215, 0],
        "preRp_m15_pair": signed_functions.get("m15") == [31221, 31220],
        "quiet_K_binding": quiet_binding.get("K") == 31205,
        "quiet_kappa_binding": quiet_binding.get("kappa") == 31204,
        "quiet_lower_binding": quiet_binding.get("lower") == 31206,
        "quiet_upper_binding": quiet_binding.get("upper") == 31162,
        "ordinary_quiet_C1_zero": signed_pieces.get("quiet", {}).get("m15", [None, None])[1] == 0,
    }
    graph_ok = graph_ok and all(graph_checks.values())
    if not graph_ok:
        gaps.append("canonical quiet graph binding failed")

    native = _mapping(reports["native_frame"].get("actual_Rp_native_frame_functions"))
    native_histories = _mapping(native.get("actual_identified_Rp_histories_C3"))
    native_m = native_histories.get("m")
    kappa_frame_node = native_m[1] if isinstance(native_m, list) and len(native_m) > 1 else None
    if kappa_frame_node is None:
        gaps.append("actual native frame kappa handle missing")

    # The native frame stores the exact turnoff integral used by the
    # canonical constructor.  Check its endpoint, Jacobian-bearing
    # integrand, and Z independence by fields only; never rebuild it.
    native_kernel_expected = {
        1: {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        2: {"operation": "exact_rational", "numerator": 40, "denominator": 1},
        25064: {"operation": "bound_variable", "name": "original_O2_turnoff_integrand_y"},
        25065: {"operation": "bound_variable", "name": "original_O2_turnoff_endpoint_y"},
        25066: {"operation": "analytic_unary", "name": "log", "argument": 25064},
        25067: {"operation": "positive_quotient", "numerator": 25066, "denominator": 2},
        25068: {"operation": "negative", "argument": 25067},
        25069: {"operation": "sum", "arguments": [1, 25068]},
        25070: {"operation": "exact_original_source_flat_sigma", "argument": 25069,
                "smooth_flat_endpoints": True, "Z_independent_argument": True},
        25071: {"operation": "negative", "argument": 25065},
        25072: {"operation": "sum", "arguments": [25064, 25071]},
        25073: {"operation": "analytic_unary", "name": "exp", "argument": 25072},
        25074: {"operation": "product", "arguments": [25073, 25070]},
        25075: {"operation": "definite_integral", "integrand": 25074,
                "variable": "original_O2_turnoff_integrand_y", "lower": 1, "upper": 25065,
                "exact_function_integral": True, "numerical_value_not_installed": True,
                "original_kernel_function": True, "coordinate_and_kernel_Z_independent": True},
    }
    native_kernel_proofs = _node_field_proofs(reports["native_frame"].get("exact_graph_nodes"), native_kernel_expected)
    native_kernel_ok = all(row.get("passed") is True for row in native_kernel_proofs.values())
    native_constant_one = native_kernel_proofs.get("1", {}).get("passed") is True
    native_constant_forty = native_kernel_proofs.get("2", {}).get("passed") is True
    native_sigma_binding = _mapping(_node(reports["native_frame"].get("exact_graph_nodes"), 25070).get("source_binding"))
    native_sigma_source_hash = _sha(source_paths["flat_sigma"])
    # The producer recorded the AST of a one-function parsed module (rather
    # than the bare FunctionDef node).  Keep both fingerprints: the local
    # guard proves the selector body and this module fingerprint proves the
    # exact raw25070 binding.
    native_sigma_module_ast_hash = _source_module_ast_sha(source_paths["flat_sigma"], "sigma_jets")
    native_sigma_ast_match = native_sigma_module_ast_hash == native_sigma_binding.get("AST_sha256")
    native_sigma_source_ok = (
        native_sigma_source_hash == native_sigma_binding.get("module_sha256")
        and native_sigma_ast_match
        and native_sigma_binding.get("function") == "sigma_jets"
    )

    # The bridge receipt is the accepted source identity for the physical
    # radius.  The pre-Rp archive is not executed here: only its saved
    # identity and its checked receipt are read.  This prevents a family
    # label, by itself, from proving the Rc pairing.
    preRp_bridge = reports["preRp_bridge"]
    preRp_bridge_check = reports["preRp_bridge_check"]
    preRp_identity = _mapping(preRp_bridge.get("source_identity"))
    preRp_check_identity = _mapping(preRp_bridge_check.get("source_identity"))
    preRp_units = _mapping(preRp_check_identity.get("original_physical_units"))
    preRp_radius_checks = {
        "bridge_receipt_passed": preRp_bridge_check.get("all_passed") is True,
        "bridge_identity_passed": preRp_identity.get("passed") is True,
        "checked_identity_passed": preRp_check_identity.get("passed") is True,
        "shared_Rc_units": preRp_identity.get("same_original_radius_units_shared_Rc_cancelled_before_boxing") is True
        and preRp_check_identity.get("same_original_radius_units_shared_Rc_cancelled_before_boxing") is True,
        "authoritative_endpoint_chain": preRp_identity.get("checked_C1_input_and_endpoint_receipt_chains_are_authoritative") is True
        and preRp_check_identity.get("checked_C1_input_and_endpoint_receipt_chains_are_authoritative") is True,
        "source_logRc_formula": preRp_units.get("source_logRc_formula") ==
        "log(110)+10*logCstar+11*exp(40)+124",
        "source_logRp_formula": preRp_units.get("source_logRp_formula") ==
        "logRc-2-60*log(mu)",
        "source_logRp_node": preRp_units.get("source_Rp_node") == native.get("original_logRp") == 26687,
        "source_native_frame": preRp_units.get("source_native_frame") == paths["native_frame"].name,
        "source_ast_bridge": ast_rows["preRp_bridge_identity"].get("passed") is True,
        "source_ast_check": ast_rows["preRp_bridge_check"].get("passed") is True,
        "source_ast_endpoint": ast_rows["preRc_geometry_recipe"].get("passed") is True,
    }
    preRp_radius_ok = all(preRp_radius_checks.values())
    axial_radius_recipe = {
        "axial_geometry_at_one":
        "log(110)+10*logCstar+10*logP+1+expm1(40)",
        "reference_logCstar": "OriginalReferenceNumeric.logC <- OriginalCoreLogAmplitude.logC",
        "reference_logP": "OriginalReferenceNumeric.logP <- OriginalCoreLogAmplitude.lp",
        "source_logP_recipe": "logP=exp(40)+11",
        "canonical_logRc_recipe": "log(110)+10*logCstar+11*exp(40)+124",
        "canonical_minus_axial_log_offset": "14",
        "canonical_Rc_equals_axial_end_exp14": True,
    }
    axial_radius_checks = {
        "axial_geometry_ast": ast_rows["axial_geometry"].get("passed") is True,
        "amplitude_logP_ast": ast_rows["core_log_amplitude"].get("passed") is True,
        "reference_logC_logP_ast": ast_rows["reference_numeric"].get("passed") is True,
        "preRp_acceptance": preRp_radius_ok,
        "canonical_Rc_offset": preRp_radius_ok and axial_radius_recipe["canonical_Rc_equals_axial_end_exp14"],
    }
    axial_radius_ok = all(axial_radius_checks.values())
    source_identity_checks = {
        "actual_kernel_node": quiet_source_identity.get("actual_kernel_node") == 25075,
        "integrand": quiet_source_identity.get("original_defining_integrand") ==
        "original_sigma(1 - log(y)/40)*exp(-u + y)",
        "cutoff": quiet_source_identity.get("exact_cutoff") == "exp(40)",
        "physical_Rc_node": quiet_source_identity.get("physical_Rc_node") == native.get("original_logRc"),
        "physical_K_source": quiet_source_identity.get("physical_K") ==
        "4*Rc*exp(-14)*(kernel_at_exp40+exp(1-exp40))",
        "Pstar_cancelled": quiet_source_identity.get("Pstar_equals_exp_logP_cancelled_before_evaluation") is True,
        "quiet_length_node": quiet_source_identity.get("exact_original_quiet_length_node") == native.get("quiet_length"),
        "flat_sigma_source": native_sigma_source_ok,
        "constant_one": native_constant_one,
        "constant_forty": native_constant_forty,
        "canonical_radius_receipt": preRp_radius_ok,
        "canonical_radius_recipe": axial_radius_ok,
    }
    source_transport_checks = {
        "forward_mean_is_transportable": ast_rows["axial_pressure_memory"].get("passed") is True
        and ast_rows["axial_pressure_fields"].get("passed") is True,
        "compact_mean_scale": ast_rows["compact_transport"].get("passed") is True,
        "quiet_mean_scale": ast_rows["quiet_transport"].get("passed") is True,
        "compact_radius_scale": ast_rows["compact_geometry"].get("passed") is True,
        "quiet_radius_scale": ast_rows["quiet_geometry"].get("passed") is True,
        "sigma_log_source": ast_rows["axial_memory_log_beta"].get("passed") is True,
        "rate_one_constructor": ast_rows["axial_rate_default"].get("passed") is True,
        "rate_one_call": ast_rows["axial_pressure_memory"].get("passed") is True,
        "axial_radius_recipe": axial_radius_ok,
        "positive_quiet_length_source": ast_rows["quiet_length_positive"].get("passed") is True,
        "physical_Rc_anchor": source_identity_checks["physical_Rc_node"]
        and ast_rows["compact_geometry"].get("passed") is True
        and ast_rows["quiet_geometry"].get("passed") is True,
    }
    source_change_of_variables = {
        "Y": "exp(40)",
        "substitution": "y=Y-q",
        "forward_q_interval": "q in [0,Y-1]",
        "canonical_y_interval": "y in [1,Y]",
        "jacobian": "dq=-dy; reversed bounds remove the minus sign",
        "forward_integrand": "exp(-q)*sigma(1-log(Y-q)/40)",
        "canonical_integrand": "exp(-Y+y)*sigma(1-log(y)/40)",
        "canonical_integrand_equals_forward": "exp(-Y+y)=exp(-q)",
        "homogeneous_endpoint": "exp(1-Y)",
        "source_graph_node": 25075,
        "native_kernel_graph": native_kernel_proofs,
        "native_sigma_binding": {
            "module": native_sigma_binding.get("module"),
            "module_sha256": native_sigma_binding.get("module_sha256"),
            "function": native_sigma_binding.get("function"),
            "AST_sha256": native_sigma_binding.get("AST_sha256"),
            "module_hash_actual": native_sigma_source_hash,
            "module_AST_sha256_actual": native_sigma_module_ast_hash,
            "function_AST_sha256_actual": ast_rows["flat_sigma"].get("function_AST_sha256"),
            "AST_hash_match": native_sigma_ast_match,
            "passed": native_sigma_source_ok,
        },
        "source_identity_checks": source_identity_checks,
        "native_constants": {"node1": native_constant_one, "node2": native_constant_forty},
        "preRp_radius_checks": preRp_radius_checks,
        "axial_radius_recipe": axial_radius_recipe,
        "axial_radius_checks": axial_radius_checks,
        "transport_checks": source_transport_checks,
        "passed": native_kernel_ok and all(source_identity_checks.values()) and all(source_transport_checks.values())
        and ast_rows["canonical_constructor"].get("passed") is True
        and ast_rows["axial_rate_evaluate"].get("passed") is True
        and source_pins_ok,
        "source_executed": False,
    }
    if not source_change_of_variables["passed"]:
        gaps.append("cross-producer quiet source equality guards failed")
    source_mean_equality = {
        "passed": source_change_of_variables["passed"]
        and ast_rows["axial_pressure_fields"].get("passed") is True
        and ast_rows["compact_geometry"].get("passed") is True
        and ast_rows["quiet_geometry"].get("passed") is True
        and axial_radius_ok
        and preRp_radius_ok
        and native_sigma_source_ok
        and native_constant_one
        and native_constant_forty
        and ast_rows["axial_rate_default"].get("passed") is True
        and ast_rows["axial_pressure_memory"].get("passed") is True,
        "mstar_equals_F": source_change_of_variables["passed"],
        "Ktrue_equals_Kcan": source_change_of_variables["passed"],
        "R_axial_pairing": source_transport_checks["forward_mean_is_transportable"]
        and source_transport_checks["compact_radius_scale"]
        and source_transport_checks["quiet_radius_scale"]
        and axial_radius_ok
        and preRp_radius_ok
        and native_sigma_source_ok
        and native_constant_one
        and native_constant_forty
        and ast_rows["axial_rate_default"].get("passed") is True
        and ast_rows["axial_pressure_memory"].get("passed") is True,
        "source_change_of_variables": source_change_of_variables,
        "independent_saved_interval_difference_unresolved": True,
        "saved_decimal_mean_or_F_selected": False,
        "old_original_quiet_K_identity_used_metadata_used_as_proof": False,
    }

    canonical_quiet = {
        "quiet_pair": [31215, 0],
        "preRp_signed_pair": [31221, 31220],
        "ordinary_C1": 0,
        "handles": {"kappa": quiet_binding.get("kappa"), "Kcan": quiet_binding.get("K"),
                     "Rp": quiet_binding.get("upper"), "R2": quiet_binding.get("lower"),
                     "scale": 25203, "Lq": 31213, "frame_kappa": kappa_frame_node},
        "sum_chain": {
            "m15": {"node": 31221, "children": [31219, 31215], "active_plus_compact": 31219,
                    "quiet": 31215, "active_plus_compact_children": [30707, 31155]},
            "C1": {"node": 31220, "children": [30708, 31156]},
        },
        "scalar_bindings": {
            "Kcan_node": quiet_binding.get("K"),
            "kappa_node": quiet_binding.get("kappa"),
            "Rp_node": quiet_binding.get("upper"),
            "R2_node": quiet_binding.get("lower"),
            "scale_node": 25203,
            "Lq_node": 31213,
            "Lq_expression": "log(Rp/(2Rc))",
            "kappa_source": {
                "report": _relative(paths["native_frame"], repo),
                "path": "actual_Rp_native_frame_functions.actual_identified_Rp_histories_C3.m[1]",
                "graph_node": kappa_frame_node,
            },
            "canonical_log_F": {
                "report": _relative(paths["quiet_source"], repo),
                "path": "actual_quiet_source.kernel.log_F",
                "saved_record": _mapping(quiet_source.get("kernel")).get("log_F"),
                "selected_as_exact_value": False,
            },
            "canonical_quiet_length": {
                "report": _relative(paths["quiet_source"], repo),
                "path": "actual_quiet_source.original_quiet_length",
                "saved_record": quiet_source.get("original_quiet_length"),
                "selected_as_exact_value": False,
            },
            "canonical_K_log_coefficient": {
                "report": _relative(paths["quiet_source"], repo),
                "path": "actual_quiet_source.K.log_coefficient",
                "saved_record": _mapping(quiet_source.get("K")).get("log_coefficient"),
                "selected_as_exact_value": False,
            },
        },
        "graph_nodes": graph_proofs,
        "graph_checks": graph_checks,
        "quiet_source_identity": {
            "actual_kernel_node": quiet_source_identity.get("actual_kernel_node"),
            "original_defining_integrand": quiet_source_identity.get("original_defining_integrand"),
            "exact_cutoff": quiet_source_identity.get("exact_cutoff"),
            "physical_Rc_node": quiet_source_identity.get("physical_Rc_node"),
            "physical_K": quiet_source_identity.get("physical_K"),
            "exact_original_quiet_length_node": quiet_source_identity.get("exact_original_quiet_length_node"),
        },
        "saved_scalar_subgraph_only": True,
    }

    # The source graph gives the exact algebraic difference, while the
    # equality of the two K scales remains an unresolved theorem obligation.
    identity = {
        "formula": "Delta_quiet=-(Lq/4)*(Ktrue^2-Kcan^2)",
        "Ktrue": "4*R_axial_end*mstar_true",
        "Kcan": "Rp*scale*kappa",
        "Lq": "log(Rp/(2Rc))",
        "ordinary_C1": 0,
        "source_graph_pair": [31215, 0],
        "Ktrue_equals_Kcan_proved": source_mean_equality["passed"],
        "equality_source_theorem": "source_mean_equality",
        "source_mean_equality": source_mean_equality,
        "old_original_quiet_K_identity_used_metadata_is_not_proof": True,
        "enclosure_or_midpoint_selected_as_equality": False,
        "numerical_difference_unresolved": True,
    }

    forward_quiet = {
        "source_report": _relative(paths["compact_quiet_m15"], repo),
        "mean_summary": dict(mean_summary),
        "mean_source": mean_source,
        "mean_summary_path": "original_axial_m15_numeric.json.pressure_row_samples[native_coordinate=1].leading_scalars.mean",
        "anchor_coordinate": "1",
        "anchor_path": "OriginalCompactQuietM15.evaluate -> self.rows.axial.geometry('1')",
        "anchor_is_2Rc": False,
        "raw_term_count": 211,
        "reduced_term_count": 181,
        "discrepancy_view_term_count": 180,
        "reduced_count_source": "parent joint discrepancy view 180 plus this quiet scalar slot",
        "target_consumed_once": _mapping(parent.get("contract")).get("independent_target_consumed_once") is True
        or all(_mapping(view).get("independent_target_consumed_once") is True
               for view in parent.get("joint_discrepancy_views", ()) if isinstance(view, Mapping)),
        "saved_incoming_term_count": 211,
        "saved_output_term_count": 212,
        "increment": "-4*R_axial_end^2*mstar_true^2*Lq",
        "Ktrue": "4*R_axial_end*mstar_true",
        "forward_wf": "4*mstar_true",
        "forward_Ktrue_factorization": "R_axial_end*forward_wf",
        "coefficient": {"sign": -4, "radial_power": 2, "mean_source": "saved mean_point/mean_truth enclosure"},
        "ordinary_C1": 0,
        "sample_receipts": forward_receipts,
        "physical_log_length_source": {
            "path": "samples[chart=quiet,native_coordinate=1].geometry.physical_local_log_distance",
            "selected_as_exact_value": False,
            "enclosure_only": True,
        },
        "original_quiet_log_scale_terms": {
            "path": "samples[chart=quiet,native_coordinate=1].m15_C0_C1_terms[-1].log_scale_terms",
            "by_Z": {key: value.get("log_scale_terms", []) for key, value in forward_receipts.items()},
            "selected_as_exact_value": False,
        },
        "original_2logRaxial_atoms": {
            "path": "samples[chart=quiet,native_coordinate=1].geometry.atomic_log_radius_terms",
            "by_Z": {key: value.get("atomic_log_radius_terms", []) for key, value in forward_receipts.items()},
            "parent_coefficient": 2,
            "selected_as_exact_value": False,
        },
        "source_function_binding": {
            "exact_primitive": "-4*R_axial_end^2*mstar_true^2*physical_local_log_length",
            "common_source_function": "original_positive_axial_terminal_mean",
            "half_S1_already_included": True,
            "global_radius_and_global_mean_paired": True,
        },
        "transport_source": {
            "compact_mean": "mean_source*exp(-14)/x",
            "quiet_mean": "mean_source*exp(-14)/2*exp(-s)",
            "compact_radius": "R_axial_end*exp(14)*x",
            "quiet_radius": "R_axial_end*exp(14)*2*exp(s)",
            "Ktrue_after_pairing": "4*R_axial_end*mean_source",
            "axial_log_radius_at_one": "log(110)+10*logCstar+10*logP+1+expm1(40)",
            "canonical_log_Rc": "log(110)+10*logCstar+11*exp(40)+124",
            "canonical_Rc_pair": "Rc=R_axial_end*exp(14)",
            "reference_logCstar_source": "OriginalReferenceNumeric.logC <- OriginalCoreLogAmplitude.logC",
            "reference_logP_source": "OriginalReferenceNumeric.logP <- OriginalCoreLogAmplitude.lp",
            "source_radius_checks": {"preRp": preRp_radius_checks, "axial": axial_radius_checks},
            "source_AST": {
                "axial_pressure_binding": ast_rows["axial_pressure_binding"],
                "axial_pressure_memory": ast_rows["axial_pressure_memory"],
                "axial_pressure_fields": ast_rows["axial_pressure_fields"],
                "axial_geometry": ast_rows["axial_geometry"],
                "axial_rate_default": ast_rows["axial_rate_default"],
                "compact_transport": ast_rows["compact_transport"],
                "quiet_transport": ast_rows["quiet_transport"],
                "compact_geometry": ast_rows["compact_geometry"],
                "quiet_geometry": ast_rows["quiet_geometry"],
            },
            "source_only": True,
        },
        "physical_scale_materialized": False,
        "source_only": True,
    }

    physical_scale_guards = {
        "forward_anchor_coordinate_one": ast_rows["forward_evaluate"].get("passed") is True,
        "mean_selector_source_only": ast_rows["mean_selection"].get("passed") is True,
        "quiet_source_Z1_exact_zero": quiet_source_guards["source_scaled_Z1_exact_zero"],
        "canonical_ordinary_C1_zero": graph_checks["ordinary_quiet_C1_zero"],
        "canonical_Lq_graph_bound": graph_proofs.get("31213", {}).get("passed") is True,
        "source_mean_equality": source_mean_equality["passed"],
        "no_physical_scale_materialization": True,
    }

    transitive_ok = input_passed == input_total and impl_passed == impl_total
    passed = bool(
        parent_pin.get("passed") and parent_pin.get("mapping_match")
        and transitive_ok and all(family_checks.values()) and ast_passed
        and mean_equal and forward_ok and all(quiet_source_guards.values())
        and graph_ok and all(physical_scale_guards.values())
        and source_pins_ok and not gaps
    )
    return {
        "passed": passed,
        "source_only": True,
        "source_executed": False,
        "parent_report": _PARENT_RELATIVE,
        "parent_report_sha256": _PARENT_SHA256,
        "parent_report_pin": parent_pin,
        "source_family": dict(_FAMILY),
        "input_hashes": selected_inputs,
        "implementation_hashes": selected_impl,
        "source_pins": source_pins,
        "source_pins_ok": source_pins_ok,
        "added_pins": {},
        "transitive_hash_validation": {
            "passed": transitive_ok,
            "input_passed": input_passed,
            "input_total": input_total,
            "implementation_passed": impl_passed,
            "implementation_total": impl_total,
        },
        "AST": ast_rows,
        "ast_function_hashes": {key: row.get("function_AST_sha256") for key, row in ast_rows.items()},
        "graph_identity": canonical_quiet,
        "canonical_quiet": canonical_quiet,
        "canonical": canonical_quiet,
        "forward_quiet": forward_quiet,
        "forward": forward_quiet,
        "mean_source": mean_source,
        "mean_source_summary": mean_source,
        "quiet_identity": identity,
        "quiet_scalar_identity": identity,
        "source_mean_equality": source_mean_equality,
        "source_change_of_variables": source_change_of_variables,
        "native_source_graph": native_kernel_proofs,
        "quiet_source_identity": quiet_source_guards,
        "physical_scale_guards": physical_scale_guards,
        "family_checks": family_checks,
        "numeric_value_available": False,
        "numeric_value_unavailable_reason": "saved producer enclosures are not function values; no scalar arithmetic performed",
        "limits": {
            "Ktrue_equals_Kcan": source_mean_equality["passed"],
            "mean_point_is_exact": False,
            "mean_truth_is_exact_function_value": False,
            "cap_or_midpoint_selected": False,
            "physical_scale_materialized": False,
            "accepted_computation_replayed": False,
        },
        "root": str(repo),
        "directory": str(work),
        "gaps": gaps,
    }


__all__ = ["build_quiet_scalar_source"]
