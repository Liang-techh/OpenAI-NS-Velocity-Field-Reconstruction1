"""Source-bound native second-micro to frozen-macro seam witness.

This module records the endpoint theorem needed before the native lower port
can be coupled to the frozen-macro source port.  It reads the already saved
macro-weighted report and its inherited physical seam, then parses the four
source locations that define the two micro prefixes, the macro endpoint
kernels, and the pressure recovery map.  It deliberately does not import or
run any owner, producer, constructor, checker, quadrature, or field evaluator.

The result is a formal coefficient/port identity.  It does not select a
saved interval endpoint or identify native transported ``Uy`` with the
macro ``source_Uy`` port.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_modal_source_point import _function_ast_guard


LATEST_REPORT = (
    "lei_ren_part1_paper_compliant_current_original_macro_weighted_m15_numeric.json"
)
LATEST_REPORT_SHA256 = "6835984804c4db15884b1832b4d95b16ccf1e394b2fd9145f612509164679cb7"


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except (OSError, ValueError):
        return None


def _digest(value: Any) -> str:
    """Use the repository's saved-witness JSON digest convention."""

    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode(
            "utf-8"
        )
    ).hexdigest()


def _repo(root: Any, directory: Any) -> Path:
    for candidate in (root, directory):
        if candidate is None:
            continue
        path = Path(candidate)
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _directory(root: Path, directory: Any) -> Path:
    if directory is None:
        return root / "experiments/root_st073"
    path = Path(directory)
    if path.is_file():
        path = path.parent
    return path


def _load_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _relative(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def _guard(
    path: Path,
    root: Path,
    function: str,
    fragments: tuple[str, ...],
    *,
    class_name: str | None = None,
) -> dict[str, Any]:
    """Run the no-execution AST helper and attach a stable source path."""

    try:
        value = _function_ast_guard(path, function, fragments, class_name=class_name)
    except (OSError, SyntaxError, UnicodeError, ValueError, StopIteration) as exc:
        return {
            "function": function,
            "class": class_name,
            "required_AST": {fragment: False for fragment in fragments},
            "function_AST_sha256": None,
            "source_path": _relative(path, root),
            "source_executed": False,
            "passed": False,
            "error": str(exc),
        }
    value = dict(value)
    value["source_path"] = _relative(path, root)
    value["source_executed"] = False
    return value


def _branch_guard(
    path: Path,
    root: Path,
    function: str,
    test: str,
    *,
    class_name: str | None = None,
) -> dict[str, Any]:
    """Confirm an exact branch test exists, using AST only."""

    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        scope: ast.AST = tree
        if class_name is not None:
            scope = next(
                node
                for node in tree.body
                if isinstance(node, ast.ClassDef) and node.name == class_name
            )
        functions = [
            node
            for node in ast.walk(scope)
            if isinstance(node, ast.FunctionDef) and node.name == function
        ]
        if len(functions) != 1:
            raise ValueError(f"Unique pinned source function required: {path}:{function}")
        node = functions[0]
        wanted = ast.parse(test, mode="eval").body
        passed = any(
            isinstance(part, ast.If)
            and ast.dump(part.test, include_attributes=False)
            == ast.dump(wanted, include_attributes=False)
            for part in ast.walk(node)
        )
        return {
            "function": function,
            "class": class_name,
            "test": test,
            "function_AST_sha256": hashlib.sha256(
                ast.dump(node, include_attributes=False).encode("utf-8")
            ).hexdigest(),
            "source_path": _relative(path, root),
            "source_executed": False,
            "passed": passed,
        }
    except (OSError, SyntaxError, UnicodeError, ValueError, StopIteration) as exc:
        return {
            "function": function,
            "class": class_name,
            "test": test,
            "function_AST_sha256": None,
            "source_path": _relative(path, root),
            "source_executed": False,
            "passed": False,
            "error": str(exc),
        }


def _source_pin(
    path: Path,
    root: Path,
    *,
    expected: str | None,
    name: str,
    kind: str,
) -> dict[str, Any]:
    actual = _sha256(path)
    passed = actual is not None and expected is not None and actual == expected
    return {
        "name": name,
        "kind": kind,
        "path": _relative(path, root),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": passed,
        "source_executed": False,
    }


def _empty_failure(gaps: list[str], *, root: Path, directory: Path) -> dict[str, Any]:
    return {
        "passed": False,
        "input_hashes": {},
        "implementation_hashes": {},
        "AST": {},
        "source_pins": {},
        "same_physical_R": {"passed": False, "source_values_selected": False},
        "nonzero_epsilon_cancellation_formal": {
            "passed": False,
            "source_values_selected": False,
        },
        "normalized_ordinary_Z_orders_0_to_2": {
            "passed": False,
            "ordinary_Z_orders": [0, 1, 2],
            "source_values_selected": False,
        },
        "derivative_port_defect_identity": {
            "passed": False,
            "source_values_selected": False,
        },
        "source_values_selected": False,
        "source_executed": False,
        "root": str(root),
        "directory": str(directory),
        "gaps": gaps,
    }


def _scope_ok(scope: Any, chart: str, query: str) -> bool:
    return (
        isinstance(scope, Mapping)
        and scope.get("chart") == chart
        and str(scope.get("native_query_value")) == query
        and scope.get("ordinary_parameter") == "Z"
        and all(order in scope.get("ordinary_Z_orders", ()) for order in (0, 1, 2))
    )


def _orders(rows: Any, quantity: str) -> list[int]:
    if not isinstance(rows, list):
        return []
    return sorted(
        {
            int(row["ordinary_Z_order"])
            for row in rows
            if isinstance(row, Mapping)
            and row.get("quantity") == quantity
            and isinstance(row.get("ordinary_Z_order"), int)
        }
    )


def _all_passed(rows: Any) -> bool:
    return isinstance(rows, list) and bool(rows) and all(
        isinstance(row, Mapping) and row.get("passed") is True for row in rows
    )


def _source_value_not_selected(row: Mapping[str, Any]) -> bool:
    """Reject selected values at either a row or direct-identity level."""

    if row.get("source_values_selected") is True:
        return False
    identity = row.get("actual_function_identity")
    return not (isinstance(identity, Mapping) and identity.get("source_values_selected") is True)


def _leading_qy_rows_pass(rows: Any, left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    """Check the saved leading common-port rows, including endpoint scopes."""

    if not isinstance(rows, list) or not rows:
        return False
    for row in rows:
        if not isinstance(row, Mapping) or row.get("passed") is not True:
            return False
        if row.get("independent_of_modified_Qy") is not True or not _source_value_not_selected(row):
            return False
        identity = row.get("actual_function_identity")
        if not isinstance(identity, Mapping) or identity.get("passed") is not True:
            return False
        if identity.get("direct_scope_match") is not True or identity.get("direct_row_nodes_match") is not True:
            return False
        if identity.get("exact_endpoint_fingerprints") is not True:
            return False
        if row.get("left_scope") != left or row.get("right_scope") != right:
            return False
    return True


def _leading_source_rows_pass(rows: Any) -> bool:
    """Check direct leading V/m source rows without using modified traces."""

    if not isinstance(rows, list) or not rows:
        return False
    return all(
        isinstance(row, Mapping)
        and row.get("passed") is True
        and _source_value_not_selected(row)
        and isinstance(row.get("source_row_sha256"), str)
        and isinstance(row.get("left_node"), int)
        and isinstance(row.get("right_node"), int)
        for row in rows
    )


def _source_guards(repo: Path, parent: Mapping[str, Any]) -> tuple[dict[str, Any], list[tuple[str, Path, str, str]]]:
    """Bind the requested source fragments without importing source owners."""

    source = repo / "src/openai_ns_reconstruction"
    base = repo / "experiments/root_st073"
    specs: tuple[tuple[str, Path, str, str, tuple[str, ...]], ...] = (
        (
            "macro_increments",
            source / "original_macro_numeric.py",
            "OriginalMacroNumeric",
            "_macro_increments",
            (
                "with mp.workdps(280):angular,axial=o.bridge.contributions(prepared,'macro',q)",
                "row=J(iv,[o._scaled(v,term) for v in term['normalized_axial_coefficients'][:order+1]])",
                "row=J(iv,[o._scaled(v,term,divide_eps=True) for v in term['normalized_axial_coefficients'][:order+1]])",
            ),
        ),
        (
            "macro_background",
            source / "original_macro_numeric.py",
            "OriginalMacroNumeric",
            "_background",
            (
                "bp,micro,_=OriginalSecondMicroNumeric._background(self,Z,iv.mpf(2))",
                "truth['Mpsi']=micro['Mpsi']*theta+a['Psi']*mass+o._norm_error_jet(self._uniform_Psi_increment(Z,prepared)*mass,order)",
            ),
        ),
        (
            "macro_geometry",
            source / "original_macro_numeric.py",
            "OriginalMacroNumeric",
            "_geometry",
            (
                "if not q:\n    theta=iv.mpf(1)",
                "mass=1-theta",
            ),
        ),
        (
            "second_increments",
            source / "original_second_micro_numeric.py",
            "OriginalSecondMicroNumeric",
            "_second_increments",
            (
                "with mp.workdps(280):angular,axial=o.bridge.contributions(prepared,'second',s)",
                "row=J(iv,[o._scaled(v,term) for v in term['normalized_axial_coefficients'][:order+1]])",
                "row=J(iv,[o._scaled(v,term,divide_eps=True) for v in term['normalized_axial_coefficients'][:order+1]])",
            ),
        ),
        (
            "bridge_contributions",
            base / "lei_ren_part1_paper_compliant_actual_bridge_integrals.py",
            "CompliantActualBridgeIntegrals",
            "contributions",
            (
                "first_phase=s if chart=='first' else c.mpf(1)",
                "second_phase=s-1 if chart=='second' else c.mpf(1 if chart=='macro' else 0)",
                "length=p['L']*s",
                "R1=c.mpf(100) if endpoints(s)==(mp.mpf(1),mp.mpf(1)) else self.bridge.r*c.exp(self.Y*s+c.mpf([0,endpoints(2*self.cap*(1-s))[1]]))",
                "kernels1=[c.mpf(0)]*3 if endpoints(s)==(mp.mpf(0),mp.mpf(0)) else radius_kernels(c,p['R0'],R1,length,1)",
                "kernels=[c.mpf(0)]*3 if endpoints(s)==(mp.mpf(0),mp.mpf(0)) else radius_kernels(c,p['R0'],R1,length,powerR[name])",
            ),
        ),
        (
            "bridge_radius_kernels",
            base / "lei_ren_part1_paper_compliant_actual_bridge_integrals.py",
            "",
            "radius_kernels",
            (
                "if power==1:\n    return [radius1-radius0,radius0*length,radius0-radius0*radius0/radius1]",
                "if power==2:\n    return [(radius1*radius1-radius0*radius0)/2,radius0*radius1-radius0*radius0,radius0*radius0*length]",
            ),
        ),
        (
            "pressure_inputs",
            source / "original_macro_pressure_rows.py",
            "OriginalMacroPressureRows",
            "inputs",
            (
                "truth['Mpsi']=incoming['actual']['Mpsi']*geo['theta']+a['Psi']*geo['mass']+o._norm_error_jet(self.micro._uniform_Psi_increment(Z,prepared)*geo['mass'],order)",
                "uy=self.source_Uy(prepared,geo,order)",
            ),
        ),
        (
            "pressure_background",
            source / "original_macro_pressure_rows.py",
            "OriginalMacroPressureRows",
            "background",
            (
                "d,L=1-z*z,1-z*z*delta",
                "qpsi=(psi*(2*z)-mean*(z*(1-delta))-self.o._derivative(mean)*d)/L",
                "my=psi-mean",
                "qpsiy=(psiy*(2*z)-my*(z*(1-delta))-self.o._derivative(my)*d)/L",
                "u0+psi*eps",
                "u0+mean*eps",
            ),
        ),
    )
    ast_rows: dict[str, Any] = {}
    paths: list[tuple[str, Path, str, str]] = []
    for name, path, cls, function, fragments in specs:
        ast_rows[name] = _guard(
            path, repo, function, fragments, class_name=cls or None
        )
        paths.append((name, path, "implementation" if path.is_relative_to(source) else "input", function))
    # The helper above checks the exact kernel branches.  This separate branch
    # witness confirms that those branches are nested under chart == 'macro'.
    ast_rows["bridge_macro_branch"] = _branch_guard(
        base / "lei_ren_part1_paper_compliant_actual_bridge_integrals.py",
        repo,
        "contributions",
        "chart == 'macro'",
        class_name="CompliantActualBridgeIntegrals",
    )
    return ast_rows, paths


def build_native_macro_seam(
    parent_report: Any, *, root: Any, directory: Any
) -> dict[str, Any]:
    """Build a compact source-bound native ``s=2`` / macro ``q=0`` witness.

    ``parent_report`` is normally the loaded latest macro-weighted report.
    A path is accepted only as a convenience for read-only loading.  Any
    missing, changed, or insufficient source witness returns ``passed=False``
    with concrete gaps and no synthesized rows.
    """

    repo = _repo(root, directory)
    base = _directory(repo, directory)
    gaps: list[str] = []
    if isinstance(parent_report, (str, Path)):
        parent_path = Path(parent_report)
        if not parent_path.is_absolute():
            parent_path = base / parent_path
        try:
            parent_report = _load_json(parent_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            return _empty_failure([f"cannot load parent report: {exc}"], root=repo, directory=base)
    if not isinstance(parent_report, Mapping):
        return _empty_failure(["loaded latest macro-weighted parent mapping is required"], root=repo, directory=base)

    latest_path = base / LATEST_REPORT
    latest_sha = _sha256(latest_path)
    latest_pin_ok = latest_sha == LATEST_REPORT_SHA256
    latest_report_digest: str | None = None
    parent_report_digest: str | None = None
    if not latest_pin_ok:
        gaps.append(
            f"latest parent report pin mismatch: expected {LATEST_REPORT_SHA256}, got {latest_sha}"
        )
    try:
        latest_mapping = _load_json(latest_path)
        latest_report_digest = _digest(latest_mapping)
        parent_report_digest = _digest(parent_report)
        if parent_report_digest != latest_report_digest:
            gaps.append("passed parent mapping does not match the pinned latest report mapping")
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        gaps.append(f"cannot bind parent mapping to latest saved report: {exc}")
    parent_checks = parent_report.get("checks")
    if not isinstance(parent_checks, Mapping) or parent_checks.get("passed") is not True:
        gaps.append("latest macro-weighted parent checks.passed is not true")
    if parent_report.get("source_primitive_only") is not True:
        gaps.append("parent source_primitive_only witness is missing")
    if parent_report.get("signed183_views_unchanged") is not True:
        gaps.append("parent signed183_views_unchanged witness is missing")

    input_hashes_parent = parent_report.get("input_hashes")
    implementation_hashes_parent = parent_report.get("implementation_hashes")
    if not isinstance(input_hashes_parent, Mapping):
        gaps.append("parent input_hashes mapping is missing")
        input_hashes_parent = {}
    if not isinstance(implementation_hashes_parent, Mapping):
        gaps.append("parent implementation_hashes mapping is missing")
        implementation_hashes_parent = {}

    input_hashes: dict[str, str | None] = {LATEST_REPORT: latest_sha}
    implementation_hashes: dict[str, str | None] = {}
    source_pins: dict[str, Any] = {}

    # The inherited endpoint ledger is the saved physical seam witness.  Its
    # bytes and the parent's input pin must both agree before reading it.
    packet = parent_report.get("inherited_m15_packet")
    ledger: Mapping[str, Any] | None = None
    seam: Mapping[str, Any] | None = None
    if not isinstance(packet, Mapping):
        gaps.append("inherited_m15_packet mapping is missing")
    else:
        ledger_name = packet.get("report")
        ledger_expected = packet.get("report_sha256")
        if not isinstance(ledger_name, str) or not isinstance(ledger_expected, str):
            gaps.append("inherited ledger report name and SHA are required")
        else:
            ledger_path = base / ledger_name
            ledger_sha = _sha256(ledger_path)
            input_hashes[ledger_name] = ledger_sha
            source_pins["inherited_endpoint_ledger"] = _source_pin(
                ledger_path,
                repo,
                expected=ledger_expected,
                name=ledger_name,
                kind="input",
            )
            if ledger_sha != ledger_expected:
                gaps.append(
                    f"inherited ledger pin mismatch: expected {ledger_expected}, got {ledger_sha}"
                )
            if input_hashes_parent.get(ledger_name) != ledger_expected:
                gaps.append("parent input_hashes does not retain the inherited ledger pin")
            try:
                ledger = _load_json(ledger_path)
                seam = ledger["source_seam_proofs"]["joins"]["micro_macro"]
            except (KeyError, OSError, ValueError, json.JSONDecodeError, TypeError) as exc:
                gaps.append(f"inherited physical micro_macro seam is unavailable: {exc}")

    partial = parent_report.get("native_macro0_partial_seam_evidence")
    if not isinstance(partial, Mapping):
        gaps.append("parent native_macro0_partial_seam_evidence is missing")
        partial = {}

    seam_digest_ok = False
    scope_ok = False
    qy_ok = False
    leading_source_rows: Any = []
    leading_qy_rows: Any = []
    leading_source_ok = False
    leading_qy_ok = False
    radius_identity: Mapping[str, Any] = {}
    pressure_identity: Mapping[str, Any] = {}
    if isinstance(seam, Mapping):
        seam_digest_ok = partial.get("saved_physical_seam_digest") == _digest(seam)
        if not seam_digest_ok:
            gaps.append("saved physical seam digest does not match inherited ledger")
        if partial.get("saved_physical_seam_passed") is not True or seam.get("passed") is not True:
            gaps.append("saved physical micro_macro seam is not passed")
        if partial.get("source_values_selected") is not False or seam.get("source_values_selected") is not False:
            gaps.append("physical seam selected a source value")
        if seam.get("numeric_caps_used") is not False:
            gaps.append("physical seam used numeric caps")
        scope_ok = (
            _scope_ok(seam.get("left_scope"), "second_micro", "2")
            and _scope_ok(seam.get("right_scope"), "frozen_macro", "0")
            and partial.get("left_scope") == seam.get("left_scope")
            and partial.get("right_scope") == seam.get("right_scope")
        )
        if not scope_ok:
            gaps.append("native s=2 and frozen-macro q=0 endpoint scopes are not pinned")
        qy = seam.get("leading_Qy_common")
        leading_source_rows = seam.get("leading_source_rows", [])
        leading_qy_rows = qy.get("input_rows", []) if isinstance(qy, Mapping) else []
        leading_source_ok = _leading_source_rows_pass(leading_source_rows)
        leading_qy_ok = _leading_qy_rows_pass(
            leading_qy_rows,
            seam.get("left_scope", {}),
            seam.get("right_scope", {}),
        )
        qy_ok = (
            isinstance(qy, Mapping)
            and qy.get("passed") is True
            and qy.get("source_values_selected") is False
            and qy.get("leading_output_alias_assigned") is False
            and leading_qy_ok
            and partial.get("common_Qy_operator") == qy
        )
        if not qy_ok:
            gaps.append("saved common leading Qy operator witness is incomplete")
        if not leading_source_ok:
            gaps.append("leading V/m source rows are incomplete or selected")
        radius_identity = seam.get("radius_phase_identity", {})
        pressure_identity = seam.get("pressure_algebra_witness", {})

    ast_rows: dict[str, Any] = {}
    try:
        ast_rows, source_paths = _source_guards(repo, parent_report)
    except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
        source_paths = []
        gaps.append(f"source AST inventory failed closed: {exc}")
    for name, path, kind, _function in source_paths:
        relative = _relative(path, repo)
        lookup = relative if kind == "implementation" else path.name
        expected = (
            implementation_hashes_parent.get(lookup)
            if kind == "implementation"
            else input_hashes_parent.get(lookup)
        )
        pin = _source_pin(path, repo, expected=expected, name=lookup, kind=kind)
        source_pins[name] = pin
        if kind == "implementation":
            implementation_hashes[relative] = pin["sha256"]
        else:
            input_hashes[lookup] = pin["sha256"]
        if expected is None:
            gaps.append(f"parent source pin missing for {lookup}")
        elif pin["sha256"] != expected:
            gaps.append(f"source pin mismatch for {lookup}: expected {expected}, got {pin['sha256']}")

    ast_passed = all(
        isinstance(row, Mapping) and row.get("passed") is True for row in ast_rows.values()
    ) and bool(ast_rows)
    if not ast_passed:
        gaps.extend(
            f"AST source guard failed: {name}"
            for name, row in ast_rows.items()
            if not isinstance(row, Mapping) or row.get("passed") is not True
        )

    # The inherited seam proves the common physical radius and keeps the
    # endpoint factor R attached exactly once in the D rows.
    radius_rows = radius_identity.get("rows", []) if isinstance(radius_identity, Mapping) else []
    defect_rows = seam.get("D_equality_rows", []) if isinstance(seam, Mapping) else []
    same_physical_R_passed = (
        scope_ok
        and isinstance(radius_identity, Mapping)
        and radius_identity.get("passed") is True
        and radius_identity.get("common_radius_source_identity_proved") is True
        and radius_identity.get("global_phase_transfer") is True
        and {row.get("role") for row in radius_rows if isinstance(row, Mapping)}
        >= {"radius", "global_offset", "actual_phase"}
        and all(
            row.get("left_source_node_present") is True
            and row.get("right_source_node_present") is True
            for row in radius_rows
            if isinstance(row, Mapping)
        )
        and _all_passed(defect_rows)
        and all(row.get("physical_R_once") is True for row in defect_rows)
        and all(row.get("signed_endpoint_value_evaluated") is False for row in defect_rows)
    )
    if not same_physical_R_passed:
        gaps.append("same-physical-radius source identity is incomplete")
    same_physical_R = {
        "passed": same_physical_R_passed,
        "left_scope": seam.get("left_scope") if isinstance(seam, Mapping) else None,
        "right_scope": seam.get("right_scope") if isinstance(seam, Mapping) else None,
        "radius_source_identity_sha256": seam.get("radius_source_identity_sha256") if isinstance(seam, Mapping) else None,
        "radius_phase_identity_digest": _digest(radius_identity) if isinstance(radius_identity, Mapping) else None,
        "roles": ["radius", "global_offset", "actual_phase"],
        "physical_R_once_in_defect": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "signed_endpoint_value_evaluated": False,
    }

    # U and m are affine in epsilon times the normalized Psi/Mpsi ports.  The
    # cancellation is stated symbolically under epsilon != 0; no epsilon or
    # source coefficient is materialized here.
    leading_port_orders = {
        "leading_source_V": _orders(leading_source_rows, "V"),
        "leading_source_m": _orders(leading_source_rows, "m"),
        "leading_source_Vy": _orders(leading_source_rows, "Vy"),
        "leading_Qy_input_V": _orders(leading_qy_rows, "V"),
        "leading_Qy_input_m": _orders(leading_qy_rows, "m"),
        "leading_Qy_input_V_y": _orders(leading_qy_rows, "V_y"),
    }
    normalized_orders = {
        "Psi": [0, 1, 2],
        "Mpsi": [0, 1, 2],
        "ordinary_Z_orders": [0, 1, 2],
        "leading_common_port_orders": leading_port_orders,
        "leading_source_rows_passed": leading_source_ok,
        "leading_Qy_input_rows_passed": leading_qy_ok,
        "coefficient_convention": "ordinary rows carry D_Z^k; Taylor coefficients carry D_Z^k/k!; k=0,1,2",
        "native_scope": seam.get("left_scope") if isinstance(seam, Mapping) else None,
        "macro0_scope": seam.get("right_scope") if isinstance(seam, Mapping) else None,
        "actual_leading_Uy_ports_need_not_match": True,
        "source_values_selected": False,
    }
    zero_endpoint_kernel = {
        "passed": bool(
            isinstance(ast_rows.get("bridge_contributions"), Mapping)
            and ast_rows["bridge_contributions"].get("required_AST", {}).get(
                "kernels1=[c.mpf(0)]*3 if endpoints(s)==(mp.mpf(0),mp.mpf(0)) else radius_kernels(c,p['R0'],R1,length,1)"
            ) is True
            and ast_rows["bridge_contributions"].get("required_AST", {}).get(
                "kernels=[c.mpf(0)]*3 if endpoints(s)==(mp.mpf(0),mp.mpf(0)) else radius_kernels(c,p['R0'],R1,length,powerR[name])"
            ) is True
        ),
        "condition": "macro endpoint s=q=0 uses endpoints(s)==(0,0)",
        "kernel_definition": "kernels1=kernels=[0,0,0] at q=0; otherwise radius_kernels(c,R0,R1,length,power)",
        "radius_kernels_at_zero_bypassed": True,
        "source_values_selected": False,
    }
    normalized_passed = (
        ast_passed
        and seam_digest_ok
        and scope_ok
        and qy_ok
        and zero_endpoint_kernel["passed"]
        and leading_source_ok
        and leading_qy_ok
        and leading_port_orders["leading_source_V"] == [0, 1, 2]
        and leading_port_orders["leading_source_m"] == [0, 1, 2]
        and leading_port_orders["leading_Qy_input_V"] == [0, 1, 2]
        and leading_port_orders["leading_Qy_input_m"] == [0, 1, 2]
        and same_physical_R_passed
    )
    normalized_orders["passed"] = normalized_passed
    if not normalized_passed:
        gaps.append("normalized Psi/Mpsi coefficient seam through ordinary Z2 is not source-bound")

    nonzero_epsilon = {
        "passed": normalized_passed,
        "premise": "epsilon != 0 (formal source parameter; no value selected)",
        "affine_ports": [
            "U=U0+epsilon*Psi",
            "m=U0+epsilon*Mpsi",
            "Uy=epsilon*Psi_y",
            "Qy=epsilon*Qpsi_y",
        ],
        "cancellation": "DeltaU_k=Delta m_k=0 with common U0 and epsilon!=0 imply DeltaPsi_k=DeltaMpsi_k=0 for k=0,1,2",
        "physical_to_normalized_direction": True,
        "epsilon_factor_cancelled_formally": True,
        "epsilon_materialized": False,
        "source_values_selected": False,
    }

    pressure_formula_ok = (
        isinstance(pressure_identity, Mapping)
        and pressure_identity.get("passed") is True
        and pressure_identity.get("numeric_caps_used") is False
        and pressure_identity.get("source_values_selected") is False
        and pressure_identity.get("method") == "exact rational polynomial ordinary-Z quotient differentiation"
    )
    derivative_passed = normalized_passed and pressure_formula_ok and same_physical_R_passed
    if not pressure_formula_ok:
        gaps.append("saved pressure Q/Qy algebra witness is incomplete")
    pressure_source_summary = {
        "passed": pressure_identity.get("passed") if isinstance(pressure_identity, Mapping) else False,
        "method": pressure_identity.get("method") if isinstance(pressure_identity, Mapping) else None,
        "formula_AST_sha256": pressure_identity.get("formula_AST_sha256") if isinstance(pressure_identity, Mapping) else None,
        "ordinary_Z_orders": sorted(
            {
                row.get("ordinary_Z_order")
                for row in (pressure_identity.get("rows", []) if isinstance(pressure_identity, Mapping) else [])
                if isinstance(row, Mapping) and isinstance(row.get("ordinary_Z_order"), int)
            }
        ),
        "quantities": sorted(
            {
                row.get("quantity")
                for row in (pressure_identity.get("rows", []) if isinstance(pressure_identity, Mapping) else [])
                if isinstance(row, Mapping) and isinstance(row.get("quantity"), str)
            }
        ),
        "reduced_difference_monomials_zero": all(
            row.get("reduced_difference_monomials") == 0
            for row in (pressure_identity.get("rows", []) if isinstance(pressure_identity, Mapping) else [])
            if isinstance(row, Mapping)
        ),
        "ordinary_derivative_product_rule": all(
            row.get("ordinary_derivative_product_rule") is True
            for row in (pressure_identity.get("rows", []) if isinstance(pressure_identity, Mapping) else [])
            if isinstance(row, Mapping)
        ),
        "positive_original_L_shared": all(
            row.get("positive_original_L_shared") is True
            for row in (pressure_identity.get("rows", []) if isinstance(pressure_identity, Mapping) else [])
            if isinstance(row, Mapping)
        ),
        "source_values_selected": False,
        "numeric_caps_used": False,
    }
    derivative_identity = {
        "passed": derivative_passed,
        "identity": "Hnative-Hmacro0=Enative-Emacro0",
        "defect_H": "H=R*Qy",
        "defect_E": "E=2*Z*R*Uy/(1-delta*Z^2)",
        "Qy": "[2*Z*Uy-(1-delta)*Z*(U-m)-(1-Z^2)*(Uz-mz)]/(1-delta*Z^2)",
        "ordinary_Z_orders": [0, 1],
        "requires_normalized_jet_orders": [0, 1, 2],
        "same_physical_R_required": True,
        "pressure_algebra_source": pressure_source_summary,
        "native_derivative_port": "second_micro native lower derivative at s=2",
        "macro_derivative_port": "frozen_macro source_Uy at q=0",
        "native_lower_equals_macro_source_Uy": False,
        "actual_port_value_equality_proved": False,
        "source_values_selected": False,
        "numeric_caps_used": False,
    }

    if not derivative_passed:
        gaps.append("derivative port defect identity lacks its normalized-jet/source premises")

    # Keep the output compact and auditable.  The copied scopes are graph
    # provenance, not evaluated source values.
    result = {
        "passed": not gaps,
        "input_hashes": input_hashes,
        "implementation_hashes": implementation_hashes,
        "AST": ast_rows,
        "source_pins": source_pins,
        "parent_report": LATEST_REPORT,
        "parent_report_sha256": LATEST_REPORT_SHA256,
        "parent_report_digest": parent_report_digest,
        "latest_report_digest": latest_report_digest,
        "inherited_ledger": packet.get("report") if isinstance(packet, Mapping) else None,
        "inherited_ledger_sha256": packet.get("report_sha256") if isinstance(packet, Mapping) else None,
        "saved_physical_seam_digest": partial.get("saved_physical_seam_digest"),
        "same_physical_R": same_physical_R,
        "nonzero_epsilon_cancellation_formal": nonzero_epsilon,
        "normalized_ordinary_Z_orders_0_to_2": normalized_orders,
        "macro_q0_zero_endpoint_kernel": zero_endpoint_kernel,
        "derivative_port_defect_identity": derivative_identity,
        # Alias without the prose suffix for consumers that use the shorter
        # theorem name.
        "normalized_ordinary_Z_orders": normalized_orders,
        "derivative_port_defect": derivative_identity,
        "native_Uy_and_macro_source_Uy_are_distinct_ports": True,
        "source_values_selected": False,
        "source_executed": False,
        "accepted_owner_or_producer_called": False,
        "numerical_quadrature_called": False,
        "gaps": gaps,
    }
    return result


__all__ = ["build_native_macro_seam"]
