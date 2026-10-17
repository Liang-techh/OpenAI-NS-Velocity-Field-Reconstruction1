"""Source-bound native core axial direction identity.

This packet binds the native prescribed-shear direction to the same actual
core functions used by the canonical fixed point.  The native source uses
the physical combination

``Gtotal = R*Ghydro + R*Pstar^2*Gpressure + R^2*F0^2*Gswirl``.

On the common core, the accepted integrated stress-free equations imply
``Gtotal = -R*d_R(U_core)``.  The collar factor ``q_source`` is kept as a
separate contract: the full requested relation follows when the collar
source proves ``q_source == ratio``.  This module does not obtain that factor
from interval overlap or a selected numerical point.

Only saved reports, source ASTs, and fresh symbolic algebra are inspected.
No accepted producer, owner, checker, recurrence, norm calculation,
quadrature, or numerical source is executed.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import sympy as sp


PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_core_collar_pressure_identity.json"
)
PARENT_SHA256 = "ae91b0d72d0415e8146f7ee9926d1a35ddbedffff6b89eac1f26558dbb5b44c7"

COMMON_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json"
)
COMMON_REPORT_SHA256 = "a1dee74791b2f9c4a9cafcad62e573baf788c7b6f6b39cef7cde90d23a2ff672"

CANONICAL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
)
CANONICAL_SOURCE_SHA256 = "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59"
CANONICAL_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json"
)
CANONICAL_CHECK_SHA256 = "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead"

INNER_BRIDGE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_inner_bridge_profiles.py"
)
INNER_BRIDGE_SOURCE_SHA256 = "142945ac1dd8a5705c380e79da29ea6bd89374c18cfa50f0412f1ba70c520840"

BRIDGE_MIXED_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_bridge_mixed_C4.py"
)
BRIDGE_MIXED_SOURCE_SHA256 = "55569c8679beb43eaf4ba9681393a2604d34dfed1a83d66d53c6219cfa7e9187"

FIRST_INTERFACE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_first_interface_check.py"
)
FIRST_INTERFACE_SOURCE_SHA256 = "5bc09c79d5d0105bd42756da12ef5a94992b552a1413ef8e32f18f9d65075751"
FIRST_INTERFACE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_first_interface_check.json"
)
FIRST_INTERFACE_REPORT_SHA256 = "e3a631d85da4cbeb62b9f6bbd7345e1109349935185edab177d2252a9d010400"

FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(value, Mapping), f"Saved mapping required: {path}")
    return value


def _family(value: Any) -> dict[str, Any]:
    value = value if isinstance(value, Mapping) else {}
    nested = value.get("source_family")
    if isinstance(nested, Mapping):
        value = nested
    return {key: value.get(key) for key in FAMILY_KEYS}


def _all_true(value: Any) -> bool:
    return isinstance(value, Mapping) and bool(value) and all(
        item is True or (isinstance(item, Mapping) and _all_true(item))
        for item in value.values()
    )


def _digest(value: Any) -> str:
    def normalize(item: Any) -> Any:
        if item is None or isinstance(item, bool) or isinstance(item, str):
            return item
        if isinstance(item, int):
            return str(item)
        if isinstance(item, float):
            return repr(item)
        if isinstance(item, Mapping):
            return {str(key): normalize(item[key]) for key in sorted(item, key=str)}
        if isinstance(item, (list, tuple)):
            return [normalize(entry) for entry in item]
        return repr(item)

    payload = json.dumps(normalize(value), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _pin(root: Path, relative: str, expected: str) -> dict[str, Any]:
    actual = _sha(root / relative)
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
        "source_executed": False,
    }


def _find_function(tree: ast.AST, name: str, cls: str | None = None) -> ast.AST:
    scope: ast.AST = tree
    if cls is not None:
        classes = [node for node in ast.walk(tree)
                   if isinstance(node, ast.ClassDef) and node.name == cls]
        _require(len(classes) == 1, f"Unique class required: {cls}")
        scope = classes[0]
    functions = [node for node in ast.walk(scope)
                 if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                 and node.name == name]
    _require(len(functions) == 1, f"Unique function required: {name}")
    return functions[0]


def _ast_guard(root: Path, relative: str, name: str,
               fragments: tuple[str, ...], cls: str | None = None) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_path": str(root / relative),
        "function": name,
        "class": cls,
        "required_fragments": {fragment: False for fragment in fragments},
        "function_AST_sha256": None,
        "source_executed": False,
        "passed": False,
    }
    try:
        source = (root / relative).read_text(encoding="utf-8")
        fn = _find_function(ast.parse(source, filename=relative), name, cls)
        segment = ast.get_source_segment(source, fn) or ""
        compact = "".join(segment.split())
        row["required_fragments"] = {
            fragment: "".join(fragment.split()) in compact for fragment in fragments
        }
        row["function_AST_sha256"] = hashlib.sha256(
            ast.dump(fn, include_attributes=False).encode("utf-8")
        ).hexdigest()
        row["passed"] = all(row["required_fragments"].values())
    except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
        row["error"] = str(exc)
    return row


def _source_ast(root: Path) -> tuple[dict[str, Any], bool]:
    guards = {
        "native_direction": _ast_guard(
            root, INNER_BRIDGE_SOURCE, "direction", (
                "return dict(D_over_R=(-W+angular/(2*f))/L",
                "drive_hydro=(-W*v+(m-z*derivative(m))*((1-delta)/2)+(z*a)*(2*delta)-d*derivative(a))/(2*L)",
                "drive_pressure=((z*pressure)*(2*(1+delta))-d*derivative(pressure))/(2*L)",
                "drive_swirl=(-(z*b)*(2*delta)+d*derivative(b)+(z*p)*(2*(1+delta))-d*derivative(p))/(2*L)",
            ),
        ),
        "native_comparison_source": _ast_guard(
            root, BRIDGE_MIXED_SOURCE, "comparison_directions", (
                "f=true(phi,ratios);h=true(rows[MTH],ratios);k=true(rows[MTHZ],ratios)",
                "m=nested(rows[MZ]);a=nested(rows['A']);vv=nested(V);pressure=lift(p0)",
                "values=dict(D_over_R=(-W+angular/(f*2))/L",
            ),
        ),
        "native_physical_scales": _ast_guard(
            root, BRIDGE_MIXED_SOURCE, "evaluate", (
                "Dbar=[rate_rows(directions['D_over_R'],scale,k)*R for k in range(4)]",
                "drive=[rate_rows(directions['hydro'],scale,k)*R+algebra.shift(rate_rows(directions['pressure'],scale,k)*R,(0,1,0,0))+algebra.shift(rate_rows(directions['swirl'],scale*2,k)*R**2,(0,0,1,0)) for k in range(4)]",
            ),
            "CompliantBridgeMixedC4",
        ),
        "canonical_core_fields": _ast_guard(
            root, CANONICAL_SOURCE, "__init__", (
                "U=a.add(U0,a.scale(Psi,epsilon))",
                "self.functions=dict(Phi=Phi,Psi=Psi,mean_Psi=M,U0=U0",
                "Lambda,epsilon=g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda))",
            ),
            "CurrentOriginalN1CoreFixedPointFunctions",
        ),
        "integrated_core_axial_direction": _ast_guard(
            root, FIRST_INTERFACE_SOURCE, "integrated_core_equations_proof", (
                "amplitude=F0/Fbase;amplitude2=F0*F0/(Fbase*Fbase)",
                "D=directions['D_over_R']",
                "drive=directions['drive_hydro']+Pstar*Pstar*directions['drive_pressure']+r*Fbase*Fbase*directions['drive_swirl']",
                "exact(axial_integrated/r,2*L*(s.diff(V,r)+drive),'Original three-scale axial drive is not integrated stress-free core recovery')",
            ),
        ),
    }
    return guards, all(row.get("passed") is True for row in guards.values())


def _report_checks(root: Path, parent: Mapping[str, Any], family: Mapping[str, Any]) -> dict[str, Any]:
    common = _load(root / COMMON_REPORT)
    canonical = _load(root / CANONICAL_CHECK)
    interface = _load(root / FIRST_INTERFACE_REPORT)
    common_proof = common.get("actual_Banach_and_common_coefficient_proof", {})
    common_models = common.get("source_models_and_domains", {})
    residual = interface.get("original_integrated_stress_free_core_equations", {})
    checks = {
        "parent_report_passed": parent.get("checks", {}).get("passed") is True,
        "parent_scope_did_not_replay_ancestors": parent.get("checks", {}).get(
            "accepted_ancestors_replayed") is False,
        "common_family": _family(common) == dict(family),
        "canonical_family": _family(canonical.get("source_family")) == dict(family),
        "interface_family": _family(interface) == dict(family),
        "common_report_passed": common.get("all_passed") is True and
        common.get("all_scoped_checks_passed") is True,
        "common_fixed_point_admitted": common.get(
            "current_core_actual_operator_Banach_fixed_point_admitted") is True,
        "same_actual_source_graph": _all_true(
            common_models.get("same_current_original_source_graph", {})),
        "same_source_functions_not_boxes": common_proof.get(
            "same_source_analytic_functions_not_just_shared_receipt_ids") is True,
        "finite_rows_bound_to_same_solution": common.get(
            "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point") is True,
        "canonical_report_passed": canonical.get("all_passed") is True and canonical.get(
            "current_original_n1_true_core_fixed_point_functions_and_enclosure_bridge_installed") is True,
        "interface_report_passed": interface.get("all_passed") is True and interface.get(
            "all_scoped_checks_passed") is True,
        "integrated_residual_source_bound": residual.get("passed") is True and
        residual.get("original_core_equations_AST_bound") is True and
        residual.get("differentiated_integrated_angular_equation_identity") is True and
        residual.get("differentiated_integrated_axial_equation_identity") is True,
        "direction_and_drives_are_source_functions": residual.get(
            "integrated_core_means_and_direction_are_source_functions_not_caps") is True,
        "physical_scales_retained": residual.get(
            "pressure_and_swirl_physical_scales_not_dropped") is True,
        "ordinary_Z_orders_available": residual.get(
            "original_nested_direction_base_identities") == 4,
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "common_report": common,
        "canonical_report": canonical,
        "interface_report": interface,
    }


def _fresh_symbolic_identity() -> dict[str, Any]:
    R, rho, Z, Lambda, epsilon, Pstar = sp.symbols(
        "R rho Z Lambda epsilon Pstar", positive=True)
    # These are source functions, not values at a selected ordinary-Z point.
    # In particular, F0(Z) and all three physical drive components retain
    # their Z derivatives before the residual identities are differentiated.
    F0 = sp.Function("F0")(Z)
    Phi = sp.Function("Phi")(rho, Z)
    Psi = sp.Function("Psi")(rho, Z)
    Ghydro = sp.Function("Ghydro")(Z)
    Gpressure = sp.Function("Gpressure")(Z)
    Gswirl = sp.Function("Gswirl")(Z)
    ratio, q_source = sp.symbols("ratio q_source")
    native_drive = Ghydro + Pstar**2 * Gpressure + R * F0**2 * Gswirl
    Gtotal = R * Ghydro + R * Pstar**2 * Gpressure + R**2 * F0**2 * Gswirl
    scale_check = sp.expand(Gtotal - R * native_drive) == 0
    dR_logF = Lambda * sp.diff(Phi, rho) / Phi
    D_over_R = -2 * dR_logF
    Dbar = sp.expand(R * D_over_R)
    U_r = epsilon * Lambda * sp.diff(Psi, rho)
    rho_scaling = sp.simplify(
        Dbar.subs(R, rho / Lambda) + 2 * rho * sp.diff(Phi, rho) / Phi
    ) == 0
    # The integrated first-interface receipt supplies native_drive = -U_r.
    # Keep that source equation as an explicit alias here: differentiating it
    # through Z orders 0..3 proves the ordinary-Z residual without replaying
    # the accepted integrated proof.
    native_drive_alias = sp.Symbol("native_drive_alias")
    pde_residual = native_drive_alias + U_r
    source_residual_under_pde = sp.simplify(
        (R * pde_residual).subs(native_drive_alias, -U_r)
    )
    ordinary_Z_source_residuals = [
        sp.simplify(sp.diff(source_residual_under_pde, Z, k)) == 0
        for k in range(4)
    ]
    ordinary_Z_scale_residuals = [
        sp.simplify(sp.diff(Gtotal - R * native_drive, Z, k)) == 0
        for k in range(4)
    ]
    ordinary_Z_residuals = all(
        ordinary_Z_source_residuals + ordinary_Z_scale_residuals
    )
    amplitude_derivative_kept = all(
        sp.diff(term, Z).has(sp.Derivative(F0, Z))
        for term in (R**2 * F0**2 * Gswirl,)
    )
    drive_derivatives_kept = all(
        sp.diff(term, Z).has(sp.Derivative(component, Z))
        for term, component in (
            (R * Ghydro, Ghydro),
            (R * Pstar**2 * Gpressure, Gpressure),
            (R**2 * F0**2 * Gswirl, Gswirl),
        )
    )
    core_Gtotal = sp.expand(-R * U_r)
    drive_residual = sp.simplify(core_Gtotal + R * U_r) == 0
    normalized_U_r = sp.simplify(U_r.subs(epsilon, 1 / Lambda) - sp.diff(Psi, rho)) == 0
    full_residual = sp.expand(q_source * core_Gtotal + ratio * R * U_r)
    conditional_full = sp.simplify(full_residual.subs(q_source, ratio)) == 0
    unconstrained_full = sp.factor(full_residual)
    return {
        "passed": all((
            scale_check,
            drive_residual,
            rho_scaling,
            normalized_U_r,
            ordinary_Z_residuals,
            amplitude_derivative_kept,
            drive_derivatives_kept,
            conditional_full,
        )),
        "checks": {
            "Gtotal_equals_R_times_native_drive": scale_check,
            "same_core_Gtotal_plus_R_dR_U_zero": drive_residual,
            "rho_equals_Lambda_R_scaling": rho_scaling,
            "Lambda_epsilon_normalization": normalized_U_r,
            "ordinary_Z0_to_Z3_residuals": ordinary_Z_residuals,
            "F0_Z_derivative_retained": amplitude_derivative_kept,
            "physical_drive_Z_derivatives_retained": drive_derivatives_kept,
            "full_relation_when_q_source_equals_ratio": conditional_full,
        },
        "ordinary_Z_residual_proof": {
            "orders": [0, 1, 2, 3],
            "source_residuals_after_integrated_PDE_identity": ordinary_Z_source_residuals,
            "physical_scale_residuals": ordinary_Z_scale_residuals,
            "F0_and_drive_are_Z_functions": True,
            "epsilon_and_Lambda_are_Z_independent_symbols": True,
            "accepted_PDE_identity_consumed": "native_drive + d_R U_core = 0",
            "accepted_PDE_identity_replayed": False,
        },
        "canonical_handles": {
            "F_core": "F0(Z)*Phi(rho,Z), rho=Lambda*R",
            "U_core": "U0(Z)+epsilon*Psi(rho,Z)",
            "d_R_log_F_core": "Lambda*Phi_rho/Phi",
            "d_R_U_core": "epsilon*Lambda*Psi_rho=Psi_rho when Lambda*epsilon=1",
        },
        "native_scales": {
            "native_drive": "Ghydro+Pstar^2*Gpressure+R*F0^2*Gswirl",
            "Gtotal": "R*Ghydro+R*Pstar^2*Gpressure+R^2*F0^2*Gswirl",
        },
        "fresh_residual_without_qsource_contract": str(unconstrained_full),
        "source_executed": False,
    }


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Build the source-only native core axial-direction packet."""

    repo = Path(root)
    parent_path = repo / PARENT
    _require(_sha(parent_path) == PARENT_SHA256, "Pinned parent report changed")
    saved_parent = _load(parent_path)
    _require(_digest(parent) == _digest(saved_parent),
             "Supplied parent mapping differs from pinned report")
    family = _family(saved_parent.get("source_family"))
    _require(all(value is not None for value in family.values()),
             "Parent source family is incomplete")

    report_specs = {
        "parent": (PARENT, PARENT_SHA256),
        "common": (COMMON_REPORT, COMMON_REPORT_SHA256),
        "canonical_source": (CANONICAL_SOURCE, CANONICAL_SOURCE_SHA256),
        "canonical_check": (CANONICAL_CHECK, CANONICAL_CHECK_SHA256),
        "inner_bridge_source": (INNER_BRIDGE_SOURCE, INNER_BRIDGE_SOURCE_SHA256),
        "bridge_mixed_source": (BRIDGE_MIXED_SOURCE, BRIDGE_MIXED_SOURCE_SHA256),
        "first_interface_source": (FIRST_INTERFACE_SOURCE, FIRST_INTERFACE_SOURCE_SHA256),
        "first_interface_report": (FIRST_INTERFACE_REPORT, FIRST_INTERFACE_REPORT_SHA256),
    }
    source_pins = {name: _pin(repo, relative, expected)
                   for name, (relative, expected) in report_specs.items()}
    source_pin_passed = all(pin["passed"] for pin in source_pins.values())
    reports = _report_checks(repo, saved_parent, family)
    source_ast, source_ast_passed = _source_ast(repo)
    symbolic = _fresh_symbolic_identity()

    same_core_drive = bool(
        source_pin_passed and reports["passed"] and source_ast_passed and
        symbolic["checks"]["Gtotal_equals_R_times_native_drive"] and
        symbolic["checks"]["same_core_Gtotal_plus_R_dR_U_zero"] and
        symbolic["checks"]["Lambda_epsilon_normalization"]
    )
    required_contract = {
        "ordinary_Z_orders": [0, 1, 2, 3],
        "domain": "same actual common-core source function, not interval overlap",
        "angular": "D_over_R=-2*d_R log(F_core)=-2*Lambda*Phi_core_rho/Phi_core",
        "Dbar": "Dbar=R*D_over_R=-2*rho*Phi_core_rho/Phi_core",
        "axial_drive": "Gtotal=-R*d_R U_core",
        "Gtotal": "R*Ghydro+R*Pstar^2*Gpressure+R^2*F0^2*Gswirl",
        "q_source": "q_source=(phi_a/barphi)*exp(ell)=phi_evolved/barphi; phi_a is the core-exit inlet",
        "full_requested_relation": "q_source*Gtotal+ratio*R*d_R U_core=0",
        "full_relation_status": "conditional on the separate q_source=ratio collar transport identity",
        "q_source_ratio_proof_scope": "external collar transport contract; not inferred from interval comparison",
    }
    checks = {
        "parent_hash_and_mapping": True,
        "source_pins": source_pin_passed,
        "saved_reports_and_family": reports["passed"],
        "source_AST": source_ast_passed,
        "fresh_symbolic_scale_and_drive": symbolic["passed"],
        "same_core_axial_drive_identity": same_core_drive,
        "no_caps_as_function_values": reports["checks"]["same_source_functions_not_boxes"],
        "no_accepted_ancestor_replay": reports["checks"]["parent_scope_did_not_replay_ancestors"],
    }
    return {
        "passed": all(checks.values()),
        "source_only": True,
        "source_executed": False,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": family,
        "checks": checks,
        "source_pins": source_pins,
        "saved_report_checks": reports["checks"],
        "source_AST": source_ast,
        "fresh_symbolic_identity": symbolic,
        "actual_core_handles": {
            "F_core": "F0(Z)*Phi(rho,Z), rho=Lambda*R",
            "U_core": "U0(Z)+epsilon*Psi(rho,Z)",
            "Phi_Psi_mean_Psi": "canonical unique actual fixed-point components",
            "ordinary_Z_orders": [0, 1, 2, 3],
            "same_source_function_identity": reports["checks"]["same_source_functions_not_boxes"],
        },
        "axial_drive_identity": {
            "passed": same_core_drive,
            "identity": "Gtotal=-R*d_R U_core",
            "Gtotal": required_contract["Gtotal"],
            "native_drive": "Ghydro+Pstar^2*Gpressure+R*F0^2*Gswirl",
            "physical_scales_retained": reports["checks"]["physical_scales_retained"],
            "ordinary_Z_orders": [0, 1, 2, 3],
            "source_function_not_interval_comparison": True,
            "six_signed_terms_and_errors_untouched": True,
        },
        "required_contract": required_contract,
        "full_qsource_relation": {
            "formula": required_contract["full_requested_relation"],
            "conditional_symbolic_check": symbolic["checks"][
                "full_relation_when_q_source_equals_ratio"],
            "status": "conditional",
            "missing_input": "q_source=ratio for the native collar transport, through ordinary Z orders 0..3",
            "not_claimed_as_unconditional": True,
        },
        "remaining": [
            "The same-core native drive identity is source-bound and proved through the accepted integrated PDE residual.",
            "The full collar relation still requires the separate same-source transport proof q_source=ratio; no interval overlap or selected point is used here.",
            "All six signed terms, source errors, and physical scales remain retained; no force or correction is changed.",
        ],
    }
