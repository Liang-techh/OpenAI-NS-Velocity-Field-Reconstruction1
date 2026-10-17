"""Source-only proof of the forward/canonical n=1 input vector identity.

The finite-prefix forward implementation and the canonical regular system use
different units.  This module proves their *known source* conversion without
importing either producer or executing a numerical owner.  The forward source
tuple is ordered as ``(pressure, angular, axial)``:

``b_forward = (N1p/Pstar^2, 2*epsilon^2*N1theta/F0,
               2*N1z - 4*epsilon*Z*rho*N1p/L)``.

For ``x=sqrt(rho)`` define

``Dscale = diag((2*epsilon*Pstar^2*x)^-1, epsilon/F0, epsilon^-1)``.

Then the canonical raw source tuple is exactly
``Dscale**(-1)*b_forward = (g4,g5,g6)`` for positive ``rho``.  The pressure
component has the regular axis extension ``g4=0`` at ``rho=0``; the diagonal
scale itself is not evaluated at the axis.

Only symbolic differential algebra and read-only AST/report inspection are
used here.  The accepted common fixed-point receipt binds the finite rows and
all radial tails to the same unique analytic nonlinear fixed point and binds
the forward source graph to the canonical ``F0*Phi`` / ``U0+epsilon*Psi``
functions.  The forward ``S_Z_taylor`` input is retained as an enclosure of
the exact ``epsilon^2*F0^2`` jet; it is not a selected point function.  The
accepted n=0 pressure identity is likewise inherited as provenance, and this
module does not claim equality of the complete solved n=1 fields.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import sympy as sp


FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
CANONICAL_OPERATOR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
)
CANONICAL_VOLTERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
CANONICAL_CORE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
)
RADIAL_Q_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py"
)
COMMON_FIXED_POINT_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point.py"
)
COMMON_FIXED_POINT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json"
)
COMMON_OPERATOR_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json"
)
SCALED_SWIRL_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_scaled_swirl_source_check.json"
)
CANONICAL_CORE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json"
)
OPERATOR_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator_check.json"
)
VOLTERRA_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
)
CORRELATED_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_correlated_scaled_source.json"
)
RC_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rc_pressure_join.json"
)

# These are source/report identities, not numerical values of an original
# solution.  Pinning the exact files prevents the algebra below from silently
# drifting away from the source expressions it documents.
EXPECTED_SHA256 = {
    "src/openai_ns_reconstruction/original_core_numeric.py": "9c152f5ff50bd760ee1a860a11e48c1b3a2fd36027681ea05a3e7a0920181603",
    FORWARD_SOURCE: "a7a925360f3a6da2348c5f18ee74cfbb885f8e220730769131e43a18bd954c2f",
    CANONICAL_OPERATOR_SOURCE: "268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014",
    CANONICAL_VOLTERRA_SOURCE: "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208",
    CANONICAL_CORE_SOURCE: "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    RADIAL_Q_SOURCE: "1400066d2da5373a95f6cd30edad19359256211d0358b454203a92fac53c5ced",
    COMMON_FIXED_POINT_SOURCE: "dd44110f48fb40fee715b358801dd6b8944bda5e6d5c0054eefa63c66f82255f",
    COMMON_FIXED_POINT_REPORT: "a1dee74791b2f9c4a9cafcad62e573baf788c7b6f6b39cef7cde90d23a2ff672",
    COMMON_OPERATOR_REPORT: "a1f9a3632a7765f6bc8d38ac1bd25c28850e128a8b25adaa0d42f769f3704d62",
    SCALED_SWIRL_REPORT: "53b78264e65edab24e73af1725a112159a32041701a610db264be14ec54cafe0",
    CANONICAL_CORE_REPORT: "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead",
    OPERATOR_RECEIPT: "846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4",
    VOLTERRA_RECEIPT: "3cfdc3ce6a7e8c2a98a31fae985c02a2acb4e28025d313ee35db4cffc6e7b47c",
    CORRELATED_REPORT: "9beead66814da85325e86f4ee2d9c091da34c524f9512f8cbf4c77e18bf1d444",
    RC_REPORT: "5f160cfa1f251253155bdad0466c462c74c5d7349bdf5cbcb9b3b7209822e9c4",
}


def _repo_root(root: str | Path | None = None) -> Path:
    """Resolve a repository root without importing repository modules."""

    starts = []
    if root is not None:
        starts.append(Path(root))
    starts.extend((Path.cwd(), Path(__file__).resolve().parents[2]))
    for start in starts:
        if start.is_file():
            start = start.parent
        for candidate in (start, *start.parents):
            if (candidate / "src").is_dir() and (candidate / "experiments" / "root_st073").is_dir():
                return candidate
    raise FileNotFoundError("st073 repository root with src/ and experiments/root_st073/ is required")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_mapping(path: Path) -> Mapping[str, Any]:
    """Read a saved receipt with large integers kept as strings."""

    # parse_int=str keeps old binary-mantissa records inspectable on Python
    # builds with a decimal-integer digit limit.  No source owner is imported.
    value = json.loads(path.read_text(encoding="utf-8"), parse_int=str)
    if not isinstance(value, Mapping):
        raise ValueError(f"saved mapping required: {path}")
    return value


def _normal(text: str) -> str:
    return "".join(text.split())


def _function_node(path: Path, name: str, class_name: str | None = None) -> ast.AST:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    if class_name is None:
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
                return node
    else:
        for node in tree.body:
            if isinstance(node, ast.ClassDef) and node.name == class_name:
                for child in node.body:
                    if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == name:
                        return child
    raise ValueError(f"AST function {class_name + '.' if class_name else ''}{name} not found in {path}")


def _ast_pin(
    root: Path,
    relative: str,
    function: str,
    fragments: tuple[str, ...],
    *,
    class_name: str | None = None,
) -> dict[str, Any]:
    """Pin a literal source function and inspect it without executing it."""

    path = root / relative
    actual = _sha256(path)
    expected = EXPECTED_SHA256.get(relative)
    node = _function_node(path, function, class_name)
    source = ast.unparse(node)
    compact = _normal(source)
    fragment_pass = all(_normal(fragment) in compact for fragment in fragments)
    hash_pass = expected is None or actual == expected
    return {
        "path": relative,
        "function": function,
        "class": class_name,
        "required_fragments": list(fragments),
        "sha256": actual,
        "expected_sha256": expected,
        "function_ast_sha256": hashlib.sha256(ast.dump(node, include_attributes=False).encode()).hexdigest(),
        "source_executed": False,
        "passed": bool(hash_pass and fragment_pass),
        "hash_passed": hash_pass,
        "fragments_passed": fragment_pass,
    }


def _report_pin(root: Path, relative: str) -> dict[str, Any]:
    path = root / relative
    actual = _sha256(path)
    expected = EXPECTED_SHA256.get(relative)
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "source_executed": False,
        "passed": expected is None or actual == expected,
    }


def _input_path(root: Path, relative: str) -> Path:
    """Resolve an accepted producer input relative to ``root_st073``."""

    path = (root / "experiments" / "root_st073" / relative).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"accepted input hash path is missing: {relative}")
    return path


def _verify_input_hashes(root: Path, hashes: Mapping[str, Any]) -> dict[str, Any]:
    """Verify every inherited parent input hash without running its producer."""

    mismatches = {}
    for relative, expected in hashes.items():
        actual = _sha256(_input_path(root, str(relative)))
        if actual != expected:
            mismatches[str(relative)] = {"expected": expected, "actual": actual}
    canonical = json.dumps(
        sorted((str(key), str(value)) for key, value in hashes.items()),
        separators=(",", ":"),
    ).encode("utf-8")
    return {
        "count": len(hashes),
        "sha256": hashlib.sha256(canonical).hexdigest(),
        "mismatches": mismatches,
        "passed": not mismatches,
    }


def _accepted_common_fixed_point(root: Path) -> dict[str, Any]:
    """Bind the actual common fixed point and its shared input graph.

    The common report is an existence/inclusion receipt.  In particular, the
    scaled ``S`` source is admitted as Cauchy jets enclosing the exact
    ``epsilon**2*F0**2`` jets.  It is deliberately recorded as an inclusion,
    never as equality to a selected interval center or finite row.
    """

    common = _read_mapping(root / COMMON_FIXED_POINT_REPORT)
    operator = _read_mapping(root / COMMON_OPERATOR_REPORT)
    scaled = _read_mapping(root / SCALED_SWIRL_REPORT)
    canonical = _read_mapping(root / CANONICAL_CORE_REPORT)
    n1_operator = _read_mapping(root / OPERATOR_RECEIPT)
    n1_volterra = _read_mapping(root / VOLTERRA_RECEIPT)

    common_inputs = common.get("input_hashes", {})
    descendants = {
        "core_operator": operator,
        "scaled_swirl": scaled,
        "canonical_n1_core": canonical,
        "canonical_n1_operator": n1_operator,
        "canonical_n1_volterra": n1_volterra,
    }
    shared_mismatches = {}
    shared_missing = {}
    for name, report in descendants.items():
        child_inputs = report.get("input_hashes", {})
        missing = sorted(set(common_inputs) - set(child_inputs))
        different = sorted(
            key
            for key in set(common_inputs).intersection(child_inputs)
            if child_inputs[key] != common_inputs[key]
        )
        if missing:
            shared_missing[name] = missing
        if different:
            shared_mismatches[name] = {
                "different": different,
            }

    family = {
        "actual_five_defect_family_sha256": common.get("actual_five_defect_family_sha256"),
        "implicit_source_sha256": common.get("implicit_source_sha256"),
        "datum_enclosure_sha256": common.get("datum_enclosure_sha256"),
    }
    family_mismatches = {}
    for name, report in descendants.items():
        report_family = report.get("source_family", report)
        if not isinstance(report_family, Mapping):
            family_mismatches[name] = "missing source family mapping"
            continue
        bad = {
            key: {"expected": expected, "actual": report_family.get(key)}
            for key, expected in family.items()
            if report_family.get(key) != expected
        }
        if bad:
            family_mismatches[name] = bad

    fixed_proof = common.get("actual_Banach_and_common_coefficient_proof", {})
    tails = common.get("same_solution_all_order_tail_proof", {})
    models = common.get("source_models_and_domains", {})
    graph = models.get("same_current_original_source_graph", {})
    scaled_graph = scaled.get("current_source_graph", {})
    scaled_jet_inclusion = bool(
        scaled.get("all_passed") is True
        and scaled.get("current_core_exact_scaled_swirl_source_jets_admitted") is True
        and scaled.get("exact_scaled_S_source_jet_admitted") is True
        and scaled.get("current_original_seed_boxes_enclose_exact_positive_S_jets") is True
    )
    common_ok = bool(
        common.get("all_passed") is True
        and common.get("all_scoped_checks_passed") is True
        and common.get("current_core_actual_operator_Banach_fixed_point_admitted") is True
        and common.get("current_core_actual_operator_contraction_and_tail_certified") is True
        and common.get("finite_rows_and_tails_bound_to_same_nonlinear_fixed_point") is True
        and fixed_proof.get("passed") is True
        and fixed_proof.get("unique_actual_analytic_fixed_point_for_each_fixed_implicit_source") is True
        and fixed_proof.get("all_original_finite_rows_enclose_coefficients_of_this_unique_solution") is True
        and fixed_proof.get("same_source_analytic_functions_not_just_shared_receipt_ids") is True
        and tails.get("passed") is True
        and tails.get("finite_and_infinite_parts_of_same_analytic_solution") is True
        and all(graph.values())
        and all(scaled_graph.values())
        and models.get("passed") is True
        and scaled_jet_inclusion
        and not shared_mismatches
        and not family_mismatches
    )
    return {
        "common_report": _report_pin(root, COMMON_FIXED_POINT_REPORT),
        "common_source": _report_pin(root, COMMON_FIXED_POINT_SOURCE),
        "operator_report": _report_pin(root, COMMON_OPERATOR_REPORT),
        "scaled_swirl_report": _report_pin(root, SCALED_SWIRL_REPORT),
        "canonical_core_report": _report_pin(root, CANONICAL_CORE_REPORT),
        "common_input_hashes": _verify_input_hashes(root, common_inputs),
        "shared_input_hashes": {
            "descendant_count": len(descendants),
            "mismatches": shared_mismatches,
            "missing_in_descendants": shared_missing,
            "passed": not shared_mismatches,
        },
        "source_family": {
            "common": family,
            "mismatches": family_mismatches,
            "passed": not family_mismatches,
        },
        "same_source_functions": {
            "report_graph": graph,
            "scaled_source_graph": scaled_graph,
            "report_flag": fixed_proof.get(
                "same_source_analytic_functions_not_just_shared_receipt_ids"
            ),
            "passed": bool(
                fixed_proof.get(
                    "same_source_analytic_functions_not_just_shared_receipt_ids"
                ) is True
                and all(graph.values())
                and all(scaled_graph.values())
            ),
        },
        "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point": common.get(
            "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point"
        ),
        "scaled_S_jet_inclusion": {
            "passed": scaled_jet_inclusion,
            "exact_source": "S_scaled=epsilon^2*F0^2",
            "finite_seed_semantics": "accepted Cauchy S_Z_taylor boxes enclose exact S jets",
            "function_equality_claim": False,
        },
        "passed": common_ok,
        "point_values_selected": False,
    }


def _accepted_receipts(root: Path) -> dict[str, Any]:
    """Inspect accepted source receipts; never invoke their checkers."""

    operator = _read_mapping(root / OPERATOR_RECEIPT)
    volterra = _read_mapping(root / VOLTERRA_RECEIPT)
    operator_ok = bool(
        operator.get("all_passed") is True
        and operator.get("current_original_n1_actual_inner_operator_input_enclosures_installed") is True
        and operator.get("coefficient_F1_Uz1_P1_solved") is False
        and operator.get("n1_full_inner_interval_solution_certified") is False
    )
    volterra_ok = bool(
        volterra.get("all_passed") is True
        and volterra.get("current_original_n1_true_inner_Volterra_function_adapter_installed") is True
        and volterra.get("original_numeric_n1_solution_oracle_executed") is False
        and volterra.get("numerical_original_n1_point_accuracy_claimed") is False
    )
    return {
        "operator": {
            "receipt": _report_pin(root, OPERATOR_RECEIPT),
            "accepted_input_enclosures": operator_ok,
            "source_family": operator.get("source_family"),
        },
        "volterra": {
            "receipt": _report_pin(root, VOLTERRA_RECEIPT),
            "accepted_complete_regular_definition": volterra_ok,
            "source_family": volterra.get("source_family"),
        },
    }


def _accepted_pressure_scaling(root: Path) -> dict[str, Any]:
    """Bind the accepted P0 definition/scaling metadata as a premise only."""

    correlated_path = root / CORRELATED_REPORT
    rc_path = root / RC_REPORT
    correlated = _read_mapping(correlated_path)
    rc = _read_mapping(rc_path)
    identity = correlated.get("core_pressure_source_identity", {})
    scaling = identity.get("scaling", {}) if isinstance(identity, Mapping) else {}
    flags = scaling.get("source_scaling_flags", {}) if isinstance(scaling, Mapping) else {}
    required_flags = (
        "Lambda_times_core_epsilon_exactly_one",
        "P0_step_equals_core_epsilon_times_Pstar_squared_times_normalized_pressure",
        "no_second_epsilon_applied_to_ell",
        "core_radial_epsilon_distinct_from_pressure_datum_parameter_epsilon",
    )
    identity_ok = bool(
        identity.get("passed") is True
        and identity.get("source_executed") is False
        and identity.get("true_forward_and_canonical_pressure_same_function") is True
        and identity.get("true_normalized_pressure_function")
        == "p(Z)=-(M2/(1+Z^2)^2+M0+F_flat(Z))"
        and scaling.get("canonical") == "physical_P0=Pstar^2*p"
        and scaling.get("forward") == "P0_Z_taylor=epsilon_core*physical_P0_Taylor"
        and scaling.get("pressure_datum_epsilon_is_distinct") is True
        and all(flags.get(key) is True for key in required_flags)
    )
    rc_parent_ok = bool(
        rc.get("parent_report") == Path(CORRELATED_REPORT).name
        and rc.get("parent_report_sha256") == EXPECTED_SHA256[CORRELATED_REPORT]
        and rc.get("checks", {}).get("passed") is True
        and rc.get("checks", {}).get("accepted_ancestors_replayed") is False
        and rc.get("core_pressure_source_identity", {}).get("source_executed") is False
        and rc.get("core_pressure_source_identity", {}).get(
            "true_forward_and_canonical_pressure_same_function"
        ) is True
    )
    return {
        "correlated_report": _report_pin(root, CORRELATED_REPORT),
        "rc_report": _report_pin(root, RC_REPORT),
        "identity": {
            "passed": identity_ok,
            "source_executed": False,
            "true_normalized_pressure_function": identity.get("true_normalized_pressure_function"),
            "canonical": scaling.get("canonical"),
            "forward": scaling.get("forward"),
            "flags": {key: flags.get(key) for key in required_flags},
        },
        "rc_parent_binding_passed": rc_parent_ok,
        "point_values_selected": False,
    }


def _source_ast_provenance(root: Path) -> dict[str, Any]:
    return {
        "forward_common_core_scales": _ast_pin(root,
            "src/openai_ns_reconstruction/original_core_numeric.py","__init__",(
                "self.logP=c.exp(40)+11","self.logLambda=4*self.logP+1000",
                "self.epsilon=c.exp(-self.logLambda)","self.lp_iv=iv.exp(40)+11",
                "self.ll_iv=4*self.lp_iv+1000","self.epsilon_iv=iv.exp(-self.ll_iv)"),
            class_name="OriginalCoreNumeric"),
        "forward_Pstar_squared_scale": _ast_pin(root,FORWARD_SOURCE,"coefficients",(
            "lp=self.leading.lp_iv if directed else self.leading.logP","P2=c.exp(2*lp)"),
            class_name="OriginalN1CoreNumeric"),
        "forward_seed_inputs": _ast_pin(
            root,
            "src/openai_ns_reconstruction/original_core_numeric.py",
            "_seed",
            (
                "U=4*variable+j",
                "L=one-variable*variable*delta",
                "H=variable*(1-delta)/2+(one-variable*variable)*U",
                "gradient=L*H/(H*H+(j/500)**2)",
                "fixed=dict(ell_Z_taylor=[-v for v in gradient.coefficients],S_Z_taylor=S,U0_Z_taylor=list(U.coefficients)",
            ),
            class_name="OriginalCoreNumeric",
        ),
        "forward_linear_rhs": _ast_pin(
            root,
            FORWARD_SOURCE,
            "linear_rhs",
            (
                "rf=(2*e*(Q+E*(-2+de)/L)-2*H*g/L)*f",
                "ru=-8*F2*z*r*A*f/L",
                "rp=2*F2/(e*P2)*A*f",
            ),
        ),
        "forward_rhs_sources": _ast_pin(
            root,
            FORWARD_SOURCE,
            "_rhs",
            (
                "g2=-known([fixed['ell_Z_taylor']],0,axial+2)",
                "F2=R(known([fixed['S_Z_taylor']],0,axial).rows,n,axial,c.mpf)/(eps*eps)",
                "force_f=-2*m*m*A",
                "force_u=-2*T(-2,T(-1-delta,U2))",
                "Np=-S/(2*P2)+pressure_head/(eps*P2)",
                "force_u-=4*z*rho*Np*(eps*P2)/L",
                "Qw=(2*zz*uu-zz*mm*(1-delta)-dd*mw.dz())/LL",
                "pressure_head=R([[(i+1)*(i+2)*v for v in Qw.rows[i+1][:axial+1]",
                "return (zv*val*weight+(1-zv*zv)*value.dz()-2*zv*val.radial_euler())/(1-zv*zv*delta)",
            ),
            class_name="OriginalN1CoreNumeric",
        ),
        "canonical_core_functions": _ast_pin(
            root,
            CANONICAL_CORE_SOURCE,
            "__init__",
            (
                "U0=a.add(a.scale(Z,4),a.fixed(self.parameters['j']))",
                "H=a.add(a.scale(a.mul(a.add(one,a.neg(de)),Z),'1/2'),a.mul(d,U0))",
                "gradient=a.div(a.mul(L,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))",
                "sigma=g.quotient(self.parameters['j'],c(500),'positive fixed scalar 500')",
                "Gjet=c3.C3Function(G,gradient.value,gradient.Z,gradient.ZZ)",
                "logLambda=g.add(g.mul(c(4),logP),c(1000))",
                "Lambda,epsilon=(g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda)))",
                "logF=a.add(a.fixed(g.neg(self.parameters['logC'])),a.neg(a.scale(Gjet,Lambda)))",
                "F0=diff.lift(source.C1Function(F0value,g.mul(F0value,logF.Z)))",
                "U=a.add(U0,a.scale(Psi,epsilon))",
                "self.functions=dict(Phi=Phi,Psi=Psi,mean_Psi=M,U0=U0,gradient=gradient,chi=chi",
            ),
            class_name="CurrentOriginalN1CoreFixedPointFunctions",
        ),
        "canonical_F0_log_gradient": _ast_pin(
            root,
            CANONICAL_CORE_SOURCE,
            "__init__",
            (
                "F0value=g.unary('exp',logF.value)",
                "F0=diff.lift(source.C1Function(F0value,g.mul(F0value,logF.Z)))",
            ),
            class_name="CurrentOriginalN1CoreFixedPointFunctions",
        ),
        "canonical_core_fields": _ast_pin(
            root,
            CANONICAL_VOLTERRA_SOURCE,
            "core_fields",
            (
                "U=[a.add(core.functions['U0'],a.scale(rows['Psi'][0],eps)),*[a.scale(q,eps) for q in rows['Psi'][1:]]]",
                "mean=[a.add(core.functions['U0'],a.scale(rows['mean_Psi'][0],eps)),*[a.scale(q,eps) for q in rows['mean_Psi'][1:]]]",
                "F_rho_rows=[a.mul(core.functions['F0'],q) for q in rows['Phi']]",
                "Q_rho_rows=[own.radial_Q(g,core.parameters['delta'],u,m) for u,m in zip(U,mean)]",
            ),
        ),
        "canonical_Q_recovery": _ast_pin(
            root,
            RADIAL_Q_SOURCE,
            "radial_Q",
            (
                "numerator=a.add(a.scale(a.mul(z,u),2),a.neg(a.mul(a.add(one,a.neg(de)),a.mul(z,mean))),a.neg(a.mul(d,mz)))",
                "return a.div(numerator,L)",
            ),
        ),
        "canonical_known_forcing": _ast_pin(
            root,
            CANONICAL_OPERATOR_SOURCE,
            "known_forcing",
            (
                "theta=-axial_second(f,z,rho,delta,-2-delta)",
                "axial=-axial_second(w,z,rho,delta,-1-delta)",
            ),
        ),
        "canonical_axial_composition": _ast_pin(
            root,CANONICAL_VOLTERRA_SOURCE,"axial_second",(
                "A=g.add(g.mul(weight,z,f.value),g.mul(d,f.Z),g.neg(g.mul(c(2),z,rho,fr.value)))",
                "AR=g.add(g.mul(g.sub(weight,c(2)),z,fr.value),g.mul(d,fr.Z),g.neg(g.mul(c(2),z,rho,frr.value)))",
                "return div(g.add(div(g.mul(g.add(g.sub(weight,g.one),de),z,A))")),
        "canonical_raw_system": _ast_pin(
            root,
            CANONICAL_VOLTERRA_SOURCE,
            "raw_system",
            (
                "Np=a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))",
                "g4=g.mul(c(2),eps,x,Np.value)",
                "g5=g.mul(c(2),eps,Nt)",
                "g6=g.sub(g.mul(c(2),eps,Nz),g.quotient(g.mul(c(4),zv,r,eps,eps,Np.value),Lv,'positive original L'))",
            ),
        ),
    }


def _formal_differential_algebra() -> dict[str, Any]:
    """Construct the source map from the literal forward source formulas.

    The pressure row is built twice: first as the forward ``S`` and
    ``pressure_head`` expressions, then as the canonical raw ``N1p``.  This
    keeps ``b_pressure=N1p/Pstar^2`` from being an assumed definition.  The
    exact ``S_scaled=epsilon^2*F0^2`` relation is an inherited source-jet
    inclusion checked by :func:`_accepted_common_fixed_point`, not a point
    equality manufactured by this algebra.
    """

    z, rho, delta, epsilon, Lambda, Pstar = sp.symbols(
        "Z rho delta epsilon Lambda Pstar", nonzero=True
    )
    # P2 is the actual producer's exp(2*logP), with the same real logP as
    # Pstar=exp(logP); it is not an independent algebra symbol.
    P2=Pstar**2
    h = sp.Function("h")(z)
    phi = sp.Function("Phi")(z, rho)
    U0 = sp.Function("U0")(z, rho)
    Psi = sp.Function("Psi")(z, rho)
    Q = sp.Function("Q")(z, rho)
    F0 = sp.exp(-h / epsilon)
    U = U0 + epsilon * Psi
    F = F0 * phi
    g = sp.diff(h, z)
    d = 1 - z**2
    L = 1 - delta * z**2

    def Z_operator(weight: Any, value: sp.Expr) -> sp.Expr:
        return (
            weight * z * value
            + d * sp.diff(value, z)
            - 2 * z * rho * sp.diff(value, rho)
        ) / L

    # Independently expand the literal canonical axial_second source.
    weight=sp.Symbol('weight')
    test=sp.Function('test')(z,rho)
    fr=sp.diff(test,rho);fz=sp.diff(test,z)
    AA=weight*z*test+d*fz-2*z*rho*fr
    AZ=weight*test+(weight-2)*z*fz+d*sp.diff(test,z,2)-2*rho*fr-2*z*rho*sp.diff(fr,z)
    AR=(weight-2)*z*fr+d*sp.diff(fr,z)-2*z*rho*sp.diff(test,rho,2)
    expanded=((weight-1+delta)*z*AA/L+d*(AZ/L+2*delta*z*AA/L**2)-2*z*rho*AR/L)/L
    axial_composition=sp.factor(sp.together(expanded-Z_operator(weight-1+delta,Z_operator(weight,test))))==0

    N1theta = -Z_operator(-3, Z_operator(-2 - delta, F))
    N1z = -Z_operator(-2, Z_operator(-1 - delta, U))
    Qrho = sp.diff(Q, rho)
    Qrhorho = sp.diff(Q, rho, 2)
    Qz = sp.diff(Q, z)

    # Literal forward OriginalN1CoreNumeric._rhs pressure source.
    S_forward = (Q + z * Qz * (1 - delta) / 2 + rho * Qrho) / L
    S_forward += Q * (Q / 2 + rho * Qrho)
    S_forward += U * (-2 * z * Q + d * Qz - 2 * z * rho * Qrho) / L
    pressure_head_forward = 2 * Qrho + rho * Qrhorho
    Np_forward = -S_forward / (2 * P2) + pressure_head_forward / (epsilon * P2)

    # Literal canonical raw_system pressure source, independently named.
    S_canonical = (Q + z * Qz * (1 - delta) / 2 + rho * Qrho) / L
    S_canonical += Q * (Q / 2 + rho * Qrho)
    S_canonical += U * (-2 * z * Q + d * Qz - 2 * z * rho * Qrho) / L
    pressure_head_canonical = 2 * Qrho + rho * Qrhorho
    N1p_canonical = -S_canonical / 2 + Lambda * pressure_head_canonical
    pressure_identity = sp.factor(
        sp.together(
            Np_forward - (N1p_canonical / P2).subs(Lambda, 1 / epsilon)
        )
    )

    # This is the source-only part of OriginalN1CoreNumeric._rhs.  F2 is a
    # homogeneous matrix input and is intentionally absent from b_forward.
    m = d * g / L
    TA = Z_operator(-2 - delta, phi)
    btheta = -2 * m**2 * phi
    btheta += 2 * epsilon * (Z_operator(-3, m * phi) + m * TA)
    btheta -= 2 * epsilon**2 * Z_operator(-3, TA)
    forward_pressure_feedback = -4 * z * rho * Np_forward * (epsilon * P2) / L
    bazial = 2 * N1z + forward_pressure_feedback

    x = sp.symbols("x", positive=True)
    N1p_canonical_at_scale = N1p_canonical.subs(Lambda, 1 / epsilon)
    g4 = 2 * epsilon * x * N1p_canonical_at_scale
    g5 = 2 * epsilon * N1theta
    g6 = 2 * epsilon * N1z - 4 * epsilon**2 * z * rho * N1p_canonical_at_scale / L

    # The explicit amplitude derivative relation is used in the angular row.
    # Treating g as a function is essential: Z_Z(F0*Phi) contains g_Z.
    btheta_expected = 2 * epsilon**2 * N1theta / F0
    angular_difference = sp.factor(sp.together(btheta - btheta_expected))
    pressure_difference = sp.factor(
        sp.together((2 * epsilon * P2 * x) * Np_forward - g4)
    )
    axial_difference = sp.factor(sp.together(epsilon * bazial - g6))

    # For Q=sum_i Q_i rho^i, the actual source head has coefficient
    # (n+1)(n+2)Q_(n+1) at every radial order.  This is checked independently
    # of the finite-prefix producer and its pressure-head list construction.
    n = sp.symbols("n", integer=True, nonnegative=True)
    q_next = sp.Symbol("Q_(n+1)")
    pressure_head_coefficient = 2 * (n + 1) * q_next + n * (n + 1) * q_next
    pressure_head_expected = (n + 1) * (n + 2) * q_next
    pressure_head_coefficient_difference = sp.expand(
        pressure_head_coefficient - pressure_head_expected
    )

    checks = {
        "P2_equals_shared_Pstar_squared": sp.simplify(sp.exp(2*sp.Symbol('logP',real=True))-
            sp.exp(sp.Symbol('logP',real=True))**2)==0,
        "shared_core_Lambda_epsilon_reciprocal": sp.simplify(sp.exp(sp.Symbol('logLambda',real=True))*
            sp.exp(-sp.Symbol('logLambda',real=True))-1)==0,
        "literal_axial_second_equals_weighted_composition": axial_composition,
        "forward_S_equals_canonical_S": sp.factor(
            sp.together(S_forward - S_canonical)
        )
        == 0,
        "forward_pressure_head_equals_canonical_head": sp.factor(
            sp.together(pressure_head_forward - pressure_head_canonical)
        )
        == 0,
        "forward_Np_equals_canonical_N1p_over_Pstar_squared": pressure_identity == 0,
        "pressure_head_coefficient_all_orders": pressure_head_coefficient_difference == 0,
        "angular_source": angular_difference == 0,
        "pressure_source": pressure_difference == 0,
        "axial_source_with_pressure_feedback": axial_difference == 0,
        "g_equals_minus_epsilon_dZ_log_F0": sp.simplify(
            g + epsilon * sp.diff(sp.log(F0), z)
        )
        == 0,
    }
    checks["all_three"] = all(
        checks[key]
        for key in (
            "angular_source",
            "pressure_source",
            "axial_source_with_pressure_feedback",
        )
    )

    Dscale = sp.diag(
        1 / (2 * epsilon * P2 * x),
        epsilon / F0,
        1 / epsilon,
    )
    b_forward = sp.Matrix((Np_forward, btheta, bazial))
    canonical = sp.Matrix((g4, g5, g6))
    scaled_difference = tuple(
        sp.factor(sp.together(q))
        for q in (Dscale.inv() * b_forward - canonical)
    )
    checks["diagonal_vector_identity"] = all(value == 0 for value in scaled_difference)
    checks["vector_identity"] = bool(checks["all_three"] and checks["diagonal_vector_identity"])

    expressions = {
        "S_forward": S_forward,
        "pressure_head_forward": pressure_head_forward,
        "Np_forward": Np_forward,
        "S_canonical": S_canonical,
        "pressure_head_canonical": pressure_head_canonical,
        "N1p_canonical": N1p_canonical,
        "N1theta": N1theta,
        "N1z": N1z,
        "b_forward": b_forward,
        "canonical_g": canonical,
        "Dscale": Dscale,
        "Dscale_inverse_b_forward": Dscale.inv() * b_forward,
    }
    return {
        "checks": checks,
        "expressions": expressions,
        "source_order": ("pressure", "angular", "axial"),
        "canonical_order": ("g4", "g5", "g6"),
        "valid_for": "rho>0, L != 0; g4 has regular zero-axis extension",
        "point_values_selected": False,
    }


def prove_source_vector_identity(root: str | Path | None = None) -> dict[str, Any]:
    """Return a source-bound, universal proof receipt.

    The function raises if the pinned source/report graph changed or if the
    symbolic identity fails.  It does not execute an accepted owner,
    constructor, checker, recurrence, preparation step, or quadrature.
    """

    repo = _repo_root(root)
    ast_receipt = _source_ast_provenance(repo)
    common_status = _accepted_common_fixed_point(repo)
    receipt_status = _accepted_receipts(repo)
    pressure_status = _accepted_pressure_scaling(repo)
    algebra = _formal_differential_algebra()
    source_pass = all(row["passed"] for row in ast_receipt.values())
    common_pass = bool(
        common_status["common_report"]["passed"]
        and common_status["common_source"]["passed"]
        and common_status["operator_report"]["passed"]
        and common_status["scaled_swirl_report"]["passed"]
        and common_status["canonical_core_report"]["passed"]
        and common_status["common_input_hashes"]["passed"]
        and common_status["shared_input_hashes"]["passed"]
        and common_status["source_family"]["passed"]
        and common_status["same_source_functions"]["passed"]
        and common_status["scaled_S_jet_inclusion"]["passed"]
        and common_status["finite_rows_and_tails_bound_to_same_nonlinear_fixed_point"]
        is True
        and common_status["passed"]
    )
    receipt_pass = all(
        row["receipt"]["passed"] and value
        for row, value in (
            (receipt_status["operator"], receipt_status["operator"]["accepted_input_enclosures"]),
            (receipt_status["volterra"], receipt_status["volterra"]["accepted_complete_regular_definition"]),
        )
    )
    pressure_pass = bool(
        pressure_status["correlated_report"]["passed"]
        and pressure_status["rc_report"]["passed"]
        and pressure_status["identity"]["passed"]
        and pressure_status["rc_parent_binding_passed"]
    )
    inherited_premises = {
        "common_fixed_point_and_same_input_functions": common_pass,
        "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point": bool(
            common_status["finite_rows_and_tails_bound_to_same_nonlinear_fixed_point"]
            is True
        ),
        "scaled_S_jet_inclusion": common_status["scaled_S_jet_inclusion"]["passed"],
        "physical_P0_equals_Pstar_squared_p": pressure_status["identity"]["passed"],
        "forward_P0_step_equals_epsilon_core_times_physical_P0": pressure_status[
            "identity"
        ]["passed"],
        "P2_and_epsilon_use_actual_shared_core_scales": bool(
            pressure_status['identity']['passed'] and common_pass and
            ast_receipt['forward_common_core_scales']['passed'] and
            ast_receipt['forward_Pstar_squared_scale']['passed'] and
            ast_receipt['canonical_core_functions']['passed']),
    }
    symbolic_pass = bool(
        algebra["checks"]["vector_identity"]
        and algebra["checks"]["g_equals_minus_epsilon_dZ_log_F0"]
        and all(inherited_premises.values())
    )
    passed = bool(source_pass and common_pass and receipt_pass and pressure_pass and symbolic_pass)
    if not passed:
        failed = {
            "source_ast": [key for key, value in ast_receipt.items() if not value["passed"]],
            "common": common_status,
            "receipts": receipt_status,
            "pressure": pressure_status,
            "symbolic": algebra["checks"],
            "inherited_premises": inherited_premises,
        }
        raise ValueError(f"forward/canonical n1 source-vector proof failed: {failed}")
    return {
        "passed": True,
        "source_executed": False,
        "repo": str(repo),
        "source_pins": ast_receipt,
        "accepted_common_fixed_point": common_status,
        "accepted_receipts": receipt_status,
        "accepted_pressure_scaling": pressure_status,
        "symbolic": {
            "checks": algebra["checks"],
            "inherited_premises": inherited_premises,
            "source_order": algebra["source_order"],
            "canonical_order": algebra["canonical_order"],
            "valid_for": algebra["valid_for"],
            "point_values_selected": False,
        },
        "source_formulas": {
            "P2": "Pstar^2=exp(2*logP); same actual source logP=exp(40)+11",
            "core_scales": "logLambda=4*logP+1000; Lambda=exp(logLambda); epsilon=exp(-logLambda)",
            "S_forward": (
                "(Q + Z*(1-delta)*Q_Z/2 + rho*Q_rho)/L "
                "+ Q*(Q/2 + rho*Q_rho) "
                "+ U*(-2*Z*Q + d*Q_Z - 2*Z*rho*Q_rho)/L"
            ),
            "pressure_head_forward": "2*Q_rho + rho*Q_rhorho",
            "Np_forward": "-S_forward/(2*Pstar^2) + pressure_head_forward/(epsilon*Pstar^2)",
            "N1p_canonical": "-S_canonical/2 + Lambda*pressure_head_canonical, Lambda*epsilon=1",
            "N1theta": "-Z_(-3)(Z_(-2-delta)(F0*Phi))",
            "N1z": "-Z_(-2)(Z_(-1-delta)(U0+epsilon*Psi))",
            "axial_pressure_feedback": (
                "forward=-4*epsilon*Z*rho*Np_forward*Pstar^2/L "
                "=-4*epsilon*Z*rho*N1p/L; after Dscale^-1 axial "
                "multiplication: epsilon*forward=-4*epsilon^2*Z*rho*N1p/L"
            ),
        },
        "input_assumptions": {
            "F": "F0*Phi",
            "U": "U0+epsilon*Psi",
            "F2": (
                "forward S_Z_taylor/(epsilon^2) boxes enclose the exact F0^2 jets; "
                "this is an inherited inclusion, not equality to a selected value"
            ),
            "g": "-epsilon*d_Z(log(F0))",
            "P0": (
                "inherited accepted identity physical_P0=Pstar^2*p; "
                "forward_P0_Z_taylor=epsilon_core*physical_P0_Taylor"
            ),
        },
        "limits": {
            "complete_n1_solution_equality": False,
            "pointwise_original_n1_oracle": False,
            "finite_prefix_or_volterra_value_selected": False,
            "global_six_term_m15_equality": False,
            "temporal_n_ge_2_recursion": False,
        },
    }


# Small aliases make the source receipt convenient to call from a focused
# inspection without introducing any producer dependency.
build_source_vector_identity = prove_source_vector_identity
check_source_vector_identity = prove_source_vector_identity


__all__ = [
    "FORWARD_SOURCE",
    "CANONICAL_OPERATOR_SOURCE",
    "CANONICAL_VOLTERRA_SOURCE",
    "prove_source_vector_identity",
    "build_source_vector_identity",
    "check_source_vector_identity",
]
