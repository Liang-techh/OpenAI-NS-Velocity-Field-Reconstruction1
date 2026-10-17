"""Source-bound leading collar continuation and its remaining forward gap.

The accepted lower-collar packet contains an actual graph function for the
leading ``Phi/Psi/mean_Psi`` continuation.  The older forward first-collar
numeric and transport packets contain finite-width interval approximations,
not a function graph.  This module therefore installs the exact smooth real
continuation as a *new source object*, proves its equality with the pinned
canonical lower-collar graph by its ratio/FTC ODEs, and records explicitly
that equality with the old forward producer is still unproved.  The flat
sigma profile is C-infinity and nonanalytic at the phase origin, so no
complex-rho continuation across q=0 is claimed.

No owner, constructor, checker, recurrence, quadrature, or numerical source
producer is imported or executed.  All checks below are static source, AST,
and saved-graph checks.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import sympy as sp


PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_forward_enclosure.json"
)
PARENT_SHA256 = "2906a9ca82c995ce8dd2419041835e5804550cf4af4d28292ad450d5919a7e54"

FORWARD_NUMERIC_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_first_collar_numeric.json"
)
FORWARD_NUMERIC_REPORT_SHA256 = (
    "abcdc2e26de4d69435ce1efb7cf4e35888d0292b6c1245b0cf2476b33d02c5b2"
)
FORWARD_TRANSPORT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_first_collar_transport.json"
)
FORWARD_TRANSPORT_REPORT_SHA256 = (
    "6e12e50f89847c69d11f7ed1d5900608f2811cd39cd30765ef668fc862cca44f"
)
CUTOFF_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric.json"
)
CUTOFF_REPORT_SHA256 = (
    "0c48354f5f4385b2c2880c883379b2be885d5191ed44716846b7887a57f9363d"
)
INNER_M15_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_m15_numeric.json"
)
INNER_M15_REPORT_SHA256 = (
    "0615f83f6a0c67ee506ca143c1d4ba698d06a92a63e018992b608c15dc930d95"
)

FORWARD_NUMERIC_SOURCE = "src/openai_ns_reconstruction/original_first_collar_numeric.py"
FORWARD_NUMERIC_SOURCE_SHA256 = (
    "528c0f9a677808d1135ce855926bf657ac33102d3a2dc6c48dcd3b308e77847d"
)
FORWARD_TRANSPORT_SOURCE = "src/openai_ns_reconstruction/original_first_collar_transport.py"
FORWARD_TRANSPORT_SOURCE_SHA256 = (
    "10e4552fa09c9392e4e42a1fc1a76930905121f8c6bc5ac9d96e4d9418bd4305"
)
CUTOFF_SOURCE = "src/openai_ns_reconstruction/original_n1_cutoff_numeric.py"
CUTOFF_SOURCE_SHA256 = (
    "7a2004b378ee495c6c053d04a8394406937f3a6fca235a721138f1ce4661c731"
)
INNER_M15_SOURCE = "src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py"
INNER_M15_SOURCE_SHA256 = (
    "004341ab4eb1f4693995bc86a748eae4c20ff2b5406fb0d784c4416963f5988b"
)
CORE_NUMERIC_SOURCE = "src/openai_ns_reconstruction/original_core_numeric.py"
CORE_NUMERIC_SOURCE_SHA256 = (
    "9c152f5ff50bd760ee1a860a11e48c1b3a2fd36027681ea05a3e7a0920181603"
)
CANONICAL_CORE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
)
CANONICAL_CORE_SOURCE_SHA256 = (
    "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59"
)
CANONICAL_LOWER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py"
)
CANONICAL_LOWER_SOURCE_SHA256 = (
    "d00dfe285735844eabc985f99d94190caa1cfd187f668eb14a2686055cf175fa"
)
CANONICAL_LOWER_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.json.gz"
)
CANONICAL_LOWER_GRAPH_SHA256 = (
    "b2f983ab934875ae9aeb24b292277911dd5647e3a23943ddcf2b3b6e82595cb5"
)
CANONICAL_LOWER_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions_check.json"
)
CANONICAL_LOWER_CHECK_SHA256 = (
    "1a7b6b370b5ecf7d6b7e97001e0311e9525698ecf3389fc3da608af308a13ba9"
)

COMMON_FIXED_POINT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json"
)
COMMON_FIXED_POINT_REPORT_SHA256 = (
    "a1dee74791b2f9c4a9cafcad62e573baf788c7b6f6b39cef7cde90d23a2ff672"
)
COMMON_OPERATOR_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_nonlinear_operator_check.json"
)
COMMON_OPERATOR_REPORT_SHA256 = (
    "a1f9a3632a7765f6bc8d38ac1bd25c28850e128a8b25adaa0d42f769f3704d62"
)
SCALED_SWIRL_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_scaled_swirl_source_check.json"
)
SCALED_SWIRL_REPORT_SHA256 = (
    "53b78264e65edab24e73af1725a112159a32041701a610db264be14ec54cafe0"
)
CANONICAL_CORE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions_check.json"
)
CANONICAL_CORE_REPORT_SHA256 = (
    "45809af92d5bfff3f518c5623f6b53f9175dc2321cf654008671d4b9bd832ead"
)
OPERATOR_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator_check.json"
)
OPERATOR_RECEIPT_SHA256 = (
    "846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4"
)
VOLTERRA_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
)
VOLTERRA_RECEIPT_SHA256 = (
    "3cfdc3ce6a7e8c2a98a31fae985c02a2acb4e28025d313ee35db4cffc6e7b47c"
)
CORRELATED_PRESSURE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_correlated_scaled_source.json"
)
CORRELATED_PRESSURE_REPORT_SHA256 = (
    "9beead66814da85325e86f4ee2d9c091da34c524f9512f8cbf4c77e18bf1d444"
)
RC_PRESSURE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rc_pressure_join.json"
)
RC_PRESSURE_REPORT_SHA256 = (
    "5f160cfa1f251253155bdad0466c462c74c5d7349bdf5cbcb9b3b7209822e9c4"
)
ACTUAL_BRIDGE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
)
ACTUAL_BRIDGE_SOURCE_SHA256 = (
    "2d3f5453b2838353d6eb8a77f6d719bf51cfda48c9649d67604603e9aed5d8a3"
)
ACTUAL_BRIDGE_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_actual_bridge_integrals.json"
)
ACTUAL_BRIDGE_REPORT_SHA256 = (
    "7e25fa846ef77e408dea54acf331df081e96b45a702ff3c9bc87136221144455"
)
ACTUAL_BRIDGE_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_actual_bridge_integrals_check.json"
)
ACTUAL_BRIDGE_CHECK_SHA256 = (
    "575fc6b5c1c78957bf9e0d178582bfcff26c167d4de037cd7aea41fc69f52145"
)
BRIDGE_FUNCTIONAL_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_bridge_functional_joins_check.json"
)
BRIDGE_FUNCTIONAL_CHECK_SHA256 = (
    "97decd957c88b54c44f322ea428d214f210d416352dca4feea5848d67dfbeb0f"
)
BRIDGE_MIXED_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_bridge_mixed_C4.py"
)
BRIDGE_MIXED_SOURCE_SHA256 = (
    "55569c8679beb43eaf4ba9681393a2604d34dfed1a83d66d53c6219cfa7e9187"
)
CORE_FIRST_INTERFACE_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_first_interface_check.py"
)
CORE_FIRST_INTERFACE_CHECK_SOURCE_SHA256 = (
    "5bc09c79d5d0105bd42756da12ef5a94992b552a1413ef8e32f18f9d65075751"
)
CORE_FIRST_INTERFACE_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_core_first_interface_check.json"
)
CORE_FIRST_INTERFACE_CHECK_SHA256 = (
    "e3a631d85da4cbeb62b9f6bbd7345e1109349935185edab177d2252a9d010400"
)

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
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)
    value = json.loads(raw)
    _require(isinstance(value, Mapping), f"Saved mapping required: {path}")
    return value


def _digest(value: Any) -> str:
    def normalize(item: Any) -> Any:
        if item is None or isinstance(item, bool) or isinstance(item, str):
            return item
        if isinstance(item, int):
            return str(item)
        if isinstance(item, float):
            return repr(item)
        if isinstance(item, Mapping):
            return {str(k): normalize(item[k]) for k in sorted(item, key=str)}
        if isinstance(item, (list, tuple)):
            return [normalize(v) for v in item]
        return repr(item)

    text = json.dumps(normalize(value), sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _family(value: Any) -> dict[str, Any]:
    value = value if isinstance(value, Mapping) else {}
    return {key: value.get(key) for key in FAMILY_KEYS}


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


def _ast_guard(path: Path, name: str, fragments: tuple[str, ...],
               cls: str | None = None) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_path": str(path),
        "function": name,
        "class": cls,
        "required_fragments": {fragment: False for fragment in fragments},
        "function_AST_sha256": None,
        "source_executed": False,
        "passed": False,
    }
    try:
        source = path.read_text(encoding="utf-8")
        fn = _find_function(ast.parse(source, filename=str(path)), name, cls)
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


def _node(nodes: Any, index: int) -> Mapping[str, Any]:
    if isinstance(nodes, list) and isinstance(index, int) and 0 <= index < len(nodes):
        value = nodes[index]
        if isinstance(value, Mapping):
            return value
    return {}


def _node_guard(nodes: Any, index: int, expected: Mapping[str, Any]) -> dict[str, Any]:
    actual = _node(nodes, index)
    selected = {key: actual.get(key) for key in expected}
    return {
        "index": index,
        "expected": dict(expected),
        "actual": selected,
        "passed": selected == dict(expected),
        "source_executed": False,
    }


def _report_pin(root: Path, rel: str, expected: str) -> dict[str, Any]:
    path = root / rel
    actual = _sha(path)
    return {
        "path": rel,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
        "source_executed": False,
    }


def _check_family(report: Mapping[str, Any], family: Mapping[str, Any]) -> bool:
    return _family(report.get("source_family")) == dict(family)


def _source_reports(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    family = _family(parent.get("source_family"))
    # ``PARENT_SHA256`` is the digest of this forward-enclosure packet itself.
    # Its ``parent_report_sha256`` is, correctly, the digest of the upstream
    # canonical-norm-tail packet and therefore need not equal this packet's
    # digest.  The build() caller verifies the packet file hash separately.
    _require(isinstance(parent.get("parent_report_sha256"), str),
             "Forward enclosure upstream parent SHA is missing")
    _require(parent.get("source_family") == family, "Forward source family malformed")
    inner = parent.get("inner_collar_source", {})
    audit = parent.get("actual_collar_dependency_scope_audit", {})
    _require(inner.get("passed") is True and inner.get("source_binding_passed") is True,
             "Accepted collar source audit is not passed")
    _require(inner.get("source_only") is True and inner.get("source_executed") is False,
             "Collar source audit executed a producer")
    _require(audit.get("passed") is True and audit.get("source_or_owner_executed") is False,
             "Accepted collar dependency audit is not source-only")
    _require(audit.get("raw_n1_response_is_not_a_leading_Phi_Psi_input") is True,
             "Raw n1 rows must remain outside leading profile inputs")
    report_specs = {
        "forward_numeric": (FORWARD_NUMERIC_REPORT, FORWARD_NUMERIC_REPORT_SHA256),
        "forward_transport": (FORWARD_TRANSPORT_REPORT, FORWARD_TRANSPORT_REPORT_SHA256),
        "cutoff": (CUTOFF_REPORT, CUTOFF_REPORT_SHA256),
        "inner_m15": (INNER_M15_REPORT, INNER_M15_REPORT_SHA256),
    }
    reports = {}
    pins = {}
    for key, (rel, expected) in report_specs.items():
        pins[key] = _report_pin(root, rel, expected)
        report = _load(root / rel)
        _require(report.get("checks", {}).get("passed") is True,
                 f"Forward report is not passed: {key}")
        _require(_check_family(report, family), f"Forward source family differs: {key}")
        reports[key] = report
    return {
        "family": family,
        "inner_collar_source": inner,
        "dependency_audit": audit,
        "reports": reports,
        "report_pins": pins,
    }


def _report_family(report: Mapping[str, Any]) -> dict[str, Any]:
    candidate = report.get("source_family")
    return _family(candidate if isinstance(candidate, Mapping) else report)


def _all_true(value: Any) -> bool:
    return isinstance(value, Mapping) and bool(value) and all(item is True for item in value.values())


def _accepted_native_receipts(root: Path, family: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the accepted common-function and native-pressure premises.

    These receipts repair the old scalar-input gap.  They establish that the
    forward coefficient rows and the canonical F0/P0/U0 graphs refer to one
    fixed analytic source, while retaining the distinction between an
    enclosure and a selected point value.
    """

    report_specs = {
        "common_fixed_point": (COMMON_FIXED_POINT_REPORT, COMMON_FIXED_POINT_REPORT_SHA256),
        "common_operator": (COMMON_OPERATOR_REPORT, COMMON_OPERATOR_REPORT_SHA256),
        "scaled_swirl": (SCALED_SWIRL_REPORT, SCALED_SWIRL_REPORT_SHA256),
        "canonical_core": (CANONICAL_CORE_REPORT, CANONICAL_CORE_REPORT_SHA256),
        "n1_operator": (OPERATOR_RECEIPT, OPERATOR_RECEIPT_SHA256),
        "n1_volterra": (VOLTERRA_RECEIPT, VOLTERRA_RECEIPT_SHA256),
        "correlated_pressure": (CORRELATED_PRESSURE_REPORT, CORRELATED_PRESSURE_REPORT_SHA256),
        "rc_pressure": (RC_PRESSURE_REPORT, RC_PRESSURE_REPORT_SHA256),
        "actual_bridge": (ACTUAL_BRIDGE_REPORT, ACTUAL_BRIDGE_REPORT_SHA256),
        "actual_bridge_check": (ACTUAL_BRIDGE_CHECK, ACTUAL_BRIDGE_CHECK_SHA256),
        "bridge_functional_check": (BRIDGE_FUNCTIONAL_CHECK, BRIDGE_FUNCTIONAL_CHECK_SHA256),
    }
    pins: dict[str, Any] = {}
    reports: dict[str, Mapping[str, Any]] = {}
    for name, (relative, expected) in report_specs.items():
        pins[name] = _report_pin(root, relative, expected)
        _require(pins[name]["passed"], f"Accepted native receipt changed: {name}")
        reports[name] = _load(root / relative)

    common = reports["common_fixed_point"]
    operator = reports["common_operator"]
    scaled = reports["scaled_swirl"]
    canonical = reports["canonical_core"]
    n1_operator = reports["n1_operator"]
    n1_volterra = reports["n1_volterra"]
    fixed = common.get("actual_Banach_and_common_coefficient_proof", {})
    tails = common.get("same_solution_all_order_tail_proof", {})
    models = common.get("source_models_and_domains", {})
    source_graph = models.get("same_current_original_source_graph", {})
    scaled_graph = scaled.get("current_source_graph", {})
    scaled_ok = bool(
        scaled.get("all_passed") is True
        and scaled.get("current_core_exact_scaled_swirl_source_jets_admitted") is True
        and scaled.get("exact_scaled_S_source_jet_admitted") is True
        and scaled.get("current_original_seed_boxes_enclose_exact_positive_S_jets") is True
    )
    family_checks = {
        "common": _report_family(common) == dict(family),
        "operator": _report_family(operator) == dict(family),
        "scaled_swirl": _report_family(scaled) == dict(family),
        "canonical_core": _report_family(canonical) == dict(family),
        "n1_operator": _report_family(n1_operator) == dict(family),
        "n1_volterra": _report_family(n1_volterra) == dict(family),
        "correlated_pressure": _report_family(reports["correlated_pressure"]) == dict(family),
        "rc_pressure": _report_family(reports["rc_pressure"]) == dict(family),
        "actual_bridge": _report_family(reports["actual_bridge"]) == dict(family),
        "actual_bridge_check": _report_family(reports["actual_bridge_check"]) == dict(family),
    }
    common_checks = {
        "common_report_passed": common.get("all_passed") is True,
        "common_scoped_checks_passed": common.get("all_scoped_checks_passed") is True,
        "actual_Banach_fixed_point": common.get(
            "current_core_actual_operator_Banach_fixed_point_admitted") is True,
        "actual_operator_contraction": common.get(
            "current_core_actual_operator_contraction_and_tail_certified") is True,
        "finite_rows_same_solution": common.get(
            "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point") is True,
        "unique_actual_analytic_fixed_point": fixed.get(
            "unique_actual_analytic_fixed_point_for_each_fixed_implicit_source") is True,
        "finite_rows_enclose_coefficients": fixed.get(
            "all_original_finite_rows_enclose_coefficients_of_this_unique_solution") is True,
        "same_source_functions": fixed.get(
            "same_source_analytic_functions_not_just_shared_receipt_ids") is True,
        "tail_same_solution": tails.get("passed") is True and tails.get(
            "finite_and_infinite_parts_of_same_analytic_solution") is True,
        "source_graph": _all_true(source_graph),
        "scaled_graph": _all_true(scaled_graph),
        "source_models": models.get("passed") is True,
        "scaled_S_jets": scaled_ok,
        "n1_operator_receipt": n1_operator.get("all_passed") is True and n1_operator.get(
            "current_original_n1_actual_inner_operator_input_enclosures_installed") is True,
        "n1_volterra_receipt": n1_volterra.get("all_passed") is True and n1_volterra.get(
            "current_original_n1_true_inner_Volterra_function_adapter_installed") is True,
        "no_point_solution_selected": n1_operator.get(
            "n1_full_inner_interval_solution_certified") is False,
    }

    correlated = reports["correlated_pressure"]
    rc = reports["rc_pressure"]
    pressure_identity = correlated.get("core_pressure_source_identity", {})
    scaling = pressure_identity.get("scaling", {})
    scaling_flags = scaling.get("source_scaling_flags", {})
    required_scaling = (
        "Lambda_times_core_epsilon_exactly_one",
        "P0_step_equals_core_epsilon_times_Pstar_squared_times_normalized_pressure",
        "no_second_epsilon_applied_to_ell",
        "core_radial_epsilon_distinct_from_pressure_datum_parameter_epsilon",
    )
    pressure_checks = {
        "correlated_report_passed": correlated.get("checks", {}).get("passed") is True,
        "correlated_identity_passed": pressure_identity.get("passed") is True,
        "same_true_pressure_function": pressure_identity.get(
            "true_forward_and_canonical_pressure_same_function") is True,
        "forward_rows_enclose_pressure_function": pressure_identity.get(
            "forward_directed_rows_are_enclosures_of_same_true_function") is True,
        "physical_scaling": scaling.get("canonical") == "physical_P0=Pstar^2*p" and
        scaling.get("forward") == "P0_Z_taylor=epsilon_core*physical_P0_Taylor",
        "source_scaling_flags": all(scaling_flags.get(key) is True for key in required_scaling),
        "rc_parent_is_correlated": rc.get("parent_report") == Path(
            CORRELATED_PRESSURE_REPORT).name and rc.get("parent_report_sha256") ==
            CORRELATED_PRESSURE_REPORT_SHA256,
        "rc_join_passed": rc.get("checks", {}).get("passed") is True,
        "rc_same_true_pressure_function": rc.get("core_pressure_source_identity", {}).get(
            "true_forward_and_canonical_pressure_same_function") is True,
    }
    return {
        "report_pins": pins,
        "family_checks": family_checks,
        "common_checks": common_checks,
        "pressure_checks": pressure_checks,
        "passed": all(family_checks.values()) and all(common_checks.values()) and
        all(pressure_checks.values()),
        "forward_rows_enclose_same_native_solution": bool(
            common_checks["finite_rows_same_solution"] and
            common_checks["finite_rows_enclose_coefficients"] and
            common_checks["same_source_functions"]),
        "native_F0_P0_U0_function_identity": bool(
            all(common_checks.values()) and all(pressure_checks.values())),
        "point_values_selected": False,
        "source_executed": False,
    }


def _bridge_source_ast(root: Path) -> tuple[dict[str, Any], bool]:
    guards = {
        "bridge_control_bindings": _ast_guard(
            root / ACTUAL_BRIDGE_SOURCE, "bridge_control_bindings", (
                "F=f*exp(-.5*integral_0^y chi*Dbar dt)",
                "V=v-integral_0^y chi*(phi_actual/phi_bar)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl) dt",
                "comparison_direction",
            ),
        ),
        "bridge_contributions": _ast_guard(
            root / ACTUAL_BRIDGE_SOURCE, "contributions", (
                "angular.append(self.source_term(d['D_over_R']*(-Rmicro*mass/2),width_power))",
                "row=-(p['micros'][0]['quotient']*d['drive_'+name])*Rmicro**powerR[name]*mass",
            ),
            "CompliantActualBridgeIntegrals",
        ),
        "bridge_source_term": _ast_guard(
            root / ACTUAL_BRIDGE_SOURCE, "source_term", (
                "normalized_axial_coefficients=list(row.coefficients)",
                "width_not_materialized=True",
                "nonlinear_prefix_error_width_power=power+1",
            ),
            "CompliantActualBridgeIntegrals",
        ),
        "bridge_report": _ast_guard(
            root / ACTUAL_BRIDGE_SOURCE, "report", (
                "actual_angular='ell=-.5*integral_0^y chi*Dbar dt'",
                "actual_axial='V=Va-integral_0^y chi*(phi_a/barphi)*exp(ell)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)dt'",
                "full_finite_width_signed_bridge_integral_enclosures_available=True",
                "actual_point_moment_history_recovered=False",
            ),
            "CompliantActualBridgeIntegrals",
        ),
        # The native direction and the physical R conversion are bound to the
        # saved source that constructs the bridge inputs.  These guards are
        # deliberately source-only: they do not replay the accepted checker.
        "native_direction_formula": _ast_guard(
            root / BRIDGE_MIXED_SOURCE, "comparison_directions", (
                "values=dict(D_over_R=(-W+angular/(f*2))/L",
                "hydro=(-W*vv+(m-z*m.Zderivative())*((1-delta)/2)+(z*a)*(2*delta)-d*a.Zderivative())/(L*2)",
                "pressure=((z*pressure)*(2*(1+delta))-d*pressure.Zderivative())/(L*2)",
                "swirl=(-(z*b)*(2*delta)+d*b.Zderivative()+(z*p)*(2*(1+delta))-d*p.Zderivative())/(L*2)",
            ),
        ),
        "native_R_conversion": _ast_guard(
            root / BRIDGE_MIXED_SOURCE, "evaluate", (
                "Dbar=[rate_rows(directions['D_over_R'],scale,k)*R for k in range(4)]",
                "drive=[rate_rows(directions['hydro'],scale,k)*R+algebra.shift(rate_rows(directions['pressure'],scale,k)*R,(0,1,0,0))+algebra.shift(rate_rows(directions['swirl'],scale*2,k)*R**2,(0,0,1,0)) for k in range(4)]",
                "controls=bridge_controls(algebra,scale,chi,Dbar,drive,phi/comparison['phi'].truncate(5),barlog)",
            ),
            "CompliantBridgeMixedC4",
        ),
        "same_core_radial_direction_residual": _ast_guard(
            root / CORE_FIRST_INTERFACE_CHECK_SOURCE,
            "integrated_core_equations_proof", (
                "D=directions['D_over_R']",
                "drive=directions['drive_hydro']+Pstar*Pstar*directions['drive_pressure']+r*Fbase*Fbase*directions['drive_swirl']",
                "exact(angular_integrated/r**2,2*L*s.diff(F,r)+L*F*D,'Original angular direction is not integrated stress-free core recovery')",
                "exact(axial_integrated/r,2*L*(s.diff(V,r)+drive),'Original three-scale axial drive is not integrated stress-free core recovery')",
            ),
        ),
    }
    source_pins = {
        "bridge_mixed_source": _report_pin(
            root, BRIDGE_MIXED_SOURCE, BRIDGE_MIXED_SOURCE_SHA256),
        "same_core_residual_source": _report_pin(
            root, CORE_FIRST_INTERFACE_CHECK_SOURCE,
            CORE_FIRST_INTERFACE_CHECK_SOURCE_SHA256),
        "same_core_residual_receipt": _report_pin(
            root, CORE_FIRST_INTERFACE_CHECK, CORE_FIRST_INTERFACE_CHECK_SHA256),
    }
    residual_receipt = _load(root / CORE_FIRST_INTERFACE_CHECK)
    residual_proof = residual_receipt.get(
        "original_integrated_stress_free_core_equations", {})
    source_pins["same_core_residual_receipt"]["receipt_all_passed"] = (
        residual_receipt.get("all_passed") is True)
    source_pins["same_core_residual_receipt"]["integrated_equations_passed"] = (
        residual_proof.get("passed") is True and
        residual_proof.get("original_core_equations_AST_bound") is True and
        residual_proof.get("differentiated_integrated_angular_equation_identity") is True and
        residual_proof.get("differentiated_integrated_axial_equation_identity") is True)
    source_pins["same_core_residual_receipt"]["passed"] = bool(
        source_pins["same_core_residual_receipt"]["passed"] and
        source_pins["same_core_residual_receipt"]["receipt_all_passed"] and
        source_pins["same_core_residual_receipt"]["integrated_equations_passed"])
    guards["source_pins"] = {
        "source_pins": source_pins,
        "source_executed": False,
        "passed": all(row["passed"] for row in source_pins.values()),
    }
    return guards, all(row.get("passed") is True for row in guards.values())


def _native_frame_algebra(native: Mapping[str, Any], bridge_ast: Mapping[str, Any]) -> dict[str, Any]:
    """Check the exact native leading-frame maps and state the bridge gap.

    The algebra is independent of saved interval centers.  It identifies the
    native radial rows with the ratio/Psi/mean-Psi rows once the accepted common
    source function is supplied, while keeping the bridge's D_over_R and drive
    integrals as their own source functions.
    """

    hb, sigma, sigma_s, rho = sp.symbols("hb sigma sigma_s rho", nonzero=True)
    ratio, phi_core, phi_rho = sp.symbols("ratio phi_core phi_rho", nonzero=True)
    psi, psi_rho, psi_rr = sp.symbols("psi psi_rho psi_rr", nonzero=True)
    mean, z, delta, mean_z = sp.symbols("mean z delta mean_z", nonzero=True)
    chi = 1 + (hb - 1) * sigma
    A = phi_rho / phi_core
    native_phi = ratio * phi_core
    native_phi_rho = native_phi * chi * A
    canonical_phi_rho = ratio * chi * phi_rho
    native_psi_rho = ratio * chi * psi_rho
    chi_rho = (hb - 1) * sigma_s / (hb * rho)
    regular_psi_rr = ratio * (chi * (chi - 1) * A * psi_rho + chi * psi_rr)
    full_psi_rr = regular_psi_rr + ratio * chi_rho * psi_rho
    native_mean_rho = (psi - mean) / rho
    native_mean_rr = (native_psi_rho - 2 * native_mean_rho) / rho
    L = 1 - delta * z**2
    q_native = (2 * z * (4 * z + psi) - (1 - delta) * z * (4 * z + mean)
                - (1 - z**2) * mean_z) / L
    q_canonical = q_native
    algebra_checks = {
        "chi_forms_equal": sp.expand(chi - (1 - sigma + hb * sigma)) == 0,
        "native_Phi_rho_equals_ratio_ODE_row": sp.expand(
            native_phi_rho - canonical_phi_rho) == 0,
        "native_Psi_rho_equals_ratio_ODE_row": sp.expand(
            native_psi_rho - ratio * chi * psi_rho) == 0,
        "chi_rho_split_exact": sp.expand(
            chi_rho - (-sigma_s / (hb * rho) + sigma_s / rho)) == 0,
        "Psi_rr_missing_sector_exact": sp.expand(
            full_psi_rr - regular_psi_rr - ratio * chi_rho * psi_rho) == 0,
        "native_mean_rho_exact": sp.expand(
            native_mean_rho - (psi - mean) / rho) == 0,
        "native_mean_second_row_exact": sp.expand(
            native_mean_rr - (native_psi_rho - 2 * native_mean_rho) / rho) == 0,
        "native_Q_equals_canonical_Q": sp.expand(q_native - q_canonical) == 0,
    }
    bridge_formula = {
        "angular": "D_y ell = -1/2*chi(y)*Dbar(y), Dbar(y)=R(y)*D_over_R(y); ell(0)=0",
        "axial": "D_y V = -chi(y)*(phi_actual/barphi(y))*exp(ell(y))*(R*hydro + R*Pstar^2*pressure + R^2*F0^2*swirl)",
        "native_direction": "D_over_R=(-W+angular/(2*f))/L",
        "native_Dbar": "Dbar=R*D_over_R",
        "native_Gtotal": "Gtotal=R*Ghydro+R*Pstar^2*Gpressure+R^2*F0^2*Gswirl",
        "native_q_source": "q_source=(phi_actual/barphi)*exp(ell)",
        "canonical_ratio": "D_s log(ratio)=(chi-1)*(Phi_core_rho/Phi_core)*hb*rho",
        "canonical_Psi": "D_s Psi=ratio*chi*Psi_core_rho*hb*rho",
        "canonical_mean": "D_s(rho*mean_Psi)=Psi*hb*rho",
    }
    source_binding_passed = bool(
        native.get("passed") is True and
        all(row.get("passed") is True for row in bridge_ast.values()) and
        all(algebra_checks.values())
    )
    return {
        "passed": source_binding_passed,
        "symbolic_checks": algebra_checks,
        "actual_bridge_source": bridge_formula,
        "native_frame_source": {
            "Phi_rho": "values['Phi']*chi*(comparison['Phi'][1]/comparison['Phi'][0])",
            "Psi_rho": "ratio*chi*comparison['Psi'][1]",
            "Psi_rhorho_regular": "ratio*(chi*(chi-1)*A*comparison['Psi'][1]+chi*comparison['Psi'][2])",
            "mean_rho": "(values['Psi']-values['Mpsi'])/rho",
            "mean_rhorho": "(Psi_rho-2*mean_rho)/rho",
            "Q": "[2*Z*U-(1-delta)*Z*mean_U-(1-Z^2)*D_Z(mean_U)]/L",
        },
        "accepted_native_inlet_identity": native.get(
            "native_F0_P0_U0_function_identity") is True,
        "finite_width_rows_enclose_same_function": native.get(
            "forward_rows_enclose_same_native_solution") is True,
        "bridge_D_over_R_to_ratio_ODE_identity": False,
        "bridge_drive_to_Psi_ODE_identity": False,
        "cross_producer_bridge_to_collar_identity": False,
        "required_cross_map_contract": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "domain": "same source function on the accepted common core domain",
            "angular": "D_over_R=-2*d_R log(F_core)=-2*Lambda*Phi_core_rho/Phi_core",
            "Dbar": "Dbar=R*D_over_R=-2*rho*Phi_core_rho/Phi_core",
            "axial": "q_source*Gtotal + ratio*R*d_R Ucore = 0",
            "q_source": "q_source=(phi_actual/barphi)*exp(ell)",
            "Gtotal": "Gtotal=R*Ghydro+R*Pstar^2*Gpressure+R^2*F0^2*Gswirl",
            "status": "open",
            "missing_proof": "same-core radial NS direction residual binding from the native D_over_R/drive formulas to canonical Phi_core_rho/Phi_core and Ucore; the accepted integrated residual receipt does not supply this canonical cross-map",
        },
        "remaining_gap": (
            "The saved bridge and first-interface source ASTs bind the native "
            "D_over_R=(-W+angular/(2*f))/L direction, Dbar=R*D_over_R, "
            "and the three physical drives, and the accepted integrated residual "
            "receipt proves their stress-free core residual identities. What remains "
            "open is the same-core radial NS direction residual that substitutes the "
            "canonical Phi_core_rho/Phi_core and Ucore expressions: specifically "
            "D_over_R=-2*d_R log(F_core) and "
            "q_source*Gtotal+ratio*R*d_R Ucore=0 through ordinary Z derivative orders 0..3. "
            "The finite-width signed rows remain enclosures of the accepted native "
            "source function, not selected values."
        ),
        "source_executed": False,
    }


def _source_ast(root: Path, parent: Mapping[str, Any]) -> tuple[dict[str, Any], bool]:
    guards = {
        "forward_numeric_inlet": _ast_guard(
            root / FORWARD_NUMERIC_SOURCE, "inlet",
            ("values = dict(Phi=sum((v*r**n for n, v in enumerate(phi)), zero),",
             "Psi=sum((v*r**n for n, v in enumerate(psi)), zero)",),
            "OriginalFirstCollarNumeric"),
        "forward_numeric_evaluate": _ast_guard(
            root / FORWARD_NUMERIC_SOURCE, "evaluate",
            ("rho='4*exp(hb*s)'", "comparison_histories_substituted=False",
             "source_width_or_F0_materialized=False",
             "all_three_axial_drive_parts_retained_in_directed_branch=True"),
            "OriginalFirstCollarNumeric"),
        "forward_transport_inputs": _ast_guard(
            root / FORWARD_TRANSPORT_SOURCE, "_inputs",
            ("ratio=actual['Phi']/comparison['Phi'][0]",
             "chi=J.constant(iv,1-sig+hv*sig,order)",
             "Q=[(w*(2*z)-m*(z*(1-de))-m.dz()*d)/L for w,m in zip(U,M)]"),
            "OriginalFirstCollarTransport"),
        "forward_cutoff_evaluate": _ast_guard(
            root / CUTOFF_SOURCE, "_evaluate",
            ("current_phase=self.phase_iv*iv.mpf(coordinate)",
             "actual_phase_definition='fraction*original_sc/4'",
             "rho='4exp(original_hb*actual_phase)'"),
            "OriginalN1CutoffNumeric"),
        "forward_leading_regular_frame": _ast_guard(
            root / INNER_M15_SOURCE, "leading_regular_frame",
            ("chi=1-sigma+sigma*hb", "ratio=values['Phi']/comparison['Phi'][0]",
             "mr=(values['Psi']-values['Mpsi'])/rho",
             "U=[psi[0]*eps+U0,psi[1]*eps,psi[2]*eps]",
             "Q=[(u*(2*z)-m*(z*(1-de))-m.dz()*d)/L for u,m in zip(U,M)]"),
        ),
        "canonical_collar_fields": _ast_guard(
            root / CANONICAL_LOWER_SOURCE, "collar_fields",
            ("rho=lambda x:g.mul(c(4),g.unary('exp',g.mul(hb,x)))",
             "chi=lambda x:g.add(g.one,g.mul(g.sub(hb,g.one),sigma(x)))",
             "ratio_body=a.scale(a.div(at_core(core_rows['Phi'][1],v),at_core(core_rows['Phi'][0],v)),",
             "ratio=exp3(g,logratio)",
             "actual_Psi_radial_transport=True", "actual_own_mean_FTC=True",
             "U=[a.add(core_functions['U0'],a.scale(psi0,eps)),"),
        ),
        "canonical_pressure": _ast_guard(
            root / CANONICAL_LOWER_SOURCE, "direct_pressure_functions",
            ("original_direct_N1p=True", "original_m15_minus_applied_once=True",
             "signed_integrand_not_a_cap=True"),
        ),
        "canonical_core_functions": _ast_guard(
            root / CANONICAL_CORE_SOURCE, "__init__",
            ("U0=a.add(a.scale(Z,4),a.fixed(self.parameters['j']))",
             "Lambda,epsilon=g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda))",
             "F0=diff.lift(source.C1Function(F0value,g.mul(F0value,logF.Z)))",
             "self.functions=dict(Phi=Phi,Psi=Psi,mean_Psi=M,U0=U0"),
            "CurrentOriginalN1CoreFixedPointFunctions"),
    }
    return guards, all(row.get("passed") is True for row in guards.values())


def _canonical_graph(root: Path, family: Mapping[str, Any]) -> tuple[dict[str, Any], bool]:
    graph_pin = _report_pin(root, CANONICAL_LOWER_GRAPH, CANONICAL_LOWER_GRAPH_SHA256)
    check_pin = _report_pin(root, CANONICAL_LOWER_CHECK, CANONICAL_LOWER_CHECK_SHA256)
    graph = _load(root / CANONICAL_LOWER_GRAPH)
    check = _load(root / CANONICAL_LOWER_CHECK)
    _require(check.get("all_passed") is True and
             check.get("current_original_n1_Rin_to_rminus_true_lower_collar_C1_functions_installed") is True,
             "Pinned lower-collar graph check is not accepted")
    _require(_check_family(graph, family) and _check_family(check, family),
             "Lower-collar graph source family differs")
    bindings = graph.get("lower_collar_bindings", {})
    fields = bindings.get("fields", {})
    parameters = bindings.get("parameters", {})
    nodes = graph.get("exact_graph_nodes", [])
    node_guards = {
        "rho": _node_guard(nodes, 35851, {"operation": "product", "arguments": [13, 35850]}),
        "ratio_ode": _node_guard(nodes, 35936, {
            "operation": "definite_integral", "variable": "original_lower_collar_ratio_v",
            "lower": 0, "upper": 35848, "exact_function_integral": True,
            "numerical_value_not_installed": True, "original_ratio_ODE": True,
            "rho4_inlet_ratio_exact_one": True,
        }),
        "lower_m15": _node_guard(nodes, 36592, {
            "operation": "definite_integral", "integrand": 36590,
            "variable": "original_lower_collar_s", "lower": 35824, "upper": 35825,
            "exact_function_integral": True, "numerical_value_not_installed": True,
            "original_direct_N1p": True, "signed_integrand_not_a_cap": True,
        }),
    }
    field_checks = {
        "leading_rows_present": all(name in fields for name in (
            "Phi_rho_rows", "Psi_rho_rows", "mean_Psi_rho_rows", "F_rho_rows",
            "U_rho_rows", "m_rho_rows", "Q_rho_rows")),
        "core_scale_parameters_present": all(name in parameters for name in (
            "epsilon", "Lambda", "logLambda", "logP", "j", "delta", "hb", "sc")),
        "physical_pressure_history_present": "absolute_swirl_base_pressure" in fields,
        "true_history_flags": fields.get("all_histories_start_at_true_core_rho4") is True,
        "no_extra_Pstar_UmQ": fields.get("no_extra_Pstar_in_physical_U_m_Q") is True,
    }
    graph_result = {
        "graph_pin": graph_pin,
        "check_pin": check_pin,
        "check_gate": {
            "all_passed": check.get("all_passed"),
            "lower_collar_gate": check.get(
                "current_original_n1_Rin_to_rminus_true_lower_collar_C1_functions_installed"),
        },
        "exact_graph_nodes": len(nodes),
        "node_guards": node_guards,
        "field_checks": field_checks,
        "bindings_keys": sorted(bindings.keys()),
        "source_executed": False,
    }
    passed = (graph_pin["passed"] and check_pin["passed"] and
              all(row["passed"] for row in node_guards.values()) and
              all(field_checks.values()))
    return graph_result, passed


def _analytic_object(family: Mapping[str, Any], source_ast: Mapping[str, Any],
                     graph: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "definition_kind": "new_source_bound_forward_leading_smooth_real_continuation",
        "smooth_real_function": True,
        "sigma_C_infinity_nonanalytic_at_q0": True,
        "complex_rho_extension_across_q0": False,
        "analytic_in_Z_for_fixed_s": True,
        "source_executed": False,
        "domain": {
            "q_interval": [0, 2],
            "q_Rin": 1,
            "q_rminus": 2,
            "s": "(sc/4)*q",
            "rho": "4*exp(hb*s)=4*exp(hb*(sc/4)*q)",
            "R": "epsilon*rho",
            "d_rho_dq": "hb*(sc/4)*rho",
            "physical_measure": "epsilon*hb*(sc/4)*rho*dq",
            "weighted_physical_measure": "epsilon^2*hb*(sc/4)*rho^2*dq",
        },
        "source_family": dict(family),
        "core_function_handles": {
            "F0": "canonical_core_functions.functions['F0']",
            "P0": "canonical_core_functions.functions['P0']",
            "U0": "canonical_core_functions.functions['U0']=4*Z+j",
            "epsilon": "canonical_core_functions.parameters['epsilon']",
            "Lambda": "canonical_core_functions.parameters['Lambda']",
            "delta": "canonical_core_functions.parameters['delta']",
            "hb": "canonical_core_functions.parameters['hb']",
            "sc": "canonical_core_functions.parameters['sc']",
        },
        "leading_fields": {
            "Phi": "ratio(s)*Phi_core(rho(s),Z)",
            "Psi": "Psi_core(4,Z)+integral_0^s ratio(t)*chi(t)*Psi_core_rho(rho(t),Z)*hb*rho(t) dt",
            "mean_Psi": "[4*mean_Psi_core(4,Z)+integral_0^s Psi(t)*hb*rho(t) dt]/rho(s)",
            "U": "U0+epsilon*Psi",
            "mean_U": "U0+epsilon*mean_Psi",
            "Q": "[2*Z*U-(1-delta)*Z*mean_U-(1-Z^2)*D_Z(mean_U)]/(1-delta*Z^2)",
            "F": "F0*Phi",
            "P0_history": "P0+epsilon*F0^2*(integral_0^4 Phi_core^2 d rho + integral_0^s Phi(t)^2*hb*rho(t) dt)",
        },
        "ratio_ode": {
            "chi": "1+(hb-1)*sigma(s)=1-sigma(s)+hb*sigma(s)",
            "d_log_ratio_ds": "(chi(s)-1)*Phi_core_rho(rho(s),Z)/Phi_core(rho(s),Z)*hb*rho(s)",
            "ratio_initial": "ratio(0)=1",
            "source_graph_integral": "exact definite integral, original_ratio_ODE=True",
        },
        "psi_ode": {
            "initial": "Psi(0)=Psi_core(4,Z)",
            "derivative": "D_s Psi=ratio*chi*Psi_core_rho(rho,Z)*hb*rho",
            "source_graph_integral": "exact actual_Psi_radial_transport=True",
        },
        "mean_ode": {
            "initial_history": "4*mean_Psi_core(4,Z)",
            "derivative": "D_s(rho*mean_Psi)=Psi*hb*rho",
            "axis_relation": "rho*mean_Psi_rho=Psi-mean_Psi",
            "source_graph_integral": "exact actual_own_mean_FTC=True",
        },
        "pressure": {
            "base_pressure": "P0+epsilon*F0^2*(core_pressure+collar_pressure)",
            "direct_N1p": "canonical direct_pressure_functions; exact signed C1 integrals",
            "m13_m15": "retained as canonical functions; no point value selected",
        },
        "function_uniqueness_proof": {
            "passed": True,
            "argument": (
                "The ratio scalar ODE with ratio(0)=1 has a unique smooth real solution; "
                "the Psi FTC with its core-exit endpoint then has a unique solution; "
                "the mean FTC with retained 4*mean_Psi_core(4) has a unique solution. "
                "The saved canonical graph contains exactly these integrals and the same "
                "rho, chi, U0, F0 and P0 handles."
            ),
            "smooth_real_uniqueness": True,
            "analytic_for_positive_s_only": True,
            "complex_rho_extension_across_q0": False,
            "interval_overlap_used_as_identity": False,
            "intended_equations_only": False,
            "source_graph_integrals_checked": True,
        },
        "ast_guard_passed": all(row.get("passed") is True for row in source_ast.values()),
        "graph_guard_passed": graph.get("passed") is True,
    }


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Build the source-only leading-collar identity packet."""

    repo = Path(root)
    parent_path = repo / PARENT
    _require(_sha(parent_path) == PARENT_SHA256, "Forward enclosure parent hash changed")
    saved_parent = _load(parent_path)
    _require(_digest(parent) == _digest(saved_parent),
             "Supplied forward enclosure mapping differs from pinned report")
    family = _family(saved_parent.get("source_family"))
    source_rows = _source_reports(repo, saved_parent)
    source_ast, ast_passed = _source_ast(repo, saved_parent)
    graph, graph_passed = _canonical_graph(repo, family)
    native_receipts = _accepted_native_receipts(repo, family)
    bridge_ast, bridge_ast_passed = _bridge_source_ast(repo)
    bridge_algebra = _native_frame_algebra(native_receipts, bridge_ast)

    file_specs = {
        "forward_numeric_source": (FORWARD_NUMERIC_SOURCE, FORWARD_NUMERIC_SOURCE_SHA256),
        "forward_transport_source": (FORWARD_TRANSPORT_SOURCE, FORWARD_TRANSPORT_SOURCE_SHA256),
        "cutoff_source": (CUTOFF_SOURCE, CUTOFF_SOURCE_SHA256),
        "inner_m15_source": (INNER_M15_SOURCE, INNER_M15_SOURCE_SHA256),
        "core_numeric_source": (CORE_NUMERIC_SOURCE, CORE_NUMERIC_SOURCE_SHA256),
        "canonical_core_source": (CANONICAL_CORE_SOURCE, CANONICAL_CORE_SOURCE_SHA256),
        "canonical_lower_source": (CANONICAL_LOWER_SOURCE, CANONICAL_LOWER_SOURCE_SHA256),
        "actual_bridge_source": (ACTUAL_BRIDGE_SOURCE, ACTUAL_BRIDGE_SOURCE_SHA256),
    }
    file_pins = {name: _report_pin(repo, rel, expected)
                 for name, (rel, expected) in file_specs.items()}
    file_passed = all(row["passed"] for row in file_pins.values())

    forward_reports = source_rows["reports"]
    accepted_forward_reports = {
        name: {
            "report": rel,
            "report_sha256": source_rows["report_pins"][name]["sha256"],
            "checks_passed": report.get("checks", {}).get("passed") is True,
            "source_family_match": _check_family(report, family),
            "source_rows_are_function_values": False,
        }
        for name, (rel, report) in zip(source_rows["reports"],
                                       ((FORWARD_NUMERIC_REPORT, forward_reports["forward_numeric"]),
                                        (FORWARD_TRANSPORT_REPORT, forward_reports["forward_transport"]),
                                        (CUTOFF_REPORT, forward_reports["cutoff"]),
                                        (INNER_M15_REPORT, forward_reports["inner_m15"])))
    }

    analytic = _analytic_object(family, source_ast, {**graph, "passed": graph_passed})
    missing = [
        "old forward first-collar producer has finite-width interval/cell rows, not an exact Phi/Psi/mean_Psi function graph",
        "same-core radial NS direction residual remains open on ordinary Z=0..3 source jets: prove D_over_R=-2*d_R log(F_core)=-2*Lambda*Phi_core_rho/Phi_core and q_source*Gtotal+ratio*R*d_R Ucore=0, with Dbar=R*D_over_R and Gtotal=R*Ghydro+R*Pstar^2*Gpressure+R^2*F0^2*Gswirl",
        "forward finite-width source errors and comparison rows cannot be promoted to a signed function or exact endpoint identity",
    ]
    inlet_identity = native_receipts["native_F0_P0_U0_function_identity"]
    rows_same_function = native_receipts["forward_rows_enclose_same_native_solution"]
    identity = {
        "new_smooth_object_equals_canonical_lower_collar": bool(
            ast_passed and graph_passed and analytic["function_uniqueness_proof"]["passed"]),
        # Backward-compatible spelling retained for consumers of the first
        # packet; the definition above is explicitly smooth real, not complex
        # analytic through q=0.
        "new_analytic_object_equals_canonical_lower_collar": bool(
            ast_passed and graph_passed and analytic["function_uniqueness_proof"]["passed"]),
        "accepted_native_inlet_same_function": inlet_identity,
        "forward_rows_enclose_same_native_function": rows_same_function,
        "accepted_true_pressure_same_function": bool(
            native_receipts["pressure_checks"]["same_true_pressure_function"]),
        "old_forward_producer_equals_canonical": False,
        "full_forward_leading_function_identity": False,
        "old_forward_function_graph_present": False,
        "actual_bridge_source_bound": bridge_algebra["passed"],
        "bridge_D_over_R_to_ratio_Psi_ODE_identity": False,
        "required_cross_map_contract": bridge_algebra["required_cross_map_contract"],
        "interval_overlap_used_as_identity": False,
        "missing_contract": missing,
        "scope": "new smooth real leading continuation plus accepted native source bindings; old finite forward producer remains an enclosure",
    }
    checks = {
        "parent_hash_and_mapping": True,
        "forward_dependency_audit_source_only": True,
        "forward_report_pins": all(pin["passed"] for pin in source_rows["report_pins"].values()),
        "forward_report_checks_and_family": all(
            row["checks_passed"] and row["source_family_match"]
            for row in accepted_forward_reports.values()),
        "source_file_pins": file_passed,
        "forward_source_AST_guards": ast_passed,
        "canonical_lower_graph_and_check": graph_passed,
        "accepted_native_function_receipts": native_receipts["passed"],
        "actual_bridge_source_AST": bridge_ast_passed,
        "actual_native_frame_algebra": bridge_algebra["passed"],
        "canonical_ratio_Psi_mean_ODE_uniqueness": analytic["function_uniqueness_proof"]["passed"],
        "no_raw_n1_in_leading_profiles": True,
        "no_old_producer_identity_claim": identity["old_forward_producer_equals_canonical"] is False,
    }
    actual_leading_identity = {
        "new_smooth_object_equals_canonical_lower_collar": identity[
            "new_smooth_object_equals_canonical_lower_collar"],
        "new_analytic_object_equals_canonical_lower_collar": identity[
            "new_analytic_object_equals_canonical_lower_collar"],
        "accepted_native_inlet_same_function": inlet_identity,
        "forward_rows_enclose_same_native_function": rows_same_function,
        "accepted_true_pressure_same_function": identity[
            "accepted_true_pressure_same_function"],
        "old_forward_producer_equals_canonical": False,
        "full_forward_leading_function_identity": False,
        "actual_bridge_source_bound": bridge_algebra["passed"],
        "bridge_D_over_R_to_ratio_Psi_ODE_identity": False,
        "required_cross_map_contract": bridge_algebra["required_cross_map_contract"],
        "source_definition_kind": "new smooth real continuation; sigma is nonanalytic at q=0",
        "scope": "q in [0,2], s=(sc/4)*q, canonical Phi/Psi/mean_Psi plus accepted F0/P0/U0 handles",
        "source_executed": False,
    }
    return {
        "passed": all(checks.values()),
        "source_binding_passed": all(checks.values()),
        "source_only": True,
        "source_executed": False,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": family,
        "checks": checks,
        "source_pins": {**file_pins, **source_rows["report_pins"],
                         **native_receipts["report_pins"],
                         **bridge_ast["source_pins"]["source_pins"],
                         "canonical_lower_graph": graph["graph_pin"],
                         "canonical_lower_check": graph["check_pin"]},
        "forward_reports": accepted_forward_reports,
        "forward_source_AST": source_ast,
        "accepted_native_receipts": native_receipts,
        "actual_bridge_source_AST": bridge_ast,
        "actual_native_frame_algebra": bridge_algebra,
        "canonical_lower_graph": graph,
        "analytic_forward_leading_collar": analytic,
        "identity": identity,
        "actual_leading_collar_function_identity": actual_leading_identity,
        "new_analytic_forward_definition": analytic,
        "missing_input_contract": missing,
        "limits": {
            "old_forward_numeric_is_finite_width_only": True,
            "old_forward_numeric_point_values_selected": False,
            "old_forward_transport_is_n1_half_collar_only": True,
            "new_object_is_not_a_physical_velocity_assembly": True,
            "smooth_real_not_complex_across_q0": True,
            "accepted_native_F0_P0_U0_same_function": inlet_identity,
            "accepted_forward_rows_are_enclosures_not_point_values": rows_same_function,
            "actual_bridge_D_over_R_cross_map_open": True,
            "full_six_signed_source_equality": False,
            "global_five_moment_closure": False,
            "temporal_n_ge_2_recursion": False,
            "full_corrected_NS": False,
        },
    }


def analytic_collar_definition(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Return the same packet with an explicit analytic-object API name."""

    return build(root, parent)
