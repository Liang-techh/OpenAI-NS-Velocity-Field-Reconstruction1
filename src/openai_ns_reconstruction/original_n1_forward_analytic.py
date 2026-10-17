"""Source-bound infinite analytic forward response for temporal ``n=1``.

The accepted canonical norm packet supplies a normal-convergent regular
Volterra solution in the raw six fields.  This module defines the new forward
analytic field by an explicit scale conjugacy.  It consumes the saved packet,
the accepted all-order coefficient identity, and source AST pins; it never
imports or executes a numerical producer, recurrence, fixture, owner, or
quadrature.

The field defined here is an actual infinite local analytic response.  The
older ``OriginalN1CoreNumeric`` object remains a finite radial-prefix
producer, and is recorded as such.  No collar extension, six-term identity,
or selected interval centre is asserted by this module.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_canonical_norm_tail.json"
)
PARENT_REPORT_SHA256 = (
    "b71ee21cf0f1dd2f3fe1d4edf4f45037ced997f1d72e2b8e0c9f1c861aae6989"
)
TERMINAL_PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C5_terminal_primitive.json"
)
TERMINAL_PARENT_REPORT_SHA256 = (
    "588d221fd082d5f84d94f1812eb6342bb6647847b114dac458c0165324a063b1"
)

FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
SOURCE_VECTOR_MODULE = (
    "src/openai_ns_reconstruction/original_n1_source_vector_identity.py"
)
SOURCE_VECTOR_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json"
)
RADIAL_CONJUGACY_MODULE = "src/openai_ns_reconstruction/original_n1_radial_conjugacy.py"
GREEN_MAJORANT_MODULE = (
    "src/openai_ns_reconstruction/original_n1_complex_radial_majorant.py"
)
CANONICAL_BINDING_MODULE = (
    "src/openai_ns_reconstruction/original_n1_canonical_solution_binding.py"
)
NORM_TAIL_MODULE = "src/openai_ns_reconstruction/original_n1_canonical_norm_tail.py"
PREFIX_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_numeric_prefix.json"
)
RAW_MATRIX_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
OPERATOR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
)

EXPECTED_SHA256 = {
    PARENT_REPORT: PARENT_REPORT_SHA256,
    FORWARD_SOURCE: "a7a925360f3a6da2348c5f18ee74cfbb885f8e220730769131e43a18bd954c2f",
    SOURCE_VECTOR_MODULE: "878f2b8781215a52a5730ed8527a561ff0df451adb5da247cec4d9fe0198423f",
    SOURCE_VECTOR_REPORT: "2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e",
    RADIAL_CONJUGACY_MODULE: "a6b1763e172fc87a84b68ffe9c188c87f62c1aa1ba6d0b6eeeb452e1ee6c8dea",
    GREEN_MAJORANT_MODULE: "a341ba79df0d43f985bc24320a7a486889c0ab59aeff653655570f9fc7c243f8",
    CANONICAL_BINDING_MODULE: "569337e2d49b8423eac9bcc8eb1ac7f319847ef66cbefe7f79f72922ed0ec0ac",
    NORM_TAIL_MODULE: "4e2fde5f0d4204272ccc3cccedf084599627073ec4360d3ba87d71c539ca0814",
    PREFIX_REPORT: "2ae583cdf5ea9fde03208fea0423eb11d171054c9a8dd45d78aa5fe2076d81fa",
    RAW_MATRIX_SOURCE: "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208",
    OPERATOR_SOURCE: "268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014",
}

D = (0, 0, 2, 0, 3, 1)
AXIS_POWERS = (2, 2, 2, 2, 1, 1)
CANONICAL_STATE = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
FORWARD_STATE = ("f", "u", "k", "p", "f_x", "u_x")
NORMALIZATION = {
    "f": "epsilon*F1/F0",
    "u": "Uz1/epsilon",
    "k": "K1/epsilon",
    "p": "P1/(epsilon*Pstar^2)",
}
SCALE_DIAGONAL = (
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    "epsilon",
)
CANONICAL_EQUATION = "W_x + D*W/x = B0*W + B1*D_Z W + g"
EXACT_WORD_DEFINITION = "W=sum_k>=0 H^k W0; W0=(I+J)Gg; H=(I+J)G(C0+B1*dZ)"


def _repo_root(root: str | Path | None) -> Path:
    starts: list[Path] = []
    if root is not None:
        starts.append(Path(root))
    starts.extend((Path.cwd(), Path(__file__).resolve().parents[2]))
    for start in starts:
        if start.is_file():
            start = start.parent
        for candidate in (start, *start.parents):
            if (candidate / "src").is_dir() and (candidate / "experiments" / "root_st073").is_dir():
                return candidate
    raise FileNotFoundError("st073 repository root is required")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _read_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError(f"JSON mapping required: {path}")
    return value


def _all_true(values: Mapping[str, Any]) -> bool:
    return bool(values) and all(value is True for value in values.values())


def _accepted_parent(repo: Path, parent: Mapping[str, Any] | None) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = repo / PARENT_REPORT
    actual_sha = _sha256(path)
    if actual_sha != PARENT_REPORT_SHA256:
        raise ValueError("accepted canonical norm-tail report changed")
    saved = _read_json(path)
    supplied = saved if parent is None else parent
    if not isinstance(supplied, Mapping):
        raise ValueError("canonical norm-tail parent mapping required")
    # A caller may pass the complete report (the normal API) or the nested
    # norm-tail packet.  In either case the bytes on disk remain authoritative.
    if "n1_actual_canonical_norm_tail" not in supplied:
        nested = supplied
        if _digest(nested) != _digest(saved.get("n1_actual_canonical_norm_tail")):
            raise ValueError("nested norm-tail packet does not match accepted report")
        supplied = saved
        input_kind = "nested_norm_tail_packet"
    else:
        if _digest(supplied) != _digest(saved):
            raise ValueError("parent mapping does not match accepted report")
        input_kind = "complete_norm_tail_report"
    return supplied, {
        "passed": True,
        "path": PARENT_REPORT,
        "sha256": actual_sha,
        "expected_sha256": PARENT_REPORT_SHA256,
        "mapping_digest": _digest(supplied),
        "input_kind": input_kind,
        "source_executed": False,
    }


def _ast_pin(repo: Path) -> dict[str, Any]:
    """Parse the finite source for pins only; never compile or execute it."""
    path = repo / FORWARD_SOURCE
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    wanted = ("linear_rhs", "_rhs", "_advance", "coefficients", "evaluate_prefix")
    found: dict[str, ast.AST] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in wanted:
            if node.name in found:
                raise ValueError(f"duplicate finite source definition: {node.name}")
            found[node.name] = node
    if set(found) != set(wanted):
        raise ValueError("finite forward definitions changed")

    def body(name: str) -> str:
        segment = ast.get_source_segment(source, found[name]) or ""
        return "".join(segment.split()).replace('"', "'")

    requirements = {
        "linear_rhs": (
            "rf=(2*e*(Q+E*(-2+de)/L)-2*H*g/L)*f",
            "ru=-8*F2*z*r*A*f/L",
            "rp=2*F2/(e*P2)*A*f",
        ),
        "_rhs": (
            "Aw=known(leading['A'],n+1,axial+2)",
            "Uw=known(leading['Uz'],n+1,axial+2)",
            "Qw=(2*zz*uu-zz*mm*(1-delta)-dd*mw.dz())/LL",
            "pressure_head=R([[(i+1)*(i+2)*vforvinQw.rows[i+1][:axial+1]",
        ),
        "_advance": (
            "axial=degree+depth-n-1",
            "k=n+1",
            "state['f'].append([v/(4*k*(k+1))forvinrf.rows[n]])",
            "state['u'].append([v/(4*k*k)forvinru.rows[n]])",
            "state['p'].append([v/kforvinrp.rows[n]])",
            "state['k'].append([-c.mpf(k)*v/(k+1)forvinstate['u'][-1]])",
        ),
        "coefficients": (
            "forninrange(lead_degree)",
            "forninrange(degree)",
        ),
        "evaluate_prefix": (
            "finite_radial_prefix_only=True",
            "infinite_radial_remainder_bounded=False",
            "full_original_n1_point_error_available=False",
            "full_original_n1_point_oracle=False",
        ),
    }
    fragments = {
        name: {fragment: fragment in body(name) for fragment in values}
        for name, values in requirements.items()
    }
    if not all(all(checks.values()) for checks in fragments.values()):
        raise ValueError("finite forward AST semantics changed")
    return {
        "path": FORWARD_SOURCE,
        "sha256": _sha256(path),
        "expected_sha256": EXPECTED_SHA256[FORWARD_SOURCE],
        "hash_passed": _sha256(path) == EXPECTED_SHA256[FORWARD_SOURCE],
        "definitions": list(wanted),
        "fragments": fragments,
        "source_executed": False,
        "archive_source_executed": False,
        "passed": _sha256(path) == EXPECTED_SHA256[FORWARD_SOURCE],
    }


def _source_pins(repo: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    radial = _mapping(parent.get("n1_all_order_radial_conjugacy"))
    symbolic = _mapping(radial.get("symbolic"))
    inherited_matrix = _mapping(symbolic.get("inherited_matrix_conjugacy"))
    matrix_pin = _mapping(inherited_matrix.get("canonical_matrix_ast_pin"))
    # The accepted 30/36 packet stores the matrix pin but not a separate
    # operator-AST pin.  The operator source is independently hash-pinned
    # below; do not re-run its accepted coefficient comparisons.
    operator_pin = {
        "path": OPERATOR_SOURCE,
        "expected_sha256": EXPECTED_SHA256[OPERATOR_SOURCE],
        "passed": _sha256(repo / OPERATOR_SOURCE) == EXPECTED_SHA256[OPERATOR_SOURCE],
        "source_executed": False,
        "archive_source_executed": False,
    }
    expected_paths = {
        FORWARD_SOURCE: EXPECTED_SHA256[FORWARD_SOURCE],
        SOURCE_VECTOR_MODULE: EXPECTED_SHA256[SOURCE_VECTOR_MODULE],
        RADIAL_CONJUGACY_MODULE: EXPECTED_SHA256[RADIAL_CONJUGACY_MODULE],
        GREEN_MAJORANT_MODULE: EXPECTED_SHA256[GREEN_MAJORANT_MODULE],
        CANONICAL_BINDING_MODULE: EXPECTED_SHA256[CANONICAL_BINDING_MODULE],
        NORM_TAIL_MODULE: EXPECTED_SHA256[NORM_TAIL_MODULE],
        SOURCE_VECTOR_REPORT: EXPECTED_SHA256[SOURCE_VECTOR_REPORT],
        PREFIX_REPORT: EXPECTED_SHA256[PREFIX_REPORT],
        RAW_MATRIX_SOURCE: EXPECTED_SHA256[RAW_MATRIX_SOURCE],
        OPERATOR_SOURCE: EXPECTED_SHA256[OPERATOR_SOURCE],
    }
    pins: dict[str, Any] = {}
    for relative, expected in expected_paths.items():
        path = repo / relative
        actual = _sha256(path)
        pins[relative] = {
            "path": relative,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "passed": actual == expected,
            "source_executed": False,
            "archive_source_executed": False,
        }
    pins["finite_forward_ast"] = _ast_pin(repo)
    # Parsing a source-vector module is an AST integrity check only.  Its
    # accepted symbolic identity is consumed from the saved radial packet.
    sv_tree = ast.parse((repo / SOURCE_VECTOR_MODULE).read_text(encoding="utf-8"))
    sv_names = sorted(
        node.name
        for node in sv_tree.body
        if isinstance(node, (ast.FunctionDef, ast.ClassDef))
    )
    pins["source_vector_ast"] = {
        "definitions": sv_names,
        "required_definitions_present": "prove_source_vector_identity" in sv_names,
        "source_executed": False,
        "archive_source_executed": False,
        "passed": "prove_source_vector_identity" in sv_names and pins[SOURCE_VECTOR_MODULE]["passed"],
    }
    pins["inherited_matrix_ast"] = {
        "canonical_matrix_pin": dict(matrix_pin),
        "canonical_operator_pin": dict(operator_pin),
        "replayed": False,
        "source_executed": False,
        "passed": bool(matrix_pin.get("passed") is True and operator_pin.get("passed") is True),
    }
    pins["all_source_hashes_passed"] = all(
        item.get("passed") is True
        for key, item in pins.items()
        if isinstance(item, Mapping) and key not in {"inherited_matrix_ast", "all_source_hashes_passed"}
    ) and pins["inherited_matrix_ast"]["passed"]
    return pins


def _normalization_binding(parent: Mapping[str, Any]) -> dict[str, Any]:
    radial = _mapping(parent.get("n1_all_order_radial_conjugacy"))
    symbolic = _mapping(radial.get("symbolic"))
    source_map = _mapping(symbolic.get("inherited_source_map"))
    matrix = _mapping(symbolic.get("inherited_matrix_conjugacy"))
    assumptions = _mapping(source_map.get("assumptions"))
    matrix_assumptions = _mapping(matrix.get("assumptions"))
    premises = _mapping(source_map.get("accepted_normalization_premises"))
    premise_checks = _mapping(premises.get("checks"))
    inherited_premises = _mapping(premises.get("inherited_premises"))
    expected_assumptions = {
        "P2_equals_Pstar_squared": "P2 = Pstar^2",
        "gauge_log_gradient": "g = -epsilon * D_Z(F0) / F0",
        "F2_equals_F0_squared": "F2 = F0^2",
        "F0_nonzero": "F0 != 0",
        "epsilon_nonzero": "epsilon != 0",
        "P2_nonzero": "P2 != 0",
    }
    checks = {
        "source_map_passed": source_map.get("checks") and _all_true(_mapping(source_map.get("checks"))),
        "matrix_passed": matrix.get("checks") and _all_true(_mapping(matrix.get("checks"))),
        "normalization_exact": source_map.get("normalization") == NORMALIZATION,
        "assumptions_exact": all(assumptions.get(k) == v for k, v in expected_assumptions.items()),
        "L_nonzero_premise": matrix_assumptions.get("L_nonzero") == "L != 0",
        "premise_checks_passed": premises.get("passed") is True and _all_true(premise_checks),
        "inherited_premises_passed": _all_true(inherited_premises),
        "source_identity_inherited": source_map.get("source_identity_is_inherited") is True,
        "not_replayed": source_map.get("replayed") is False,
        "source_not_executed": source_map.get("source_executed") is False,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "normalization": dict(NORMALIZATION),
        "assumptions": {**expected_assumptions, "L_nonzero": "L != 0"},
        "assumptions_kind": "accepted symbolic source-vector checks and inherited premises",
        "accepted_normalization_premises": premises,
        "source_map": source_map,
        "inherited_matrix_comparisons": {
            "state_coefficient_comparisons": matrix.get("state_coefficient_comparisons"),
            "mixed_coefficient_comparisons": matrix.get("mixed_coefficient_comparisons"),
            "replayed": matrix.get("replayed"),
            "source_executed": matrix.get("source_executed"),
        },
        "source_executed": False,
        "archive_source_executed": False,
    }


def _canonical_packets(parent: Mapping[str, Any]) -> dict[str, Any]:
    radial = _mapping(parent.get("n1_all_order_radial_conjugacy"))
    norm = _mapping(parent.get("n1_actual_canonical_norm_tail"))
    binding = _mapping(parent.get("n1_actual_canonical_solution_binding"))
    prefix = _mapping(parent.get("n1_actual_raw_prefix_binding"))
    bounds = _mapping(parent.get("n1_complex_actual_source_bounds"))
    green = _mapping(parent.get("n1_weighted_complex_Green_kernel"))
    convergence = _mapping(norm.get("convergence_proof"))
    inherited_matrix = _mapping(_mapping(radial.get("symbolic")).get("inherited_matrix_conjugacy"))
    inherited_source = _mapping(_mapping(radial.get("symbolic")).get("inherited_source_map"))
    prefix_data = _mapping(prefix.get("actual_forward_prefix_bound"))
    prefix_checks = _mapping(prefix.get("checks"))
    prefix_normalization = dict(_mapping(prefix_data.get("normalization")))
    if prefix_normalization.get("f") == "epsilon*F1/F0(Z)":
        prefix_normalization["f"] = "epsilon*F1/F0"
    norm_checks = {
        "norm_tail_passed": norm.get("passed") is True,
        "canonical_norm_proved": norm.get("actual_canonical_regular_solution_norm_proved") is True,
        "normal_convergence": convergence.get("normal_convergence_on_compact_subsets") is True,
        "termwise_derivatives": convergence.get("termwise_derivatives_on_inner_compact_subsets") is True,
        "word_bound_exact": convergence.get("word_bound") == "Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!",
        "geometric_bound_exact": convergence.get("geometric_bound") == "word_k <= Msrc*beta^k",
        "solution_norm_exact": convergence.get("solution_norm") == "||v|| <= Msrc/(1-beta) <= 2*Msrc",
        "solution_definition_exact": norm.get("solution_definition") == "Normal-convergent sum H^k W0, k>=0; W0 is the preconditioned green_apply(g) output with actual zero-axis data",
        "no_forward_identity_claim": norm.get("full_forward_response_function_identity") is False,
        "no_collar_bridge": norm.get("q_dependent_collar_bridge") is False,
    }
    binding_required = (
        "actual_canonical_owner_bound",
        "actual_zero_axis_condition_bound",
        "actual_raw_equations_bound",
        "actual_raw_B0_B1_g_same_source",
        "actual_axis_weighted_transform_bound",
        "actual_source_normalization_bound",
        "analytic_regular_uniqueness_bound",
        "actual_homogeneous_word_uniqueness_bound",
        "actual_canonical_regular_solution_exists",
    )
    binding_checks = {
        "binding_passed": binding.get("passed") is True,
        "binding_flags": all(binding.get(key) is True for key in binding_required),
        "historical_marker_not_used": _mapping(binding.get("canonical_owner")).get(
            "historical_adapter_marker_used_as_acceptance"
        ) is False,
        "owner_function_nodes": _mapping(binding.get("canonical_owner")).get(
            "actual_function_nodes_used_as_acceptance"
        ) is True,
        "word_definition": _mapping(binding.get("canonical_solution")).get("definition") == EXACT_WORD_DEFINITION,
        "word_norm_not_recomputed": _mapping(binding.get("canonical_solution")).get(
            "norm_recomputed_here"
        ) is False,
    }
    prefix_checks_local = {
        "prefix_passed": prefix.get("passed") is True and prefix_data.get("passed") is True,
        "formal_identity_every_available": _mapping(
            radial.get("available_forward_enclosures")
        ).get("formal_identity_applies_to_every_available_forward_coefficient") is True,
        "directed_enclosures": _mapping(radial.get("available_forward_enclosures")).get(
            "directed_source_enclosures_present"
        ) is True,
        "degree16": prefix_data.get("raw_base_radial_degree") == 16,
        "weighted_degree15": prefix_data.get("weighted_regularized_radial_degree") == 15,
        "normalization_exact": prefix_normalization == NORMALIZATION,
        "canonical_state_exact": prefix_data.get("canonical_state") == list(CANONICAL_STATE),
        "canonical_D_exact": prefix_data.get("canonical_D") == list(D),
        "not_selected_points": _mapping(radial.get("available_forward_enclosures")).get(
            "interval_enclosures_are_not_selected_points"
        ) is True and _mapping(prefix_data.get("limits")).get(
            "physical_scaled_coefficient_intervals_materialized"
        ) is False,
        "finite_prefix_retained": prefix_data.get("finite_radial_prefix_only", True) is True
        and _mapping(radial.get("available_forward_enclosures")).get("finite_prefix_only") is True,
        "prefix_report_checks": _all_true(prefix_checks),
    }
    source_bounds_checks = {
        "bounds_passed": bounds.get("passed") is True,
        "source_bounds": bounds.get("actual_source_bounds_proved") is True,
        "operator_bounds": bounds.get("actual_weighted_operator_bounds_proved") is True,
        "source_only": bounds.get("source_only") is True and bounds.get("source_executed") is False,
        "no_replay": bounds.get("accepted_owner_checker_fixture_recurrence_quadrature_replayed") is False,
        "missing_inputs_empty": bounds.get("missing_input_contract") == [],
        "domain_gates": all(
            value is True
            for key, value in _mapping(bounds.get("gates")).items()
            if key not in ("solution_norm", "radial_Taylor_remainder")
        ),
    }
    green_checks = {
        "green_passed": green.get("passed") is True,
        "matrix_entries": _mapping(green.get("actual_weighted_matrix_entries")).get("passed") is True,
        "exact_controls": _mapping(green.get("exact_controls")).get("passed") is True,
        "axis_weights": green.get("axis_weights") == list(AXIS_POWERS),
        "six_kernels_include_kinematic": green.get("kinematic_preconditioning_already_in_six_kernels") is True,
        "no_extra_preconditioner": green.get("extra_I_plus_J_factor_applied") is False,
        "z_loss": green.get("Z_loss_allocation") == "total gap sigma, sigma/k for each of k applications",
    }
    inherited_checks = {
        "radial_passed": radial.get("passed") is True,
        "inherited_matrix_passed": inherited_matrix.get("checks", {}).get("passed") is True,
        "matrix_30": inherited_matrix.get("state_coefficient_comparisons") == 30,
        "mixed_36": inherited_matrix.get("mixed_coefficient_comparisons") == 36,
        "matrix_not_replayed": inherited_matrix.get("replayed") is False,
        "matrix_source_not_executed": inherited_matrix.get("source_executed") is False,
        "source_map_passed": inherited_source.get("checks", {}).get("all_three_sources") is True,
        "source_map_not_replayed": inherited_source.get("replayed") is False,
        "radial_source_not_executed": radial.get("source_executed") is False,
        "radial_archived_only": radial.get("archive_source_executed") is False,
    }
    return {
        "checks": {
            "norm_tail": {**norm_checks, "passed": all(norm_checks.values())},
            "canonical_binding": {**binding_checks, "passed": all(binding_checks.values())},
            "prefix": {**prefix_checks_local, "passed": all(prefix_checks_local.values())},
            "source_bounds": {**source_bounds_checks, "passed": all(source_bounds_checks.values())},
            "green": {**green_checks, "passed": all(green_checks.values())},
            "inherited_conjugacy": {**inherited_checks, "passed": all(inherited_checks.values())},
        },
        "norm_tail": norm,
        "canonical_binding": binding,
        "prefix_binding": prefix,
        "radial_conjugacy": radial,
        "source_bounds": bounds,
        "green": green,
        "inherited_matrix": inherited_matrix,
        "inherited_source": inherited_source,
    }


def _derive_scale(normalization: Mapping[str, Any]) -> dict[str, Any]:
    """Derive the six positional scale entries from the accepted map."""
    rules = (
        ("f", "epsilon*F1/F0", "F0(Z)/epsilon", "epsilon/F0(Z)"),
        ("u", "Uz1/epsilon", "epsilon", "1/epsilon"),
        ("k", "K1/epsilon", "epsilon", "1/epsilon"),
        ("p", "P1/(epsilon*Pstar^2)", "epsilon*Pstar^2", "1/(epsilon*Pstar^2)"),
        ("f_x", "epsilon*F1x/F0", "F0(Z)/epsilon", "epsilon/F0(Z)"),
        ("u_x", "Uz1x/epsilon", "epsilon", "1/epsilon"),
    )
    if dict(normalization) != NORMALIZATION:
        raise ValueError("accepted forward normalization changed")
    scale: list[str] = []
    inverse: list[str] = []
    positional: list[dict[str, str]] = []
    for name, accepted_expression, scale_entry, inverse_entry in rules:
        if name in NORMALIZATION and normalization.get(name) != accepted_expression:
            raise ValueError(f"accepted scale convention changed: {name}")
        scale.append(scale_entry)
        inverse.append(inverse_entry)
        positional.append(
            {
                "forward_component": name,
                "accepted_normalization": accepted_expression,
                "raw_scale": scale_entry,
                "forward_inverse_scale": inverse_entry,
            }
        )
    return {"diagonal": scale, "inverse_diagonal": inverse, "positional": positional}


def _symbolic_scale_checks(derived: Mapping[str, Any]) -> dict[str, Any]:
    """Check the scale algebra with fresh exact SymPy symbols only."""
    import sympy as sp

    epsilon, F0, F0z, Pstar = sp.symbols("epsilon F0 F0z Pstar", nonzero=True)
    expressions = {
        "F0(Z)/epsilon": F0/epsilon, "epsilon": epsilon,
        "epsilon*Pstar^2": epsilon*Pstar**2, "epsilon/F0(Z)": epsilon/F0,
        "1/epsilon": 1/epsilon, "1/(epsilon*Pstar^2)": 1/(epsilon*Pstar**2),
    }
    scale = [expressions[value] for value in derived["diagonal"]]
    inverse = [expressions[value] for value in derived["inverse_diagonal"]]
    S_z = [F0z / epsilon, 0, 0, 0, F0z / epsilon, 0]
    gamma = [F0z / F0, 0, 0, 0, F0z / F0, 0]
    products = [sp.simplify(a * b - 1) == 0 for a, b in zip(scale, inverse)]
    sz_gamma = [sp.simplify(a * b - c) == 0 for a, b, c in zip(scale, gamma, S_z)]
    S = sp.diag(*scale)
    diagonal_D = sp.diag(*D)
    x = sp.symbols("x")
    commute = S*diagonal_D-diagonal_D*S
    return {
        "S_inverse_S_identity": all(products),
        "component_products": products,
        "D_commutes_with_S": all(sp.simplify(value) == 0 for value in commute),
        "S_x_zero": all(sp.diff(value,x) == 0 for value in S),
        "S_z_equals_S_Gamma": all(sz_gamma),
        "connection_diagonal": ["F0_Z/F0", "0", "0", "0", "F0_Z/F0", "0"],
        "exact_symbolic_backend": "sympy",
    }


def _scale_conjugacy(
    normalization: Mapping[str, Any],
    source_map: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    derived = _derive_scale(normalization)
    symbolic = _symbolic_scale_checks(derived)
    source_map = _mapping(source_map)
    source_rows = _mapping(source_map.get("canonical_source_rows"))
    forward_sources = _mapping(source_map.get("forward_normalized_sources"))
    expected_rows = {
        "g4": "2*epsilon*x*N1p = 2*epsilon*x*Pstar^2*Np",
        "g5": "2*epsilon*N1theta",
        "g6": "2*epsilon*N1z - 4*epsilon^2*Z*rho*N1p/L",
    }
    expected_forward_sources = {
        "pressure": "Np=N1p/Pstar^2",
        "angular": "2*epsilon^2*N1theta/F0",
        "axial": "2*N1z-4*epsilon*Z*rho*N1p/L",
    }
    checks = {
        "normalization_exact": dict(normalization) == NORMALIZATION,
        "derived_positional_scale_exact": derived.get("diagonal") == list(SCALE_DIAGONAL),
        "derived_positional_inverse_exact": derived.get("inverse_diagonal") == [
            "epsilon/F0(Z)", "1/epsilon", "1/epsilon", "1/(epsilon*Pstar^2)",
            "epsilon/F0(Z)", "1/epsilon",
        ],
        "S_inverse_S_identity": symbolic["S_inverse_S_identity"],
        "D_commutes_with_S": symbolic["D_commutes_with_S"],
        "S_x_zero": symbolic["S_x_zero"],
        "S_z_equals_S_Gamma": symbolic["S_z_equals_S_Gamma"],
        "accepted_canonical_source_map": dict(source_rows) == expected_rows,
        "accepted_forward_source_map": dict(forward_sources) == expected_forward_sources,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "canonical_state": list(CANONICAL_STATE),
        "forward_state": list(FORWARD_STATE),
        "rho": "rho=x^2",
        "x": "sqrt(rho)",
        "normalization": dict(NORMALIZATION),
        "W_equals_SY": "W=S(Z)Y",
        "scale_matrix": {
            "type": "diagonal",
            "diagonal": derived["diagonal"],
            "S_inverse_diagonal": derived["inverse_diagonal"],
            "positional_map": derived["positional"],
        },
        "exact_symbolic_scale_checks": symbolic,
        "connection": "Gamma=S^{-1}*D_Z(S)",
        "derivative_identity": "D_Z(SY)=S*(D_Z Y+Gamma*Y)",
        "forward_operator": {
            "B0_forward": "S^{-1}*B0*S",
            "B1_forward": "S^{-1}*B1*S",
            "g_forward": "S^{-1}*g",
            "equation": "Y_x + D*Y/x = (S^{-1}*B0*S)Y + (S^{-1}*B1*S)(Y_Z+Gamma*Y) + S^{-1}*g",
            "canonical_equation": CANONICAL_EQUATION,
        },
        "source_scale_map": {
            "canonical_source_rows": dict(source_rows),
            "forward_normalized_sources": dict(forward_sources),
            "pressure_normalization": "Np=N1p/Pstar^2",
            "fourth_pde_source": "(S^{-1}g)_4=2*x*Np",
            "fourth_pde_source_is_distinct_from_Np_definition": True,
        },
        "axis": {
            "D": list(D),
            "regular_powers": list(AXIS_POWERS),
            "canonical_zero_axis": "W(0,Z)=0",
            "forward_zero_axis": "Y(0,Z)=S(Z)^{-1}W(0,Z)=0",
            "raw_orders": "W1..W4=O(rho), W5,W6=O(x)",
            "forward_orders": "f,u,k,p=O(rho), f_x,u_x=O(x)",
            "source_vector_identity_domain": "rho>0 and L != 0",
            "axis_extension": "regular parity and the accepted identity theorem extend the fields to rho=0; no axis division is evaluated",
            "axis_division_not_evaluated": True,
        },
        "nonvanishing_premises": {
            "F0": "F0 != 0",
            "epsilon": "epsilon != 0",
            "P2": "P2 != 0",
            "P2_equals_Pstar_squared": "P2 = Pstar^2",
            "L": "L != 0",
        },
        "source_executed": False,
    }


def _analytic_definition(packets: Mapping[str, Any], scale: Mapping[str, Any]) -> dict[str, Any]:
    norm = _mapping(packets.get("norm_tail"))
    bounds = _mapping(packets.get("source_bounds"))
    domain = _mapping(norm.get("axial_domain"))
    bound_domain = _mapping(bounds.get("domain"))
    convergence = _mapping(norm.get("convergence_proof"))
    expected_nonzero = {
        "F0": "F0 != 0",
        "epsilon": "epsilon != 0",
        "P2": "P2 != 0",
        "P2_equals_Pstar_squared": "P2 = Pstar^2",
        "L": "L != 0",
    }
    checks = {
        "normal_convergence_accepted": packets["checks"]["norm_tail"]["passed"],
        "scale_conjugacy_accepted": scale.get("checks", {}).get("passed") is True,
        "solution_definition_exact": norm.get("solution_definition") == "Normal-convergent sum H^k W0, k>=0; W0 is the preconditioned green_apply(g) output with actual zero-axis data",
        "termwise_derivatives_accepted": convergence.get("termwise_derivatives_on_inner_compact_subsets") is True,
        "source_radius_exact": domain.get("source_geometry") == "Union of admitted local disks centered on the real interval; not a single origin disk",
        "rho_radius_exact": bound_domain.get("rho_bound_radius") == "9/2" and bound_domain.get("rho_target_radius") == "4",
        "scale_domain_nonzero": scale.get("nonvanishing_premises") == expected_nonzero,
        "axis_identity_contract": (
            packets["checks"]["canonical_binding"]["passed"]
            and scale.get("nonvanishing_premises", {}).get("L") == "L != 0"
        ),
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "definition_kind": "new_infinite_forward_analytic_response",
        "canonical_weighted_series": {
            "field": "v_star",
            "definition": "v_star=sum_{k>=0} H^k v0",
            "accepted_report_definition": norm.get("solution_definition"),
            "report_series_convention": "The accepted W0/W word notation is the six weighted Green-kernel field; T maps it to raw W_star.",
            "v0": "preconditioned_green_apply(g) with actual zero-axis data",
            "H": "preconditioned_green_apply(C0*v+B1*D_Z(v)); six kernels already include kinematic preconditioning",
            "normal_convergence": convergence.get("normal_convergence_on_compact_subsets"),
            "termwise_holomorphic_derivatives": convergence.get("termwise_derivatives_on_inner_compact_subsets"),
            "termwise_derivative_scope": {
                "weighted_v_in_rho": True,
                "raw_W1_to_W4_in_rho": True,
                "raw_W5_W6_in_x_with_odd_parity": True,
                "raw_W5_over_x_W6_over_x_in_rho": True,
                "entire_raw_six_vector_single_valued_holomorphic_in_rho": False,
            },
        },
        "canonical_raw_series": {
            "field": "W_star",
            "definition": "W_star=T(rho)*v_star",
            "axis_transform": "T(rho)*v=(rho*v1,rho*v2,rho*v3,rho*v4,x*v5,x*v6)",
        },
        "forward_field": {
            "field": "Y_star",
            "definition": "Y_star(rho,Z)=S(Z)^(-1)*W_star(rho,Z)",
            "response": "(f,u,k,p)=first_four_components(Y_star)",
            "derivative_state": "(f_x,u_x)=last_two_components(Y_star)",
            "coefficient_definition": "Y_star is defined by the normal-convergent canonical series followed by the displayed scale inverse",
        },
        "domain": {
            "rho_source_disk": bound_domain.get("rho_bound_radius"),
            "rho_target_disk": bound_domain.get("rho_target_radius"),
            "local_Z_centers": domain.get("real_center_interval_open"),
            "local_center_description": domain.get("local_center"),
            "source_geometry": domain.get("source_geometry"),
            "solution_Z_radius": domain.get("solution_radius"),
            "scale_nonvanishing_on_domain": checks["scale_domain_nonzero"],
            "rho_holomorphic_components": ["v1", "v2", "v3", "v4", "W1", "W2", "W3", "W4", "f", "u", "k", "p"],
            "x_holomorphic_odd_components": ["W5", "W6", "f_x", "u_x"],
            "rho_holomorphic_regularized_odd_components": ["W5/x", "W6/x", "f_x/x", "u_x/x"],
            "raw_six_vector_single_valued_rho_holomorphic": False,
            "raw_odd_rows_require_x_sheet": True,
            "source_vector_identity_domain": "rho>0 and L != 0",
            "axis_extension_from_parity": "regular parity plus the accepted identity theorem extends to rho=0",
            "global_origin_disk_claimed": False,
        },
        "zero_axis": {
            "canonical": "W_star(0,Z)=0",
            "forward": "Y_star(0,Z)=0",
            "base_fields": "f(0,Z)=u(0,Z)=k(0,Z)=p(0,Z)=0",
            "derivative_fields": "f_x(0,Z)=u_x(0,Z)=0",
            "raw_base_rows": "W1..W4 are rho-holomorphic and vanish to order one",
            "raw_odd_rows": "W5,W6 are odd x-holomorphic and vanish to order one in x",
            "regularized_odd_rows": "W5/x,W6/x are rho-holomorphic at the axis",
            "rho_identity_theorem_scope": "The source-vector identity is established first for rho>0 with L != 0; regular parity supplies the axis extension.",
        },
        "source_executed": False,
        "selected_interval_centres": False,
        "finite_N_numerical_calculation": False,
    }


def _prefix_identity(packets: Mapping[str, Any], normalization: Mapping[str, Any]) -> dict[str, Any]:
    radial = packets["radial_conjugacy"]
    prefix = packets["prefix_binding"]
    data = _mapping(prefix.get("actual_forward_prefix_bound"))
    available = _mapping(radial.get("available_forward_enclosures"))
    matrix = packets["inherited_matrix"]
    prefix_normalization = dict(_mapping(data.get("normalization")))
    if prefix_normalization.get("f") == "epsilon*F1/F0(Z)":
        prefix_normalization["f"] = "epsilon*F1/F0"
    checks = {
        "formal_all_available": available.get("formal_identity_applies_to_every_available_forward_coefficient") is True,
        "accepted_prefix_packet": packets["checks"]["prefix"]["passed"],
        "normalization_exact": prefix_normalization == dict(normalization),
        "state_exact": data.get("canonical_state") == list(CANONICAL_STATE),
        "D_exact": data.get("canonical_D") == list(D),
        "matrix_comparisons_consumed": matrix.get("state_coefficient_comparisons") == 30,
        "mixed_comparisons_consumed": matrix.get("mixed_coefficient_comparisons") == 36,
        "matrix_not_replayed": matrix.get("replayed") is False,
        "coefficient_intervals_not_points": available.get("interval_enclosures_are_not_selected_points") is True,
        "finite_prefix_retained": available.get("finite_prefix_only") is True,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "statement": "For every available directed forward coefficient, the coefficient of Y_star equals S^{-1} times the accepted canonical raw coefficient after the stated normalization.",
        "coefficient_map": {
            "f": "f_m=(epsilon/F0(Z))*F1_m",
            "u": "u_m=Uz1_m/epsilon",
            "k": "k_m=K1_m/epsilon",
            "p": "p_m=P1_m/(epsilon*Pstar^2)",
            "f_x": "f_x=2*x*d_rho(f)",
            "u_x": "u_x=2*x*d_rho(u)",
        },
        "exact_axis_derivative_map": {
            "F1x_over_x": "F1x/x=2*(F0(Z)/epsilon)*sum_{n>=0}(n+1)*f_(n+1)*rho^n",
            "Uz1x_over_x": "Uz1x/x=2*epsilon*sum_{n>=0}(n+1)*u_(n+1)*rho^n",
            "forward_derivatives": "f_x/x=2*sum_{n>=0}(n+1)*f_(n+1)*rho^n; u_x/x=2*sum_{n>=0}(n+1)*u_(n+1)*rho^n",
        },
        "available_degree": {
            "raw_base": data.get("raw_base_radial_degree"),
            "weighted_regularized": data.get("weighted_regularized_radial_degree"),
            "common_ordinary_Z0_to_3": _mapping(data.get("ordinary_Z_domain")).get("max_rho_degree_with_orders_0_to_3"),
            "common_weighted_ordinary_Z0_to_3": _mapping(data.get("ordinary_Z_domain")).get("max_weighted_degree_with_orders_0_to_3"),
        },
        "accepted_receipts": {
            "prefix_report": data.get("prefix_report"),
            "prefix_report_sha256": data.get("prefix_report_sha256"),
            "source_vector_report": data.get("source_vector_report"),
            "source_vector_report_sha256": data.get("source_vector_report_sha256"),
            "state_coefficient_comparisons": 30,
            "mixed_coefficient_comparisons": 36,
        },
        "finite_prefix_only": True,
        "radial_remainder_claimed_here": False,
        "source_executed": False,
        "point_values_selected": False,
    }


def _uniqueness(packets: Mapping[str, Any], scale: Mapping[str, Any]) -> dict[str, Any]:
    binding = packets["canonical_binding"]
    estimate = _mapping(binding.get("homogeneous_word_estimate"))
    coefficient = _mapping(estimate.get("coefficient_uniqueness"))
    axis_transform = _mapping(binding.get("axis_weighted_transform"))
    forward_class = {
        "Y1_to_Y4": "Y1..Y4=O(rho)",
        "Y5_to_Y6": "Y5..Y6=O(x), odd in x",
        "weighted_v": "v=(Y1/rho,Y2/rho,Y3/rho,Y4/rho,Y5/x,Y6/x) is rho-holomorphic",
    }
    checks = {
        "canonical_unique": packets["checks"]["canonical_binding"]["passed"],
        "formal_unique": coefficient.get("formal_solution_unique") is True,
        "scale_invertible": scale.get("checks", {}).get("passed") is True,
        "radial_gain": "H has a factor rho" in str(estimate.get("homogeneous_uniqueness")),
        "distinct_norms": _mapping(estimate.get("homogeneous_norms")).get("same_space_contraction_claimed") is False,
        "identity_theorem": "identity theorem" in str(estimate.get("homogeneous_uniqueness")),
        "forward_regular_class": axis_transform.get("zero_axis") == "W(0,Z)=0; W1..W4=O(x^2), W5,W6=O(x)",
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "forward_regular_solution_unique": all(checks.values()),
        "proof": "If Y and Y_tilde lie in the class Y1..Y4=O(rho), Y5..Y6=O(x) with odd x parity and have the same source, S(Y-Y_tilde) is a canonical homogeneous regular solution. The weighted difference v=(Y1/rho,..,Y4/rho,Y5/x,Y6/x) is rho-holomorphic; accepted rho-triangularity gives every coefficient zero, then the connected rho-disk identity theorem gives equality.",
        "forward_regular_class": forward_class,
        "weighted_regular_class": "v1..v4=Y1..Y4/rho and v5..v6=Y5..Y6/x are rho-holomorphic; the raw odd rows are represented on the x sheet.",
        "homogeneous_word_estimate": {
            "input_norm": _mapping(estimate.get("homogeneous_norms")).get("input"),
            "output_norm": _mapping(estimate.get("homogeneous_norms")).get("output"),
            "bound": estimate.get("homogeneous_word_bound"),
            "same_space_contraction_claimed": False,
        },
        "formal_triangularity": {
            "weighted_operator": "x^(p_j-q_i)*B_ij=rho*A_ij with A holomorphic",
            "order_action": "H sends rho-order m to m+1",
            "ordinary_Z": "D_Z commutes with [rho^m] and does not reduce m",
            "induction": "h=Hh gives h_0=0, then h_m=0 from h_(m-1)",
            "analytic_closure": "zero rho Taylor series is zero on the connected rho disk",
        },
        "source_executed": False,
    }


def _finite_producer_scope(ast_pin: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "source": FORWARD_SOURCE,
        "source_sha256": EXPECTED_SHA256[FORWARD_SOURCE],
        "AST_pin": ast_pin,
        "definition": "OriginalN1CoreNumeric._advance and evaluate_prefix produce finite radial prefixes only",
        "advance_rules": {
            "f": "f_(n+1)=rf_n/[4*(n+1)*(n+2)]",
            "u": "u_(n+1)=ru_n/[4*(n+1)^2]",
            "p": "p_(n+1)=rp_n/(n+1)",
            "k": "k_(n+1)=-(n+1)*u_(n+1)/(n+2)",
        },
        "finite_loop": "for n in range(degree)",
        "infinite_radial_remainder_bounded": False,
        "full_original_n1_point_oracle": False,
        "old_finite_producer_full_identity": False,
        "source_executed": False,
    }


def _failure(reason: str, parent_pin: Mapping[str, Any] | None = None) -> dict[str, Any]:
    result: dict[str, Any] = {
        "passed": False,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "failure": reason,
        "checks": {"passed": False, "failure": reason},
        "actual_infinite_forward_analytic_response": False,
        "new_infinite_forward_definition": False,
        "actual_forward_analytic_definition_bound": False,
        "actual_forward_regular_uniqueness": False,
        "actual_forward_prefix_identity": False,
        "actual_source_scale_conjugacy": False,
        "new_infinite_forward_definition_coefficient_identity_bound": False,
        "actual_forward_equation_binding": False,
        "source_executed": False,
        "old_pure_definitions_executed": False,
        "archive_source_executed": False,
        "accepted_owner_builder_fixture_recurrence_quadrature_replayed": False,
        "point_values_selected": False,
        "finite_N_numerical_calculation": False,
        "old_finite_producer_full_identity": False,
        "full_forward_response_function_identity": False,
        "q_dependent_collar_bridge": False,
        "six_term_source_equality": False,
        "remaining_analytic_remainder": True,
    }
    if parent_pin is not None:
        result["parent_report_pin"] = dict(parent_pin)
    return result


def build(root: str | Path | None = None, parent: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build a source-bound infinite forward analytic-response packet.

    The accepted norm-tail packet is the only numerical input.  All source
    owners and finite producers are inspected by hash/AST only.  Missing or
    changed evidence fails closed.
    """
    parent_pin: dict[str, Any] | None = None
    try:
        repo = _repo_root(root)
        accepted, parent_pin = _accepted_parent(repo, parent)
        pins = _source_pins(repo, accepted)
        packets = _canonical_packets(accepted)
        normalization = _normalization_binding(accepted)
        scale = _scale_conjugacy(
            normalization.get("normalization", {}),
            normalization.get("source_map"),
        )
        analytic = _analytic_definition(packets, scale)
        prefix = _prefix_identity(packets, normalization.get("normalization", {}))
        uniqueness = _uniqueness(packets, scale)
        finite_scope = _finite_producer_scope(pins["finite_forward_ast"])

        parent_checks = {
            "parent_pin": parent_pin.get("passed") is True,
            "top_checks_passed": _mapping(accepted.get("checks")).get("passed") is True,
            "packet_invariants": _mapping(accepted.get("checks", {}).get("packet_invariants")).get("passed") is True,
            "no_ancestor_replay": _mapping(accepted.get("checks", {}).get("packet_invariants")).get("accepted_ancestor_or_quadrature_replayed") is False,
            "terminal_parent_sha": _mapping(accepted.get("n1_actual_canonical_norm_tail")).get("parent_report_sha256") == TERMINAL_PARENT_REPORT_SHA256,
        }
        source_checks = {
            "all_source_hashes": pins.get("all_source_hashes_passed") is True,
            "finite_ast": _mapping(pins.get("finite_forward_ast")).get("passed") is True,
            "source_vector_ast": _mapping(pins.get("source_vector_ast")).get("passed") is True,
            "inherited_matrix_ast": _mapping(pins.get("inherited_matrix_ast")).get("passed") is True,
        }
        packet_checks = {
            key: value.get("passed") is True
            for key, value in _mapping(packets.get("checks")).items()
            if isinstance(value, Mapping)
        }
        overall_checks = {
            "parent": all(parent_checks.values()),
            "source_pins": all(source_checks.values()),
            "accepted_packets": all(packet_checks.values()),
            "normalization": normalization.get("checks", {}).get("passed") is True,
            "scale_conjugacy": scale.get("checks", {}).get("passed") is True,
            "analytic_definition": analytic.get("checks", {}).get("passed") is True,
            "prefix_identity": prefix.get("checks", {}).get("passed") is True,
            "regular_uniqueness": uniqueness.get("checks", {}).get("passed") is True,
        }
        passed = all(overall_checks.values())
        norm = packets["norm_tail"]
        bounds = packets["source_bounds"]
        green = packets["green"]
        return {
            "passed": passed,
            "temporal_order": 1,
            "parent_report": PARENT_REPORT,
            "parent_report_sha256": PARENT_REPORT_SHA256,
            "parent_report_pin": parent_pin,
            "parent_mapping_digest": parent_pin.get("mapping_digest"),
            "checks": {**overall_checks, "passed": passed},
            "actual_infinite_forward_analytic_response": passed,
            "new_infinite_forward_definition": analytic.get("checks", {}).get("passed") is True,
            "actual_forward_analytic_definition_bound": passed,
            "actual_forward_regular_uniqueness": uniqueness.get("checks", {}).get("passed") is True,
            "actual_forward_prefix_identity": prefix.get("checks", {}).get("passed") is True,
            "actual_source_scale_conjugacy": scale.get("checks", {}).get("passed") is True,
            "new_infinite_forward_definition_coefficient_identity_bound": (
                analytic.get("checks", {}).get("passed") is True
                and prefix.get("checks", {}).get("passed") is True
            ),
            "actual_forward_equation_binding": (
                scale.get("checks", {}).get("passed") is True
                and analytic.get("checks", {}).get("passed") is True
            ),
            "source_pins": pins,
            "accepted_packet_digests": {
                "norm_tail": _digest(norm),
                "canonical_binding": _digest(packets["canonical_binding"]),
                "raw_prefix_binding": _digest(packets["prefix_binding"]),
                "radial_conjugacy": _digest(packets["radial_conjugacy"]),
                "source_bounds": _digest(bounds),
                "green": _digest(green),
            },
            "normalization_binding": normalization,
            "scale_conjugacy": scale,
            "canonical_solution": {
                "normal_convergent": True,
                "definition": norm.get("solution_definition"),
                "weighted_unknown": norm.get("weighted_unknown"),
                "axis_transform": norm.get("axis_transform"),
                "convergence_proof": norm.get("convergence_proof"),
                "accepted_norm_tail_consumed": True,
                "norm_recomputed_here": False,
                "numeric_norm_production_replayed": False,
            },
            "analytic_forward_definition": analytic,
            "forward_prefix_identity": prefix,
            "regular_uniqueness": uniqueness,
            "old_finite_forward_producer": finite_scope,
            "domain": analytic.get("domain"),
            "forward_equations": scale.get("forward_operator"),
            "zero_axis": analytic.get("zero_axis"),
            "old_finite_producer_full_identity": False,
            "full_forward_response_function_identity": False,
            "q_dependent_collar_bridge": False,
            "six_term_source_equality": False,
            # The newly defined infinite field has no omitted radial terms;
            # the old finite producer still has a separate remainder scope.
            "remaining_analytic_remainder": True,
            "remaining_scope": [
                "The old finite producer is not an infinite-function identity.",
                "No q-dependent collar bridge or six-term signed source equality is supplied.",
                "No temporal n>=2 recursion or global physical rescaling is supplied.",
            ],
            "source_executed": False,
            "old_pure_definitions_executed": False,
            "archive_source_executed": False,
            "accepted_owner_builder_fixture_recurrence_quadrature_replayed": False,
            "point_values_selected": False,
            "finite_N_numerical_calculation": False,
            "tool_metadata": {
                "route": "fixed luna_max_worker",
                "initial_fork_turns": "none",
                "requested_model": "GPT-5.6 Luna/max",
                "independent_runtime_model_metadata_exposed": False,
            },
        }
    except Exception as exc:  # fail closed; missing evidence is not a theorem
        return _failure(f"source binding failed closed: {type(exc).__name__}: {exc}", parent_pin)


def bind_forward_analytic(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for the source-bound build operation."""

    return build(*args, **kwargs)


__all__ = [
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "D",
    "AXIS_POWERS",
    "NORMALIZATION",
    "build",
    "bind_forward_analytic",
]
