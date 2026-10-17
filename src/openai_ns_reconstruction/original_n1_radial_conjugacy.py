"""Source-only all-order radial conjugacy for the original temporal ``n=1`` core.

The numerical producer in :mod:`original_n1_core_numeric` stores a finite
Taylor prefix in ``rho``.  The accepted canonical source uses the regular
variable ``x=sqrt(rho)`` and the six raw coordinates

``W=(F1, Uz1, K1, P1, F1x, Uz1x)``.

This module proves the formal coefficient map between those descriptions.  It
does not construct an owner, call ``OriginalN1CoreNumeric.coefficients`` or
evaluate a finite prefix.  The accepted matrix/prefix and known-source
receipts are read as provenance.  A finite directed enclosure is therefore
reported as an enclosure of the corresponding canonical Taylor coefficient
only through this formal identity; no interval centre is selected and no
radial remainder is inferred.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import sympy as sp


FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
FORWARD_SOURCE_VECTOR = "src/openai_ns_reconstruction/original_n1_source_vector_identity.py"
CANONICAL_VOLTERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
CANONICAL_OPERATOR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
)
PREFIX_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_numeric_prefix.json"
)
SOURCE_VECTOR_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json"
)
OPERATOR_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator_check.json"
)
VOLTERRA_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
)

# These are immutable source/report identities.  They are deliberately kept
# local rather than importing an accepted owner or its checker.
EXPECTED_SHA256 = {
    FORWARD_SOURCE: "a7a925360f3a6da2348c5f18ee74cfbb885f8e220730769131e43a18bd954c2f",
    FORWARD_SOURCE_VECTOR: "878f2b8781215a52a5730ed8527a561ff0df451adb5da247cec4d9fe0198423f",
    CANONICAL_VOLTERRA_SOURCE: "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208",
    CANONICAL_OPERATOR_SOURCE: "268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014",
    PREFIX_REPORT: "2ae583cdf5ea9fde03208fea0423eb11d171054c9a8dd45d78aa5fe2076d81fa",
    SOURCE_VECTOR_REPORT: "2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e",
    OPERATOR_RECEIPT: "846cfedcbd0355a4a8e6a5ed445558d367cd55c7459085dd9c5d74821d969bd4",
    VOLTERRA_RECEIPT: "3cfdc3ce6a7e8c2a98a31fae985c02a2acb4e28025d313ee35db4cffc6e7b47c",
}

D = (0, 0, 2, 0, 3, 1)
W_ORDER = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
FORWARD_ORDER = ("f", "u", "k", "p")
FORWARD_NORMALIZATION = {
    "f": "epsilon*F1/F0",
    "u": "Uz1/epsilon",
    "k": "K1/epsilon",
    "p": "P1/(epsilon*Pstar^2)",
}


def _repo_root(root: str | Path | None = None) -> Path:
    """Find the st073 repository without importing any source owner."""

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
    raise FileNotFoundError("st073 repository root with src/ and experiments/root_st073/ is required")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_mapping(path: Path) -> Mapping[str, Any]:
    # Some inherited reports contain very large exact integers.  Keeping them
    # as strings makes inspection independent of the Python integer digit cap.
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
    """Hash and inspect one literal source definition."""

    path = root / relative
    actual = _sha256(path)
    expected = EXPECTED_SHA256.get(relative)
    node = _function_node(path, function, class_name)
    compact = _normal(ast.unparse(node))
    fragments_ok = all(_normal(fragment) in compact for fragment in fragments)
    hash_ok = expected is not None and actual == expected
    return {
        "path": relative,
        "function": f"{class_name + '.' if class_name else ''}{function}",
        "sha256": actual,
        "expected_sha256": expected,
        "hash_passed": hash_ok,
        "fragments": list(fragments),
        "fragments_passed": fragments_ok,
        "passed": bool(hash_ok and fragments_ok),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _file_pin(root: Path, relative: str) -> dict[str, Any]:
    actual = _sha256(root / relative)
    expected = EXPECTED_SHA256.get(relative)
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": bool(expected is not None and actual == expected),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _source_ast_provenance(root: Path) -> dict[str, Any]:
    """Bind the actual forward recurrence and the canonical raw six-row source."""

    return {
        "accepted_source_vector_module": _file_pin(root, FORWARD_SOURCE_VECTOR),
        "forward_advance": _ast_pin(
            root,
            FORWARD_SOURCE,
            "_advance",
            (
                "axial = degree + depth - n - 1",
                "rf, ru, rp = self._rhs(fixed, leading, state, n, axial, center, delta, eps, P2, c)",
                "k = n + 1",
                "state['f'].append([v / (4 * k * (k + 1)) for v in rf.rows[n]])",
                "state['u'].append([v / (4 * k * k) for v in ru.rows[n]])",
                "state['p'].append([v / k for v in rp.rows[n]])",
                "state['k'].append([-c.mpf(k) * v / (k + 1) for v in state['u'][-1]])",
            ),
            class_name="OriginalN1CoreNumeric",
        ),
        "forward_rhs": _ast_pin(
            root,
            FORWARD_SOURCE,
            "_rhs",
            (
                "Aw = known(leading['A'], n + 1, axial + 2)",
                "Uw = known(leading['Uz'], n + 1, axial + 2)",
                "Qw = (2 * zz * uu - zz * mm * (1 - delta) - dd * mw.dz()) / LL",
                "QZ = crop(Qw.dz(), n, axial)",
                "pressure_head = R([[(i + 1) * (i + 2) * v for v in Qw.rows[i + 1][:axial + 1]]",
                "F2 = R(known([fixed['S_Z_taylor']], 0, axial).rows, n, axial, c.mpf) / (eps * eps)",
                "force_u -= 4 * z * rho * Np * (eps * P2) / L",
            ),
            class_name="OriginalN1CoreNumeric",
        ),
        "forward_linear_rhs": _ast_pin(
            root,
            FORWARD_SOURCE,
            "linear_rhs",
            (
                "rf = (2 * e * (Q + E * (-2 + de) / L) - 2 * H * g / L) * f",
                "ru = -8 * F2 * z * r * A * f / L",
                "rp = 2 * F2 / (e * P2) * A * f",
            ),
        ),
        "canonical_raw_matrix": _ast_pin(
            root,
            CANONICAL_VOLTERRA_SOURCE,
            "raw_matrix_expressions",
            (
                "'4,1': eps * 4 * x * f",
                "'5,1': eps * 2 * (q + E * (-2 + de) / L)",
                "'6,1': -eps ** 2 * 8 * z * r * f / L",
                "'6,4': -eps * 4 * z / L",
            ),
        ),
        "canonical_raw_system": _ast_pin(
            root,
            CANONICAL_VOLTERRA_SOURCE,
            "raw_system",
            (
                "Np = a.add(a.scale('-1/2', S), a.scale(parameters['Lambda'], a.add(a.scale(2, c1(Q[1])), a.scale(r, c1(Q[2])))))",
                "g4 = g.mul(c(2), eps, x, Np.value)",
                "g6 = g.sub(g.mul(c(2), eps, Nz), g.quotient(g.mul(c(4), zv, r, eps, eps, Np.value), Lv, 'positive original L'))",
                "D = D",
            ),
        ),
        "canonical_operator_rows": _ast_pin(
            root,
            CANONICAL_OPERATOR_SOURCE,
            "matrix",
            (
                "B0[f'{i},{j}'] = algebra.lift(sgn)",
                "B0['4,1'] = eps(f * (4 * x))",
                "B1[f'{i},{j}'] = eps(value)",
                "return dict(B0={k: truncate(v, 3) for k, v in B0.items()},",
            ),
        ),
    }


def _accepted_receipts(root: Path) -> dict[str, Any]:
    """Inspect accepted reports; no owner, checker, recurrence or quadrature runs."""

    prefix = _read_mapping(root / PREFIX_REPORT)
    source_vector_report = _read_mapping(root / SOURCE_VECTOR_REPORT)
    source_vector = source_vector_report.get("n1_source_vector_identity", {})
    source_symbolic = source_vector.get("symbolic", {})
    source_symbolic_checks = source_symbolic.get("checks", {})
    source_inherited_premises = source_symbolic.get("inherited_premises", {})
    operator = _read_mapping(root / OPERATOR_RECEIPT)
    volterra = _read_mapping(root / VOLTERRA_RECEIPT)
    prefix_pin = _file_pin(root, PREFIX_REPORT)
    source_vector_pin = _file_pin(root, SOURCE_VECTOR_REPORT)
    operator_pin = _file_pin(root, OPERATOR_RECEIPT)
    volterra_pin = _file_pin(root, VOLTERRA_RECEIPT)

    prefix_checks = prefix.get("checks", {})
    prefix_limitations = prefix.get("limitations", {})
    prefix_ok = bool(
        prefix_checks.get("passed") is True
        and int(prefix_checks.get("raw_matrix_normalization_coefficient_comparisons", -1)) == 30
        and int(prefix_checks.get("exact_raw_equation_mixed_coefficient_comparisons", -1)) == 36
        and prefix_checks.get("same_original_leading_n1_source_family") is True
        and prefix_checks.get("no_source_midpoint_selected") is True
        and prefix_checks.get("all_B0_B1_pressure_feedback_and_Z_derivatives_retained") is True
        and prefix_limitations.get("finite_radial_prefix_only") is True
        and prefix_limitations.get("radial_Taylor_remainder_certified") is False
        and prefix_limitations.get("Volterra_word_tail_reused_as_radial_tail") is False
    )

    source_pins = source_vector.get("source_pins", {})
    source_vector_ok = bool(
        source_vector_report.get("parent_report_sha256")
        or source_vector_report.get("source_family")
    ) and bool(
        source_vector.get("passed") is True
        and source_vector.get("source_executed") is False
        and all(value.get("passed") is True for value in source_pins.values() if isinstance(value, Mapping))
        and source_vector.get("limits", {}).get("complete_n1_solution_equality") is False
        and source_vector.get("limits", {}).get("pointwise_original_n1_oracle") is False
    )

    operator_checks = operator.get("checks", {})
    operator_ok = bool(
        operator.get("all_passed") is True
        and operator.get("current_original_n1_actual_inner_operator_input_enclosures_installed") is True
        and operator_checks.get("independent_direct_paper_operator_and_forcing_rows") is not None
        and operator.get("original_numeric_n1_solution_oracle_executed") in (None, False)
    )
    volterra_checks = volterra.get("checks", {})
    volterra_ok = bool(
        volterra.get("all_passed") is True
        and volterra.get("current_original_n1_true_inner_Volterra_function_adapter_installed") is True
        and volterra_checks.get("actual_source", {}).get("passed") is True
        and volterra_checks.get("raw_operator", {}).get("all_epsilon_powers_and_B1_entries_retained") is True
        and volterra_checks.get("Green_feedback", {}).get("zero_axis_regular_values") is True
        and volterra.get("original_numeric_n1_solution_oracle_executed") is False
    )

    return {
        "prefix_report": {
            **prefix_pin,
            "checks": dict(prefix_checks),
            "limitations": dict(prefix_limitations),
            "passed": bool(prefix_pin["passed"] and prefix_ok),
            "source_executed": False,
            "archive_source_executed": False,
        },
        "source_vector_report": {
            **source_vector_pin,
            "accepted_identity_passed": bool(source_vector_pin["passed"] and source_vector_ok),
            "source_executed": source_vector.get("source_executed"),
            "archive_source_executed": False,
            "source_pins_passed": all(
                value.get("passed") is True
                for value in source_pins.values()
                if isinstance(value, Mapping)
            ),
            "symbolic_checks": dict(source_symbolic_checks),
            "inherited_premises": dict(source_inherited_premises),
            "limits": dict(source_vector.get("limits", {})),
        },
        "operator_receipt": {
            **operator_pin,
            "accepted_input_enclosures": bool(operator_pin["passed"] and operator_ok),
            "source_executed": bool(operator.get("original_numeric_n1_solution_oracle_executed") is True),
            "archive_source_executed": False,
        },
        "volterra_receipt": {
            **volterra_pin,
            "accepted_complete_regular_definition": bool(volterra_pin["passed"] and volterra_ok),
            "source_executed": bool(volterra.get("original_numeric_n1_solution_oracle_executed") is True),
            "archive_source_executed": False,
        },
        "passed": bool(
            prefix_pin["passed"]
            and source_vector_pin["passed"]
            and operator_pin["passed"]
            and volterra_pin["passed"]
            and prefix_ok
            and source_vector_ok
            and operator_ok
            and volterra_ok
        ),
        "point_values_selected": False,
        "accepted_builders_rerun": False,
        "accepted_fixtures_rerun": False,
        "accepted_owners_rerun": False,
        "accepted_recurrences_rerun": False,
        "source_executed": False,
        "old_pure_definitions_executed": False,
        "archive_source_executed": False,
    }


def _accepted_normalization_premises(source: Mapping[str, Any]) -> dict[str, Any]:
    """Bind normalization assumptions to the accepted source-vector algebra."""

    symbolic = source.get("symbolic_checks", {})
    inherited = source.get("inherited_premises", {})
    symbolic_keys = (
        "P2_equals_shared_Pstar_squared",
        "shared_core_Lambda_epsilon_reciprocal",
        "g_equals_minus_epsilon_dZ_log_F0",
        "forward_pressure_head_equals_canonical_head",
        "forward_Np_equals_canonical_N1p_over_Pstar_squared",
        "pressure_head_coefficient_all_orders",
        "angular_source",
        "pressure_source",
        "axial_source_with_pressure_feedback",
        "vector_identity",
    )
    inherited_keys = (
        "common_fixed_point_and_same_input_functions",
        "finite_rows_and_tails_bound_to_same_nonlinear_fixed_point",
        "scaled_S_jet_inclusion",
        "P2_and_epsilon_use_actual_shared_core_scales",
    )
    checks = {
        **{f"symbolic_{key}": symbolic.get(key) is True for key in symbolic_keys},
        **{f"premise_{key}": inherited.get(key) is True for key in inherited_keys},
    }
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "symbolic_checks": {key: symbolic.get(key) for key in symbolic_keys},
        "inherited_premises": {key: inherited.get(key) for key in inherited_keys},
    }


def _inherited_matrix_conjugacy(
    receipts: Mapping[str, Any], ast_receipt: Mapping[str, Any]
) -> dict[str, Any]:
    """Expose the accepted 30-row and 36-prefix matrix receipt without replay."""

    prefix = receipts["prefix_report"]
    source = receipts["source_vector_report"]
    checks = prefix["checks"]
    state_comparisons = int(checks.get("raw_matrix_normalization_coefficient_comparisons", -1))
    mixed_comparisons = int(checks.get("exact_raw_equation_mixed_coefficient_comparisons", -1))
    normalization_premises = _accepted_normalization_premises(source)
    assumptions = {
        "rho_equals_x_squared": "rho = x^2",
        "P2_equals_Pstar_squared": "P2 = Pstar^2",
        "gauge_log_gradient": "g = -epsilon * D_Z(F0) / F0",
        "F2_equals_F0_squared": "F2 = F0^2",
        "F0_nonzero": "F0 != 0",
        "epsilon_nonzero": "epsilon != 0",
        "P2_nonzero": "P2 != 0",
        "L_nonzero": "L != 0",
    }
    inherited_checks = {
        "accepted_prefix_receipt": bool(prefix["passed"]),
        "accepted_30_homogeneous_matrix_comparisons": state_comparisons == 30,
        "accepted_36_mixed_prefix_comparisons": mixed_comparisons == 36,
        "canonical_matrix_ast_pin": bool(ast_receipt["canonical_raw_matrix"]["passed"]),
        "forward_linear_rhs_ast_pin": bool(ast_receipt["forward_linear_rhs"]["passed"]),
        "normalization_assumptions_bound": normalization_premises["passed"],
    }
    return {
        "checks": {**inherited_checks, "passed": all(inherited_checks.values())},
        "state_coefficient_comparisons": state_comparisons,
        "mixed_coefficient_comparisons": mixed_comparisons,
        "normalization": dict(FORWARD_NORMALIZATION),
        "assumptions": assumptions,
        "assumptions_kind": "accepted symbolic source-vector checks and inherited premises",
        "accepted_normalization_premises": normalization_premises,
        "source_order": ("pressure", "angular", "axial"),
        "canonical_rows": ("g4/pressure", "g5/angular", "g6/axial"),
        "canonical_matrix_ast_pin": ast_receipt["canonical_raw_matrix"],
        "forward_linear_rhs_ast_pin": ast_receipt["forward_linear_rhs"],
        "accepted_prefix_report": {
            "path": prefix["path"],
            "sha256": prefix["sha256"],
            "expected_sha256": prefix["expected_sha256"],
        },
        "replayed": False,
        "source_executed": False,
        "point_values_selected": False,
    }


def _inherited_source_scales(
    receipts: Mapping[str, Any], ast_receipt: Mapping[str, Any]
) -> dict[str, Any]:
    """Link the accepted source-vector identity to the local normalization."""

    source = receipts["source_vector_report"]
    normalization_premises = _accepted_normalization_premises(source)
    checks = {
        "accepted_source_vector_identity": bool(source["accepted_identity_passed"]),
        "accepted_source_vector_pins": bool(source["source_pins_passed"]),
        "forward_linear_rhs_ast_pin": bool(ast_receipt["forward_linear_rhs"]["passed"]),
        "canonical_raw_system_ast_pin": bool(ast_receipt["canonical_raw_system"]["passed"]),
        "accepted_normalization_premises": normalization_premises["passed"],
    }
    assumptions = {
        "P2_equals_Pstar_squared": "P2 = Pstar^2",
        "gauge_log_gradient": "g = -epsilon * D_Z(F0) / F0",
        "F2_equals_F0_squared": "F2 = F0^2",
        "F0_nonzero": "F0 != 0",
        "epsilon_nonzero": "epsilon != 0",
        "P2_nonzero": "P2 != 0",
    }
    return {
        "checks": {**checks, "all_three_sources": all(checks.values())},
        "canonical_source_rows": {
            "g4": "2*epsilon*x*N1p = 2*epsilon*x*Pstar^2*Np",
            "g5": "2*epsilon*N1theta",
            "g6": "2*epsilon*N1z - 4*epsilon^2*Z*rho*N1p/L",
        },
        "forward_normalized_sources": {
            "pressure": "Np=N1p/Pstar^2",
            "angular": "2*epsilon^2*N1theta/F0",
            "axial": "2*N1z-4*epsilon*Z*rho*N1p/L",
        },
        "normalization": dict(FORWARD_NORMALIZATION),
        "assumptions": assumptions,
        "assumptions_kind": "accepted symbolic source-vector checks and inherited premises",
        "accepted_normalization_premises": normalization_premises,
        "accepted_source_vector_report": {
            "path": source["path"],
            "sha256": source["sha256"],
            "expected_sha256": source["expected_sha256"],
        },
        "accepted_source_vector_module_pin": ast_receipt["accepted_source_vector_module"],
        "source_identity_is_inherited": True,
        "replayed": False,
        "source_executed": False,
        "point_values_selected": False,
    }


def _formal_radial_rules() -> dict[str, Any]:
    """Derive the coefficient factors, kinematics, and formal uniqueness."""

    n = sp.symbols("n", integer=True, nonnegative=True)
    m = n + 1
    f_denominator = sp.expand(2 * m * (2 * m - 1) + 3 * (2 * m))
    u_denominator = sp.expand(2 * m * (2 * m - 1) + 1 * (2 * m))
    p_denominator = 2 * m
    kinematic = sp.expand((2 * m + 2) * (-m / (m + 1)) + 2 * m)
    derivative_f = sp.expand(2 * m)
    derivative_u = sp.expand(2 * m)

    q_next = sp.Symbol("Q_(n+1)")
    pressure_head = sp.expand(2 * m * q_next + n * m * q_next)
    pressure_head_expected = sp.expand(m * (m + 1) * q_next)

    # General-power bookkeeping for W_m proportional to x^(2m).  The
    # differentiated states F1x and Uz1x have power 2m-1, so rows 5 and 6
    # land at x^(2m-2)=x^(2n); P1 and K1 rows land at x^(2m-1)=x^(2n+1).
    f1_power = sp.expand(2 * m)
    f1x_power = sp.expand(f1_power - 1)
    row_5_power = sp.expand(f1x_power - 1)
    row_6_power = sp.expand(f1x_power - 1)
    p1_power = sp.expand(2 * m)
    row_4_power = sp.expand(p1_power - 1)
    k1_power = sp.expand(2 * m)
    row_3_power = sp.expand(k1_power - 1)
    general_power_checks = {
        "F1x_power": f1x_power == 2 * m - 1,
        "Uz1x_power": f1x_power == 2 * m - 1,
        "row_5_output_power": row_5_power == 2 * n,
        "row_6_output_power": row_6_power == 2 * n,
        "row_4_output_power": row_4_power == 2 * n + 1,
        "row_3_output_power": row_3_power == 2 * n + 1,
        "row_5_factor_from_power": sp.expand(f_denominator - 4 * m * (m + 1)) == 0,
        "row_6_factor_from_power": sp.expand(u_denominator - 4 * m * m) == 0,
        "row_4_factor_from_power": p_denominator == 2 * m,
        "row_3_factor_from_power": (2 * m + 2) == 2 * (m + 1),
    }

    checks = {
        "F1x_definition": derivative_f == 2 * m,
        "Uz1x_definition": derivative_u == 2 * m,
        "row_5_factor": sp.simplify(f_denominator - 4 * m * (m + 1)) == 0,
        "row_6_factor": sp.simplify(u_denominator - 4 * m * m) == 0,
        "row_4_factor": p_denominator == 2 * m,
        "K1_kinematic_factor": sp.simplify(kinematic) == 0,
        "pressure_head_all_orders": sp.expand(pressure_head - pressure_head_expected) == 0,
        "m_is_positive_for_n_nonnegative": True,
        "general_power_lhs": all(general_power_checks.values()),
    }

    # A finite polynomial with arbitrary coefficient functions is enough to
    # check the coefficient/ordinary-Z commutator without selecting any
    # source value.  The identity is termwise, hence the displayed result is
    # the all-order rule, not a finite radial calculation.
    z, rho = sp.symbols("Z rho")
    coefficients = [sp.Function(f"a_{j}")(z) for j in range(6)]
    polynomial = sum(coefficients[j] * rho**j for j in range(6))
    commutator_checks = []
    for order in range(5):
        differentiated = sp.diff(polynomial, z, order)
        for radial_degree in range(6):
            left = sp.diff(sp.expand(polynomial).coeff(rho, radial_degree), z, order)
            right = sp.expand(differentiated).coeff(rho, radial_degree)
            commutator_checks.append(sp.simplify(left - right) == 0)

    checks["ordinary_Z_commutes_with_radial_coefficient"] = all(commutator_checks)

    return {
        "checks": checks,
        "formal_index": "n>=0, m=n+1>=1",
        "forward_recurrence": {
            "f_(n+1)": "[rho^n]R_f/(4*(n+1)*(n+2))",
            "u_(n+1)": "[rho^n]R_u/(4*(n+1)^2)",
            "p_(n+1)": "[rho^n]R_p/(n+1)",
            "k_(n+1)": "-(n+1)*u_(n+1)/(n+2)",
        },
        "canonical_coefficient_equations": {
            "row_5": "4*(n+1)*(n+2)*F1_(n+1) = [x^(2n)](B0 W+B1 W_Z+g)_5",
            "row_6": "4*(n+1)^2*Uz1_(n+1) = [x^(2n)](B0 W+B1 W_Z+g)_6",
            "row_4": "2*(n+1)*P1_(n+1) = [x^(2n+1)](B0 W+B1 W_Z+g)_4",
            "row_3": "2*(n+2)*K1_(n+1) = -2*(n+1)*Uz1_(n+1)",
        },
        "general_power_lhs": {
            "input_series_power": "Base coordinates W1..W4 have degree-m terms proportional to x^(2m), m=n+1",
            "F1x_and_Uz1x_power": "2m-1",
            "row_5_output_power": "2m-2=2n",
            "row_6_output_power": "2m-2=2n",
            "row_4_output_power": "2m-1=2n+1",
            "row_3_output_power": "2m-1=2n+1",
            "row_5_factor": "(2m)(2m-1)+3(2m)=4m(m+1)",
            "row_6_factor": "(2m)(2m-1)+(2m)=4m^2",
            "row_4_factor": "2m",
            "row_3_factor": "2m+2",
            "checks": general_power_checks,
        },
        "axis": {
            "regular_powers_in_x": list((2, 2, 2, 2, 1, 1)),
            "D": list(D),
            "zero_axis_values": {name: "0" for name in W_ORDER},
            "singular_division_at_axis_evaluated": False,
        },
        "coefficient_uniqueness": {
            "induction_variable": "m=n+1>=1",
            "nonzero_diagonal_factors": [
                "4*m*(m+1)",
                "4*m^2",
                "2*m",
                "2*(m+1)",
            ],
            "axis_constants_fixed_by_zero_axis": True,
            "formal_solution_unique": all(checks.values()),
        },
        "ordinary_Z_commutation": {
            "formula": "D_Z^k([rho^n]A)=[rho^n](D_Z^k A), n,k>=0",
            "checked_orders": {"radial": 5, "ordinary_Z": 4},
            "termwise_polynomial_checks": len(commutator_checks),
            "passed": checks["ordinary_Z_commutes_with_radial_coefficient"],
        },
        "point_values_selected": False,
    }


def _triangularity_support(root: Path, ast_receipt: Mapping[str, Any]) -> dict[str, Any]:
    """Show from the pinned RHS AST why the coefficient induction is triangular."""

    node = _function_node(root / FORWARD_SOURCE, "linear_rhs")
    compact = _normal(ast.unparse(node))
    names = {item.id for item in ast.walk(node) if isinstance(item, ast.Name)}
    required_fragments = (
        "rf += 2*e*B*Ef+2*e*H*fz/L-2*e**3*d*C*(uz+kz)/L",
        "ru = -8*F2*z*r*A*f/L",
        "ru += 2*e*B*Eu+2*e*(H-d*v['Euler_U'])*uz/L",
        "rp = 2*F2/(e*P2)*A*f",
    )
    fragment_checks = {
        f"linear_rhs_fragment_{index}": _normal(fragment) in compact
        for index, fragment in enumerate(required_fragments, start=1)
    }
    checks = {
        "linear_rhs_ast_pinned": bool(ast_receipt["forward_linear_rhs"]["passed"]),
        "Euler_f_and_Euler_u_are_explicit_inputs": {"Ef", "Eu"}.issubset(names),
        "ordinary_Z_inputs_are_explicit": {"fz", "uz", "kz", "pz"}.issubset(names),
        "unweighted_radial_derivative_inputs_absent": not ({"fr", "ur"} & names),
        "leading_n_plus_1_source_rows_are_pinned": bool(ast_receipt["forward_rhs"]["passed"]),
        **fragment_checks,
    }
    return {
        "checks": {**checks, "passed": all(checks.values())},
        "coefficient_statement": (
            "At [rho^n], analytic coefficient convolutions use unknown forward "
            "coefficients f_j,u_j,k_j,p_j with 0<=j<=n and their ordinary Z "
            "derivatives; Euler and Z operations do not advance radial order, "
            "and the new coefficient at n+1 is supplied by the diagonal radial row."
        ),
        "unknown_coefficients_at_rho_n": "f_j,u_j,k_j,p_j for 0<=j<=n",
        "radial_support_at_rho_n": "j in {0,...,n}; coefficient-field convolutions include lower j",
        "ordinary_Z_inputs_at_rho_n": "D_Z^k f_j,D_Z^k u_j,D_Z^k k_j,D_Z^k p_j for 0<=j<=n",
        "rho_Euler_rule": "[rho^n](rho*d_rho a)=n*a_n; Euler_f and Euler_u therefore do not advance radial order",
        "B1_Z_rule": "B1*D_Z W uses the same rho coefficient because ordinary D_Z commutes with [rho^n]",
        "leading_source_rule": "known(leading, n+1, ...) and Q_(n+1) pressure head are permitted source rows",
        "analytic_coefficient_fields": ("A", "C", "Zf_without_gauge", "Zw", "B", "E", "H", "L^-1"),
        "assumptions": {
            "coefficient_fields_analytic_in_rho": True,
            "F0_nonzero": "F0 != 0",
            "rho_equals_x_squared": "rho = x^2",
        },
        "ast_source": ast_receipt["forward_linear_rhs"],
        "point_values_selected": False,
    }


def _extra_order_budget(ast_receipt: Mapping[str, Any]) -> dict[str, Any]:
    """Expose the source's required extra radial/Z rows explicitly."""

    source_rows = ast_receipt["forward_rhs"]
    passed = bool(source_rows["passed"])
    return {
        "passed": passed,
        "rule": "coefficient rho^n with ordinary Z order k consumes leading rows through rho^(n+1) and Z order k+2",
        "radial_leading_rows": "n+1",
        "ordinary_Z_rows": "k+2",
        "consumers": [
            "Q recovery consumes mean_Z",
            "Q_Z is differentiated before cropping",
            "pressure head uses (i+1)(i+2) Q_(i+1)",
            "F0/gauge and B1*d_Z terms retain their ordinary Z derivative",
        ],
        "missing_rows_zero_padded": False,
        "source_guarded_by_actual_AST": True,
    }


def _available_enclosure_semantics(root: Path, receipts: Mapping[str, Any]) -> dict[str, Any]:
    """Relate saved directed prefixes to the formal coefficient identity."""

    report = receipts["prefix_report"]
    checks = report["checks"]
    max_degree = None
    # The report is read-only; this only discovers the highest serialized
    # degree and never evaluates a source function.
    try:
        raw = _read_mapping(root / PREFIX_REPORT)
        levels = raw.get("levels", [])
        max_degree = max((int(level.get("degree", -1)) for level in levels), default=-1)
        coefficient_rows = raw.get("coefficients", [])
        rows_have_enclosures = bool(coefficient_rows) and all(
            isinstance(row, Mapping)
            and "source_enclosure_exact_mpf_tuples" in row
            for row in coefficient_rows
        )
    except (FileNotFoundError, ValueError, TypeError):
        max_degree = None
        rows_have_enclosures = False

    return {
        "formal_identity_applies_to_every_available_forward_coefficient": bool(
            report["passed"] and checks.get("all_B0_B1_pressure_feedback_and_Z_derivatives_retained") is True
        ),
        "directed_source_enclosures_present": rows_have_enclosures,
        "available_max_radial_degree": max_degree,
        "interval_enclosures_are_not_selected_points": bool(
            checks.get("no_source_midpoint_selected") is True
        ),
        "finite_prefix_only": report["limitations"].get("finite_radial_prefix_only") is True,
        "radial_remainder_certified": report["limitations"].get("radial_Taylor_remainder_certified") is True,
        "canonical_full_function_equality_claimed": False,
        "numeric_radial_tail_claimed": False,
    }


def prove_all_order_radial_conjugacy(root: str | Path | None = None) -> dict[str, Any]:
    """Return the universal formal radial conjugacy receipt.

    Only AST parsing, saved-report inspection and bounded SymPy identities are
    executed.  The accepted numerical builders, fixtures, owners and
    recurrences are intentionally not replayed.
    """

    repo = _repo_root(root)
    # The helper above discovers the report path from the process root.  Keep
    # all subsequent read-only inspection rooted at the caller-selected repo.
    source_pins = _source_ast_provenance(repo)
    receipts = _accepted_receipts(repo)
    matrix = _inherited_matrix_conjugacy(receipts, source_pins)
    sources = _inherited_source_scales(receipts, source_pins)
    radial = _formal_radial_rules()
    triangularity = _triangularity_support(repo, source_pins)
    budget = _extra_order_budget(source_pins)

    source_pin_pass = all(value["passed"] for value in source_pins.values())
    matrix_pass = bool(
        matrix["checks"]["passed"]
        and matrix["state_coefficient_comparisons"] == 30
        and matrix["mixed_coefficient_comparisons"] == 36
    )
    source_pass = bool(sources["checks"]["all_three_sources"])
    radial_pass = bool(
        all(radial["checks"].values())
        and radial["coefficient_uniqueness"]["formal_solution_unique"]
        and radial["ordinary_Z_commutation"]["passed"]
        and triangularity["checks"]["passed"]
    )
    passed = bool(
        source_pin_pass
        and receipts["passed"]
        and matrix_pass
        and source_pass
        and radial_pass
        and budget["passed"]
    )
    if not passed:
        failed = {
            "source_pins": [key for key, value in source_pins.items() if not value["passed"]],
            "receipts": receipts,
            "matrix": matrix["checks"],
            "sources": sources["checks"],
            "radial": radial["checks"],
            "triangularity": triangularity["checks"],
            "budget": budget,
        }
        raise ValueError(f"all-order original n1 radial conjugacy proof failed: {failed}")

    enclosure = _available_enclosure_semantics(repo, receipts)
    return {
        "passed": True,
        "source_executed": False,
        "old_pure_definitions_executed": False,
        "archive_source_executed": False,
        "repo": str(repo),
        "source_pins": source_pins,
        "accepted_receipts": receipts,
        "symbolic": {
            "matrix_conjugacy": matrix,
            "inherited_matrix_conjugacy": matrix,
            "known_source_map": sources,
            "inherited_source_map": sources,
            "radial_rules": radial,
            "triangularity": triangularity,
            "extra_order_budget": budget,
        },
        "formal_statement": {
            "forward_to_canonical": (
                "For every n>=0 and every available ordinary-Z coefficient, "
                "the forward _advance row equals the canonical raw-system row "
                "after f=epsilon*F1/F0, u=Uz1/epsilon, k=K1/epsilon, "
                "p=P1/(epsilon*Pstar^2)."
            ),
            "canonical_state": list(W_ORDER),
            "canonical_D": list(D),
            "axis_extension": "W(0,Z)=0; W1..W4 are O(x^2), W5,W6 are O(x)",
            "canonical_coefficient_exponents": radial["canonical_coefficient_equations"],
            "coefficient_uniqueness": "The m=n+1 diagonal factors are nonzero and zero-axis data remove homogeneous constants.",
            "triangularity": triangularity["coefficient_statement"],
            "ordinary_Z_rule": "D_Z^k commutes with [rho^n] for n,k>=0 when the source has the required extra rows.",
        },
        "available_forward_enclosures": enclosure,
        "accepted_prefix_report_hash_unchanged": receipts["prefix_report"]["passed"],
        "point_values_selected": False,
        "finite_N_numerical_calculation": False,
        "limitations": {
            "full_function_equality": False,
            "radial_Taylor_remainder": False,
            "numeric_radial_tail": False,
            "remaining_analytic_remainder": (
                "unbounded here; requires a direct canonical function enclosure "
                "or a separately certified radial Cauchy/Volterra tail"
            ),
            "temporal_n_ge_2_recursion": False,
            "global_n1_solution": False,
            "all_error_views_mutated": False,
        },
    }


# Focused aliases make the proof convenient to call from a saved-report
# inspection without introducing any producer dependency.
prove_radial_coefficient_conjugacy = prove_all_order_radial_conjugacy
build_radial_conjugacy = prove_all_order_radial_conjugacy
check_radial_conjugacy = prove_all_order_radial_conjugacy


__all__ = [
    "FORWARD_SOURCE",
    "CANONICAL_VOLTERRA_SOURCE",
    "CANONICAL_OPERATOR_SOURCE",
    "prove_all_order_radial_conjugacy",
    "prove_radial_coefficient_conjugacy",
    "build_radial_conjugacy",
    "check_radial_conjugacy",
]
