"""Source-only compact terminal algebra and endpoint-support audit.

The producer sources are read as text and parsed with :mod:`ast`; they are
never imported or evaluated here.  The algebra retains all five independent
target defects and nonlinear terms with exact :class:`fractions.Fraction`
coefficients.  Any terminal identity remains conditional on the same source
branch, inverse, endpoint integral identity, and raw bump semantics.
"""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping

COMPONENT = "compact_terminal_algebra_v1"
SCHEMA = COMPONENT
COMPACT_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py"
)
CONTROL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_Rc_functional_controls.py"
)
COMPACT_SOURCE_SHA256 = "f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b"
CONTROL_SOURCE_SHA256 = "cb6f11d4a62d59a45f05d0961236d5c6236bd544dca8b18499f68b2d8d98d009"
CENTERS = (Fraction(1, 5), Fraction(1, 2), Fraction(4, 5))
RAW_BETA_DEFINITION = "exp(-1/(1-t^2)) on abs(t)<1, zero otherwise"


class SourceContractError(ValueError):
    """The pinned source no longer has the admitted AST shape."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceContractError(message)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_pinned(root: Path, relative: str, expected: str) -> tuple[str, str]:
    path = root / relative
    _require(path.is_file(), f"missing pinned source: {relative}")
    actual = _sha256(path)
    _require(actual == expected, f"changed pinned source: {relative}")
    text = path.read_text(encoding="utf-8")
    _require(bool(text.strip()), f"empty pinned source: {relative}")
    return text, actual


def _dump(node: ast.AST) -> str:
    return ast.dump(node, include_attributes=False)


def _expected_nodes(text: str) -> list[ast.AST]:
    parsed = ast.parse(text, mode="exec")
    _require(parsed.body, "expected source fragment must contain a statement")
    return list(parsed.body)


def _function(tree: ast.AST, name: str, owner: str | None = None) -> ast.AST:
    found: list[ast.AST] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) or node.name != name:
            continue
        if owner is None or getattr(node, "_compact_owner", None) == owner:
            found.append(node)
    _require(len(found) == 1, f"ambiguous or missing source function {owner + '.' if owner else ''}{name}")
    return found[0]


def _annotate_class_methods(tree: ast.Module) -> None:
    for cls in (node for node in ast.walk(tree) if isinstance(node, ast.ClassDef)):
        for child in cls.body:
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                child._compact_owner = cls.name  # type: ignore[attr-defined]


def _parse(text: str, filename: str) -> ast.Module:
    try:
        tree = ast.parse(text, filename=filename, mode="exec")
    except SyntaxError as exc:  # pragma: no cover - pinned source is expected valid
        raise SourceContractError(f"source AST failed for {filename}: {exc.msg}") from exc
    _annotate_class_methods(tree)
    return tree


def _contains_statement(node: ast.AST, source: str, label: str, *, count: int = 1) -> ast.AST:
    expected = _expected_nodes(source)
    if len(expected) == 1:
        matches = [item for item in ast.walk(node) if _dump(item) == _dump(expected[0])]
    else:
        matches = []
        for item in ast.walk(node):
            body = getattr(item, "body", None)
            if not isinstance(body, list) or len(body) < len(expected):
                continue
            for start in range(len(body) - len(expected) + 1):
                if all(_dump(body[start + offset]) == _dump(expected[offset]) for offset in range(len(expected))):
                    matches.append(body[start])
    _require(len(matches) == count, f"{label}: expected {count} exact AST occurrence(s), found {len(matches)}")
    return matches[0]


def _require_args(function: ast.AST, expected: Iterable[str], label: str) -> None:
    args = getattr(function, "args", None)
    actual = [item.arg for item in getattr(args, "args", ())]
    _require(actual == list(expected), f"{label}: argument contract changed")


def _check_normalized_targets(tree: ast.Module) -> dict[str, bool]:
    fn = _function(tree, "normalized_targets")
    _require_args(fn, ("g", "correction", "A", "mu", "N"), "normalized_targets")
    expected = (
        ("trim_amplitude", "A=a.trim(A);AA=a.mul(A,A)"),
        ("scale_correction", "scaled={k:a.scale(v,N) for k,v in correction.items()}"),
        ("five_target_rows", "d={row:a.div(scaled[k],A if degree==1 else AA,'original positive Rc amplitude') for row,k,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2))}"),
        ("joint_target", "joint=a.sub(scaled['k'],a.mul(A,scaled['m']))"),
        ("joint_row_normalization", "d[compact.ROWS[1]]=a.div(joint,a.scale(AA,mu),'original positive mu times Rc amplitude squared')"),
        ("return_target_order", "return {key:d[key] for key in compact.ROWS},scaled,joint"),
    )
    for label, fragment in expected:
        _contains_statement(fn, fragment, f"normalized_targets {label}")
    return {label: True for label, _ in expected}


def _check_nonlinear(tree: ast.Module) -> dict[str, bool]:
    fn = _function(tree, "nonlinear")
    _require_args(fn, ("a", "W", "h"), "nonlinear")
    expected = (
        "zero=a.fixed(0);cross=a.add(a.scale(a.mul(h[0],h[2]),W['cross'][0]),a.scale(a.mul(h[1],h[4]),W['cross'][1]))",
        "axial=a.add(a.scale(a.mul(h[0],h[0]),W['energy'][0]),a.scale(a.mul(h[1],h[1]),W['energy'][2]))",
        "swirl=a.add(*(a.scale(a.mul(h[i+2],h[i+2]),W['energy'][i]) for i in range(3)))",
        "pressure=a.add(*(a.scale(a.mul(h[i+2],h[i+2]),W['pressure'][i]) for i in range(3)))",
        "return [zero,a.scale(cross,a.g.quotient(a.g.one,W['mu'],'same original positive mu')),zero,a.sub(axial,a.scale(swirl,'1/2')),a.scale(pressure,'1/2')]",
    )
    for index, fragment in enumerate(expected):
        _contains_statement(fn, fragment, f"nonlinear statement {index}")
    return {f"statement_{index}": True for index in range(len(expected))}


def _check_implicit_rows(tree: ast.Module) -> dict[str, bool]:
    fn = _function(tree, "implicit_rows")
    _require_args(fn, ("g", "W", "N", "d", "values", "order"), "implicit_rows")
    expected = (
        "invN=g.quotient(g.one,N,'same original exact positive integer N')",
        "columns=[compact.quadratic_action(g,W,values,[g.one if i==j else g.zero for i in range(5)]) for j in range(5)]",
        "J=[[g.add(W['B'][i][j],g.mul(invN,columns[j][i])) for j in range(5)] for i in range(5)]",
        "inverse,det=leading.inverse_matrix(g,J)",
        "out=[[v] for v in values];recurrence=[]",
        "rhs=[g.add(rows(d[i])[n],g.mul(invN,rows(Q[i])[n])) for i in range(5)]",
        "solved=[g.neg(g.add(*(g.mul(inverse[i][j],rhs[j]) for j in range(5)))) for i in range(5)]",
        "for i in range(5):out[i].append(solved[i])",
        "return [FunctionJet(tuple(q)) for q in out],dict(Jacobian=J,inverse=inverse,determinant=det,recurrence=recurrence)",
    )
    for index, fragment in enumerate(expected):
        _contains_statement(fn, fragment, f"implicit_rows statement {index}")
    return {f"statement_{index}": True for index in range(len(expected))}


def _check_terminal_constructor(tree: ast.Module) -> dict[str, bool]:
    fn = _function(tree, "__init__", owner="CurrentOriginalTrueCompactC5Functions")
    expected = (
        ("generic_five_controls", "generic=[g.symbol('true_compact_branch_h%d'%i) for i in range(5)]"),
        ("nonlinear_Q", "Q=compact.quadratic(g,self.W,generic)"),
        ("target_plus_Q_over_N", "vec=[g.add(self.targets[key].value,g.quotient(q,self.N,'same original exact positive integer N')) for key,q in zip(compact.ROWS,Q)]"),
        ("negative_B_inverse_action", "control_map=[g.neg(v) for v in compact.inverse_action(g,self.W,vec)]"),
        ("same_branch_system", "self.system=g.node('original_compact_unique_Banach_limit_system',variables=[q.node for q in generic],exact_map=[q.node for q in control_map],actual_target_values=[q.value.node for q in self.targets.values()],actual_selected_N=self.N.node,original_weights=encode(self.W),initial_vector=[g.zero.node]*5,source_family=self.family,original_branch_transfer=self.branch_transfer,definition='same original unique C1 Banach limit; no finite iterate or enclosure chosen as its value')"),
        ("same_branch_h", "self.h,self.implicit=implicit_rows(g,self.W,self.N,list(self.targets.values()),values)"),
        ("residual_definition", "self.residual={key:a.add(self.targets[key],a.add(*(a.scale(q,b) for q,b in zip(self.h,row))),a.scale(qq,invN)) for key,row,qq in zip(compact.ROWS,self.W['B'],Q)}"),
        (
            "terminal_rhs_rows",
            "for key,row,amplitude in (('m','M',a.trim(self.A)),('h','I',a.trim(self.A)),('e','S',AA),('p','Cp',AA)):\n    decay=g.unary('exp',g.mul(g.constant(-RATES[key]),g.unary('log',g.constant(2))))\n    self.terminal_rhs[key]=a.scale(a.mul(amplitude,self.residual[row]),g.mul(decay,invN))",
        ),
        ("terminal_rhs_decay", "decay=g.unary('exp',g.mul(g.constant(-RATES[key]),g.unary('log',g.constant(2))))"),
        ("terminal_rhs_formula", "self.terminal_rhs[key]=a.scale(a.mul(amplitude,self.residual[row]),g.mul(decay,invN))"),
        ("terminal_joint", "joint=a.add(self.residual['M'],a.scale(self.residual[compact.ROWS[1]],self.mu))"),
        ("terminal_k_decay", "decay=g.unary('exp',g.mul(g.constant('-3/2'),g.unary('log',g.constant(2))))"),
        ("terminal_k_rhs", "self.terminal_rhs['k']=a.scale(a.mul(AA,joint),g.mul(decay,invN))"),
    )
    for label, fragment in expected:
        _contains_statement(fn, fragment, f"compact terminal {label}")
    return {label: True for label, _ in expected}


def _check_bump_support(compact_tree: ast.Module, controls_tree: ast.Module) -> dict[str, bool]:
    bump = _function(compact_tree, "bump_rows")
    profile = _function(compact_tree, "profiles")
    weights = _function(controls_tree, "exact_weights")
    _require_args(bump, ("g", "W", "x"), "bump_rows")
    _require_args(profile, ("g", "W", "N", "A", "h", "x"), "profiles")
    _require_args(weights, ("g", "mu"), "exact_weights")
    expected = (
        ("bump_log_coordinate", bump, "logx=g.unary('log',x)"),
        ("bump_inverse_ell", bump, "invell=g.quotient(g.one,W['ell'],'original positive log bump width')"),
        ("bump_argument", bump, "argument=g.mul(g.sub(logx,center),invell)"),
        ("raw_beta_call", bump, "raw=leading.raw_beta(g,argument)"),
        ("profile_F", profile, "F=a.scale(a.add(*(a.scale(h[i+2],bumps[i][j]) for i in range(3))),invN)"),
        ("profile_delta_E", profile, "dE=a.mul(a.trim(A),F);dV=a.mul(a.trim(A),G)"),
        ("profile_radial_delta_E", profile, "radial.append(dict(F=F,G=G,original_E=E,delta_E=dE,delta_V=dV,E=a.add(E,dE),V=dV))"),
        ("controls_log_length", weights, "L=g.unary('log',g.constant(2));ell=g.mul(L,g.constant('1/40'))"),
        ("controls_centers", weights, "centers=[g.mul(L,g.constant(q)) for q in ('1/5','1/2','4/5')]"),
        ("controls_raw_beta_definition", weights, "beta=g.node('compact_raw_beta',argument=t.node,definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',defining_module=PREFIX+'outer_pulse_map.py',defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))"),
    )
    for label, node, fragment in expected:
        _contains_statement(node, fragment, f"endpoint support {label}")
    return {label: True for label, _, _ in expected}


def _monomial(*atoms: str) -> tuple[str, ...]:
    return tuple(sorted(atoms))


def _poly_add(*values: Mapping[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    out: dict[tuple[str, ...], Fraction] = {}
    for value in values:
        for monomial, coefficient in value.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
    return {key: value for key, value in out.items() if value}


def _poly_scale(value: Mapping[tuple[str, ...], Fraction], factor: Fraction) -> dict[tuple[str, ...], Fraction]:
    return {key: factor * value for key, value in value.items() if factor * value}


def _poly_mul(left: Mapping[tuple[str, ...], Fraction], right: Mapping[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    out: dict[tuple[str, ...], Fraction] = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            key = _monomial(*(left_monomial + right_monomial))
            out[key] = out.get(key, Fraction(0)) + left_coefficient * right_coefficient
    return {key: value for key, value in out.items() if value}


def _atom(name: str) -> dict[tuple[str, ...], Fraction]:
    return {(_monomial(name)): Fraction(1)}


def _poly_serializable(value: Mapping[tuple[str, ...], Fraction]) -> list[dict[str, Any]]:
    return [
        {"monomial": list(monomial), "coefficient": f"{coefficient.numerator}/{coefficient.denominator}"}
        for monomial, coefficient in sorted(value.items())
    ]


def _fraction_terminal_algebra() -> dict[str, Any]:
    """Prove d + B h + Q/N = 0 conditional on B B_inv = I."""
    n_inv = _atom("N_inv")
    d = [_atom(f"d_{index}") for index in range(5)]
    q = [_atom(f"Q_{index}") for index in range(5)]
    h: list[dict[tuple[str, ...], Fraction]] = []
    for row in range(5):
        terms: list[dict[tuple[str, ...], Fraction]] = []
        for column in range(5):
            b_inv = _atom(f"B_inv_{row}_{column}")
            terms.append(_poly_scale(_poly_mul(b_inv, d[column]), Fraction(-1)))
            terms.append(_poly_scale(_poly_mul(_poly_mul(b_inv, n_inv), q[column]), Fraction(-1)))
        h.append(_poly_add(*terms))
    residual: list[dict[tuple[str, ...], Fraction]] = []
    for row in range(5):
        terms = [d[row], _poly_mul(n_inv, q[row])]
        for column in range(5):
            terms.append(_poly_mul(_atom(f"B_{row}_{column}"), h[column]))
        residual.append(_poly_add(*terms))
    independent = [f"d_{index}" for index in range(5)] + [f"Q_{index}" for index in range(5)]
    independent_presence = {
        name: any(any(name in monomial for monomial in row) for row in h)
        for name in independent
    }
    _require(all(independent_presence.values()), "Fraction algebra dropped an independent d_i or Q_i")
    inverse_relation_checks: dict[str, bool] = {}
    coefficient_witnesses: dict[str, dict[str, Any]] = {}
    for row in range(5):
        for column in range(5):
            key = f"{row},{column}"
            direct = Fraction(1) if row == column else Fraction(0)
            summands = [
                {"coefficient": "-1", "term": f"B_{row}_{inner}*B_inv_{inner}_{column}"}
                for inner in range(5)
            ]
            coefficient_witnesses[f"d_{key}"] = {
                "direct_fraction": f"{direct.numerator}/{direct.denominator}",
                "inverse_product": summands,
                "conditional_reduction": "delta[row,column] - sum_inner B[row,inner] B_inv[inner,column] = 0",
            }
            coefficient_witnesses[f"Q_{key}"] = {
                "direct_fraction": f"{direct.numerator}/{direct.denominator}",
                "inverse_product": summands,
                "common_factor": "N_inv",
                "conditional_reduction": "N_inv*(delta[row,column] - sum_inner B[row,inner] B_inv[inner,column]) = 0",
            }
            expected_d = {(): direct} if direct else {}
            for inner in range(5):
                expected_d[_monomial(f"B_{row}_{inner}", f"B_inv_{inner}_{column}")] = Fraction(-1)
            for name, extra in ((f"d_{column}", ()), (f"Q_{column}", ("N_inv",))):
                actual = {}
                for monomial, value in residual[row].items():
                    if name not in monomial:
                        continue
                    parts = list(monomial)
                    parts.remove(name)
                    for atom in extra:
                        _require(atom in parts, "Q term lost N_inv")
                        parts.remove(atom)
                    actual[tuple(parts)] = value
                inverse_relation_checks[f"{'d' if name.startswith('d_') else 'Q'}_{key}"] = actual == expected_d
    return {
        "dimension": 5,
        "independent_coordinates": independent,
        "fraction_coefficient_domain": True,
        "h_formula": "h_i = -sum_j B_inv[i,j] * (d_j + Q_j/N)",
        "residual_formula": "residual_i = d_i + sum_j B[i,j] h_j + Q_i/N",
        "h_polynomials": [_poly_serializable(value) for value in h],
        "raw_residual_polynomials": [_poly_serializable(value) for value in residual],
        "independent_coordinate_presence": independent_presence,
        "coefficient_witnesses": coefficient_witnesses,
        "conditional_on_B_inverse": True,
        "B_inverse_relation_checks": inverse_relation_checks,
        "all_fraction_checks_passed": all(inverse_relation_checks.values()),
        "terminal_h_rhs": "2^(-3/2) * A * residual_I * N_inv",
        "terminal_h_exponent": "-3/2",
        "terminal_h_row": "I",
    }


def _endpoint_support_proof() -> dict[str, Any]:
    arguments = {str(center): (Fraction(1) - center) / Fraction(1, 40) for center in CENTERS}
    _require(all(value > 1 for value in arguments.values()), "x=2 bump argument did not leave every support")
    return {
        "centers": [str(value) for value in CENTERS],
        "ell_over_log2": "1/40",
        "x": "2",
        "arguments_at_x2": {key: str(value) for key, value in arguments.items()},
        "arguments_strictly_outside_unit_support": all(value > 1 for value in arguments.values()),
        "raw_beta_definition": RAW_BETA_DEFINITION,
        "delta_E_support_formula": "delta_E = A * F, F = sum_i h[i+2] * bump_i / N",
        "conditional_on_raw_beta_zero_outside_abs_t_lt_1": True,
        "conditional_on_same_source_beta_binding": True,
        "endpoint_delta_E_zero_from_support_conditional": True,
        "actual_delta_E_zero_proved": False,
        "actual_endpoint_integral_identity_proved": False,
    }


def recover(root: str | Path) -> dict[str, Any]:
    """Audit pinned compact source and return conditional exact evidence."""
    root = Path(root)
    compact_text, compact_sha = _read_pinned(root, COMPACT_SOURCE, COMPACT_SOURCE_SHA256)
    controls_text, controls_sha = _read_pinned(root, CONTROL_SOURCE, CONTROL_SOURCE_SHA256)
    compact_tree = _parse(compact_text, COMPACT_SOURCE)
    controls_tree = _parse(controls_text, CONTROL_SOURCE)
    source_checks: dict[str, bool] = {}
    source_checks.update({f"targets.{key}": value for key, value in _check_normalized_targets(compact_tree).items()})
    source_checks.update({f"nonlinear.{key}": value for key, value in _check_nonlinear(compact_tree).items()})
    source_checks.update({f"implicit_rows.{key}": value for key, value in _check_implicit_rows(compact_tree).items()})
    source_checks.update({f"terminal_constructor.{key}": value for key, value in _check_terminal_constructor(compact_tree).items()})
    source_checks.update({f"endpoint_support.{key}": value for key, value in _check_bump_support(compact_tree, controls_tree).items()})
    algebra = _fraction_terminal_algebra()
    endpoint = _endpoint_support_proof()
    checks_passed = bool(all(source_checks.values()) and algebra["all_fraction_checks_passed"] and endpoint["arguments_strictly_outside_unit_support"])
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "source_ast_and_fraction_conditional" if checks_passed else "failed",
        "source_only": True,
        "source_paths": [COMPACT_SOURCE, CONTROL_SOURCE],
        "source_sha256": {COMPACT_SOURCE: compact_sha, CONTROL_SOURCE: controls_sha},
        "source_ast_checks": source_checks,
        "checks_passed": checks_passed,
        "fraction_algebra": algebra,
        "endpoint_support": endpoint,
        "conditions": [
            "B_inv is the inverse of the same source B on the same exact compact branch",
            "the source endpoint integral identity applies to the same h and residual rows",
            "RATES['h'] is the source rate 3/2 for the h/I terminal row",
            "raw beta is zero outside abs(t)<1 and controls and profiles use the same source binding",
        ],
        "actual_graph_verified": False,
        "actual_function_identity_proved": False,
        "actual_derivative_identity_proved": False,
        "actual_terminal_h_zero_proved": False,
        "actual_terminal_rhs_zero_proved": False,
        "terminal_h_unconditional": False,
        "actual_source_values_computed": False,
        "NS_residual": False,
        "runtime_model_metadata": None,
    }


__all__ = [
    "COMPONENT",
    "COMPACT_SOURCE",
    "CONTROL_SOURCE",
    "SCHEMA",
    "SourceContractError",
    "recover",
]
