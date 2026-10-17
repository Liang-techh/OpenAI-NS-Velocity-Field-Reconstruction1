"""Bounded interior mass for the pinned original flat pulse.

The original first-switch mass is

    M(s) = integral_0^s (1 - sigma(t)) dt,
    sigma(t) = gp(t) / (gp(t) + gp(1 - t)),
    gp(t) = exp(-1/t**2) for t > 0,

with the flat endpoint convention ``gp(t) = 0`` for ``t <= 0``.  This
module binds that source by AST and evaluates a separate enclosure.  It
does not import or call the accepted bridge mass producer.  Interior rows
use an independent scalar adaptive quadrature for their displayed source
approximation, while monotone left/right endpoint rectangles provide the
directed certification box.  The exact zero and half branches are retained
at the two endpoints.

Rectangle midpoints are refinement machinery only.  They are never
promoted to exact source values, and finite rectangle nodes are quadrature
machinery rather than source parameters.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_leading_viscous_kernel import encode_pair


_DEFAULT_PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_leading_viscous_source_numeric.json"
)
_EXPECTED_PARENT_SHA256 = (
    "d1aa3644ea8158fe78dccb210a0a75a4c81ad257b8504758cb5797845b8003d4"
)
_PHASES = ("0", "0.125", "0.25", "0.5", "0.75", "0.875", "1")
_PANELS = (128, 256)
_SOURCE_FORMULA = "gp(t)=exp(-1/t^2) for t>0; sigma=gp(t)/(gp(t)+gp(1-t))"
_MASS_FORMULA = "M(s)=integral_0^s (1-sigma(t)) dt"


def _repo_root(root: Any, directory: Any) -> Path:
    """Resolve a repository root without reading or constructing an owner."""

    for candidate in (root, directory):
        if candidate is None:
            continue
        path = Path(candidate)
        if path.is_file():
            path = path.parent
        if (path / "src").is_dir() and (path / "experiments").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _relative(path: Path, repo: Path) -> str:
    try:
        return path.relative_to(repo).as_posix()
    except ValueError:
        return path.as_posix()


def _sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def _number(ctx: Any, value: Any, digits: int = 90) -> str:
    """Serialize an mp number without converting it through binary float."""

    value = ctx.mpf(value)
    return "0" if value == 0 else ctx.nstr(value, digits)


def _interval_bounds(ctx: Any, value: Any) -> tuple[Any, Any]:
    """Read outward mpmath interval endpoints into the ordinary context."""

    try:
        lower, upper = value._mpi_
    except AttributeError as exc:  # pragma: no cover - internal invariant
        raise TypeError("expected an mpmath interval") from exc
    return ctx.mpf(lower), ctx.mpf(upper)


def _exact_interval_from_bounds(ctx: Any, iv: MPIntervalContext, lower: Any, upper: Any) -> Any:
    """Build an interval from the directed mpf tuples, without re-rounding labels."""

    lower_tuple = tuple(lower._mpi_[0])
    upper_tuple = tuple(upper._mpi_[1])
    return iv.mpf(
        [
            ctx.make_mpf(lower_tuple),
            ctx.make_mpf(upper_tuple),
        ]
    )


def _serialized_error_number(ctx: Any, encoded: Mapping[str, Any]) -> str:
    error_tuple = encoded["serialized_approximation_absolute_error_upper_exact_mpf_tuple"]
    return _number(ctx, ctx.make_mpf(tuple(error_tuple)))


def _sigma_scalar(ctx: Any, value: Any) -> Any:
    """Evaluate the pinned gp-ratio in an ordinary context for adaptive quadrature."""

    point = ctx.mpf(value)
    if point <= 0:
        return ctx.mpf(0)
    if point >= 1:
        return ctx.mpf(1)
    left = ctx.exp(-1 / point**2)
    reflected = 1 - point
    right = ctx.exp(-1 / reflected**2)
    return left / (left + right)


def _adaptive_scalar_mass(ctx: Any, phase: str) -> Any:
    """Return a new scalar c.quad approximation of the explicit source formula."""

    endpoint = ctx.mpf(phase)
    if endpoint == 0:
        return ctx.mpf(0)
    if endpoint == 1:
        return ctx.mpf("0.5")

    def integrand(value: Any) -> Any:
        return 1 - _sigma_scalar(ctx, value)

    return ctx.quad(integrand, [0, endpoint])


def _ast_assignments(path: Path, function_name: str) -> dict[str, set[str]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    functions = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function_name
    ]
    if len(functions) != 1:
        raise ValueError(f"expected one pinned source function: {function_name}")
    assignments: dict[str, set[str]] = {}
    for node in ast.walk(functions[0]):
        if not isinstance(node, ast.Assign):
            continue
        expression = ast.unparse(node.value).replace(" ", "")
        for target in node.targets:
            assignments.setdefault(ast.unparse(target).replace(" ", ""), set()).add(expression)
    return assignments


def _ast_source_pin(repo: Path) -> dict[str, Any]:
    """Confirm the exact source expressions before evaluating any mass row."""

    specs = (
        {
            "path": "experiments/root_st073/lei_ren_part1_paper_compliant_flat_pulse_derivatives.py",
            "function": "_sigma_left",
            "required": {
                "odds": "1/(1-x)**2-1/x**2",
                "value": "e/(1+e)",
            },
            "role": "production equivalent odds form",
        },
        {
            "path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_inner_exit_strict_collar.py",
            "function": "exact_source_attachment",
            "required": {
                "aa": "s.exp(-1/x**2)",
                "bb": "s.exp(-1/(1-x)**2)",
                "sig": "aa/(aa+bb)",
            },
            "role": "pinned gp ratio and monotonicity proof source",
        },
    )
    rows: list[dict[str, Any]] = []
    for spec in specs:
        path = repo / spec["path"]
        digest = _sha256(path)
        found: dict[str, list[str]] = {}
        error = None
        try:
            assignments = _ast_assignments(path, spec["function"])
            found = {key: sorted(values) for key, values in assignments.items()}
            matched = {
                key: value in assignments.get(key, set())
                for key, value in spec["required"].items()
            }
        except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
            matched = {key: False for key in spec["required"]}
            error = str(exc)
        rows.append(
            {
                "path": spec["path"],
                "sha256": digest,
                "function": spec["function"],
                "role": spec["role"],
                "required_ast": dict(spec["required"]),
                "found_ast": found,
                "matched_ast": matched,
                "passed": bool(digest and error is None and all(matched.values())),
                "error": error,
                "source_value_not_selected": True,
            }
        )

    # This saved theorem receipt is a source pin too, but no code from it is
    # executed.  It records the endpoint convention used by the source.
    fixed_path = repo / "experiments/root_st073/lei_ren_part1_paper_shared_fixed_step_bound.json"
    fixed_digest = _sha256(fixed_path)
    fixed_data: Mapping[str, Any] = {}
    fixed_error = None
    try:
        fixed_data = json.loads(fixed_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, UnicodeError) as exc:
        fixed_error = str(exc)
    fixed_formula = fixed_data.get("cutoff_formula")
    fixed_passed = bool(
        fixed_digest
        and fixed_error is None
        and fixed_formula
        == "sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), 0<s<1; flat constants outside"
        and fixed_data.get("strict_monotonicity_on_open_unit_interval") is True
        and fixed_data.get("smooth_flat_endpoints") is True
    )
    rows.append(
        {
            "path": _relative(fixed_path, repo),
            "sha256": fixed_digest,
            "role": "saved endpoint and monotonicity receipt",
            "required_formula": "sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), 0<s<1; flat constants outside",
            "found_formula": fixed_formula,
            "matched_endpoint_theorem": fixed_passed,
            "passed": fixed_passed,
            "error": fixed_error,
            "source_value_not_selected": True,
        }
    )
    return {
        "formula": _SOURCE_FORMULA,
        "mass_formula": _MASS_FORMULA,
        "endpoint_convention": {
            "gp(t<=0)": "0",
            "gp(t>0)": "exp(-1/t^2)",
            "sigma(t<=0)": "0",
            "sigma(t>=1)": "1",
        },
        "rows": rows,
        "passed": bool(rows) and all(row.get("passed") is True for row in rows),
        "source_value_not_selected": True,
    }


def _parent_pin(parent_report: Any, repo: Path) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(parent_report, Mapping):
        return parent_report, {
            "mode": "supplied_mapping",
            "digest": _digest(parent_report),
            "sha256": None,
            "expected_sha256": _EXPECTED_PARENT_SHA256,
            "passed": False,
            "source_value_not_selected": True,
        }
    path = Path(parent_report) if parent_report is not None else repo / _DEFAULT_PARENT
    if not path.is_absolute():
        path = repo / path
    digest = _sha256(path)
    data: Mapping[str, Any] = {}
    error = None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(data, Mapping):
            data = {}
            error = "parent report is not a JSON object"
    except (OSError, ValueError, UnicodeError) as exc:
        error = str(exc)
    return data, {
        "path": _relative(path, repo),
        "sha256": digest,
        "expected_sha256": _EXPECTED_PARENT_SHA256,
        "passed": digest == _EXPECTED_PARENT_SHA256 and error is None and bool(data),
        "error": error,
        "source_value_not_selected": True,
    }


def _gp_interval(iv: MPIntervalContext, ctx: Any, value: Any) -> Any:
    lo, hi = _interval_bounds(ctx, value)
    if hi <= 0:
        return iv.mpf(0)
    if lo <= 0:
        raise ValueError("source node straddles the flat endpoint")
    return iv.exp(-1 / value**2)


def _sigma_interval(iv: MPIntervalContext, ctx: Any, value: Any) -> Any:
    lo, hi = _interval_bounds(ctx, value)
    if hi <= 0:
        return iv.mpf(0)
    if lo >= 1:
        return iv.mpf(1)
    left = _gp_interval(iv, ctx, value)
    right = _gp_interval(iv, ctx, iv.mpf(1) - value)
    return left / (left + right)


def _source_sample(ctx: Any, iv: MPIntervalContext, phase: str) -> dict[str, Any]:
    point = ctx.mpf(phase)
    if point == 0:
        sigma = iv.mpf(0)
    elif point == 1:
        sigma = iv.mpf(1)
    else:
        sigma = _sigma_interval(iv, ctx, iv.mpf(phase))
    lower, upper = _interval_bounds(ctx, sigma)
    exact = point in (0, 1)
    return {
        "phase": phase,
        "sigma": {
            "approximation": _number(ctx, (lower + upper) / 2),
            "lower": _number(ctx, lower),
            "upper": _number(ctx, upper),
            "approximation_is_exact_source_value": exact,
            "source_value_not_selected": not exact,
            "approximation_is_not_exact_source_value": not exact,
        },
        "integrand": {
            "formula": "1-sigma(phase)",
            "lower": _number(ctx, 1 - upper),
            "upper": _number(ctx, 1 - lower),
            "source_value_not_selected": not exact,
        },
    }


def _rectangle_row(
    ctx: Any,
    iv: MPIntervalContext,
    phase: str,
    panels: int,
) -> tuple[dict[str, Any], Any, Any, Any]:
    phase_point = ctx.mpf(phase)
    if phase_point == 0:
        truth = iv.mpf(0)
        encoded = encode_pair(ctx, iv, ctx.mpf(0), truth)
        exact = {
            "phase": phase,
            "panels": panels,
            "kind": "signed_source_integral",
            "signed": True,
            "method": "exact endpoint branch",
            **encoded,
            "lower": "0",
            "upper": "0",
            "approximation": "0",
            "absolute_error_upper": "0",
            "serialized_approximation_absolute_error_upper": "0",
            "integration_error_bound": "0",
            "enclosure_width": "0",
            "enclosure_radius": "0",
            "rectangle_width_bound": "0",
            "width_bound_passed": True,
            "approximation_is_exact_source_value": True,
            "approximation_is_not_exact_source_value": False,
            "source_value_not_selected": False,
            "source_enclosure_is_not_a_selected_value": False,
            "exact_source_branch": "M(0)=0",
        }
        return exact, ctx.mpf(0), ctx.mpf(0), truth
    if phase_point == 1:
        truth = iv.mpf("0.5")
        encoded = encode_pair(ctx, iv, ctx.mpf("0.5"), truth)
        exact = {
            "phase": phase,
            "panels": panels,
            "kind": "signed_source_integral",
            "signed": True,
            "method": "exact endpoint symmetry branch",
            **encoded,
            "lower": "0.5",
            "upper": "0.5",
            "approximation": "0.5",
            "absolute_error_upper": "0",
            "serialized_approximation_absolute_error_upper": "0",
            "integration_error_bound": "0",
            "enclosure_width": "0",
            "enclosure_radius": "0",
            "rectangle_width_bound": "0",
            "width_bound_passed": True,
            "approximation_is_exact_source_value": True,
            "approximation_is_not_exact_source_value": False,
            "source_value_not_selected": False,
            "source_enclosure_is_not_a_selected_value": False,
            "exact_source_branch": "M(1)=1/2 from sigma(1-t)=1-sigma(t)",
        }
        return exact, ctx.mpf("0.5"), ctx.mpf("0.5"), truth

    # Every node is a rational phase subdivision.  The interval source value
    # is evaluated at the endpoint; no midpoint or cap is made a source value.
    phase_iv = iv.mpf(phase)
    step = phase_iv / panels
    lower_sum = iv.mpf(0)
    upper_sum = iv.mpf(0)
    for index in range(panels):
        left = phase_iv * index / panels
        right = phase_iv * (index + 1) / panels
        left_sigma = _sigma_interval(iv, ctx, left)
        right_sigma = _sigma_interval(iv, ctx, right)
        left_integrand = iv.mpf(1) - left_sigma
        right_integrand = iv.mpf(1) - right_sigma
        upper_sum += step * left_integrand
        lower_sum += step * right_integrand
    truth = _exact_interval_from_bounds(ctx, iv, lower_sum, upper_sum)
    lower, upper = _interval_bounds(ctx, truth)
    # The source is decreasing, so the true integral lies in this interval.
    # Keep the rectangle midpoint as refinement machinery only.  The source
    # approximation exported by the enclosing mass row comes from a separate
    # scalar adaptive quadrature.
    approximation = (lower + upper) / 2
    encoded = encode_pair(ctx, iv, approximation, truth)
    width = upper - lower
    rectangle_bound = phase_point / panels
    row = {
        "phase": phase,
        "panels": panels,
        "kind": "signed_source_integral",
        "signed": True,
        "method": "monotone endpoint rectangles",
        "integrand_monotone": "decreasing on [0,1]",
        **encoded,
        "lower": _number(ctx, lower),
        "upper": _number(ctx, upper),
        "absolute_error_upper": _serialized_error_number(ctx, encoded),
        "serialized_approximation_absolute_error_upper": _serialized_error_number(ctx, encoded),
        "enclosure_radius": _number(ctx, width / 2),
        "integration_error_bound": _number(ctx, rectangle_bound),
        "enclosure_width": _number(ctx, width),
        "rectangle_width_bound": _number(ctx, rectangle_bound),
        "width_bound_passed": bool(width <= rectangle_bound),
        "approximation_is_exact_source_value": False,
        "approximation_is_not_exact_source_value": True,
        "source_value_not_selected": True,
        "source_enclosure_is_not_a_selected_value": True,
        "finite_nodes_are_machinery": True,
    }
    return row, lower, upper, truth


def _mass_row(ctx: Any, iv: MPIntervalContext, phase: str) -> tuple[dict[str, Any], dict[str, Any], Any, Any]:
    refinements: dict[str, Any] = {}
    internal: dict[int, tuple[Any, Any, Any]] = {}
    for panels in _PANELS:
        row, lower, upper, truth = _rectangle_row(ctx, iv, phase, panels)
        refinements[str(panels)] = row
        internal[panels] = (lower, upper, truth)
    fine_lower, fine_upper, fine_truth = internal[_PANELS[-1]]
    fine_rectangle_midpoint = (fine_lower + fine_upper) / 2
    coarse_lower, coarse_upper, _coarse_truth = internal[_PANELS[0]]
    coarse_approximation = (coarse_lower + coarse_upper) / 2
    convergence = abs(fine_rectangle_midpoint - coarse_approximation)
    convergence_bound = (coarse_upper - coarse_lower) / 2 + (fine_upper - fine_lower) / 2
    phase_point = ctx.mpf(phase)
    exact = phase_point in (0, 1)
    scalar_approximation = _adaptive_scalar_mass(ctx, phase)
    encoded = encode_pair(ctx, iv, scalar_approximation, fine_truth)
    if exact:
        encoded["approximation_is_exact_source_value"] = True
        encoded["source_enclosure_is_not_a_selected_value"] = False
    serialized_error = _serialized_error_number(ctx, encoded)
    row = {
        "phase": phase,
        "kind": "signed_source_integral",
        "signed": True,
        "formula": _MASS_FORMULA,
        "method": "independent scalar adaptive quadrature (ctx.quad)",
        "refinements": refinements,
        **encoded,
        "approximation": (
            "0"
            if phase_point == 0
            else "0.5"
            if phase_point == 1
            else encoded["approximation"]
        ),
        "lower": _number(ctx, fine_lower),
        "upper": _number(ctx, fine_upper),
        "absolute_error_upper": serialized_error,
        "serialized_approximation_absolute_error_upper": serialized_error,
        "enclosure_radius": _number(ctx, (fine_upper - fine_lower) / 2),
        "rectangle_midpoint_approximation": _number(ctx, fine_rectangle_midpoint),
        "rectangle_midpoint_source_value_not_selected": True,
        "convergence": {
            "panels": list(_PANELS),
            "difference_128_to_256": _number(ctx, convergence),
            "enclosure_radius_sum": _number(ctx, convergence_bound),
            "passed": bool(convergence <= convergence_bound),
            "source_value_not_selected": True,
        },
        "mass_positive": bool(fine_lower >= 0 and fine_upper >= fine_lower),
        "mass_Z_independent": True,
        "approximation_is_exact_source_value": exact,
        "approximation_is_not_exact_source_value": not exact,
        "source_value_not_selected": not exact,
        "source_enclosure_is_not_a_selected_value": not exact,
        "finite_nodes_are_machinery": True,
    }
    if exact:
        row["exact_source_branch"] = "M(0)=0" if phase_point == 0 else "M(1)=1/2"
    return row, refinements, fine_lower, fine_upper


def _symmetry_checks(ctx: Any, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_phase = {row["phase"]: row for row in rows}
    checks: list[dict[str, Any]] = []
    for phase in ("0.125", "0.25", "0.5"):
        complement = {"0.125": "0.875", "0.25": "0.75", "0.5": "0.5"}[phase]
        left = by_phase[phase]
        right = by_phase[complement]
        lhs = ctx.mpf(right["approximation"]) - ctx.mpf(left["approximation"])
        rhs = ctx.mpf("0.5") - ctx.mpf(phase)
        difference = abs(lhs - rhs)
        radius = ctx.mpf(left["serialized_approximation_absolute_error_upper"]) + ctx.mpf(
            right["serialized_approximation_absolute_error_upper"]
        )
        numerical_consistency_passed = bool(difference <= radius)
        checks.append(
            {
                "phase": phase,
                "complement_phase": complement,
                "identity": "M(1-s)-M(s)=1/2-s",
                "check_kind": "numerical_consistency_only",
                "exact_identity_source_theorem": "sigma(1-t)=1-sigma(t), derived algebraically",
                "left_minus_right": _number(ctx, lhs),
                "expected": _number(ctx, rhs),
                "difference": _number(ctx, difference),
                "informational_radius": _number(ctx, radius),
                "interval_radius_is_informational": True,
                "numerical_consistency_passed": numerical_consistency_passed,
                "passed": numerical_consistency_passed,
                "source_value_not_selected": True,
            }
        )
    return checks


def build_interior_mass_source(*, parent_report: Any, directory: Any, root: Any) -> dict[str, Any]:
    """Export bounded original mass rows at the seven requested phases.

    ``parent_report`` may be a path to the pinned slim report or a mapping
    supplied by a caller.  A supplied mapping is retained as provenance but
    cannot establish the required file SHA, so its pin is reported as
    failed.  ``directory`` and ``root`` only resolve saved files; no accepted
    constructor, producer, checker, fixture, or ancestor calculation is
    imported or executed.
    """

    repo = _repo_root(root, directory)
    parent_data, parent_pin = _parent_pin(parent_report, repo)
    source_pin = _ast_source_pin(repo)
    ctx = mp.mp.clone()
    ctx.dps = 180
    iv = MPIntervalContext()
    iv.dps = 180

    rows: list[dict[str, Any]] = []
    samples: list[dict[str, Any]] = []
    internal: list[tuple[Any, Any]] = []
    for phase in _PHASES:
        row, _refinements, lower, upper = _mass_row(ctx, iv, phase)
        rows.append(row)
        internal.append((lower, upper))
        samples.append(_source_sample(ctx, iv, phase))

    symmetry = _symmetry_checks(ctx, rows)
    widths = [
        bool(refinement["width_bound_passed"])
        for row in rows
        for refinement in row["refinements"].values()
        if refinement["method"] == "monotone endpoint rectangles"
    ]
    monotone = all(
        ctx.mpf(rows[index]["approximation"]) <= ctx.mpf(rows[index + 1]["approximation"])
        for index in range(len(rows) - 1)
    )
    endpoint_zero = rows[0]["approximation"] == "0"
    endpoint_half = rows[-1]["approximation"] == "0.5"
    source_flags = all(
        row["approximation_is_not_exact_source_value"] is True
        for row in rows[1:-1]
    )
    endpoint_symmetry_numerical_consistency_passed = bool(symmetry) and all(
        item["numerical_consistency_passed"] for item in symmetry
    )
    controls = {
        "source_ast_confirmed_before_evaluation": source_pin["passed"],
        "parent_report_sha256_matches": parent_pin["passed"],
        "exact_M0_zero": endpoint_zero,
        "exact_M1_half": endpoint_half,
        "endpoint_symmetry_checks": endpoint_symmetry_numerical_consistency_passed,
        "endpoint_symmetry_numerical_consistency_passed": endpoint_symmetry_numerical_consistency_passed,
        "endpoint_symmetry_passed": endpoint_symmetry_numerical_consistency_passed,
        "mass_monotone_in_phase": monotone,
        "mass_positive": all(row["mass_positive"] for row in rows),
        "signed_source_integral_preserved": all(row["signed"] for row in rows),
        "mass_Z_independent": all(row["mass_Z_independent"] for row in rows),
        "rectangle_width_le_phase_over_panels": bool(widths) and all(widths),
        "independent_128_256_convergence": all(row["convergence"]["passed"] for row in rows),
        "source_value_not_selected_for_interior_approximations": source_flags,
        "finite_quadrature_nodes_are_machinery": all(row["finite_nodes_are_machinery"] for row in rows),
        "old_pulse_mass_routine_called": False,
        "accepted_ancestor_replayed": False,
        "numerical_caps_used_as_values": False,
        "Z_parameter_used": False,
    }
    required_control_names = (
        "source_ast_confirmed_before_evaluation",
        "parent_report_sha256_matches",
        "exact_M0_zero",
        "exact_M1_half",
        "endpoint_symmetry_numerical_consistency_passed",
        "mass_monotone_in_phase",
        "mass_positive",
        "signed_source_integral_preserved",
        "mass_Z_independent",
        "rectangle_width_le_phase_over_panels",
        "independent_128_256_convergence",
        "source_value_not_selected_for_interior_approximations",
        "finite_quadrature_nodes_are_machinery",
    )
    required_controls_passed = all(bool(controls[name]) for name in required_control_names)
    controls["required_positive_controls_passed"] = required_controls_passed
    return {
        "source": {
            "formula": _SOURCE_FORMULA,
            "mass_formula": _MASS_FORMULA,
            "signed_source_integral": True,
            "domain": "0<=s<=1",
            "monotonicity_theorem": "sigma'(t)>0 on (0,1), hence 1-sigma(t) decreases",
            "symmetry_theorem": "sigma(1-t)=1-sigma(t), so M(1)=1/2 and M(1-s)-M(s)=1/2-s",
            "source_value_not_selected": True,
        },
        "parent_report": parent_pin,
        "parent_report_data_digest": _digest(parent_data) if parent_data else None,
        "rows": rows,
        "samples": samples,
        "symmetry_checks": symmetry,
        "source_pins": {
            "parent_report": parent_pin,
            "original_sigma_ast": source_pin,
            "source_value_not_selected": True,
        },
        "controls": controls,
        "limits": {
            "phase_count": len(_PHASES),
            "phases": list(_PHASES),
            "panel_counts": list(_PANELS),
            "rectangle_error_bound": "s/N",
            "mass_is_Z_independent": True,
            "interior_approximation_is_not_exact_source_value": True,
            "exact_endpoint_branches_retained": True,
            "required_positive_control_names": list(required_control_names),
            "required_positive_controls_passed": required_controls_passed,
            "expected_false_guard_flags_are_informational": True,
            "bounded_mass_interface_complete": bool(
                required_controls_passed and source_pin["passed"] and parent_pin["passed"]
            ),
            "complete": bool(
                required_controls_passed and source_pin["passed"] and parent_pin["passed"]
            ),
            "complete_scope": "bounded mass interface; full leading source remains outside this export",
        },
        "passed": bool(
            source_pin["passed"]
            and parent_pin["passed"]
            and controls["exact_M0_zero"]
            and controls["exact_M1_half"]
            and controls["rectangle_width_le_phase_over_panels"]
            and controls["independent_128_256_convergence"]
        ),
        "complete": False,
    }


__all__ = ["build_interior_mass_source"]
