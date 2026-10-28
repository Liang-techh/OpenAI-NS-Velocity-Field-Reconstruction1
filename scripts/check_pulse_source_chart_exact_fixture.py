"""Static-source fixture for the exact three-chart pulse coordinate kernel.

The fixture is intentionally independent of the older pulse modules and all
source-expression consumers.  It uses an independent positive rational
synthetic ``mu`` for exact symbolic and numeric checks, then separately checks
that the real nested source ``mu`` remains a SymPy expression rather than a
floating value.  It does not evaluate physical fields, compare graph nodes,
claim pulse equality, or claim error elimination.

This file is prepared for a later single execution by the parent.  The current
task deliberately does not execute, import, compile, or test it.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from openai_ns_reconstruction.pulse_source_chart_exact import (  # noqa: E402
    MU_SYMBOL,
    build_exact_chart_system,
    source_mu_expression,
)


def _symbolic_zero(expression: sp.Expr) -> bool:
    return sp.simplify(expression) == sp.Integer(0)


def _numeric_zero(expression: sp.Expr) -> bool:
    # The synthetic substitution is rational.  ``N`` is used only for an
    # independent numeric view of that exact rational check; it is never used
    # for the real source mu.
    return sp.N(expression, 60) == sp.Integer(0)


def _require(controls: dict[str, bool], name: str, value: bool) -> None:
    controls[name] = bool(value)
    if not value:
        raise AssertionError(f"source chart fixture control failed: {name}")


def run_fixture() -> dict[str, Any]:
    # An independent exact positive rational, unrelated to the real source
    # expression.  No Python float is introduced anywhere in this fixture.
    synthetic_mu = sp.Rational(1, 4)
    synthetic_log_rp = sp.Rational(7, 5)
    system = build_exact_chart_system(
        mu=synthetic_mu,
        log_rp=synthetic_log_rp,
        ordinary_z_independent=True,
    )
    charts = system.chart_map
    gap = charts["gap"]
    end = charts["end"]
    xi = gap.coordinate
    s = gap.reparameterized_coordinate
    assert s is not None
    assert gap.xi_from_s is not None
    assert gap.s_from_xi is not None
    assert gap.reparameterized_to_y is not None
    assert gap.reparameterized_bounds is not None
    assert gap.reparameterized_theta_exponent_original is not None
    assert gap.reparameterized_theta_exponent_log_rewrite is not None

    controls: dict[str, bool] = {}

    # The source-coordinate endpoint expression remains exactly 13 - 4 * mu.
    _require(
        controls,
        "gap_upper_retains_exact_13_minus_4mu",
        _symbolic_zero(
            gap.source_bounds.upper
            - (sp.Integer(13) - sp.Integer(4) * synthetic_mu)
        ),
    )

    # xi = 13 + mu*s and y = xi/mu = 13/mu + s are exact identities.
    gap_y_from_xi = gap.coordinate_to_y
    gap_y_from_s = gap.reparameterized_to_y
    _require(
        controls,
        "gap_coordinate_map_symbolic",
        _symbolic_zero(
            gap_y_from_xi.subs(xi, gap.xi_from_s) - gap_y_from_s
        ),
    )
    _require(
        controls,
        "gap_coordinate_map_numeric_rational",
        _numeric_zero(
            (gap_y_from_xi.subs(xi, gap.xi_from_s) - gap_y_from_s).subs(
                {MU_SYMBOL: synthetic_mu}
            )
        ),
    )

    # The lower/upper endpoint map is [11, 13 - 4 * mu] -> [-2/mu, -4].
    s_lower = gap.reparameterized_bounds.lower
    s_upper = gap.reparameterized_bounds.upper
    _require(
        controls,
        "gap_lower_endpoint_map",
        _symbolic_zero(gap.xi_from_s.subs(s, s_lower) - sp.Integer(11)),
    )
    _require(
        controls,
        "gap_upper_endpoint_map",
        _symbolic_zero(
            gap.xi_from_s.subs(s, s_upper)
            - (sp.Integer(13) - sp.Integer(4) * synthetic_mu)
        ),
    )
    _require(
        controls,
        "gap_reparameterized_bounds_exact",
        s_lower == -sp.Integer(2) / synthetic_mu and s_upper == -sp.Integer(4),
    )

    # dxi/ds = mu, so dxi/mu = ds exactly.
    _require(
        controls,
        "gap_measure_symbolic",
        _symbolic_zero(
            sp.diff(gap.xi_from_s, s) / synthetic_mu
            - gap.reparameterized_measure_density
        ),
    )
    _require(
        controls,
        "gap_measure_numeric_rational",
        _numeric_zero(
            (
                sp.diff(gap.xi_from_s, s) / synthetic_mu
                - gap.reparameterized_measure_density
            ).subs({MU_SYMBOL: synthetic_mu})
        ),
    )

    # beta = 1/2 + mu and the log-space theta exponent are identical before
    # and after rewriting, both in the reparameterized gap and on end.
    _require(
        controls,
        "gap_theta_rewrite_symbolic",
        _symbolic_zero(
            gap.reparameterized_theta_exponent_original
            - gap.reparameterized_theta_exponent_log_rewrite
        ),
    )
    _require(
        controls,
        "gap_theta_rewrite_numeric_rational",
        _numeric_zero(
            (
                gap.reparameterized_theta_exponent_original
                - gap.reparameterized_theta_exponent_log_rewrite
            ).subs({MU_SYMBOL: synthetic_mu})
        ),
    )
    _require(
        controls,
        "end_theta_rewrite_symbolic",
        _symbolic_zero(
            end.theta_exponent_original - end.theta_exponent_log_rewrite
        ),
    )
    _require(
        controls,
        "end_theta_rewrite_numeric_rational",
        _numeric_zero(
            (
                end.theta_exponent_original - end.theta_exponent_log_rewrite
            ).subs({MU_SYMBOL: synthetic_mu})
        ),
    )

    # Gap reparameterization and end have the same y, log-radius and physical
    # measure formulas.  This is a coordinate identity, not field equality.
    _require(
        controls,
        "gap_end_y_map_same",
        _symbolic_zero(gap.reparameterized_to_y - end.coordinate_to_y),
    )
    _require(
        controls,
        "gap_end_log_radius_map_same",
        _symbolic_zero(
            gap.reparameterized_log_radius - end.log_radius
        ),
    )
    _require(
        controls,
        "gap_end_measure_map_same",
        gap.reparameterized_measure_density == end.physical_measure_density
        and gap.reparameterized_measure_identity == "dy = ds"
        and end.physical_measure_identity == "dy = ds",
    )

    _require(
        controls,
        "all_three_charts_returned",
        tuple(system.chart_map) == ("main", "gap", "end"),
    )
    _require(
        controls,
        "ordinary_z_declared_explicitly",
        system.ordinary_z_independence.declared
        and all(chart.ordinary_z_independence.declared for chart in system.charts),
    )
    _require(
        controls,
        "ordinary_z_certificate_not_automatic",
        not system.ordinary_z_independence.automatic_certificate
        and all(
            not chart.ordinary_z_independence.automatic_certificate
            for chart in system.charts
        ),
    )
    _require(
        controls,
        "radius_value_never_materialized",
        not system.radius_value_materialized
        and all(not chart.radius_value_materialized for chart in system.charts),
    )

    # The actual source expression is retained as exact SymPy syntax.  This
    # check is deliberately independent of the synthetic rational path above.
    source_mu = source_mu_expression()
    _require(
        controls,
        "real_source_mu_is_sympy_expression",
        isinstance(source_mu, sp.Expr) and not isinstance(source_mu, float),
    )
    _require(
        controls,
        "real_source_mu_not_sympy_float",
        not source_mu.is_Float,
    )
    _require(
        controls,
        "real_source_mu_retains_exp40",
        source_mu.has(sp.exp(sp.Integer(40))),
    )
    _require(
        controls,
        "real_source_substitution_is_exact",
        system.source_substitution[MU_SYMBOL] == source_mu,
    )
    source_system = build_exact_chart_system(ordinary_z_independent=True)
    source_gap = source_system.chart_map["gap"].source_substituted()
    _require(
        controls,
        "source_gap_upper_retains_mu_term",
        source_gap["source_upper"] == sp.Integer(13) - sp.Integer(4) * source_mu,
    )
    _require(
        controls,
        "source_gap_reparameterized_endpoint_exact",
        source_gap["reparameterized_lower"] == -sp.Integer(2) / source_mu
        and source_gap["reparameterized_upper"] == -sp.Integer(4),
    )
    _require(
        controls,
        "real_source_substitution_is_not_materialized",
        not system.as_dict()["radius_value_materialized"]
        and "Float" not in system.as_dict()["source_mu"],
    )

    report: dict[str, Any] = {
        "passed": all(controls.values()),
        "synthetic_only": True,
        "synthetic_mu": sp.sstr(synthetic_mu),
        "source_values_selected": False,
        "physical_fields_bound": False,
        "pulse_equality_proved": False,
        "error_elimination_claimed": False,
        "radius_value_materialized": system.radius_value_materialized,
        "ordinary_z_independence": system.ordinary_z_independence.as_dict(),
        "source_mu_recipe": system.source_mu_recipe,
        "source_mu": sp.sstr(source_mu),
        "synthetic_charts": system.as_dict(),
        "source_charts": source_system.as_dict(include_source_view=True),
        "controls": controls,
    }
    if not report["passed"]:
        raise AssertionError(f"exact source chart fixture failed: {report}")
    return report


if __name__ == "__main__":
    print(json.dumps(run_fixture(), indent=2, sort_keys=True))
