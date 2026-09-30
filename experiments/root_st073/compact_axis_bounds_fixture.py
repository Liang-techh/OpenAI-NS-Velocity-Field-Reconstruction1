"""Resolved compact-interval fixture for ``axis_norm_bounds``.

The fixture checks the real derivatives directly with ``mp.diff`` at two
compact radii.  It deliberately disables endpoint ``G`` evaluations because
those calls perform adaptive MP quadrature and are not needed to validate the
radius-weighted majorants.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets  # noqa: E402
from lei_ren_part1_paper_axis_norm_bounds import (  # noqa: E402
    _direct_gradient,
    _nstr,
    axis_norm_bounds,
)


def _compact_points(radius: mp.mpf, z0: mp.mpf) -> list[mp.mpf]:
    points = [
        -mp.mpf("0.8") * radius,
        -mp.mpf("0.2") * radius,
        mp.mpf("0.2") * radius,
        mp.mpf("0.8") * radius,
    ]
    if abs(z0) <= radius:
        points.insert(2, z0)
    return points


def _check_radius(axis: Any, radius: Any) -> dict[str, Any]:
    radius_mp = mp.mpf(str(radius))
    result = axis_norm_bounds(
        axis,
        axial_radius=radius_mp,
        include_endpoint_values=False,
    )
    rows = []
    covered = True
    max_ratio = mp.mpf(0)
    for point in _compact_points(radius_mp, mp.mpf(axis.Z0)):
        derivatives = []
        for order, bound in enumerate(result["G_derivative_bounds"]):
            actual = mp.diff(lambda value: _direct_gradient(axis, value), point, order)
            ratio = abs(actual) / bound if bound else mp.mpf(0)
            max_ratio = max(max_ratio, ratio)
            item = {
                "order": order,
                "value": _nstr(actual),
                "bound": _nstr(bound),
                "ratio": _nstr(ratio),
                "bound_covers": bool(abs(actual) <= bound),
            }
            covered = covered and item["bound_covers"]
            derivatives.append(item)
        rows.append({"Z": _nstr(point), "derivatives": derivatives})

    anchor_expected = radius_mp + abs(mp.mpf(axis.Z0))
    anchor_derivative_expected = max(radius_mp, abs(mp.mpf(axis.Z0)))
    anchor_path_matches = result["G_anchor_path_length"] == anchor_expected
    anchor_derivative_radius_covers_path = (
        result["G_anchor_derivative_radius"] == anchor_derivative_expected
        and result["G_anchor_derivative_radius"] >= radius_mp
        and result["G_anchor_derivative_radius"] >= abs(mp.mpf(axis.Z0))
    )
    return {
        "axial_radius": _nstr(radius_mp),
        "interval": [_nstr(-radius_mp), _nstr(radius_mp)],
        "endpoint_values_evaluated": result["endpoint_values_evaluated"],
        "G_anchor_path_length": _nstr(result["G_anchor_path_length"]),
        "G_anchor_path_matches_radius_plus_abs_Z0": bool(anchor_path_matches),
        "G_anchor_derivative_radius": _nstr(
            result["G_anchor_derivative_radius"]
        ),
        "G_anchor_derivative_radius_covers_path": bool(
            anchor_derivative_radius_covers_path
        ),
        "G_value_upper_covers_anchor_bound": bool(
            result["G_value_upper"] >= result["G_anchor_value_upper"]
        ),
        "derivative_checks": rows,
        "all_derivative_checks_covered": covered,
        "max_derivative_ratio": _nstr(max_ratio),
        "scalar_rounding_enclosed": result["scalar_rounding_enclosed"],
        "source_error_enclosed": result["source_error_enclosed"],
    }


def run_fixture() -> dict[str, Any]:
    with mp.workdps(180):
        axis = RegularCoreAxisJets(
            j="1e-14",
            Lambda="1e36",
            delta="1e-200",
            precision=180,
        )
        radii = [_check_radius(axis, "0.8"), _check_radius(axis, "1e-30")]
        default = axis_norm_bounds(
            axis,
            axial_radius=1,
            include_endpoint_values=False,
        )
        default_legacy_value_bound = 2 * default["G_derivative_bounds"][0]
        return {
            "source": "experiments/root_st073/lei_ren_part1_paper_axis_norm_bounds.py",
            "parameters": {
                "j": _nstr(axis.j),
                "Lambda": _nstr(axis.Lambda),
                "delta": _nstr(axis.delta),
                "precision": axis.precision,
            },
            "Z0": _nstr(axis.Z0),
            "radii": radii,
            "default_radius_preserves_legacy_G_value_bound": bool(
                default["G_value_upper"] == default_legacy_value_bound
            ),
            "default_endpoint_values_evaluated": default["endpoint_values_evaluated"],
            "all_checks_passed": bool(
                all(row["all_derivative_checks_covered"] for row in radii)
                and all(row["G_anchor_path_matches_radius_plus_abs_Z0"] for row in radii)
                and all(row["G_anchor_derivative_radius_covers_path"] for row in radii)
                and all(row["G_value_upper_covers_anchor_bound"] for row in radii)
            ),
            "scope": (
                "resolved real-axis derivative replay only; scalar rounding, source "
                "quadrature, mixed radial, complex A_Omega, pressure, nonlinear core, "
                "cone, PDE, and endpoint-strip errors are not enclosed"
            ),
        }


def _json_value(value: Any) -> Any:
    if isinstance(value, mp.mpf):
        return _nstr(value)
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


if __name__ == "__main__":
    receipt = run_fixture()
    output = HERE / "compact_axis_bounds_fixture.json"
    output.write_text(json.dumps(_json_value(receipt), indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "all_checks_passed": receipt["all_checks_passed"],
                "Z0": receipt["Z0"],
                "radii": [row["axial_radius"] for row in receipt["radii"]],
            }
        ),
        flush=True,
    )


__all__ = ["run_fixture"]
