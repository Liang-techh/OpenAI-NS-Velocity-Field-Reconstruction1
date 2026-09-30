"""Conditional C1 majorants for the finite five-bump response.

This module bounds the fixed finite ``FiveBumpMomentMap`` using a supplied
uniform C1 defect bound ``e_bound``.  It uses the paper's conditional estimate
``||A_m**-2||_C1 <= 40 / Pstar**2`` and the scalar inverse norm of the fixed
linear matrix.  The finite quadrature weights are data, so the result is a
conditional majorant rather than an analytic enclosure of the bump family or
of an actual source.

No large response roots are subtracted.  The contraction radius, formal
Taylor remainder, and tail are computed from positive majorant quantities.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
import operator
from typing import Any

import mpmath as mp


def _integer(value: Any, name: str) -> int:
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    return int(result)


def _mp(value: Any, name: str) -> mp.mpf:
    try:
        result = value if isinstance(value, mp.mpf) else mp.mpf(str(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real scalar") from exc
    if not mp.isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


def _nstr(value: Any, digits: int = 60) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (int, str)):
        return str(value)
    return mp.nstr(value, digits)


def _catalan(index: int) -> int:
    return int(mp.binomial(2 * index, index) / (index + 1))


def _matrix_value(matrix: Any, row: int, column: int) -> mp.mpf:
    try:
        value = matrix[row, column]
    except (TypeError, IndexError, KeyError):
        value = matrix[row][column]
    return _mp(value, "quadratic weight")


def _quadratic_pair_bounds(moment_map: Any, pstar: mp.mpf, amplitude_bound=None) -> tuple[list[list[mp.mpf]], dict[str, Any]]:
    """Build the symmetric column-pair l1 tensor bound for Q."""

    weights = getattr(moment_map, "quadratic_weights", None)
    if not isinstance(weights, Mapping):
        raise TypeError("moment_map must expose quadratic_weights mapping")
    required = ("fg_sqrt", "gg", "ff", "ff_over_x")
    if any(name not in weights for name in required):
        raise KeyError("quadratic_weights must contain fg_sqrt, gg, ff, ff_over_x")

    pair = [[mp.mpf(0) for _ in range(5)] for _ in range(5)]
    fg = weights["fg_sqrt"]
    gg = weights["gg"]
    ff = weights["ff"]
    ff_over_x = weights["ff_over_x"]
    axial_angular: list[list[mp.mpf]] = []
    axial_axial: list[list[mp.mpf]] = []
    angular_angular: list[list[mp.mpf]] = []

    # The symmetric polarization of c^T F xi has .5 F in each orientation.
    for axial in range(2):
        row: list[mp.mpf] = []
        for angular in range(3):
            coefficient = abs(_matrix_value(fg, axial, angular)) / 2
            pair[axial][2 + angular] += coefficient
            pair[2 + angular][axial] += coefficient
            row.append(coefficient)
        axial_angular.append(row)

    # ||Am^-2||_C1 <= 40 / Pstar^2 is the supplied paper bound.
    am_inverse_square_bound = mp.mpf(40) / (pstar * pstar) if amplitude_bound is None else amplitude_bound
    for left in range(2):
        row = []
        for right in range(2):
            coefficient = abs(_matrix_value(gg, left, right)) * am_inverse_square_bound
            pair[left][right] += coefficient
            row.append(coefficient)
        axial_axial.append(row)

    for left in range(3):
        row = []
        for right in range(3):
            coefficient = (
                abs(_matrix_value(ff, left, right)) / 2
                + abs(_matrix_value(ff_over_x, left, right)) / 2
            )
            pair[2 + left][2 + right] += coefficient
            row.append(coefficient)
        angular_angular.append(row)

    contribution_receipt = {
        "axial_angular_half_fg": axial_angular,
        "axial_axial_gg_times_amplitude_bound": axial_axial,
        "angular_angular_half_abs_ff_plus_ff_over_x": angular_angular,
        "Am_inverse_square_norm_bound": am_inverse_square_bound,
    }
    return pair, contribution_receipt


def compute_majorant(
    moment_map: Any,
    e_bound: Any,
    Pstar: Any,
    degree: int = 3,
    *, amplitude_inverse_square_bound=None, norm_label="C1",
) -> dict[str, Any]:
    """Compute a conditional uniform C1 response and tail majorant.

    ``e_bound`` must be supplied as a declared uniform C1 norm of the five
    defects.  The function does not infer or verify that norm from a local
    ``Z`` sample.  Raw MP values are returned under ``raw``; ``metadata`` is
    string/bool/list based and can be serialized directly.
    """

    if not hasattr(moment_map, "inverse_l1_norm"):
        raise TypeError("moment_map must expose inverse_l1_norm")
    precision = int(getattr(moment_map, "precision", 100))
    degree = _integer(degree, "degree")
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    with mp.workdps(max(80, precision + 20)):
        e = _mp(e_bound, "e_bound")
        pstar = _mp(Pstar, "Pstar")
        if e < 0:
            raise ValueError("e_bound must be nonnegative")
        if pstar < 1:
            raise ValueError("Pstar must be at least 1")
        ca = max(mp.mpf(1), abs(_mp(moment_map.inverse_l1_norm, "inverse_l1_norm")))
        amplitude_bound = mp.mpf(40)/(pstar*pstar) if amplitude_inverse_square_bound is None else _mp(amplitude_inverse_square_bound, 'amplitude_inverse_square_bound')
        if amplitude_bound < 0:
            raise ValueError('amplitude_inverse_square_bound must be nonnegative')
        pair_bounds, contribution_receipt = _quadratic_pair_bounds(moment_map, pstar, amplitude_bound)
        # Each pair entry already sums the absolute output-row coefficients
        # for one ordered (j,k) coefficient pair.  The maximum pair entry is
        # the l1 bilinear operator bound.  Keep a floor of one so all
        # downstream inequalities have a positive scale.
        cq_raw = max(
            (value for row in pair_bounds for value in row),
            default=mp.mpf(0),
        )
        cq = max(mp.mpf(1), cq_raw)
        threshold = 1 / (8 * ca * ca * cq)
        condition_passed = e <= threshold
        radius = 2 * ca * e
        lipschitz = 4 * ca * ca * cq * e
        response_scale = ca * e
        quadratic_scale = ca * ca * cq * e
        next_term = mp.mpf(_catalan(degree)) * response_scale * quadratic_scale ** degree
        ratio_bound = 4 * quadratic_scale
        if ratio_bound < 1:
            tail_bound = next_term / (1 - ratio_bound)
            tail_finite = True
        else:
            tail_bound = mp.inf
            tail_finite = False

        raw = {
            "e_bound": e,
            "Pstar": pstar,
            "degree": degree,
            "CA": ca,
            "CQ_raw_max_pair_sum": cq_raw,
            "CQ": cq,
            "Am_inverse_square_norm_bound": amplitude_bound,
            "condition_threshold": threshold,
            "contraction_condition_passed": condition_passed,
            "contraction_radius": radius,
            "lipschitz_bound": lipschitz,
            "response_scale_CA_e": response_scale,
            "quadratic_scale_CA2_CQ_e": quadratic_scale,
            "formal_next_term_bound": next_term,
            "tail_ratio_bound": ratio_bound,
            "tail_bound": tail_bound,
            "tail_bound_finite": tail_finite,
            "pair_bounds": pair_bounds,
        }
        metadata = {
            "norm_label": norm_label,
            "degree": degree,
            "CA": _nstr(ca),
            "CQ_raw_max_pair_sum": _nstr(cq_raw),
            "CQ": _nstr(cq),
            "Pstar": _nstr(pstar),
            "e_bound": _nstr(e),
            "condition_threshold": _nstr(threshold),
            "contraction_condition_passed": bool(condition_passed),
            "contraction_radius": _nstr(radius),
            "lipschitz_bound": _nstr(lipschitz),
            "formal_next_term_bound": _nstr(next_term),
            "tail_ratio_bound": _nstr(ratio_bound),
            "tail_bound": _nstr(tail_bound),
            "tail_bound_finite": bool(tail_finite),
            "pair_bounds": [[_nstr(value) for value in row] for row in pair_bounds],
            "contribution_receipt": {
                key: (
                    _nstr(value)
                    if isinstance(value, mp.mpf)
                    else [[_nstr(item) for item in row] for row in value]
                )
                for key, value in contribution_receipt.items()
            },
            "supplied_uniform_e_not_verified": True,
            "analytic_bump_norm_enclosed": False,
            "quadrature_enclosed": bool(getattr(moment_map, "quadrature_enclosed", False)),
            "actual_input_closure": False,
            "global_field_installed": False,
            "convergent_response_certified": False,
            "temporal_recursion_certified": False,
            "smallness_condition_is_conditional": True,
            "paper_pressure_bound": "||Am^-2||_C1 <= 40/Pstar^2" if amplitude_inverse_square_bound is None else "explicit supplied Banach algebra norm bound",
        }
        # ``raw_mp`` is an explicit alias for callers that distinguish the
        # arbitrary-precision payload from the serializable receipt.
        if amplitude_inverse_square_bound is None:
            raw['Am_inverse_square_C1_bound'] = amplitude_bound
            contribution_receipt['Am_inverse_square_C1_bound'] = amplitude_bound
            metadata['contribution_receipt']['Am_inverse_square_C1_bound'] = _nstr(amplitude_bound)
            metadata['contribution_receipt']['axial_axial_gg_times_40_over_Pstar2'] = metadata['contribution_receipt']['axial_axial_gg_times_amplitude_bound']
        return {"raw": raw, "raw_mp": raw, "metadata": metadata}


def _smoke() -> dict[str, Any]:
    try:
        from .lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
    except (ImportError, ValueError):
        from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap  # type: ignore
    with mp.workdps(90):
        moment_map = FiveBumpMomentMap(precision=90, order=48)
        first = compute_majorant(moment_map, 0, 1, degree=3)
        bound = first["raw"]["condition_threshold"]
        boundary = compute_majorant(moment_map, bound, 1, degree=3)
        failing = compute_majorant(moment_map, 2 * bound, 1, degree=3)
        return {
            "zero_passed": first["raw"]["contraction_condition_passed"],
            "boundary_passed": boundary["raw"]["contraction_condition_passed"],
            "boundary_lipschitz": mp.nstr(boundary["raw"]["lipschitz_bound"], 20),
            "failing_condition": failing["raw"]["contraction_condition_passed"],
            "zero_tail": mp.nstr(first["raw"]["tail_bound"], 10),
            "conditional_only": first["metadata"]["supplied_uniform_e_not_verified"],
        }


if __name__ == "__main__":
    print(json.dumps(_smoke(), indent=2))


__all__ = ["compute_majorant"]
