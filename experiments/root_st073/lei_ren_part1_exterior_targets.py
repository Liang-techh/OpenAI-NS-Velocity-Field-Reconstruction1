"""Transport Part I terminal moments inward through a pure-swirl exterior.

Source: Lei--Ren, arXiv:2609.35406v1, (1.2), (2.21)--(2.22).
The supplied F must be the actual exterior, with Uz=0, and must agree with
the supplied heat exterior from R_outer to infinity. This module does not
construct the unknown interior or certify its compatibility.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss

from lei_ren_part1_heat_moments import heat_tail_moments


def exterior_targets(F: Callable, heat, R_inner: float, R_outer: float,
                     Z: float, *, n: int = 64, P0: float | None = None) -> dict:
    """Return axis-to-inner moments required by the actual terminal data.

    Pressure is returned as P(R_inner). Its cumulative moment is
    P(R_inner)-P0, so a numerical pressure-moment target is provided only
    when the SAME candidate's axis pressure P0 has been supplied.
    Quadratic target is minus the entire exterior quadratic integral.
    """
    ri, ro, z = float(R_inner), float(R_outer), float(Z)
    if not all(math.isfinite(v) for v in (ri, ro, z)):
        raise ValueError("coordinates must be finite")
    if not 0 < ri <= ro or abs(z) > 1:
        raise ValueError("require 0 < R_inner <= R_outer and |Z| <= 1")
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 8:
        raise ValueError("n must be an integer at least 8")
    if P0 is not None and not math.isfinite(float(P0)):
        raise ValueError("P0 must be finite")
    nodes, weights = leggauss(int(n))
    width = math.log(ro / ri)
    ys = 0.5 * width * (nodes + 1)
    radii = ro * np.exp(-ys)
    fs = np.array([float(F(float(r), z)) for r in radii])
    if not np.all(np.isfinite(fs)):
        raise ValueError("exterior F returned a nonfinite value")
    # dR=R dy for the inward logarithmic integration.
    w = 0.5 * width * weights
    angular = float(np.dot(w, 2 * radii**2 * fs))
    quadratic = float(np.dot(w, -radii**2 * fs**2))
    pressure = float(np.dot(w, radii * fs**2))
    tail = heat_tail_moments(ro, z, heat, series_order=12)
    p_inner = -float(tail["pressure_tail"]) - pressure
    return {
        "R_inner": ri, "R_outer": ro, "Z": z,
        "angular_target": float(tail["angular_target"]) - angular,
        "axial_target": 0.0, "mixed_target": 0.0,
        "quadratic_target": -float(tail["quadratic_tail"]) - quadratic,
        "pressure_at_inner": p_inner,
        "pressure_moment_target": None if P0 is None else p_inner - float(P0),
        "axis_pressure": None if P0 is None else float(P0),
        "finite_exterior_angular": angular,
        "finite_exterior_quadratic": quadratic,
        "finite_exterior_pressure": pressure,
        "pure_swirl_exterior_assumed": True,
        "interior_closure_validated": False,
    }
